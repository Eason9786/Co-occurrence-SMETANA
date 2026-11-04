# Statistical method and framework review

## Sources and intended fidelity

This implementation follows the higher-order search strategy in **Computing co-occurrence**, Machado et al., *Polarization of microbial communities between competitive and cooperative metabolism*, Nature Ecology & Evolution 5, 195-203 (2021), [DOI](https://doi.org/10.1038/s41559-020-01353-4), [open manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC7610595/). In the supplied journal PDF, this section spans pages 201-202 (PDF pages 7-8).

We also inspected [HiOrCo's method](https://github.com/cdanielmachado/HiOrCo/blob/ee9af1f1911f6313376247352051f25e09b2fd0a/hiorco/method.py), pinned to commit `ee9af1f1911f6313376247352051f25e09b2fd0a`. Workflow code here is independently implemented. Paper principles and reference-code conventions are distinguished below; no claim is made to recreate published communities or ecological conclusions.

## 1. Define the observations and units

An observation is an ecological sample, not a model file. Rows identify species and columns identify samples. A model folder alone supplies neither marginal occurrence probabilities nor evidence of enrichment. The retained filename `Computing co-occurrence_Single.py` denotes one multi-sample dataset, not one ecological sample.

Define the sample population before looking at significance. Mixing habitats, treatments or sequencing batches can create associations driven by shared environments. Repeated observations from the same subject/site also undermine independent-sample assumptions. The script does not fit covariate-adjusted, repeated-measures or habitat-stratified nulls. A dataset of 48 columns need not contain 48 independent biological replicates.

Counts are normalized by each **complete supplied column total**, before model filtering. An incomplete count table has an incomplete denominator. Supply the full relevant table, or already computed fractions with `--input-kind relative`. Relative inputs are not rescaled, so sums below one are allowed for an already selected species subset. Percentages must first be converted to fractions.

For `--input-kind tpm`, values are divided by one million without renormalization. Each column must total at most 1,000,000 TPM, allowing a tolerance of 1 TPM for rounding. This mode preserves the original denominator when only selected taxa appear in the table. It assumes that the supplied values are valid species-level TPM on a consistent denominator; the script does not infer how gene-level TPM was aggregated into taxa or correct an upstream aggregation error. Using `counts` for a partial TPM table would incorrectly rescale those values.

Abundance-based presence is `relative_abundance > abundance_cutoff`. The default `0.001` is 0.1%, or strictly more than 1,000 TPM. `0.0001` is 0.01%, or strictly more than 100 TPM. Presence-absence input accepts only 0/1 values and assumes thresholding was already performed upstream. These cutoffs do not themselves classify communities as competitive or cooperative.

All sample columns are retained, including zero-total samples and columns with fewer than two eligible present species. They contribute absences to the denominator. Failed samples should be removed using a documented prior quality-control rule. Species with zero prevalence remain in the specified tested universe. This differs from the inspected HiOrCo preprocessing and makes the sample denominator and hypothesis family explicit.

`--models-dir` optionally restricts species to exact lowercase `.xml` filename stems, after abundance conversion and thresholding. Exclusions are recorded. This defines a different hypothesis family and can change adjusted p-values; it does not recover organisms lacking models or validate model contents. Apply model coverage and taxonomic mapping consistently when comparing runs.

## 2. Count joint occurrence and specify the null

For `n` retained samples, let `p_i` be the fraction in which species `i` is present. For group `G`:

```text
observed(G) = samples where all members of G are present
p0(G)       = product of p_i over members of G
expected(G) = n * p0(G)
X_G         ~ Binomial(n, p0(G)) under independent species occurrence
```

The default enrichment p-value includes the observed count:

```text
p_value(G) = P(X_G >= observed(G))
           = scipy.stats.binom.sf(observed(G) - 1, n, p0(G))
```

A survival function avoids subtracting a cumulative probability close to one. Under this definition a group with zero observations has a p-value of one. Marginal prevalences are estimated from the data; this is a plug-in independence model, not a mechanistic ecological null or a permutation model preserving all sample characteristics.

`fold` is observed divided by expected. Zero observed and expected counts do not supply evidence of enrichment; the support and significance requirements still fail. A large fold is not sufficient by itself.

In the included artificial 48-sample dataset, A/B/C each occur in the same 12 samples. Each pair therefore has expected count `48*(12/48)^2 = 3`, observed count 12 and fold 4. Their triple has expected count `48*(12/48)^3 = 0.75`, observed count 12 and fold 16. These are deliberately simple numerical checks, not estimated effects from the user's data.

## 3. Correct within each evaluated order

Size 2 evaluates every unique pair in the defined species universe. At later sizes, the tested family is every unique candidate generated from the previous propagated population. Deduplication occurs before testing.

Benjamini-Hochberg adjustment uses **all p-values in that evaluated family**, including zero-support groups. Support and fold filtering happen afterwards; removing unpromising groups before correction would shrink the correction family and alter the result.

A group must satisfy all three conditions:

- `observed >= min_samples`, default 10;
- `observed >= min_fold * expected`, default 2;
- `q_value <= fdr_cutoff`, default 0.05.

Here `q_value` is the Benjamini-Hochberg adjusted p-value, not a separate Storey q-value estimate. `passes_fdr` is reported independently of `significant`, because adjusted significance can pass while the support or fold requirement fails. With 48 samples, support 10 corresponds to at least 20.8% prevalence for a whole group. Results with fewer observations are not automatically rescued by lowering the threshold.

## 4. Rank, propagate and extend

At each size, significant groups are ranked by observed joint-occurrence count, with stable lexical ordering for ties. Up to `top_n=1000` groups **per size** are saved. There is no shared total budget and no random filling if fewer groups pass.

The population for the next order is selected separately from this saved top list. At most `population_size=10000` significant groups are drawn without replacement, using roulette-style weights proportional to observed frequency. When fewer significant groups exist, all propagate. A fixed seed controls this selection with deterministic candidate ordering.

Every propagated group is extended by every other eligible species; equivalent results are deduplicated. A higher-order candidate must have a propagated significant subset, but does **not** require every pair within it to have passed. Requiring a fully significant pairwise clique would implement a different search.

Search stops at the maximum size, when no significant groups remain, or when no further extensions exist. The default maximum size 5 is a configurable engineering starting point, not the paper's largest explored group. Higher-order search is heuristic and can miss significant groups without a propagated subset. Failure to find a group is not proof of no association.

The candidate cap, default 1,000,000, raises an error instead of silently reducing the tested family. Candidates and per-candidate scores are materialized, with results retained across completed orders; memory use therefore grows with evaluated candidates across the run. The cap does not make the implementation streaming or guarantee feasibility at its maximum. Raising it increases resource use. Restricting a predefined species universe or reducing the propagated population changes the inferential/search scope and should be recorded. Separate runs with different seeds can assess search stability; record all choices rather than choosing only the most favorable run.

## 5. Differences from the inspected HiOrCo source

The paper specifies a binomial test and frequency-proportional selection but leaves some implementation conventions unstated. The inspected source contains these details:

- Its `1 - binom.cdf(observed, n, p0)` is `P(X > observed)`. This repository defaults to the standard inclusive enrichment tail `P(X >= observed)`. `--pvalue-tail strict` exposes the reference tail using a stable survival function. Strict tails can give zero probabilities at boundary observations; they are not a general improvement over the inclusive test.
- It uses `sqrt(observed)` for propagation weights. This repository defaults to raw frequency, following the paper's text; `--selection-weight sqrt-frequency` exposes the source convention. Both rank final top lists identically, but may propagate different groups.
- It filters species by marginal support and removes low-richness samples during preprocessing. This implementation retains all supplied samples and eligible species. A fixed denominator avoids making the observation population depend on model-covered richness.
- Its callable and command-line defaults differ from the paper's stated 1,000 saved and 10,000 propagated groups. This release uses those paper-level values, with a modest configurable maximum size.

Both reference-convention flags reproduce **only those two conventions**. They do not recreate upstream preprocessing, random-number sequences, tie ordering, EMP data processing or published outputs. Compare methods only with their input populations, software versions and settings specified.

## 6. Output and inferential limits

`size_N_statistics.tsv` contains every evaluated candidate. Columns are `community_id` (blank unless saved), `size`, `species` (JSON array), `observed`, `expected`, `fold`, `p_value`, `q_value`, `passes_fdr`, `significant` and `selected`. `selected` identifies saved top results, not all groups eligible for propagation.

`communities.tsv` combines selected groups from completed orders into SMETANA's headerless two-column definition. `size_N_communities.tsv` is the corresponding per-order file. `run_summary.json` records the input and implementation hashes, Python/NumPy/SciPy versions, sample and analyzed species IDs, settings, processing metadata, order trace and stopping reason. No hits is a valid successful result, including when support cannot reach the configured minimum. Existing nonempty output directories are rejected to prevent stale/new result mixtures.

Adjusted p-values summarize the tested family at each stage. Adaptive construction of later families, dependence among overlapping groups, estimated marginal prevalences and ecological sampling assumptions mean this workflow does **not** establish a global FDR bound across all possible groups or stages. Causal or confirmatory ecological claims need additional study design and validation.

## 7. Connect statistical groups to metabolic models

The coherent default framework is to infer co-occurrence once across a specified sample population, then simulate selected groups using a consistent GEM library. Each abundance ID must denote the same biological entity as its model filename. Resolve taxonomic aggregation, multiple strains and many-to-one mapping upstream; silent renaming would obscure those assumptions.

With fixed GEMs and unchanged SMETANA settings, each group need only be simulated once. Repeating the same model set in different sample folders does not make predictions sample-specific. Sample-dependent GEMs or media are possible, but require matched metadata and separate outputs. A globally selected group need not occur in every original sample, so never blindly assign the global file to sample folders with missing members.

`Smetana_FWQ.py --single` accepts a shared model library and an explicit group file. Global mode requests MIP/MRO, and detailed mode requests interaction scores. These describe model-based metabolic potentials rather than measured ecological competition/cooperation. To study polarization as in the paper, compare MIP/MRO under consistent conditions and use appropriately matched same-size random controls. The archived random-combination script alone does not supply this full control design or infer a cooperative/competitive label.

The synthetic example and unit tests verify computational behavior, not the user's biological data, model quality, solver environment, media assumptions or published metabolic landscape.
