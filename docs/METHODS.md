# Higher-order co-occurrence analysis

## Input and observation unit

`Computing co-occurrence_Single.py` analyzes one species-by-sample abundance table. Each column is an ecological sample and each row is a species with a unique identifier. At least two species and two samples are required. Species identifiers must match model filename stems when the results are used for SMETANA.

Define sample inclusion and quality-control rules before analysis. Habitat mixtures, sequencing batches and repeated measurements can affect associations. A table with 48 columns does not necessarily represent 48 independent biological replicates. Failed samples should be excluded upstream using a documented rule.

## Abundance units and presence

- `counts`: normalize each complete supplied sample column by its total, before any model-based species restriction. An incomplete count table supplies an incomplete denominator.
- `tpm`: divide species-level TPM by 1,000,000. Subsets are not renormalized. Each value must be at most 1,000,000 and each column total at most 1,000,001, allowing 1 TPM for rounding.
- `relative`: use fractions in [0, 1], with column totals at most 1, allowing a tolerance of 1e-8. Subsets are not renormalized. Convert percentages to fractions before input.
- `presence-absence`: accept only 0 and 1; the abundance cutoff is not reapplied.

For abundance input, a species is present only when its relative abundance is **strictly greater than** `--abundance-cutoff`. The default 0.001 is 0.1%, corresponding to more than 1,000 TPM; 0.0001 corresponds to more than 100 TPM. These thresholds define detection and do not classify communities as competitive or cooperative.

TPM input assumes species-level values with a consistent denominator. Taxonomic aggregation and gene-to-species mapping must be resolved before analysis. Do not normalize a partial TPM table as counts.

All sample columns remain in the analysis denominator, including zero-total and low-richness samples. All species in the defined analysis universe remain eligible, including species with zero prevalence. `--models-dir` optionally restricts that universe to exact lowercase `.xml` filename stems, after abundance conversion and thresholding. Excluded IDs are recorded. The restriction changes the hypothesis family and can affect adjusted p-values.

## Joint occurrence and independence null

Let `n` be the retained sample count, `p_i` the prevalence of species `i`, and `G` a group of species:

```text
observed(G) = number of samples containing every member of G
p0(G)       = product of member prevalences p_i
expected(G) = n * p0(G)
X_G         ~ Binomial(n, p0(G))
```

The default one-sided enrichment test includes the observed count:

```text
p_value(G) = P(X_G >= observed(G))
           = scipy.stats.binom.sf(observed(G) - 1, n, p0(G))
```

The survival function avoids subtracting a cumulative probability close to one. With the inclusive definition, zero observations give a p-value of one. `--pvalue-tail strict` instead selects `P(X_G > observed(G))`; record the selected definition when reporting a run.

Marginal prevalences are estimated from the same dataset. This is a plug-in independence null, not a covariate-adjusted ecological model or a permutation test preserving sample-specific richness. `fold` is observed divided by expected; support and significance requirements still apply when counts are zero or enrichment is large.

## Multiple testing and selection

The search begins with every unique species pair. Each subsequent order evaluates every unique candidate generated from the propagated population. Benjamini-Hochberg adjustment uses **all candidates evaluated at that order**, including groups with zero support. Support and fold filtering occur after adjustment.

A qualifying group satisfies all three criteria:

- observed count at least `--min-samples`, default 10;
- observed count at least `--min-fold` times the expected count, default 2;
- BH-adjusted p-value at most `--fdr-cutoff`, default 0.05.

`q_value` denotes the BH-adjusted p-value. `passes_fdr` records only its threshold, while `significant` records all three criteria. With 48 samples, support 10 requires joint presence in at least 10/48, approximately 20.8%, of samples. Floating-point comparisons allow a relative tolerance of 1e-12 at the inclusive fold boundary.

## Ranking and higher-order search

At each order, qualifying groups are ranked by descending joint-occurrence count. Lexical species order resolves ties. At most `--top-n` groups, default 1,000 **per size**, are exported.

The population for expansion is selected separately from the exported top list. Up to `--population-size`, default 10,000, qualifying groups are sampled without replacement. The default `--selection-weight frequency` uses weights proportional to observed frequency; `sqrt-frequency` uses its square root. If all qualifying groups fit within the population limit, all propagate. A fixed seed and deterministic candidate order control reproducibility.

Each propagated group is extended by every other eligible species, and equivalent groups are deduplicated before testing. A candidate requires a propagated significant subset; it does not require all member pairs to be significant.

The default maximum group size is 5 and is configurable with `--max-combo-size`. Search stops at that size, when no significant groups remain, or when no further extensions exist. This heuristic can miss groups without a propagated significant subset. A failure to discover a group does not establish absence of an association.

`--max-candidates`, default 1,000,000 per order, raises an error when exceeded. Candidates and statistics are held in memory across orders; the limit does not guarantee that its maximum is feasible on every computer. Restrict a species universe based on a defined study scope, or reduce the propagated population with its search implications documented. Record seed sensitivity analyses without selecting only favorable runs.

## Outputs and reproducibility

- `communities.tsv`: selected groups across completed orders, in headerless two-column SMETANA format.
- `size_N_communities.tsv`: selected groups of size N in the same format.
- `size_N_statistics.tsv`: all tested candidates, with columns `community_id`, `size`, `species`, `observed`, `expected`, `fold`, `p_value`, `q_value`, `passes_fdr`, `significant` and `selected`. The species field is a JSON array; the community ID is blank unless the group is exported.
- `run_summary.json`: input and implementation hashes, software versions, sample/species identifiers, settings, preprocessing metadata, order counts and stopping reason.

Use a new or empty output directory. An empty selected-community file is a valid result when no groups pass; no random replacement or automatic threshold relaxation is performed. Do not pass an empty community file to SMETANA.

The synthetic 48-sample example has A, B and C present together in 12 samples. Each pair has observed count 12, expected count 3 and fold 4. The triple has expected count 0.75 and fold 16. These known values test computation; they are not biological observations.

## Metabolic analysis and interpretation

Infer groups across the defined sample population, then simulate the selected groups against a consistent genome-scale metabolic model library. Each abundance identifier must refer to the same biological entity as its model. Resolve strains, taxonomic aggregation and many-to-one mappings explicitly.

Use `Smetana_FWQ.py` with `--single` for one model library and `--communities` for the inferred community file. The community path is relative to that library unless it is absolute. Global mode requests MIP/MRO; detailed mode requests donor-receiver-metabolite scores. Exact commands and directory layout are in [README.md](../README.md).

With fixed models and environmental settings, each unique group needs one simulation. Sample-specific models or media require corresponding metadata and separate outputs. A group selected across a dataset need not occur in every sample; all referenced models must exist wherever it is simulated.

Co-occurrence is an association, not evidence of direct interaction. SMETANA scores describe model-based metabolic potentials. Comparing cooperative and competitive tendencies requires consistent simulation settings and appropriate same-size random controls. The optional random-community utility does not automatically construct a matched control analysis.

Per-order BH adjustment does not establish a global FDR bound for the entire adaptive search. Overlapping groups, fitted marginal probabilities, habitat differences and non-independent samples require care in interpretation. Synthetic and mocked tests do not validate biological data, model quality, solver configuration or media assumptions.

## References

The higher-order search strategy is based on Machado et al. (2021), *Polarization of microbial communities between competitive and cooperative metabolism*, Nature Ecology & Evolution 5, 195-203. [Article](https://doi.org/10.1038/s41559-020-01353-4), [open manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC7610595/).

[HiOrCo](https://github.com/cdanielmachado/HiOrCo/tree/ee9af1f1911f6313376247352051f25e09b2fd0a) provides a reference method. The statistical definitions, preprocessing rules and configurable parameters documented here specify this implementation; they do not guarantee reproduction of the paper's datasets or results.
