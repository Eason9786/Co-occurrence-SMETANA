# Higher-order co-occurrence analysis and SMETANA

Infer enriched microbial species groups across samples, then evaluate their metabolic interaction potential with [SMETANA](https://github.com/cdanielmachado/smetana).

[中文说明](README_zh-CN.md) · [Installation](INSTALL.md) · [Method and limitations](docs/METHODS.md)

## Workflow and input

The workflow is: define and quality-control the sample population, threshold abundances, infer higher-order co-occurrence across samples, match species to metabolic models, and simulate the selected groups. `Computing co-occurrence_Single.py` processes one multi-sample abundance table using an iterative strategy based on [Machado et al. (2021)](https://doi.org/10.1038/s41559-020-01353-4). `Smetana_FWQ.py` runs the selected groups against a consistent SBML model library. Co-occurrence describes association; metabolic models provide estimates of interaction potential.

The abundance file has a header, species as rows and samples as columns. The first column contains unique species IDs:

```text
species_id	sample_01	sample_02	sample_03
species_A	100	0	50
species_B	120	0	60
species_C	0	200	90
```

This three-sample illustration explains the format, but cannot satisfy the default minimum of 10 co-occurring samples. Values must be complete, numeric, finite and nonnegative. Duplicate or empty species/sample IDs are rejected. TSV is the default; use `--sep comma` for CSV.

- `--input-kind counts` (default): divide each sample column by its full supplied total, then threshold. Supply the complete count table relevant to that denominator.
- `--input-kind tpm`: divide values by 1,000,000 without renormalizing a species subset. Column totals must not exceed 1,000,000, apart from a small rounding tolerance.
- `--input-kind relative`: supply fractions from 0 to 1, with column sums at most 1; a partial table is not renormalized. Convert percentages to fractions first.
- `--input-kind presence-absence`: supply exactly 0/1 values using a detection threshold defined upstream.

Abundance must be **strictly greater than** the cutoff to count as present. All sample columns are retained, including zero-total columns; remove failed samples upstream using a documented quality-control rule. Use `tpm` for a partial TPM table to preserve its original abundance denominator.

## Installation

Use Python 3.10 or newer. For inference and SMETANA with the supplied Conda environment:

```bash
conda env create -f environment_smetana.yml
conda activate smetana_env
```

Configure a supported solver before metabolic simulation; see [INSTALL.md](INSTALL.md). For inference alone in an existing Python environment, install `requirements.txt` with `python -m pip install -r requirements.txt`.

## Quick start with 48 synthetic TPM samples

Install dependencies using [INSTALL.md](INSTALL.md), then run from the repository root:

```bash
python "Computing co-occurrence_Single.py" "examples/synthetic_tpm_48.tsv" --input-kind tpm --output "results/demo"
```

This artificial dataset is expected to identify A/B, A/C, B/C and A/B/C. Their observed count is 12 of 48 samples; expected counts are 3 for each pair and 0.75 for the triple. It tests computational behavior only and includes no usable SBML models.

For your own TPM data and model library:

```text
data/
  abundance_tpm.tsv
  models/
    species_A.xml
    species_B.xml
    ...
```

```bash
python "Computing co-occurrence_Single.py" "data/abundance_tpm.tsv" --input-kind tpm --models-dir "data/models" --abundance-cutoff 0.001 --max-combo-size 5 --top-n 1000 --population-size 10000 --seed 42 --output "results/tpm_0p1"
```

`--models-dir` is optional. It limits the tested species universe to exact model filename stems, after abundance conversion and thresholding. It changes the hypothesis family and reports excluded species. Use lowercase `.xml` filenames. Without it, all input species are tested; ensure every exported member has a suitable model before simulation.

Use a new or empty output directory for each run; the default is `cooccurrence_results` if `--output` is omitted. A cutoff of `0.001` is 0.1%, or strictly more than 1,000 TPM. For sensitivity to 0.01%, rerun with `--abundance-cutoff 0.0001` (strictly more than 100 TPM) and a different output directory. Those cutoffs do not automatically label groups competitive or cooperative.

## Statistical defaults

For each evaluated group, count samples where every member is present. Under an independent-occurrence null, the expected count is the number of samples multiplied by the product of member prevalences. A one-sided binomial enrichment test is followed by Benjamini-Hochberg correction across **all evaluated groups of that size**, before support or fold filtering.

A significant group must have at least 10 joint observations, at least twice its expected count, and an adjusted p-value at most 0.05. With 48 samples, the minimum support requires presence in at least 20.8% of the samples. The top 1,000 significant groups by frequency are exported **per size**. Up to 10,000 significant groups, sampled without replacement with weights proportional to observed frequency, seed the next order; each is extended by every eligible additional species.

Change criteria with `--min-samples`, `--min-fold`, `--fdr-cutoff`, `--top-n` and `--population-size`. The default maximum size is 5 and can be increased. Higher-order search is heuristic, not exhaustive. `--max-candidates` defaults to 1,000,000 and fails explicitly if a family exceeds the limit, rather than silently subsampling its tests. Candidates and statistics are held in memory across orders; the cap is a resource guard, not a streaming/scalability guarantee. No hits is a valid result; thresholds are not automatically relaxed and random communities are not substituted.

The default binomial enrichment tail is `P(X >= observed)` and the default propagation weight is observed frequency. `--pvalue-tail strict` selects `P(X > observed)`; `--selection-weight sqrt-frequency` selects square-root frequency weights. The exact definitions, assumptions and relation to [HiOrCo](https://github.com/cdanielmachado/HiOrCo) are documented in [METHODS.md](docs/METHODS.md).

## Outputs

- `communities.tsv`: selected groups across completed sizes, in SMETANA's two-column, **headerless** format.
- `size_N_communities.tsv`: selected groups for one size, in the same format.
- `size_N_statistics.tsv`: every evaluated candidate, JSON member IDs, observed/expected counts, enrichment fold, raw and adjusted p-values, and selection flags.
- `run_summary.json`: input hash, settings, processing metadata, per-order counts and stopping reason.

Example community definition:

```text
community_2_0001	species_A
community_2_0001	species_B
```

An empty `communities.tsv` with a completed summary means no groups passed. Inspect the statistics before simulation; the runner rejects an empty community file. Inference outputs use TSV/JSON.

## Run SMETANA

Use the inferred group file once with a model library containing all referenced IDs. For the directory layout above:

```bash
conda activate smetana_env
python Smetana_FWQ.py "data/models" --single --communities "../../results/tpm_0p1/communities.tsv" --solver gurobi --mode global --dry-run
python Smetana_FWQ.py "data/models" --single --communities "../../results/tpm_0p1/communities.tsv" --solver gurobi --mode global
```

`--communities` is resolved relative to the model folder, hence `../../results/...`; an absolute path also works. `--single` means **one model library**, not one ecological sample. Global mode requests MIP/MRO; use `--mode detailed` for donor-receiver-metabolite scores. Alternatively add `--conda-env smetana_env` to launch SMETANA through Conda.

`global.tsv` or `detailed.tsv` and `smetana_run.log` are written in the model folder. Repeated runs can overwrite these files; retain independent model-folder copies for comparisons. Preflight checks model names and community membership, not SBML validity or solver licensing.

Inference requires multiple ecological samples. With fixed GEMs and simulation settings, simulate each selected group once. Sample-specific model variants or media require an explicit design and separate outputs. The runner supports sample-folder and group-folder batches for such prepared inputs; see `--help`.

## Validation, limitations and citation

```bash
python -m unittest discover -s tests -v
```

Tests cover statistical calculations, search behavior, input checks, reproducibility and command handling. Synthetic and mocked tests do not establish biological validity or solver compatibility. Shared habitats, repeated measurements, sequencing batches and model coverage can affect associations. Per-order BH correction is not a global FDR guarantee for the complete adaptive search. Claims of polarization require consistent MIP/MRO comparisons and suitable same-size controls.

`legacy/` contains optional dataset-specific preprocessing helpers and `random_subcommunities.py` for artificial control candidates. Match control group sizes and counts explicitly when designing a comparison. Optional model reconstruction and solver configuration are described in [INSTALL.md](INSTALL.md). Environment files specify Python 3.10 and SMETANA 1.2.1; export the resolved environment after installation for reproducibility.

Workflow code uses the [MIT License](LICENSE). External tools retain their own licenses; no solver key or private dataset is included. Cite this repository and the methods/tools used:

- Machado et al. (2021). *Polarization of microbial communities between competitive and cooperative metabolism.* Nature Ecology & Evolution 5, 195-203. [Article](https://doi.org/10.1038/s41559-020-01353-4), [HiOrCo](https://github.com/cdanielmachado/HiOrCo).
- Machado et al. (2018). *Fast automated reconstruction of genome-scale metabolic models for microbial species and communities.* Nucleic Acids Research 46, 7542-7553. [CarveMe](https://github.com/cdanielmachado/carveme).
- Zelezniak et al. (2015). *Metabolic dependencies drive species co-occurrence in diverse microbial communities.* PNAS 112, 6449-6454. [SMETANA](https://github.com/cdanielmachado/smetana).

Questions and contributions: [GitHub issues](https://github.com/Eason9786/Co-occurrence-SMETANA/issues).
