# MetaCoNet: Co-occurrence–SMETANA workflow

Scripts for preparing candidate microbial subcommunities from genome-scale metabolic models and running [SMETANA](https://github.com/cdanielmachado/smetana) in batches.

[中文说明](README_zh-CN.md) · [Installation](INSTALL.md) · [Changes from the local scripts](CHANGELOG.md)

## Scope

This repository shares the practical model-combination and simulation workflow of MetaCoNet (Metabolic Co-occurrence Network). Starting with one SBML XML model per species in a sample or module folder, it samples candidate communities of 2–5 species, writes SMETANA-compatible community definitions, and runs global or detailed metabolic interaction analysis.

**The supplied community generator samples combinations within a single folder. It does not calculate cross-sample co-occurrence frequencies, binomial p-values, FDR/q-values, abundance thresholds, or statistically significant associations.** Any ecological network/module selection and abundance filtering must be performed upstream. Generated combinations are candidate communities, not evidence of mutualism or observed co-occurrence.

## Main files

- `Computing co-occurrence_Single.py`: reproducible sampling across community sizes, with a configurable total budget; exports TSV and optional Excel.
- `Smetana_FWQ.py`: validates community definitions, runs SMETANA per sample, saves logs, and reports failures.
- `environment_smetana.yml`: proposed Python 3.10 / SMETANA 1.2.1 environment, including the Gurobi Python interface.
- `environment_carveme.yml`: optional model-reconstruction environment with CarveMe and DIAMOND.
- `legacy/`: original dataset-specific FASTA preparation helpers, with limitations documented separately.
- `tests/`: checks for community sampling and batch execution; no solver required.

## Input layout

Provide your own valid SBML models; biological datasets, model files and solver licenses are not bundled.

```text
data/
  sample_01/
    species_A.xml
    species_B.xml
    species_C.xml
    communities.tsv       # generated in step 1
  sample_02/
    ...
```

The model filename without `.xml` is the species identifier. Use unique identifiers and lowercase `.xml` extensions. The scripts read XML files directly inside each sample folder; compressed models are not supported by these wrappers.

## Quick start

First follow [INSTALL.md](INSTALL.md) to install SMETANA and configure a solver. Run commands from the repository root. Quote the generator filename because it contains a space.

### 1. Generate communities for one sample

```bash
conda activate smetana_env
python "Computing co-occurrence_Single.py" "data/sample_01" --max-combo-size 5 --total-budget 1000 --seed 42 --output "data/sample_01/communities.tsv" --excel "data/sample_01/communities.xlsx"
```

Repeat for each sample to be simulated. The generator distributes the budget as evenly as possible across sizes 2 to `min(max_combo_size, species_count)`, redistributing capacity left over from smaller sizes. It returns at most the available number of distinct combinations. The same seed, model identifiers, arguments and Python version give reproducible output.

TSV files contain **two tab-separated columns without a header**, one community/member pair per row:

```text
community1	species_A
community1	species_B
community2	species_A
community2	species_C
```

Species IDs must match the model filenames. Excel output includes a `Community` / `Species` header for readability; do not pass the Excel file or a headered TSV directly to SMETANA.

### 2. Check inputs and commands

```bash
python Smetana_FWQ.py "data" --solver gurobi --dry-run
```

This checks model filenames and community membership and prints the command for each sample. It does not parse SBML, check solver licensing, or run simulations.

### 3. Run SMETANA

```bash
python Smetana_FWQ.py "data" --solver gurobi --mode global
```

Global mode requests SMETANA's MIP/MRO analysis. Detailed mode requests donor–receiver–metabolite interaction scores:

```bash
python Smetana_FWQ.py "data" --solver gurobi --mode detailed
```

Outputs are written inside each sample folder: `global.tsv` or `detailed.tsv`, plus `smetana_run.log`. Running again can overwrite these files, so use a separate copy of the input folders for runs you want to retain. Scores reflect the supplied models and SMETANA's environmental assumptions; they require biological interpretation.

For a single model folder, use `--single`:

```bash
python Smetana_FWQ.py "data/sample_01" --single --solver gurobi
```

For the original group/sample layout (`data/CK_XML/sample_01`, etc.):

```bash
python Smetana_FWQ.py "data" --groups CK_XML DL_XML SL_XML UL_XML --solver gurobi
```

By default the runner uses `smetana` from the active environment. Alternatively, launch SMETANA with Conda explicitly:

```bash
python Smetana_FWQ.py "data" --conda-env smetana_env --solver gurobi
```

`--communities` changes the default `communities.tsv` filename; an absolute path can be used for a shared file only when all referenced models exist in every sample. A failed sample does not prevent later samples from running, but the overall command exits with an error. See `--help` on either script for all options.

## Optional upstream reconstruction

Protein FASTA files can be reconstructed with CarveMe after installing its own environment and a supported solver:

```bash
conda activate carveme_env
carve -r "species_faa/*.faa" -o "data/sample_01" --solver gurobi
```

Use `--dna` for nucleotide input, not protein FASTA. Verify model quality, growth feasibility, identifier consistency, and biological suitability before downstream analysis. The optional `legacy/` scripts demonstrate dataset-specific sequence preparation and file selection; they do not implement the abundance filtering discussed in the original project outline.

## Validation and reproducibility

```bash
python -m unittest discover -s tests -v
```

Tests exercise sampling, TSV formatting, input validation, and subprocess handling. Batch execution tests use mocked commands; they do not establish biological validity or solver compatibility. The environment files are starting specifications, not fully locked or end-to-end benchmarked environments. Historical local notes used SMETANA 1.2.0 and Python 3.6; this release modernizes the wrappers (Python >=3.8) and proposes SMETANA 1.2.1 with Python 3.10. Existing validated environments can be selected with `--conda-env`; revalidate results when changing software or solver versions.


## License and citation

These workflow scripts are distributed under the [MIT License](LICENSE), following the supplied project README. SMETANA, CarveMe, solvers and their dependencies retain their own licenses. Gurobi requires an appropriate license for your usage; no key is included here.

A MetaCoNet publication citation has not yet been supplied. Until available, cite the repository and the tools used:

- Machado D., Andrejev S., Tramontano M., Patil K. R. (2018). *Fast automated reconstruction of genome-scale metabolic models for microbial species and communities.* Nucleic Acids Research, 46(15), 7542–7553. [CarveMe](https://github.com/cdanielmachado/carveme).
- Zelezniak A. et al. (2015). *Metabolic dependencies drive species co-occurrence in diverse microbial communities.* PNAS, 112(20), 6449–6454. [SMETANA citation](https://github.com/cdanielmachado/smetana#citing-smetana).

Questions and contributions: [GitHub issues](https://github.com/Eason9786/Co-occurrence-SMETANA/issues).
