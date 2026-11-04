# Installation / 安装说明

## Statistical inference only

Use Python 3.10 or newer with NumPy and SciPy. The inference script needs no SBML parser or optimization solver:

```bash
python -m pip install -r requirements.txt
python "Computing co-occurrence_Single.py" --help
python "Computing co-occurrence_Single.py" "examples/synthetic_tpm_48.tsv" --input-kind tpm --output "results/demo"
python -m unittest discover -s tests -v
```

`requirements.txt` also includes openpyxl for optional Excel output from the random-community utility. The statistical inference script writes TSV and JSON.

The example is synthetic and verifies computation only; it contains no usable SBML models. Use a new or empty output directory for a repeated run. TPM values require `--input-kind tpm`; they are divided by one million without rescaling a species subset. A relative-abundance cutoff of `0.001` therefore means strictly more than 1,000 TPM.

## Inference and SMETANA environment

With Conda installed and initialized:

```bash
conda env create -f environment_smetana.yml
conda activate smetana_env
python --version
python -c "import numpy, scipy; print(numpy.__version__, scipy.__version__)"
smetana --help
```

The environment specifies Python 3.10, inference dependencies, SMETANA 1.2.1 and Gurobi's Python interface. It is an installation specification rather than a fully locked environment. Record the resolved package and solver versions used for an analysis.

The wrappers can be used on Windows when `smetana` is a working command in the selected environment. Linux/WSL is convenient for the optional DIAMOND/CarveMe workflow. WSL needs Linux paths and environments installed within WSL, rather than Windows Conda paths.

## Configure a solver

SMETANA requires a supported solver. Configure Gurobi using the [vendor's instructions](https://docs.gurobi.com/projects/optimizer/en/current/reference/misc/licensing.html), then use `--solver gurobi`. Installing `gurobipy` does not grant an unrestricted license. Keep keys and license files outside the repository.

For an established CPLEX setup, install its Python interface in the SMETANA environment and pass `--solver cplex`. Follow the selected solver's version-compatibility requirements.

A basic Gurobi initialization check is:

```bash
python -c "import gurobipy as gp; m = gp.Model(); m.optimize()"
```

This checks only a trivial problem. Test a small group of valid metabolic models before a large run; a restricted license may not support the full community model.

## Connect inference to simulations

Prepare a multi-sample abundance table and model library with IDs matching model filename stems. For TPM input:

```bash
python "Computing co-occurrence_Single.py" "data/abundance_tpm.tsv" --input-kind tpm --models-dir "data/models" --output "results/tpm_0p1"
```

Check `run_summary.json` and the per-size statistics, and confirm that `communities.tsv` has selected groups. With the repository layout in [README.md](README.md):

```bash
python Smetana_FWQ.py "data/models" --single --communities "../../results/tpm_0p1/communities.tsv" --solver gurobi --dry-run
python Smetana_FWQ.py "data/models" --single --communities "../../results/tpm_0p1/communities.tsv" --solver gurobi --mode global
```

The community path is relative to the model folder, not the terminal's working directory; an absolute path avoids ambiguity. `--single` selects one model library. Inference still requires multiple ecological samples.

Activate the environment containing SMETANA first, or add `--conda-env ENV_NAME`. The latter uses `conda run`, so Conda must be available in the terminal. Preflight `--dry-run` checks names and memberships without parsing SBML or invoking the solver.

## Optional CarveMe reconstruction

```bash
conda env create -f environment_carveme.yml
conda activate carveme_env
carve --help
diamond version
```

CarveMe is separate to keep reconstruction dependencies independent. Its version is not locked; follow the [official installation guide](https://carveme.readthedocs.io/en/latest/installation.html), configure a supported solver and record installed versions. DIAMOND is installed through Bioconda, not the similarly named Python package.

Example protein FASTA input:

```bash
carve -r "species_faa/*.faa" -o "data/models" --solver gurobi
```

Use `--dna` only for nucleotide input. Validate model quality and taxonomic/model-ID mapping. Model-file presence does not establish species presence in an ecological sample.

## Preserve a reproducible run

Keep the abundance input, sample inclusion rules, units/cutoff, ID mapping, parameters, seed, summary and statistics. For simulations, retain exact models, SMETANA/solver versions, media settings and logs. After installation succeeds, export resolved environments:

```bash
conda env export -n smetana_env > environment_smetana_resolved.yml
conda env export -n carveme_env > environment_carveme_resolved.yml
```

Inspect exports before sharing for local paths and private package sources. Model coverage, thresholds, search conventions, pool size, media and software versions can change findings. Unit tests and artificial data do not establish biological or solver validation.
