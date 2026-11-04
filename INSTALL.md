# Installation / 安装说明

## Recommended starting environment

Use Conda on Linux/WSL for the optional DIAMOND/CarveMe workflow. The two main Python wrappers can also run on Windows with a working SMETANA command and solver. Commands below assume Conda is installed and initialized in your terminal.

```bash
conda env create -f environment_smetana.yml
conda activate smetana_env
python --version
smetana --help
```

This proposes Python 3.10, SMETANA 1.2.1 and the Gurobi Python interface. It is a starting specification, not a tested lockfile. SMETANA 1.2.1 advertises Python 3.8–3.12 support in its [release metadata](https://pypi.org/project/smetana/1.2.1/). The wrappers need Python >=3.8; do not run the sampler under the historical Python 3.6 environment because it uses `math.comb`.

For community generation alone, only the Python standard library is needed. Optional Excel export needs:

```bash
python -m pip install -r requirements.txt
```

## Solver

SMETANA requires a supported optimization solver. Install and configure a suitable Gurobi license following the [official instructions](https://docs.gurobi.com/projects/optimizer/en/current/reference/misc/licensing.html). Installing `gurobipy` does not by itself grant an unrestricted license. Keep license files and activation keys outside this repository.

Use `--solver gurobi` in the batch runner. If you already use CPLEX, install its Python interface in your SMETANA environment and use `--solver cplex` instead. Follow the solver vendor's compatibility requirements; solver installations are not managed by these wrappers.

For Gurobi, a basic environment check is:

```bash
python -c "import gurobipy as gp; m = gp.Model(); m.optimize()"
```

This tests solver initialization only; a successful trivial model does not prove your license supports a full community model. Run a small community of your own valid models before a large batch.

## Optional CarveMe reconstruction

```bash
conda env create -f environment_carveme.yml
conda activate carveme_env
carve --help
diamond version
```

The CarveMe environment is separate to keep reconstruction dependencies independent. CarveMe's version is not locked here; consult its [installation documentation](https://carveme.readthedocs.io/en/latest/installation.html) and record the installed version for reproducibility. DIAMOND is installed through Bioconda, not the unrelated similarly named Python package. Configure a supported solver in this environment as well.

Example protein input:

```bash
carve -r "species_faa/*.faa" -o "data/sample_01" --solver gurobi
```

Use `--dna` only for nucleotide input. Model reconstruction and solving have not been executed as part of this repository's unit tests.

## Verify the wrappers

```bash
python "Computing co-occurrence_Single.py" --help
python Smetana_FWQ.py --help
python -m unittest discover -s tests -v
```

Generate community definitions for your sample, then run `Smetana_FWQ.py` with `--dry-run` before simulation; see [README.md](README.md). All paths can be supplied on the command line; edit no machine-specific script paths.

If SMETANA is installed in another environment, use `--conda-env ENV_NAME`, or activate that environment first. On Windows, run from a terminal where Conda is initialized. On WSL, use Linux paths (for example `/mnt/d/models`) and environments installed inside WSL.

## Record an actual run

Keep a copy of the models, generated communities, sampling seed/options, software versions, solver version, media choices and logs used for the analysis. Export your resolved environments after successful installation:

```bash
conda env export -n smetana_env > environment_smetana_resolved.yml
conda env export -n carveme_env > environment_carveme_resolved.yml
```

Inspect exports before sharing for local paths or private package sources. Changing SMETANA, CarveMe, solver or media configuration can change results. The original notes used SMETANA 1.2.0 / Python 3.6 and CarveMe 1.6.1 / Python 3.12; those historical environments are not reproduced by the supplied environment files.
