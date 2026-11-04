# Publication notes

## Initial public release

Based on the local `Computing co-occurrence_Single.py` and `Smetana_FWQ.py` scripts.

- Retained the main script filenames and the original workflow: generate candidate subcommunities from model filenames, then run SMETANA per sample.
- Replaced machine-specific paths with command-line arguments.
- Corrected community-size budget allocation and bounded it by available distinct combinations.
- Made sampling reproducible across processes, without changing global random state.
- Added headerless community TSV export for SMETANA and kept optional human-readable Excel export.
- Replaced shell-based environment activation with an active-environment command or explicit `conda run`; expanded model files in Python and checked subprocess results.
- Added per-sample logs, input validation and dry-run mode.
- Added installation guidance, environment templates, English/Chinese documentation and focused unit tests.
- Kept selected original preprocessing helpers in `legacy/`, with their dataset-specific assumptions documented.
- Corrected the README to describe candidate-combination sampling. Statistical significance testing, abundance filtering, frequency ranking, benchmark results and automatic network/flux-map generation are not implemented or established by these scripts.

The original local files were not edited. Raw sequence datasets, private license information, local installation notes, destructive cleanup helpers and superseded workflow variants are not included in this release.
