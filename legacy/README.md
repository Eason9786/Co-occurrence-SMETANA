# Optional dataset-specific preparation scripts

## Random community generator

`random_subcommunities.py` reads a folder of XML model filenames and samples synthetic community combinations. It can supply random candidates for a defined control design; it does not calculate statistical co-occurrence. Run from the repository root:

```bash
python legacy/random_subcommunities.py data/models --max-combo-size 5 --total-budget 1000 --seed 42 --output results/random_candidates.tsv
```

The budget is balanced across feasible sizes. Size-matched controls require matching the inferred groups' counts separately for each size; this utility does not do that automatically. Optional `--excel` output requires openpyxl.

## Preprocessing helpers

These helpers use the Python standard library. Set their input/output path constants before running; paths are relative to the working directory, not to each script. Check the dataset-specific identifier rules before use.

- `1_tsv_to_faa/p02_split_faa_file_from_CK.py` reads a tab-separated annotation file with `Species` and `#Query` columns and a protein FASTA file, then writes one `.faa` file per species. Rows with an empty species or `-` are excluded. Query IDs lose their final underscore-delimited component; FASTA IDs lose one or two trailing numeric components. These dataset-specific rules must be checked against your actual headers. Repeated normalized FASTA IDs overwrite earlier entries, and missing IDs are reported and skipped. Species output names replace spaces and `/` with `_`.
- `2_faa_TOP1/filter_files_single.py` selects one `.faa` file from each filename group in one input folder. The group is the text after `s__` and before the next `_` (a genus-like label). It copies the largest file in **bytes** per group to the output folder; it does not rank species by abundance. Files without the expected naming pattern are skipped.
- `2_faa_TOP1/filter_files_all.py` applies the same filename grouping and largest-file selection separately to each immediate subfolder, creating an output folder with the `_TOP1` suffix for each sample.

Provide your own annotation and sequence data. Existing files with matching output names can be overwritten. Neither relative-abundance filtering, model reconstruction, statistical co-occurrence tests, nor SMETANA execution is performed by these helpers.
