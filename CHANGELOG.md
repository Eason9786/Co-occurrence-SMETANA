# Publication notes

## 2.0 - Statistical higher-order co-occurrence

This release replaces the primary random sampler with cross-sample inference inspired by the Methods of Machado et al. (2021), *Polarization of microbial communities between competitive and cooperative metabolism*, doi:10.1038/s41559-020-01353-4.

- **Breaking input change:** `Computing co-occurrence_Single.py` now reads one multi-sample species-abundance table. It cannot infer co-occurrence from a single model folder. `--output` now names a new/empty results directory; `--top-n` is a per-size export limit. The old `--total-budget` and `--excel` options belong to the archived sampler.
- Supports counts, fractional relative abundance, species-level TPM and presence/absence input. TPM is divided by 1,000,000 without renormalizing partial tables. Counts use full input column totals before model restriction.
- Starts with all pairs, counts joint presence across all samples, computes the product-of-prevalences binomial null, applies BH adjustment to the complete evaluated family per order, then applies minimum support, enrichment and q-value thresholds.
- Exports up to 1,000 qualifying communities per size and propagates up to 10,000 qualifying groups by weighted sampling without replacement before extension. Limits are configurable; the population limit is separate from the exported top results.
- Records complete tested statistics, selected headerless community files, settings, input/script hashes, software versions, species/sample identifiers, search trace and stopping reason. Empty results are valid, with no random replacement.
- Defaults to inclusive enrichment tails and frequency-proportional propagation, following the stated statistical question and paper description. Optional strict tails and square-root frequency weights expose two differences in the inspected HiOrCo source; they do not make this an exact historical reproduction. See `docs/METHODS.md` for all differences and limitations.
- Preserves the former portable random sampler as `legacy/random_subcommunities.py`. It is a control-generation utility, not statistical inference or an automatically matched paper-control workflow.
- Adds a fabricated 48-sample TPM example and numerical/input/integration tests; actual user data and biological solver simulations were not used.
- Uses a global, consistent species model library with the existing SMETANA runner rather than treating each single sample as an independent discovery dataset. No runner command-line interface change is required.

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
