#!/usr/bin/env python3
"""Infer higher-order co-occurrence from a species-by-sample abundance table.

Independent implementation of the iterative strategy described by Machado et al.
(2021), doi:10.1038/s41559-020-01353-4. Requires Python >=3.10, NumPy and SciPy.
Default binomial tails are inclusive and propagation weights are frequencies;
see docs/METHODS.md for differences from the public HiOrCo implementation.
"""

import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import random
import sys

import numpy as np
import scipy
from scipy.stats import binom


METHOD_VERSION = "2.0"
PAPER_URL = "https://doi.org/10.1038/s41559-020-01353-4"
REFERENCE_CODE = ("https://github.com/cdanielmachado/HiOrCo/blob/"
                  "ee9af1f1911f6313376247352051f25e09b2fd0a/hiorco/method.py")


def _validate_ids(values, label):
    if any(not isinstance(value, str) or not value.strip()
           or any(char in value for char in "\t\r\n") for value in values):
        raise ValueError("{} IDs must be nonempty strings without tabs or newlines".format(label))
    if len(set(values)) != len(values):
        raise ValueError("Duplicate {} IDs are not allowed".format(label.lower()))


def _positive_integer(value, name, minimum=1):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError("{} must be an integer >= {}".format(name, minimum))


def _finite_number(value, name, minimum=None, maximum=None):
    try:
        valid = (not isinstance(value, (bool, np.bool_)) and math.isfinite(value)
                 and (minimum is None or value >= minimum)
                 and (maximum is None or value <= maximum))
    except TypeError:
        valid = False
    if not valid:
        raise ValueError("{} must be finite and between {} and {}".format(name, minimum, maximum))


def load_presence_table(path, input_kind="counts", abundance_cutoff=0.001,
                        sep="tab", models_dir=None):
    """Load species rows/sample columns and apply a strict presence threshold.

    Counts normalize FULL input columns before model restriction. Relative data
    are fractions in [0,1], while TPM values are divided by 1,000,000. Neither
    relative data nor TPM are renormalized after input/model subsetting.
    All sample columns remain in the denominator, including empty samples.
    """
    if input_kind not in ("counts", "relative", "tpm", "presence-absence"):
        raise ValueError("input_kind must be counts, relative, tpm, or presence-absence")
    _finite_number(abundance_cutoff, "abundance_cutoff", 0, 1)
    if sep not in ("tab", "comma"):
        raise ValueError("sep must be tab or comma")
    source = Path(path)
    if not source.is_file():
        raise ValueError("Abundance table does not exist: {}".format(source))
    species_ids, values = [], []
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t" if sep == "tab" else ",")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("Abundance table is empty") from None
        if len(header) < 3:
            raise ValueError("Table requires a header with at least two sample columns")
        sample_ids = header[1:]
        _validate_ids(sample_ids, "Sample")
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError("Line {}: expected {} columns, found {}".format(line_number, len(header), len(row)))
            species_ids.append(row[0])
            if any(not value.strip() for value in row[1:]):
                raise ValueError("Line {}: missing abundance value".format(line_number))
            try:
                values.append([float(value) for value in row[1:]])
            except ValueError:
                raise ValueError("Line {}: all abundance values must be numeric".format(line_number)) from None
    _validate_ids(species_ids, "Species")
    if len(species_ids) < 2:
        raise ValueError("At least two species rows are required")
    data = np.asarray(values, dtype=float)
    if not np.isfinite(data).all() or np.any(data < 0):
        raise ValueError("Abundance values must be finite and nonnegative")
    with np.errstate(over="ignore", invalid="ignore"):
        column_totals = data.sum(axis=0)
    if not np.isfinite(column_totals).all():
        raise ValueError("Abundance column totals overflow; rescale the input values")
    if input_kind == "counts":
        relative = np.divide(data, column_totals, out=np.zeros_like(data), where=column_totals > 0)
        presence = relative > abundance_cutoff
    elif input_kind == "relative":
        if np.any(data > 1) or np.any(column_totals > 1 + 1e-8):
            raise ValueError("Relative abundances must be fractions in [0,1] with column sums <= 1; convert percentages first")
        presence = data > abundance_cutoff
    elif input_kind == "tpm":
        if np.any(data > 1e6) or np.any(column_totals > 1e6 + 1):
            raise ValueError("TPM values must be <= 1,000,000 with column sums <= 1,000,001; supply species-level TPM, not unaggregated gene TPM")
        presence = (data / 1e6) > abundance_cutoff
    else:
        if not np.isin(data, [0, 1]).all():
            raise ValueError("presence-absence input must contain only 0 and 1")
        presence = data.astype(bool)
    normalization = {"counts": "full input columns before model restriction",
                     "relative": "none (already fractions)",
                     "tpm": "divide by 1,000,000; no column renormalization",
                     "presence-absence": "none"}[input_kind]
    metadata = {
        "path": str(source.resolve()), "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "input_kind": input_kind, "delimiter": sep,
        "abundance_cutoff": abundance_cutoff if input_kind != "presence-absence" else None,
        "threshold_raw_units": abundance_cutoff * 1e6 if input_kind == "tpm" else None,
        "presence_rule": "value == 1" if input_kind == "presence-absence" else "relative_abundance > cutoff",
        "normalization": normalization, "input_species": len(species_ids), "input_samples": len(sample_ids),
        "sample_ids": sample_ids,
        "zero_total_sample_ids": [sample_ids[i] for i in np.flatnonzero(column_totals == 0)],
        "model_restriction": None, "unmatched_species_ids": [],
    }
    if models_dir is not None:
        model_folder = Path(models_dir)
        if not model_folder.is_dir():
            raise ValueError("Model directory does not exist: {}".format(model_folder))
        model_ids = {item.stem for item in model_folder.iterdir()
                     if item.is_file() and item.name.endswith(".xml")}
        matched = [i for i, name in enumerate(species_ids) if name in model_ids]
        metadata["model_restriction"] = str(model_folder.resolve())
        metadata["unmatched_species_ids"] = [name for name in species_ids if name not in model_ids]
        if len(matched) < 2:
            raise ValueError("Fewer than two table species match .xml filenames in the model directory")
        species_ids = [species_ids[i] for i in matched]
        presence = presence[matched, :]
    metadata["analyzed_species"] = len(species_ids)
    metadata["retained_samples"] = len(sample_ids)
    metadata["samples_with_fewer_than_two_present_species"] = int(np.count_nonzero(presence.sum(axis=0) < 2))
    return {"presence": presence, "species_ids": species_ids, "sample_ids": sample_ids, "metadata": metadata}


def adjust_bh(p_values):
    """Benjamini-Hochberg adjusted p-values for the complete supplied family."""
    p_values = np.asarray(p_values, dtype=float)
    if p_values.ndim != 1 or not np.isfinite(p_values).all() or np.any((p_values < 0) | (p_values > 1)):
        raise ValueError("p_values must be a one-dimensional finite array in [0, 1]")
    count = len(p_values)
    if not count:
        return np.array([], dtype=float)
    order = np.argsort(p_values, kind="stable")
    ranked = p_values[order] * count / np.arange(1, count + 1)
    adjusted = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    result = np.empty(count, dtype=float)
    result[order] = adjusted
    return result


def _weighted_population(scores, population_size, rng, selection_weight):
    """Weighted sampling without replacement through exponential-race keys."""
    if len(scores) <= population_size:
        return [score["species"] for score in scores]
    priorities = []
    for score in sorted(scores, key=lambda item: item["species"]):
        weight = score["observed"] if selection_weight == "frequency" else math.sqrt(score["observed"])
        priority = -math.log1p(-rng.random()) / weight
        priorities.append((priority, score["species"]))
    priorities.sort()
    return sorted(species for _, species in priorities[:population_size])


def _expand_population(population, species_ids, max_candidates):
    candidates = set()
    for community in sorted(population):
        present = set(community)
        for species in species_ids:
            if species not in present:
                candidates.add(tuple(sorted(community + (species,))))
                if len(candidates) > max_candidates:
                    raise ValueError("Candidate count exceeds --max-candidates {}; increase the cap with adequate memory, explicitly restrict species, or reduce --population-size".format(max_candidates))
    return sorted(candidates)


def infer_cooccurrence(presence, species_ids, max_combo_size=5, top_n=1000,
                       population_size=10000, min_samples=10, min_fold=2,
                       fdr_cutoff=0.05, seed=42, selection_weight="frequency",
                       pvalue_tail="inclusive", max_candidates=1000000):
    """Infer enriched communities, preserving all species and sample columns.

    BH adjustment is per order over ALL generated candidates, before support/fold
    filtering. Higher orders are an adaptive heuristic, not exhaustive testing or
    globally FDR-controlled discoveries. ``significant`` means all three criteria.
    """
    for name, value, minimum in (("max_combo_size", max_combo_size, 2), ("top_n", top_n, 1),
                                ("population_size", population_size, 1), ("min_samples", min_samples, 1),
                                ("max_candidates", max_candidates, 1)):
        _positive_integer(value, name, minimum)
    _finite_number(min_fold, "min_fold", 1)
    _finite_number(fdr_cutoff, "fdr_cutoff", 0, 1)
    if fdr_cutoff == 0:
        raise ValueError("fdr_cutoff must be greater than zero")
    if selection_weight not in ("frequency", "sqrt-frequency"):
        raise ValueError("selection_weight must be frequency or sqrt-frequency")
    if pvalue_tail not in ("inclusive", "strict"):
        raise ValueError("pvalue_tail must be inclusive or strict")
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)):
        raise ValueError("seed must be an integer")
    species_ids = list(species_ids)
    _validate_ids(species_ids, "Species")
    data = np.asarray(presence)
    if (data.ndim != 2 or data.shape[0] != len(species_ids)
            or len(species_ids) < 2 or data.shape[1] < 2):
        raise ValueError("presence must have one row per species, >=2 species and >=2 samples")
    if not np.isin(data, [False, True]).all():
        raise ValueError("presence must contain only booleans or numeric 0/1 values")
    data = data.astype(bool)
    n_samples = data.shape[1]
    row_order = sorted(range(len(species_ids)), key=lambda index: species_ids[index])
    data = data[row_order, :]
    species_ids = [species_ids[i] for i in row_order]
    row_index = {species: i for i, species in enumerate(species_ids)}
    frequencies = data.sum(axis=1) / n_samples
    bitsets = [int.from_bytes(np.packbits(row, bitorder="little").tobytes(), "little") for row in data]
    rng = random.Random(int(seed))
    settings = {
        "max_combo_size": int(max_combo_size), "top_n": int(top_n), "population_size": int(population_size),
        "min_samples": int(min_samples), "min_fold": float(min_fold), "fdr_cutoff": float(fdr_cutoff),
        "seed": int(seed), "selection_weight": selection_weight, "pvalue_tail": pvalue_tail,
        "max_candidates": int(max_candidates),
    }
    pair_count = math.comb(len(species_ids), 2)
    if pair_count > max_candidates:
        raise ValueError("{} pairwise candidates exceed --max-candidates {}; increase the cap with adequate memory or explicitly restrict the species universe".format(pair_count, max_candidates))
    candidates = list(itertools.combinations(species_ids, 2))
    orders, trace = [], []
    stop_reason = "max_size_reached"
    maximum_size = min(max_combo_size, len(species_ids))
    for size in range(2, maximum_size + 1):
        observed = np.empty(len(candidates), dtype=np.int64)
        probabilities = np.empty(len(candidates), dtype=float)
        for index, community in enumerate(candidates):
            indices = [row_index[species] for species in community]
            intersection = bitsets[indices[0]]
            for member in indices[1:]:
                intersection &= bitsets[member]
            observed[index] = intersection.bit_count()
            probabilities[index] = math.prod(frequencies[member] for member in indices)
        expected = n_samples * probabilities
        tails = observed - 1 if pvalue_tail == "inclusive" else observed
        p_values = np.asarray(binom.sf(tails, n_samples, probabilities), dtype=float)
        q_values = adjust_bh(p_values)
        scores = []
        for index, community in enumerate(candidates):
            total, expectation = int(observed[index]), float(expected[index])
            passes_fdr = bool(q_values[index] <= fdr_cutoff)
            fold_threshold = min_fold * expectation
            # Preserve inclusive effect-size boundaries despite roundoff, e.g.
            # 14 observations versus a mathematically exact threshold of 14.
            passes_fold = (total >= fold_threshold or
                           math.isclose(total, fold_threshold, rel_tol=1e-12, abs_tol=0.0))
            significant = bool(passes_fdr and total >= min_samples and passes_fold)
            scores.append({"species": community, "observed": total, "expected": expectation,
                           "fold": total / expectation if expectation > 0 else (math.inf if total > 0 else 0.0),
                           "p_value": float(p_values[index]), "q_value": float(q_values[index]),
                           "passes_fdr": passes_fdr, "significant": significant,
                           "selected": False, "community_id": ""})
        qualifying = sorted((score for score in scores if score["significant"]),
                            key=lambda score: (-score["observed"], score["species"]))
        selected = qualifying[:top_n]
        for rank, score in enumerate(selected, start=1):
            score["selected"] = True
            score["community_id"] = "community_{}_{:04d}".format(size, rank)
        propagated = (_weighted_population(qualifying, population_size, rng, selection_weight)
                      if qualifying and size < maximum_size else [])
        orders.append({"size": size, "scores": scores, "selected": selected, "propagated": propagated})
        trace.append({"size": size, "evaluated_candidates": len(candidates),
                      "passes_fdr": int(np.count_nonzero(q_values <= fdr_cutoff)),
                      "significant": len(qualifying), "selected": len(selected), "propagated": len(propagated)})
        if not qualifying:
            stop_reason = "no_significant_communities_at_size_{}".format(size)
            break
        if size < maximum_size:
            candidates = _expand_population(propagated, species_ids, max_candidates)
            if not candidates:
                stop_reason = "no_extensions_at_size_{}".format(size)
                break
    else:
        if maximum_size < max_combo_size:
            stop_reason = "all_species_used"
    return {"orders": orders, "trace": trace, "stop_reason": stop_reason, "settings": settings,
            "n_species": len(species_ids), "n_samples": int(n_samples), "species_ids": species_ids}


def _write_communities(path, selected):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for score in selected:
            for species in score["species"]:
                writer.writerow([score["community_id"], species])


def write_results(result, output, input_metadata=None):
    """Write complete statistics, selected headerless TSVs, and run provenance.

    Reject existing nonempty directories to prevent mixing new and stale output.
    No XML models or biological input files are copied into the output.
    """
    output = Path(output)
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("Output must be a new or empty directory: {}".format(output))
    output.mkdir(parents=True, exist_ok=True)
    all_selected = []
    columns = ["community_id", "size", "species", "observed", "expected", "fold",
               "p_value", "q_value", "passes_fdr", "significant", "selected"]
    for order in result["orders"]:
        size = order["size"]
        with (output / "size_{}_statistics.tsv".format(size)).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for score in order["scores"]:
                row = dict(score, size=size)
                row["species"] = json.dumps(score["species"], ensure_ascii=False)
                writer.writerow(row)
        _write_communities(output / "size_{}_communities.tsv".format(size), order["selected"])
        all_selected.extend(order["selected"])
    _write_communities(output / "communities.tsv", all_selected)
    summary = {
        "method": "iterative higher-order co-occurrence inference", "method_version": METHOD_VERSION,
        "software": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "paper": PAPER_URL, "reference_implementation_inspected": REFERENCE_CODE,
        "implementation": "independent implementation, not exact reproduction of HiOrCo",
        "input": input_metadata or {}, "settings": result["settings"],
        "n_species": result["n_species"], "n_samples": result["n_samples"],
        "analyzed_species_ids": result["species_ids"],
        "n_selected_communities": len(all_selected), "trace": result["trace"], "stop_reason": result["stop_reason"],
        "definitions": {
            "expected": "n_samples * product(individual_species_prevalences)",
            "binomial_tail": "P(X >= observed)" if result["settings"]["pvalue_tail"] == "inclusive" else "P(X > observed)",
            "multiple_testing": "Benjamini-Hochberg per order over all evaluated candidates before support/fold filtering",
            "tested_universe": "all input species after optional model restriction, including zero-prevalence species",
            "sample_denominator": "all input sample columns, including zero-total and low-richness samples",
            "significant": "passes FDR, minimum observed support, and minimum observed/expected fold",
            "ranking": "observed descending; lexical species tuple breaks ties",
            "propagation": "weighted without replacement among all significant sets, independent of top-n export",
        },
        "limitations": [
            "Adaptive higher-order search can miss communities without a propagated significant subset.",
            "Per-order BH values describe the evaluated family; adaptive search and dependent hypotheses do not guarantee global FDR control.",
            "The independent-species binomial null does not adjust for habitat, batch, repeated measures, or compositional confounding.",
            "Significant co-occurrence is an association, not evidence of direct interaction or cooperation.",
        ],
    }
    with (output / "run_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("abundance_table", help="Headered table: species rows and sample columns; first column contains species IDs")
    parser.add_argument("--input-kind", choices=("counts", "relative", "tpm", "presence-absence"), default="counts")
    parser.add_argument("--sep", choices=("tab", "comma"), default="tab")
    parser.add_argument("--abundance-cutoff", type=float, default=0.001,
                        help="Strict relative-abundance presence threshold (default 0.001 = 0.1%%); ignored for presence-absence input")
    parser.add_argument("--models-dir", help="Optional flat XML library; restrict exact filename-stem species IDs after normalization")
    parser.add_argument("--max-combo-size", type=int, default=5)
    parser.add_argument("--top-n", type=int, default=1000, help="Maximum exported communities PER SIZE, not across all sizes")
    parser.add_argument("--population-size", type=int, default=10000, help="Maximum significant sets propagated to the next order")
    parser.add_argument("--min-samples", type=int, default=10)
    parser.add_argument("--min-fold", type=float, default=2)
    parser.add_argument("--fdr-cutoff", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--selection-weight", choices=("frequency", "sqrt-frequency"), default="frequency")
    parser.add_argument("--pvalue-tail", choices=("inclusive", "strict"), default="inclusive")
    parser.add_argument("--max-candidates", type=int, default=1000000,
                        help="Hard per-order cap; exceedance fails rather than silently truncating tests")
    parser.add_argument("--output", default="cooccurrence_results", help="New or empty output DIRECTORY (default: cooccurrence_results)")
    args = parser.parse_args(argv)
    output = Path(args.output)
    try:
        if output.exists() and (not output.is_dir() or any(output.iterdir())):
            raise ValueError("Output must be a new or empty directory: {}".format(output))
        loaded = load_presence_table(args.abundance_table, args.input_kind, args.abundance_cutoff, args.sep, args.models_dir)
        result = infer_cooccurrence(
            loaded["presence"], loaded["species_ids"], args.max_combo_size, args.top_n,
            args.population_size, args.min_samples, args.min_fold, args.fdr_cutoff, args.seed,
            args.selection_weight, args.pvalue_tail, args.max_candidates)
        summary = write_results(result, output, loaded["metadata"])
    except (OSError, ValueError, UnicodeError, csv.Error) as error:
        parser.exit(1, "Error: {}\n".format(error))
    for row in result["trace"]:
        print("Size {size}: evaluated {evaluated_candidates}; significant {significant}; exported {selected}; propagated {propagated}".format(**row))
    print("Saved {} selected communities to {}".format(summary["n_selected_communities"], output / "communities.tsv"))
    print("Stop reason: {}".format(result["stop_reason"]))
    if not summary["n_selected_communities"]:
        print("No communities passed all criteria. Do not run SMETANA on the empty communities file.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
