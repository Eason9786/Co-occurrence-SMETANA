"""Numerical and integration tests for statistical higher-order co-occurrence."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
from scipy.stats import binom

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "Computing co-occurrence_Single.py"
SPEC = importlib.util.spec_from_file_location("statistical_cooccurrence", SCRIPT)
cooccurrence = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = cooccurrence
SPEC.loader.exec_module(cooccurrence)


class CooccurrenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT)
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)

    def table(self, rows, samples=None):
        path = self.folder / "abundances.tsv"
        if samples is None and len(rows[0]) == 2:
            rows = [[row[0], row[1], row[1]] for row in rows]
        if samples is None:
            samples = ["S{}".format(i) for i in range(len(rows[0]) - 1)]
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["species"] + samples)
            writer.writerows(rows)
        return path

    def call_cli(self, *args):
        return subprocess.run([sys.executable, "-X", "utf8", str(SCRIPT)] + list(map(str, args)),
                              capture_output=True, text=True, encoding="utf-8")

    def test_bh_known_values(self):
        np.testing.assert_allclose(cooccurrence.adjust_bh([0.01, 0.04, 0.03, 0.002]),
                                   [0.02, 0.04, 0.04, 0.008], atol=1e-15)
        np.testing.assert_allclose(cooccurrence.adjust_bh([0.01, 0.04, 0.2]), [0.03, 0.06, 0.2])

    def test_synthetic_48_tpm_recovers_known_pairs_and_triple(self):
        loaded = cooccurrence.load_presence_table(ROOT / "examples/synthetic_tpm_48.tsv", input_kind="tpm")
        self.assertEqual(loaded["presence"].shape, (4, 48))
        result = cooccurrence.infer_cooccurrence(loaded["presence"], loaded["species_ids"])
        selected = [score for order in result["orders"] for score in order["selected"]]
        self.assertEqual({tuple(score["species"]) for score in selected},
                         {("species_A", "species_B"), ("species_A", "species_C"),
                          ("species_B", "species_C"), ("species_A", "species_B", "species_C")})
        pair = next(score for score in selected if tuple(score["species"]) == ("species_A", "species_B"))
        self.assertEqual(pair["observed"], 12)
        self.assertAlmostEqual(pair["expected"], 3)
        self.assertAlmostEqual(pair["fold"], 4)
        triple = next(score for score in selected if len(score["species"]) == 3)
        self.assertAlmostEqual(triple["expected"], 0.75)
        self.assertAlmostEqual(triple["fold"], 16)

    def test_inclusive_binomial_and_bh_over_all_pairs(self):
        presence = np.zeros((3, 100), dtype=bool)
        presence[:2, :10] = True
        result = cooccurrence.infer_cooccurrence(presence, ["A", "B", "C"], max_combo_size=2)
        scores = result["orders"][0]["scores"]
        self.assertEqual(len(scores), 3)
        pair = next(score for score in scores if tuple(score["species"]) == ("A", "B"))
        self.assertAlmostEqual(pair["expected"], 1)
        p = binom.sf(9, 100, 0.01)
        self.assertAlmostEqual(pair["p_value"], p, delta=p * 1e-12)
        self.assertAlmostEqual(pair["q_value"], 3 * p, delta=p * 1e-11)
        for score in scores:
            if "C" in score["species"]:
                self.assertEqual(score["p_value"], 1)
                self.assertEqual(score["q_value"], 1)

    def test_strict_tail_is_explicit_reference_option(self):
        presence = np.zeros((2, 100), dtype=bool)
        presence[:, :10] = True
        result = cooccurrence.infer_cooccurrence(presence, ["A", "B"], pvalue_tail="strict")
        self.assertAlmostEqual(result["orders"][0]["scores"][0]["p_value"],
                               binom.sf(10, 100, 0.01), delta=1e-18)

    def test_full_counts_denominator_precedes_model_restriction(self):
        path = self.table([["A", 1, 1], ["B", 1, 1], ["unmodeled", 98, 98]])
        models = self.folder / "models"
        models.mkdir()
        for name in ("A.xml", "B.xml"):
            (models / name).touch()
        loaded = cooccurrence.load_presence_table(path, input_kind="counts", abundance_cutoff=0.02, models_dir=models)
        self.assertEqual(loaded["species_ids"], ["A", "B"])
        self.assertFalse(loaded["presence"].any())

    def test_partial_tpm_is_not_renormalized_and_cutoff_is_strict(self):
        loaded = cooccurrence.load_presence_table(self.table([["A", 1000, 1001], ["B", 1, 999]]), input_kind="tpm")
        np.testing.assert_array_equal(loaded["presence"], [[False, True], [False, False]])

    def test_partial_relative_and_zero_samples_are_retained(self):
        loaded = cooccurrence.load_presence_table(self.table([["A", 0.001, 0], ["B", 0.002, 0]]), input_kind="relative")
        np.testing.assert_array_equal(loaded["presence"], [[False, False], [True, False]])
        self.assertEqual(len(loaded["sample_ids"]), 2)

    def test_presence_absence_preserves_literal_na_id(self):
        loaded = cooccurrence.load_presence_table(self.table([["NA", 1, 0], ["B", 0, 1]]),
                                                  input_kind="presence-absence", abundance_cutoff=0.9)
        self.assertIn("NA", loaded["species_ids"])
        self.assertEqual(int(loaded["presence"].sum()), 2)

    def test_invalid_values_and_duplicate_identifiers_fail(self):
        cases = [([["A", ""], ["B", 1]], None), ([["A", "nan"], ["B", 1]], None),
                 ([["A", "inf"], ["B", 1]], None), ([["A", -1], ["B", 1]], None),
                 ([["A", "text"], ["B", 1]], None), ([["A", 1], ["A", 2]], None),
                 ([["", 1], ["B", 2]], None), ([["A", 1, 2], ["B", 3, 4]], ["S1", "S1"])]
        for rows, samples in cases:
            with self.subTest(rows=rows, samples=samples), self.assertRaises(ValueError):
                cooccurrence.load_presence_table(self.table(rows, samples))

    def test_unit_validation(self):
        cases = [("relative", [["A", 50], ["B", 50]]), ("relative", [["A", 0.8], ["B", 0.8]]),
                 ("tpm", [["A", 800000], ["B", 800000]]), ("presence-absence", [["A", 2], ["B", 0]])]
        for kind, rows in cases:
            with self.subTest(kind=kind), self.assertRaises(ValueError):
                cooccurrence.load_presence_table(self.table(rows), input_kind=kind)

    def test_extension_does_not_require_a_pairwise_clique(self):
        presence = np.zeros((3, 200), dtype=bool)
        presence[:2, :40] = True
        presence[2, :20] = True
        presence[2, 40:120] = True
        result = cooccurrence.infer_cooccurrence(presence, ["A", "B", "C"])
        self.assertEqual({tuple(score["species"]) for score in result["orders"][0]["selected"]}, {("A", "B")})
        triples = [score for order in result["orders"] if order["size"] == 3 for score in order["selected"]]
        self.assertEqual(len(triples), 1)
        self.assertEqual(triples[0]["observed"], 20)
        self.assertAlmostEqual(triples[0]["expected"], 4)

    def test_heuristic_misses_pure_three_way_association(self):
        patterns = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0], [1, 1, 1]], dtype=bool)
        result = cooccurrence.infer_cooccurrence(np.repeat(patterns, 25, axis=0).T, ["A", "B", "C"])
        self.assertEqual([order["size"] for order in result["orders"]], [2])
        self.assertFalse(result["orders"][0]["selected"])

    def test_no_hits_has_no_random_fallback(self):
        result = cooccurrence.infer_cooccurrence(np.ones((4, 48), dtype=bool), list("ABCD"))
        self.assertFalse(any(order["selected"] for order in result["orders"]))
        self.assertTrue(result["stop_reason"])

    def test_exact_twofold_boundary_is_inclusive_despite_roundoff(self):
        presence = np.zeros((2, 50), dtype=bool)
        presence[0, :14] = True
        presence[1, :25] = True
        result = cooccurrence.infer_cooccurrence(presence, ["A", "B"])
        self.assertEqual(len(result["orders"][0]["selected"]), 1)

    def test_count_sum_overflow_and_too_few_samples(self):
        with self.assertRaisesRegex(ValueError, "overflow"):
            cooccurrence.load_presence_table(self.table([["A", 1e308, 1], ["B", 1e308, 1]]))
        path = self.table([["A", 1], ["B", 1]], samples=["S1"])
        with self.assertRaisesRegex(ValueError, "two sample"):
            cooccurrence.load_presence_table(path)

    def test_sample_shortfall_produces_no_hits(self):
        presence = np.zeros((2, 9), dtype=bool)
        presence[:, :4] = True
        result = cooccurrence.infer_cooccurrence(presence, ["A", "B"])
        self.assertFalse(result["orders"][0]["selected"])

    def test_top_n_is_per_order_and_distinct_from_population(self):
        presence = np.zeros((4, 100), dtype=bool)
        presence[:, :10] = True
        result = cooccurrence.infer_cooccurrence(presence, list("ABCD"), max_combo_size=3, top_n=1, population_size=6)
        self.assertEqual([len(order["selected"]) for order in result["orders"]], [1, 1])
        self.assertEqual(len(result["orders"][0]["propagated"]), 6)

    def test_population_reproducibility_and_candidate_deduplication(self):
        presence = np.zeros((5, 100), dtype=bool)
        presence[:, :10] = True
        kwargs = dict(max_combo_size=4, top_n=2, population_size=2, seed=99)
        first = cooccurrence.infer_cooccurrence(presence, list("ABCDE"), **kwargs)
        self.assertEqual(first, cooccurrence.infer_cooccurrence(presence, list("ABCDE"), **kwargs))
        for order in first["orders"]:
            candidates = [tuple(score["species"]) for score in order["scores"]]
            self.assertEqual(len(candidates), len(set(candidates)))

    def test_candidate_guard_fails_instead_of_truncating(self):
        with self.assertRaises((ValueError, RuntimeError)):
            cooccurrence.infer_cooccurrence(np.ones((5, 48), dtype=bool), list("ABCDE"), max_candidates=9)

    def test_cli_output_and_runner_contract(self):
        output = self.folder / "run output"
        result = self.call_cli(ROOT / "examples/synthetic_tpm_48.tsv", "--input-kind", "tpm", "--output", output)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIsInstance(json.loads((output / "run_summary.json").read_text(encoding="utf-8")), dict)
        with (output / "communities.tsv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 9)
        self.assertTrue(all(len(row) == 2 and row[1].startswith("species_") for row in rows))
        self.assertEqual(len({row[0] for row in rows}), 4)
        for size in (2, 3):
            self.assertTrue((output / ("size_{}_statistics.tsv".format(size))).is_file())
        models = self.folder / "model library"
        models.mkdir()
        for species in "ABCD":
            (models / ("species_" + species + ".xml")).write_text("<placeholder/>", encoding="utf-8")
        dryrun = subprocess.run([sys.executable, "-X", "utf8", str(ROOT / "Smetana_FWQ.py"), str(models),
                                 "--single", "--communities", str(output / "communities.tsv"), "--dry-run"],
                                capture_output=True, text=True, encoding="utf-8")
        self.assertEqual(dryrun.returncode, 0, dryrun.stderr)
        self.assertFalse((models / "smetana_run.log").exists())

    def test_nonempty_output_is_not_overwritten(self):
        output = self.folder / "output"
        output.mkdir()
        marker = output / "keep.txt"
        marker.write_text("keep", encoding="utf-8")
        result = self.call_cli(ROOT / "examples/synthetic_tpm_48.tsv", "--input-kind", "tpm", "--output", output)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(marker.read_text(encoding="utf-8"), "keep")


if __name__ == "__main__":
    unittest.main()
