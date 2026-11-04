"""Sampling and export checks that do not require metabolic models or a solver."""

import csv
import importlib.util
import itertools
import math
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "legacy" / "random_subcommunities.py"
SPEC = importlib.util.spec_from_file_location("community_sampling", SCRIPT)
sampling = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sampling)


class SamplingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=SCRIPT.parent)
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name) / "model files"
        self.folder.mkdir()

    def models(self, count):
        for index in range(count):
            (self.folder / "species_{}.xml".format(index)).touch()

    def test_small_population_returns_all_without_overallocation(self):
        for count in range(2, 7):
            with self.subTest(count=count):
                self.models(count)
                communities = sampling.generate_balanced_sub_communities(self.folder)
                expected = sum(math.comb(count, size)
                               for size in range(2, min(count, 5) + 1))
                self.assertEqual(len(communities), expected)
                self.assertEqual(len(set(communities)), expected)

    def test_quota_is_balanced_and_seed_is_local(self):
        self.models(15)
        before = random.getstate()
        communities = sampling.generate_balanced_sub_communities(
            self.folder, total_budget=101, seed=17)
        self.assertEqual(random.getstate(), before)
        self.assertEqual(len(communities), 101)
        self.assertEqual(len(set(communities)), 101)
        counts = [sum(len(combo) == size for combo in communities)
                  for size in range(2, 6)]
        self.assertLessEqual(max(counts) - min(counts), 1)
        self.assertEqual(communities, sampling.generate_balanced_sub_communities(
            self.folder, total_budget=101, seed=17))

    def test_one_community_and_capacity_redistribution(self):
        self.models(4)
        self.assertEqual(len(sampling.generate_balanced_sub_communities(
            self.folder, total_budget=1)), 1)
        communities = sampling.generate_balanced_sub_communities(
            self.folder, total_budget=10)
        self.assertEqual(len(communities), 10)
        self.assertEqual(sorted(sum(len(combo) == k for combo in communities)
                                for k in (2, 3, 4)), [1, 4, 5])

    def test_unranking_matches_all_combinations(self):
        species = [str(index) for index in range(8)]
        for size in range(2, 7):
            expected = [frozenset(combo) for combo in itertools.combinations(species, size)]
            self.assertEqual([sampling._unrank_combination(species, size, rank)
                              for rank in range(len(expected))], expected)

    def test_rank_sampling_handles_huge_and_dense_populations(self):
        huge = 10 ** 30
        ranks = sampling._sample_ranks(random.Random(42), huge, 100)
        self.assertEqual(len(set(ranks)), 100)
        self.assertTrue(all(0 <= rank < huge for rank in ranks))
        self.assertEqual(sampling._sample_ranks(random.Random(42), 100, 100),
                         list(range(100)))

    def test_only_regular_xml_files_are_models(self):
        self.models(1)
        (self.folder / "ignored.xml").mkdir()
        (self.folder / "ignored.txt").touch()
        self.assertEqual(sampling.generate_balanced_sub_communities(self.folder), [])

    def test_invalid_parameters_fail(self):
        self.models(2)
        for kwargs in ({"max_combo_size": 1}, {"total_budget": 0}, {"total_budget": -1}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                sampling.generate_balanced_sub_communities(self.folder, **kwargs)
        with self.assertRaises(ValueError):
            sampling.generate_balanced_sub_communities(self.folder / "missing")

    def test_tsv_has_no_header_and_uses_final_xml_suffix(self):
        (self.folder / "species.xml.extra.xml").touch()
        (self.folder / "other.xml").touch()
        communities = sampling.generate_balanced_sub_communities(self.folder)
        output = Path(self.temp.name) / "result.tsv"
        sampling.generate_tsv_output(communities, output)
        with output.open(encoding="utf-8", newline="") as handle:
            self.assertEqual(list(csv.reader(handle, delimiter="\t")),
                             [["community1", "other"], ["community1", "species.xml.extra"]])

    def test_cli_output_is_identical_across_hash_seeds(self):
        self.models(25)
        outputs = []
        for seed in ("1", "54321"):
            output = Path(self.temp.name) / (seed + ".tsv")
            env = dict(os.environ, PYTHONHASHSEED=seed)
            result = subprocess.run(
                [sys.executable, str(SCRIPT), str(self.folder),
                 "--total-budget", "37", "--output", str(output)],
                env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            outputs.append(output.read_bytes())
        self.assertEqual(outputs[0], outputs[1])


if __name__ == "__main__":
    unittest.main()
