"""Runner integration checks use model placeholders and mock the SMETANA process."""

import contextlib
import importlib.util
import io
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


SCRIPT = Path(__file__).resolve().parents[1] / "Smetana_FWQ.py"
SPEC = importlib.util.spec_from_file_location("smetana_runner", SCRIPT)
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=SCRIPT.parent)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def sample(self, name="sample with spaces", content="community1\tmodel A\ncommunity1\tmodel_B\n"):
        folder = self.root / name
        folder.mkdir(parents=True)
        (folder / "model A.xml").write_text("<placeholder/>", encoding="utf-8")
        (folder / "model_B.xml").write_text("<placeholder/>", encoding="utf-8")
        (folder / "communities.tsv").write_text(content, encoding="utf-8")
        return folder

    def run_main(self, args):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return runner.main([str(item) for item in args])

    def test_single_paths_modes_and_conda_are_separate_arguments(self):
        folder = self.sample()
        with patch.object(runner.subprocess, "run", return_value=SimpleNamespace(returncode=0)) as run:
            result = self.run_main([folder, "--single", "--mode", "detailed", "--solver", "gurobi", "--conda-env", "smetana_env"])
        self.assertEqual(result, 0)
        command = run.call_args.args[0]
        self.assertEqual(command[:7], ["conda", "run", "--no-capture-output", "-n", "smetana_env", "smetana", str(folder / "model A.xml")])
        self.assertEqual(command[-5:], ["-c", str(folder / "communities.tsv"), "-d", "--solver", "gurobi"])
        self.assertEqual(run.call_args.kwargs["cwd"], str(folder))
        self.assertNotIn("shell", run.call_args.kwargs)
        self.assertEqual(run.call_args.kwargs["stderr"], subprocess.STDOUT)
        self.assertTrue((folder / "smetana_run.log").is_file())

    def test_dry_run_does_not_execute_or_create_logs(self):
        folder = self.sample()
        with patch.object(runner.subprocess, "run") as run:
            self.assertEqual(self.run_main([folder, "--single", "--dry-run"]), 0)
        run.assert_not_called()
        self.assertFalse((folder / "smetana_run.log").exists())

    def test_batch_continues_after_process_failure(self):
        self.sample("a")
        self.sample("b")
        with patch.object(runner.subprocess, "run", side_effect=[SimpleNamespace(returncode=3), SimpleNamespace(returncode=0)]) as run:
            self.assertEqual(self.run_main([self.root]), 1)
        self.assertEqual(run.call_count, 2)

    def test_batch_continues_after_invalid_sample(self):
        self.sample("a", "community1\tunknown\ncommunity1\tmodel_B\n")
        self.sample("b")
        with patch.object(runner.subprocess, "run", return_value=SimpleNamespace(returncode=0)) as run:
            self.assertEqual(self.run_main([self.root]), 1)
        self.assertEqual(run.call_count, 1)

    def test_reject_header(self):
        folder = self.sample(content="Community\tSpecies\ncommunity1\tmodel A\ncommunity1\tmodel_B\n")
        with self.assertRaisesRegex(ValueError, "headerless"):
            runner.validate_sample(folder)

    def test_reject_unknown_member(self):
        folder = self.sample(content="community1\tmissing\ncommunity1\tmodel_B\n")
        with self.assertRaisesRegex(ValueError, "no matching"):
            runner.validate_sample(folder)

    def test_reject_duplicate_and_single_member_and_extra_column(self):
        for content, message in [("c\tmodel A\nc\tmodel A\n", "duplicate"),
                                 ("c\tmodel A\n", "at least two"),
                                 ("c\tmodel A\textra\n", "two nonempty")]:
            with self.subTest(content=content):
                folder = self.sample(name=str(len(list(self.root.iterdir()))), content=content)
                with self.assertRaisesRegex(ValueError, message):
                    runner.validate_sample(folder)

    def test_groups_and_shared_absolute_community_file(self):
        a = self.sample("CK_XML/a")
        b = self.sample("DL_XML/b")
        self.assertEqual(runner.find_sample_folders(self.root, groups=["CK_XML", "DL_XML"]), [a, b])
        (b / "communities.tsv").unlink()
        with patch.object(runner.subprocess, "run") as run:
            self.assertEqual(self.run_main([self.root, "--groups", "CK_XML", "DL_XML", "--communities", a / "communities.tsv", "--dry-run"]), 0)
        run.assert_not_called()

    def test_missing_root_empty_root_and_missing_group_fail(self):
        self.assertEqual(self.run_main([self.root / "missing"]), 1)
        self.assertEqual(self.run_main([self.root]), 1)
        self.assertEqual(self.run_main([self.root, "--groups", "absent"]), 1)

    def test_single_and_groups_conflict(self):
        self.assertEqual(self.run_main([self.root, "--single", "--groups", "CK_XML"]), 1)


if __name__ == "__main__":
    unittest.main()
