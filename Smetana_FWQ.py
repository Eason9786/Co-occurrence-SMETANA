#!/usr/bin/env python3
"""Validate community inputs and run SMETANA in one or more model folders.

Activate the environment containing SMETANA before running this script, or use
--conda-env. Each sample contains .xml models and a headerless communities.tsv.
"""

import argparse
import csv
import shlex
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def find_sample_folders(root, single=False, groups=None):
    """Return immediate sample folders, optionally within named group folders."""
    root = Path(root).resolve()
    if not root.is_dir():
        raise ValueError("Root directory does not exist: {}".format(root))
    if single:
        if groups:
            raise ValueError("--single and --groups cannot be used together")
        return [root]

    parents = [root / group for group in groups] if groups else [root]
    samples = []
    for parent in parents:
        if not parent.is_dir():
            raise ValueError("Group directory does not exist: {}".format(parent))
        children = sorted(path for path in parent.iterdir() if path.is_dir())
        if not children:
            raise ValueError("No sample subdirectories found in: {}".format(parent))
        samples.extend(children)
    # Repeated group names must not cause duplicate executions.
    return list(dict.fromkeys(samples))


def validate_sample(sample, communities="communities.tsv"):
    """Validate the two-column TSV and return model and community file paths.

    This checks names and membership, not SBML content or model feasibility.
    SMETANA reads community files with header=None, so no header is accepted.
    """
    sample = Path(sample).resolve()
    models = sorted(path for path in sample.glob("*.xml") if path.is_file())
    if len(models) < 2:
        raise ValueError("At least two .xml model files are required in {}".format(sample))
    model_ids = {path.stem for path in models}
    community_file = Path(communities)
    if not community_file.is_absolute():
        community_file = sample / community_file
    community_file = community_file.resolve()
    if not community_file.is_file():
        raise ValueError("Community file does not exist: {}".format(community_file))

    members = defaultdict(set)
    with community_file.open("r", encoding="utf-8-sig", newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(row) != 2 or any(not value.strip() for value in row):
                raise ValueError("{}:{}: expected two nonempty tab-separated columns".format(
                    community_file, line_number))
            community, species = row
            if line_number == 1 and [value.lower() for value in row] == ["community", "species"]:
                raise ValueError("{}: remove the Community/Species header; SMETANA requires a headerless TSV".format(
                    community_file))
            if species not in model_ids:
                raise ValueError("{}:{}: no matching .xml model for species {!r}".format(
                    community_file, line_number, species))
            if species in members[community]:
                raise ValueError("{}:{}: duplicate species {!r} in community {!r}".format(
                    community_file, line_number, species, community))
            members[community].add(species)
    if not members:
        raise ValueError("Community file is empty: {}".format(community_file))
    for community, species in members.items():
        if len(species) < 2:
            raise ValueError("Community {!r} must contain at least two distinct species".format(community))
    return models, community_file


def build_command(models, community_file, mode="global", solver=None, conda_env=None):
    """Build an argument list; do not rely on shell expansion or activation."""
    if mode not in ("global", "detailed"):
        raise ValueError("Unsupported SMETANA mode: {}".format(mode))
    if solver not in (None, "gurobi", "cplex"):
        raise ValueError("Unsupported solver: {}".format(solver))
    command = (["conda", "run", "--no-capture-output", "-n", conda_env, "smetana"]
               if conda_env else ["smetana"])
    command.extend(str(path) for path in models)
    command.extend(["-c", str(community_file), "-g" if mode == "global" else "-d"])
    if solver:
        command.extend(["--solver", solver])
    return command


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, help="Root containing sample subfolders (or the model folder with --single)")
    parser.add_argument("--single", action="store_true", help="Run directly in root instead of its sample subfolders")
    parser.add_argument("--groups", nargs="+", help="Group subfolders to search, for example CK_XML DL_XML")
    parser.add_argument("--communities", default="communities.tsv", help="TSV path relative to each sample, or an absolute shared path")
    parser.add_argument("--conda-env", help="Run SMETANA using conda run -n ENV")
    parser.add_argument("--mode", choices=("global", "detailed"), default="global")
    parser.add_argument("--solver", choices=("gurobi", "cplex"), help="Forward solver selection to SMETANA")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs and print commands without running SMETANA")
    args = parser.parse_args(argv)

    try:
        samples = find_sample_folders(args.root, args.single, args.groups)
    except (ValueError, OSError) as exc:
        print("Error: {}".format(exc), file=sys.stderr)
        return 1

    failed = 0
    for sample in samples:
        try:
            models, community_file = validate_sample(sample, args.communities)
            command = build_command(models, community_file, args.mode, args.solver, args.conda_env)
            print("Working directory: {}".format(sample))
            print("Command (display only): {}".format(shlex.join(command)))
            if args.dry_run:
                continue
            log_path = sample / "smetana_run.log"
            with log_path.open("w", encoding="utf-8") as log:
                result = subprocess.run(command, cwd=str(sample), stdout=log,
                                        stderr=subprocess.STDOUT, check=False)
            if result.returncode:
                raise RuntimeError("SMETANA exited with code {}; see {}".format(result.returncode, log_path))
            print("Completed: {} (log: {})".format(sample, log_path))
        except (ValueError, OSError, UnicodeError, csv.Error, RuntimeError) as exc:
            failed += 1
            print("Failed [{}]: {}".format(sample, exc), file=sys.stderr)
    print("{}: {} sample(s), {} failed".format(
        "Validation complete" if args.dry_run else "Batch complete", len(samples), failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
