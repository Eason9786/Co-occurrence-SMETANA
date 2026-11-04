"""Sample balanced subcommunities from one folder of SBML XML model files.

This script samples combinations; it does not perform statistical co-occurrence
tests. Model identifiers are the filenames without the final ``.xml`` suffix.
"""

import argparse
import csv
import math
from pathlib import Path
import random
import sys


def _allocate_sizes(capacities, total_budget):
    """Distribute the budget equally, redistributing any exhausted capacity."""
    remaining = min(total_budget, sum(capacities.values()))
    active = sorted(capacities, key=lambda k: (capacities[k], k))
    allocated = {}
    while active:
        share, extra = divmod(remaining, len(active))
        size = active[0]
        if capacities[size] <= share:
            allocated[size] = capacities[size]
            remaining -= capacities[size]
            active.pop(0)
        else:
            for index, size in enumerate(sorted(active)):
                allocated[size] = share + (index < extra)
            break
    return allocated


def _sample_ranks(rng, population_size, sample_size):
    """Floyd's algorithm: unique ranks, without materializing the population.

    Unlike sampling a range through random.sample, this also supports population
    sizes larger than the platform's maximum sequence length.
    """
    selected = set()
    for upper in range(population_size - sample_size, population_size):
        candidate = rng.randrange(upper + 1)
        selected.add(upper if candidate in selected else candidate)
    return sorted(selected)


def _unrank_combination(species, size, rank):
    """Decode a zero-based lexicographic combination rank."""
    members = []
    start = 0
    for remaining in range(size, 0, -1):
        for index in range(start, len(species) - remaining + 1):
            block_size = math.comb(len(species) - index - 1, remaining - 1)
            if rank < block_size:
                members.append(species[index])
                start = index + 1
                break
            rank -= block_size
    return frozenset(members)


def generate_balanced_sub_communities(folder_path, max_combo_size=5,
                                     total_budget=1000, seed=42):
    """Return unique subcommunities, balanced across feasible sizes 2..max.

    At most ``total_budget`` communities are returned. If fewer combinations
    exist, all are returned. The seed controls sampling and output order without
    changing Python's global random state. The folder is read nonrecursively;
    only regular files ending in lowercase ``.xml`` are included.
    """
    if not isinstance(max_combo_size, int) or max_combo_size < 2:
        raise ValueError("max_combo_size must be an integer of at least 2")
    if not isinstance(total_budget, int) or total_budget < 1:
        raise ValueError("total_budget must be a positive integer")
    folder = Path(folder_path)
    if not folder.is_dir():
        raise ValueError("Model folder does not exist: {}".format(folder))

    all_species = sorted(path.stem for path in folder.iterdir()
                         if path.is_file() and path.name.endswith(".xml"))
    if any(not name or any(char in name for char in "\t\r\n")
           for name in all_species):
        raise ValueError("Model filenames must have nonempty IDs without tabs or newlines")
    if len(all_species) < 2:
        return []

    capacities = {size: math.comb(len(all_species), size)
                  for size in range(2, min(max_combo_size, len(all_species)) + 1)}
    allocated = _allocate_sizes(capacities, total_budget)
    rng = random.Random(seed)
    communities = []
    for size in sorted(allocated):
        for rank in _sample_ranks(rng, capacities[size], allocated[size]):
            communities.append(_unrank_combination(all_species, size, rank))
    rng.shuffle(communities)
    return communities


def _community_rows(communities):
    for index, combo in enumerate(communities, start=1):
        for species in sorted(combo):
            yield ("community{}".format(index), species)


def generate_tsv_output(co_occurrences, output_filename="communities.tsv"):
    """Write SMETANA's two-column community definition, WITHOUT a header."""
    output = Path(output_filename)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(_community_rows(co_occurrences))


def generate_excel_output_simplified(co_occurrences,
                                     output_filename="cooccurrence_communities.xlsx"):
    """Write an optional human-readable Excel table with a column header."""
    import openpyxl

    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.title = "Co-occurring Communities"
    sheet.append(["Community", "Species"])
    for row in _community_rows(co_occurrences):
        sheet.append(row)
    output = Path(output_filename)
    output.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(output)
    workbook.close()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model_folder", help="Folder containing species .xml model files")
    parser.add_argument("--max-combo-size", type=int, default=5)
    parser.add_argument("--total-budget", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", default="communities.tsv",
                        help="Headerless SMETANA TSV output (default: communities.tsv)")
    parser.add_argument("--excel", help="Optional Excel output path (requires openpyxl)")
    args = parser.parse_args(argv)
    if args.excel and Path(args.output).resolve() == Path(args.excel).resolve():
        parser.error("TSV and Excel output paths must differ")
    try:
        communities = generate_balanced_sub_communities(
            args.model_folder, args.max_combo_size, args.total_budget, args.seed)
        if not communities:
            raise ValueError("At least two .xml model files are required")
        generate_tsv_output(communities, args.output)
        if args.excel:
            generate_excel_output_simplified(communities, args.excel)
    except (OSError, ValueError, ImportError) as error:
        parser.exit(1, "Error: {}\n".format(error))
    print("Saved {} communities to {}".format(len(communities), args.output))
    if args.excel:
        print("Saved Excel table to {}".format(args.excel))
    return 0


if __name__ == "__main__":
    sys.exit(main())
