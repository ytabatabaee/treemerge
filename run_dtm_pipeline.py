#!/usr/bin/env python3

"""
===========================================================
Divide-and-Conquer Species Tree Pipeline
(ASTRID + FastME + ASTRAL4 + TreeMerge)
===========================================================

Pipeline:

1. Compute AGID matrix using ASTRID
2. Infer starting tree using FastME
3. Decompose taxa into subsets
4. Prune gene trees to subsets
5. Run ASTRAL4 on each subset
6. Merge subset trees using TreeMerge
7. Score final species tree using ASTRAL4

Input:
    Single file containing estimated gene trees
    (one Newick tree per line)

Output:
    Final scored species tree

===========================================================
"""

import argparse
import multiprocessing as mp
import os
import subprocess
import sys
from pathlib import Path

from ete3 import Tree


# ===========================================================
# Utilities
# ===========================================================

def run(cmd):

    print("\n[RUN]")
    print(" ".join(map(str, cmd)))
    sys.stdout.flush()

    result = subprocess.run(cmd)

    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed:\n{' '.join(map(str, cmd))}"
        )


def mkdir(path):
    os.makedirs(path, exist_ok=True)


# ===========================================================
# ASTRID
# ===========================================================

def run_astrid(
    astrid_bin,
    gene_trees,
    agid_matrix
):

    cmd = [
        astrid_bin,
        "-i", gene_trees,
        "-c", agid_matrix
    ]

    run(cmd)


# ===========================================================
# FastME
# ===========================================================

def infer_starting_tree(
    agid_matrix,
    starting_tree
):

    cmd = [
        "fastme",
        "-i", agid_matrix,
        "-o", starting_tree
    ]

    run(cmd)


# ===========================================================
# Taxon decomposition
# ===========================================================

def recursive_split(
    tree,
    max_subset_size,
    subsets
):

    taxa = tree.get_leaf_names()

    if len(taxa) <= max_subset_size:

        subsets.append(set(taxa))
        return

    children = sorted(
        tree.children,
        key=lambda x: len(x.get_leaf_names()),
        reverse=True
    )

    if len(children) < 2:

        subsets.append(set(taxa))
        return

    recursive_split(
        children[0],
        max_subset_size,
        subsets
    )

    recursive_split(
        children[1],
        max_subset_size,
        subsets
    )


def decompose_tree(
    starting_tree,
    max_subset_size
):

    t = Tree(starting_tree)

    subsets = []

    recursive_split(
        t,
        max_subset_size,
        subsets
    )

    return subsets


# ===========================================================
# Subset gene trees
# ===========================================================

def prune_tree_to_subset(
    tree_str,
    subset_taxa
):

    t = Tree(tree_str)

    keep = [
        x for x in t.get_leaf_names()
        if x in subset_taxa
    ]

    if len(keep) < 4:
        return None

    t.prune(keep)

    return t.write(format=1)


def build_subset_gene_tree_file(
    gene_trees,
    subset_taxa,
    outfile
):

    count = 0

    with open(gene_trees) as infile, \
         open(outfile, "w") as out:

        for line in infile:

            line = line.strip()

            if not line:
                continue

            pruned = prune_tree_to_subset(
                line,
                subset_taxa
            )

            if pruned is not None:

                out.write(pruned + "\n")
                count += 1

    return count


# ===========================================================
# ASTRAL4
# ===========================================================

def run_astral4_subset(job):

    astral4_bin, subset_gene_trees, output_tree = job

    cmd = [
        astral4_bin,
        "-i", subset_gene_trees,
        "-o", output_tree
    ]

    run(cmd)

    return output_tree


# ===========================================================
# TreeMerge
# ===========================================================

def run_treemerge(
    treemerge_script,
    paup_path,
    starting_tree,
    subset_species_trees,
    agid_matrix,
    taxlist,
    output_tree,
    workdir
):

    cmd = [
        "python",
        treemerge_script,

        "-s", starting_tree,

        "-t"
    ]

    cmd.extend(subset_species_trees)

    cmd.extend([
        "-m", agid_matrix,
        "-x", taxlist,
        "-o", output_tree,
        "-w", workdir,
        "-p", paup_path
    ])

    run(cmd)


# ===========================================================
# Final ASTRAL4 scoring
# ===========================================================

def score_species_tree(
    astral4_bin,
    species_tree,
    gene_trees,
    output_tree
):

    cmd = [
        astral4_bin,
        "-q", species_tree,
        "-i", gene_trees,
        "-o", output_tree
    ]

    run(cmd)


# ===========================================================
# Main
# ===========================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--gene_trees",
        required=True,
        help="Single file containing gene trees"
    )

    parser.add_argument(
        "--outdir",
        required=True
    )

    parser.add_argument(
        "--astrid_bin",
        required=True,
        help="Path to ASTRID executable"
    )

    parser.add_argument(
        "--astral4_bin",
        required=True,
        help="Path to astral4 executable"
    )

    parser.add_argument(
        "--treemerge_script",
        required=True,
        help="Path to treemerge.py"
    )

    parser.add_argument(
        "--paup",
        required=True,
        help="Path to PAUP* binary"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=4
    )

    parser.add_argument(
        "--max_subset_size",
        type=int,
        default=100
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    mkdir(outdir)

    intermediate = outdir / "intermediate"
    mkdir(intermediate)

    # =======================================================
    # Step 1: AGID matrix
    # =======================================================

    agid_matrix = intermediate / "agid_matrix.txt"

    run_astrid(
        args.astrid_bin,
        args.gene_trees,
        str(agid_matrix)
    )

    taxlist = str(agid_matrix) + "_taxlist"

    # =======================================================
    # Step 2: Starting tree
    # =======================================================

    starting_tree = intermediate / "starting_tree.nwk"

    infer_starting_tree(
        str(agid_matrix),
        str(starting_tree)
    )

    # =======================================================
    # Step 3: Taxon decomposition
    # =======================================================

    subsets = decompose_tree(
        str(starting_tree),
        args.max_subset_size
    )

    print(f"\nGenerated {len(subsets)} subsets")

    subset_gene_tree_dir = \
        intermediate / "subset_gene_trees"

    mkdir(subset_gene_tree_dir)

    subset_gene_tree_files = []

    for i, subset_taxa in enumerate(subsets):

        outfile = (
            subset_gene_tree_dir /
            f"subset_{i+1}.trees"
        )

        count = build_subset_gene_tree_file(
            args.gene_trees,
            subset_taxa,
            outfile
        )

        print(
            f"[Subset {i+1}] "
            f"{count} pruned gene trees"
        )

        subset_gene_tree_files.append(str(outfile))

    # =======================================================
    # Step 4: Subset ASTRAL4
    # =======================================================

    subset_species_dir = \
        intermediate / "subset_species_trees"

    mkdir(subset_species_dir)

    jobs = []

    subset_species_trees = []

    for i, subset_gt in enumerate(
        subset_gene_tree_files
    ):

        output_tree = (
            subset_species_dir /
            f"subset_{i+1}.nwk"
        )

        jobs.append(
            (
                args.astral4_bin,
                subset_gt,
                str(output_tree)
            )
        )

        subset_species_trees.append(
            str(output_tree)
        )

    pool = mp.Pool(args.threads)

    pool.map(run_astral4_subset, jobs)

    pool.close()
    pool.join()

    # =======================================================
    # Step 5: TreeMerge
    # =======================================================

    merged_tree = (
        outdir /
        "treemerge_species_tree.nwk"
    )

    run_treemerge(
        args.treemerge_script,
        args.paup,
        str(starting_tree),
        subset_species_trees,
        str(agid_matrix),
        taxlist,
        str(merged_tree),
        str(intermediate)
    )

    # =======================================================
    # Step 6: Final ASTRAL4 scoring
    # =======================================================

    final_tree = (
        outdir /
        "final_species_tree_scored.nwk"
    )

    score_species_tree(
        args.astral4_bin,
        str(merged_tree),
        args.gene_trees,
        str(final_tree)
    )

    # =======================================================
    # Finished
    # =======================================================

    print("\n===================================")
    print("PIPELINE COMPLETED")
    print("===================================")

    print("\nFinal species tree:")
    print(final_tree)


if __name__ == "__main__":
    main()
