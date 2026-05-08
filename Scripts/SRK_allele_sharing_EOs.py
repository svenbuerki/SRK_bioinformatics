#!/usr/bin/env python3
"""
SRK_allele_sharing_EOs.py - Step 18 of Canu_amplicon pipeline

Compares S-allele composition between groups using:
  1. UpSet plot - all intersection sizes, sorted by set size
  2. Pairwise sharing heatmap - N shared alleles between each group pair

Two stratification levels are produced in parallel:
  - Element Occurrences (EOs) with N >= 5 individuals, sorted by parent BL
  - Bottleneck Lineages (BL1-BL5), all aggregated

Both levels share the locked Set1 BL palette used by Steps 14-17 and the
sibling LEPA_EO_spatial_clustering project.

Inputs:
  SRK_individual_allele_genotypes.tsv      - allele count matrix
  SRK_individual_BL_assignments.tsv        - from Step 13 (EO/BL labels)

Outputs (PDF in project root, PNG in figures/):
  SRK_allele_upset_EOs.{pdf,png}           - EO-level UpSet
  SRK_allele_sharing_heatmap_EOs.{pdf,png} - EO-level pairwise heatmap
  SRK_allele_upset_BLs.{pdf,png}           - BL-level UpSet (NEW)
  SRK_allele_sharing_heatmap_BLs.{pdf,png} - BL-level pairwise heatmap (NEW)
"""

import argparse
import itertools
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Locked BL color palette (Set1; matches Steps 14-17 and spatial clustering)
# ---------------------------------------------------------------------------

BL_PALETTE = {
    "BL1": "#E41A1C",  # red
    "BL2": "#377EB8",  # blue
    "BL3": "#4DAF4A",  # green
    "BL4": "#984EA3",  # purple
    "BL5": "#FF7F00",  # orange
}
BL_UNASSIGNED_COLOR = "#999999"


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--genotypes",
        default="SRK_individual_allele_genotypes.tsv",
        help="Allele count matrix TSV",
    )
    p.add_argument(
        "--bl-assignments",
        default="SRK_individual_BL_assignments.tsv",
        help="Per-individual EO/BL bridge table from Step 13",
    )
    p.add_argument(
        "--min-samples",
        type=int,
        default=5,
        help="Minimum individuals to include an EO (default: 5)",
    )
    p.add_argument(
        "--pdf-dir",
        default=".",
        help="Output directory for PDFs (default: project root)",
    )
    p.add_argument(
        "--png-dir",
        default="figures",
        help="Output directory for PNGs (default: figures/)",
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_data(genotypes_path, bl_path):
    geno = pd.read_csv(
        genotypes_path, sep="\t", index_col=0, encoding="utf-8-sig"
    )
    bl = pd.read_csv(bl_path, sep="\t", encoding="utf-8-sig")
    bl = bl[bl["BL_status"].isin(["Assigned", "Inferred"])].copy()
    bl = bl.set_index("Individual")
    return geno, bl


# ---------------------------------------------------------------------------
# Build per-group allele sets
# ---------------------------------------------------------------------------


def build_allele_sets(geno, bl_df, group_col, min_samples=1):
    """Return dict {group_label: set(allele_ids)}.
    Only groups with >= min_samples individuals are kept.
    Order: by parent BL, then alphabetical within BL.
    """
    allele_cols = [c for c in geno.columns if c.startswith("Allele_")]
    sets = {}
    sizes = {}

    # Determine ordering
    if group_col == "EO":
        # Group EOs by parent BL, then sort alphabetically within each BL
        # so the UpSet column order reflects BL nesting.
        eo_to_bl = bl_df.groupby("EO")["BL"].first()
        order = sorted(
            eo_to_bl.index, key=lambda eo: (eo_to_bl[eo], eo)
        )
    else:  # BL
        order = sorted(bl_df["BL"].unique())

    for group in order:
        individuals = bl_df.index[bl_df[group_col] == group].tolist()
        individuals = [i for i in individuals if i in geno.index]
        if len(individuals) < min_samples:
            print(
                f"  Skipping {group}: only {len(individuals)} individuals "
                f"(< {min_samples})"
            )
            continue
        sub = geno.loc[individuals, allele_cols]
        col_sums = sub.values.sum(axis=0)
        present = {c for c, s in zip(allele_cols, col_sums) if s > 0}
        sets[group] = present
        sizes[group] = len(individuals)
        print(
            f"  {group}: {len(individuals)} individuals, "
            f"{len(present)} alleles"
        )

    return sets, sizes


# ---------------------------------------------------------------------------
# Resolve parent BL for an EO label (used to color EO bars/dots by BL)
# ---------------------------------------------------------------------------


def eo_to_bl_color(label, bl_df):
    if label.startswith("BL"):
        return BL_PALETTE.get(label, BL_UNASSIGNED_COLOR)
    bl = bl_df.loc[bl_df["EO"] == label, "BL"]
    if bl.empty:
        return BL_UNASSIGNED_COLOR
    return BL_PALETTE.get(bl.iloc[0], BL_UNASSIGNED_COLOR)


# ---------------------------------------------------------------------------
# Intersections
# ---------------------------------------------------------------------------


def compute_intersections(group_sets):
    """Return list of (frozenset_of_groups, allele_set) for every non-empty
    intersection where alleles are EXCLUSIVE to exactly that combo of groups.
    """
    groups = list(group_sets.keys())  # preserve insertion order
    intersections = []
    for r in range(1, len(groups) + 1):
        for combo in itertools.combinations(groups, r):
            combo_set = frozenset(combo)
            shared = set.intersection(*(group_sets[g] for g in combo))
            others = set(groups) - combo_set
            if others:
                in_others = set.union(*(group_sets[g] for g in others))
                exclusive = shared - in_others
            else:
                exclusive = shared
            if exclusive:
                intersections.append((combo_set, exclusive))
    intersections.sort(key=lambda x: len(x[1]), reverse=True)
    return intersections, groups


# ---------------------------------------------------------------------------
# UpSet plot
# ---------------------------------------------------------------------------


def plot_upset(group_sets, sizes, bl_df, title, level_name, save_to):
    intersections, groups = compute_intersections(group_sets)
    n_inter = len(intersections)
    n_grp = len(groups)

    counts = [len(alleles) for _, alleles in intersections]
    combos = [combo for combo, _ in intersections]

    bar_height_ratio = 3
    dot_height_ratio = n_grp
    fig_height = 4 + n_grp * 0.45
    fig_width = max(8, n_inter * 0.55 + 2)

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(
        2, 1,
        height_ratios=[bar_height_ratio, dot_height_ratio],
        hspace=0.05,
    )
    ax_bar = fig.add_subplot(gs[0])
    ax_dot = fig.add_subplot(gs[1], sharex=ax_bar)

    x = np.arange(n_inter)

    # Group color = parent BL (for EOs) or own BL (for BLs)
    group_color = {g: eo_to_bl_color(g, bl_df) for g in groups}

    # Bar color: single-group bars get that group's BL color; multi-group
    # bars stay neutral grey.
    bar_colors = []
    for combo in combos:
        if len(combo) == 1:
            bar_colors.append(group_color[next(iter(combo))])
        else:
            bar_colors.append("#555555")

    bars = ax_bar.bar(
        x, counts, color=bar_colors, edgecolor="white",
        linewidth=0.5, zorder=3,
    )
    for bar, count in zip(bars, counts):
        ax_bar.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.05,
            str(count),
            ha="center", va="bottom", fontsize=8,
        )

    ax_bar.set_ylabel("Allele count", fontsize=10)
    ax_bar.set_ylim(0, max(counts) * 1.18)
    ax_bar.spines[["right", "top", "bottom"]].set_visible(False)
    ax_bar.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_bar.yaxis.grid(True, linestyle="--", alpha=0.5, zorder=0)
    ax_bar.set_axisbelow(True)
    ax_bar.set_title(title, fontsize=12, pad=10)

    # Dot matrix (groups in input order, but plotted top-to-bottom)
    y_pos = {g: i for i, g in enumerate(reversed(groups))}

    ax_dot.set_xlim(-0.5, n_inter - 0.5)
    ax_dot.set_ylim(-0.5, n_grp - 0.5)

    for i, g in enumerate(reversed(groups)):
        color = "#f0f0f0" if i % 2 == 0 else "white"
        ax_dot.axhspan(i - 0.5, i + 0.5, facecolor=color, zorder=0)

    dot_size = 120
    for xi, combo in enumerate(combos):
        active_ys = []
        for g in groups:
            yi = y_pos[g]
            if g in combo:
                ax_dot.scatter(
                    xi, yi, s=dot_size, color=group_color[g],
                    zorder=4, linewidths=0,
                )
                active_ys.append(yi)
            else:
                ax_dot.scatter(
                    xi, yi, s=dot_size, color="#dddddd",
                    zorder=3, linewidths=0,
                )
        if len(active_ys) > 1:
            ax_dot.plot(
                [xi, xi], [min(active_ys), max(active_ys)],
                color="#333333", linewidth=2, zorder=3.5,
            )

    # Y labels include sample size
    ylabels = [f"{g}  (n={sizes.get(g, 0)})" for g in reversed(groups)]
    ax_dot.set_yticks(range(n_grp))
    ax_dot.set_yticklabels(ylabels, fontsize=10)
    ax_dot.spines[["right", "top", "left", "bottom"]].set_visible(False)
    ax_dot.tick_params(axis="both", which="both", length=0)
    ax_dot.set_xticks([])

    # Legend: BL colors used in this plot
    bls_in_plot = sorted(
        {bl_df.loc[bl_df["EO"] == g, "BL"].iloc[0]
         if g in bl_df["EO"].values else g
         for g in groups
         if g.startswith("BL") or g in bl_df["EO"].values}
    )
    legend_patches = [
        mpatches.Patch(color=BL_PALETTE.get(b, BL_UNASSIGNED_COLOR), label=b)
        for b in bls_in_plot
    ]
    ax_bar.legend(
        handles=legend_patches, loc="upper right",
        fontsize=8, framealpha=0.9, ncol=min(3, len(bls_in_plot)),
        title="Bottleneck lineage",
    )

    for path in save_to:
        fig.savefig(path, dpi=200, bbox_inches="tight")
        print(f"  Saved: {path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Pairwise sharing heatmap
# ---------------------------------------------------------------------------


def plot_heatmap(group_sets, title, level_name, save_to):
    groups = list(group_sets.keys())
    n = len(groups)
    matrix = np.zeros((n, n), dtype=int)

    for i, g1 in enumerate(groups):
        for j, g2 in enumerate(groups):
            matrix[i, j] = len(group_sets[g1] & group_sets[g2])

    fig, ax = plt.subplots(figsize=(max(5, 0.85 * n + 2), max(4, 0.85 * n + 1)))
    im = ax.imshow(matrix, cmap="YlOrRd", aspect="equal")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(groups, fontsize=11)
    ax.set_yticklabels(groups, fontsize=11)

    threshold = matrix.max() * 0.6 if matrix.max() > 0 else 0
    for i in range(n):
        for j in range(n):
            val = matrix[i, j]
            color = "white" if val > threshold else "black"
            ax.text(
                j, i, str(val), ha="center", va="center",
                fontsize=13, fontweight="bold", color=color,
            )

    cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("Shared alleles (N)", fontsize=10)

    ax.set_title(title, fontsize=12, pad=12)
    ax.set_xlabel(level_name, fontsize=10)
    ax.set_ylabel(level_name, fontsize=10)

    plt.tight_layout()
    for path in save_to:
        fig.savefig(path, dpi=200, bbox_inches="tight")
        print(f"  Saved: {path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------


def print_summary(group_sets, intersections, level_name):
    groups = list(group_sets.keys())
    print(f"\n--- Allele set sizes ({level_name}) ---")
    for g in groups:
        print(f"  {g}: {len(group_sets[g])} alleles")

    print(f"\n--- Private alleles (exclusive to one {level_name}) ---")
    for combo, alleles in intersections:
        if len(combo) == 1:
            g = next(iter(combo))
            print(
                f"  {g}: {len(alleles)} private alleles  "
                f"({sorted(alleles)})"
            )

    print(f"\n--- Alleles shared across all {level_name}s ---")
    if group_sets:
        shared_all = set.intersection(*group_sets.values())
        if shared_all:
            print(
                f"  {len(shared_all)} alleles shared by all "
                f"{len(groups)} {level_name}s: {sorted(shared_all)}"
            )
        else:
            print("  None")

    print(f"\n--- Top intersections ({level_name}) ---")
    for combo, alleles in intersections[:10]:
        label = " ∩ ".join(combo)
        print(f"  {label}: {len(alleles)}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    args = parse_args()
    os.makedirs(args.pdf_dir, exist_ok=True)
    os.makedirs(args.png_dir, exist_ok=True)

    print("Loading data (Step 18)...")
    geno, bl_df = load_data(args.genotypes, args.bl_assignments)
    print(
        f"  {len(geno)} individuals in genotype matrix, "
        f"{len(bl_df)} BL-assigned individuals"
    )

    # ---- EO level ----
    print("\nBuilding EO allele sets (sorted by parent BL, N >= "
          f"{args.min_samples})...")
    eo_sets, eo_sizes = build_allele_sets(
        geno, bl_df, "EO", args.min_samples
    )
    if len(eo_sets) >= 2:
        eo_intersections, _ = compute_intersections(eo_sets)
        print_summary(eo_sets, eo_intersections, "EO")

        print("\nGenerating EO figures...")
        plot_upset(
            eo_sets, eo_sizes, bl_df,
            title=("S-allele sharing among Element Occurrences "
                   "(UpSet, EOs sorted by parent BL)"),
            level_name="EO",
            save_to=[
                os.path.join(args.pdf_dir, "SRK_allele_upset_EOs.pdf"),
                os.path.join(args.png_dir, "SRK_allele_upset_EOs.png"),
            ],
        )
        plot_heatmap(
            eo_sets,
            title="Pairwise S-allele sharing between Element Occurrences",
            level_name="Element Occurrence",
            save_to=[
                os.path.join(args.pdf_dir, "SRK_allele_sharing_heatmap_EOs.pdf"),
                os.path.join(args.png_dir, "SRK_allele_sharing_heatmap_EOs.png"),
            ],
        )
    else:
        print("\nFewer than 2 EOs with sufficient data - skipping EO plots.")

    # ---- BL level ----
    print("\nBuilding BL allele sets (all 5 BLs, no minimum)...")
    bl_sets, bl_sizes = build_allele_sets(
        geno, bl_df, "BL", min_samples=1
    )
    if len(bl_sets) >= 2:
        bl_intersections, _ = compute_intersections(bl_sets)
        print_summary(bl_sets, bl_intersections, "BL")

        print("\nGenerating BL figures...")
        plot_upset(
            bl_sets, bl_sizes, bl_df,
            title="S-allele sharing among Bottleneck Lineages (UpSet)",
            level_name="BL",
            save_to=[
                os.path.join(args.pdf_dir, "SRK_allele_upset_BLs.pdf"),
                os.path.join(args.png_dir, "SRK_allele_upset_BLs.png"),
            ],
        )
        plot_heatmap(
            bl_sets,
            title="Pairwise S-allele sharing between Bottleneck Lineages",
            level_name="Bottleneck Lineage",
            save_to=[
                os.path.join(args.pdf_dir, "SRK_allele_sharing_heatmap_BLs.pdf"),
                os.path.join(args.png_dir, "SRK_allele_sharing_heatmap_BLs.png"),
            ],
        )
    else:
        print("\nFewer than 2 BLs - skipping BL plots.")

    print("\nStep 18 complete.")


if __name__ == "__main__":
    main()
