#!/usr/bin/env python3
"""
SRK_allele_sharing_EOs.py
Phase 3 — Allele Composition Comparison Across Element Occurrences

Compares S-allele composition between Element Occurrences using:
  1. UpSet plot  — all intersection sizes, sorted by set size
  2. Pairwise sharing heatmap — N shared alleles between each EO pair

Inputs (default: current directory):
  SRK_individual_allele_genotypes.tsv  — allele count matrix (individuals x allele bins)
  sampling_metadata.csv                — sample metadata with Pop column

Outputs (--outdir, default: current directory):
  SRK_allele_upset_EOs.pdf
  SRK_allele_sharing_heatmap_EOs.pdf

Usage:
  python SRK_allele_sharing_EOs.py
  python SRK_allele_sharing_EOs.py --genotypes path/to/genotypes.tsv \\
      --metadata path/to/metadata.csv --pops 25 27 67 70 76 --outdir figures/
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
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genotypes", default="SRK_individual_allele_genotypes.tsv",
                        help="Allele count matrix TSV (default: SRK_individual_allele_genotypes.tsv)")
    parser.add_argument("--metadata", default="sampling_metadata.csv",
                        help="Sampling metadata CSV with Pop column (default: sampling_metadata.csv)")
    parser.add_argument("--pops", nargs="+", default=["25", "27", "67", "70", "76"],
                        help="EO population IDs to include (default: 25 27 67 70 76)")
    parser.add_argument("--min-samples", type=int, default=5,
                        help="Minimum individuals to include an EO (default: 5)")
    parser.add_argument("--outdir", default=".",
                        help="Output directory for figures (default: current directory)")
    return parser.parse_args()

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_data(genotypes_path, metadata_path):
    geno = pd.read_csv(genotypes_path, sep="\t", index_col=0, encoding="utf-8-sig")
    meta = pd.read_csv(metadata_path, encoding="utf-8-sig")
    meta["Pop"] = meta["Pop"].astype(str).str.strip()
    # Use SampleID or first column that matches genotype index
    id_col = "SampleID" if "SampleID" in meta.columns else meta.columns[0]
    meta = meta.set_index(id_col)
    return geno, meta

# ---------------------------------------------------------------------------
# Build per-EO allele sets
# ---------------------------------------------------------------------------

def build_eo_allele_sets(geno, meta, pops, min_samples):
    allele_cols = [c for c in geno.columns if c.startswith("Allele_")]
    eo_sets = {}
    eo_sizes = {}
    for pop in pops:
        individuals = meta.index[meta["Pop"] == pop].tolist()
        # Keep only individuals present in genotype matrix
        individuals = [i for i in individuals if i in geno.index]
        if "Ingroup" in meta.columns:
            individuals = [i for i in individuals
                           if str(meta.loc[i, "Ingroup"]).strip() in ("1", "True", "yes")]
        if len(individuals) < min_samples:
            print(f"  Skipping EO{pop}: only {len(individuals)} individuals (< {min_samples})")
            continue
        sub = geno.loc[individuals, allele_cols]
        # Allele present in an EO if at least one individual has count > 0
        col_sums = sub.values.sum(axis=0)
        present = set(c for c, s in zip(allele_cols, col_sums) if s > 0)
        eo_sets[f"EO{pop}"] = present
        eo_sizes[f"EO{pop}"] = len(individuals)
        print(f"  EO{pop}: {len(individuals)} individuals, {len(present)} alleles")
    return eo_sets

# ---------------------------------------------------------------------------
# UpSet plot (pure matplotlib)
# ---------------------------------------------------------------------------

EO_COLORS = {
    "EO25": "#4C72B0",
    "EO27": "#DD8452",
    "EO67": "#55A868",
    "EO70": "#C44E52",
    "EO76": "#8172B2",
}

def compute_intersections(eo_sets):
    """Return list of (frozenset_of_eos, allele_set) for every non-empty intersection,
    where the intersection is the set of alleles EXCLUSIVE to exactly those EOs."""
    eos = sorted(eo_sets.keys())
    intersections = []
    for r in range(1, len(eos) + 1):
        for combo in itertools.combinations(eos, r):
            combo_set = frozenset(combo)
            # Alleles in ALL members of combo
            shared = set.intersection(*(eo_sets[e] for e in combo))
            # Alleles NOT in any EO outside combo
            others = set(eos) - combo_set
            if others:
                in_others = set.union(*(eo_sets[e] for e in others))
                exclusive = shared - in_others
            else:
                exclusive = shared
            if exclusive:
                intersections.append((combo_set, exclusive))
    # Sort by size descending
    intersections.sort(key=lambda x: len(x[1]), reverse=True)
    return intersections, eos


def plot_upset(eo_sets, outpath):
    intersections, eos = compute_intersections(eo_sets)
    n_inter = len(intersections)
    n_eos   = len(eos)

    counts  = [len(alleles) for _, alleles in intersections]
    combos  = [combo for combo, _ in intersections]

    # Layout: top bar chart + bottom dot matrix
    bar_height_ratio  = 3
    dot_height_ratio  = n_eos
    fig_height = 4 + n_eos * 0.45
    fig_width  = max(8, n_inter * 0.55 + 2)

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs  = fig.add_gridspec(2, 1,
                           height_ratios=[bar_height_ratio, dot_height_ratio],
                           hspace=0.05)
    ax_bar = fig.add_subplot(gs[0])
    ax_dot = fig.add_subplot(gs[1], sharex=ax_bar)

    x = np.arange(n_inter)

    # --- Bar chart ---
    bar_colors = []
    for combo in combos:
        if len(combo) == 1:
            eo = list(combo)[0]
            bar_colors.append(EO_COLORS.get(eo, "#888888"))
        else:
            bar_colors.append("#555555")

    bars = ax_bar.bar(x, counts, color=bar_colors, edgecolor="white", linewidth=0.5, zorder=3)
    for bar, count in zip(bars, counts):
        ax_bar.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05,
                    str(count), ha="center", va="bottom", fontsize=8)

    ax_bar.set_ylabel("Allele count", fontsize=10)
    ax_bar.set_ylim(0, max(counts) * 1.18)
    ax_bar.spines[["right", "top", "bottom"]].set_visible(False)
    ax_bar.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_bar.yaxis.grid(True, linestyle="--", alpha=0.5, zorder=0)
    ax_bar.set_axisbelow(True)
    ax_bar.set_title("S-allele sharing among Element Occurrences (UpSet plot)",
                     fontsize=12, pad=10)

    # --- Dot matrix ---
    y_pos = {eo: i for i, eo in enumerate(reversed(eos))}

    ax_dot.set_xlim(-0.5, n_inter - 0.5)
    ax_dot.set_ylim(-0.5, n_eos - 0.5)

    # Background alternating rows
    for i, eo in enumerate(reversed(eos)):
        color = "#f0f0f0" if i % 2 == 0 else "white"
        ax_dot.axhspan(i - 0.5, i + 0.5, facecolor=color, zorder=0)

    # Dots and connecting lines
    dot_size = 120
    for xi, combo in enumerate(combos):
        active_ys = []
        for eo in eos:
            yi = y_pos[eo]
            if eo in combo:
                fc = EO_COLORS.get(eo, "#555555")
                ax_dot.scatter(xi, yi, s=dot_size, color=fc, zorder=4, linewidths=0)
                active_ys.append(yi)
            else:
                ax_dot.scatter(xi, yi, s=dot_size, color="#dddddd", zorder=3, linewidths=0)
        # Vertical line connecting active dots
        if len(active_ys) > 1:
            ax_dot.plot([xi, xi], [min(active_ys), max(active_ys)],
                        color="#333333", linewidth=2, zorder=3.5)

    ax_dot.set_yticks(range(n_eos))
    ax_dot.set_yticklabels(list(reversed(eos)), fontsize=10)
    ax_dot.spines[["right", "top", "left", "bottom"]].set_visible(False)
    ax_dot.tick_params(axis="both", which="both", length=0)
    ax_dot.set_xticks([])

    # Legend patches
    legend_patches = [mpatches.Patch(color=EO_COLORS.get(eo, "#888888"), label=eo)
                      for eo in eos]
    ax_bar.legend(handles=legend_patches, loc="upper right",
                  fontsize=8, framealpha=0.9, ncol=min(3, len(eos)))

    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")

# ---------------------------------------------------------------------------
# Pairwise sharing heatmap
# ---------------------------------------------------------------------------

def plot_heatmap(eo_sets, outpath):
    eos    = sorted(eo_sets.keys())
    n      = len(eos)
    matrix = np.zeros((n, n), dtype=int)

    for i, e1 in enumerate(eos):
        for j, e2 in enumerate(eos):
            matrix[i, j] = len(eo_sets[e1] & eo_sets[e2])

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(matrix, cmap="YlOrRd", aspect="equal")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(eos, fontsize=11)
    ax.set_yticklabels(eos, fontsize=11)

    for i in range(n):
        for j in range(n):
            val = matrix[i, j]
            # Choose text colour for contrast
            threshold = matrix.max() * 0.6
            color = "white" if val > threshold else "black"
            ax.text(j, i, str(val), ha="center", va="center",
                    fontsize=13, fontweight="bold", color=color)

    cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("Shared alleles (N)", fontsize=10)

    ax.set_title("Pairwise S-allele sharing between Element Occurrences",
                 fontsize=12, pad=12)

    # Diagonal label
    ax.set_xlabel("Element Occurrence", fontsize=10)
    ax.set_ylabel("Element Occurrence", fontsize=10)

    plt.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def print_summary(eo_sets, intersections):
    eos = sorted(eo_sets.keys())
    print("\n--- Allele set sizes ---")
    for eo in eos:
        print(f"  {eo}: {len(eo_sets[eo])} alleles")

    print("\n--- Private alleles (exclusive to one EO) ---")
    for combo, alleles in intersections:
        if len(combo) == 1:
            eo = list(combo)[0]
            print(f"  {eo}: {len(alleles)} private alleles  ({sorted(alleles)})")

    print("\n--- Alleles shared across all EOs ---")
    shared_all = set.intersection(*eo_sets.values())
    if shared_all:
        print(f"  {len(shared_all)} alleles shared by all {len(eos)} EOs: {sorted(shared_all)}")
    else:
        print("  None")

    print("\n--- Top intersections ---")
    for combo, alleles in intersections[:10]:
        label = " ∩ ".join(sorted(combo))
        print(f"  {label}: {len(alleles)}")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    pops = [str(p) for p in args.pops]

    print("Loading data...")
    geno, meta = load_data(args.genotypes, args.metadata)

    print(f"\nBuilding allele sets for EOs: {pops}")
    eo_sets = build_eo_allele_sets(geno, meta, pops, args.min_samples)

    if len(eo_sets) < 2:
        sys.exit("ERROR: Fewer than 2 EOs with sufficient data. Check --pops and --min-samples.")

    intersections, _ = compute_intersections(eo_sets)
    print_summary(eo_sets, intersections)

    print("\nGenerating figures...")
    upset_path   = os.path.join(args.outdir, "SRK_allele_upset_EOs.pdf")
    heatmap_path = os.path.join(args.outdir, "SRK_allele_sharing_heatmap_EOs.pdf")

    plot_upset(eo_sets, upset_path)
    plot_heatmap(eo_sets, heatmap_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
