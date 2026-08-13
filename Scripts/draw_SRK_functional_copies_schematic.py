#!/usr/bin/env python3
"""Schematic showing the mechanistic link between SRK functional copy number
(molecular Step 22b) and breeding-system classification (ISI phenotypic).

LEPA is an allopolyploid tetraploid — two subgenomes (A + B), each contributing
two chromosomes carrying the S-locus. This figure walks the five possible
functional-copy states (4 → 3 → 2 → 1 → 0) as five side-by-side panels.

Each panel:
    - 4 chromosomes (2 A + 2 B) drawn as rounded grey bars
    - S-locus on each chromosome shown as a coloured circle:
        green = functional SRK, red = non-functional
    - Category label below: SI / pSI / SC in the ISI-axis palette
      (SI = red, Partial-SI = amber, SC = green)

Two colour scales are used on purpose:
    - Loci follow the MOLECULAR convention (green = working).
    - Category labels follow the ISI-axis convention (green = SC band).
A caption note at the bottom of the figure spells this out.

Output:
    figures/SRK_functional_copies_schematic.{pdf,png}
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch

DEFAULT_FIGURES = "figures"

# --- Locus colours (deliberately outside the ISI red/amber/green palette
#     to avoid semantic collision with the category labels below) ---
COL_FUNCTIONAL = "#2166ac"       # medium blue — functional SRK allele
COL_BROKEN     = "#4d4d4d"       # dark grey  — non-functional SRK allele
COL_LOCUS_EDGE = "black"

# --- Chromosome colours ---
COL_CHROM_A = "#B9CBDB"          # pale blue-grey for subgenome A
COL_CHROM_B = "#D8CBE0"          # pale mauve for subgenome B
COL_CHROM_EDGE = "#555555"

# --- Category-label palette (ISI axis) ---
COL_SI_LABEL  = "#b2182b"        # red — SI band
COL_PSI_LABEL = "#e08214"        # amber — Partial-SI band
COL_SC_LABEL  = "#1b7837"        # green — SC band

# --- Five scenarios, arranged left-to-right to align with the ISI axis
#     (SC on the left, SI on the right — matches the ISI density figure) ---
#     Category labels match the ISI figure band labels exactly: SC / Partial-SI / SI.
SCENARIOS = [
    {"n_func": 0, "loci": [0, 0, 0, 0], "category": "SC",
     "cat_color": COL_SC_LABEL,  "sublabel": "4 broken copies"},
    {"n_func": 1, "loci": [0, 0, 0, 1], "category": "Partial-SI",
     "cat_color": COL_PSI_LABEL, "sublabel": "1 functional  |  3 broken"},
    {"n_func": 2, "loci": [0, 0, 1, 1], "category": "Partial-SI",
     "cat_color": COL_PSI_LABEL, "sublabel": "2 functional  |  2 broken"},
    {"n_func": 3, "loci": [0, 1, 1, 1], "category": "Partial-SI",
     "cat_color": COL_PSI_LABEL, "sublabel": "3 functional  |  1 broken"},
    {"n_func": 4, "loci": [1, 1, 1, 1], "category": "SI",
     "cat_color": COL_SI_LABEL,  "sublabel": "4 functional copies"},
]

CHROM_POSITIONS = [(1.15, "A"), (1.85, "A"), (3.15, "B"), (3.85, "B")]

def draw_chromosome(ax, cx, cy, height, width, body_color, locus_functional):
    """One chromosome: rounded grey-ish bar with a coloured S-locus circle."""
    ax.add_patch(FancyBboxPatch(
        (cx - width / 2, cy - height / 2), width, height,
        boxstyle="round,pad=0.01,rounding_size=0.14",
        facecolor=body_color, edgecolor=COL_CHROM_EDGE, linewidth=1.6,
        zorder=2,
    ))
    locus_color = COL_FUNCTIONAL if locus_functional else COL_BROKEN
    ax.add_patch(Circle(
        (cx, cy + height * 0.18), width * 0.55,
        facecolor=locus_color, edgecolor=COL_LOCUS_EDGE, linewidth=1.6,
        zorder=3,
    ))

def build_panel(ax, scenario):
    ax.set_xlim(0, 5); ax.set_ylim(0, 6)
    ax.set_aspect("equal"); ax.axis("off")
    # Panel header — big number of functional copies
    ax.text(2.5, 5.60, f"{scenario['n_func']}", ha="center", va="center",
            fontsize=32, fontweight="bold", color="#2c3e50")
    ax.text(2.5, 5.00, "functional copies", ha="center", va="center",
            fontsize=10, fontstyle="italic", color="#555555")
    # Chromosomes
    for (x, sub), func_flag in zip(CHROM_POSITIONS, scenario["loci"]):
        body = COL_CHROM_A if sub == "A" else COL_CHROM_B
        draw_chromosome(ax, cx=x, cy=3.30, height=1.90, width=0.55,
                        body_color=body, locus_functional=bool(func_flag))
    # Subgenome group labels underneath
    ax.text(1.50, 2.05, "Subgenome A", ha="center", va="top",
            fontsize=10, fontweight="bold", color="#2c3e50")
    ax.text(3.50, 2.05, "Subgenome B", ha="center", va="top",
            fontsize=10, fontweight="bold", color="#2c3e50")
    # Category label (large, ISI-axis-coloured) + sublabel (italic count)
    ax.text(2.5, 1.05, scenario["category"], ha="center", va="center",
            fontsize=22, fontweight="bold", color=scenario["cat_color"])
    ax.text(2.5, 0.35, scenario["sublabel"], ha="center", va="center",
            fontsize=10, fontstyle="italic", color="#555555")

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--figures-dir", type=Path, default=Path(DEFAULT_FIGURES))
    args = parser.parse_args()
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 5, figsize=(16, 6.5))
    for ax, scenario in zip(axes, SCENARIOS):
        build_panel(ax, scenario)

    fig.suptitle(
        "SRK functional copy number  →  breeding-system classification",
        fontsize=14, y=0.98,
    )
    # Caption note — two colour scales in one figure
    fig.text(0.5, 0.03,
             "S-locus colour  =  molecular function  (blue = functional SRK; "
             "grey = non-functional).      "
             "Category label colour  =  ISI-axis band  (red = SI; amber = Partial-SI; "
             "green = SC).",
             ha="center", fontsize=10, fontstyle="italic", color="#555555")
    fig.tight_layout(rect=[0, 0.06, 1, 0.94])

    out_pdf = args.figures_dir / "SRK_functional_copies_schematic.pdf"
    out_png = args.figures_dir / "SRK_functional_copies_schematic.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    print(f"[chromo] {out_pdf}")
    print(f"[chromo] {out_png}")

if __name__ == "__main__":
    main()
