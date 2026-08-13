#!/usr/bin/env python3
"""Pie chart of per-individual SI status (Step 22b, robust subset), coloured to
match the ISI breeding-system figures. Bridges the phenotypic ISI story
(species-level) with the molecular per-individual genotyping (Step 22).

Source: Tables/Phase4/step22b_individual_SI_status.tsv (497 ingroup rows).
Robust subset excludes Insufficient_data (247 SI + 15 pSI + 1 SC = 263 rows).
The 3-tier pSI severity ladder (1/2/3 non-functional copies) is merged into a
single pSI category — the point of the pie is to show the ISI-aligned three
breeding-system categories, not the pSI severity structure.

Colour mapping — matches the ISI breeding-system figures, NOT the existing
Step 22b figures (which use the opposite red/green semantics for
"machinery-intact = good" framing). The colour swap is deliberate for the
talk: it lets the pie sit visually alongside the ISI density figure without
a colour-legend clash.

    SC  → #1b7837  green   (matches ISI SC band, ISI < 0.2)
    pSI → #e08214  amber   (matches ISI Partial-SI band, 0.2 ≤ ISI < 0.8)
    SI  → #b2182b  red     (matches ISI SI band, ISI ≥ 0.8)

Output:
    figures/SRK_genotypic_SI_status_pie.{pdf,png}
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

DEFAULT_INPUT = "Tables/Phase4/step22b_individual_SI_status.tsv"
DEFAULT_FIGURES = "figures"

# Match ISI breeding-system figure palette (analyze_ISI_breeding_system.py)
COL_SC  = "#1b7837"
COL_PSI = "#e08214"
COL_SI  = "#b2182b"

CATEGORY_ORDER = ["SC", "pSI", "SI"]        # pie order: left-to-right on ISI axis
LEGEND_ORDER   = ["SI", "pSI", "SC"]        # legend order: dominant result first
CATEGORY_COLOURS = {"SC": COL_SC, "pSI": COL_PSI, "SI": COL_SI}
CATEGORY_LONG = {
    "SC":  "Self-compatible",
    "pSI": "Partial-SI",
    "SI":  "Self-incompatible",
}
PIE_ALPHA = 0.40                            # match soft feel of ISI band shading


def load_and_classify(path: Path) -> pd.DataFrame:
    """Read Step 22b TSV, drop Insufficient_data, merge pSI subcategories."""
    warnings.filterwarnings("ignore")
    df = pd.read_csv(path, sep="\t", encoding="utf-8-sig")
    df = df[df["SI_status"] != "Insufficient_data"].copy()
    df["Category"] = df["SI_status"].where(
        df["SI_status"].isin(["SI", "SC"]), "pSI"
    )
    return df


def plot_pie(df: pd.DataFrame, out_pdf: Path, out_png: Path) -> None:
    counts = {c: int((df["Category"] == c).sum()) for c in CATEGORY_ORDER}
    n_total = sum(counts.values())

    sizes = [counts[c] for c in CATEGORY_ORDER]
    colours = [CATEGORY_COLOURS[c] for c in CATEGORY_ORDER]
    # SC is 1 individual (0.4 %); pSI is 15 (5.7 %) — explode both a touch so
    # the small slices are visible against the dominant SI wedge
    explode = {"SC": 0.14, "pSI": 0.03, "SI": 0.0}
    explode_list = [explode[c] for c in CATEGORY_ORDER]

    fig, (ax_pie, ax_legend) = plt.subplots(
        1, 2, figsize=(11.0, 6.5),
        gridspec_kw={"width_ratios": [1.0, 1.15]},
    )
    wedges, _ = ax_pie.pie(
        sizes,
        colors=colours,
        labels=None,
        startangle=90,
        counterclock=False,
        explode=explode_list,
        wedgeprops=dict(edgecolor="white", linewidth=2.5, alpha=PIE_ALPHA),
    )
    ax_pie.set_aspect("equal")

    # Legend panel on the right: colour swatch + short + long + count + percent
    ax_legend.set_xlim(0, 1); ax_legend.set_ylim(0, 1)
    ax_legend.axis("off")
    ys = [0.72, 0.50, 0.28]
    swatch_x, text_x = 0.02, 0.14
    for y, cat in zip(ys, LEGEND_ORDER):
        colour = CATEGORY_COLOURS[cat]
        ax_legend.add_patch(plt.Rectangle(
            (swatch_x, y - 0.045), 0.08, 0.09,
            facecolor=colour, edgecolor="white", linewidth=1.5,
            alpha=PIE_ALPHA,
            transform=ax_legend.transAxes,
        ))
        pct = counts[cat] / n_total * 100
        ax_legend.text(
            text_x, y + 0.020, f"{cat}   —   {CATEGORY_LONG[cat]}",
            fontsize=15, fontweight="bold", color=colour,
            va="center", ha="left", transform=ax_legend.transAxes,
        )
        ax_legend.text(
            text_x, y - 0.035, f"n = {counts[cat]}   ({pct:.1f} %)",
            fontsize=13, color="#333333",
            va="center", ha="left", transform=ax_legend.transAxes,
        )

    fig.suptitle(
        f"Per-individual SI status  —  molecular genotyping  (robust n = {n_total})",
        fontsize=13, y=0.97,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", type=Path, default=Path(DEFAULT_INPUT))
    parser.add_argument("--figures-dir", type=Path, default=Path(DEFAULT_FIGURES))
    args = parser.parse_args()
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    df = load_and_classify(args.input)
    out_pdf = args.figures_dir / "SRK_genotypic_SI_status_pie.pdf"
    out_png = args.figures_dir / "SRK_genotypic_SI_status_pie.png"
    plot_pie(df, out_pdf, out_png)

    print(f"[pie] Robust n = {len(df)}")
    print(df["Category"].value_counts().reindex(CATEGORY_ORDER))
    print(f"[pie] Figure: {out_pdf}")
    print(f"[pie] Figure: {out_png}")


if __name__ == "__main__":
    main()
