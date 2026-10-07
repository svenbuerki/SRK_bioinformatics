"""Step 29a — per-population per-year deme-size figure.

Visualises the data in `Tables/Phase5/step29a_demes_per_population_year.tsv`
(258 rows, one per (population, year, deme)) as a per-BL faceted figure:
one row per population, two years overlaid as the open/filled circle
convention used across Phase 5 (see Figures 8 and 10 of the compact
doc / Figures 3 and 5 of the long doc).

Design
------
- One subplot per BL (BL1 → BL5, area DESC → connectivity DESC).
- One row per population (ordered by the mean of each population's
  per-year total `component_N_fertile`, ascending).
- Each row holds up to **two series**: 2025 (open circles) and 2026
  (filled circles); vertically offset by ± 0.14 so a row that has
  both years shows two parallel strips of dots.
- Each dot is a single deme at x = `component_N_fertile` (log scale).
- Reference lines: red at x = 1 (single-plant SI floor); grey at
  x = 8 (species-pool floor — 4 alleles × 8 plants = 32 copies);
  dashed at x = 32 (species ceiling).
- Right margin: `2025 Kd/N  ·  2026 Kd/N` using the compact doc's
  standard annotation convention.
- BL label in the top-left corner of each subplot (coloured box).

Output
------
    figures/Phase5/step29a_demes_per_population_year.png + .pdf

The `populationID → locationCode(s)` lookup is in
`step30g_populations_classified.tsv` (compact doc § Populations has
the inline table).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from srk_bl_constants import BL_COLORS

TABLES  = Path("Tables/Phase5")
FIGURES = Path("figures/Phase5")

NEW_BL_ORDER = ["BL1", "BL2", "BL3", "BL4", "BL5"]
Y_OFFSET     = 0.14   # vertical half-offset between year strips

REF_SI_FLOOR          = 1
REF_SPECIES_POOL_8    = 8
REF_SPECIES_CEILING   = 32


def main() -> None:
    df = pd.read_csv(TABLES / "step29a_demes_per_population_year.tsv",
                     sep="\t", encoding="utf-8-sig")

    # Per-population totals across years → used only for row ordering
    totals_per_pop_year = (df.groupby(["populationID", "BL", "year"])
                              ["component_N_fertile"].sum().reset_index())
    pop_order = (totals_per_pop_year.groupby(["populationID", "BL"])
                   ["component_N_fertile"].mean()
                   .reset_index()
                   .rename(columns={"component_N_fertile": "mean_total"}))

    # Counts per (population, year) for the right-margin annotation
    deme_counts = (df.groupby(["populationID", "year"])
                      .agg(n_demes=("demeID", "nunique"),
                           N_fert=("component_N_fertile", "sum"))
                      .reset_index())

    bls = [b for b in NEW_BL_ORDER if (pop_order["BL"] == b).any()]
    heights = [max(int((pop_order["BL"] == b).sum()), 1) for b in bls]

    fig, axes = plt.subplots(
        len(bls), 1,
        figsize=(11.5, max(6.5, 0.30 * sum(heights) + 1.8)),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )
    if len(bls) == 1:
        axes = [axes]

    for ax, bl in zip(axes, bls):
        sub = (pop_order[pop_order["BL"] == bl]
                 .sort_values("mean_total", ascending=True)
                 .reset_index(drop=True))
        y = np.arange(len(sub))
        colour = BL_COLORS.get(bl, "#777777")

        # Reference lines
        ax.axvline(REF_SI_FLOOR, color="#b2182b", ls=":", lw=0.9,
                   alpha=0.6)
        ax.axvline(REF_SPECIES_POOL_8, color="#888", ls=":", lw=0.9,
                   alpha=0.7)
        ax.axvline(REF_SPECIES_CEILING, color="#444", ls="--", lw=0.9,
                   alpha=0.65)

        for i, row in sub.iterrows():
            pid = int(row["populationID"])
            yi = y[i]
            for yr in (2025, 2026):
                g = df[(df["populationID"] == pid) & (df["year"] == yr)]
                if g.empty:
                    continue
                sizes = g["component_N_fertile"].to_numpy().astype(float)
                sizes = np.clip(sizes, 1, None)
                yoff  = -Y_OFFSET if yr == 2025 else +Y_OFFSET
                mfc   = ("white" if yr == 2025 else colour)
                mew   = (1.4 if yr == 2025 else 0.5)
                ax.scatter(sizes,
                           np.full_like(sizes, yi + yoff, dtype=float),
                           s=42, facecolor=mfc, edgecolor=colour,
                           linewidth=mew, zorder=3)

        # Left margin: short P{N}
        ax.set_yticks(y)
        ax.set_yticklabels([f"P{int(r['populationID']):>2}"
                             for _, r in sub.iterrows()],
                            fontsize=9)

        # Right margin via twinx (per-year K-demes / N adults)
        ax2 = ax.twinx()
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks(y)
        right = []
        for _, r in sub.iterrows():
            pid = int(r["populationID"])
            parts = []
            for yr in (2025, 2026):
                row2 = deme_counts[(deme_counts["populationID"] == pid)
                                     & (deme_counts["year"] == yr)]
                if not row2.empty:
                    k = int(row2["n_demes"].iat[0])
                    n = int(row2["N_fert"].iat[0])
                    parts.append(f"{yr} {k}d/{n}")
            right.append("  ·  ".join(parts))
        ax2.set_yticklabels(right, fontsize=7.5, color="#444")
        ax2.tick_params(axis="y", length=0, pad=2)
        for s in ("top", "right", "left"):
            ax2.spines[s].set_visible(False)

        ax.set_xscale("log")
        ax.set_xlim(0.8, max(df["component_N_fertile"].max() * 1.3, 50))
        ax.set_ylim(-0.7, len(sub) - 0.3)
        ax.invert_yaxis()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="x", which="major", alpha=0.3)

        # BL label in top-left corner, coloured box
        ax.text(0.008, 0.97, bl,
                transform=ax.transAxes,
                fontsize=12, fontweight="bold", color=colour,
                va="top", ha="left",
                bbox=dict(facecolor="white", edgecolor=colour,
                           boxstyle="round,pad=0.25", alpha=0.9,
                           linewidth=1.0))

    axes[-1].set_xlabel("Deme size — component_N_fertile (log scale)",
                         fontsize=11)

    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="",
                markerfacecolor="white", markeredgecolor="#444",
                markeredgewidth=1.4, markersize=8, label="2025 deme (open)"),
        Line2D([0], [0], marker="o", linestyle="",
                markerfacecolor="#444", markeredgecolor="white",
                markeredgewidth=0.5, markersize=8, label="2026 deme (filled)"),
        Line2D([0], [0], color="#b2182b", ls=":", lw=0.9,
                label="N = 1 (single-plant SI floor)"),
        Line2D([0], [0], color="#888", ls=":", lw=0.9,
                label="N = 8 (species-pool floor, 4 × 8 = 32 copies)"),
        Line2D([0], [0], color="#444", ls="--", lw=0.9,
                label="N = 32 (species ceiling of distinct Fgs)"),
    ]
    fig.legend(handles=legend_handles,
                loc="upper center", bbox_to_anchor=(0.5, 0.965),
                fontsize=9, frameon=True, ncol=3)

    fig.suptitle(
        "Phase 5 — deme-size distribution per population per year  "
        "(2025 open · 2026 filled, log x-scale)",
        fontsize=13, y=0.998,
    )
    fig.tight_layout(rect=[0, 0, 0.90, 0.93])
    out_png = FIGURES / "step29a_demes_per_population_year.png"
    out_pdf = FIGURES / "step29a_demes_per_population_year.pdf"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[step29a-demes] Wrote {out_png.name} + .pdf")


if __name__ == "__main__":
    main()
