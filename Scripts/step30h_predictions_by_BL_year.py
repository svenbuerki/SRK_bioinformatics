"""Step 30h / Phase IV — predictions split by year, grouped by BL.

Phase 5 § A.4.5. The core "re-run predictions split by year" view
the user asked for: each population's predicted SRK diversity and
pollen compatibility shown as parallel 2025 vs 2026 bars with 95 %
CIs, with populations grouped by new BL (BL1 → BL5, each a
sub-panel). No simulation re-run — just reads the step30g per-
(populationID, year) prediction TSV that already has new population
IDs and BL after Phase III.

Outputs
-------
figures/Phase5/step30h_pred_diversity_by_BL_year.png / .pdf
figures/Phase5/step30h_pred_pcompat_by_BL_year.png / .pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from srk_bl_constants import make_location_label

TABLES  = Path("Tables/Phase5")
FIGURES = Path("figures/Phase5")

BAR_2025 = "#0072B2"   # Okabe-Ito blue
BAR_2026 = "#E69F00"   # Okabe-Ito orange


def population_labels() -> dict[int, str]:
    """Build {new populationID → '{EO_locID}, ...'} by walking the
    (now remapped) crosswalk."""
    cw = pd.read_csv(TABLES / "step29a_population_crosswalk.tsv",
                      sep="\t", encoding="utf-8-sig")
    pairs = (cw[["populationID", "locationCode", "locationID"]]
                .drop_duplicates()
                .sort_values(["populationID", "locationID"]))
    out: dict[int, str] = {}
    for pop_id, g in pairs.groupby("populationID"):
        parts = []
        for _, r in g.iterrows():
            code = str(r["locationCode"]).strip()
            lid = int(r["locationID"])
            if code and code != "nan":
                parts.append(make_location_label(code, lid))
            else:
                parts.append(f"locID_{lid}")
        out[int(pop_id)] = ", ".join(sorted(set(parts)))
    return out


def plot_metric(pred: pd.DataFrame, labels: dict[int, str],
                  metric: str,
                  display_name: str,
                  x_axis_label: str,
                  out_png: Path, out_pdf: Path) -> None:
    """One sub-panel per BL; within each panel one row per population
    sorted by populationID; two bars per population (2025 blue + 2026
    orange) with 95 % CI error bars."""
    bls = sorted(pred["BL"].dropna().unique())
    heights = [max(int((pred["BL"] == b).sum() / 2), 1)   # /2 because 2 years
               for b in bls]
    x_max = max(pred[f"pred_{metric}_hi95"].max(),
                 pred[f"pred_{metric}_mean"].max())

    fig, axes = plt.subplots(
        len(bls), 1,
        figsize=(10.5, max(9, 0.35 * (len(pred) // 2) + 3)),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )
    if len(bls) == 1:
        axes = [axes]

    width = 0.4
    for ax, bl in zip(axes, bls):
        sub = (pred[pred["BL"] == bl]
                  .sort_values(["populationID", "year"])
                  .reset_index(drop=True))
        pops = sub["populationID"].unique()
        y_pos = np.arange(len(pops))
        for yr, color, offset in [(2025, BAR_2025, -width/2),
                                     (2026, BAR_2026, +width/2)]:
            yr_sub = sub[sub["year"] == yr].set_index("populationID")
            means = [yr_sub.loc[p, f"pred_{metric}_mean"]
                        if p in yr_sub.index else np.nan
                     for p in pops]
            los   = [yr_sub.loc[p, f"pred_{metric}_lo95"]
                        if p in yr_sub.index else np.nan
                     for p in pops]
            his   = [yr_sub.loc[p, f"pred_{metric}_hi95"]
                        if p in yr_sub.index else np.nan
                     for p in pops]
            means_arr = np.array(means, dtype=float)
            errs_lo = means_arr - np.array(los, dtype=float)
            errs_hi = np.array(his, dtype=float) - means_arr
            errs = np.vstack([np.clip(np.nan_to_num(errs_lo, nan=0),
                                        0, None),
                               np.clip(np.nan_to_num(errs_hi, nan=0),
                                        0, None)])
            ax.barh(y_pos + offset, means_arr, height=width,
                    color=color, alpha=0.85, label=str(yr),
                    xerr=errs, ecolor="#333",
                    error_kw={"elinewidth": 0.7, "capsize": 2})

        row_labels = [f"P{int(p):>2}  ({labels.get(int(p), '')})"
                       for p in pops]
        ax.set_yticks(y_pos)
        ax.set_yticklabels(row_labels, fontsize=8)
        ax.set_xlim(0, x_max * 1.05)
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(1.012, 0.5, f"{bl}\n(n = {len(pops)})",
                transform=ax.transAxes,
                color="#5e3c99", fontsize=11, fontweight="bold",
                va="center", ha="left")

    axes[-1].set_xlabel(x_axis_label)
    axes[0].legend(loc="lower right", fontsize=10, frameon=True)
    fig.suptitle(
        f"Phase 5 predicted {display_name} per population — "
        "2025 vs 2026, grouped by new BL",
        fontsize=12, y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 0.92, 0.985])
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[step30h-IV] Wrote {out_png.name} + .pdf")


def main() -> None:
    pred = pd.read_csv(TABLES / "step30g_prediction_population_year.tsv",
                        sep="\t", encoding="utf-8-sig")
    labels = population_labels()

    plot_metric(
        pred, labels,
        metric="srk_diversity",
        display_name="SRK allele diversity (distinct Fgs)",
        x_axis_label="Predicted distinct SRK alleles",
        out_png=FIGURES / "step30h_pred_diversity_by_BL_year.png",
        out_pdf=FIGURES / "step30h_pred_diversity_by_BL_year.pdf",
    )
    plot_metric(
        pred, labels,
        metric="pcompat",
        display_name="pollen compatibility",
        x_axis_label="Predicted pollen compatibility",
        out_png=FIGURES / "step30h_pred_pcompat_by_BL_year.png",
        out_pdf=FIGURES / "step30h_pred_pcompat_by_BL_year.pdf",
    )


if __name__ == "__main__":
    main()
