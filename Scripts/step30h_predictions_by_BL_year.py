"""Step 30h / Phase IV — predictions in the fecundation style, per year.

Phase 5 § A.4.5. Rebuilt to match `step30_A_prediction_fecundation.png`:
dots + 95 % CI error bars, traffic-light background bands (pollen
compatibility only), species-mean reference line, BL-coloured points,
BL label on the right in bold BL colour, one row per population,
panels stacked by BL (BL_ORDER from srk_bl_constants).

User's request: PREDICTIONS SPLIT BY YEAR → one figure per year per
metric, i.e. four PNGs total:

    step30h_pred_pcompat_2025.png / .pdf
    step30h_pred_pcompat_2026.png / .pdf
    step30h_pred_diversity_2025.png / .pdf
    step30h_pred_diversity_2026.png / .pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from srk_bl_constants import BL_COLORS, make_location_label

TABLES  = Path("Tables/Phase5")
FIGURES = Path("figures/Phase5")

# BL_ORDER for the NEW population-based BL framework = BL1 → BL5
# (area DESC → connectivity DESC rule, locked in step30h Phase III).
NEW_BL_ORDER = ["BL1", "BL2", "BL3", "BL4", "BL5"]

# Traffic-light thresholds — same as step30_A_prediction_fecundation
T_FAILED_HI     = 0.259
T_STRUGGLING_HI = 0.519
SPECIES_MEAN_PC = 0.778
X_UPPER_PC      = 1.0

K_FG_CEILING    = 32  # SRK diversity species-wide ceiling


def population_labels() -> dict[int, str]:
    """Build {populationID (new) → '{EO_locID}, ...'} by walking the
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


def plot_metric_year(pred_y: pd.DataFrame,
                       labels: dict[int, str],
                       year: int,
                       metric: str,                 # "pcompat" | "srk_diversity"
                       x_label: str,
                       suptitle: str,
                       x_upper: float,
                       traffic_light: bool,
                       species_mean_x: float | None,
                       out_png: Path, out_pdf: Path) -> None:

    bls = [b for b in NEW_BL_ORDER if b in pred_y["BL"].values]
    heights = [max(int((pred_y["BL"] == b).sum()), 1) for b in bls]

    fig, axes = plt.subplots(
        len(bls), 1,
        figsize=(10, max(6.5, 0.30 * sum(heights) + 1.8)),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )
    if len(bls) == 1:
        axes = [axes]

    for ax, bl in zip(axes, bls):
        sub = (pred_y[pred_y["BL"] == bl]
                   .sort_values(f"pred_{metric}_mean", ascending=True)
                   .reset_index(drop=True))
        y = np.arange(len(sub))
        colour = BL_COLORS.get(bl, "#777777")

        # Traffic-light bands (pollen compatibility only)
        if traffic_light:
            band_alpha = 0.12
            ax.axvspan(0.0, T_FAILED_HI, color="#b2182b",
                       alpha=band_alpha, zorder=0)
            ax.axvspan(T_FAILED_HI, T_STRUGGLING_HI, color="#e08214",
                       alpha=band_alpha, zorder=0)
            ax.axvspan(T_STRUGGLING_HI, x_upper, color="#1b7837",
                       alpha=band_alpha, zorder=0)
        # Species-mean reference line
        if species_mean_x is not None:
            ax.axvline(species_mean_x, color="#1b7837", ls=":",
                       lw=1.0, alpha=0.75)

        # Error bars + dots
        means = sub[f"pred_{metric}_mean"].to_numpy()
        xerr_lo = np.clip(means - sub[f"pred_{metric}_lo95"].to_numpy(),
                           0, None)
        xerr_hi = np.clip(sub[f"pred_{metric}_hi95"].to_numpy() - means,
                           0, None)
        ax.errorbar(means, y, xerr=[xerr_lo, xerr_hi],
                     fmt="none", ecolor=colour, alpha=0.5,
                     elinewidth=1.2, capsize=2.5, zorder=1)
        sizes = 30 + 8 * np.sqrt(np.clip(sub["N_fertile_total"], 1, None))
        ax.scatter(means, y, s=sizes, c=colour, edgecolor="white",
                   linewidth=0.6, zorder=2)

        # Row labels: 'P{N}  ({EO_locID})  (K demes, N adults)'
        row_lbls = []
        for _, r in sub.iterrows():
            pid = int(r["populationID"])
            loc_lbl = labels.get(pid, "")
            k = int(r["n_demes"])
            n = int(r["N_fertile_total"])
            row_lbls.append(
                f"P{pid:>2}  ({loc_lbl})  "
                f"({k} deme{'s' if k != 1 else ''}, {n} adults)"
            )
        ax.set_yticks(y)
        ax.set_yticklabels(row_lbls, fontsize=8)
        ax.set_xlim(0.0, x_upper)
        ax.set_ylim(-0.7, len(sub) - 0.3)

        # BL label on the right
        ax.text(1.01, 0.5, bl,
                transform=ax.transAxes,
                fontsize=13, fontweight="bold", color=colour,
                va="center", ha="left")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Legend on the first panel
    legend_handles: list = []
    if traffic_light:
        legend_handles += [
            Patch(facecolor="#b2182b", alpha=0.35,
                  label=f"failed  (< {T_FAILED_HI:.3f})"),
            Patch(facecolor="#e08214", alpha=0.35,
                  label=f"struggling  "
                        f"({T_FAILED_HI:.3f}–{T_STRUGGLING_HI:.3f})"),
            Patch(facecolor="#1b7837", alpha=0.35,
                  label=f"sustainable  (≥ {T_STRUGGLING_HI:.3f})"),
        ]
    if species_mean_x is not None:
        tag = ("sporophytic species mean"
               if traffic_light else "species-wide ceiling")
        legend_handles.append(
            Line2D([0], [0], color="#1b7837", ls=":", lw=1.2,
                   label=f"{tag}  ({species_mean_x:.3f})"),
        )
    if legend_handles:
        axes[0].legend(handles=legend_handles, loc="upper left",
                        fontsize=9, frameon=True)

    axes[-1].set_xlabel(x_label, fontsize=11)
    fig.suptitle(suptitle, fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0, 0.94, 0.97])
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"[step30h-IV] Wrote {out_png.name} + .pdf")


def main() -> None:
    pred = pd.read_csv(TABLES / "step30g_prediction_population_year.tsv",
                        sep="\t", encoding="utf-8-sig")
    labels = population_labels()

    for yr in (2025, 2026):
        pred_y = pred[pred["year"] == yr].copy()

        # Pollen compatibility
        plot_metric_year(
            pred_y, labels, year=yr,
            metric="pcompat",
            x_label=("Predicted pollen compatibility under random "
                      "mating  (mean per population; 95 % credible interval)"),
            suptitle=f"Predicted per-mother pollen compatibility — "
                      f"Snake River Plain populations — {yr}",
            x_upper=X_UPPER_PC,
            traffic_light=True,
            species_mean_x=SPECIES_MEAN_PC,
            out_png=FIGURES / f"step30h_pred_pcompat_{yr}.png",
            out_pdf=FIGURES / f"step30h_pred_pcompat_{yr}.pdf",
        )

        # SRK diversity
        plot_metric_year(
            pred_y, labels, year=yr,
            metric="srk_diversity",
            x_label=("Predicted distinct SRK alleles at the population  "
                      "(mean; 95 % credible interval)"),
            suptitle=f"Predicted SRK allele diversity — Snake River Plain "
                      f"populations — {yr}",
            x_upper=K_FG_CEILING + 1,
            traffic_light=False,
            species_mean_x=K_FG_CEILING,
            out_png=FIGURES / f"step30h_pred_diversity_{yr}.png",
            out_pdf=FIGURES / f"step30h_pred_diversity_{yr}.pdf",
        )


if __name__ == "__main__":
    main()
