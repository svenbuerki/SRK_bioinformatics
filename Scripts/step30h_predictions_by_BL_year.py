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


def plot_metric_both_years(pred: pd.DataFrame,
                             labels: dict[int, str],
                             metric: str,
                             x_label: str, suptitle: str,
                             x_upper: float, traffic_light: bool,
                             species_mean_x: float | None,
                             out_png: Path, out_pdf: Path) -> None:
    """Two-year overlay. One row per population; the 2025 draw is
    shown as an **open** circle, the 2026 draw as a **filled**
    circle in the same BL colour; a thin connector line joins the
    two when both years are present (so the 2025→2026 trajectory
    is read at a glance). Single-year populations show only the
    applicable year's marker, no connector."""

    # Union of populations present in either year, keeping their BL + legacy label.
    rep = (pred.sort_values(["populationID", "year"])
                .drop_duplicates("populationID", keep="first")
                [["populationID", "BL"]].copy())

    # Per-population year × metric lookup
    def _by(year: int) -> pd.DataFrame:
        return pred[pred["year"] == year].set_index("populationID")

    p25, p26 = _by(2025), _by(2026)

    # Row sort key within BL: mean of available-year metric means (ascending)
    def _row_key(pid: int) -> float:
        vals = []
        if pid in p25.index:
            vals.append(p25.loc[pid, f"pred_{metric}_mean"])
        if pid in p26.index:
            vals.append(p26.loc[pid, f"pred_{metric}_mean"])
        return float(np.mean(vals)) if vals else 0.0

    rep["sort_key"] = rep["populationID"].map(_row_key)

    bls = [b for b in NEW_BL_ORDER if b in rep["BL"].values]
    heights = [max(int((rep["BL"] == b).sum()), 1) for b in bls]

    fig, axes = plt.subplots(
        len(bls), 1,
        figsize=(11.5, max(6.5, 0.30 * sum(heights) + 1.8)),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )
    if len(bls) == 1:
        axes = [axes]

    for ax, bl in zip(axes, bls):
        sub = (rep[rep["BL"] == bl]
                   .sort_values("sort_key", ascending=True)
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
        # Species-mean / ceiling reference line
        if species_mean_x is not None:
            ax.axvline(species_mean_x, color="#1b7837", ls=":",
                       lw=1.0, alpha=0.75)

        # Per-row: fetch both years' values and plot
        for i, row in sub.iterrows():
            pid = int(row["populationID"])
            yi = y[i]
            has25 = pid in p25.index
            has26 = pid in p26.index

            m25 = p25.loc[pid, f"pred_{metric}_mean"] if has25 else None
            m26 = p26.loc[pid, f"pred_{metric}_mean"] if has26 else None

            # Connector line between the two years' means
            if has25 and has26:
                ax.plot([m25, m26], [yi, yi],
                         color=colour, lw=1.0, alpha=0.55, zorder=1)

            # 2025 draw — open circle
            if has25:
                lo = p25.loc[pid, f"pred_{metric}_lo95"]
                hi = p25.loc[pid, f"pred_{metric}_hi95"]
                xerr = [[max(m25 - lo, 0)], [max(hi - m25, 0)]]
                ax.errorbar([m25], [yi], xerr=xerr, fmt="none",
                             ecolor=colour, alpha=0.4, elinewidth=1.0,
                             capsize=2.0, zorder=1)
                n25 = p25.loc[pid, "N_fertile_total"]
                sz25 = 30 + 8 * np.sqrt(max(float(n25), 1.0))
                ax.scatter([m25], [yi], s=sz25,
                            facecolor="white", edgecolor=colour,
                            linewidth=1.5, zorder=2)
            # 2026 draw — filled circle
            if has26:
                lo = p26.loc[pid, f"pred_{metric}_lo95"]
                hi = p26.loc[pid, f"pred_{metric}_hi95"]
                xerr = [[max(m26 - lo, 0)], [max(hi - m26, 0)]]
                ax.errorbar([m26], [yi], xerr=xerr, fmt="none",
                             ecolor=colour, alpha=0.55, elinewidth=1.0,
                             capsize=2.0, zorder=1)
                n26 = p26.loc[pid, "N_fertile_total"]
                sz26 = 30 + 8 * np.sqrt(max(float(n26), 1.0))
                ax.scatter([m26], [yi], s=sz26, c=colour,
                            edgecolor="white", linewidth=0.6, zorder=3)

        # LEFT margin: short population IDs only. The locationCode
        # list is in the lookup TSV cited in the figure caption, so
        # the figure stays legible when population labels are long.
        ax.set_yticks(y)
        ax.set_yticklabels([f"P{int(r['populationID']):>2}"
                             for _, r in sub.iterrows()],
                            fontsize=9)
        ax.set_xlim(0.0, x_upper)
        ax.set_ylim(-0.7, len(sub) - 0.3)

        # RIGHT margin (secondary y-axis): per-year K-demes / N-adults
        # annotation, flush just outside the plot's right edge.
        ax2 = ax.twinx()
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks(y)
        right_lbls = []
        for _, r in sub.iterrows():
            pid = int(r["populationID"])
            parts = []
            if pid in p25.index:
                k25 = int(p25.loc[pid, "n_demes"])
                n25 = int(p25.loc[pid, "N_fertile_total"])
                parts.append(f"2025 {k25}d/{n25}")
            if pid in p26.index:
                k26 = int(p26.loc[pid, "n_demes"])
                n26 = int(p26.loc[pid, "N_fertile_total"])
                parts.append(f"2026 {k26}d/{n26}")
            right_lbls.append("  ·  ".join(parts))
        ax2.set_yticklabels(right_lbls, fontsize=7.5, color="#444")
        ax2.tick_params(axis="y", length=0, pad=2)
        for s in ("top", "right", "left"):
            ax2.spines[s].set_visible(False)

        # BL label — top-left corner of the subplot, in BL colour,
        # with a translucent BL-tinted box so it reads as a facet strip
        # without overlapping the right-margin per-year annotations.
        ax.text(0.008, 0.97, bl,
                transform=ax.transAxes,
                fontsize=12, fontweight="bold", color=colour,
                va="top", ha="left",
                bbox=dict(facecolor="white", edgecolor=colour,
                           boxstyle="round,pad=0.25", alpha=0.9,
                           linewidth=1.0))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Figure-level legend at the top, below the suptitle — frees up
    # the subplots' upper-left corners for the per-BL facet labels.
    legend_handles: list = [
        Line2D([0], [0], marker="o", linestyle="",
                markerfacecolor="white", markeredgecolor="#444",
                markeredgewidth=1.4, markersize=8, label="2025 (open)"),
        Line2D([0], [0], marker="o", linestyle="",
                markerfacecolor="#444", markeredgecolor="white",
                markeredgewidth=0.5, markersize=8, label="2026 (filled)"),
    ]
    if traffic_light:
        legend_handles += [
            Patch(facecolor="#b2182b", alpha=0.35,
                  label=f"failed  (< {T_FAILED_HI:.3f})"),
            Patch(facecolor="#e08214", alpha=0.35,
                  label=f"struggling"),
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
    fig.legend(handles=legend_handles,
                loc="upper center",
                bbox_to_anchor=(0.5, 0.965),
                fontsize=9, frameon=True,
                ncol=min(len(legend_handles), 4))

    axes[-1].set_xlabel(x_label, fontsize=11)
    fig.suptitle(suptitle, fontsize=13, y=0.998)
    fig.tight_layout(rect=[0, 0, 0.90, 0.93])
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[step30h-IV] Wrote {out_png.name} + .pdf")


def main() -> None:
    pred = pd.read_csv(TABLES / "step30g_prediction_population_year.tsv",
                        sep="\t", encoding="utf-8-sig")
    labels = population_labels()

    # SRK diversity — merged 2025+2026 overlay (user request 2026-10-07)
    plot_metric_both_years(
        pred, labels,
        metric="srk_diversity",
        x_label=("Predicted distinct SRK alleles at the population  "
                  "(mean; 95 % credible interval)"),
        suptitle=("Predicted SRK allele diversity — Snake River Plain "
                   "populations — 2025 (open) vs 2026 (filled)"),
        x_upper=K_FG_CEILING + 1,
        traffic_light=False,
        species_mean_x=K_FG_CEILING,
        out_png=FIGURES / "step30h_pred_diversity.png",
        out_pdf=FIGURES / "step30h_pred_diversity.pdf",
    )

    # Pollen compatibility — merged 2025+2026 overlay (user request 2026-10-07)
    plot_metric_both_years(
        pred, labels,
        metric="pcompat",
        x_label=("Predicted pollen compatibility under random mating  "
                  "(mean per population; 95 % credible interval)"),
        suptitle=("Predicted per-mother pollen compatibility — Snake "
                   "River Plain populations — 2025 (open) vs 2026 (filled)"),
        x_upper=X_UPPER_PC,
        traffic_light=True,
        species_mean_x=SPECIES_MEAN_PC,
        out_png=FIGURES / "step30h_pred_pcompat.png",
        out_pdf=FIGURES / "step30h_pred_pcompat.pdf",
    )


if __name__ == "__main__":
    main()
