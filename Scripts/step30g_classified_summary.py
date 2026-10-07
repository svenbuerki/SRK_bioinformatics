"""Step 30g companion — classified population summary.

Combines step29a_population_summary.tsv + step30g_across_year_comparison.tsv
into one master table that lists every population with its:

  - legacy locationCodes (comma-separated)
  - occupancy pattern (both_years / 2025_only / 2026_only)
  - trend classification across years:
      stable       — both CIs overlap on both metrics (both-year only)
      growth       — N_fert_2026 / N_fert_2025 >= 4 with positive Δdiv
      crash        — N_fert_2026 / N_fert_2025 <= 0.25 with negative Δdiv
      ambiguous    — both-year but doesn't meet stable/growth/crash
      single_year  — present only in 2025 or 2026 above-ground
  - N_fertile per year + total
  - n_events per year
  - deme counts per year
  - prediction summary (pred diversity + pred pcompat per year)

Outputs
-------
Tables/Phase5/step30g_populations_classified.tsv
figures/Phase5/step30g_populations_classified.png / .pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from srk_bl_constants import make_location_label

TABLES  = Path("Tables/Phase5")
FIGURES = Path("figures/Phase5")

CRASH_RATIO  = 0.25
GROWTH_RATIO = 4.0
STABLE_FLOOR = 20

CLASS_ORDER = ["crash", "stable", "growth", "ambiguous",
               "2025_only", "2026_only"]
CLASS_COLORS = {
    "crash":     "#D55E00",   # vermillion
    "stable":    "#5e3c99",   # purple
    "growth":    "#009E73",   # teal green
    "ambiguous": "#aaaaaa",   # grey
    "2025_only": "#0072B2",   # blue
    "2026_only": "#E69F00",   # orange
}
BAR_2025 = "#0072B2"
BAR_2026 = "#E69F00"


def classify(row: pd.Series) -> str:
    if row["occupancy"] == "2025_only":
        return "2025_only"
    if row["occupancy"] == "2026_only":
        return "2026_only"
    # both-year populations: use the across-year comparison outputs
    if pd.isna(row.get("stable_across_years")):
        return "ambiguous"
    if row["stable_across_years"] and row["n_fertile_total"] >= STABLE_FLOOR:
        return "stable"
    nfert_25 = row["n_fertile_2025"]
    nfert_26 = row["n_fertile_2026"]
    if nfert_25 > 0 and (nfert_26 / nfert_25) <= CRASH_RATIO \
            and row.get("delta_diversity", 0) < 0:
        return "crash"
    if nfert_25 > 0 and (nfert_26 / nfert_25) >= GROWTH_RATIO:
        return "growth"
    return "ambiguous"


def build_population_labels() -> pd.Series:
    """From the crosswalk, build a project-standard per-population
    label: comma-separated `{locationCode}_{locationID}` strings in
    sorted order, one row per unique (locationCode, locationID) pair."""
    cw = pd.read_csv(TABLES / "step29a_population_crosswalk.tsv",
                      sep="\t", encoding="utf-8-sig")
    cw["locationCode"] = cw["locationCode"].astype(str).replace("nan", "")
    pairs = cw[["populationID", "locationCode", "locationID"]].drop_duplicates()
    labels: dict[int, str] = {}
    for pop_id, sub in pairs.groupby("populationID"):
        parts = []
        for _, r in sub.iterrows():
            code = str(r["locationCode"]).strip()
            loc_id = int(r["locationID"])
            if code and code != "nan":
                parts.append(make_location_label(code, loc_id))
            else:
                # Historical locationID with no Phase 5 locationCode
                parts.append(f"locID_{loc_id}")
        parts = sorted(set(parts))
        labels[int(pop_id)] = ", ".join(parts)
    return pd.Series(labels, name="population_label")


def build_master() -> pd.DataFrame:
    pop  = pd.read_csv(TABLES / "step29a_population_summary.tsv",
                        sep="\t", encoding="utf-8-sig")
    comp = pd.read_csv(TABLES / "step30g_across_year_comparison.tsv",
                        sep="\t", encoding="utf-8-sig")
    pop = pop.merge(build_population_labels(),
                     left_on="populationID", right_index=True, how="left")

    # Rename comp n_fertile columns to match pop (avoid conflict)
    comp = comp.rename(columns={
        "N_fert_2025": "n_fertile_2025_comp",
        "N_fert_2026": "n_fertile_2026_comp",
        "N_fert_total": "n_fertile_total_comp",
    })
    # Keep only the across-year-specific columns from comp (drop
    # n_demes columns because they already exist in pop)
    keep_comp = [
        "populationID",
        "pred_diversity_2025", "pred_diversity_2026",
        "pred_pcompat_2025", "pred_pcompat_2026",
        "delta_diversity", "delta_pcompat",
        "diversity_CI_overlap", "pcompat_CI_overlap",
        "stable_across_years",
        "mean_nfert_per_event_2025", "mean_nfert_per_event_2026",
        "crash_ratio_nfert",
    ]
    comp = comp[[c for c in keep_comp if c in comp.columns]]
    df = pop.merge(comp, on="populationID", how="left")

    df["trend_class"] = df.apply(classify, axis=1)
    df["class_rank"] = df["trend_class"].map(
        {c: i for i, c in enumerate(CLASS_ORDER)}
    )
    # Within class: sort by BL (BL1 → BL5), then by populationID
    # (within-BL order already follows the dendrogram leaf position
    # from Phase III, so adjacent IDs within a BL are spatial
    # neighbours too). Rationale (team feedback 2026-10-06): grouping
    # by BL within each trend class exposes which BLs carry each
    # dynamic, and makes selection of pair candidates per BL direct.
    df["bl_rank"] = (df["BL"].str.extract(r"BL(\d+)")[0]
                      .astype("Int64")
                      .fillna(99))
    df = df.sort_values(["class_rank", "bl_rank", "populationID"],
                         ascending=[True, True, True]).reset_index(drop=True)

    # Reorder key columns to front
    front = [
        "populationID", "population_label", "trend_class",
        "locationCodes", "locationIDs", "occupancy",
        "n_events_2025", "n_events_2026",
        "n_fertile_2025", "n_fertile_2026", "n_fertile_total",
        "n_slickspots", "n_slickspots_both_years",
        "n_demes_2025", "n_demes_2026",
        "pred_diversity_2025", "pred_diversity_2026",
        "pred_pcompat_2025", "pred_pcompat_2026",
    ]
    rest = [c for c in df.columns if c not in front + ["class_rank", "bl_rank"]]
    df = df[front + rest]
    return df


def plot_summary(df: pd.DataFrame,
                  out_png: Path, out_pdf: Path) -> None:
    """Horizontal grouped-bar chart, faceted vertically by trend class.
    One subplot per class; within each subplot, populations are sorted
    by BL (BL1 → BL5) then populationID. BL grouping inside each
    class panel is shown as a **facet-strip** on the left margin —
    a thin coloured vertical stripe + bold BL label — not as colour
    on the row labels themselves (user preference, 2026-10-07)."""
    classes_present = [c for c in CLASS_ORDER
                       if (df["trend_class"] == c).any()]
    heights = [max(int((df["trend_class"] == c).sum()), 1)
                for c in classes_present]
    x_max = max(df["n_fertile_2025"].max(), df["n_fertile_2026"].max(),
                 1)

    fig, axes = plt.subplots(
        len(classes_present), 1,
        figsize=(11.5, max(8, 0.33 * len(df) + 2.5)),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )
    if len(classes_present) == 1:
        axes = [axes]

    try:
        from srk_bl_constants import BL_COLORS as _BL_COL
        bl_colour = {k: v for k, v in _BL_COL.items()}
    except Exception:
        bl_colour = {}

    # BL facet strip positions (in axes-fraction x, data y via
    # ax.get_yaxis_transform()). Strip sits in the left margin,
    # outside the plot area.
    STRIP_X     = -0.085   # x-centre of coloured stripe
    STRIP_W     = 0.016    # stripe width (axes-fraction)
    LABEL_X     = -0.100   # x-centre of BL text (further left)

    width = 0.4
    for ax, cls in zip(axes, classes_present):
        sub = df[df["trend_class"] == cls].reset_index(drop=True)
        y_pos = np.arange(len(sub))
        ax.barh(y_pos - width/2, sub["n_fertile_2025"], height=width,
                color=BAR_2025, alpha=0.85, label="2025")
        ax.barh(y_pos + width/2, sub["n_fertile_2026"], height=width,
                color=BAR_2026, alpha=0.85, label="2026")

        # Row labels: short 'P{N}' only. The full locationCode list
        # is in step30g_populations_classified.tsv — the figure
        # caption cites that lookup. Keeps the figure legible when
        # populations merge several legacy locationCodes.
        ax.set_yticks(y_pos)
        ax.set_yticklabels(
            [f"P{int(r['populationID']):>2}"
             for _, r in sub.iterrows()],
            fontsize=9,
        )

        # BL facet-strip rendering within this class panel
        bls_in_class = (sub["BL"].astype(str)
                           .where(sub["BL"].notna(), "")
                           .tolist())
        if bls_in_class:
            # Walk contiguous BL blocks (sub is already BL-sorted)
            i = 0
            while i < len(bls_in_class):
                bl = bls_in_class[i]
                j = i
                while j + 1 < len(bls_in_class) and bls_in_class[j + 1] == bl:
                    j += 1
                # Block = rows i..j inclusive
                y_top = i - 0.5
                y_bot = j + 0.5
                y_center = (y_top + y_bot) / 2
                col = bl_colour.get(bl, "#777777")
                # Coloured strip (axes-fraction x, data y)
                ax.add_patch(Rectangle(
                    (STRIP_X - STRIP_W / 2, y_top),
                    STRIP_W, y_bot - y_top,
                    transform=ax.get_yaxis_transform(),
                    facecolor=col, edgecolor="none", alpha=0.95,
                    clip_on=False, zorder=5,
                ))
                # Bold BL label to the left of the stripe
                if bl:
                    ax.text(LABEL_X, y_center, bl,
                            transform=ax.get_yaxis_transform(),
                            fontsize=9.5, fontweight="bold", color=col,
                            ha="center", va="center", rotation=90,
                            clip_on=False, zorder=6)
                # Horizontal separator between BL blocks (not above the
                # first block)
                if i > 0:
                    ax.axhline(i - 0.5, color="#888888", lw=0.7, alpha=0.7)
                i = j + 1

        ax.set_xlim(0, x_max * 1.05)
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Class title block on the right margin
        ax.text(1.012, 0.5, f"{cls}\n(n = {len(sub)})",
                transform=ax.transAxes,
                color=CLASS_COLORS[cls], fontsize=12, fontweight="bold",
                va="center", ha="left")

    axes[-1].set_xlabel("N_fertile per year")
    axes[0].legend(loc="lower right", fontsize=10, frameon=True)

    fig.suptitle(
        "Phase 5 populations — 2025 vs 2026 N_fertile by across-year "
        f"trend class (n = {len(df)} populations); BL facet strip on "
        f"left",
        fontsize=12, y=0.995,
    )
    fig.tight_layout(rect=[0.05, 0, 0.92, 0.985])
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[step30g-classified] Wrote {out_png.name} + .pdf")



def main() -> None:
    df = build_master()
    df_write = df.drop(columns=["class_rank", "bl_rank"], errors="ignore")
    df_write.to_csv(TABLES / "step30g_populations_classified.tsv",
                     sep="\t", index=False)
    print("[step30g-classified] Wrote Tables/Phase5/step30g_populations_classified.tsv")
    print()
    print("Class counts:")
    print(df["trend_class"].value_counts().reindex(CLASS_ORDER)
          .fillna(0).astype(int).to_string())
    print()
    print("Full breakdown (populationID, label, N_fert 2025/2026, class):")
    for _, r in df.iterrows():
        print(f"  P{int(r['populationID']):>2}  "
              f"{str(r['population_label']):<36}  "
              f"{int(r['n_fertile_2025']):>4}/{int(r['n_fertile_2026']):>4}  "
              f"{r['trend_class']}")

    plot_summary(
        df,
        FIGURES / "step30g_populations_classified.png",
        FIGURES / "step30g_populations_classified.pdf",
    )


if __name__ == "__main__":
    main()
