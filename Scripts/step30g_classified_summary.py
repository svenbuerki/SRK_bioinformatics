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


def build_master() -> pd.DataFrame:
    pop  = pd.read_csv(TABLES / "step29a_population_summary.tsv",
                        sep="\t", encoding="utf-8-sig")
    comp = pd.read_csv(TABLES / "step30g_across_year_comparison.tsv",
                        sep="\t", encoding="utf-8-sig")

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
    # Within class: sort by N_fert_total DESC
    df = df.sort_values(["class_rank", "n_fertile_total"],
                         ascending=[True, False]).reset_index(drop=True)

    # Reorder key columns to front
    front = [
        "populationID", "trend_class", "locationCodes", "locationIDs",
        "occupancy",
        "n_events_2025", "n_events_2026",
        "n_fertile_2025", "n_fertile_2026", "n_fertile_total",
        "n_slickspots", "n_slickspots_both_years",
        "n_demes_2025", "n_demes_2026",
        "pred_diversity_2025", "pred_diversity_2026",
        "pred_pcompat_2025", "pred_pcompat_2026",
    ]
    rest = [c for c in df.columns if c not in front + ["class_rank"]]
    df = df[front + rest]
    return df


def plot_summary(df: pd.DataFrame,
                  out_png: Path, out_pdf: Path) -> None:
    """Horizontal grouped-bar chart, faceted vertically by trend class.
    One subplot per class, with populations ranked by N_fertile_total
    within class. Row labels: 'P{N}  (locationCodes)'. Each subplot
    gets a title with the class name + count."""
    classes_present = [c for c in CLASS_ORDER
                       if (df["trend_class"] == c).any()]
    heights = [max(int((df["trend_class"] == c).sum()), 1)
                for c in classes_present]
    x_max = max(df["n_fertile_2025"].max(), df["n_fertile_2026"].max(),
                 1)

    fig, axes = plt.subplots(
        len(classes_present), 1,
        figsize=(10.5, max(8, 0.33 * len(df) + 2.5)),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )
    if len(classes_present) == 1:
        axes = [axes]

    def _lbl(row: pd.Series) -> str:
        codes = row["locationCodes"]
        if pd.isna(codes) or str(codes) == "nan" or codes == "":
            codes = f"locID {row['locationIDs']}"
        return f"P{int(row['populationID']):>2}  ({codes})"

    width = 0.4
    for ax, cls in zip(axes, classes_present):
        sub = df[df["trend_class"] == cls].reset_index(drop=True)
        y_pos = np.arange(len(sub))
        ax.barh(y_pos - width/2, sub["n_fertile_2025"], height=width,
                color=BAR_2025, alpha=0.85, label="2025")
        ax.barh(y_pos + width/2, sub["n_fertile_2026"], height=width,
                color=BAR_2026, alpha=0.85, label="2026")

        ax.set_yticks(y_pos)
        ax.set_yticklabels([_lbl(r) for _, r in sub.iterrows()],
                             fontsize=8.5)
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
        f"trend class (n = {len(df)} populations)",
        fontsize=12, y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 0.92, 0.985])
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[step30g-classified] Wrote {out_png.name} + .pdf")


def main() -> None:
    df = build_master()
    df_write = df.drop(columns=["class_rank"], errors="ignore")
    df_write.to_csv(TABLES / "step30g_populations_classified.tsv",
                     sep="\t", index=False)
    print("[step30g-classified] Wrote Tables/Phase5/step30g_populations_classified.tsv")
    print()
    print("Class counts:")
    print(df["trend_class"].value_counts().reindex(CLASS_ORDER)
          .fillna(0).astype(int).to_string())
    print()
    print("Full breakdown (populationID, locationCodes, N_fert 2025/2026, class):")
    for _, r in df.iterrows():
        codes = r["locationCodes"] if pd.notna(r["locationCodes"]) \
                                   else f"locID {r['locationIDs']}"
        print(f"  P{int(r['populationID']):>2}  "
              f"{str(codes):<24}  "
              f"{int(r['n_fertile_2025']):>4}/{int(r['n_fertile_2026']):>4}  "
              f"{r['trend_class']}")

    plot_summary(
        df,
        FIGURES / "step30g_populations_classified.png",
        FIGURES / "step30g_populations_classified.pdf",
    )


if __name__ == "__main__":
    main()
