"""Step 29d — Mating pool structure per location.

Phase 5 defines a **50 m connected component** as one mating pool —
the set of adult plants that share a pollen pool at the 50 m
pollinator radius. A location may hold one mating pool (fully
connected) or many (fragmented), and each pool has its own
`component_N_fertile` adult count. This information is used by
step30's per-component simulations but has not had a dedicated
display figure.

This script fills that gap. It reads the event → 50 m component
lookup written by step29c and emits:

Figures/Phase5/step29d_mating_pool_structure.{png,pdf}
    Two-panel figure.
      Panel A — Dataset summary scatter: total adults (x, log) vs
        number of 50 m mating pools (y) per location, coloured by
        Bottleneck Lineage. Shows the fragmentation × size pattern
        across the dataset.
      Panel B — Per-location strip plot of per-mating-pool adult
        count. One row per location, grouped vertically by BL in
        canonical BL_ORDER. Each dot = one 50 m mating pool, placed
        at its `component_N_fertile` on the log x-axis. Reference
        lines at N = 1 (single-plant SI floor) and N = 8 (coupon-
        collector 32-tetraploid-allele floor) — same convention as
        Figure 2b.

Tables/Phase5/step29d_mating_pool_summary.tsv
    One row per location with the mating-pool structure (number of
    pools, largest pool, median pool size, etc.) summarised.

Dependency: run `step29c_fragmentation_aware_sampling.py` first — it
writes the `step29c_event_to_component_50m.tsv` lookup this script
reads.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from srk_bl_constants import (
    BL_COLORS, BL_ORDER, locationCode_to_bl, location_label_series,
)

DEFAULT_TABLES  = Path("Tables/Phase5")
DEFAULT_FIGURES = Path("figures/Phase5")
INPUT_LOOKUP    = DEFAULT_TABLES / "step29c_event_to_component_50m.tsv"


def build_summary(lookup: pd.DataFrame) -> pd.DataFrame:
    """Reduce the event-level lookup to one row per location, with
    mating-pool structure metrics."""
    # One row per (locationID, component_id_50m)
    pools = lookup.drop_duplicates(
        ["locationID", "component_id_50m"])[
        ["locationID", "locationCode", "component_id_50m",
         "component_N_fertile", "component_N_events"]
    ].copy()
    pools["component_N_fertile"] = pools["component_N_fertile"].astype(int)

    rows = []
    for (loc_id, loc_code), sub in pools.groupby(
            ["locationID", "locationCode"], sort=False):
        sizes = sub["component_N_fertile"].astype(int).values
        rows.append({
            "locationID":           int(loc_id),
            "locationCode":         loc_code,
            "n_mating_pools_50m":   int(len(sizes)),
            "total_adults":         int(sizes.sum()),
            "largest_pool_N":       int(sizes.max()),
            "smallest_pool_N":      int(sizes.min()),
            "median_pool_N":        float(np.median(sizes)),
            "largest_pool_share":   float(sizes.max() / sizes.sum()),
            "n_pools_below_SI_floor_1":        int((sizes <= 1).sum()),
            "n_pools_below_coupon_floor_8":    int((sizes < 8).sum()),
        })
    summary = pd.DataFrame(rows)
    summary["BL"] = locationCode_to_bl(summary["locationCode"]).values
    summary["BL"] = summary["BL"].fillna("Unassigned")

    # Project-wide label convention (2026-10-03): every per-location
    # row across every figure uses `{locationCode}_{locationID}` so
    # distinct physical sites that share a locationCode (EO8_27 /
    # EO8_28 / EO8_29) are never confused.
    summary["display_label"] = location_label_series(summary)

    return summary, pools


def _bl_order_for(summary: pd.DataFrame) -> list[str]:
    present = summary["BL"].unique().tolist()
    ordered = [b for b in BL_ORDER if b in present]
    if "Unassigned" in present:
        ordered.append("Unassigned")
    return ordered


def plot_mating_pool_structure(summary: pd.DataFrame,
                                pools: pd.DataFrame,
                                out_png: Path, out_pdf: Path) -> None:
    bls = _bl_order_for(summary)
    palette = {**BL_COLORS, "Unassigned": "#8a8a8a"}
    pools = pools.merge(summary[["locationID", "BL"]], on="locationID",
                         how="left")

    # Per-BL block heights proportional to the number of locations
    # in each BL.
    heights = [max(int((summary["BL"] == b).sum()), 1) for b in bls]
    fig = plt.figure(
        figsize=(11.0, max(6.5, 0.30 * sum(heights) + 1.8)),
    )
    gs = fig.add_gridspec(len(bls), 1, hspace=0.14,
                           height_ratios=heights)
    axesB = [fig.add_subplot(gs[i, 0]) for i in range(len(bls))]

    # -----------------------------------------------------------------
    # Per-location strip plot: component_N_fertile per pool
    # -----------------------------------------------------------------
    for ax_i, bl in zip(axesB, bls):
        sub = summary[summary["BL"] == bl].sort_values(
            "largest_pool_N", ascending=True).reset_index(drop=True)
        colour = palette[bl]
        y_labels = []
        for y_pos, (_, row) in enumerate(sub.iterrows()):
            loc_pools = pools[pools["locationID"] == row["locationID"]]
            sizes = loc_pools["component_N_fertile"].astype(int).values
            # jitter dots vertically within the row for readability
            # Keep all dots on the row centre line — no vertical jitter,
            # dot size encodes pool size instead.
            ax_i.scatter(
                sizes, np.full(len(sizes), y_pos),
                s=np.clip(np.sqrt(sizes) * 4, 25, 180),
                color=colour, edgecolor="white", linewidth=0.5,
                alpha=0.85, zorder=3,
            )
            y_labels.append(
                f"{row['display_label']}  "
                f"({row['n_mating_pools_50m']} pool"
                f"{'s' if row['n_mating_pools_50m'] != 1 else ''}, "
                f"{row['total_adults']} adults)"
            )
        ax_i.axvline(1, color="#b2182b", ls=":", lw=1.0, alpha=0.7)
        ax_i.axvline(8, color="#333333", ls="--", lw=0.9, alpha=0.5)
        ax_i.set_xscale("log", base=2)
        ax_i.set_yticks(np.arange(len(sub)))
        ax_i.set_yticklabels(y_labels, fontsize=8)
        ax_i.set_ylim(-0.6, len(sub) - 0.4)
        ax_i.text(1.01, 0.5, bl, transform=ax_i.transAxes,
                   fontsize=12, fontweight="bold", color=colour,
                   va="center", ha="left")
        ax_i.spines["top"].set_visible(False)
        ax_i.spines["right"].set_visible(False)

    xticks = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
    axesB[-1].set_xticks(xticks)
    axesB[-1].set_xticklabels([str(x) for x in xticks], fontsize=9)
    axesB[-1].set_xlabel(
        "Adult plants in each 50 m mating pool  (log₂ scale)  —  "
        "red dotted = N = 1 (single-plant SI floor); "
        "grey dashed = N = 8 (32-allele coupon-collector floor)",
        fontsize=10,
    )
    # Share x-limits
    all_sizes = pools["component_N_fertile"].astype(int)
    xmin = max(0.6, 0.9 * all_sizes.min())
    xmax = 1.15 * all_sizes.max()
    for ax_i in axesB:
        ax_i.set_xlim(xmin, xmax)
    # Only the last subpanel keeps x-tick labels; higher subpanels
    # share them visually.
    for ax_i in axesB[:-1]:
        ax_i.set_xticks(xticks)
        ax_i.set_xticklabels([])

    axesB[0].set_title(
        "Per-location mating-pool structure  (one dot per 50 m "
        "connected component)",
        fontsize=11, loc="left")

    fig.suptitle(
        "Mating-pool structure per LEPA location — 50 m connected "
        "components (= distinct mating pools)",
        fontsize=12, y=0.995,
    )
    fig.subplots_adjust(left=0.26, right=0.96, top=0.94, bottom=0.07)
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    if not INPUT_LOOKUP.exists():
        raise SystemExit(
            f"Missing {INPUT_LOOKUP} — run step29c first.")
    lookup = pd.read_csv(INPUT_LOOKUP, sep="\t", encoding="utf-8-sig")
    summary, pools = build_summary(lookup)

    DEFAULT_TABLES.mkdir(parents=True, exist_ok=True)
    DEFAULT_FIGURES.mkdir(parents=True, exist_ok=True)

    tsv_path = DEFAULT_TABLES / "step29d_mating_pool_summary.tsv"
    summary.to_csv(tsv_path, sep="\t", index=False)
    print(f"[step29d] Wrote {tsv_path} ({len(summary)} locations, "
          f"{len(pools)} mating pools total)")

    fig_png = DEFAULT_FIGURES / "step29d_mating_pool_structure.png"
    fig_pdf = DEFAULT_FIGURES / "step29d_mating_pool_structure.pdf"
    plot_mating_pool_structure(summary, pools, fig_png, fig_pdf)
    print(f"[step29d] Wrote {fig_png} + .pdf")

    # Short summary
    print()
    print("[step29d] ========== Dataset-wide mating-pool structure ==========")
    print(f"[step29d] Locations: {len(summary)}")
    print(f"[step29d] Total 50 m mating pools: {len(pools)}")
    print(f"[step29d] Mean pools per location: "
          f"{summary['n_mating_pools_50m'].mean():.1f}  "
          f"(median {int(summary['n_mating_pools_50m'].median())}, "
          f"max {summary['n_mating_pools_50m'].max()})")
    print(f"[step29d] Pool sizes: min {pools['component_N_fertile'].min()}, "
          f"median {int(pools['component_N_fertile'].median())}, "
          f"max {pools['component_N_fertile'].max()}")
    print(f"[step29d] Pools at or below SI floor (N=1): "
          f"{(pools['component_N_fertile'] <= 1).sum()} of {len(pools)}")
    print(f"[step29d] Pools below coupon-collector floor (N<8): "
          f"{(pools['component_N_fertile'] < 8).sum()} of {len(pools)}")


if __name__ == "__main__":
    main()
