#!/usr/bin/env python3
"""Step 30b — Fragmentation of pollen flow (Phase A prediction).

A **purely spatial** metric of how much pollen flow at each LEPA site
is limited by pollinator range. It does not use allele frequencies or
drift; it only asks *how many pollen-donor plants can a mother reach
within 25 m?* — the biological definition of fragmentation for a
selfing-incompatible plant with short-range pollinators.

Two scales are reported side by side so the effect can be traced from
the individual mother up to the location:

Event scale (per event)
-----------------------
    K_spatial_25m = 2 · (sum of N_fertile in events within 25 m − 1)

The number of pollen-donor SRK allele copies a mother at that event
can access under the 25 m primary radius. Small K_spatial_25m = few
donors reachable = fragmentation-limited event.

We report, per event, the raw K_spatial_25m, the "any-neighbour" flag
(has ≥ 1 other event within 25 m), and the fragmentation index at
event scale
    F_event = 1 − K_spatial_25m / K_SPECIES_FG
so 0 = pollinator range gives the mother the full 32-allele species
pool, 1 = she has no reachable donor at all.

Location scale (per location)
-----------------------------
    connected_share_25m         (from step29_location_connectivity.tsv)
    largest_component_share_25m (from same)

At the location level, fragmentation is the fraction of adults NOT in
a multi-event pollen-flow component. We report
    F_location = 1 − connected_share_25m
plus the median event-scale F within the location, so the two scales
are visible together.

Both scales are inputs to Phase C's mate-limitation regression: the
per-mother K_spatial_25m becomes the β₂ (fragmentation) predictor; the
per-location F_location becomes an explicit location-level covariate.

Outputs
-------
    tables/Phase5/step30_A_fragmentation_per_event.tsv
    tables/Phase5/step30_A_fragmentation_per_location.tsv
    figures/Phase5/step30_A_fragmentation_index.png/pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from srk_bl_constants import BL_COLORS, BL_ORDER, locationCode_to_bl
from step28_seed_sampling_per_mother import (
    PLOIDY, K_SPECIES_FG, K_pool,
)

DEFAULT_TABLES = Path("Tables/Phase5")
DEFAULT_FIGURES = Path("figures/Phase5")

CONN_TSV  = DEFAULT_TABLES / "step29_location_connectivity.tsv"
EVENT_TSV = DEFAULT_TABLES / "step28_events_spatial_neighborhood.tsv"

R_PRIMARY = 25
# Under tetraploid (PLOIDY = 4), each donor plant contributes 4 SRK allele
# copies. The reference point where a mother's 25 m neighbourhood delivers
# enough copies to potentially saturate the 32-Fg species pool is
# N_ceiling_plants such that PLOIDY * (N-1) >= K_SPECIES_FG, i.e.
# N >= K_SPECIES_FG / PLOIDY + 1 = 9. Rounded to 8 (the point at which
# PLOIDY * N = K_SPECIES_FG copies excluding the mother).
N_FLOOR_SINGLE = 1                                # single-plant SI floor
N_CEILING_SPECIES = K_SPECIES_FG // PLOIDY        # = 8 under tetraploid


def build_event_table() -> pd.DataFrame:
    ev = pd.read_csv(EVENT_TSV, sep="\t", encoding="utf-8-sig")
    conn = pd.read_csv(CONN_TSV, sep="\t", encoding="utf-8-sig")
    loc_meta = conn[["locationID", "locationCode"]].drop_duplicates()

    K_col = f"K_spatial_{R_PRIMARY}m"
    N_col = f"N_neighbor_events_{R_PRIMARY}m"
    N_reach_col = f"N_compatible_spatial_{R_PRIMARY}m"

    df = ev.merge(loc_meta, on="locationID", how="left")
    df["K_spatial_25m"] = df[K_col].astype(int)          # PLOIDY-aware
    # N_reachable_25m = donor plant count in the 25 m neighbourhood
    # (excluding the mother's own N_fertile-1). This is the axis the
    # A.7 figure now uses so readers see plants, not allele copies.
    df["N_reachable_25m"] = df[N_reach_col].astype(int)
    df["has_neighbour_25m"] = (df[N_col].astype(int) > 0).astype(int)
    df["F_event"] = 1.0 - np.minimum(df["K_spatial_25m"], K_SPECIES_FG) / K_SPECIES_FG
    df["BL"] = locationCode_to_bl(df["locationCode"]).values
    df["BL"] = df["BL"].fillna("Unassigned")
    return df[[
        "eventID", "locationID", "locationCode", "BL",
        "n_fertile", "N_reachable_25m", "K_spatial_25m",
        "has_neighbour_25m", "F_event",
    ]].sort_values(["locationCode", "F_event"], ascending=[True, False])


def build_location_table(event_df: pd.DataFrame) -> pd.DataFrame:
    conn = pd.read_csv(CONN_TSV, sep="\t", encoding="utf-8-sig")
    agg = (event_df.groupby("locationID")
           .agg(median_K_spatial_25m=("K_spatial_25m", "median"),
                min_K_spatial_25m=("K_spatial_25m", "min"),
                frac_isolated_events=("has_neighbour_25m",
                                       lambda s: float((s == 0).mean())),
                median_F_event=("F_event", "median"))
           .reset_index())
    df = conn[["locationID", "locationCode", "n_events", "total_n_fertile",
               "connected_share_25m", "largest_component_share_25m",
               "n_components_25m"]].merge(agg, on="locationID", how="left")
    df["F_location"] = 1.0 - df["connected_share_25m"]
    df["BL"] = locationCode_to_bl(df["locationCode"]).values
    df["BL"] = df["BL"].fillna("Unassigned")
    return df.sort_values(["F_location", "median_F_event"],
                          ascending=[False, False])


def plot_event_K_by_location(ev_df: pd.DataFrame, loc_df: pd.DataFrame,
                              out_png: Path, out_pdf: Path):
    """Per-location distribution of the event-scale pollen-donor plant
    count reachable within 25 m, panelled by BL.

    x-axis: N_reachable_25m = donor plants a mother at that event can
    physically reach within pollinator range (excludes the mother's own
    event's plants but includes the mother's neighbours' plants). Plotted
    as a plain plant count so readers cannot confuse it with the 32-Fg
    count of distinct species-wide SRK allele classes (a different
    quantity — allele identities vs plant copies).

    Reference lines:
      * Red dotted at N = 1 — the single-plant SI floor. Events at or
        below this line have no reachable pollen donor and no seed
        set is possible under strict self-incompatibility.
      * Grey dashed at N = N_CEILING_SPECIES (= 8 under tetraploid) —
        the coupon-collector floor for the 32-Fg species pool. A mother
        reaching this many donor plants receives PLOIDY · N = 32
        allele copies, which is the *minimum physical count* needed for
        the species pool to be reachable in principle. Whether drift
        preserved the diversity is A.5's question.

    x-axis is log2 to keep both the 1 → 8 range (below the species
    floor) and the 8 → 500 range (above it) legible.
    """
    ev = ev_df.copy()
    ev["N_plot"] = np.maximum(ev["N_reachable_25m"].astype(int), 1)

    bls = [b for b in BL_ORDER if b in loc_df["BL"].values]
    if (loc_df["BL"] == "Unassigned").any():
        bls.append("Unassigned")
    palette = {**BL_COLORS, "Unassigned": "#8a8a8a"}

    heights = [max(int((loc_df["BL"] == b).sum()), 1) for b in bls]
    fig, axes = plt.subplots(
        len(bls), 1,
        figsize=(10.5, max(6.0, 0.32 * sum(heights) + 1.5)),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )
    if len(bls) == 1:
        axes = [axes]

    for ax, bl in zip(axes, bls):
        sub_locs = (loc_df[loc_df["BL"] == bl]
                    .sort_values("median_K_spatial_25m", ascending=True)
                    .reset_index(drop=True))
        colour = palette[bl]
        for i, r in sub_locs.iterrows():
            vals = ev[ev["locationID"] == r["locationID"]]["N_plot"].values
            if len(vals) == 0:
                continue
            ax.boxplot(
                [vals], positions=[i], vert=False, widths=0.6,
                patch_artist=True,
                boxprops=dict(facecolor=colour, alpha=0.55, edgecolor=colour),
                medianprops=dict(color="white", linewidth=1.6),
                whiskerprops=dict(color=colour, linewidth=1.0),
                capprops=dict(color=colour, linewidth=1.0),
                flierprops=dict(marker="o", markersize=3,
                                markerfacecolor=colour, alpha=0.55,
                                markeredgecolor="none"),
            )
            if len(vals) >= 1:
                jitter = (np.random.default_rng(int(r["locationID"]))
                          .uniform(-0.20, 0.20, size=len(vals)))
                ax.scatter(vals, np.full(len(vals), i) + jitter,
                           s=12, color=colour, edgecolor="white",
                           linewidth=0.4, alpha=0.75, zorder=3)
        labels = [f"{r['locationCode']}   (events = {int(r['n_events'])}, "
                  f"adults = {int(r['total_n_fertile'])})"
                  for _, r in sub_locs.iterrows()]
        ax.set_yticks(range(len(sub_locs)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_ylim(-0.7, len(sub_locs) - 0.3)
        ax.axvline(N_FLOOR_SINGLE, color="#b2182b", ls=":", lw=1.2, alpha=0.8)
        ax.axvline(N_CEILING_SPECIES,
                   color="#333333", ls="--", lw=1.0, alpha=0.6)
        ax.set_xscale("log", base=2)
        ax.text(1.01, 0.5, bl, transform=ax.transAxes,
                fontsize=13, fontweight="bold", color=colour,
                va="center", ha="left")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[-1].set_xlabel(
        "Pollen-donor plants reachable within 25 m per event  "
        "(N_reachable_25m = Σ N_fertile in events within 25 m − 1)\n"
        f"Red dotted: N = {N_FLOOR_SINGLE} single-plant SI floor  ·  "
        f"Grey dashed: N = {N_CEILING_SPECIES} plants = {PLOIDY * N_CEILING_SPECIES} "
        f"tetraploid allele copies, the coupon-collector floor for the 32-Fg "
        f"species pool",
        fontsize=10,
    )
    xticks = [1, 2, 4, 8, 16, 32, 64, 128, 256]
    axes[-1].set_xticks(xticks)
    axes[-1].set_xticklabels([str(x) for x in xticks])
    fig.suptitle(
        "Event-scale pollen-donor plant count per LEPA location "
        "(purely spatial, 25 m pollinator range, LEPA tetraploid)",
        fontsize=13, y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 0.94, 0.97])
    fig.savefig(out_png, dpi=200); fig.savefig(out_pdf); plt.close(fig)


def main() -> None:
    for p in (CONN_TSV, EVENT_TSV):
        if not p.exists():
            raise SystemExit(f"Missing input {p}")
    ev = build_event_table()
    out_ev = DEFAULT_TABLES / "step30_A_fragmentation_per_event.tsv"
    ev.to_csv(out_ev, sep="\t", index=False)
    print(f"[step30b] Wrote {out_ev}  ({len(ev)} events)")

    loc = build_location_table(ev)
    out_loc = DEFAULT_TABLES / "step30_A_fragmentation_per_location.tsv"
    loc.to_csv(out_loc, sep="\t", index=False)
    print(f"[step30b] Wrote {out_loc}  ({len(loc)} locations)")

    print("[step30b] Fragmentation summary (25 m primary):")
    print(f"  Events with ≥ 1 neighbour within 25 m : "
          f"{int(ev['has_neighbour_25m'].sum())} / {len(ev)}")
    print(f"  Locations with F_location ≥ 0.5      : "
          f"{int((loc['F_location'] >= 0.5).sum())} / {len(loc)}")
    print(f"  Median event-scale F_event           : "
          f"{ev['F_event'].median():.3f}")
    print(f"  Median location-scale F_location     : "
          f"{loc['F_location'].median():.3f}")

    DEFAULT_FIGURES.mkdir(parents=True, exist_ok=True)
    plot_event_K_by_location(
        ev, loc,
        out_png=DEFAULT_FIGURES / "step30_A_fragmentation_index.png",
        out_pdf=DEFAULT_FIGURES / "step30_A_fragmentation_index.pdf",
    )
    print(f"[step30b] Figure in {DEFAULT_FIGURES}/")


if __name__ == "__main__":
    main()
