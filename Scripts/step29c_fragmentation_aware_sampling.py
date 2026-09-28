#!/usr/bin/env python3
"""Step 29c — Fragmentation-aware sampling allocation.

Uses the 25 m primary pollen-flight radius to group each location's
events into connected components, then allocates mothers per
COMPONENT (not per whole location) using the same coupon-collector
90 %-probability-of-full-detection rule as § B.4.1. Within each
component, mothers are distributed proportional to N_fertile, subject
to a maternal-genotype floor of ≥ 1 mother per event.

The component-level total is summed to give the location-level
M_frag_aware, which is compared head-to-head with the existing
per-location recommendation in
step28_mothers_for_full_detection_by_location.tsv (§ B.4.1).

Why this is different
---------------------
The current allocation (Step 29 + B.4.1) uses the location's
*aggregate* N_fertile as a single pool. Two events sitting 300 m apart
at the same locationID are treated as one mating unit even though no
pollen crosses. The fragmentation-aware version scopes the coupon-
collector maths to each 25 m component, so:

  * A location that is one well-connected cluster → unchanged.
  * A location that is many isolated events → n_events mothers
    (private-allele floor was already forcing this).
  * A location that is a big cluster + a few isolated singletons
    (e.g. EO27-1, EO26-3) → the cluster only needs its coupon-
    collector M (~6), the isolated events each need 1, and the total
    can DROP substantially below the current recommendation while
    still delivering the same 90 %-see-every-allele guarantee.

Outputs
-------
    tables/Phase5/step29c_sampling_frag_aware_per_event.tsv
    tables/Phase5/step29c_sampling_comparison_per_location.tsv
    figures/Phase5/step29c_sampling_comparison.png/pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from step28_seed_sampling_per_mother import (
    haversine_meters,
    mothers_for_full_detection,
    n_for_miss_probability,
    K_SPECIES_FG,
    K_pool,
    PLOIDY,
)
from srk_bl_constants import BL_COLORS, BL_ORDER, locationCode_to_bl

DEFAULT_TABLES = Path("Tables/Phase5")
DEFAULT_FIGURES = Path("figures/Phase5")

EVENT_TSV = DEFAULT_TABLES / "step28_events_spatial_neighborhood.tsv"
CURRENT_M_TSV = (DEFAULT_TABLES /
                 "step28_mothers_for_full_detection_by_location.tsv")
R_PRIMARY = 25.0
TARGET_PROB = 0.90


def connected_components_25m(sub: pd.DataFrame) -> list[list[int]]:
    """BFS on the 25 m adjacency graph. Returns positional-index lists,
    one per component."""
    n = len(sub)
    if n == 0:
        return []
    lat = sub["lat"].values
    lon = sub["lon"].values
    d = haversine_meters(lat[:, None], lon[:, None],
                         lat[None, :], lon[None, :])
    adj = (d > 0) & (d <= R_PRIMARY)
    seen = np.zeros(n, bool)
    comps: list[list[int]] = []
    for start in range(n):
        if seen[start]:
            continue
        stack = [start]
        comp: list[int] = []
        while stack:
            v = stack.pop()
            if seen[v]:
                continue
            seen[v] = True
            comp.append(v)
            for u in np.where(adj[v] & ~seen)[0]:
                stack.append(int(u))
        comps.append(comp)
    return comps


def allocate_component(comp_subset: pd.DataFrame,
                        target_prob: float = TARGET_PROB
                        ) -> pd.DataFrame:
    """Compute component K, coupon-collector M, then distribute across
    events proportional to N_fertile with a floor of ≥ 1 per event and
    a cap at each event's N_fertile."""
    N_c = int(comp_subset["n_fertile"].sum())
    n_events_c = len(comp_subset)
    K_c_raw = K_pool(N_c)   # PLOIDY * (N_c - 1) under tetraploid
    K_c = min(K_c_raw, K_SPECIES_FG)
    M_paternal = (1 if K_c <= 1
                  else mothers_for_full_detection(K_c, target_prob))
    M_c_final = max(M_paternal, n_events_c)

    events_out = comp_subset.copy().reset_index(drop=True)
    events_out["M_frag"] = 1  # maternal-genotype floor
    remaining = M_c_final - n_events_c
    if remaining > 0:
        weights = events_out["n_fertile"].astype(float).values
        weights = weights / weights.sum() if weights.sum() > 0 else weights
        raw = weights * remaining
        base = np.floor(raw).astype(int)
        residuals = raw - base
        extra = int(remaining - base.sum())
        if extra > 0:
            top = np.argsort(-residuals)[:extra]
            base[top] += 1
        events_out["M_frag"] = events_out["M_frag"].values + base

    # Cap at each event's N_fertile and redistribute any overflow
    n_cap = events_out["n_fertile"].astype(int).values
    while True:
        over_mask = events_out["M_frag"].astype(int).values > n_cap
        if not over_mask.any():
            break
        excess = int((events_out["M_frag"].values[over_mask]
                      - n_cap[over_mask]).sum())
        events_out.loc[over_mask, "M_frag"] = n_cap[over_mask]
        while excess > 0:
            headroom = (n_cap - events_out["M_frag"].astype(int).values)
            if (headroom > 0).sum() == 0:
                break
            i = int(np.argmax(headroom))
            events_out.at[i, "M_frag"] += 1
            excess -= 1
        if excess > 0:
            break

    events_out["component_K"] = K_c
    events_out["component_M_target"] = int(events_out["M_frag"].sum())
    events_out["component_N_events"] = n_events_c
    events_out["component_N_fertile"] = N_c
    events_out["M_paternal_only"] = M_paternal
    return events_out


def build_per_event() -> pd.DataFrame:
    ev = pd.read_csv(EVENT_TSV, sep="\t", encoding="utf-8-sig")
    parts = []
    for loc_id, sub in ev.groupby("locationID"):
        sub = sub.reset_index(drop=True)
        comps = connected_components_25m(sub)
        pieces = []
        for cid, idxs in enumerate(comps):
            comp_sub = sub.iloc[idxs].copy()
            allocated = allocate_component(comp_sub)
            allocated["component_id_within_loc"] = cid
            pieces.append(allocated)
        parts.append(pd.concat(pieces, ignore_index=True))
    return pd.concat(parts, ignore_index=True)


def build_per_location(event_df: pd.DataFrame) -> pd.DataFrame:
    # Per-component summary first so per-location aggregates are honest
    comp = (event_df
            .drop_duplicates(["locationID", "component_id_within_loc"])
            .groupby("locationID")
            .agg(sum_M_paternal_only=("M_paternal_only", "sum"))
            .reset_index())
    per_loc = (event_df.groupby("locationID")
               .agg(n_events=("eventID", "count"),
                    total_N_fertile=("n_fertile", "sum"),
                    n_components_25m=("component_id_within_loc",
                                       lambda s: int(s.nunique())),
                    M_frag_aware=("M_frag", "sum"))
               .reset_index()
               .merge(comp, on="locationID"))
    cur = pd.read_csv(CURRENT_M_TSV, sep="\t", encoding="utf-8-sig")
    cur_slim = cur[["locationID", "locationCode", "K_local",
                    "M_uniform_full_detection", "M_event_coverage",
                    "M_recommended"]].rename(columns={
        "M_recommended":           "M_current",
        "M_uniform_full_detection": "M_current_uniform",
        "M_event_coverage":         "M_current_event_floor",
    })
    per_loc = per_loc.merge(cur_slim, on="locationID", how="left")
    per_loc["delta"] = per_loc["M_frag_aware"] - per_loc["M_current"]
    per_loc["BL"] = locationCode_to_bl(per_loc["locationCode"]).values
    per_loc["BL"] = per_loc["BL"].fillna("Unassigned")
    return per_loc.sort_values(["BL", "delta"], ascending=[True, True])


def plot_comparison(loc_df: pd.DataFrame, out_png: Path, out_pdf: Path):
    bls = [b for b in BL_ORDER if b in loc_df["BL"].values]
    if (loc_df["BL"] == "Unassigned").any():
        bls.append("Unassigned")
    palette = {**BL_COLORS, "Unassigned": "#8a8a8a"}
    heights = [max(int((loc_df["BL"] == b).sum()), 1) for b in bls]
    fig, axes = plt.subplots(
        len(bls), 1,
        figsize=(11.0, max(6.0, 0.36 * sum(heights) + 2.0)),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )
    if len(bls) == 1:
        axes = [axes]
    for ax, bl in zip(axes, bls):
        sub = (loc_df[loc_df["BL"] == bl]
               .sort_values("delta", ascending=True)
               .reset_index(drop=True))
        y = np.arange(len(sub))
        colour = palette[bl]
        h = 0.36
        ax.barh(y - h/2, sub["M_current"],   h, color=colour, alpha=0.35,
                edgecolor=colour, linewidth=1.0,
                label="M current (§ B.4.1)")
        ax.barh(y + h/2, sub["M_frag_aware"], h, color=colour, alpha=0.9,
                edgecolor=colour, linewidth=1.0,
                label="M fragmentation-aware (§ B.4.2)")
        labels = [f"{r['locationCode']}  "
                  f"(events = {int(r['n_events'])}, "
                  f"components = {int(r['n_components_25m'])}, "
                  f"adults = {int(r['total_N_fertile'])})"
                  for _, r in sub.iterrows()]
        for i, r in sub.iterrows():
            delta = int(r["delta"])
            sign = "+" if delta > 0 else ""
            ax.text(max(r["M_current"], r["M_frag_aware"]) + 0.7, i,
                    f"Δ = {sign}{delta}", fontsize=8,
                    color=colour, va="center", ha="left")
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_ylim(-0.7, len(sub) - 0.3)
        ax.text(1.01, 0.5, bl, transform=ax.transAxes,
                fontsize=13, fontweight="bold", color=colour,
                va="center", ha="left")
        if bl == bls[0]:
            ax.legend(loc="lower right", fontsize=9, frameon=True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    xmax = float(max(loc_df["M_current"].max(),
                     loc_df["M_frag_aware"].max()))
    axes[-1].set_xlim(0, xmax + 12)
    axes[-1].set_xlabel(
        f"Number of mothers recommended per location "
        f"({n_for_miss_probability()} seeds each; tetraploid Rule 2)",
        fontsize=11)
    fig.suptitle(
        "Sampling recommendation — current (§ B.4.1) vs fragmentation-aware "
        "(§ B.4.2)",
        fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0, 0.94, 0.97])
    fig.savefig(out_png, dpi=200); fig.savefig(out_pdf); plt.close(fig)


def main() -> None:
    for p in (EVENT_TSV, CURRENT_M_TSV):
        if not p.exists():
            raise SystemExit(f"Missing input {p}")

    ev = build_per_event()
    ev_out = ev[[
        "locationID", "eventID", "component_id_within_loc",
        "n_fertile", "K_spatial_25m",
        "component_N_fertile", "component_N_events",
        "component_K", "component_M_target", "M_paternal_only",
        "M_frag"]].rename(columns={"component_id_within_loc": "component_id"})
    out_ev = DEFAULT_TABLES / "step29c_sampling_frag_aware_per_event.tsv"
    ev_out.to_csv(out_ev, sep="\t", index=False)
    print(f"[step29c] Wrote {out_ev}  ({len(ev_out)} events)")

    loc = build_per_location(ev)
    tidy = loc[[
        "locationID", "locationCode", "BL",
        "n_events", "n_components_25m", "total_N_fertile", "K_local",
        "M_current_uniform", "M_current_event_floor", "M_current",
        "sum_M_paternal_only", "M_frag_aware", "delta"]]
    out_loc = DEFAULT_TABLES / "step29c_sampling_comparison_per_location.tsv"
    tidy.to_csv(out_loc, sep="\t", index=False)
    print(f"[step29c] Wrote {out_loc}  ({len(tidy)} locations)")

    gained = int((loc["delta"] > 0).sum())
    saved  = int((loc["delta"] < 0).sum())
    same   = int((loc["delta"] == 0).sum())
    print(f"[step29c] Delta summary (fragmentation-aware − current):")
    print(f"  Locations MORE sampling  (Δ > 0): {gained}")
    print(f"  Locations LESS sampling  (Δ < 0): {saved}")
    print(f"  Locations unchanged      (Δ = 0): {same}")
    print(f"  Total mothers, current    : {int(loc['M_current'].sum())}")
    print(f"  Total mothers, frag-aware : {int(loc['M_frag_aware'].sum())}")
    print(f"  Median Δ per location    : {int(loc['delta'].median())}")
    print(f"  Range Δ per location     : "
          f"{int(loc['delta'].min())} to {int(loc['delta'].max())}")

    DEFAULT_FIGURES.mkdir(parents=True, exist_ok=True)
    plot_comparison(
        loc,
        out_png=DEFAULT_FIGURES / "step29c_sampling_comparison.png",
        out_pdf=DEFAULT_FIGURES / "step29c_sampling_comparison.pdf",
    )
    print(f"[step29c] Figure in {DEFAULT_FIGURES}/")


if __name__ == "__main__":
    main()
