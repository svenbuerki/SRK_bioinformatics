#!/usr/bin/env python3
"""Step 29c — Fragmentation-aware sampling allocation.

Uses the **50 m primary pollen-flight radius** (see
`step29a_pollinator_radius_sensitivity.py` for the biological
justification) to group each location's events into connected
components, then allocates mothers per COMPONENT (not per whole
location) using the same coupon-collector 90 %-probability-of-full-
detection rule as § B.4.1. Within each component, mothers are
distributed proportional to N_fertile, subject to a maternal-
genotype floor of ≥ 1 mother per event.

Dependency order: run **Step 28 → Step 29 → Step 29b → Step 29c**.

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
collector maths to each 50 m component, so:

  * A location that is one well-connected cluster → unchanged.
  * A location that is many isolated events → n_events mothers
    (private-allele floor was already forcing this).
  * A location that is a big cluster + a few isolated singletons
    (e.g. EO27-1, EO26-3) → the cluster only needs its coupon-
    collector M (~6), the isolated events each need 1, and the total
    can DROP substantially below the current recommendation while
    still delivering the same 90 %-see-every-allele guarantee.

Coupon-collector maths are scoped to each 50 m component rather than
the pooled location, so a location that is many isolated events is
allocated fewer mothers than the pooled Step 29 version would ask for,
and a location that is one well-connected cluster is unchanged.

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
from srk_bl_constants import (
    BL_COLORS, BL_ORDER, locationCode_to_bl, base_eo,
)

# Field → lab yield assumption used for Part C sampling design.
# Approximate LEPA germination rate under greenhouse conditions.
# Updating this constant propagates through n_seeds_to_germinate /
# n_seedlings_expected in step29c_partC_germplasmID_selection.tsv.
GERMINATION_RATE = 0.60

DEFAULT_TABLES = Path("Tables/Phase5")
DEFAULT_FIGURES = Path("figures/Phase5")

EVENT_TSV = DEFAULT_TABLES / "step28_events_spatial_neighborhood.tsv"
CURRENT_M_TSV = (DEFAULT_TABLES /
                 "step28_mothers_for_full_detection_by_location.tsv")
R_PRIMARY = 50.0
TARGET_PROB = 0.90


def connected_components_50m(sub: pd.DataFrame) -> list[list[int]]:
    """BFS on the 50 m adjacency graph. Returns positional-index lists,
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
        comps = connected_components_50m(sub)
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
                    n_components_50m=("component_id_within_loc",
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
    # Merge step29b connectivity so the figure label can show N_fert_eff
    # (drift-relevant mating pool) alongside the raw census.
    conn_tsv = DEFAULT_TABLES / "step29_location_connectivity.tsv"
    if conn_tsv.exists():
        conn = pd.read_csv(conn_tsv, sep="\t", encoding="utf-8-sig")
        conn_slim = conn[["locationID", "largest_component_share_50m"]]
        per_loc = per_loc.merge(conn_slim, on="locationID", how="left")
    # Count mothers already collected in the LEPA DB per location, so the
    # figure can show 'mothers available' vs 'mothers needed by design'
    # side-by-side. Shortage = target − available, floored at 0.
    field_recipe_tsv = DEFAULT_TABLES / "step29_field_team_sampling_recipe.tsv"
    if field_recipe_tsv.exists():
        rec = pd.read_csv(field_recipe_tsv, sep="\t", encoding="utf-8-sig")
        avail = (rec.groupby("locationID")
                    .size()
                    .rename("n_mothers_available_in_DB")
                    .reset_index())
        per_loc = per_loc.merge(avail, on="locationID", how="left")
        per_loc["n_mothers_available_in_DB"] = (
            per_loc["n_mothers_available_in_DB"].fillna(0).astype(int))
    else:
        per_loc["n_mothers_available_in_DB"] = 0
    per_loc["shortage"] = (per_loc["M_frag_aware"]
                            - per_loc["n_mothers_available_in_DB"]).clip(lower=0)
    per_loc["BL"] = locationCode_to_bl(per_loc["locationCode"]).values
    per_loc["BL"] = per_loc["BL"].fillna("Unassigned")
    # Sort within each BL by M_frag_aware ascending (small → large) so
    # every panel reads consistently top-to-bottom.
    return per_loc.sort_values(["BL", "M_frag_aware"], ascending=[True, True])


def select_germplasm_for_partC(per_event_df: pd.DataFrame,
                                field_recipe_tsv: Path,
                                seeds_per_mother_target: int = 15
                                ) -> pd.DataFrame:
    """Pick specific germplasmIDs from the LEPA DB to satisfy step29c's
    fragmentation-aware allocation for Part C testing — **at the 50 m
    component scale**, not per event.

    Rationale. Events inside the same 50 m connected component share the
    same pollen pool (that is the biological definition of connectivity
    at the primary pollinator radius). A mother sampled at event A
    therefore samples the same pool as a mother at event B when both
    belong to the same component. So the coupon-collector coverage
    guarantee for the component travels freely across its events, and
    a shortage at one event can be absorbed by picking extra mothers
    at a connected event in the same component. Only the ≥ 1-mother-
    per-event *maternal-genotype floor* is a strict per-event
    requirement, and even that is waived at events that hold no
    germplasm in the DB at all.

    Algorithm (per 50 m component):
      1. Collect every germplasmID in the DB across every event of the
         component.
      2. Enforce the maternal-genotype floor: at each event that has at
         least one germplasmID, reserve the top-seeded mother.
      3. Fill the remaining `component_M_target − (# events with
         germplasm)` slots from the leftover mothers in the component,
         sorted by `seeds_available` (descending) then germplasmID
         (ascending) for a deterministic tie-break.
      4. If the component's total germplasm is smaller than
         `component_M_target`, take every available mother and record
         the shortfall in `component_gap`.

    Returns one row per SELECTED germplasmID with columns:
        germplasmID, occurrenceID, eventID, locationID, locationCode,
        seeds_available, n_seeds_to_genotype (= min(15, seeds_available)),
        component_id_within_loc, component_M_target,
        component_n_available_in_DB, component_gap,
        selection_reason ('floor' or 'top_seeds'),
        selection_priority_within_component,
        event_M_frag (for reference — the per-event allocation upstream).
    """
    if not field_recipe_tsv.exists():
        raise SystemExit(
            f"Missing {field_recipe_tsv} — run step29 first to build "
            "the per-germplasmID recipe.")
    recipe = pd.read_csv(field_recipe_tsv, sep="\t", encoding="utf-8-sig")

    # Bring per-event M_frag AND component identity + component target
    # from step29c's per-event output.
    event_meta = per_event_df[[
        "locationID", "eventID", "M_frag",
        "component_id_within_loc", "component_M_target",
    ]].copy()
    merged = recipe.merge(event_meta, on=["locationID", "eventID"], how="left")
    merged["M_frag"] = merged["M_frag"].fillna(0).astype(int)
    # A germplasmID whose event has no component info (data cleanup
    # gaps) is assigned its own singleton component.
    missing = merged["component_id_within_loc"].isna()
    merged.loc[missing, "component_id_within_loc"] = -1
    merged.loc[missing, "component_M_target"] = merged.loc[missing, "M_frag"]
    merged["component_id_within_loc"] = merged["component_id_within_loc"].astype(int)
    merged["component_M_target"] = merged["component_M_target"].astype(int)

    # Component-level DB availability (all germplasmIDs at all events
    # in the same component).
    comp_avail = (merged.groupby(["locationID", "component_id_within_loc"])
                        .size().rename("component_n_available_in_DB")
                        .reset_index())
    merged = merged.merge(
        comp_avail, on=["locationID", "component_id_within_loc"], how="left")

    picked_rows = []
    # Iterate over components.
    for (loc_id, comp_id), grp in merged.groupby(
            ["locationID", "component_id_within_loc"]):
        target = int(grp["component_M_target"].iloc[0])
        n_avail = int(len(grp))
        # Sort component-wide by seeds_available DESC, then germplasmID ASC.
        grp = grp.sort_values(
            ["seeds_available", "germplasmID"],
            ascending=[False, True]).reset_index(drop=True)

        # Step 1 — floor: one top-seeded mother per event with germplasm.
        floor_ids = (grp.groupby("eventID", sort=False)
                        .head(1)["germplasmID"].tolist())
        floor_set = set(floor_ids)
        floor_picked = grp[grp["germplasmID"].isin(floor_set)].copy()
        floor_picked["selection_reason"] = "floor"

        # Step 2 — fill remaining coverage slots.
        remaining_target = max(target - len(floor_picked), 0)
        leftover = grp[~grp["germplasmID"].isin(floor_set)]
        extra_picked = leftover.head(remaining_target).copy()
        extra_picked["selection_reason"] = "top_seeds"

        picked = pd.concat([floor_picked, extra_picked], ignore_index=True)
        picked["component_gap"] = max(target - n_avail, 0)
        picked_rows.append(picked)

    if picked_rows:
        selected = pd.concat(picked_rows, ignore_index=True)
    else:
        selected = merged.iloc[0:0].copy()
        selected["selection_reason"] = pd.Series(dtype="object")
        selected["component_gap"] = pd.Series(dtype=int)

    # Rank within component by seeds_available desc for the output.
    selected = selected.sort_values(
        ["locationID", "component_id_within_loc",
         "selection_reason", "seeds_available", "germplasmID"],
        ascending=[True, True, True, False, True]).reset_index(drop=True)
    selected["selection_priority_within_component"] = (
        selected.groupby(["locationID", "component_id_within_loc"]).cumcount() + 1)

    # Part C genotyping is on SEEDLINGS, not seeds (user-confirmed
    # 2026-10-03). The tetraploid Rule 2 coupon-collector floor remains
    # 15 **seedlings** per mother; given a ~60 % germination rate we
    # need to germinate ~ceil(15 / 0.60) = 25 seeds per mother to hit
    # that seedling target in expectation. Columns:
    #   seedlings_target          — Rule 2 target (15 seedlings)
    #   n_seeds_to_germinate      — min(25, seeds_available)
    #   n_seedlings_expected      — round(n_seeds_to_germinate * 0.60)
    #   n_seedlings_to_genotype   — min(15, n_seedlings_expected)
    seeds_to_germinate_target = int(np.ceil(
        seeds_per_mother_target / GERMINATION_RATE))  # 25
    selected["seedlings_target"] = seeds_per_mother_target
    selected["n_seeds_to_germinate"] = np.minimum(
        selected["seeds_available"].astype(int),
        seeds_to_germinate_target,
    )
    selected["n_seedlings_expected"] = np.rint(
        selected["n_seeds_to_germinate"] * GERMINATION_RATE
    ).astype(int)
    selected["n_seedlings_to_genotype"] = np.minimum(
        selected["n_seedlings_expected"],
        seeds_per_mother_target,
    )
    # Backward-compat column (deprecated — kept for one release).
    selected["n_seeds_to_genotype"] = selected["n_seedlings_to_genotype"]

    # Project-wide EOID = base EO code (strips Phase 5 dash suffixes:
    # EO18-7 → EO18, EO27RT → EO27, EO27-3 → EO27). The TSV's primary
    # sort key is EOID; locationID disambiguates within-EO splits.
    selected["EOID"] = selected["locationCode"].astype(str).map(base_eo)

    keep = [
        "EOID", "locationCode", "locationID",
        "component_id_within_loc", "germplasmID", "occurrenceID",
        "eventID", "seeds_available",
        "seedlings_target", "n_seeds_to_germinate",
        "n_seedlings_expected", "n_seedlings_to_genotype",
        "n_seeds_to_genotype",
        "component_M_target", "component_n_available_in_DB",
        "component_gap", "selection_reason",
        "selection_priority_within_component",
        "M_frag",
    ]
    keep = [c for c in keep if c in selected.columns]
    out = (selected[keep]
           .rename(columns={"M_frag": "event_M_frag"})
           .sort_values(
               ["EOID", "locationID", "component_id_within_loc",
                "germplasmID"])
           .reset_index(drop=True))
    return out


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
               .sort_values("M_frag_aware", ascending=True)
               .reset_index(drop=True))
        y = np.arange(len(sub))
        colour = palette[bl]
        h = 0.36
        # Solid bar — fragmentation-aware design target (M_frag).
        # Visually appears BELOW the light bar because barh with y - h/2
        # sits on the lower half of each row.
        ax.barh(y - h/2, sub["M_frag_aware"], h, color=colour, alpha=0.9,
                edgecolor=colour, linewidth=1.0)
        # Light bar — mothers already collected and stored in the LEPA DB.
        # Appears ABOVE the solid bar in each row.
        ax.barh(y + h/2, sub["n_mothers_available_in_DB"], h,
                color=colour, alpha=0.35, edgecolor=colour, linewidth=1.0)
        # Row label: locationCode + spatial context matching Figure 2b.
        share = sub.get("largest_component_share_50m",
                        pd.Series([1.0] * len(sub))).fillna(1.0)
        eff = (sub["total_N_fertile"].astype(float) * share).round().astype(int)
        labels = [f"{r['locationCode']}  "
                  f"(events = {int(r['n_events'])}, "
                  f"50 m components = {int(r['n_components_50m'])}, "
                  f"census = {int(r['total_N_fertile'])}, "
                  f"effective = {int(e)})"
                  for (_, r), e in zip(sub.iterrows(), eff)]
        # Status annotation to the right of each row:
        #   short by N  → red, DB lacks mothers to hit the target
        #   covered     → grey, DB has enough (or more) already
        for i, r in sub.iterrows():
            short = int(r["shortage"])
            end_x = max(int(r["M_frag_aware"]),
                        int(r["n_mothers_available_in_DB"]))
            if short > 0:
                ax.text(end_x + 0.7, i, f"short by {short}",
                        fontsize=8, color="#c94b4b",
                        va="center", ha="left", fontweight="bold")
            else:
                ax.text(end_x + 0.7, i, "covered",
                        fontsize=8, color="#3c8f4c",
                        va="center", ha="left")
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_ylim(-0.7, len(sub) - 0.3)
        ax.text(1.01, 0.5, bl, transform=ax.transAxes,
                fontsize=13, fontweight="bold", color=colour,
                va="center", ha="left")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Figure-level legend with neutral grey handles so it doesn't look
    # BL-specific. Order matches the visual top-to-bottom layout in
    # each row: light 'available' bar on top, solid 'needed' bar below.
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="#555555", alpha=0.35, edgecolor="#333",
              label="Mothers already in the LEPA DB (light bar, on top)"),
        Patch(facecolor="#555555", alpha=0.9, edgecolor="#333",
              label="Mothers needed by § B.4.2 target (solid bar, below)"),
    ]
    xmax = float(max(loc_df["M_frag_aware"].max(),
                     loc_df["n_mothers_available_in_DB"].max()))
    axes[-1].set_xlim(0, xmax + 18)
    axes[-1].set_xlabel(
        "Mothers per location  —  target (§ B.4.2) vs already collected "
        "in the LEPA DB",
        fontsize=11)
    n_short_loc  = int((loc_df["shortage"] > 0).sum())
    total_short  = int(loc_df["shortage"].sum())
    total_target = int(loc_df["M_frag_aware"].sum())
    total_avail  = int(loc_df["n_mothers_available_in_DB"].sum())
    fig.suptitle(
        f"Mothers per LEPA location — do we have what the design asks for?\n"
        f"Target: {total_target}. In LEPA DB: {total_avail}. "
        f"Short: {n_short_loc} locations ({total_short} mothers, "
        f"2026 top-up).",
        fontsize=11, y=0.998)
    fig.legend(handles=legend_handles, loc="upper center",
                bbox_to_anchor=(0.5, 0.94), ncol=2,
                fontsize=10, frameon=True)
    fig.tight_layout(rect=[0, 0, 0.94, 0.92])
    fig.savefig(out_png, dpi=200); fig.savefig(out_pdf); plt.close(fig)


def main() -> None:
    for p in (EVENT_TSV, CURRENT_M_TSV):
        if not p.exists():
            raise SystemExit(f"Missing input {p}")

    ev = build_per_event()
    ev_out = ev[[
        "locationID", "eventID", "component_id_within_loc",
        "n_fertile", "K_spatial_50m",
        "component_N_fertile", "component_N_events",
        "component_K", "component_M_target", "M_paternal_only",
        "M_frag"]].rename(columns={"component_id_within_loc": "component_id"})
    out_ev = DEFAULT_TABLES / "step29c_sampling_frag_aware_per_event.tsv"
    ev_out.to_csv(out_ev, sep="\t", index=False)
    print(f"[step29c] Wrote {out_ev}  ({len(ev_out)} events)")

    # ---- Dedicated event → 50 m component lookup ----
    # A lean TSV joining (locationID, eventID) → component_id_50m so
    # any downstream query 'which events share a pollen pool?' has a
    # single canonical source. One row per event.
    field_recipe_tsv_for_lookup = (
        DEFAULT_TABLES / "step29_field_team_sampling_recipe.tsv")
    lookup = ev[[
        "locationID", "eventID", "component_id_within_loc",
        "n_fertile", "component_N_fertile", "component_N_events",
        "component_K", "component_M_target"]].rename(columns={
            "component_id_within_loc": "component_id_50m",
            "n_fertile":                "event_n_fertile",
        }).drop_duplicates(subset=["locationID", "eventID"]).copy()
    # locationCode attached from the field recipe. Deduplicate to one
    # (locationID, eventID, locationCode) row before merging so the
    # lookup doesn't blow up to per-mother rows.
    if field_recipe_tsv_for_lookup.exists():
        rec = pd.read_csv(field_recipe_tsv_for_lookup, sep="\t",
                          encoding="utf-8-sig")
        loc_key = (rec[["locationID", "locationCode", "eventID"]]
                     .drop_duplicates(subset=["locationID", "eventID"]))
        lookup = lookup.merge(loc_key, on=["locationID", "eventID"], how="left")
    lookup = lookup[[
        c for c in ["locationID", "locationCode", "eventID",
                    "component_id_50m", "event_n_fertile",
                    "component_N_fertile", "component_N_events",
                    "component_K", "component_M_target"]
        if c in lookup.columns
    ]].sort_values(["locationCode", "component_id_50m", "eventID"] \
                    if "locationCode" in lookup.columns \
                    else ["locationID", "component_id_50m", "eventID"]
                    ).reset_index(drop=True)
    out_lookup = DEFAULT_TABLES / "step29c_event_to_component_50m.tsv"
    lookup.to_csv(out_lookup, sep="\t", index=False)
    print(f"[step29c] Wrote {out_lookup}  "
          f"({len(lookup)} events across "
          f"{lookup.groupby(['locationID','component_id_50m']).ngroups} components)")

    loc = build_per_location(ev)
    tidy_cols = [
        "locationID", "locationCode", "BL",
        "n_events", "n_components_50m", "total_N_fertile",
        "largest_component_share_50m",
        "K_local",
        "M_current_uniform", "M_current_event_floor", "M_current",
        "sum_M_paternal_only", "M_frag_aware",
        "n_mothers_available_in_DB", "shortage",
    ]
    tidy = loc[[c for c in tidy_cols if c in loc.columns]]
    out_loc = DEFAULT_TABLES / "step29c_sampling_comparison_per_location.tsv"
    tidy.to_csv(out_loc, sep="\t", index=False)
    print(f"[step29c] Wrote {out_loc}  ({len(tidy)} locations)")

    n_target  = int(loc["M_frag_aware"].sum())
    n_avail   = int(loc["n_mothers_available_in_DB"].sum())
    n_short_loc = int((loc["shortage"] > 0).sum())
    n_covered   = int((loc["shortage"] == 0).sum())
    n_short_m   = int(loc["shortage"].sum())
    print(f"[step29c] Fragmentation-aware sampling summary:")
    print(f"  Target (M_frag_aware)               : {n_target} mothers")
    print(f"  Already in LEPA DB (n_available)    : {n_avail} mothers")
    print(f"  Locations fully covered by DB       : {n_covered} / {len(loc)}")
    print(f"  Locations short of target           : {n_short_loc} / {len(loc)}")
    print(f"  Total mothers short (2026 top-up)   : {n_short_m}")

    # ---- Part C germplasmID selection ----
    # Pick specific germplasmIDs from the LEPA DB to satisfy the
    # fragmentation-aware allocation at the 50 m COMPONENT scale
    # (mothers within a component pool their coverage; only the
    # per-event maternal-genotype floor is a strict per-event rule).
    field_recipe_tsv = DEFAULT_TABLES / "step29_field_team_sampling_recipe.tsv"
    partC = select_germplasm_for_partC(ev, field_recipe_tsv)
    out_partC = DEFAULT_TABLES / "step29c_partC_germplasmID_selection.tsv"
    partC.to_csv(out_partC, sep="\t", index=False)
    n_selected   = len(partC)
    total_seeds  = int(partC["n_seeds_to_germinate"].sum())
    total_slings = int(partC["n_seedlings_to_genotype"].sum())
    per_comp_gap = (partC.groupby(["locationID", "component_id_within_loc"])
                          ["component_gap"].first())
    n_comps_short  = int((per_comp_gap > 0).sum())
    total_shortage = int(per_comp_gap.sum())
    m_frag_target  = int(ev["M_frag"].sum())
    print(f"[step29c] Wrote {out_partC}  "
          f"({n_selected} germplasmIDs selected of {m_frag_target} "
          f"M_frag target; {total_seeds} seeds to germinate at "
          f"{GERMINATION_RATE:.0%} germination → "
          f"{total_slings} seedlings to genotype; "
          f"{n_comps_short} components short by {total_shortage} mothers total)")

    DEFAULT_FIGURES.mkdir(parents=True, exist_ok=True)
    plot_comparison(
        loc,
        out_png=DEFAULT_FIGURES / "step29c_sampling_comparison.png",
        out_pdf=DEFAULT_FIGURES / "step29c_sampling_comparison.pdf",
    )
    print(f"[step29c] Figure in {DEFAULT_FIGURES}/")


if __name__ == "__main__":
    main()
