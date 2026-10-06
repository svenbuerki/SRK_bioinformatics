"""Step 29a — build the Phase 5 population + slickspot crosswalk.

Phase 5 § A.5.0. Three-level spatial hierarchy for an annual plant:

    Species > Population (500 m, stable across years) > Deme (50 m,
    within-year) > Event.

L. papilliferum is annual and its above-ground footprint can switch
on/off between years with seed-bank germination; eventIDs are also
NOT stable across years (fresh barcodes each season — but locationID
is stable). So before we build the population graph we match
slickspots across years with a 10 m coordinate buffer that absorbs
typical consumer-GPS drift.

Algorithm
---------
1. Load 2025 + 2026 LEPA events (coords + n_fertile; same data-scope
   filters as § A.2) tagged with their source year.
2. Load the historical centroid of every DB locationID.
3. Cross-year slickspot match: within each locationID, 10 m
   haversine graph on 2025 ∪ 2026 events; components = slickspots.
4. Pool slickspot centroids + locationID centroids; build 500 m
   haversine graph; components = populations. populationID sorted
   by total n_fertile DESC then by westmost lon ASC.
5. Per (populationID, year), rebuild the 50 m deme partition on
   that year's events only.
6. Seed-cleaning priority ranking: populations present in BOTH 2025
   and 2026 above-ground, with ≥ 2 re-visited slickspots, ranked by
   total n_fertile.

Outputs (Tables/Phase5/)
------------------------
    step29a_population_crosswalk.tsv
    step29a_slickspot_summary.tsv
    step29a_population_summary.tsv
    step29a_seedcleaning_priority.tsv
    step29a_demes_per_population_year.tsv
"""
from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

from step28_seed_sampling_per_mother import (
    DEFAULT_DB, load_all_events, haversine_meters,
)
from srk_bl_constants import locationCode_to_bl

# ---- Tunables ----
SLICKSPOT_BUFFER_M = 10.0     # cross-year slickspot match radius
POPULATION_GAP_M   = 500.0    # population-level separation threshold
DEME_RADIUS_M      = 50.0     # within-year deme radius (same as Phase 5 primary)
YEARS              = (2025, 2026)

TABLES = Path("Tables/Phase5")


# ---------------------------------------------------------------------------
# 1. Data loaders
# ---------------------------------------------------------------------------
def load_all_year_events() -> pd.DataFrame:
    """Load 2025 + 2026 events into one frame with an `event_year` column."""
    frames = []
    for yr in YEARS:
        ev = load_all_events(DEFAULT_DB, year=yr)
        if len(ev):
            frames.append(ev.assign(event_year=yr))
    if not frames:
        raise SystemExit("[step29a] No events loaded for 2025 or 2026.")
    return pd.concat(frames, ignore_index=True)


def load_locationID_centroids(db_path: Path) -> pd.DataFrame:
    """Mean lat/lon across every event ever recorded at each DB locationID."""
    con = sqlite3.connect(str(db_path))
    try:
        df = pd.read_sql_query(
            """
            SELECT DISTINCT
                e.locationID              AS locationID,
                e.eventDecimalLatitude    AS lat,
                e.eventDecimalLongitude   AS lon
            FROM Events e
            JOIN Occurrences o ON o.eventID = e.eventID
            JOIN Taxonomy    t ON o.taxonID = t.taxonID
            WHERE t.genus='Lepidium' AND t.specificEpithet='papilliferum'
              AND e.eventDecimalLatitude  IS NOT NULL
              AND e.eventDecimalLongitude IS NOT NULL
              AND (o.provenance IS NULL OR o.provenance = 'in situ')
            """,
            con,
        )
    finally:
        con.close()
    df["lat"] = pd.to_numeric(df["lat"], errors="coerce")
    df["lon"] = pd.to_numeric(df["lon"], errors="coerce")
    df = df.dropna(subset=["lat", "lon"])
    df = df[(df["lat"].between(30, 55)) & (df["lon"].between(-125, -100))]
    return (df.groupby("locationID", as_index=False)
             .agg(lat=("lat", "mean"), lon=("lon", "mean")))


def load_locationID_to_locationCode() -> pd.Series:
    """Build locationID → locationCode map from the Phase 5 step29d
    summary (which has one row per locationID with its locationCode)."""
    tsv = TABLES / "step29d_mating_pool_summary.tsv"
    if not tsv.exists():
        return pd.Series(dtype=object)
    df = pd.read_csv(tsv, sep="\t", encoding="utf-8-sig")
    return df.set_index("locationID")["locationCode"]


# ---------------------------------------------------------------------------
# 2. Haversine connectivity → components (reusable)
# ---------------------------------------------------------------------------
def connected_components(lat: np.ndarray, lon: np.ndarray,
                          radius_m: float) -> list[list[int]]:
    """List of positional-index lists, one per connected component in
    the <= radius_m haversine adjacency graph. Isolated vertices form
    singleton components."""
    n = len(lat)
    if n == 0:
        return []
    d = haversine_meters(lat[:, None], lon[:, None], lat[None, :], lon[None, :])
    adj = (d > 0) & (d <= radius_m)
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
        comps.append(sorted(comp))
    return comps


# ---------------------------------------------------------------------------
# 3. Cross-year slickspot matching (within each locationID)
# ---------------------------------------------------------------------------
def assign_slickspots(events: pd.DataFrame) -> pd.DataFrame:
    """Add a `slickspotID` column. Within each DB locationID, build
    the 10 m haversine graph on 2025 ∪ 2026 events; connected
    components = slickspots. SlickspotIDs are globally unique
    integers."""
    out = events.copy().reset_index(drop=True)
    out["slickspotID"] = pd.NA
    next_id = 1
    for loc_id, sub in events.groupby("locationID"):
        idx = sub.index.to_numpy()
        lat = sub["lat"].to_numpy()
        lon = sub["lon"].to_numpy()
        comps = connected_components(lat, lon, SLICKSPOT_BUFFER_M)
        for comp in comps:
            for k in comp:
                out.at[idx[k], "slickspotID"] = next_id
            next_id += 1
    out["slickspotID"] = out["slickspotID"].astype("Int64")
    return out


def build_slickspot_summary(events_ss: pd.DataFrame) -> pd.DataFrame:
    """One row per slickspotID with centroid + per-year n_fertile +
    per-year event count + parent locationID."""
    def _agg_year(sub, yr):
        sel = sub[sub["event_year"] == yr]
        return int(sel["n_fertile"].sum()), int(len(sel))

    rows = []
    for ss_id, sub in events_ss.groupby("slickspotID"):
        n_f_25, n_e_25 = _agg_year(sub, 2025)
        n_f_26, n_e_26 = _agg_year(sub, 2026)
        rows.append({
            "slickspotID":      int(ss_id),
            "locationID":       int(sub["locationID"].iloc[0]),
            "lat":              float(sub["lat"].mean()),
            "lon":              float(sub["lon"].mean()),
            "n_fertile_2025":   n_f_25,
            "n_fertile_2026":   n_f_26,
            "n_events_2025":    n_e_25,
            "n_events_2026":    n_e_26,
            "both_years":       (n_e_25 > 0) and (n_e_26 > 0),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 4. 500 m population graph
# ---------------------------------------------------------------------------
def build_populations(slickspots: pd.DataFrame,
                       loc_centroids: pd.DataFrame
                       ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pool slickspot centroids + locationID centroids and run the
    500 m haversine graph. Returns (slickspots_with_populationID,
    pooled_points_with_populationID)."""
    # Pooled point set: slickspot centroids first (so we can map back
    # by positional index), then locationID centroids.
    ss_pts = slickspots[["lat", "lon"]].assign(
        source="slickspot_centroid",
        slickspotID=slickspots["slickspotID"].astype(int),
        locationID=slickspots["locationID"].astype(int),
        n_fert=slickspots["n_fertile_2025"] + slickspots["n_fertile_2026"],
    )
    loc_pts = loc_centroids[["lat", "lon"]].assign(
        source="locationID_centroid",
        slickspotID=pd.NA,
        locationID=loc_centroids["locationID"].astype(int),
        n_fert=0,
    )
    pts = pd.concat([ss_pts, loc_pts], ignore_index=True)

    comps = connected_components(pts["lat"].to_numpy(),
                                   pts["lon"].to_numpy(),
                                   POPULATION_GAP_M)
    pts["component_idx"] = -1
    for cidx, comp in enumerate(comps):
        pts.loc[comp, "component_idx"] = cidx

    # Rank components by total n_fert DESC then westmost lon ASC
    comp_rank = (pts.groupby("component_idx", as_index=False)
                   .agg(n_fert_total=("n_fert", "sum"),
                        westmost_lon=("lon", "min")))
    comp_rank = comp_rank.sort_values(
        ["n_fert_total", "westmost_lon"], ascending=[False, True]
    ).reset_index(drop=True)
    comp_rank["populationID"] = np.arange(1, len(comp_rank) + 1)
    comp_to_pop = dict(zip(comp_rank["component_idx"], comp_rank["populationID"]))

    pts["populationID"] = pts["component_idx"].map(comp_to_pop).astype(int)
    ss_with_pop = slickspots.copy()
    ss_with_pop["populationID"] = pts.loc[:len(ss_with_pop) - 1,
                                           "populationID"].to_numpy()
    return ss_with_pop, pts


# ---------------------------------------------------------------------------
# 5. Per-year deme rebuild inside each population
# ---------------------------------------------------------------------------
def rebuild_demes(events_ss: pd.DataFrame,
                   slickspot_to_pop: pd.Series) -> pd.DataFrame:
    out = events_ss.copy()
    out["populationID"] = out["slickspotID"].map(slickspot_to_pop).astype("Int64")
    rows: list[dict] = []
    for (pop_id, yr), sub in out.groupby(["populationID", "event_year"]):
        if pd.isna(pop_id):
            continue
        comps = connected_components(sub["lat"].to_numpy(),
                                       sub["lon"].to_numpy(),
                                       DEME_RADIUS_M)
        n_fert = sub["n_fertile"].to_numpy()
        for d_idx, comp in enumerate(comps, start=1):
            rows.append({
                "populationID":        int(pop_id),
                "year":                int(yr),
                "demeID":              d_idx,
                "n_events":            len(comp),
                "component_N_fertile": int(n_fert[comp].sum()),
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 6. Population summary + seed-cleaning priority
# ---------------------------------------------------------------------------
def population_summary(slickspots: pd.DataFrame,
                        demes: pd.DataFrame,
                        locid_to_loccode: pd.Series) -> pd.DataFrame:
    agg = (slickspots.groupby("populationID", as_index=False)
                     .agg(n_slickspots=("slickspotID", "nunique"),
                          n_slickspots_both_years=("both_years", "sum"),
                          n_fertile_2025=("n_fertile_2025", "sum"),
                          n_fertile_2026=("n_fertile_2026", "sum"),
                          n_events_2025=("n_events_2025", "sum"),
                          n_events_2026=("n_events_2026", "sum"),
                          locationIDs=("locationID",
                                        lambda s: ",".join(sorted(map(str, set(s)))))))
    agg["n_fertile_total"] = agg["n_fertile_2025"] + agg["n_fertile_2026"]
    agg["present_both_years"] = (
        (agg["n_fertile_2025"] > 0) & (agg["n_fertile_2026"] > 0)
    )
    agg["occupancy"] = np.where(
        agg["present_both_years"], "both_years",
        np.where(agg["n_fertile_2025"] > 0, "2025_only", "2026_only"),
    )
    agg["ratio_2026_over_2025"] = np.where(
        agg["n_fertile_2025"] > 0,
        agg["n_fertile_2026"] / agg["n_fertile_2025"].replace(0, np.nan),
        np.nan,
    )

    # Deme counts per year
    for yr in YEARS:
        sub = demes[demes["year"] == yr]
        counts = sub.groupby("populationID").size().rename(f"n_demes_{yr}")
        agg = agg.merge(counts, on="populationID", how="left")
        agg[f"n_demes_{yr}"] = agg[f"n_demes_{yr}"].fillna(0).astype(int)

    # locationCode list (legacy)
    def _codes(id_str: str) -> str:
        ids = [int(x) for x in id_str.split(",") if x]
        codes = sorted({locid_to_loccode.get(i, f"loc{i}") for i in ids})
        return ",".join(codes)
    agg["locationCodes"] = agg["locationIDs"].apply(_codes)

    agg = agg.sort_values("n_fertile_total", ascending=False).reset_index(drop=True)
    return agg


def seedcleaning_priority(pop_summary: pd.DataFrame) -> pd.DataFrame:
    cands = pop_summary[
        pop_summary["present_both_years"]
        & (pop_summary["n_slickspots_both_years"] >= 2)
    ].copy()
    cands = cands.sort_values("n_fertile_total", ascending=False).reset_index(drop=True)
    return cands


# ---------------------------------------------------------------------------
# 7. Driver
# ---------------------------------------------------------------------------
def main() -> None:
    events = load_all_year_events()
    print(f"[step29a] Loaded {len(events)} events "
          f"({(events['event_year'] == 2025).sum()} in 2025, "
          f"{(events['event_year'] == 2026).sum()} in 2026)")

    loc_centroids = load_locationID_centroids(DEFAULT_DB)
    print(f"[step29a] Loaded {len(loc_centroids)} historical locationID centroids")

    locid_to_loccode = load_locationID_to_locationCode()
    print(f"[step29a] locationID → locationCode map: "
          f"{len(locid_to_loccode)} entries from step29d_mating_pool_summary.tsv")

    events_ss = assign_slickspots(events)
    n_slickspots = int(events_ss["slickspotID"].nunique())
    print(f"[step29a] Cross-year slickspot match ({SLICKSPOT_BUFFER_M:.0f} m): "
          f"{n_slickspots} slickspots from {len(events)} events")

    slickspots = build_slickspot_summary(events_ss)
    slickspots, pts = build_populations(slickspots, loc_centroids)
    n_pop = int(slickspots["populationID"].nunique())
    print(f"[step29a] Population graph ({POPULATION_GAP_M:.0f} m): "
          f"{n_pop} populations from {n_slickspots} slickspots + "
          f"{len(loc_centroids)} locationID centroids")

    slickspot_to_pop = slickspots.set_index("slickspotID")["populationID"]
    demes = rebuild_demes(events_ss, slickspot_to_pop)
    print(f"[step29a] Deme partition ({DEME_RADIUS_M:.0f} m): "
          f"{len(demes)} per-year demes across {n_pop} populations")

    # ---- Crosswalk ----
    cw = events_ss.copy()
    cw["populationID"] = cw["slickspotID"].map(slickspot_to_pop).astype("Int64")
    cw["locationCode"] = cw["locationID"].map(locid_to_loccode)
    cw["BL"] = locationCode_to_bl(cw["locationCode"]).values
    cw = cw[[
        "populationID", "slickspotID", "event_year", "eventID",
        "locationID", "locationCode", "BL", "lat", "lon", "n_fertile",
    ]].sort_values(
        ["populationID", "slickspotID", "event_year"]
    ).reset_index(drop=True)
    (TABLES / "step29a_population_crosswalk.tsv").write_text("")
    cw.to_csv(TABLES / "step29a_population_crosswalk.tsv", sep="\t", index=False)
    print(f"[step29a] Wrote Tables/Phase5/step29a_population_crosswalk.tsv")

    # ---- Slickspot summary ----
    slickspots_out = slickspots.copy()
    slickspots_out["locationCode"] = slickspots_out["locationID"].map(locid_to_loccode)
    slickspots_out.to_csv(TABLES / "step29a_slickspot_summary.tsv",
                           sep="\t", index=False)
    print(f"[step29a] Wrote Tables/Phase5/step29a_slickspot_summary.tsv")

    # ---- Population summary ----
    pop_sum = population_summary(slickspots, demes, locid_to_loccode)
    pop_sum.to_csv(TABLES / "step29a_population_summary.tsv",
                    sep="\t", index=False)
    print(f"[step29a] Wrote Tables/Phase5/step29a_population_summary.tsv")

    # ---- Seed-cleaning priority ----
    cands = seedcleaning_priority(pop_sum)
    cands.to_csv(TABLES / "step29a_seedcleaning_priority.tsv",
                  sep="\t", index=False)
    print(f"[step29a] Wrote Tables/Phase5/step29a_seedcleaning_priority.tsv  "
          f"({len(cands)} candidate populations)")

    # ---- Per-population-year demes ----
    demes.to_csv(TABLES / "step29a_demes_per_population_year.tsv",
                  sep="\t", index=False)
    print(f"[step29a] Wrote Tables/Phase5/step29a_demes_per_population_year.tsv")

    # ---- Headline numbers ----
    occ = pop_sum["occupancy"].value_counts().to_dict()
    print()
    print(f"[step29a] ========== Population summary ==========")
    print(f"[step29a] {n_pop} populations total")
    print(f"[step29a]   both years above-ground   : {occ.get('both_years', 0)}")
    print(f"[step29a]   2025 only                 : {occ.get('2025_only', 0)}")
    print(f"[step29a]   2026 only                 : {occ.get('2026_only', 0)}")
    print(f"[step29a]   candidate (≥2 re-visited) : {len(cands)}")

    # ---- Legacy locationCode → populationID aggregation ----
    code_counts = (cw.dropna(subset=["locationCode"])
                     .groupby("locationCode")["populationID"].nunique())
    merged = code_counts[code_counts == 1]
    split = code_counts[code_counts > 1]
    print()
    print(f"[step29a] Legacy locationCode mapping: "
          f"{len(merged)} codes map to a single populationID; "
          f"{len(split)} codes span multiple populations.")

    # Any locationCodes that got MERGED into the same populationID?
    pop_counts = (cw.dropna(subset=["locationCode"])
                    .groupby("populationID")["locationCode"].nunique())
    merged_pops = pop_counts[pop_counts > 1]
    print(f"[step29a] Populations aggregating multiple locationCodes: "
          f"{len(merged_pops)} populations (merge ≥ 2 Phase 5 locationCodes)")


if __name__ == "__main__":
    main()
