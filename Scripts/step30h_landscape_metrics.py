"""Step 30h / Phase II — landscape + within-population fragmentation metrics.

Phase 5 § A.4.5. The main fragmentation signal in LEPA lives at the
LANDSCAPE scale (between populations), not within a population at
the deme level. This script computes the two layers of fragmentation
metrics that feed the BL definition:

LANDSCAPE (between-population) — primary input to BL definition
    - nearest-neighbor population distance (per population)
    - cluster membership at k_optimal (per population, from Phase I)
    - mean within-BL inter-population distance (per BL)
    - within-BL connectivity at threshold distances (per BL)
    - BL convex-hull area + mean centroid (per BL)
    - min between-BL distance = isolation (per BL)

WITHIN-POPULATION (deme-level) — secondary diagnostic
    - n_demes per year (already in step29a output)
    - mean inter-deme distance per year (computed here)

Reads
-----
Tables/Phase5/step30h_cluster_assignments.tsv       (Phase I)
Tables/Phase5/step30h_population_distances.tsv      (Phase I)
Tables/Phase5/step29a_demes_per_population_year.tsv
Tables/Phase5/step29a_population_crosswalk.tsv

Outputs
-------
Tables/Phase5/step30h_landscape_per_population.tsv
Tables/Phase5/step30h_landscape_per_BL.tsv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

from step28_seed_sampling_per_mother import haversine_meters

TABLES = Path("Tables/Phase5")

CONNECTIVITY_THRESHOLDS_M = [1_000, 2_000, 5_000, 10_000]


# ---------------------------------------------------------------------------
# 1. Load inputs
# ---------------------------------------------------------------------------
def load_inputs():
    cents = pd.read_csv(TABLES / "step30h_cluster_assignments.tsv",
                         sep="\t", encoding="utf-8-sig")
    long_D = pd.read_csv(TABLES / "step30h_population_distances.tsv",
                          sep="\t", encoding="utf-8-sig")
    demes_py = pd.read_csv(TABLES / "step29a_demes_per_population_year.tsv",
                            sep="\t", encoding="utf-8-sig")
    cw = pd.read_csv(TABLES / "step29a_population_crosswalk.tsv",
                      sep="\t", encoding="utf-8-sig")
    return cents, long_D, demes_py, cw


# ---------------------------------------------------------------------------
# 2. Per-population landscape metrics (nearest-neighbor etc.)
# ---------------------------------------------------------------------------
def per_population_landscape(cents: pd.DataFrame,
                                long_D: pd.DataFrame) -> pd.DataFrame:
    """Nearest-neighbor population distance (across all other
    populations regardless of BL) + nearest same-BL population
    distance. BL labels from the k_optimal column."""
    cl = cents[["populationID", "cluster_k_optimal"]].copy()

    # Pairwise distance → nearest neighbor
    # long_D has one row per unordered pair (A < B); expand to both
    # directions so each population sees all others.
    both = pd.concat([
        long_D.rename(columns={"populationA": "pop",
                                  "populationB": "other"}),
        long_D.rename(columns={"populationB": "pop",
                                  "populationA": "other"}),
    ])
    both = both.merge(cl.rename(columns={"populationID": "pop",
                                           "cluster_k_optimal": "pop_cluster"}),
                       on="pop", how="left")
    both = both.merge(cl.rename(columns={"populationID": "other",
                                           "cluster_k_optimal": "other_cluster"}),
                       on="other", how="left")
    both["same_cluster"] = both["pop_cluster"] == both["other_cluster"]

    agg_all = (both.groupby("pop", as_index=False)
                    .agg(nn_distance_m=("distance_m", "min"),
                         nn_populationID=("distance_m", lambda s: int(
                             both.loc[s.idxmin(), "other"]))))
    same = both[both["same_cluster"]]
    agg_same = (same.groupby("pop", as_index=False)
                     .agg(nn_same_BL_distance_m=("distance_m", "min")))
    out = (cl.rename(columns={"populationID": "pop"})
             .merge(agg_all, on="pop", how="left")
             .merge(agg_same, on="pop", how="left"))
    out = out.rename(columns={"pop": "populationID"})
    return out


# ---------------------------------------------------------------------------
# 3. Per-BL landscape metrics (connectivity, area, isolation)
# ---------------------------------------------------------------------------
def per_bl_landscape(cents: pd.DataFrame,
                       long_D: pd.DataFrame) -> pd.DataFrame:
    cl = cents[["populationID", "lat", "lon", "n_fert_total",
                 "cluster_k_optimal"]].copy()

    # Attach cluster labels to the pairwise distance table
    both = long_D.copy()
    both = both.merge(cl[["populationID", "cluster_k_optimal"]]
                        .rename(columns={"populationID": "populationA",
                                           "cluster_k_optimal": "clA"}),
                       on="populationA")
    both = both.merge(cl[["populationID", "cluster_k_optimal"]]
                        .rename(columns={"populationID": "populationB",
                                           "cluster_k_optimal": "clB"}),
                       on="populationB")

    rows: list[dict] = []
    for cl_id, sub in cl.groupby("cluster_k_optimal"):
        n_pop = len(sub)
        pop_ids = set(sub["populationID"])

        # Within-BL distances
        within = both[(both["clA"] == cl_id) & (both["clB"] == cl_id)]
        if len(within):
            d_within = within["distance_m"].to_numpy()
        else:
            d_within = np.array([])

        # Between-BL distances (to any other BL, min)
        between = both[(both["clA"] == cl_id) ^ (both["clB"] == cl_id)]
        min_between = float(between["distance_m"].min()) if len(between) else float("nan")

        # Convex-hull area of population centroids (approximate — use
        # Euclidean area in degrees then multiply by local scaling)
        hull_area_km2 = float("nan")
        if n_pop >= 3:
            pts = sub[["lon", "lat"]].to_numpy()
            try:
                hull = ConvexHull(pts)
                deg2_area = hull.volume   # 2D area in deg^2
                # Approx km^2 at this latitude
                mean_lat = np.deg2rad(float(sub["lat"].mean()))
                km_per_deg_lat = 110.574
                km_per_deg_lon = 111.320 * np.cos(mean_lat)
                hull_area_km2 = float(deg2_area * km_per_deg_lat * km_per_deg_lon)
            except Exception:
                pass

        row = {
            "BL_cluster_id":            int(cl_id),
            "n_populations":            int(n_pop),
            "total_N_fert":             int(sub["n_fert_total"].sum()),
            "centroid_lat":             float(sub["lat"].mean()),
            "centroid_lon":             float(sub["lon"].mean()),
            "convex_hull_area_km2":     round(hull_area_km2, 3)
                                          if hull_area_km2 == hull_area_km2
                                          else float("nan"),
            "mean_within_BL_dist_m":    (float(np.mean(d_within))
                                          if len(d_within) else float("nan")),
            "median_within_BL_dist_m":  (float(np.median(d_within))
                                          if len(d_within) else float("nan")),
            "max_within_BL_dist_m":     (float(np.max(d_within))
                                          if len(d_within) else float("nan")),
            "min_between_BL_dist_m":    min_between,
        }
        # Within-BL connectivity at thresholds
        if len(d_within):
            total_pairs = len(d_within)
            for t in CONNECTIVITY_THRESHOLDS_M:
                frac = (d_within <= t).sum() / total_pairs
                row[f"frac_within_BL_pairs_le_{t}m"] = round(float(frac), 3)
        else:
            for t in CONNECTIVITY_THRESHOLDS_M:
                row[f"frac_within_BL_pairs_le_{t}m"] = float("nan")

        rows.append(row)

    return pd.DataFrame(rows).sort_values("BL_cluster_id").reset_index(drop=True)


# ---------------------------------------------------------------------------
# 4. Within-population deme-level fragmentation (secondary diagnostic)
# ---------------------------------------------------------------------------
def within_population_deme_frag(demes_py: pd.DataFrame,
                                   cw: pd.DataFrame) -> pd.DataFrame:
    """For each (populationID, year), compute n_demes AND the mean
    pairwise deme-centroid distance. Deme centroids are per-year
    aggregates of event coordinates in that deme.

    `step29a_demes_per_population_year.tsv` has one row per (pop, year,
    demeID) with `component_N_fertile` but no coordinates — we recover
    deme centroids by re-running the 50 m connected-component logic
    on the crosswalk."""
    from srk_bl_constants import make_location_label  # noqa: F401  (dependency check)
    # Reuse the connected_components function from step29a without
    # importing matplotlib etc.
    from step29a_populations import connected_components, DEME_RADIUS_M

    cw_both = cw[cw["event_year"].isin((2025, 2026))].copy()
    cw_both = cw_both.dropna(subset=["lat", "lon"])

    rows: list[dict] = []
    for (pop_id, yr), g in cw_both.groupby(["populationID", "event_year"]):
        lat = g["lat"].to_numpy()
        lon = g["lon"].to_numpy()
        n_fert = g["n_fertile"].to_numpy()
        comps = connected_components(lat, lon, DEME_RADIUS_M)
        if not comps:
            continue
        deme_centroids = []
        for comp in comps:
            w = n_fert[comp]
            w = w if w.sum() > 0 else np.ones_like(w)
            deme_centroids.append((
                float(np.average(lat[comp], weights=w)),
                float(np.average(lon[comp], weights=w)),
            ))
        n_demes = len(deme_centroids)
        if n_demes >= 2:
            d_lat = np.array([c[0] for c in deme_centroids])
            d_lon = np.array([c[1] for c in deme_centroids])
            D = haversine_meters(d_lat[:, None], d_lon[:, None],
                                   d_lat[None, :], d_lon[None, :])
            # upper triangle only
            iu = np.triu_indices(n_demes, k=1)
            mean_d = float(np.mean(D[iu]))
            max_d  = float(np.max(D[iu]))
        else:
            mean_d = 0.0
            max_d  = 0.0
        rows.append({
            "populationID":            int(pop_id),
            "year":                    int(yr),
            "n_demes":                 n_demes,
            "mean_inter_deme_dist_m":  round(mean_d, 1),
            "max_inter_deme_dist_m":   round(max_d, 1),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 5. Driver
# ---------------------------------------------------------------------------
def main() -> None:
    cents, long_D, demes_py, cw = load_inputs()

    # ---- Per-population landscape ----
    pop_land = per_population_landscape(cents, long_D)

    # ---- Within-population deme fragmentation ----
    deme_frag = within_population_deme_frag(demes_py, cw)
    # Pivot per-year deme frag to wide
    wide = deme_frag.pivot(index="populationID", columns="year",
                             values=["n_demes",
                                      "mean_inter_deme_dist_m",
                                      "max_inter_deme_dist_m"])
    wide.columns = [f"{a}_{int(b)}" for a, b in wide.columns]
    wide = wide.reset_index()
    pop_land = pop_land.merge(wide, on="populationID", how="left")

    # Add N_fert / event counts for context
    pop_land = pop_land.merge(
        cents[["populationID", "lat", "lon", "n_fert_total", "n_events"]],
        on="populationID", how="left",
    )
    pop_land = pop_land.sort_values(
        ["cluster_k_optimal", "populationID"]
    ).reset_index(drop=True)
    pop_land.to_csv(TABLES / "step30h_landscape_per_population.tsv",
                     sep="\t", index=False)
    print(f"[step30h-L] Wrote step30h_landscape_per_population.tsv "
          f"({len(pop_land)} rows)")

    # ---- Per-BL landscape ----
    bl_land = per_bl_landscape(cents, long_D)
    bl_land.to_csv(TABLES / "step30h_landscape_per_BL.tsv",
                    sep="\t", index=False)
    print(f"[step30h-L] Wrote step30h_landscape_per_BL.tsv "
          f"({len(bl_land)} BLs)")

    # ---- Headline ----
    print()
    print("Per-BL landscape metrics (k_optimal):")
    cols_show = ["BL_cluster_id", "n_populations", "total_N_fert",
                  "convex_hull_area_km2",
                  "mean_within_BL_dist_m", "min_between_BL_dist_m",
                  "frac_within_BL_pairs_le_5000m"]
    print(bl_land[cols_show].round(1).to_string(index=False))


if __name__ == "__main__":
    main()
