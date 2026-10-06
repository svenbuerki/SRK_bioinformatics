"""Step 30h / Phase III — define new BL + renumber populationIDs.

Phase 5 § A.4.5. Takes the Ward's-D2 cluster assignments from Phase I
and the landscape metrics from Phase II, applies the ordering rule
chosen by the user, and produces:

  1. **BL ordering rule: area DESC → connectivity DESC** (same spirit
     as the external LEPA_EO_spatial_clustering repo). The 5
     Ward clusters are relabeled BL1..BL5 under this rule.

  2. **Within-BL population ordering: latitude DESC (N → S).** New
     populationID = sequential across (BL, within-BL order).

Then rewrites every Phase 5 population-level TSV with the new
populationIDs and BL labels, preserving the old populationID (and
old cluster id) in a crosswalk for traceability. Pure remap of
existing output — no simulation is re-run.

Outputs
-------
Tables/Phase5/step30h_bl_definition.tsv
    one row per BL with ordering keys (area, connectivity),
    centroid, population count, total N_fert.
Tables/Phase5/step30h_populationID_crosswalk.tsv
    old → new populationID + old cluster → new BL mapping.

Rewrites (new populationID + BL columns):
    step29a_population_crosswalk.tsv
    step29a_population_summary.tsv
    step29a_slickspot_summary.tsv
    step29a_demes_per_population_year.tsv
    step29a_seedcleaning_priority.tsv
    step30g_prediction_population_year.tsv
    step30g_across_year_comparison.tsv
    step30g_stable_candidate_pair.tsv
    step30g_crash_candidates.tsv
    step30g_populations_classified.tsv
    step30h_cluster_assignments.tsv
    step30h_landscape_per_population.tsv
    step30h_landscape_per_BL.tsv

A backup copy of every rewritten TSV is first written to
Tables/Phase5/_bl_renumber_backup/ so the pre-Phase-III state is
recoverable.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

TABLES = Path("Tables/Phase5")
BACKUP_DIR = TABLES / "_bl_renumber_backup"

# Which column in each file carries the OLD populationID that we
# need to remap. "cluster_k_optimal" or "BL_cluster_id" are the OLD
# cluster labels remapped to NEW BL names.
POP_ID_TSVS = {
    "step29a_population_crosswalk.tsv":        "populationID",
    "step29a_population_summary.tsv":          "populationID",
    "step29a_slickspot_summary.tsv":           "populationID",
    "step29a_demes_per_population_year.tsv":   "populationID",
    "step29a_seedcleaning_priority.tsv":       "populationID",
    "step30g_prediction_population_year.tsv":  "populationID",
    "step30g_across_year_comparison.tsv":      "populationID",
    "step30g_stable_candidate_pair.tsv":       "populationID",
    "step30g_crash_candidates.tsv":            "populationID",
    "step30g_populations_classified.tsv":      "populationID",
    "step30h_cluster_assignments.tsv":         "populationID",
    "step30h_landscape_per_population.tsv":    "populationID",
}


# ---------------------------------------------------------------------------
# 1. BL definition — ordering rule
# ---------------------------------------------------------------------------
def define_bl(bl_land: pd.DataFrame) -> pd.DataFrame:
    """Apply BL ordering: area DESC → connectivity DESC → cluster id ASC.

    'connectivity' = frac_within_BL_pairs_le_5000m.
    """
    df = bl_land.copy()
    df["sort_connectivity"] = df["frac_within_BL_pairs_le_5000m"].fillna(0)
    df = df.sort_values(
        ["convex_hull_area_km2", "sort_connectivity", "BL_cluster_id"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    df["BL"] = [f"BL{i + 1}" for i in range(len(df))]
    df["BL_rank"] = range(1, len(df) + 1)
    return df


# ---------------------------------------------------------------------------
# 2. Population renumbering
# ---------------------------------------------------------------------------
def renumber_populations(cl: pd.DataFrame,
                           bl_def: pd.DataFrame) -> pd.DataFrame:
    """Within each new BL, sort populations by lat DESC (N → S) then
    lon ASC (W → E) as tiebreak. New populationID = running index."""
    cluster_to_bl = dict(zip(bl_def["BL_cluster_id"], bl_def["BL"]))
    bl_to_rank    = dict(zip(bl_def["BL"], bl_def["BL_rank"]))

    df = cl.copy()
    df["BL"] = df["cluster_k_optimal"].map(cluster_to_bl)
    df["BL_rank"] = df["BL"].map(bl_to_rank)
    df = df.sort_values(
        ["BL_rank", "lat", "lon"],
        ascending=[True, False, True],
    ).reset_index(drop=True)
    df["populationID_new"] = range(1, len(df) + 1)
    df = df.rename(columns={"populationID": "populationID_old"})
    return df[[
        "populationID_new", "populationID_old", "BL",
        "cluster_k_optimal", "lat", "lon",
    ]]


# ---------------------------------------------------------------------------
# 3. TSV rewriters
# ---------------------------------------------------------------------------
def backup_tsvs() -> None:
    BACKUP_DIR.mkdir(exist_ok=True)
    for name in POP_ID_TSVS:
        src = TABLES / name
        if src.exists():
            shutil.copy2(src, BACKUP_DIR / name)
    # also backup the per-BL landscape TSV (gets rewritten but has
    # BL_cluster_id as its key, not populationID — handled separately)
    for extra in ("step30h_landscape_per_BL.tsv",):
        src = TABLES / extra
        if src.exists():
            shutil.copy2(src, BACKUP_DIR / extra)
    print(f"[step30h-III] Backed up {len(POP_ID_TSVS) + 1} TSVs → {BACKUP_DIR.name}/")


def remap_population_tsvs(crosswalk: pd.DataFrame) -> None:
    """Remap the populationID column + add a BL column in every TSV
    that carries populationID."""
    pop_map = dict(zip(crosswalk["populationID_old"],
                        crosswalk["populationID_new"]))
    bl_map  = dict(zip(crosswalk["populationID_old"],
                        crosswalk["BL"]))
    for name, col in POP_ID_TSVS.items():
        path = TABLES / name
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t", encoding="utf-8-sig")
        if col not in df.columns:
            continue
        df["populationID_old"] = df[col]
        df[col] = df[col].map(pop_map).astype("Int64")
        # Add BL column if not already present
        if "BL" not in df.columns:
            df.insert(1, "BL", df["populationID_old"].map(bl_map))
        # Reorder: put the renamed populationID + BL + populationID_old at the front
        front = [col, "BL"]
        if "populationID_old" in df.columns and "populationID_old" not in front:
            front.append("populationID_old")
        rest = [c for c in df.columns if c not in front]
        df = df[front + rest].sort_values(col).reset_index(drop=True)
        df.to_csv(path, sep="\t", index=False)
    print(f"[step30h-III] Remapped {len(POP_ID_TSVS)} population-keyed TSVs")


def rewrite_landscape_per_BL(bl_def: pd.DataFrame) -> None:
    """Rewrite step30h_landscape_per_BL.tsv so each row has the NEW
    BL name as the primary key (BL_cluster_id becomes the old key)."""
    path = TABLES / "step30h_landscape_per_BL.tsv"
    if not path.exists():
        return
    df = pd.read_csv(path, sep="\t", encoding="utf-8-sig")
    df = df.merge(bl_def[["BL_cluster_id", "BL", "BL_rank"]],
                   on="BL_cluster_id", how="left")
    front = ["BL", "BL_rank", "BL_cluster_id"]
    rest = [c for c in df.columns if c not in front]
    df = df[front + rest].sort_values("BL_rank").reset_index(drop=True)
    df.to_csv(path, sep="\t", index=False)
    print(f"[step30h-III] Rewrote step30h_landscape_per_BL.tsv with new BL names")


# ---------------------------------------------------------------------------
# 4. Driver
# ---------------------------------------------------------------------------
def main() -> None:
    bl_land = pd.read_csv(TABLES / "step30h_landscape_per_BL.tsv",
                           sep="\t", encoding="utf-8-sig")
    cl      = pd.read_csv(TABLES / "step30h_cluster_assignments.tsv",
                           sep="\t", encoding="utf-8-sig")

    # ---- Define BLs under the area DESC → connectivity DESC rule ----
    bl_def = define_bl(bl_land)
    print("[step30h-III] BL definition (area DESC → connectivity DESC):")
    for _, r in bl_def.iterrows():
        print(f"  {r['BL']}  = Ward cluster {int(r['BL_cluster_id'])}  "
              f"({int(r['n_populations']):>2} pops, "
              f"{r['convex_hull_area_km2']:>6.1f} km², "
              f"frac≤5km = {r['frac_within_BL_pairs_le_5000m']:.2f}, "
              f"total N_fert = {int(r['total_N_fert'])})")

    # ---- Renumber populations ----
    crosswalk = renumber_populations(cl, bl_def)
    crosswalk.to_csv(TABLES / "step30h_populationID_crosswalk.tsv",
                      sep="\t", index=False)
    print()
    print(f"[step30h-III] Wrote step30h_populationID_crosswalk.tsv "
          f"({len(crosswalk)} populations)")

    # ---- Write BL definition TSV ----
    bl_def.drop(columns=["sort_connectivity"]).to_csv(
        TABLES / "step30h_bl_definition.tsv", sep="\t", index=False,
    )
    print(f"[step30h-III] Wrote step30h_bl_definition.tsv")

    # ---- Backup + remap every population-keyed TSV ----
    backup_tsvs()
    remap_population_tsvs(crosswalk)
    rewrite_landscape_per_BL(bl_def)

    # ---- Headline: new population → BL mapping summary ----
    print()
    print("[step30h-III] New populationID assignment (first 5 of each BL):")
    merged = crosswalk.sort_values("populationID_new").head(50)
    for bl in bl_def["BL"]:
        sub = merged[merged["BL"] == bl].head(5)
        print(f"  {bl}:")
        for _, r in sub.iterrows():
            print(f"    P{int(r['populationID_new']):>2} "
                  f"← old P{int(r['populationID_old']):>2}  "
                  f"(lat {r['lat']:.4f}, lon {r['lon']:.4f})")


if __name__ == "__main__":
    main()
