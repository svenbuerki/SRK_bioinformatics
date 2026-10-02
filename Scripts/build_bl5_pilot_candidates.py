"""Build the BL5 pilot-candidate TSV for Part C preliminary analysis.

For each of three within-BL5 opposite-pair options (A1, A2, A3) the
TSV lists the two candidate locations with their mating-pool
structure (from step29d_mating_pool_summary.tsv), the mother count
already in the LEPA DB (from step29c_sampling_comparison_per_location.tsv),
and the design rationale. A single `recommended` column flags the
option Phase 5 defaults to (user-selected A2 on 2026-10-03).

Output: Tables/Phase5/step29c_partC_BL5_pilot_candidates.tsv
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from srk_bl_constants import make_location_label, base_eo

TABLES = Path("Tables/Phase5")

OPTIONS = [
    # (option, recommended, role, locationCode, locationID, rationale)
    ("A1", False, "predicted_robust",
     "EO32", 6,
     "Large and well-populated: 466 adults across 5 mating pools, "
     "all above the N=8 coupon-collector floor. Expected sustainable "
     "pollen compatibility, close to species mean."),
    ("A1", False, "predicted_drift_sensitive",
     "EO25-B", 21,
     "Small and drift-collapsed: 10 adults across 2 mating pools, "
     "both pools at or below the N=8 coupon-collector floor. "
     "Expected sustainable mean but wide credible interval entering "
     "the struggling band."),
    ("A2", True, "predicted_robust",
     "EO48", 7,
     "Medium and fully connected: 98 adults in 1 mating pool, well "
     "above N=8. The purest 'drift-safe' BL5 location — tests whether "
     "the model's sustainable prediction holds when fragmentation is "
     "absent."),
    ("A2", True, "predicted_drift_sensitive",
     "EO18-7", 19,
     "Medium and fragmented: 34 adults across 3 mating pools. Similar "
     "order-of-magnitude total census as the robust partner but very "
     "different pool structure — isolates the fragmentation × drift "
     "effect from the raw-size effect. **User-recommended pair "
     "(2026-10-03).**"),
    ("A3", False, "predicted_robust",
     "EO18-7", 17,
     "Largest fragmented BL5 location: 242 adults across 5 mating "
     "pools. Big enough that per-pool drift should be mild, but "
     "fragmentation is in play."),
    ("A3", False, "predicted_drift_sensitive",
     "EO24-7", 25,
     "Near-extinction floor: 3 adults in 1 mating pool. Maximum "
     "contrast within BL5 but low genotyping statistical power "
     "(n ≤ 3 adults means every seedling observation has large CI)."),
]


def main() -> None:
    pools = pd.read_csv(TABLES / "step29d_mating_pool_summary.tsv",
                         sep="\t", encoding="utf-8-sig")
    mothers = pd.read_csv(TABLES / "step29c_sampling_comparison_per_location.tsv",
                          sep="\t", encoding="utf-8-sig")
    mothers_by_loc = mothers.set_index("locationID")[
        ["M_frag_aware", "n_mothers_available_in_DB"]
    ]

    rows = []
    for option, recommended, role, code, loc_id, rationale in OPTIONS:
        p = pools[(pools["locationCode"] == code)
                   & (pools["locationID"] == loc_id)]
        if len(p) == 0:
            raise SystemExit(f"[bl5-candidates] {code} loc {loc_id} not in "
                              f"step29d_mating_pool_summary.tsv")
        p = p.iloc[0]
        m = (mothers_by_loc.loc[loc_id]
             if loc_id in mothers_by_loc.index
             else pd.Series(
                 {"M_frag_aware": pd.NA,
                  "n_mothers_available_in_DB": pd.NA}))
        rows.append({
            "option":                        option,
            "recommended":                   recommended,
            "role":                          role,
            "EOID":                          base_eo(code),
            "locationCode":                  code,
            "locationID":                    loc_id,
            "display_label":                 make_location_label(code, loc_id),
            "BL":                            p["BL"],
            "n_mating_pools_50m":            int(p["n_mating_pools_50m"]),
            "total_adults":                  int(p["total_adults"]),
            "largest_pool_N":                int(p["largest_pool_N"]),
            "smallest_pool_N":               int(p["smallest_pool_N"]),
            "median_pool_N":                 float(p["median_pool_N"]),
            "n_pools_below_coupon_floor_8":  int(p["n_pools_below_coupon_floor_8"]),
            "n_pools_below_SI_floor_1":      int(p["n_pools_below_SI_floor_1"]),
            "M_frag_aware_target":           (int(m["M_frag_aware"])
                                               if pd.notna(m["M_frag_aware"])
                                               else pd.NA),
            "n_mothers_available_in_DB":     (int(m["n_mothers_available_in_DB"])
                                               if pd.notna(m["n_mothers_available_in_DB"])
                                               else pd.NA),
            "design_rationale":              rationale,
        })

    df = pd.DataFrame(rows)
    out = TABLES / "step29c_partC_BL5_pilot_candidates.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"[bl5-candidates] Wrote {out}  ({len(df)} rows, "
          f"{df['option'].nunique()} options, "
          f"recommended = {df[df['recommended']]['display_label'].tolist()})")

    # Short summary
    print()
    for option, sub in df.groupby("option"):
        tag = "  ★" if sub["recommended"].any() else ""
        robust = sub[sub["role"] == "predicted_robust"].iloc[0]
        drift  = sub[sub["role"] == "predicted_drift_sensitive"].iloc[0]
        print(f"  {option}{tag}")
        print(f"    robust      : {robust['display_label']}  "
              f"({robust['n_mating_pools_50m']} pools, "
              f"{robust['total_adults']} adults, "
              f"{robust['n_mothers_available_in_DB']} mothers in DB)")
        print(f"    drift-sens. : {drift['display_label']}  "
              f"({drift['n_mating_pools_50m']} pools, "
              f"{drift['total_adults']} adults, "
              f"{drift['n_mothers_available_in_DB']} mothers in DB)")


if __name__ == "__main__":
    main()
