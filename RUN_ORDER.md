# Phase 5 pipeline — running order

Phase 5 has one strict linear dependency chain (Steps 28 → 29 → 30) with
three branches that plug in at defined points (Steps 29a, 29b, 29c, 30b,
30c). Every script writes its outputs into `Tables/Phase5/` and
`Figures/Phase5/` and reads its inputs from the same folders + the
LEPA SQL database.

## The dependency graph

```
LEPA_SQL.db
    │
    ▼
step28  ─────►  Tables/Phase5/step28_events_spatial_neighborhood.tsv
                Tables/Phase5/step28_seed_sampling_per_mother.tsv
                Tables/Phase5/step28_mothers_for_full_detection_by_location.tsv
                Tables/Phase5/step28_coverage_curves_by_Nfertile.tsv
                Tables/Phase5/step29_sampling_per_location.tsv  (seed)
    │
    ▼
step29  ─────►  Tables/Phase5/step29_sampling_per_event.tsv
    │           Tables/Phase5/step29_sampling_per_location.tsv  (extends)
    │           Tables/Phase5/step29_field_team_sampling_recipe.tsv
    │           Tables/Phase5/step29_location_coverage_curves.tsv
    │
    ├─► step29a (sensitivity — one-time radius validation)
    │       reads  step28 spatial, step29 sampling, step26i P1
    │       writes step30_A_radius_sensitivity_*.tsv, .png
    │
    ├─► step29b (within-location connectivity)
    │       reads  step28 spatial, step29 sampling
    │       writes step29_location_connectivity.tsv,
    │              step29_location_connectivity_*.png
    │
    └─► step29c (fragmentation-aware sampling)
            reads  step28 spatial, step28 mothers-for-detection
            writes step29c_sampling_frag_aware_per_event.tsv,
                   step29c_sampling_comparison_per_location.tsv,
                   step29c_sampling_comparison.png
                (this is the authoritative field-team recipe)

step29b + step29 → feed step30 + step30b
    │
    ▼
step30   ────► All step30_A_prediction_* tables and figures
    (Phase A per-location predictions: diversity, P_compat,
     fragmentation-aware N_fertile_effective, cross-plot)

    step30b (fragmentation-index scatter, purely spatial)
        reads  step28 spatial, step29b connectivity
        writes step30_A_fragmentation_per_event.tsv,
               step30_A_fragmentation_per_location.tsv,
               step30_A_fragmentation_index.png

    step30c (§ C.0 empirical validation)
        reads  Tables/Phase4/step23_individual_allele_genotypes_with_nulls.tsv,
               step26i_L1_carrier_inventory.tsv,
               step30_A_prediction_prior_frequencies.tsv
        writes step30_C_pcompat_validation_at_eo.tsv,
               step30_C_pcompat_observed_vs_predicted.png
```

## The canonical run order

Run these in order for a complete Phase A rebuild:

```bash
python step28_seed_sampling_per_mother.py       --year 2025
python step29_event_location_sampling.py         --year 2025
python step29b_location_connectivity.py          --year 2025
python step29c_fragmentation_aware_sampling.py                # no --year flag
python step30_srk_diversity_prediction_vs_observed.py --year 2025
python step30b_fragmentation_index.py                          # no --year flag
python step30c_srk_validation_at_eo_level.py                   # no --year flag
```

## Optional / one-time

```bash
# Rerun only if the 50 m primary pollinator radius needs revisiting.
python step29a_pollinator_radius_sensitivity.py               # no --year flag
```

`step29a` is a one-time validation of the 50 m primary pollinator radius
adopted by every downstream script. Its output justifies the choice; it
does not need to be rerun on every pipeline iteration.

## Phase B (once seed genotypes are in hand)

`step30` takes real seed data via CLI flags:

```bash
python step30_srk_diversity_prediction_vs_observed.py \
    --year 2025 \
    --seed-genotypes  real_seed_genotypes.tsv \
    --mother-genotypes real_mother_genotypes.tsv \
    --match-seed-count
```

This adds the `step30_B_*` output family (mate-limitation regression +
Phase B comparison tables and figures).

## Notes on `--year`

Only Steps 28, 29, 29b, and 30 carry a `--year` flag — the year is used
to filter events / germplasm records before analysis. Steps 29a, 29c,
30b, and 30c are year-agnostic (they consume the year-filtered tables
that the year-aware scripts produced).
