#!/usr/bin/env python3
"""Step 30 — SRK diversity prediction vs observed (two-generation).

Given the preliminary-study 32-Fg carrier inventory as an *a priori*
expectation, predict per-location SRK diversity and per-mother fecundation
failure under a random-mating null. If seed genotypes are supplied (or with
`--demo`, simulated from the prior), the script also produces the
*comparison* outputs — the two-generation test in which each seed lot
contributes one maternal allele (mother's own genotype) and one paternal
allele (pollen donor). The two spectra are then compared to prediction.

Output families are kept separate:

    Tables/Phase5/step30_A_prediction_*      ← Phase A (prior only, no seed data)
    Tables/Phase5/step30_B_comparison_*      ← Phase B (needs observed seed data)
    Tables/Phase5/step30_B_mate_limitation_* ← Phase B (needs seed + mother genotypes)
    Tables/Phase5/step30_B_si_escape_*       ← Phase B (needs seed + mother genotypes)
    figures/Phase5/step30_A_*                ← Phase A figures
    figures/Phase5/step30_B_*                ← Phase B figures

Inputs (read-only)
    /Users/sven/Documents/Current_projects/LEPA_fieldwork_protocol/SQL_DB/LEPA_SQL.db
    Tables/Phase5/step26i_L1_carrier_inventory.tsv  (prior — species-wide Fg counts)
    Tables/Phase5/step29_sampling_per_location.tsv  (locations + M recommended)

CLI
    python step30_srk_diversity_prediction_vs_observed.py
        Runs prediction only (writes the step30_A_prediction_* family).
    python step30_srk_diversity_prediction_vs_observed.py --demo
        Additionally simulates seed genotypes from the P1 prior and writes
        the step30_B_comparison_* family — end-to-end pipeline demo.
    python step30_srk_diversity_prediction_vs_observed.py \
        --seed-genotypes path/to/seeds.tsv
        Reads a real seed-genotype TSV (columns: locationID, eventID,
        germplasmID, seed_id, maternal_Fg, paternal_Fg) and writes the
        comparison family from it.
"""
from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
DEFAULT_DB = Path(
    "/Users/sven/Documents/Current_projects/LEPA_fieldwork_protocol/"
    "SQL_DB/LEPA_SQL.db"
)
DEFAULT_TABLES  = Path("Tables/Phase5")
DEFAULT_FIGURES = Path("figures/Phase5")
DEFAULT_PRIOR_TSV       = DEFAULT_TABLES / "step26i_L1_carrier_inventory.tsv"
DEFAULT_LOCATIONS_TSV   = DEFAULT_TABLES / "step29_sampling_per_location.tsv"
DEFAULT_CONNECTIVITY_TSV = DEFAULT_TABLES / "step29_location_connectivity.tsv"
DEFAULT_EVENT_COMPONENT_TSV = DEFAULT_TABLES / "step29c_event_to_component_50m.tsv"

PRIOR_ESS = 20.0        # effective sample size of the P1 prior; small so real
                        # seed data can shift it after modest sampling.
N_POSTERIOR_DRAWS = 4_000
RNG_SEED = 2030


# ---------------------------------------------------------------------------
# 1 — build the P1 (species-wide) prior over Fgs
# ---------------------------------------------------------------------------
def build_p1_prior(carrier_tsv: Path) -> pd.DataFrame:
    """Aggregate per-allele carrier counts to per-Fg frequency estimates.

    Returns a DataFrame indexed by Fg with columns:
        n_carriers  — total individuals carrying any allele of that Fg
        f_mean      — empirical mean frequency (normalised)
        alpha       — Dirichlet alpha (f_mean * PRIOR_ESS + a small smoother)
    """
    df = pd.read_csv(carrier_tsv, sep="\t", encoding="utf-8-sig")
    fg = (df.groupby("Fg", as_index=False)
            .agg(n_carriers=("n_carriers", "sum")))
    total = fg["n_carriers"].sum()
    fg["f_mean"] = fg["n_carriers"] / total
    fg["alpha"] = fg["f_mean"] * PRIOR_ESS + 0.5   # 0.5 = weak smoother
    fg = fg.sort_values("f_mean", ascending=False).reset_index(drop=True)
    return fg


# ---------------------------------------------------------------------------
# 2 — Prediction — from prior alone
# ---------------------------------------------------------------------------
def draw_frequencies(prior: pd.DataFrame, n_draws: int,
                     rng: np.random.Generator) -> np.ndarray:
    """Draw n_draws Dirichlet samples of Fg frequency vectors under the prior."""
    return rng.dirichlet(prior["alpha"].values, size=n_draws)


def _largest_remainder(float_counts: np.ndarray, total: int) -> np.ndarray:
    """Distribute `total` integer units across bins proportional to
    `float_counts`, using the largest-remainder method so Σ = total
    exactly (unlike `np.round`, which can drift)."""
    if total <= 0 or float_counts.sum() <= 0:
        return np.zeros(len(float_counts), dtype=int)
    weights = float_counts / float_counts.sum()
    exact = weights * total
    base = np.floor(exact).astype(int)
    remainder = total - base.sum()
    if remainder > 0:
        order = np.argsort(-(exact - base))
        base[order[:remainder]] += 1
    return base


def predicted_diversity_per_location(
        locations: pd.DataFrame,
        prior: pd.DataFrame,
        rng: np.random.Generator,
        seeds_per_mother: int = 15,
        n_local_replicates: int = 1000,
        event_component_tsv: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-50 m-component SRK diversity prediction aggregated to the
    location level by **set union** (Fgs are a set, not a continuous
    quantity, so simple size-weighted means do not apply).

    Rationale (symmetry with per-component P_compat)
    ------------------------------------------------
    The drift unit is a 50 m connected component inside a location,
    not the whole location: plants in different components do not share
    pollen, so their drift histories are independent. Phase A therefore
    simulates each component's SRK pool independently and derives the
    location-level numbers by **union** across components (what Fgs are
    present *somewhere* in the location vs. what the sampling actually
    recovers *somewhere* in the location).

    Three deliverables per replicate `k`:

      1. **Species-wide coverage** (unchanged — not drift-dependent).
         E[distinct Fgs observed | PLOIDY · M adult allele draws from
         the P1 species-wide prior].

      2. **Location-local pool size** (per-component, unioned).
         For each component *c*: draw `4 × component_N_fertile_c`
         alleles from P1, record its set of present Fgs. Location pool
         size = |union of component Fg sets|.

      3. **Location-local coverage** (per-component sampling, unioned).
         Distribute the location's mothers and seeds across components
         proportional to `component_N_fertile_c` (largest-remainder
         integer allocation). For each component: `A_delivered_c =
         PLOIDY · M_c + PATERNAL_ALLELES_PER_SEED · seeds_c` allele
         draws, record which Fgs this stochastic sample hits. Location
         detected = |union of component detected sets|; coverage =
         detected / pool_size.

    Returns
    -------
    (location_df, component_df)
        location_df : one row per location — species-wide coverage,
            union-based local pool size and coverage.
        component_df : one row per (locationID, component_id_50m) —
            per-component pool size, coverage, and the share of
            mothers/seeds the proportional allocation gave it.
    """
    from step28_seed_sampling_per_mother import PLOIDY, PATERNAL_ALLELES_PER_SEED

    freqs = draw_frequencies(prior, N_POSTERIOR_DRAWS, rng)   # (draws, K_fg)
    f_mean = prior["f_mean"].values
    K_fg = freqs.shape[1]

    # Load component inventory from step29c (same source P_compat uses).
    if event_component_tsv is None:
        event_component_tsv = DEFAULT_EVENT_COMPONENT_TSV
    comp_by_loc: dict[int, list[tuple[int, int]]] = {}
    if Path(event_component_tsv).exists():
        comp_df = pd.read_csv(event_component_tsv, sep="\t",
                              encoding="utf-8-sig")
        comp_unique = (
            comp_df.drop_duplicates(["locationID", "component_id_50m"])
                   [["locationID", "component_id_50m", "component_N_fertile"]]
        )
        for _, r in comp_unique.iterrows():
            comp_by_loc.setdefault(int(r["locationID"]), []).append(
                (int(r["component_id_50m"]), int(r["component_N_fertile"]))
            )

    loc_rows: list[dict] = []
    comp_rows: list[dict] = []
    has_actual = "M_actual_in_step28" in locations.columns

    for _, row in locations.iterrows():
        loc_id = int(row["locationID"])
        loc_code = row["locationCode"]
        M_ceiling = int(row["M_achievable_location"]) \
            if "M_achievable_location" in row else 0
        if has_actual and pd.notna(row.get("M_actual_in_step28")):
            M_used = int(row["M_actual_in_step28"])
        else:
            M_used = M_ceiling
        M_used = max(M_used, 0)

        realised_seeds = row.get("total_n_seeds_realised_exp")
        total_seeds = (int(realised_seeds)
                       if pd.notna(realised_seeds) and realised_seeds > 0
                       else 0)

        # Component inventory — fall back to a single component equal
        # to the whole location if the lookup is missing.
        comps = comp_by_loc.get(loc_id)
        if not comps:
            fallback_N = row.get("N_fertile_effective_50m",
                                  row.get("total_n_fertile", 1))
            try:
                fallback_N = int(fallback_N) if pd.notna(fallback_N) else 1
            except (TypeError, ValueError):
                fallback_N = 1
            comps = [(0, max(fallback_N, max(M_used, 1)))]

        comp_sizes = np.array([max(N, 1) for _, N in comps], dtype=float)
        N_fert_eff = int(comp_sizes.sum())

        # Allocate mothers and realised seeds across components
        # proportional to component_N_fertile (largest remainder keeps
        # the per-component integer total consistent with the whole-
        # location M_used and total_seeds).
        M_c_each = _largest_remainder(comp_sizes, M_used)
        seeds_c_each = _largest_remainder(comp_sizes, total_seeds)

        # Design A_delivered per component (unused in the prediction,
        # kept as a per-component column for reference).
        A_design_each = M_c_each * (PLOIDY
                                     + PATERNAL_ALLELES_PER_SEED * seeds_per_mother)
        A_actual_each = M_c_each * PLOIDY + PATERNAL_ALLELES_PER_SEED * seeds_c_each

        # --- 1. Species-wide coverage (tetraploid: PLOIDY · M draws) ---
        alleles_drawn = PLOIDY * max(M_used, 1)
        exp_distinct = (1.0 - (1.0 - freqs) ** alleles_drawn).sum(axis=1)

        # --- 2 + 3. Per-component drift + sampling, union at location ---
        n_comps = len(comps)
        comp_pool_mat = np.zeros((n_comps, n_local_replicates))
        comp_detected_mat = np.zeros((n_comps, n_local_replicates))
        comp_coverage_mat = np.zeros((n_comps, n_local_replicates))

        loc_pool_sizes = np.empty(n_local_replicates)
        loc_detected = np.empty(n_local_replicates)
        loc_coverages = np.empty(n_local_replicates)

        for k in range(n_local_replicates):
            present_at_loc = np.zeros(K_fg, dtype=bool)
            detected_at_loc = np.zeros(K_fg, dtype=bool)
            for c_idx, (comp_id, N_c) in enumerate(comps):
                pool_size_alleles = PLOIDY * max(int(N_c), 1)
                local = rng.choice(K_fg, size=pool_size_alleles, p=f_mean)
                counts = np.bincount(local, minlength=K_fg)
                present_c = counts > 0
                pool_size_c = int(present_c.sum())
                A_c = int(A_actual_each[c_idx])
                if A_c > 0 and counts.sum() > 0:
                    f_local = counts / counts.sum()
                    sampled = rng.choice(K_fg, size=A_c, p=f_local)
                    detected_c = np.zeros(K_fg, dtype=bool)
                    detected_c[np.unique(sampled)] = True
                    # (sampled Fgs are necessarily in present_c)
                else:
                    detected_c = np.zeros(K_fg, dtype=bool)

                present_at_loc |= present_c
                detected_at_loc |= detected_c

                comp_pool_mat[c_idx, k] = pool_size_c
                comp_detected_mat[c_idx, k] = int(detected_c.sum())
                comp_coverage_mat[c_idx, k] = (
                    detected_c.sum() / pool_size_c if pool_size_c > 0 else 0.0
                )

            loc_pool_sizes[k] = int(present_at_loc.sum())
            loc_detected[k] = int(detected_at_loc.sum())
            loc_coverages[k] = (
                detected_at_loc.sum() / present_at_loc.sum()
                if present_at_loc.any() else 0.0
            )

        # Emit per-component rows
        for c_idx, (comp_id, N_c) in enumerate(comps):
            pool_arr = comp_pool_mat[c_idx]
            cov_arr = comp_coverage_mat[c_idx]
            det_arr = comp_detected_mat[c_idx]
            comp_rows.append({
                "locationID":                 loc_id,
                "locationCode":               loc_code,
                "component_id_50m":           int(comp_id),
                "component_N_fertile":        int(N_c),
                "component_M_allocated":      int(M_c_each[c_idx]),
                "component_seeds_allocated":  int(seeds_c_each[c_idx]),
                "component_A_delivered_actual": int(A_actual_each[c_idx]),
                "component_A_delivered_design_at_15seeds":
                                              int(A_design_each[c_idx]),
                "component_pool_size_mean":   float(pool_arr.mean()),
                "component_pool_size_lo":     float(np.quantile(pool_arr, 0.025)),
                "component_pool_size_hi":     float(np.quantile(pool_arr, 0.975)),
                "component_Fgs_detected_mean": float(det_arr.mean()),
                "component_Fgs_detected_lo":   float(np.quantile(det_arr, 0.025)),
                "component_Fgs_detected_hi":   float(np.quantile(det_arr, 0.975)),
                "component_coverage_mean":    float(cov_arr.mean()),
                "component_coverage_lo":      float(np.quantile(cov_arr, 0.025)),
                "component_coverage_hi":      float(np.quantile(cov_arr, 0.975)),
            })

        A_delivered_design = M_used * (PLOIDY
                                        + PATERNAL_ALLELES_PER_SEED * seeds_per_mother)
        A_delivered_actual = PLOIDY * M_used + PATERNAL_ALLELES_PER_SEED * total_seeds

        loc_rows.append({
            "locationID":                    loc_id,
            "locationCode":                  loc_code,
            "M_mothers_in_db":               M_used,
            "M_achievable_ceiling":          M_ceiling,
            "N_fertile_effective":           N_fert_eff,
            "n_components_50m":              n_comps,
            "A_delivered_actual":            A_delivered_actual,
            "A_delivered_design_at_15seeds": A_delivered_design,
            "total_n_seeds_realised_exp":    total_seeds,
            # Species-wide (of the 32 P1 alleles):
            "predicted_distinct_Fgs":        float(np.mean(exp_distinct)),
            "predicted_distinct_Fgs_lo":     float(np.quantile(exp_distinct, 0.025)),
            "predicted_distinct_Fgs_hi":     float(np.quantile(exp_distinct, 0.975)),
            "predicted_species_coverage_mean": float(np.mean(exp_distinct) / K_fg),
            "n_fg_species_wide":             K_fg,
            # Location-local — UNBIASED (union across components):
            "predicted_local_pool_size_mean": float(loc_pool_sizes.mean()),
            "predicted_local_pool_size_lo":   float(np.quantile(loc_pool_sizes, 0.025)),
            "predicted_local_pool_size_hi":   float(np.quantile(loc_pool_sizes, 0.975)),
            # Location-local — SAMPLING (union across components):
            "predicted_local_Fgs_detected_from_sampling_mean": float(loc_detected.mean()),
            "predicted_local_Fgs_detected_from_sampling_lo":   float(np.quantile(loc_detected, 0.025)),
            "predicted_local_Fgs_detected_from_sampling_hi":   float(np.quantile(loc_detected, 0.975)),
            # Location-local — COVERAGE (detected union / present union):
            "predicted_local_coverage_mean": float(loc_coverages.mean()),
            "predicted_local_coverage_lo":   float(np.quantile(loc_coverages, 0.025)),
            "predicted_local_coverage_hi":   float(np.quantile(loc_coverages, 0.975)),
        })
    return pd.DataFrame(loc_rows), pd.DataFrame(comp_rows)


def predicted_diversity_matched_to_seeds(seeds: pd.DataFrame,
                                         prior: pd.DataFrame,
                                         rng: np.random.Generator) -> pd.DataFrame:
    """Per-location prediction whose sample size matches the actual number of
    allele observations in the seed data. Under tetraploid sporophytic each
    seed contributes PLOIDY // 2 + PATERNAL_ALLELES_PER_SEED = 2 + 2 = 4
    alleles (2 maternal from the mother's gamete + 2 paternal from the
    pollen donor's gamete). Use this in comparisons so the y-axis
    (observed) and x-axis (predicted) are on the same detection scale.
    """
    from step28_seed_sampling_per_mother import PLOIDY, PATERNAL_ALLELES_PER_SEED
    alleles_per_seed = PLOIDY // 2 + PATERNAL_ALLELES_PER_SEED    # 2 + 2 = 4
    freqs = draw_frequencies(prior, N_POSTERIOR_DRAWS, rng)
    rows = []
    for loc_id, sub in seeds.groupby("locationID"):
        n_alleles = alleles_per_seed * len(sub)   # tetraploid: 4 per seed
        exp_distinct = (1.0 - (1.0 - freqs) ** n_alleles).sum(axis=1)
        rows.append({
            "locationID":                        loc_id,
            "n_seed_alleles":                    int(n_alleles),
            "predicted_distinct_Fgs_matched":    float(np.mean(exp_distinct)),
            "predicted_distinct_Fgs_matched_lo": float(np.quantile(exp_distinct, 0.025)),
            "predicted_distinct_Fgs_matched_hi": float(np.quantile(exp_distinct, 0.975)),
        })
    return pd.DataFrame(rows)


def predicted_pcompat_distribution(prior: pd.DataFrame,
                                    n_mothers: int,
                                    rng: np.random.Generator) -> pd.DataFrame:
    """Species-wide reference distribution of per-mother P_compat under
    the **sporophytic tetraploid** Class I / II model (see
    `srk_si_model.py`). Simulates `n_mothers` random tetraploid mothers
    from the P1 prior — one row per mother with her genotype (4 alleles)
    and her analytical sporophytic P_compat.
    """
    from step28_seed_sampling_per_mother import PLOIDY
    from srk_si_model import (
        load_class_map, build_class_i_mask,
        load_zygosity_dist, sample_genotypes_empirical,
        p_compat_sporophytic_empirical,
    )
    f_mean = prior["f_mean"].values
    fg_ids = prior["Fg"].astype(str).values
    class_map = load_class_map()
    class_i_mask = build_class_i_mask(fg_ids.tolist(), class_map)
    zygosity_probs = load_zygosity_dist()
    # Sample mother tetraploid genotypes under EMPIRICAL LEPA zygosity.
    mothers = sample_genotypes_empirical(n_mothers, f_mean, zygosity_probs, rng)
    p_compat = p_compat_sporophytic_empirical(
        mothers, f_mean, class_i_mask, zygosity_probs,
        n_fathers=2_000, rng=rng)
    has_class_i = class_i_mask[mothers].any(axis=1)
    # Distinct-identity count per mother
    n_distinct = np.array([len(set(m.tolist())) for m in mothers])
    return pd.DataFrame({
        "sim_mother_id":  np.arange(n_mothers),
        **{f"mother_Fg_{i+1}": fg_ids[mothers[:, i]] for i in range(PLOIDY)},
        "n_distinct_functional_alleles": n_distinct,
        "expresses_class_I": has_class_i,
        "P_compat":       p_compat,
        "fecundation_failure_rate": 1.0 - p_compat,
    })


def predicted_pcompat_per_location(
        locations: pd.DataFrame,
        prior: pd.DataFrame,
        rng: np.random.Generator,
        n_draws: int = 400,
        n_sim_mothers: int = 200,
        event_component_tsv: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Phase-A **per-component** prediction of random-mating P_compat
    under the sporophytic tetraploid SI model with Class I / Class II
    dominance (see `srk_si_model.py` and § A.6 of the Phase 5 doc), then
    aggregated to the location level by component-size weighting.

    Rationale (per-component, not location-wide proxy)
    --------------------------------------------------
    Effective deme size is a property of each 50 m connected
    component inside a location, not of the location as a whole: plants
    inside a component share pollen; plants in a different component
    (same location but no pollen link) do not. Earlier revisions used
    `largest_component_share_50m` as a single "effective-N" proxy — this
    hid intra-location fragmentation. The current implementation
    simulates each component's local Fg pool independently and reports
    the location-level P_compat as the size-weighted mean across its
    components, matching the mating biology.

    Per component c, on each of `n_draws` simulation replicates:

      1. Simulate the component pool: draw PLOIDY · N_c = 4 · N_c SRK
         alleles i.i.d. from the species-wide P1 prior. Small components
         drift from P1 by ~1/√(4 N_c).
      2. Compute local Fg frequencies from that pool.
      3. Sample `n_sim_mothers` mothers under EMPIRICAL LEPA zygosity
         (66 % single-identity, 32 % 2-distinct, 2 % 3-distinct — see
         § A.6.3a) drawing identities from the local Fg frequencies.
      4. Per mother, compute sporophytic P_compat against
         empirically-zygotic candidate fathers drawn from the same
         local pool, under Class I / II dominance.
      5. Component mean = mean over `n_sim_mothers` sampled mothers.
      6. Component-level posterior CI across the `n_draws` replicates.
      7. Location mean on replicate k =
            Σ_c  P_compat_{c,k} · N_c / Σ_c N_c
         Posterior CI across the resulting n_draws location means.

    Returns
    -------
    (location_df, component_df)
        location_df : one row per location — the headline Phase-A
            prediction with the component-weighted mean and 95 % CI.
        component_df : one row per (locationID, component_id_50m) —
            exposes the finer-grained within-location variability so
            readers can see which components inside a fragmented
            location are predicted to struggle.
    """
    from step28_seed_sampling_per_mother import PLOIDY
    from srk_si_model import (
        load_class_map,
        build_class_i_mask,
        load_zygosity_dist,
        sample_genotypes_empirical,
        p_compat_sporophytic_empirical,
    )

    f_mean = prior["f_mean"].values
    K_fg = len(f_mean)
    fg_labels = prior["Fg"].astype(str).tolist() if "Fg" in prior.columns else \
                [f"FG{i+1:03d}" for i in range(K_fg)]
    class_map = load_class_map()
    class_i_mask = build_class_i_mask(fg_labels, class_map)
    zygosity_probs = load_zygosity_dist()

    # ------------------------------------------------------------------
    # Load the event → 50 m component lookup written by step29c; derive
    # a per-location component inventory (one row per unique
    # (locationID, component_id_50m) with component_N_fertile).
    # ------------------------------------------------------------------
    if event_component_tsv is None:
        event_component_tsv = DEFAULT_EVENT_COMPONENT_TSV
    comp_by_loc: dict[int, list[tuple[int, int]]] = {}
    if Path(event_component_tsv).exists():
        comp_df = pd.read_csv(event_component_tsv, sep="\t",
                              encoding="utf-8-sig")
        comp_unique = (
            comp_df.drop_duplicates(["locationID", "component_id_50m"])
                   [["locationID", "component_id_50m", "component_N_fertile"]]
        )
        for _, r in comp_unique.iterrows():
            comp_by_loc.setdefault(int(r["locationID"]), []).append(
                (int(r["component_id_50m"]), int(r["component_N_fertile"]))
            )
    else:
        print(f"[step30] {event_component_tsv} not found — falling back to "
              "location-level N_fertile_effective_50m proxy.")

    loc_rows: list[dict] = []
    comp_rows: list[dict] = []

    for _, row in locations.iterrows():
        loc_id = int(row["locationID"])
        loc_code = row["locationCode"]
        # Mothers actually in DB — permit-realistic sample (only used
        # for the location-level M label; prediction itself uses a
        # fixed n_sim_mothers so the prediction's precision does not
        # vary with observed sample size).
        raw = (row.get("M_actual_in_step28")
               or row.get("M_mothers_in_db")
               or row.get("M_achievable_location", 0))
        try:
            M_db = int(raw) if pd.notna(raw) else 0
        except (TypeError, ValueError):
            M_db = 0

        # Component inventory — fall back to a single component equal
        # to the whole location if the lookup is missing.
        comps = comp_by_loc.get(loc_id)
        if not comps:
            fallback_N = row.get("N_fertile_effective_50m",
                                  row.get("total_n_fertile", 1))
            try:
                fallback_N = int(fallback_N) if pd.notna(fallback_N) else 1
            except (TypeError, ValueError):
                fallback_N = 1
            comps = [(0, max(fallback_N, 1))]

        # Simulate each component independently.
        comp_sizes = np.array([max(N, 1) for _, N in comps], dtype=float)
        comp_weights = comp_sizes / comp_sizes.sum()
        # (n_components × n_draws)
        comp_mean_mat = np.zeros((len(comps), n_draws), dtype=float)
        for c_idx, (comp_id, N_c) in enumerate(comps):
            pool_size = PLOIDY * max(N_c, 1)
            for k in range(n_draws):
                local_alleles = rng.choice(K_fg, size=pool_size, p=f_mean)
                local_f = (np.bincount(local_alleles, minlength=K_fg)
                           / pool_size)
                if not (local_f > 0).any():
                    comp_mean_mat[c_idx, k] = 0.0
                    continue
                mother_genotypes = sample_genotypes_empirical(
                    n_sim_mothers, local_f, zygosity_probs, rng)
                pc = p_compat_sporophytic_empirical(
                    mother_genotypes, local_f, class_i_mask,
                    zygosity_probs, n_fathers=300, rng=rng)
                comp_mean_mat[c_idx, k] = pc.mean()

            arr = comp_mean_mat[c_idx]
            comp_rows.append({
                "locationID":                 loc_id,
                "locationCode":               loc_code,
                "component_id_50m":           int(comp_id),
                "component_N_fertile":        int(N_c),
                "component_weight":           float(comp_weights[c_idx]),
                "component_P_compat_mean":    float(arr.mean()),
                "component_P_compat_lo":      float(np.quantile(arr, 0.025)),
                "component_P_compat_hi":      float(np.quantile(arr, 0.975)),
            })

        # Size-weighted location mean per replicate.
        loc_means_per_draw = (comp_weights[:, None] * comp_mean_mat).sum(axis=0)
        N_fert_eff = int(comp_sizes.sum())
        loc_rows.append({
            "locationID":              loc_id,
            "locationCode":            loc_code,
            "N_fertile_effective":     N_fert_eff,
            "n_components_50m":        len(comps),
            "M_mothers_in_db":         M_db,
            "predicted_P_compat_mean": float(loc_means_per_draw.mean()),
            "predicted_P_compat_lo":   float(np.quantile(loc_means_per_draw, 0.025)),
            "predicted_P_compat_hi":   float(np.quantile(loc_means_per_draw, 0.975)),
        })

    return pd.DataFrame(loc_rows), pd.DataFrame(comp_rows)


# ---------------------------------------------------------------------------
# 3 — Seed genotype ingestion (real or simulated)
# ---------------------------------------------------------------------------
def simulate_seed_genotypes(locations: pd.DataFrame, prior: pd.DataFrame,
                            seeds_per_mother: int,
                            rng: np.random.Generator,
                            si_escape_rate: float = 0.0,
                            mean_ovules_per_mother: int = 100,
                            ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Demo-mode: simulate seed-DNA-level genotypes from the P1 prior,
    respecting each location's realistic sampled mothers (`M_actual_in_step28`
    when available; falls back to `M_achievable_location`).

    Also returns a per-mother table that carries the ground-truth maternal
    Fg pair, K_spatial (fake, sampled from a broad range), and a simulated
    observed seed count = Poisson(mean_ovules × P_compat) — used by the
    mate-limitation regression (§ 2.2.1). A user-supplied `si_escape_rate`
    (fraction of seeds where SI is bypassed) mimics partial-SI populations
    and lets us demonstrate the SI-escape permutation test (§ 2.2.2).
    """
    f_mean = prior["f_mean"].values
    fg_ids = prior["Fg"].values

    seed_rows: list[dict] = []
    mother_rows: list[dict] = []
    seed_counter = 0
    has_actual = "M_actual_in_step28" in locations.columns
    for _, row in locations.iterrows():
        if has_actual and pd.notna(row.get("M_actual_in_step28")):
            M = int(row["M_actual_in_step28"])
        else:
            M = int(row.get("M_achievable_location", 0))
        for m_idx in range(M):
            a = rng.choice(len(f_mean), p=f_mean)
            b = rng.choice(len(f_mean), p=f_mean)
            p_compat = 1.0 - f_mean[a] - (f_mean[b] if a != b else 0.0)
            p_compat = float(np.clip(p_compat, 0.01, 1.0))
            # simulate a K_spatial by drawing #compatible mates ∈ [1, 20];
            # correlates loosely with N_fertile at the location.
            K_spatial = int(2 * rng.integers(
                1, max(2, int(row["total_n_fertile"]) + 1)
            ))
            # observed seeds = Poisson(mean_ovules × P_compat)
            seeds_est = int(rng.poisson(mean_ovules_per_mother * p_compat))
            germplasm_id = f"sim_L{row['locationID']}_M{m_idx}"
            mother_rows.append({
                "locationID":     row["locationID"],
                "locationCode":   row["locationCode"],
                "germplasmID":    germplasm_id,
                "mother_Fg_a":    fg_ids[a],
                "mother_Fg_b":    fg_ids[b],
                "P_compat_true":  p_compat,
                "K_spatial":      K_spatial,
                "seeds_est":      seeds_est,
            })
            for s in range(seeds_per_mother):
                mat = a if rng.random() < 0.5 else b
                pat = rng.choice(len(f_mean), p=f_mean)
                # SI filter with escape probability
                if (pat == a or pat == b) and (rng.random() >= si_escape_rate):
                    continue
                seed_rows.append({
                    "locationID":   row["locationID"],
                    "locationCode": row["locationCode"],
                    "germplasmID":  germplasm_id,
                    "seed_id":      seed_counter,
                    "maternal_Fg":  fg_ids[mat],
                    "paternal_Fg":  fg_ids[pat],
                })
                seed_counter += 1
    return pd.DataFrame(seed_rows), pd.DataFrame(mother_rows)


# ---------------------------------------------------------------------------
# 4 — Comparison — posterior update + maternal vs paternal test
# ---------------------------------------------------------------------------
def location_posterior(seeds: pd.DataFrame, prior: pd.DataFrame,
                       rng: np.random.Generator) -> pd.DataFrame:
    """For each location, compute the posterior Dirichlet(alpha + counts) over
    Fgs using ALL observed allele draws (maternal + paternal) and return the
    posterior mean + 95 % CrI on the number of distinct Fgs implied by the
    sample.
    """
    n_fg = len(prior)
    fg_index = {fg: i for i, fg in enumerate(prior["Fg"])}
    alpha = prior["alpha"].values.copy()

    rows = []
    for loc_id, sub in seeds.groupby("locationID"):
        counts = np.zeros(n_fg)
        for col in ("maternal_Fg", "paternal_Fg"):
            vc = sub[col].value_counts()
            for fg, n in vc.items():
                if fg in fg_index:
                    counts[fg_index[fg]] += n
        posterior_alpha = alpha + counts
        draws = rng.dirichlet(posterior_alpha, size=1000)
        # Observed distinct Fgs in this location's sample
        observed = set(sub["maternal_Fg"]).union(set(sub["paternal_Fg"]))
        rows.append({
            "locationID":            loc_id,
            "n_seeds":               len(sub),
            "observed_distinct_Fgs": len(observed),
            "posterior_mean_distinct_Fgs":
                float(np.mean((draws > 0).sum(axis=1))),
            "posterior_distinct_Fgs_lo":
                float(np.quantile((draws > 0).sum(axis=1), 0.025)),
            "posterior_distinct_Fgs_hi":
                float(np.quantile((draws > 0).sum(axis=1), 0.975)),
        })
    return pd.DataFrame(rows)


def maternal_vs_paternal(seeds: pd.DataFrame) -> pd.DataFrame:
    """Per-location chi-square-style comparison of maternal vs paternal Fg
    frequency spectra."""
    rows = []
    for loc_id, sub in seeds.groupby("locationID"):
        mat = sub["maternal_Fg"].value_counts()
        pat = sub["paternal_Fg"].value_counts()
        fgs = sorted(set(mat.index).union(pat.index))
        m = np.array([mat.get(fg, 0) for fg in fgs], float)
        p = np.array([pat.get(fg, 0) for fg in fgs], float)
        expected_m = (m.sum() / (m.sum() + p.sum())) * (m + p)
        expected_p = (p.sum() / (m.sum() + p.sum())) * (m + p)
        with np.errstate(divide="ignore", invalid="ignore"):
            chi2 = np.nansum(
                (m - expected_m) ** 2 / np.where(expected_m > 0, expected_m, 1) +
                (p - expected_p) ** 2 / np.where(expected_p > 0, expected_p, 1)
            )
        rows.append({
            "locationID":       loc_id,
            "n_distinct_Fgs":   len(fgs),
            "n_maternal_alleles": int(m.sum()),
            "n_paternal_alleles": int(p.sum()),
            "chi2":             float(chi2),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 4b — Mate-limitation regression (Test 1, § 2.2.1 of the design doc)
# ---------------------------------------------------------------------------
def compute_location_posterior_frequencies(seeds: pd.DataFrame,
                                            prior: pd.DataFrame
                                            ) -> dict[int, np.ndarray]:
    """Per-location posterior mean Fg frequency vector, used to compute
    P_compat for each mother at her own location."""
    fg_index = {fg: i for i, fg in enumerate(prior["Fg"])}
    alpha = prior["alpha"].values.astype(float)
    out: dict[int, np.ndarray] = {}
    for loc_id, sub in seeds.groupby("locationID"):
        counts = np.zeros(len(prior))
        for col in ("maternal_Fg", "paternal_Fg"):
            for fg, n in sub[col].value_counts().items():
                if fg in fg_index:
                    counts[fg_index[fg]] += n
        posterior_alpha = alpha + counts
        out[loc_id] = posterior_alpha / posterior_alpha.sum()
    return out


def mate_limitation_regression(seeds: pd.DataFrame,
                                mothers: pd.DataFrame,
                                prior: pd.DataFrame,
                                connectivity: pd.DataFrame | None = None,
                                ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Fit the mate-limitation regression at the LOCATION level:

        mean_seeds_per_mother ~ mean_P_compat + mean_K_spatial + connectivity

    weighted by n_mothers per location. `connectivity` (optional) is the
    per-location `largest_component_share_50m` from Step 28c — the
    fraction of the location's adults connected via 10 m pollen flow.
    When provided, it enters the regression as an explicit predictor
    so we can separate three effects:

      β₁ (P_compat)     — allele-frequency drift → mate limitation
      β₂ (K_spatial)    — mother-level spatial neighbourhood size
      β₃ (connectivity) — within-location pollen-flow fragmentation

    Rationale (design doc § 2.2.1). Mate limitation is a *population*-level
    phenomenon: whether reduced pollen-donor diversity depresses seed set
    is a question about *locations*, not individual mothers. The per-mother
    variation within a location is heterogeneity we want to average over,
    not variance we want to model.

    Returns (per_mother_annotated, per_location_summary, coefficient_table).
    per_mother has predicted_seeds / residual columns at the individual
    level (useful for looking at outlier mothers). per_location is what the
    regression is fit on.
    """
    import statsmodels.api as sm

    fg_index = {fg: i for i, fg in enumerate(prior["Fg"])}
    loc_post = compute_location_posterior_frequencies(seeds, prior)

    per_mother = mothers.copy()
    p_compat = []
    for _, m in per_mother.iterrows():
        f = loc_post.get(m["locationID"])
        if f is None:
            p_compat.append(np.nan); continue
        a_i = fg_index.get(m["mother_Fg_a"]); b_i = fg_index.get(m["mother_Fg_b"])
        if a_i is None or b_i is None:
            p_compat.append(np.nan); continue
        if a_i == b_i:
            p_compat.append(1.0 - f[a_i])
        else:
            p_compat.append(1.0 - f[a_i] - f[b_i])
    per_mother["P_compat_posterior"] = np.clip(p_compat, 0.0, 1.0)

    per_mother = per_mother.dropna(subset=[
        "P_compat_posterior", "K_spatial", "seeds_est"]).copy()
    per_mother["seeds_est"] = per_mother["seeds_est"].astype(float)

    # ---------------- Location-level aggregation ----------------
    loc_code_col = ("locationCode"
                    if "locationCode" in per_mother.columns else None)
    agg_dict = {
        "n_mothers":         ("germplasmID", "count"),
        "mean_seeds":        ("seeds_est", "mean"),
        "std_seeds":         ("seeds_est", "std"),
        "mean_P_compat":     ("P_compat_posterior", "mean"),
        "std_P_compat":      ("P_compat_posterior", "std"),
        "mean_K_spatial":    ("K_spatial", "mean"),
    }
    per_location = per_mother.groupby("locationID", as_index=False).agg(**agg_dict)
    if loc_code_col:
        code_map = per_mother.drop_duplicates("locationID").set_index(
            "locationID")["locationCode"]
        per_location["locationCode"] = per_location["locationID"].map(code_map)

    # Merge in connectivity as an explicit predictor if available.
    have_connectivity = False
    if connectivity is not None and "largest_component_share_50m" in connectivity.columns:
        per_location = per_location.merge(
            connectivity[["locationID", "largest_component_share_50m"]],
            on="locationID", how="left",
        )
        per_location["largest_component_share_50m"] = (
            per_location["largest_component_share_50m"].fillna(1.0)
        )
        have_connectivity = True

    # Fit location-level OLS weighted by n_mothers
    X_dict = {
        "const":     1.0,
        "P_compat":  per_location["mean_P_compat"].values,
        "K_spatial": per_location["mean_K_spatial"].astype(float).values,
    }
    terms = ["P_compat", "K_spatial"]
    if have_connectivity:
        X_dict["connectivity"] = per_location["largest_component_share_50m"].astype(float).values
        terms.append("connectivity")
    X = pd.DataFrame(X_dict).astype(float)
    y = per_location["mean_seeds"].values.astype(float)
    w = per_location["n_mothers"].values.astype(float)

    fit = sm.WLS(y, X, weights=w).fit()
    per_location["predicted_mean_seeds"] = fit.predict(X)
    per_location["residual_mean_seeds"]  = y - per_location["predicted_mean_seeds"]

    # Mother-level predictions from the location-level model
    per_mother = per_mother.merge(
        per_location[["locationID", "predicted_mean_seeds"]],
        on="locationID", how="left",
    ).rename(columns={"predicted_mean_seeds": "predicted_seeds_from_location_model"})
    per_mother["residual_seeds"] = (per_mother["seeds_est"]
                                    - per_mother["predicted_seeds_from_location_model"])

    coefs = pd.DataFrame({
        "term":     terms,
        "estimate": [fit.params[t] for t in terms],
        "std_err":  [fit.bse[t] for t in terms],
        "ci_lo":    [fit.conf_int().loc[t, 0] for t in terms],
        "ci_hi":    [fit.conf_int().loc[t, 1] for t in terms],
        "p_value":  [fit.pvalues[t] for t in terms],
    })
    coefs["significant_95"] = (coefs["ci_lo"] > 0) | (coefs["ci_hi"] < 0)
    return per_mother, per_location, coefs


# ---------------------------------------------------------------------------
# 4c — SI-escape rate permutation test (Test 2, § 2.2.2 of the design doc)
# ---------------------------------------------------------------------------
def si_escape_permutation(seeds: pd.DataFrame,
                          mothers: pd.DataFrame,
                          p_null: float = 0.02,
                          rng: np.random.Generator | None = None
                          ) -> pd.DataFrame:
    """Per location, test observed rate of self-matching paternal alleles
    (paternal Fg ∈ mother's Fg pair) against the strict-SI null.

    Under **strict SI**, a self-matching paternal Fg is genetically
    impossible — the expected rate is 0. In practice we allow a small
    baseline `p_null` (default 0.02 = 2 %) to accommodate genotyping error
    and rare true escapes. The test is a one-sided **exact Binomial test**:

        H0 : pi_self_match  <=  p_null       (strict SI or near-strict)
        H1 : pi_self_match  >   p_null       (partial-SI escape signature)

    p-value = P( X >= n_self_matching  |  X ~ Binomial(n_seeds, p_null) ).

    Returns per-location table with n_seeds_scored, n_self_matching,
    pi_self_match, wilson_ci_lo / hi, binomial p-value and a
    Benjamini-Hochberg FDR q-value across locations.

    The permutation-style API of the earlier implementation is retained
    (the `rng` argument is accepted but no longer used, so existing callers
    don't break).
    """
    from scipy import stats

    _ = rng  # unused now; kept for backwards-compatible signature
    seeds_with_mother = seeds.merge(
        mothers[["germplasmID", "mother_Fg_a", "mother_Fg_b"]],
        on="germplasmID", how="left"
    ).dropna(subset=["mother_Fg_a", "mother_Fg_b"])

    rows = []
    for loc_id, sub in seeds_with_mother.groupby("locationID"):
        paternal = sub["paternal_Fg"].values
        mA = sub["mother_Fg_a"].values
        mB = sub["mother_Fg_b"].values
        match = (paternal == mA) | (paternal == mB)
        n_seeds = int(len(sub))
        n_self  = int(match.sum())
        pi      = n_self / n_seeds if n_seeds > 0 else np.nan

        # One-sided exact Binomial test against p_null
        # p-value = P(X >= n_self | Binomial(n_seeds, p_null))
        if n_seeds == 0:
            p_bin = np.nan
        else:
            p_bin = float(stats.binom.sf(n_self - 1, n_seeds, p_null))

        # Wilson score CI for observed rate (95 %)
        if n_seeds > 0:
            ci = stats.binomtest(n_self, n_seeds).proportion_ci(
                confidence_level=0.95, method="wilson")
            ci_lo, ci_hi = float(ci.low), float(ci.high)
        else:
            ci_lo, ci_hi = np.nan, np.nan

        rows.append({
            "locationID":       loc_id,
            "n_seeds_scored":   n_seeds,
            "n_self_matching":  n_self,
            "pi_self_match":    pi,
            "wilson_ci_lo":     ci_lo,
            "wilson_ci_hi":     ci_hi,
            "p_null_used":      p_null,
            "p_binomial":       p_bin,
        })
    df = pd.DataFrame(rows).sort_values("p_binomial").reset_index(drop=True)
    # Benjamini-Hochberg FDR
    m = len(df)
    if m > 0:
        ranks = np.arange(1, m + 1)
        q = df["p_binomial"].values * m / ranks
        q = np.minimum.accumulate(q[::-1])[::-1]
        df["q_bh_fdr"] = np.clip(q, 0.0, 1.0)
    return df


# ---------------------------------------------------------------------------
# 5 — figures
# ---------------------------------------------------------------------------
def _year_suffix(year: int | None) -> str:
    return f"  —  {year}" if year is not None else ""


def _demo_suffix(demo: bool) -> str:
    return "   [DEMO — synthetic data]" if demo else ""


def _add_demo_watermark(fig, demo: bool):
    """Diagonal 'DEMO' watermark across the figure — obvious enough that a
    Phase B figure produced under `--demo` cannot be mistaken for real
    analysis, faint enough not to obscure the plot."""
    if not demo:
        return
    fig.text(
        0.5, 0.5, "DEMO",
        ha="center", va="center",
        fontsize=90, color="#b2182b",
        alpha=0.10, rotation=25, fontweight="bold",
        zorder=100,
    )


def plot_diversity_unbiased_vs_sampling(pred: pd.DataFrame,
                                         prior: pd.DataFrame,
                                         out_png: Path, out_pdf: Path,
                                         year: int | None = None):
    """Three-column × BL-row figure separating the UNBIASED local Fg
    diversity per location (what the population holds, driven by
    N_fertile_effective alone) from the SAMPLING-inferred detection
    (what our seed genotyping will recover, driven by M × 15 seeds)
    and the coverage ratio.

    Layout matches every other Phase 5 per-location figure: rows =
    Bottleneck Lineages in canonical BL_ORDER (BL4 → BL5 → BL3 → BL1
    → BL2, top-to-bottom); columns = Panel A / B / C. Row heights
    ∝ number of locations per BL. Within each BL row, locations are
    sorted by unbiased pool size (ascending, small → large).

    Panel A — Unbiased local Fg diversity. Depends only on
              N_fertile_effective. Answers: how many SRK alleles does
              this location actually hold under drift?
    Panel B — Sampling-inferred detection. Depends on M × 15 seeds
              against the true local pool. Each mother's 15 seeds
              contribute 2·15 = 30 paternal-allele samples of the
              local pollen donor pool; the mother's own 4-copy
              genotype adds 4 maternal-allele samples. Total
              A_delivered = M × 34 draws.
    Panel C — Coverage fraction = Panel B / Panel A. Dotted target
              line at 0.90. Answers: how well does our sampling
              recover the truth?
    """
    from srk_bl_constants import BL_COLORS, BL_ORDER, locationCode_to_bl, location_label_series

    K_species = len(prior)
    df = pred.copy()
    df["BL"] = locationCode_to_bl(df["locationCode"]).values
    df["BL"] = df["BL"].fillna("Unassigned")
    m_col = "M_mothers_in_db" if "M_mothers_in_db" in df.columns \
        else "M_achievable_location"
    df["M"] = df[m_col].astype(float)
    df["N_eff"] = df["N_fertile_effective"].astype(float)
    # Number of 50 m demes at the location — the explicit link to
    # Figure 2's per-location deme-structure display.
    df["n_pools"] = df.get("n_components_50m",
                             pd.Series(1, index=df.index)).astype(int)
    # Total seeds recorded at each location in the LEPA DB — the raw
    # per-location input that drives Panel B (A_delivered = 4·M +
    # 2·total_seeds). Showing the total, not a per-mother mean, keeps
    # the arithmetic transparent and does not hide the (often skewed)
    # per-mother distribution behind an average.
    df["total_seeds"] = df["total_n_seeds_realised_exp"].astype(float)

    # BL row order — canonical BL_ORDER first, Unassigned last (if any).
    bls = [b for b in BL_ORDER if b in df["BL"].values]
    if (df["BL"] == "Unassigned").any():
        bls.append("Unassigned")
    panel_colours = {**BL_COLORS, "Unassigned": "#8a8a8a"}

    # Row heights ∝ N locations per BL (min 1 for empty rows).
    heights = [max(int((df["BL"] == b).sum()), 1) for b in bls]
    fig, axes = plt.subplots(
        len(bls), 3,
        figsize=(15.5, max(6.5, 0.30 * sum(heights) + 2.5)),
        gridspec_kw={
            "height_ratios": heights,
            "width_ratios":  [1.0, 1.0, 0.75],
        },
        sharex="col",
    )
    if len(bls) == 1:
        axes = np.array([axes])

    def _plot_metric(ax, sub, y, mean_col, lo_col, hi_col, colour):
        # Clip to non-negative — when the empirical distribution is
        # nearly constant (e.g. coverage = 1.0 for every replicate),
        # floating-point noise can push the quantile a hair above
        # the mean, which matplotlib otherwise rejects.
        xerr_lo = (sub[mean_col] - sub[lo_col]).clip(lower=0.0)
        xerr_hi = (sub[hi_col]   - sub[mean_col]).clip(lower=0.0)
        for i in range(len(sub)):
            ax.errorbar(
                sub[mean_col].iloc[i], y[i],
                xerr=[[xerr_lo.iloc[i]], [xerr_hi.iloc[i]]],
                fmt="none", ecolor=colour, alpha=0.4,
                elinewidth=1.0, capsize=2.0, zorder=1,
            )
        ax.scatter(
            sub[mean_col], y,
            s=55, c=colour, edgecolor="white", linewidth=0.6, zorder=2,
        )

    for row_idx, bl in enumerate(bls):
        sub = df[df["BL"] == bl].sort_values(
            "predicted_local_pool_size_mean", ascending=True
        ).reset_index(drop=True)
        y = np.arange(len(sub))
        colour = panel_colours[bl]
        axA_row, axB_row, axC_row = axes[row_idx]

        # ---- Panel A: unbiased ----
        _plot_metric(axA_row, sub, y,
                     "predicted_local_pool_size_mean",
                     "predicted_local_pool_size_lo",
                     "predicted_local_pool_size_hi", colour)
        axA_row.axvline(K_species, color="#333", ls="--", lw=1.0, alpha=0.4)
        axA_row.set_xlim(0, K_species + 2)
        axA_row.set_ylim(-0.7, len(sub) - 0.3)
        site_labels = location_label_series(sub)
        labels = [
            f"{label}  ({int(k)} deme{'s' if int(k) != 1 else ''}, "
            f"{int(n)} adults, "
            f"mothers sampled = {int(m)}, "
            f"seeds recorded = {int(s)})"
            for label, k, n, m, s in zip(
                site_labels, sub["n_pools"], sub["N_eff"],
                sub["M"], sub["total_seeds"],
            )
        ]
        axA_row.set_yticks(y)
        axA_row.set_yticklabels(labels, fontsize=8)
        axA_row.spines["top"].set_visible(False)
        axA_row.spines["right"].set_visible(False)

        # ---- Panel B: sampling ----
        _plot_metric(axB_row, sub, y,
                     "predicted_local_Fgs_detected_from_sampling_mean",
                     "predicted_local_Fgs_detected_from_sampling_lo",
                     "predicted_local_Fgs_detected_from_sampling_hi", colour)
        axB_row.axvline(K_species, color="#333", ls="--", lw=1.0, alpha=0.4)
        axB_row.set_xlim(0, K_species + 2)
        axB_row.set_ylim(-0.7, len(sub) - 0.3)
        axB_row.set_yticks(y); axB_row.set_yticklabels([])
        axB_row.spines["top"].set_visible(False)
        axB_row.spines["right"].set_visible(False)

        # ---- Panel C: coverage ----
        _plot_metric(axC_row, sub, y,
                     "predicted_local_coverage_mean",
                     "predicted_local_coverage_lo",
                     "predicted_local_coverage_hi", colour)
        axC_row.axvline(0.9, color="#333", ls=":", lw=1.0, alpha=0.6)
        axC_row.axvline(1.0, color="#333", ls="--", lw=1.0, alpha=0.4)
        axC_row.set_xlim(0, 1.05)
        axC_row.set_ylim(-0.7, len(sub) - 0.3)
        axC_row.set_yticks(y); axC_row.set_yticklabels([])
        axC_row.spines["top"].set_visible(False)
        axC_row.spines["right"].set_visible(False)

        # BL label on the right of Panel C
        axC_row.text(
            1.03, 0.5, bl,
            transform=axC_row.transAxes,
            fontsize=13, fontweight="bold", color=colour,
            va="center", ha="left",
        )

    # Column headers (only on the top row)
    axes[0, 0].set_title("A. What the population holds at each location",
                          fontsize=11)
    axes[0, 1].set_title("B. What our sampling will detect",
                          fontsize=11)
    axes[0, 2].set_title("C. How well the sampling recovers the truth",
                          fontsize=11)
    # Species-wide ceiling label on the top A/B panels
    axes[0, 0].text(K_species - 0.3, axes[0, 0].get_ylim()[1] - 0.2,
                    f"species-wide ceiling = {K_species}",
                    fontsize=8, color="#333", ha="right", va="top")

    # Column x-labels — short single lines, panel-bounded to avoid overlap.
    axes[-1, 0].set_xlabel(
        "SRK allele classes present at the location",
        fontsize=10)
    axes[-1, 1].set_xlabel(
        "SRK allele classes recovered by our sampling",
        fontsize=10)
    axes[-1, 2].set_xlabel(
        "Fraction of local SRK pool recovered",
        fontsize=10)

    fig.suptitle(
        "Predicted SRK allele diversity per LEPA location"
        f"{_year_suffix(year)}\n"
        "A: what the population holds — driven by the effective breeding unit at the location.  "
        "B: what our sampling detects — from the mothers and seeds already recorded in the LEPA DB.  "
        "C: coverage = B ÷ A (dotted line = 90 % target).\n"
        "Panelled by Bottleneck Lineage in BL_ORDER.",
        fontsize=11, y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 0.96, 0.97])
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)



def plot_prediction_fecundation(pcompat_per_loc: pd.DataFrame,
                                out_png: Path, out_pdf: Path,
                                year: int | None = None,
                                bands: dict[str, float] | None = None):
    """Per-location Phase-A prediction of random-mating pollen compatibility
    under the sporophytic tetraploid Class I / II model (§ A.6, Part 2),
    grouped by BL, with the failed / struggling / sustainable traffic-light
    bands recalibrated against the sporophytic species mean.
    """
    from srk_bl_constants import BL_COLORS, BL_ORDER, locationCode_to_bl, location_label_series

    df = pcompat_per_loc.copy()
    df["BL"] = locationCode_to_bl(df["locationCode"]).values
    df["BL"] = df["BL"].fillna("Unassigned")

    bls = [b for b in BL_ORDER if b in df["BL"].values]
    if (df["BL"] == "Unassigned").any():
        bls.append("Unassigned")
    panel_colours = {**BL_COLORS, "Unassigned": "#8a8a8a"}

    heights = [max(int((df["BL"] == b).sum()), 1) for b in bls]
    fig, axes = plt.subplots(
        len(bls), 1,
        figsize=(9.5, max(6.0, 0.28 * sum(heights) + 1.5)),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )
    if len(bls) == 1:
        axes = [axes]

    # Traffic-light thresholds — recalibrated against sporophytic species mean
    T_FAILED_HI     = bands["failed_max"]     if bands else 0.20
    T_STRUGGLING_HI = bands["struggling_max"] if bands else 0.40
    SPECIES_MEAN    = bands["species_mean"]   if bands else 0.63
    x_upper = max(0.30, min(1.0, 2.0 * SPECIES_MEAN))
    band_alpha = 0.12

    for ax, bl in zip(axes, bls):
        sub = df[df["BL"] == bl].sort_values(
            "predicted_P_compat_mean", ascending=True,
        ).reset_index(drop=True)
        y = np.arange(len(sub))

        # Traffic-light bands
        ax.axvspan(0.0, T_FAILED_HI,       color="#b2182b",
                    alpha=band_alpha, zorder=0)
        ax.axvspan(T_FAILED_HI, T_STRUGGLING_HI, color="#e08214",
                    alpha=band_alpha, zorder=0)
        ax.axvspan(T_STRUGGLING_HI, x_upper, color="#1b7837",
                    alpha=band_alpha, zorder=0)
        # Species-mean reference line
        ax.axvline(SPECIES_MEAN, color="#1b7837", ls=":", lw=1.0, alpha=0.7)

        xerr_lo = sub["predicted_P_compat_mean"] - sub["predicted_P_compat_lo"]
        xerr_hi = sub["predicted_P_compat_hi"] - sub["predicted_P_compat_mean"]
        colour = panel_colours[bl]
        ax.errorbar(
            sub["predicted_P_compat_mean"], y,
            xerr=[xerr_lo, xerr_hi],
            fmt="none", ecolor=colour, alpha=0.5,
            elinewidth=1.2, capsize=2.5, zorder=1,
        )
        sizes = 30 + 8 * np.sqrt(np.clip(sub["M_mothers_in_db"], 1, None))
        ax.scatter(
            sub["predicted_P_compat_mean"], y,
            s=sizes, c=colour, edgecolor="white",
            linewidth=0.6, zorder=2,
        )
        n_pools_col = sub.get(
            "n_components_50m",
            pd.Series(1, index=sub.index)).astype(int)
        site_labels = location_label_series(sub)
        labels = [
            f"{label}  ({int(k)} deme{'s' if int(k) != 1 else ''}, "
            f"{int(n)} adults, "
            f"mothers sampled = {int(m)})"
            for label, k, n, m in zip(
                site_labels,
                n_pools_col,
                sub["N_fertile_effective"],
                sub["M_mothers_in_db"],
            )
        ]
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlim(0.0, x_upper)
        ax.set_ylim(-0.7, len(sub) - 0.3)
        ax.text(1.01, 0.5, bl,
                transform=ax.transAxes,
                fontsize=13, fontweight="bold", color=colour,
                va="center", ha="left")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Traffic-light legend at the top of the first panel
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    tl_handles = [
        Patch(facecolor="#b2182b", alpha=0.35,
              label=f"failed  (< {T_FAILED_HI:.3f})"),
        Patch(facecolor="#e08214", alpha=0.35,
              label=f"struggling  ({T_FAILED_HI:.3f}–{T_STRUGGLING_HI:.3f})"),
        Patch(facecolor="#1b7837", alpha=0.35,
              label=f"sustainable  (≥ {T_STRUGGLING_HI:.3f})"),
        Line2D([0], [0], color="#1b7837", ls=":", lw=1.2,
               label=f"sporophytic species mean  ({SPECIES_MEAN:.3f})"),
    ]
    axes[0].legend(handles=tl_handles, loc="upper left",
                    fontsize=9, frameon=True)

    axes[-1].set_xlabel(
        "Predicted pollen compatibility under random mating  "
        "(mean across sampled mothers per location; 95 % credible interval)",
        fontsize=11,
    )
    fig.suptitle(
        f"Predicted per-mother pollen compatibility across LEPA"
        f"{_year_suffix(year)}",
        fontsize=13, y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 0.94, 0.97])
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)


def plot_diversity_vs_pcompat(pred_div: pd.DataFrame,
                               pred_pc: pd.DataFrame,
                               prior: pd.DataFrame,
                               out_png: Path, out_pdf: Path,
                               year: int | None = None,
                               bands: dict[str, float] | None = None):
    """Cross-plot of the two Phase-A per-location predictions:
      x = predicted number of distinct SRK alleles at that location
      y = predicted mean pollen compatibility under random mating
    One dot per location, coloured by BL. Error bars on both axes come
    from the same Dirichlet posterior draws used to build the two
    single-quantity figures. The dashed reference line shows the
    theoretical relation E[P_compat] = 1 - 2/k_eff (with k_eff ranging
    across the panel's x-axis), which is the analytic tie between SRK
    diversity and mate compatibility under random mating; individual
    locations sit on or near this curve depending on how uneven their
    Fg-frequency distribution is.
    """
    from srk_bl_constants import BL_COLORS, BL_ORDER, locationCode_to_bl, location_label_series

    m = pred_div.merge(
        pred_pc[["locationID", "predicted_P_compat_mean",
                 "predicted_P_compat_lo", "predicted_P_compat_hi",
                 "M_mothers_in_db"]],
        on="locationID", how="inner", suffixes=("", "_pc"),
    )
    m["BL"] = locationCode_to_bl(m["locationCode"]).values
    m["BL"] = m["BL"].fillna("Unassigned")

    fig, ax = plt.subplots(figsize=(10.0, 7.0))

    # Theoretical curve: E[P_compat] = 1 - 2/k_eff, with k_eff varying
    # from 1 to the species-wide pool size (32). We plot in the observed-
    # coverage domain by using the predicted number of distinct alleles
    # as a proxy for effective k, so the curve is qualitative — a visual
    # reminder that mate compatibility rises with SRK diversity.
    ks = np.arange(1, len(prior) + 1)
    # Species-mean P_compat reference (sporophytic Part 2)
    SPECIES_MEAN = bands["species_mean"] if bands else 0.63
    ax.axhline(SPECIES_MEAN, color="#1b7837", ls=":", lw=1.2, alpha=0.7,
               label=f"sporophytic species-mean P_compat  ({SPECIES_MEAN:.3f})")

    bls = [b for b in BL_ORDER if b in m["BL"].values]
    if (m["BL"] == "Unassigned").any():
        bls.append("Unassigned")
    palette = {**BL_COLORS, "Unassigned": "#8a8a8a"}

    xerr_lo = m["predicted_distinct_Fgs"] - m["predicted_distinct_Fgs_lo"]
    xerr_hi = m["predicted_distinct_Fgs_hi"] - m["predicted_distinct_Fgs"]
    yerr_lo = m["predicted_P_compat_mean"] - m["predicted_P_compat_lo"]
    yerr_hi = m["predicted_P_compat_hi"] - m["predicted_P_compat_mean"]

    for bl in bls:
        mask = m["BL"] == bl
        colour = palette[bl]
        ax.errorbar(
            m.loc[mask, "predicted_distinct_Fgs"],
            m.loc[mask, "predicted_P_compat_mean"],
            xerr=[xerr_lo[mask], xerr_hi[mask]],
            yerr=[yerr_lo[mask], yerr_hi[mask]],
            fmt="none", ecolor=colour, alpha=0.35,
            elinewidth=1.0, capsize=2.5, zorder=1,
        )
        sizes = 30 + 6 * np.sqrt(np.clip(m.loc[mask, "M_mothers_in_db"], 1, None))
        ax.scatter(
            m.loc[mask, "predicted_distinct_Fgs"],
            m.loc[mask, "predicted_P_compat_mean"],
            s=sizes, c=colour, edgecolor="white", linewidth=0.6,
            label=f"{bl}  (n = {int(mask.sum())})",
            alpha=0.9, zorder=2,
        )

    # Label the extreme dots (top/bottom-3 by predicted diversity)
    m_sorted = m.sort_values("predicted_distinct_Fgs")
    for _, row in pd.concat([m_sorted.head(3), m_sorted.tail(3)]).iterrows():
        ax.text(row["predicted_distinct_Fgs"] + 0.3,
                row["predicted_P_compat_mean"],
                f" {row['locationCode']}",
                fontsize=8, color="#222", va="center")

    ax.set_xlabel(
        "Predicted number of distinct SRK alleles detected at the location",
        fontsize=11,
    )
    ax.set_ylabel(
        "Predicted pollen compatibility under random mating "
        "(mean per location)",
        fontsize=11,
    )
    ax.set_xlim(0, len(prior) + 2)
    y_upper = max(0.30, min(1.0, 2.0 * SPECIES_MEAN))
    ax.set_ylim(0.0, y_upper)
    ax.set_title(
        f"SRK allele diversity vs pollen compatibility across LEPA locations"
        f"{_year_suffix(year)}",
        fontsize=13,
    )
    ax.legend(loc="lower right", fontsize=9, frameon=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)


def plot_mate_limitation(per_location: pd.DataFrame, coefs: pd.DataFrame,
                         out_png: Path, out_pdf: Path,
                         year: int | None = None,
                         demo: bool = False,
                         bands: dict[str, float] | None = None):
    """Single-panel LOCATION-level scatter of mean seeds per mother vs
    location-mean random-mating compatibility. Dot size = √n_mothers;
    error bars on both axes = per-mother SEM. Background traffic-light
    tiers (failed / struggling / sustainable) match the wording used
    elsewhere in the SRK random-mating framework. Weighted OLS slope +
    coefficient p-value are printed in the legend. The coefficient
    detail (β₁, β₂ CI) lives in the accompanying tables:
        step30_B_mate_limitation_coefficients.tsv
        step30_B_mate_limitation_per_location.tsv
    """
    fig, ax_scatter = plt.subplots(1, 1, figsize=(11.0, 6.5))

    # ---- location-level scatter ----
    x = per_location["mean_P_compat"].values
    y = per_location["mean_seeds"].values
    w = per_location["n_mothers"].values
    # y-error = SEM per location (std / sqrt(n))
    y_sem = per_location["std_seeds"] / np.sqrt(per_location["n_mothers"].clip(lower=1))
    # x-error = SEM of within-location P_compat variability across mothers
    x_sem = per_location["std_P_compat"] / np.sqrt(per_location["n_mothers"].clip(lower=1))

    # Traffic-light bands: failed (red) / struggling (orange) / sustainable
    # (green). Thresholds recalibrated against the sporophytic tetraploid
    # species mean (§ A.6, Part 2). Defaults preserve backward-compatible
    # gametophytic thresholds if `bands` is not supplied.
    T_FAILED_HI     = bands["failed_max"]     if bands else 0.20
    T_STRUGGLING_HI = bands["struggling_max"] if bands else 0.40
    SPECIES_MEAN    = bands["species_mean"]   if bands else 0.63
    x_upper = max(0.30, min(1.0, 2.0 * SPECIES_MEAN))
    band_alpha  = 0.14
    ax_scatter.axvspan(0.0, T_FAILED_HI,
                        color="#b2182b", alpha=band_alpha, zorder=0)
    ax_scatter.axvspan(T_FAILED_HI, T_STRUGGLING_HI,
                        color="#e08214", alpha=band_alpha, zorder=0)
    ax_scatter.axvspan(T_STRUGGLING_HI, x_upper,
                        color="#1b7837", alpha=band_alpha, zorder=0)
    ax_scatter.axvline(SPECIES_MEAN, color="#1b7837", ls=":", lw=1.0, alpha=0.7)

    ax_scatter.errorbar(
        x, y, xerr=x_sem, yerr=y_sem,
        fmt="none", ecolor="#a6cee3", elinewidth=1.0, capsize=2.5, zorder=1,
    )
    sizes = 40 + 12 * np.sqrt(w)   # dot size ~ √n_mothers (slightly larger)
    # Dot colour by traffic-light tier at the mean compatibility
    def _tier_colour(pc: float) -> str:
        if pc < T_FAILED_HI:     return "#b2182b"
        if pc < T_STRUGGLING_HI: return "#e08214"
        return "#1b7837"
    dot_colours = [_tier_colour(pc) for pc in x]
    sc = ax_scatter.scatter(x, y, s=sizes,
                             c=dot_colours, edgecolor="white", linewidth=0.6,
                             alpha=0.9, zorder=2)

    # locationCode labels for the extreme dots (top / bottom 5)
    label_col = "locationCode" if "locationCode" in per_location.columns else "locationID"
    ranked = per_location.assign(_resid=(y - np.mean(y))).sort_values("_resid")
    for _, row in pd.concat([ranked.head(5), ranked.tail(5)]).iterrows():
        ax_scatter.text(row["mean_P_compat"] + 0.006, row["mean_seeds"],
                        f" {row[label_col]}", fontsize=8, color="#222",
                        va="center")

    # Weighted OLS fit line (visual reference; official coefficients live in
    # step30_B_mate_limitation_coefficients.tsv).
    slope_line_label = None
    if len(per_location) > 2:
        w_norm = w / w.sum()
        x_bar = np.sum(w_norm * x); y_bar = np.sum(w_norm * y)
        num = np.sum(w * (x - x_bar) * (y - y_bar))
        den = np.sum(w * (x - x_bar) ** 2)
        slope = num / den if den > 0 else 0.0
        intercept = y_bar - slope * x_bar
        xs = np.linspace(x.min(), x.max(), 100)
        # p-value on β₁ from the accompanying coefs table (P_compat term)
        p_val = float(coefs.loc[coefs["term"] == "P_compat", "p_value"].iloc[0])
        p_str = "< 1e-4" if p_val < 1e-4 else f"= {p_val:.3g}"
        slope_line_label = (f"fitted slope = {slope:.1f}   "
                            f"(p {p_str})")
        ax_scatter.plot(xs, intercept + slope * xs,
                        color="#333333", ls="--", lw=1.6, alpha=0.9,
                        label=slope_line_label, zorder=3)

    # Traffic-light legend — sporophytic bands rescaled to species mean.
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    tl_handles = [
        Patch(facecolor="#b2182b", alpha=0.35,
              label=f"failed  (< {T_FAILED_HI:.3f})"),
        Patch(facecolor="#e08214", alpha=0.35,
              label=f"struggling  ({T_FAILED_HI:.3f} – {T_STRUGGLING_HI:.3f})"),
        Patch(facecolor="#1b7837", alpha=0.35,
              label=f"sustainable  (≥ {T_STRUGGLING_HI:.3f})"),
        Line2D([0], [0], color="#1b7837", ls=":", lw=1.2,
               label=f"species mean  ({SPECIES_MEAN:.3f})"),
    ]
    ax_scatter.set_xlim(0.0, x_upper)
    lg_handles, lg_labels = ax_scatter.get_legend_handles_labels()
    ax_scatter.legend(
        lg_handles + tl_handles,
        lg_labels + [h.get_label() for h in tl_handles],
        loc="upper left", fontsize=10, frameon=True,
    )

    ax_scatter.set_xlabel(
        "Predicted pollen compatibility under random mating  "
        "(mean across sampled mothers per location)",
        fontsize=11,
    )
    ax_scatter.set_ylabel("Observed seeds per mother  (mean per location)",
                          fontsize=11)
    ax_scatter.set_xlim(0.0, 1.0)
    ax_scatter.set_title(
        f"Mate limitation across LEPA locations{_year_suffix(year)}"
        f"{_demo_suffix(demo)}",
        fontsize=13,
    )
    _add_demo_watermark(fig, demo)
    ax_scatter.spines["top"].set_visible(False)
    ax_scatter.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200); fig.savefig(out_pdf); plt.close(fig)


def plot_si_escape(perm: pd.DataFrame,
                    location_codes: dict[int, str] | None,
                    out_png: Path, out_pdf: Path,
                    q_threshold: float = 0.05,
                    year: int | None = None,
                    demo: bool = False):
    """Per-location bars of pi_self_match with FDR-significant locations
    highlighted red. Vertical line at 0 (strict-SI null).
    Y-axis uses `location_codes` mapping (locationID → human-readable code)
    when supplied, so bars are labelled by locationCode rather than opaque
    integer IDs."""
    df = perm.sort_values("pi_self_match", ascending=True).reset_index(drop=True)
    colours = np.where(
        (df["q_bh_fdr"] <= q_threshold) & (df["pi_self_match"] > 0),
        "#b2182b", "#8a8a8a",
    )
    fig, ax = plt.subplots(figsize=(8.5, max(4.0, 0.20 * len(df) + 1.5)))
    y = np.arange(len(df))
    ax.barh(y, df["pi_self_match"], color=colours, edgecolor="white", height=0.75)
    ax.axvline(0, color="#333", ls="--", lw=1.0, alpha=0.6)
    if location_codes:
        labels = [location_codes.get(int(lid), str(lid))
                  for lid in df["locationID"]]
    else:
        labels = df["locationID"].astype(str).tolist()
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Observed rate of pollen alleles matching the mother's "
                  "own alleles  (= self-incompatibility escape rate)")
    n_sig = int(((df['q_bh_fdr'] <= q_threshold) & (df['pi_self_match'] > 0)).sum())
    ax.set_title(
        f"Self-incompatibility escape signal across LEPA locations"
        f"{_year_suffix(year)}   —   {n_sig} / {len(df)} locations "
        f"reject strict self-incompatibility "
        f"({int(q_threshold*100)} % false-discovery rate)"
        f"{_demo_suffix(demo)}",
        fontsize=12,
    )
    _add_demo_watermark(fig, demo)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor="#b2182b",
              label=f"reject strict self-incompatibility "
                    f"({int(q_threshold*100)} % false-discovery rate)"),
        Patch(facecolor="#8a8a8a",
              label="consistent with strict self-incompatibility"),
    ], loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200); fig.savefig(out_pdf); plt.close(fig)


def plot_comparison_diversity(merged: pd.DataFrame,
                              out_png: Path, out_pdf: Path,
                              use_matched: bool = False,
                              year: int | None = None,
                              demo: bool = False):
    """Predicted vs observed # distinct Fgs per location.

    use_matched=False → x = predicted at M mothers (may sit below observed
                         if far more seeds were genotyped than M implies).
    use_matched=True  → x = predicted at the actual n_seeds × 2 allele draws
                         per location (like-for-like scale).
    """
    if use_matched:
        x     = merged["predicted_distinct_Fgs_matched"]
        x_lo  = merged["predicted_distinct_Fgs_matched_lo"]
        x_hi  = merged["predicted_distinct_Fgs_matched_hi"]
        xlab = "Predicted number of distinct SRK alleles"
    else:
        x     = merged["predicted_distinct_Fgs"]
        x_lo  = merged["predicted_distinct_Fgs_lo"]
        x_hi  = merged["predicted_distinct_Fgs_hi"]
        xlab = "Predicted number of distinct SRK alleles"

    fig, ax = plt.subplots(figsize=(7.5, 6.0))
    ax.errorbar(
        x, merged["observed_distinct_Fgs"],
        xerr=[x - x_lo, x_hi - x],
        fmt="o", color="#1f78b4", ecolor="#a6cee3",
        elinewidth=1.2, capsize=2.5, markersize=5, alpha=0.85,
    )
    lim = max(float(x_hi.max()),
              float(merged["observed_distinct_Fgs"].max())) + 1
    ax.plot([0, lim], [0, lim], color="#333", ls="--", lw=1.0, alpha=0.6)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel(xlab)
    ax.set_ylabel("Observed number of distinct SRK alleles in seed genotypes")
    ax.set_title(
        f"Observed vs predicted SRK allele diversity across LEPA locations"
        f"{_year_suffix(year)}{_demo_suffix(demo)}",
        fontsize=13,
    )
    _add_demo_watermark(fig, demo)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 6 — main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed-genotypes", type=Path,
                    help="TSV with columns locationID, eventID, germplasmID, "
                         "seed_id, maternal_Fg, paternal_Fg. If omitted, only "
                         "the prediction family is produced.")
    ap.add_argument("--demo", action="store_true",
                    help="Simulate seed genotypes from the P1 prior to "
                         "demonstrate the comparison pipeline.")
    ap.add_argument("--demo-seeds-per-mother", type=int, default=29,
                    help="Seeds per mother in demo mode (default: Step 28 "
                         "miss-probability rule n = 29).")
    ap.add_argument("--demo-si-escape-rate", type=float, default=0.08,
                    help="Fraction of seeds where SI is bypassed in the "
                         "demo simulator. 0 = strict SI. Default 0.08 so the "
                         "SI-escape permutation test has signal to detect.")
    ap.add_argument("--mother-genotypes", type=Path,
                    help="Optional TSV of per-mother metadata "
                         "(locationID, germplasmID, mother_Fg_a, mother_Fg_b, "
                         "K_spatial, seeds_est). Required for real-data mode "
                         "of the mate-limitation regression and SI-escape "
                         "permutation test.")
    ap.add_argument("--n-mothers-pcompat", type=int, default=5000,
                    help="Number of simulated mothers for the P_compat "
                         "distribution (default: 5000).")
    ap.add_argument("--match-seed-count", action="store_true",
                    help="Add a scale-matched prediction column to the "
                         "comparison output whose sample size equals the "
                         "actual n_seed_alleles per location (2 x #seeds). "
                         "The diversity comparison figure then uses that "
                         "column instead of the M-mothers prediction.")
    ap.add_argument("--year", type=int, default=None,
                    help="Survey year to show in figure titles (display "
                         "only). The upstream filtering by year is done "
                         "in Steps 28 and 29; this flag simply threads the "
                         "same value into the Step 30 figure titles.")
    args = ap.parse_args()

    tables_dir  = Path(DEFAULT_TABLES); tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = Path(DEFAULT_FIGURES); figures_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)

    # ---- Load prior ----
    prior = build_p1_prior(DEFAULT_PRIOR_TSV)
    print(f"[step30] P1 prior — {len(prior)} Fgs, "
          f"ESS = {PRIOR_ESS}, most common Fg = "
          f"{prior.iloc[0]['Fg']} ({prior.iloc[0]['f_mean']:.2%}).")
    prior.to_csv(tables_dir / "step30_A_prediction_prior_frequencies.tsv",
                 sep="\t", index=False)

    # ---- Load locations from Step 29 ----
    if not DEFAULT_LOCATIONS_TSV.exists():
        raise SystemExit(
            f"[step30] Missing {DEFAULT_LOCATIONS_TSV} — run step29 first.")
    locations = pd.read_csv(DEFAULT_LOCATIONS_TSV, sep="\t",
                             encoding="utf-8-sig")

    # Weave in within-location connectivity from Step 28c. Connectivity
    # is a **predictor** of both realised SRK diversity and random-mating
    # compatibility: an event that shares no mating neighbourhood with
    # the rest of its "location" behaves as a smaller population and
    # therefore sees more drift and lower compatibility. We use
    # `largest_component_share_50m` (the share of adults in the largest
    # within-location connected component at 10 m) to define the
    # effective mating N; if connectivity data are missing the fallback
    # is `total_n_fertile` (i.e. treat the location as one unit).
    if DEFAULT_CONNECTIVITY_TSV.exists():
        conn = pd.read_csv(DEFAULT_CONNECTIVITY_TSV, sep="\t",
                            encoding="utf-8-sig")
        keep = ["locationID", "n_events", "total_n_fertile"] + [
            c for c in conn.columns
            if c.startswith(("connected_share_", "largest_component_share_",
                             "n_components_"))
        ]
        conn = conn[keep].rename(columns={
            "n_events":        "conn_n_events",
            "total_n_fertile": "conn_total_n_fertile",
        })
        locations = locations.merge(conn, on="locationID", how="left")
        # Effective mating N = adults in the largest connected component
        # at the 50 m primary pollen-flight radius. Use step29b's OWN
        # census (conn_total_n_fertile) as the base — the connectivity
        # share was computed against that number, not step29's SQL
        # census, so multiplying by step29's number produces a subtly
        # wrong N_fert_eff. Also overwrite total_n_fertile so every
        # downstream consumer sees the spatially-analysable census
        # (some events may be dropped by the coord filter — those are
        # the ~2-4 % discrepancy between step29 and step29b totals).
        locations["total_n_fertile"] = locations["conn_total_n_fertile"].astype(int)
        locations["N_fertile_effective_50m"] = (
            locations["conn_total_n_fertile"].astype(float)
            * locations["largest_component_share_50m"].fillna(1.0)
        ).round().astype(int)
        print(f"[step30] Woven in within-location connectivity from "
              f"{DEFAULT_CONNECTIVITY_TSV.name} "
              f"(mean largest-component share at 50 m = "
              f"{locations['largest_component_share_50m'].mean():.2f}).")
    else:
        print(f"[step30] {DEFAULT_CONNECTIVITY_TSV} not found — "
              "connectivity predictor unavailable. Run step28c first.")
        locations["largest_component_share_50m"] = 1.0
        locations["N_fertile_effective_50m"] = locations["total_n_fertile"]

    # ---- Prediction outputs ----
    # Per-50 m-component drift + sampling simulation; location-level
    # pool size and coverage derived by set-union across components.
    pred_div, pred_div_comp = predicted_diversity_per_location(
        locations, prior, rng)
    pred_div_path = tables_dir / "step30_A_prediction_location_diversity.tsv"
    pred_div.to_csv(pred_div_path, sep="\t", index=False)
    print(f"[step30] Wrote {pred_div_path}")
    pred_div_comp_path = tables_dir / "step30_A_prediction_component_diversity.tsv"
    pred_div_comp.to_csv(pred_div_comp_path, sep="\t", index=False)
    print(f"[step30] Wrote {pred_div_comp_path} "
          f"({len(pred_div_comp)} rows across "
          f"{pred_div_comp['locationID'].nunique()} locations).")

    # Per-mother reference distribution (kept for downstream reference).
    pcompat_ref = predicted_pcompat_distribution(prior, args.n_mothers_pcompat, rng)
    pcompat_ref_path = tables_dir / "step30_A_prediction_per_mother_fecundation.tsv"
    pcompat_ref.to_csv(pcompat_ref_path, sep="\t", index=False)
    print(f"[step30] Wrote {pcompat_ref_path} "
          f"(species-wide mean P_compat = {pcompat_ref['P_compat'].mean():.3f})")

    # Per-location Phase-A prediction, built from per-component
    # P_compat simulations aggregated by component-size weighting.
    # Each 50 m connected component inside a location gets its own
    # pool, mothers, P_compat; the location-level number is the
    # size-weighted mean. The finer-grained per-component table is
    # saved alongside so readers can spot fragmented locations whose
    # headline mean hides a struggling sub-component.
    pcompat_loc, pcompat_comp = predicted_pcompat_per_location(
        locations, prior, rng)
    pcompat_loc_path = tables_dir / "step30_A_prediction_location_pcompat.tsv"
    pcompat_loc.to_csv(pcompat_loc_path, sep="\t", index=False)
    print(f"[step30] Wrote {pcompat_loc_path}")
    pcompat_comp_path = tables_dir / "step30_A_prediction_component_pcompat.tsv"
    pcompat_comp.to_csv(pcompat_comp_path, sep="\t", index=False)
    print(f"[step30] Wrote {pcompat_comp_path} "
          f"({len(pcompat_comp)} rows across "
          f"{pcompat_comp['locationID'].nunique()} locations).")

    # --- Recalibrate traffic-light bands against the sporophytic
    # + empirical-zygosity species mean (see § A.7)
    from srk_si_model import (
        load_class_map, build_class_i_mask,
        load_zygosity_dist,
        species_mean_p_compat_empirical, traffic_light_bands,
    )
    fg_labels_ordered = prior["Fg"].astype(str).tolist() \
        if "Fg" in prior.columns \
        else [f"FG{i+1:03d}" for i in range(len(prior))]
    class_map = load_class_map()
    class_i_mask = build_class_i_mask(fg_labels_ordered, class_map)
    zygosity_probs_main = load_zygosity_dist()
    species_mean_pc = species_mean_p_compat_empirical(
        prior["f_mean"].values, class_i_mask, zygosity_probs_main,
        n_mothers=10_000, n_fathers=1_500, rng=rng)
    bands = traffic_light_bands(species_mean_pc)
    bands_df = pd.DataFrame([{
        "si_model":                    "sporophytic_class_I_dominant_empirical_zygosity",
        "species_mean":                bands["species_mean"],
        "failed_max":                  bands["failed_max"],
        "struggling_max":              bands["struggling_max"],
        "class_I_count":               int(class_i_mask.sum()),
        "class_I_P1_freq":             float(prior["f_mean"].values[class_i_mask].sum()),
        "zygosity_p_1_distinct":       float(zygosity_probs_main[0]),
        "zygosity_p_2_distinct":       float(zygosity_probs_main[1]),
        "zygosity_p_3_distinct":       float(zygosity_probs_main[2]),
        "zygosity_p_4_distinct":       float(zygosity_probs_main[3]),
    }])
    bands_path = tables_dir / "step30_A_traffic_light_bands.tsv"
    bands_df.to_csv(bands_path, sep="\t", index=False)
    print(f"[step30] Sporophytic + empirical-zygosity species-mean "
          f"P_compat = {bands['species_mean']:.4f}; bands failed < "
          f"{bands['failed_max']:.4f}, struggling < "
          f"{bands['struggling_max']:.4f}, sustainable ≥ "
          f"{bands['struggling_max']:.4f}.")
    print(f"[step30] Wrote {bands_path}")

    # Legacy N_fertile_effective figure dropped (2026-10-03) — the
    # per-location deme structure is now shown directly in
    # step29d_mating_pool_structure.py, and the per-component drift
    # simulations (Figures 3 and 4) are the authoritative display.
    plot_diversity_unbiased_vs_sampling(pred_div, prior,
        out_png=figures_dir / "step30_A_diversity_unbiased_vs_sampling.png",
        out_pdf=figures_dir / "step30_A_diversity_unbiased_vs_sampling.pdf",
        year=args.year)
    plot_prediction_fecundation(pcompat_loc,
        out_png=figures_dir / "step30_A_prediction_fecundation.png",
        out_pdf=figures_dir / "step30_A_prediction_fecundation.pdf",
        year=args.year, bands=bands)
    plot_diversity_vs_pcompat(pred_div, pcompat_loc, prior,
        out_png=figures_dir / "step30_A_diversity_vs_pcompat.png",
        out_pdf=figures_dir / "step30_A_diversity_vs_pcompat.pdf",
        year=args.year, bands=bands)
    print(f"[step30] Phase A prediction figures in {figures_dir}/")

    # ---- Comparison outputs (optional) ----
    seeds = None
    mothers = None
    if args.seed_genotypes:
        seeds = pd.read_csv(args.seed_genotypes, sep="\t",
                            encoding="utf-8-sig")
        if args.mother_genotypes:
            mothers = pd.read_csv(args.mother_genotypes, sep="\t",
                                   encoding="utf-8-sig")
    elif args.demo:
        seeds, mothers = simulate_seed_genotypes(
            locations, prior,
            seeds_per_mother=args.demo_seeds_per_mother, rng=rng,
            si_escape_rate=args.demo_si_escape_rate,
        )
        demo_seeds_path   = tables_dir / "step30_B_DEMO_seed_genotypes.tsv"
        demo_mothers_path = tables_dir / "step30_B_DEMO_mother_genotypes.tsv"
        seeds.to_csv(demo_seeds_path, sep="\t", index=False)
        mothers.to_csv(demo_mothers_path, sep="\t", index=False)
        print(f"[step30] Simulated {len(seeds)} seeds from {len(mothers)} "
              f"mothers across {seeds['locationID'].nunique()} locations "
              f"(SI-escape rate = {args.demo_si_escape_rate:.2%}) "
              f"→ {demo_seeds_path}, {demo_mothers_path}")

    if seeds is not None:
        # Distinguish demo (synthetic) Phase B outputs from real ones via
        # filename infix, so a `--demo` run never overwrites a real Phase B
        # deliverable and never sits next to real files without a warning.
        # Titles get a "[DEMO — synthetic data]" suffix and figures get a
        # diagonal DEMO watermark; both make the source unambiguous.
        is_demo = bool(args.demo)
        tag = "_DEMO" if is_demo else ""

        post = location_posterior(seeds, prior, rng)
        merged = pred_div.merge(post, on="locationID", how="inner")
        if args.match_seed_count:
            matched = predicted_diversity_matched_to_seeds(seeds, prior, rng)
            merged = merged.merge(matched, on="locationID", how="left")
        comparison_path = tables_dir / f"step30_B{tag}_comparison_location_diversity.tsv"
        merged.to_csv(comparison_path, sep="\t", index=False)
        print(f"[step30] Wrote {comparison_path}"
              + ("  (scale-matched prediction included)"
                 if args.match_seed_count else ""))

        mp = maternal_vs_paternal(seeds)
        mp_path = tables_dir / f"step30_B{tag}_comparison_maternal_vs_paternal.tsv"
        mp.to_csv(mp_path, sep="\t", index=False)
        print(f"[step30] Wrote {mp_path}")

        plot_comparison_diversity(merged,
            out_png=figures_dir / f"step30_B{tag}_comparison_diversity.png",
            out_pdf=figures_dir / f"step30_B{tag}_comparison_diversity.pdf",
            use_matched=args.match_seed_count,
            year=args.year, demo=is_demo)
        print(f"[step30] Comparison figure in {figures_dir}/")

        # ---- Test 1: Mate-limitation regression (§ 2.2.1) ----
        if mothers is not None:
            conn_df = None
            if "largest_component_share_50m" in locations.columns:
                conn_df = locations[[
                    "locationID", "largest_component_share_50m"
                ]].copy()
            reg_df, per_loc_reg, coefs = mate_limitation_regression(
                seeds, mothers, prior, connectivity=conn_df)
            reg_path      = tables_dir / f"step30_B{tag}_mate_limitation_per_mother.tsv"
            per_loc_path  = tables_dir / f"step30_B{tag}_mate_limitation_per_location.tsv"
            coefs_path    = tables_dir / f"step30_B{tag}_mate_limitation_coefficients.tsv"
            reg_df.to_csv(reg_path, sep="\t", index=False)
            per_loc_reg.to_csv(per_loc_path, sep="\t", index=False)
            coefs.to_csv(coefs_path, sep="\t", index=False)
            print(f"[step30] Test 1  Mate-limitation regression (LOCATION "
                  f"level, {len(per_loc_reg)} locations) → {per_loc_path}")
            print(f"[step30]                                per-mother detail "
                  f"→ {reg_path}")
            for _, coef_row in coefs.iterrows():
                print(f"[step30]   β ({coef_row['term']}) = "
                      f"{coef_row['estimate']:.3f} "
                      f"[95 % CI {coef_row['ci_lo']:.3f}, "
                      f"{coef_row['ci_hi']:.3f}]  "
                      f"p = {coef_row['p_value']:.3g}"
                      + ("  SIGNIFICANT" if coef_row["significant_95"] else ""))
            plot_mate_limitation(per_loc_reg, coefs,
                out_png=figures_dir / f"step30_B{tag}_mate_limitation.png",
                out_pdf=figures_dir / f"step30_B{tag}_mate_limitation.pdf",
                year=args.year, demo=is_demo, bands=bands)

            # ---- Test 2: SI-escape permutation (§ 2.2.2) ----
            perm = si_escape_permutation(seeds, mothers, rng=rng)
            perm_path = tables_dir / f"step30_B{tag}_si_escape_permutation.tsv"
            perm.to_csv(perm_path, sep="\t", index=False)
            n_sig = int((perm["q_bh_fdr"] <= 0.05).sum())
            print(f"[step30] Test 2  SI-escape permutation written to "
                  f"{perm_path} — {n_sig}/{len(perm)} locations reject "
                  f"strict-SI at FDR ≤ 0.05.")
            # locationCode mapping for the SI-escape figure
            loc_code_map = None
            if "locationCode" in seeds.columns:
                loc_code_map = (seeds.drop_duplicates("locationID")
                                .set_index("locationID")["locationCode"].to_dict())
            plot_si_escape(perm, loc_code_map,
                out_png=figures_dir / f"step30_B{tag}_si_escape.png",
                out_pdf=figures_dir / f"step30_B{tag}_si_escape.pdf",
                year=args.year, demo=is_demo)
        else:
            print("[step30] Skipping Tests 1 + 2 — no --mother-genotypes "
                  "provided (or --demo not enabled).")


if __name__ == "__main__":
    main()
