"""Step 30g — Phase A predictions per (populationID, year).

Phase 5 § A.5.0 companion. Takes the per-year deme partition built
in Step 29a and runs the standard Phase A simulation (per-deme
draws from P1 → union to the population) at each (populationID,
year) combination. The point: within a population, 2025 and 2026
above-ground samples are two independent draws from the same
seed bank, so the across-year gap in predicted SRK diversity /
pollen compatibility is a direct estimate of per-population drift
beyond the species-wide prior (§ C.0.c H2).

Outputs
-------
Tables/Phase5/step30g_prediction_population_year.tsv
    one row per (populationID, year) with:
      - n_demes, N_fertile_total
      - pred_srk_diversity_mean / lo95 / hi95
      - pred_pcompat_mean / lo95 / hi95

Tables/Phase5/step30g_across_year_comparison.tsv
    one row per populationID present in BOTH years with:
      - pred_srk_diversity 2025 vs 2026 (delta + CI overlap flag)
      - pred_pcompat 2025 vs 2026 (delta + CI overlap flag)
      - a candidate flag for the preliminary large/small pair
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from step28_seed_sampling_per_mother import PLOIDY
from step30_srk_diversity_prediction_vs_observed import build_p1_prior
from srk_si_model import (
    load_class_map, build_class_i_mask, load_zygosity_dist,
    sample_genotypes_empirical, p_compat_sporophytic_empirical,
)

TABLES = Path("Tables/Phase5")

N_REPLICATES = 150               # per-population-year replicates
N_FATHERS    = 300               # candidate fathers in each p_compat call
N_MOTHERS    = 60                # mothers simulated per population-year
RNG_SEED     = 2030              # (lighter params than step30; this is a
                                 #  comparison across pop-years, not an
                                 #  absolute-CI product. Up these if the
                                 #  pcompat CIs are too wide to call
                                 #  stable vs unstable.)

# ---- Candidate-selection tunables ----
# Site visitation effort was equal across 2025 and 2026. The extra 2026
# occurrence records (mother plants collected for seed banking) are not
# used at this Phase 5 stage, so raw N_fertile sums are the appropriate
# across-year comparison axis. Mean-per-event is kept as an informational
# column in the comparison TSV but is not the primary classifier.
MIN_N_FERT_SMALL_STABLE = 20     # SMALL stable candidate floor on
                                 # N_fert_total across both years.
CRASH_RATIO_THRESHOLD   = 0.25   # raw-N_fert ratio N_fert_2026 /
                                 # N_fert_2025 <= this AND negative
                                 # delta_diversity.


def simulate_population_year(deme_sizes: list[int],
                               p1_vec: np.ndarray,
                               class_i_mask: np.ndarray,
                               zygosity_probs: np.ndarray,
                               rng: np.random.Generator,
                               ) -> tuple[np.ndarray, np.ndarray]:
    """Return (diversity_counts, pcompat_means) each of length
    N_REPLICATES.

    Diversity per replicate: for each deme of size N draw PLOIDY*N
    alleles from P1, record present Fgs; union across demes.
    P_compat per replicate: build the local Fg frequency vector as
    the union-weighted mean of the per-deme local frequencies, then
    run p_compat_sporophytic_empirical on N_MOTHERS mothers from that
    pool and take the mean."""
    K = len(p1_vec)
    diversity = np.zeros(N_REPLICATES, dtype=int)
    pcompat = np.zeros(N_REPLICATES, dtype=float)
    for i in range(N_REPLICATES):
        present = np.zeros(K, dtype=bool)
        local_counts = np.zeros(K, dtype=int)
        for size in deme_sizes:
            n_alleles = PLOIDY * int(size)
            if n_alleles <= 0:
                continue
            draws = rng.choice(K, size=n_alleles, replace=True, p=p1_vec)
            present[draws] = True
            np.add.at(local_counts, draws, 1)
        diversity[i] = int(present.sum())
        if local_counts.sum() == 0:
            pcompat[i] = np.nan
            continue
        local_f = local_counts / local_counts.sum()
        mothers = sample_genotypes_empirical(N_MOTHERS, local_f, zygosity_probs, rng)
        per_mother = p_compat_sporophytic_empirical(
            mothers, local_f, class_i_mask, zygosity_probs,
            n_fathers=N_FATHERS, rng=rng,
        )
        pcompat[i] = float(np.nanmean(per_mother))
    return diversity, pcompat


def ci_overlap(lo1: float, hi1: float, lo2: float, hi2: float) -> bool:
    """True if the two 95 % CIs overlap."""
    return not (hi1 < lo2 or hi2 < lo1)


def main() -> None:
    rng = np.random.default_rng(RNG_SEED)

    # ---- Load prerequisites ----
    prior = build_p1_prior(TABLES / "step26i_L1_carrier_inventory.tsv")
    p1_vec = prior["f_mean"].values
    fg_labels = prior["Fg"].tolist()
    class_i_mask = build_class_i_mask(fg_labels, load_class_map())
    zygosity_probs = load_zygosity_dist()
    print(f"[step30g] P1: {len(p1_vec)} Fgs; "
          f"Class I = {int(class_i_mask.sum())}, Class II = {int((~class_i_mask).sum())}; "
          f"zygosity probs {zygosity_probs}")

    demes = pd.read_csv(TABLES / "step29a_demes_per_population_year.tsv",
                         sep="\t", encoding="utf-8-sig")

    # ---- Per (populationID, year) prediction ----
    rows: list[dict] = []
    for (pop_id, yr), sub in demes.groupby(["populationID", "year"]):
        deme_sizes = sub["component_N_fertile"].astype(int).tolist()
        div, pc = simulate_population_year(
            deme_sizes, p1_vec, class_i_mask, zygosity_probs, rng,
        )
        rows.append({
            "populationID":            int(pop_id),
            "year":                    int(yr),
            "n_demes":                 len(deme_sizes),
            "N_fertile_total":         int(sum(deme_sizes)),
            "largest_deme_N":          int(max(deme_sizes)),
            "smallest_deme_N":         int(min(deme_sizes)),
            "pred_srk_diversity_mean": float(div.mean()),
            "pred_srk_diversity_lo95": float(np.quantile(div, 0.025)),
            "pred_srk_diversity_hi95": float(np.quantile(div, 0.975)),
            "pred_pcompat_mean":       float(np.nanmean(pc)),
            "pred_pcompat_lo95":       float(np.nanquantile(pc, 0.025)),
            "pred_pcompat_hi95":       float(np.nanquantile(pc, 0.975)),
        })
    pred = pd.DataFrame(rows)
    pred.to_csv(TABLES / "step30g_prediction_population_year.tsv",
                 sep="\t", index=False)
    print(f"[step30g] Wrote step30g_prediction_population_year.tsv "
          f"({len(pred)} rows across {pred['populationID'].nunique()} populations)")

    # ---- Across-year comparison ----
    wide = pred.pivot(index="populationID", columns="year",
                       values=["pred_srk_diversity_mean",
                               "pred_srk_diversity_lo95",
                               "pred_srk_diversity_hi95",
                               "pred_pcompat_mean",
                               "pred_pcompat_lo95",
                               "pred_pcompat_hi95",
                               "n_demes",
                               "N_fertile_total"])
    both = wide.dropna(subset=[("pred_pcompat_mean", 2025),
                                ("pred_pcompat_mean", 2026)]).copy()

    comp_rows: list[dict] = []
    for pop_id, r in both.iterrows():
        d25 = r[("pred_srk_diversity_mean", 2025)]
        d26 = r[("pred_srk_diversity_mean", 2026)]
        p25 = r[("pred_pcompat_mean", 2025)]
        p26 = r[("pred_pcompat_mean", 2026)]
        div_overlap = ci_overlap(
            r[("pred_srk_diversity_lo95", 2025)],
            r[("pred_srk_diversity_hi95", 2025)],
            r[("pred_srk_diversity_lo95", 2026)],
            r[("pred_srk_diversity_hi95", 2026)],
        )
        pc_overlap = ci_overlap(
            r[("pred_pcompat_lo95", 2025)],
            r[("pred_pcompat_hi95", 2025)],
            r[("pred_pcompat_lo95", 2026)],
            r[("pred_pcompat_hi95", 2026)],
        )
        comp_rows.append({
            "populationID":            int(pop_id),
            "n_demes_2025":            int(r[("n_demes", 2025)]),
            "n_demes_2026":            int(r[("n_demes", 2026)]),
            "N_fert_2025":             int(r[("N_fertile_total", 2025)]),
            "N_fert_2026":             int(r[("N_fertile_total", 2026)]),
            "N_fert_total":            int(r[("N_fertile_total", 2025)]
                                            + r[("N_fertile_total", 2026)]),
            "pred_diversity_2025":     float(d25),
            "pred_diversity_2026":     float(d26),
            "delta_diversity":         float(d26 - d25),
            "diversity_CI_overlap":    bool(div_overlap),
            "pred_pcompat_2025":       float(p25),
            "pred_pcompat_2026":       float(p26),
            "delta_pcompat":           float(p26 - p25),
            "pcompat_CI_overlap":      bool(pc_overlap),
            "stable_across_years":     bool(div_overlap and pc_overlap),
        })
    comp = (pd.DataFrame(comp_rows)
              .sort_values("N_fert_total", ascending=False)
              .reset_index(drop=True))
    comp.to_csv(TABLES / "step30g_across_year_comparison.tsv",
                 sep="\t", index=False)
    print(f"[step30g] Wrote step30g_across_year_comparison.tsv "
          f"({len(comp)} populations present in both years)")

    # ---- Candidate pair selection: stable LARGE + SMALL ----
    stable = comp[comp["stable_across_years"]].copy()
    stable_small_pool = stable[
        stable["N_fert_total"] >= MIN_N_FERT_SMALL_STABLE
    ].sort_values("N_fert_total", ascending=False)
    pair_rows: list[dict] = []
    if len(stable_small_pool) >= 2:
        large = stable_small_pool.iloc[0]
        small = stable_small_pool.iloc[-1]
        for role, row in [("LARGE", large), ("SMALL", small)]:
            pair_rows.append({"role": role, **row.to_dict()})
        pair = pd.DataFrame(pair_rows)
        pair.to_csv(TABLES / "step30g_stable_candidate_pair.tsv",
                     sep="\t", index=False)
        print(f"[step30g] Wrote step30g_stable_candidate_pair.tsv")
        print()
        print("[step30g] ========== Candidate pair (stable; N_fert_total >= "
              f"{MIN_N_FERT_SMALL_STABLE}) ==========")
        print(f"[step30g] LARGE : populationID={large['populationID']}  "
              f"N_fert 2025/2026={large['N_fert_2025']}/{large['N_fert_2026']}  "
              f"pred_pcompat 2025/2026={large['pred_pcompat_2025']:.3f}/{large['pred_pcompat_2026']:.3f}  "
              f"pred_diversity 2025/2026={large['pred_diversity_2025']:.1f}/{large['pred_diversity_2026']:.1f}")
        print(f"[step30g] SMALL : populationID={small['populationID']}  "
              f"N_fert 2025/2026={small['N_fert_2025']}/{small['N_fert_2026']}  "
              f"pred_pcompat 2025/2026={small['pred_pcompat_2025']:.3f}/{small['pred_pcompat_2026']:.3f}  "
              f"pred_diversity 2025/2026={small['pred_diversity_2025']:.1f}/{small['pred_diversity_2026']:.1f}")

    # ---- Secondary target: crash populations (clear 2026 < 2025 collapse) ----
    comp["crash_ratio"] = np.where(
        comp["N_fert_2025"] > 0,
        comp["N_fert_2026"] / comp["N_fert_2025"],
        np.nan,
    )
    crash = comp[
        (comp["crash_ratio"] <= CRASH_RATIO_THRESHOLD)
        & (comp["delta_diversity"] < 0)
        & (~comp["diversity_CI_overlap"])
    ].sort_values("crash_ratio").reset_index(drop=True)
    crash.to_csv(TABLES / "step30g_crash_candidates.tsv", sep="\t", index=False)
    print(f"[step30g] Wrote step30g_crash_candidates.tsv  "
          f"({len(crash)} crash populations)")
    if len(crash):
        print()
        print("[step30g] Crash candidates (2026 collapsed vs 2025; "
              "clearest within-pipeline H2 signal):")
        for _, r in crash.iterrows():
            print(f"[step30g]   populationID={int(r['populationID']):>3}  "
                  f"N_fert 2025/2026={int(r['N_fert_2025']):>4}/{int(r['N_fert_2026']):>4}  "
                  f"ratio={r['crash_ratio']:.3f}  "
                  f"Δdiv={r['delta_diversity']:+.1f}  "
                  f"Δpc={r['delta_pcompat']:+.3f}")

    # ---- Headlines ----
    unstable = comp[~comp["stable_across_years"]]
    print()
    print(f"[step30g] Across-year classification:")
    print(f"[step30g]   stable (both CIs overlap)   : {len(stable)}")
    print(f"[step30g]   unstable (one or both shift): {len(unstable)}")
    print(f"[step30g]   crash subset of unstable    : {len(crash)}")


if __name__ == "__main__":
    main()
