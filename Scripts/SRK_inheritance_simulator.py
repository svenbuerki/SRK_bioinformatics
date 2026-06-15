#!/usr/bin/env python3
"""
SRK_inheritance_simulator.py — Step 27

Forward-time Wright–Fisher tetraploid simulator with null-aware sporophytic
SI. Built on the empirical per-BL state from Step 26 to project the SI → SC
trajectory under drift, mutation, mating, and inbreeding depression — and to
quantify the rescue dose required to halt SC fixation via inter-BL migration.

Model summary (Sven, 2026-06-11)
────────────────────────────────
Each individual carries 4 SRK haplotypes encoded as integers; 0 = `Allele_NULL`
(non-functional), > 0 = functional allele identity (catalogue index).

Per generation (in this order):
  1. **Mutation** — each copy → NULL with probability μ.
  2. **Mating draws** — each new individual:
       a. Pick a mother uniformly at random.
       b. If mother has ≥ MIN_SC_NULL nulls AND U(0,1) < SELFING_RATE → self.
       c. Else try outcross up to MAX_ATTEMPTS times: pick a random father,
          sample one gamete (2 of 4 with double-reduction prob ALPHA).
          Reject if **any** functional father allele matches any functional
          mother allele (sporophytic SI). NULL alleles in the pollen never
          trigger rejection; NULL alleles in the stigma never reject. If no
          compatible father is found, fall back to selfing.
  3. **Gamete formation** for chosen mother + father (tetrasomic, α).
  4. **Inbreeding depression** — offspring from selfing has fitness 1 − δ.
       Implemented by weighted sampling of the N candidate offspring.
  5. **Optional inter-BL migration** — with probability m, replace one of
     the new individuals with a random individual drawn from another BL.

Outputs (per replicate × per generation × per BL):
  Tables/SRK_inheritance_trajectories.tsv
      BL, replicate, generation, scenario, mean_p_NULL, n_SI, n_pSI1, n_pSI2,
      n_pSI3, n_SC, n_individuals, mean_fitness.
  Tables/SRK_inheritance_time_to_sc.tsv
      Per-replicate first-passage time to SC_THRESHOLD per BL per scenario.

Scenarios (default ladder):
  baseline    μ=1e-4, s=0.5, δ=0.5, m=0       (current isolated state)
  rescue_low  m=0.001                          (light inter-BL flow)
  rescue_high m=0.01                           (substantial inter-BL flow)
  high_drift  Ne_scaling=0.5                   (further habitat loss)

The script keeps the per-BL N at empirical sample size unless `Ne_scaling`
is set, which multiplies each BL's N. This is the cleanest way to test
"what if Ne were smaller / larger" scenarios.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO    = Path(__file__).resolve().parent
TABLES  = REPO / "Tables"

GENO_TSV = TABLES / "Phase4" / "step23_individual_allele_genotypes_with_nulls.tsv"
OUT_TRAJ = TABLES / "Phase4" / "step24_inheritance_trajectories.tsv"
OUT_TTSC = TABLES / "Phase4" / "step24_inheritance_time_to_sc.tsv"

# ─── Defaults ─────────────────────────────────────────────────────────────────
DEFAULTS = dict(
    n_generations    = 200,
    n_replicates     = 50,
    mu               = 1e-4,
    selfing_rate     = 0.50,
    delta            = 0.50,
    alpha            = 0.10,
    max_attempts     = 20,
    min_sc_null      = 3,         # ≥ 3 NULL stigma copies → SC-like (can self)
    sc_threshold     = 0.50,      # 50 % SC frequency = "system collapsed"
    ne_scaling       = 1.0,       # multiplier on empirical BL N
    migration_rate   = 0.0,       # per-individual per-generation inter-BL replacement
    seed             = 20260611,
)

SCENARIOS = {
    "baseline":    dict(migration_rate=0.0),
    "rescue_low":  dict(migration_rate=0.001),
    "rescue_high": dict(migration_rate=0.01),
    "high_drift":  dict(ne_scaling=0.5),
}

BL_ORDER = ["BL1", "BL2", "BL3", "BL4", "BL5"]


# ─── Load empirical initial state ─────────────────────────────────────────────
def load_initial_state(geno_tsv: Path):
    df = pd.read_csv(geno_tsv, sep="\t", encoding="utf-8-sig")
    allele_cols = [c for c in df.columns if c.startswith("Allele_")]
    # Map each allele to an integer code; NULL = 0
    func_cols = [c for c in allele_cols if c != "Allele_NULL"]
    allele_to_code = {a: i + 1 for i, a in enumerate(func_cols)}
    allele_to_code["Allele_NULL"] = 0

    per_bl = {}
    for bl, sub in df.groupby("BL_inferred"):
        if bl not in BL_ORDER:
            continue
        n = len(sub)
        state = np.zeros((n, 4), dtype=np.int32)
        for i, (_, r) in enumerate(sub.iterrows()):
            # Expand counts into per-copy allele assignments
            copy_slot = 0
            for a in allele_cols:
                k = int(r[a])
                if k <= 0:
                    continue
                code = allele_to_code[a]
                for _ in range(k):
                    state[i, copy_slot] = code
                    copy_slot += 1
            assert copy_slot == 4, f"{r['Individual']}: row sums to {copy_slot}, not 4"
        per_bl[bl] = state
    return per_bl, allele_to_code


# ─── Core simulator step ──────────────────────────────────────────────────────
def sample_gamete(individual: np.ndarray, alpha: float, rng: np.random.Generator):
    """Sample 2 alleles from a tetraploid with double-reduction probability alpha."""
    if rng.random() < alpha:
        idx = rng.integers(0, 4)
        return individual[[idx, idx]]
    idx = rng.choice(4, size=2, replace=False)
    return individual[idx]


def is_compatible(pollen: np.ndarray, stigma: np.ndarray) -> bool:
    """Sporophytic SI check. Pollen = 2 alleles, stigma = 4 alleles."""
    stigma_set = set(int(x) for x in stigma if x != 0)
    for p in pollen:
        if p == 0:
            continue
        if int(p) in stigma_set:
            return False
    return True


def step_one_generation(pop: np.ndarray,
                         params: dict,
                         rng: np.random.Generator,
                         migrants: np.ndarray | None = None) -> np.ndarray:
    n = pop.shape[0]
    # 1. Mutation (functional → NULL)
    if params["mu"] > 0:
        mut_mask = rng.random((n, 4)) < params["mu"]
        if mut_mask.any():
            pop = pop.copy()
            pop[mut_mask] = 0

    # 2. Build candidate offspring + per-candidate fitness
    new_pop = np.zeros((n, 4), dtype=pop.dtype)
    fitness = np.ones(n, dtype=np.float64)

    for i in range(n):
        mother_idx = rng.integers(0, n)
        mother = pop[mother_idx]
        n_null_m = int(np.sum(mother == 0))

        selfed = False
        if (n_null_m >= params["min_sc_null"]
                and rng.random() < params["selfing_rate"]):
            father = mother
            selfed = True
        else:
            father = None
            for _ in range(params["max_attempts"]):
                cand_idx = rng.integers(0, n)
                if cand_idx == mother_idx:
                    continue
                cand = pop[cand_idx]
                p_gam = sample_gamete(cand, params["alpha"], rng)
                if is_compatible(p_gam, mother):
                    father = cand
                    father_gam = p_gam
                    break
            if father is None:
                father = mother
                father_gam = sample_gamete(mother, params["alpha"], rng)
                selfed = True
            elif father is mother:
                father_gam = sample_gamete(mother, params["alpha"], rng)

        mother_gam = sample_gamete(mother, params["alpha"], rng)
        if selfed:
            father_gam = sample_gamete(mother, params["alpha"], rng)
        new_pop[i, :2] = mother_gam
        new_pop[i, 2:] = father_gam
        if selfed:
            fitness[i] = 1.0 - params["delta"]

    # 3. Inbreeding-weighted resampling (Wright–Fisher with selection)
    if fitness.sum() > 0:
        probs = fitness / fitness.sum()
        chosen = rng.choice(n, size=n, replace=True, p=probs)
        new_pop = new_pop[chosen]

    # 4. Optional migration from other BLs
    if migrants is not None and params["migration_rate"] > 0:
        mig_mask = rng.random(n) < params["migration_rate"]
        n_mig = int(mig_mask.sum())
        if n_mig > 0:
            picks = rng.integers(0, migrants.shape[0], size=n_mig)
            new_pop[mig_mask] = migrants[picks]

    return new_pop


# ─── Per-generation summary ──────────────────────────────────────────────────
def summarise(pop: np.ndarray) -> dict:
    n = pop.shape[0]
    n_null = (pop == 0).sum(axis=1)
    return dict(
        n_individuals = n,
        mean_p_NULL   = float(np.mean((pop == 0).sum() / (4 * n))),
        n_SI          = int(np.sum(n_null == 0)),
        n_pSI1        = int(np.sum(n_null == 1)),
        n_pSI2        = int(np.sum(n_null == 2)),
        n_pSI3        = int(np.sum(n_null == 3)),
        n_SC          = int(np.sum(n_null == 4)),
    )


# ─── Run one replicate ────────────────────────────────────────────────────────
def run_replicate(initial: np.ndarray,
                  params: dict,
                  rng: np.random.Generator,
                  migrants: np.ndarray | None) -> pd.DataFrame:
    pop = initial.copy()
    # Apply Ne scaling by resampling with replacement
    n_target = max(1, int(round(pop.shape[0] * params["ne_scaling"])))
    if n_target != pop.shape[0]:
        idx = rng.integers(0, pop.shape[0], size=n_target)
        pop = pop[idx]
    rows = []
    rows.append({"generation": 0, **summarise(pop)})
    for g in range(1, params["n_generations"] + 1):
        pop = step_one_generation(pop, params, rng, migrants=migrants)
        rows.append({"generation": g, **summarise(pop)})
    return pd.DataFrame(rows)


# ─── Main driver ──────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n_generations", type=int, default=DEFAULTS["n_generations"])
    ap.add_argument("--n_replicates",  type=int, default=DEFAULTS["n_replicates"])
    ap.add_argument("--mu",            type=float, default=DEFAULTS["mu"])
    ap.add_argument("--seed",          type=int, default=DEFAULTS["seed"])
    ap.add_argument("--scenarios",     nargs="+", default=list(SCENARIOS.keys()))
    ap.add_argument("--quick", action="store_true",
                    help="Quick test: 30 gens, 5 reps, all scenarios.")
    args = ap.parse_args()

    if args.quick:
        args.n_generations = 30
        args.n_replicates  = 5

    print(f"Loading initial state from {GENO_TSV.name} ...")
    per_bl, _ = load_initial_state(GENO_TSV)
    for bl in BL_ORDER:
        if bl in per_bl:
            print(f"  {bl}: N={per_bl[bl].shape[0]}")

    base = {**DEFAULTS,
            "n_generations": args.n_generations,
            "n_replicates":  args.n_replicates,
            "mu":            args.mu}

    all_traj = []
    all_ttsc = []
    rng_master = np.random.default_rng(args.seed)

    t0 = time.time()
    for scen_name in args.scenarios:
        if scen_name not in SCENARIOS:
            print(f"  Skipping unknown scenario: {scen_name}")
            continue
        params = {**base, **SCENARIOS[scen_name]}
        print(f"\n=== Scenario: {scen_name} "
              f"(mu={params['mu']}, m={params['migration_rate']}, "
              f"Ne_scale={params['ne_scaling']}) ===")
        for bl in BL_ORDER:
            if bl not in per_bl:
                continue
            initial = per_bl[bl]
            # migration pool = all other BLs concatenated
            other = np.concatenate([per_bl[b] for b in BL_ORDER
                                    if b in per_bl and b != bl], axis=0) \
                    if params["migration_rate"] > 0 else None
            for rep in range(params["n_replicates"]):
                seed = rng_master.integers(0, 2**31)
                rng = np.random.default_rng(seed)
                df = run_replicate(initial, params, rng, other)
                df["BL"]       = bl
                df["replicate"] = rep
                df["scenario"] = scen_name
                all_traj.append(df)
                # First-passage time to SC threshold
                sc_frac = df["n_SC"] / df["n_individuals"]
                hit = sc_frac >= params["sc_threshold"]
                t_to_sc = int(df["generation"][hit].iloc[0]) if hit.any() else None
                all_ttsc.append(dict(scenario=scen_name, BL=bl, replicate=rep,
                                      time_to_sc=t_to_sc))
            print(f"  {bl}: {params['n_replicates']} reps × "
                  f"{params['n_generations']} gens done")

    traj_df = pd.concat(all_traj, ignore_index=True)
    ttsc_df = pd.DataFrame(all_ttsc)
    OUT_TRAJ.parent.mkdir(exist_ok=True)
    traj_df.to_csv(OUT_TRAJ, sep="\t", index=False)
    ttsc_df.to_csv(OUT_TTSC, sep="\t", index=False)

    elapsed = time.time() - t0
    print(f"\nTotal runtime: {elapsed:.1f} s")
    print(f"Wrote {OUT_TRAJ} ({len(traj_df)} rows)")
    print(f"Wrote {OUT_TTSC} ({len(ttsc_df)} rows)")

    # Console summary: median time-to-SC per (scenario, BL)
    print("\nMedian time-to-SC (generations) per scenario × BL:")
    summary = (ttsc_df.assign(time_to_sc=ttsc_df["time_to_sc"].fillna(args.n_generations + 1))
               .groupby(["scenario", "BL"])["time_to_sc"]
               .median().unstack())
    print(summary.to_string())


if __name__ == "__main__":
    main()
