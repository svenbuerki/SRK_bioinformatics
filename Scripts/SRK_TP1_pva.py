#!/usr/bin/env python3
"""
SRK_TP1_pva.py — Population Viability Analysis (PVA) for LEPA SI tipping
point. v1 baseline projection (no interventions).

Spec: SRK_TP1_pva_design.md — read that first.

Per-EO Monte Carlo forward simulation that couples genetic drift, SI
leakage feedback (C3 Stage 5 -> Stage 2), and a three-compartment seed
bank (non-dormant / physiologically dormant / physical dormancy).

Cohort representation
---------------------
Each seed cohort stores REPRESENTATIVE_N "representative" offspring
genotypes (4-allele tuples) plus a total seed-count `n`. The
representatives are a random sample of the offspring produced this year;
storing the full cohort (millions of seeds) is infeasible. Mortality
shrinks `n` proportionally each year; the representatives are unchanged
(no drift in storage). Germination samples with replacement from the
representatives, scaled by `n` × germination_rate.

This sample-based approach avoids the analytical class-distribution
computation in earlier drafts, which was both slow (O(k^3) per call) and
incorrect for tetraploid mating (treated 4 alleles per offspring as
independent draws, breaking tetraploid genotype structure and collapsing
AAAA frequency artificially).

Outputs
-------
Tables/SRK_TP1_pva_checkpoint_summary.tsv
Tables/SRK_TP1_pva_allele_extinction.tsv
"""

from __future__ import annotations

import csv
import math
import random
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parent
TABLES = REPO / "Tables"

METRICS_TSV = TABLES / "SRK_EO_allele_richness.tsv"
GENOTYPES_TSV = TABLES / "SRK_individual_allele_genotypes.tsv"
METADATA_CSV = TABLES / "sampling_metadata.csv"
FECUNDITY_CSV = Path(
    "/Users/sven/Documents/Current_projects/NSF24-543_Self-incompatibility"
    "/Brainstorm_NSF/Data/Data_aggregation_by_location_fecundity.csv"
)

OUT_SUMMARY = TABLES / "SRK_TP1_pva_checkpoint_summary.tsv"
OUT_ALLELES = TABLES / "SRK_TP1_pva_allele_extinction.tsv"

# ----- Constants --------------------------------------------------------------

SUBCODE_POOL = {"27pv": "27", "27rt": "27"}
GENOTYPE_CLASSES = ["AAAA", "AAAB", "AABB", "AABC", "ABCD"]
P_COLLAPSE = 0.05               # P_compat collapse threshold
REPRESENTATIVE_N = 200          # cohort representatives sampled per new cohort

COMPARTMENT_RATES = {
    # name: (germ_rate, mortality_rate, max_age_yr)
    "ND": dict(germ=0.80, mort=0.10, max_age=5),
    "PD": dict(germ=0.15, mort=0.06, max_age=20),
    "PH": dict(germ=0.02, mort=0.015, max_age=60),
}

JUV_SURVIVAL = 0.50
K_MULT_MEAN = 1.5
K_MULT_SIGMA = 0.30
SIGMA_ENV = 0.30
L_CAP = 0.50

N_REPS = 1000
HORIZON = 50
CHECKPOINTS = [10, 20, 50]
SEED = 20260526

PCTILES = [10, 50, 90]
FOCAL_EOS = ["18", "25", "27", "67", "70", "76"]

DORMANCY_SPLIT = {
    # EO: (f_ND, f_PD, f_PH, f_NV)  — values from figures/dormancy_stacked_barplot.jpeg
    "18": (0.60, 0.13, 0.27, 0.00),
    "25": (0.67, 0.08, 0.23, 0.02),
    "27": (0.25, 0.40, 0.20, 0.15),
    "67": (0.68, 0.08, 0.22, 0.02),
    "70": (0.46, 0.15, 0.36, 0.03),
    "76": (0.93, 0.00, 0.07, 0.00),
}


# ============================================================================
# Data loading
# ============================================================================


def load_metadata() -> dict[str, str]:
    out: dict[str, str] = {}
    with METADATA_CSV.open(encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            if (row.get("Ingroup") or "").strip() != "1":
                continue
            sid = row["SampleID"].strip()
            eo = (row.get("EO_w_sub") or "").strip()
            if not eo:
                continue
            out[sid] = SUBCODE_POOL.get(eo, eo)
    return out


def infer_tetraploid_genotype(raw: dict[str, int]) -> tuple[str, ...]:
    """Mirror SRK_zygosity_from_genotype.R: fill to 4 copies; return tuple
    of 4 allele ids (with repetition)."""
    present = sorted([(a, c) for a, c in raw.items() if c > 0],
                     key=lambda x: -x[1])
    n = len(present)
    if n == 0:
        return ()
    if n == 1:
        a, _ = present[0]
        return (a, a, a, a)
    if n == 2:
        (a, ca), (b, cb) = present
        if ca >= 2 * cb:
            return (a, a, a, b)
        return (a, a, b, b)
    if n == 3:
        (a, _), (b, _), (c, _) = present
        return (a, a, b, c)
    (a, _), (b, _), (c, _), (d, _) = present[:4]
    return (a, b, c, d)


def load_per_eo_individuals(meta: dict[str, str]) -> dict[str, list[tuple[str, ...]]]:
    by_eo: dict[str, list[tuple[str, ...]]] = {eo: [] for eo in FOCAL_EOS}
    with GENOTYPES_TSV.open(encoding="utf-8-sig", newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        alleles = header[1:]
        for row in rdr:
            sid = row[0]
            eo = meta.get(sid)
            if eo not in by_eo:
                continue
            raw = {a: int(v) for a, v in zip(alleles, row[1:])}
            geno = infer_tetraploid_genotype(raw)
            if geno:
                by_eo[eo].append(geno)
    return by_eo


def load_eo_census_and_fecundity() -> dict[str, dict]:
    out: dict[str, dict] = defaultdict(lambda: {
        "census_N": 0, "fec_per_plant": 0.0, "n_germ_obs": 0,
        "seeds_banked": 0.0,
    })
    with FECUNDITY_CSV.open(encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            eo_id = row["EOID"].strip()
            if eo_id.startswith("EO"):
                eo_id = eo_id[2:]
            if eo_id not in FOCAL_EOS:
                continue
            d = out[eo_id]
            try:
                d["census_N"] += int(row["total_OrganismQuantityFertile"] or 0)
                d["census_N"] += int(row["total_OrganismQuantityVegetative"] or 0)
            except (ValueError, TypeError):
                pass
            try:
                d["fec_per_plant"] += float(row["mean_number_of_seeds_per_germplasm"])
                d["n_germ_obs"] += 1
            except (ValueError, TypeError):
                pass
            try:
                d["seeds_banked"] += float(row["total_number_of_seeds_banked"])
            except (ValueError, TypeError):
                pass
    for eo, d in out.items():
        d["fec_per_plant"] = (
            d["fec_per_plant"] / d["n_germ_obs"] if d["n_germ_obs"] else 300.0
        )
        del d["n_germ_obs"]
    return dict(out)


# ============================================================================
# Core simulation primitives
# ============================================================================


def class_of(genotype: tuple[str, ...]) -> str:
    counts = sorted(Counter(genotype).values(), reverse=True)
    if counts == [4]:
        return "AAAA"
    if counts == [3, 1]:
        return "AAAB"
    if counts == [2, 2]:
        return "AABB"
    if counts == [2, 1, 1]:
        return "AABC"
    return "ABCD"


def gamete_freqs_from_adults(adults: list[tuple[str, ...]]) -> dict[str, float]:
    copies: Counter = Counter()
    for g in adults:
        for a in g:
            copies[a] += 1
    total = sum(copies.values())
    if total == 0:
        return {}
    return {a: c / total for a, c in copies.items()}


def p_compat_tetraploid(freqs: dict[str, float]) -> float:
    """Exact tetraploid sporophytic SI compatibility under multinomial sampling.
    Sums over ordered pollen pairs; stigma must contain neither.
    Matches SRK_TP1_compatibility_metrics.py."""
    ps = list(freqs.values())
    n = len(ps)
    total = 0.0
    for i in range(n):
        pi = ps[i]
        if pi == 0:
            continue
        total += pi * pi * (1 - pi) ** 4
        for j in range(n):
            if i == j:
                continue
            pj = ps[j]
            if pj == 0:
                continue
            rem = max(0.0, 1 - pi - pj)
            total += pi * pj * (rem ** 4)
    return total


def sample_outcross_offspring(
    parents: list[tuple[str, ...]], rng: random.Random
) -> tuple[str, ...]:
    """One outcrossed offspring: 2 alleles from each of 2 distinct parents,
    drawn without replacement within each parent (tetrasomic random pairing).
    SI rejection is NOT enforced at the per-pair level here — instead the
    expected count of viable outcrossed seeds is reduced upstream by
    multiplying by P_compat. This is the design-doc simplification (a
    population-level Bernoulli on cross success, not a per-pair filter)."""
    p1, p2 = rng.sample(parents, 2)
    g1 = rng.sample(p1, 2)
    g2 = rng.sample(p2, 2)
    return tuple(g1 + g2)


def sample_selfed_offspring(
    parent: tuple[str, ...], rng: random.Random
) -> tuple[str, ...]:
    """One selfed offspring under tetrasomic inheritance: two independent
    gametes from the same parent."""
    g1 = rng.sample(parent, 2)
    g2 = rng.sample(parent, 2)
    return tuple(g1 + g2)


# ============================================================================
# Seed cohort
# ============================================================================


@dataclass
class SeedCohort:
    age: int
    n: float                                     # remaining seed count
    representatives: list[tuple[str, ...]]       # REPRESENTATIVE_N sample of cohort

    def alleles_present(self) -> set[str]:
        s: set[str] = set()
        for ind in self.representatives:
            s.update(ind)
        return s


@dataclass
class CompartmentState:
    name: str
    cohorts: list[SeedCohort] = field(default_factory=list)

    def alleles_present(self) -> set[str]:
        s: set[str] = set()
        for coh in self.cohorts:
            s.update(coh.alleles_present())
        return s

    def total_seeds(self) -> float:
        return sum(c.n for c in self.cohorts)

    def age_one_year(self) -> None:
        rates = COMPARTMENT_RATES[self.name]
        out = []
        for coh in self.cohorts:
            coh.age += 1
            if coh.age > rates["max_age"]:
                continue
            coh.n *= (1.0 - rates["mort"])
            if coh.n < 0.5:
                continue
            out.append(coh)
        self.cohorts = out

    def germinate(self, rng: random.Random) -> list[tuple[str, ...]]:
        rates = COMPARTMENT_RATES[self.name]
        gr = rates["germ"]
        recruits: list[tuple[str, ...]] = []
        for coh in self.cohorts:
            n_germ = int(round(coh.n * gr))
            if n_germ <= 0 or not coh.representatives:
                continue
            recruits.extend(
                rng.choices(coh.representatives, k=n_germ)
            )
            coh.n -= n_germ
        self.cohorts = [c for c in self.cohorts if c.n >= 0.5]
        return recruits


# ============================================================================
# EO state
# ============================================================================


@dataclass
class EOState:
    eo: str
    adults: list[tuple[str, ...]]
    K: float
    fec_per_plant: float
    f_ND: float
    f_PD: float
    f_PH: float
    f_NV: float
    L: float
    nd: CompartmentState = field(default_factory=lambda: CompartmentState("ND"))
    pd: CompartmentState = field(default_factory=lambda: CompartmentState("PD"))
    ph: CompartmentState = field(default_factory=lambda: CompartmentState("PH"))
    extinction_year: dict[str, int] = field(default_factory=dict)


def initial_state(
    eo: str,
    individuals: list[tuple[str, ...]],
    census: dict,
    rng: random.Random,
) -> EOState:
    fec = census["fec_per_plant"] if census["fec_per_plant"] > 0 else 300.0
    census_N = max(census["census_N"], len(individuals))
    f_ND, f_PD, f_PH, f_NV = DORMANCY_SPLIT[eo]
    viable = f_ND + f_PD + f_PH
    if viable > 0:
        f_ND_n, f_PD_n, f_PH_n = f_ND / viable, f_PD / viable, f_PH / viable
    else:
        f_ND_n = f_PD_n = f_PH_n = 0.0
    log_K = rng.gauss(math.log(census_N * K_MULT_MEAN), K_MULT_SIGMA)
    K_val = max(census_N, math.exp(log_K))
    if census_N <= len(individuals):
        adults_init = rng.sample(individuals, census_N)
    else:
        adults_init = [rng.choice(individuals) for _ in range(census_N)]
    n_aaaa = sum(1 for g in adults_init if class_of(g) == "AAAA")
    prop_aaaa = n_aaaa / len(adults_init) if adults_init else 0.0
    L0 = min(L_CAP, prop_aaaa / 3.5)
    state = EOState(
        eo=eo, adults=adults_init, K=K_val, fec_per_plant=fec,
        f_ND=f_ND_n, f_PD=f_PD_n, f_PH=f_PH_n, f_NV=f_NV, L=L0,
    )
    # Initial bank: partition standing seeds into the three compartments,
    # representatives bootstrapped from the standing adults (steady-state).
    init_bank = max(0.0, float(census.get("seeds_banked", 0.0)))
    if init_bank > 0 and adults_init:
        rep_seeds = _sample_representatives(adults_init, REPRESENTATIVE_N, rng)
        for compartment, share in [
            (state.nd, f_ND), (state.pd, f_PD), (state.ph, f_PH),
        ]:
            n_in_pool = init_bank * share
            if n_in_pool < 0.5:
                continue
            # Use the same representative seed set for all three compartments
            compartment.cohorts.append(SeedCohort(
                age=1, n=n_in_pool, representatives=list(rep_seeds),
            ))
    return state


def _sample_representatives(
    individuals: list[tuple[str, ...]], n: int, rng: random.Random
) -> list[tuple[str, ...]]:
    """Sample n offspring representatives from a pool of individuals (treated
    as parents). Used for the initial seed bank (steady-state assumption)."""
    out: list[tuple[str, ...]] = []
    n_ind = len(individuals)
    if n_ind == 0:
        return out
    for _ in range(n):
        if n_ind == 1:
            out.append(sample_selfed_offspring(individuals[0], rng))
            continue
        # default to outcrossed for the initial bank
        out.append(sample_outcross_offspring(individuals, rng))
    return out


# ============================================================================
# Annual transition
# ============================================================================


def collapsed(state: EOState) -> bool:
    if not state.adults:
        if (state.nd.total_seeds() + state.pd.total_seeds()
                + state.ph.total_seeds()) < 1.0:
            return True
        return False
    freqs = gamete_freqs_from_adults(state.adults)
    if not freqs:
        return True
    return p_compat_tetraploid(freqs) < P_COLLAPSE


def step_year(state: EOState, year: int, rng: random.Random) -> None:
    if state.adults:
        freqs = gamete_freqs_from_adults(state.adults)
        p_compat = p_compat_tetraploid(freqs)
        N_adult = len(state.adults)

        seeds_outcr = N_adult * state.fec_per_plant * p_compat * (1 - state.L)
        seeds_self = N_adult * state.fec_per_plant * state.L
        seeds_total = seeds_outcr + seeds_self
        viable_seeds = seeds_total * (1 - state.f_NV)

        if viable_seeds >= 1:
            # Generate representatives for the new cohort
            outcr_share = seeds_outcr / max(seeds_total, 1e-9)
            n_outcr = int(round(REPRESENTATIVE_N * outcr_share))
            n_self = REPRESENTATIVE_N - n_outcr
            reps: list[tuple[str, ...]] = []
            if N_adult >= 2:
                for _ in range(n_outcr):
                    reps.append(sample_outcross_offspring(state.adults, rng))
            else:
                # Can't outcross with only one adult; all reproduction is selfing
                n_self = REPRESENTATIVE_N
            for _ in range(n_self):
                parent = rng.choice(state.adults)
                reps.append(sample_selfed_offspring(parent, rng))

            # Partition viable_seeds into compartments
            for comp, share in [
                (state.nd, state.f_ND),
                (state.pd, state.f_PD),
                (state.ph, state.f_PH),
            ]:
                n_in_pool = viable_seeds * share
                if n_in_pool >= 0.5 and share > 0:
                    comp.cohorts.append(SeedCohort(
                        age=0, n=n_in_pool, representatives=list(reps),
                    ))

    # Age and germinate (mortality first so age-0 cohorts don't lose seeds)
    for comp in (state.nd, state.pd, state.ph):
        comp.age_one_year()

    juveniles: list[tuple[str, ...]] = []
    for comp in (state.nd, state.pd, state.ph):
        juveniles.extend(comp.germinate(rng))

    # Juvenile survival
    n_alive = int(round(len(juveniles) * JUV_SURVIVAL))
    if n_alive <= 0:
        new_adults: list[tuple[str, ...]] = []
    elif n_alive < len(juveniles):
        new_adults = rng.sample(juveniles, n_alive)
    else:
        new_adults = list(juveniles)

    # Environmental noise on annual K
    K_t = state.K * math.exp(rng.gauss(0, SIGMA_ENV))
    cap = max(0, int(round(K_t)))
    if len(new_adults) > cap:
        new_adults = rng.sample(new_adults, cap)

    state.adults = new_adults

    # Update L
    if state.adults:
        n_aaaa = sum(1 for g in state.adults if class_of(g) == "AAAA")
        prop_aaaa = n_aaaa / len(state.adults)
        state.L = min(L_CAP, prop_aaaa / 3.5)


# ============================================================================
# Replicate runner
# ============================================================================


def run_replicate(
    eo: str,
    base_individuals: list[tuple[str, ...]],
    census: dict,
    rng: random.Random,
) -> dict:
    state = initial_state(eo, base_individuals, census, rng)
    ever_present: set[str] = set()
    for g in state.adults:
        ever_present.update(g)
    for comp in (state.nd, state.pd, state.ph):
        ever_present.update(comp.alleles_present())

    snapshots = {}
    collapse_year: int | None = None
    for year in range(1, HORIZON + 1):
        step_year(state, year, rng)
        present_now: set[str] = set()
        for g in state.adults:
            present_now.update(g)
        for comp in (state.nd, state.pd, state.ph):
            present_now.update(comp.alleles_present())
        for a in ever_present - present_now - set(state.extinction_year):
            state.extinction_year[a] = year
        if collapse_year is None and collapsed(state):
            collapse_year = year
        if year in CHECKPOINTS:
            freqs = gamete_freqs_from_adults(state.adults)
            p_compat = p_compat_tetraploid(freqs) if freqs else 0.0
            snapshots[year] = {
                "N_adult": len(state.adults),
                "k_alleles": len(present_now),
                "P_compat": p_compat,
                "collapsed": collapse_year is not None and collapse_year <= year,
                "L": state.L,
            }
    for cp in CHECKPOINTS:
        snapshots.setdefault(cp, {
            "N_adult": 0, "k_alleles": 0, "P_compat": 0.0,
            "collapsed": True, "L": 0.0,
        })
    return {
        "snapshots": snapshots,
        "collapse_year": collapse_year,
        "extinction_year": dict(state.extinction_year),
        "ever_present": ever_present,
    }


def percentile(values: list[float], q: float) -> float:
    if not values:
        return float("nan")
    s = sorted(values)
    k = (len(s) - 1) * q / 100.0
    lo, hi = int(math.floor(k)), int(math.ceil(k))
    if lo == hi:
        return s[lo]
    return s[lo] + (s[hi] - s[lo]) * (k - lo)


def run_eo(
    eo: str,
    base_individuals: list[tuple[str, ...]],
    census: dict,
    n_reps: int,
    rng: random.Random,
) -> tuple[list[dict], list[dict]]:
    reps = [run_replicate(eo, base_individuals, census, rng) for _ in range(n_reps)]

    summary_rows = []
    for cp in CHECKPOINTS:
        snaps = [r["snapshots"][cp] for r in reps]
        N_vals = [s["N_adult"] for s in snaps]
        k_vals = [s["k_alleles"] for s in snaps]
        pc_vals = [s["P_compat"] for s in snaps]
        L_vals = [s["L"] for s in snaps]
        p_collapsed = sum(1 for s in snaps if s["collapsed"]) / max(1, len(snaps))
        row = {
            "EO": f"EO{eo}",
            "checkpoint_year": cp,
            "P_collapsed": f"{p_collapsed:.3f}",
        }
        for q in PCTILES:
            row[f"N_adult_p{q}"] = f"{percentile(N_vals, q):.0f}"
            row[f"k_alleles_p{q}"] = f"{percentile(k_vals, q):.1f}"
            row[f"P_compat_p{q}"] = f"{percentile(pc_vals, q):.3f}"
            row[f"L_p{q}"] = f"{percentile(L_vals, q):.3f}"
        summary_rows.append(row)

    all_alleles: set[str] = set()
    for r in reps:
        all_alleles.update(r["ever_present"])
    allele_rows = []
    for a in sorted(all_alleles):
        years = [r["extinction_year"].get(a) for r in reps if a in r["ever_present"]]
        n_present = len(years)
        loss_years = [y for y in years if y is not None]
        p_lost_by = {cp: sum(1 for y in loss_years if y <= cp) / max(1, n_present)
                     for cp in CHECKPOINTS}
        allele_rows.append({
            "EO": f"EO{eo}",
            "allele": a,
            "n_reps_present_initially": n_present,
            "P_lost_by_y10": f"{p_lost_by[10]:.3f}",
            "P_lost_by_y20": f"{p_lost_by[20]:.3f}",
            "P_lost_by_y50": f"{p_lost_by[50]:.3f}",
            "median_loss_year": (
                f"{percentile(loss_years, 50):.1f}" if loss_years else "NA"
            ),
        })
    return summary_rows, allele_rows


# ============================================================================
# Entry point
# ============================================================================


def main() -> None:
    print("Loading inputs ...")
    meta = load_metadata()
    individuals_by_eo = load_per_eo_individuals(meta)
    census_by_eo = load_eo_census_and_fecundity()

    rng = random.Random(SEED)
    all_summary = []
    all_alleles = []

    for eo in FOCAL_EOS:
        inds = individuals_by_eo[eo]
        census = census_by_eo.get(eo, {"census_N": len(inds), "fec_per_plant": 300.0,
                                        "seeds_banked": 0.0})
        if not inds:
            print(f"  EO{eo}: no genotyped individuals — skipping")
            continue
        print(f"  EO{eo}: census_N={census['census_N']}, "
              f"fec/plant={census['fec_per_plant']:.0f}, "
              f"seeds_banked={census['seeds_banked']:.0f}, "
              f"n_reps={N_REPS}")
        summary, alleles = run_eo(eo, inds, census, N_REPS, rng)
        all_summary.extend(summary)
        all_alleles.extend(alleles)

    with OUT_SUMMARY.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_summary[0].keys()),
                            delimiter="\t")
        w.writeheader()
        w.writerows(all_summary)
    with OUT_ALLELES.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_alleles[0].keys()),
                            delimiter="\t")
        w.writeheader()
        w.writerows(all_alleles)

    print(f"\nWrote {OUT_SUMMARY.relative_to(REPO)} ({len(all_summary)} rows)")
    print(f"Wrote {OUT_ALLELES.relative_to(REPO)} ({len(all_alleles)} rows)")

    cols = ["EO", "checkpoint_year", "P_collapsed",
            "N_adult_p50", "k_alleles_p50", "P_compat_p50", "L_p50"]
    col_w = {k: max(len(k), max(len(str(r[k])) for r in all_summary)) for k in cols}
    line = "  ".join(k.ljust(col_w[k]) for k in cols)
    print(f"\nCheckpoint summary (medians):\n{line}\n{'-' * len(line)}")
    for r in all_summary:
        print("  ".join(str(r[k]).ljust(col_w[k]) for k in cols))


if __name__ == "__main__":
    main()
