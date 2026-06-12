#!/usr/bin/env python3
"""
SRK_TP1_pva_prediction.py — Population Viability Analysis v1.6 PREDICTION run.

Spec: SRK_TP1_pva_design.md.

Single-scenario forward projection with **empirically-grounded per-EO
parameters** derived from the IDFG Element Occurrence Records (kept in
EO_data_sensitive_data/). Replaces the generic σ_env = 0.30 + 5 % catastrophe
defaults of v1/v1.5 with the actual variability observed in each EO's
historical census record (1989-2023, depending on EO).

This is the "what will actually happen" run, distinct from v1.5 which is the
"sensitivity to assumptions" sweep.

The prediction scenario combines all empirically-justified mechanisms:
  * per-EO σ_env from HIP-transect variability (capped at 1.0)
  * per-EO catastrophe rate (from observed extreme-low years)
  * area-bounded K (where habitat area actually constrains population)
  * per-mating SI rejection — proper NFDS rare-allele advantage
  * SC + inbreeding depression for EO76 (5 SI-escape candidates from Step 12c)

NOT in this run:
  * K decline trend (no empirical evidence for systematic K decline)
  * Other speculative stressors

Inputs
------
EO_data_sensitive_data/extracted_pva_inputs.tsv  ← required for prediction

Outputs
-------
Tables/SRK_TP1_pva_prediction_summary.tsv     (per EO × checkpoint)
Tables/SRK_TP1_pva_prediction_alleles.tsv     (per EO × allele)

Sensitive-data handling
-----------------------
This script READS sensitive census-derived parameters from
EO_data_sensitive_data/ but writes ONLY aggregated/derived PVA metrics
(probabilities, percentiles, allele extinction rates) to Tables/. No raw
historical census counts are included in any published output.
"""

from __future__ import annotations

import csv
import math
import multiprocessing as mp
import random
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parent
TABLES = REPO / "Tables"
SENSITIVE = REPO / "EO_data_sensitive_data"

METRICS_TSV = TABLES / "SRK_EO_allele_richness.tsv"
GENOTYPES_TSV = TABLES / "SRK_individual_allele_genotypes.tsv"
METADATA_CSV = TABLES / "sampling_metadata.csv"
AREAS_CSV = TABLES / "EO_group_BL_summary.csv"
FECUNDITY_CSV = Path(
    "/Users/sven/Documents/Current_projects/NSF24-543_Self-incompatibility"
    "/Brainstorm_NSF/Data/Data_aggregation_by_location_fecundity.csv"
)
PVA_INPUTS_TSV = SENSITIVE / "extracted_pva_inputs.tsv"

OUT_SUMMARY = TABLES / "SRK_TP1_pva_prediction_summary.tsv"
OUT_ALLELES = TABLES / "SRK_TP1_pva_prediction_alleles.tsv"

# ----- Constants (kept in sync with SRK_TP1_pva.py) -----------------------

SUBCODE_POOL = {"27pv": "27", "27rt": "27"}
GENOTYPE_CLASSES = ["AAAA", "AAAB", "AABB", "AABC", "ABCD"]
P_COLLAPSE = 0.05
REPRESENTATIVE_N = 200

COMPARTMENT_RATES = {
    "ND": dict(germ=0.80, mort=0.10, max_age=5),
    "PD": dict(germ=0.15, mort=0.06, max_age=20),
    "PH": dict(germ=0.02, mort=0.015, max_age=60),
}

JUV_SURVIVAL = 0.50
K_MULT_MEAN = 1.5
K_MULT_SIGMA = 0.30
L_CAP = 0.50
SIGMA_ENV_FALLBACK = 0.30                 # default if no empirical value
SIGMA_ENV_CAP = 1.0                       # observed 1.18 capped for stability
CATASTROPHE_RECRUIT_MULT = 0.05
SI_MAX_ATTEMPTS = 100

# Area-cap (same as v1.5 area_capped)
AREA_MAX_DENSITY_PER_HA = 5000.0

# SC + ID
ID_SEVERITY = 0.5

N_REPS = 1000
HORIZON = 50
CHECKPOINTS = [10, 20, 50]
SEED = 20260527
PCTILES = [10, 50, 90]

FOCAL_EOS = ["18", "25", "27", "67", "70", "76"]

DORMANCY_SPLIT = {
    "18": (0.60, 0.13, 0.27, 0.00),
    "25": (0.67, 0.08, 0.23, 0.02),
    "27": (0.25, 0.40, 0.20, 0.15),
    "67": (0.68, 0.08, 0.22, 0.02),
    "70": (0.46, 0.15, 0.36, 0.03),
    "76": (0.93, 0.00, 0.07, 0.00),
}

SC_INDIVIDUALS = {eo: 0 for eo in FOCAL_EOS}
SC_INDIVIDUALS["76"] = 5


# ============================================================================
# Sensitive inputs loader
# ============================================================================


@dataclass(frozen=True)
class EmpiricalParams:
    sigma_env: float           # from HIP-transect log SD
    catastrophe_rate: float    # observed frequency of extreme-low years
    historical_extinction: bool


def load_empirical_inputs() -> dict[str, EmpiricalParams]:
    """Read per-EO sigma_env and catastrophe_rate from the sensitive TSV.
    Returns empty dict if file missing (script will fall back to defaults)."""
    out: dict[str, EmpiricalParams] = {}
    if not PVA_INPUTS_TSV.exists():
        print(f"WARNING: {PVA_INPUTS_TSV} missing — using fallback defaults",
              file=sys.stderr)
        return out
    with PVA_INPUTS_TSV.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            eo = row["EO"].strip()
            sigma = min(SIGMA_ENV_CAP, float(row["sigma_env_log"]))
            cat_years = row["catastrophe_years"].strip()
            # Catastrophe rate = #catastrophe-years / span of observations.
            # Cap at 0.20 to avoid one historically catastrophic period
            # dominating future projections.
            n_cat = len([y for y in cat_years.replace(",", " ").split() if y.isdigit()])
            n_hip = int(row["HIP_n_obs"]) if row["HIP_n_obs"].isdigit() else 5
            cat_rate = min(0.20, n_cat / max(n_hip, 1))
            ext = row["historical_extinction"].strip().lower().startswith("yes")
            out[eo] = EmpiricalParams(sigma_env=sigma,
                                       catastrophe_rate=cat_rate,
                                       historical_extinction=ext)
    return out


# ============================================================================
# Data loading (shared with v1/v1.5)
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


def load_eo_areas() -> dict[str, float]:
    out: dict[str, float] = {eo: 0.0 for eo in FOCAL_EOS}
    with AREAS_CSV.open(encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            eo_raw = row["EO"].strip()
            try:
                area = float(row["Area_ha"])
            except (ValueError, TypeError):
                continue
            parts = [p.strip().lstrip("EO") for p in eo_raw.replace(",", ";").split(";")]
            parts = [p for p in parts if p]
            if not parts:
                continue
            share = area / len(parts)
            for p in parts:
                if p in out:
                    out[p] += share
    return out


# ============================================================================
# Core simulation primitives (same as v1.5)
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


def sample_outcross_offspring_si(
    parents: list[tuple[str, ...]], rng: random.Random
) -> tuple[str, ...] | None:
    """Per-mating SI rejection (always used in prediction run)."""
    for _ in range(SI_MAX_ATTEMPTS):
        p1, p2 = rng.sample(parents, 2)
        g2 = rng.sample(p2, 2)
        stigma_set = set(p1)
        if g2[0] not in stigma_set and g2[1] not in stigma_set:
            g1 = rng.sample(p1, 2)
            return tuple(g1 + g2)
    return None


def sample_selfed_offspring(
    parent: tuple[str, ...], rng: random.Random
) -> tuple[str, ...]:
    g1 = rng.sample(parent, 2)
    g2 = rng.sample(parent, 2)
    return tuple(g1 + g2)


# ============================================================================
# Seed cohort + compartment
# ============================================================================


@dataclass
class SeedCohort:
    age: int
    n: float
    representatives: list[tuple[str, ...]]
    f_sc: float = 0.0          # cohort-average ID survival multiplier (re-purposed field)

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
        for c in self.cohorts:
            s.update(c.alleles_present())
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

    def germinate(self, rng: random.Random) -> list[tuple[tuple[str, ...], float]]:
        rates = COMPARTMENT_RATES[self.name]
        gr = rates["germ"]
        recruits: list[tuple[tuple[str, ...], float]] = []
        for coh in self.cohorts:
            n_germ = int(round(coh.n * gr))
            if n_germ <= 0 or not coh.representatives:
                continue
            picks = rng.choices(coh.representatives, k=n_germ)
            id_mult = 1.0 - coh.f_sc        # f_sc field holds the ID-derived survival reduction
            recruits.extend((p, id_mult) for p in picks)
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
    K0: float
    fec_per_plant: float
    f_ND: float
    f_PD: float
    f_PH: float
    f_NV: float
    L: float
    f_SC: float
    sigma_env: float
    catastrophe_rate: float
    nd: CompartmentState = field(default_factory=lambda: CompartmentState("ND"))
    pd: CompartmentState = field(default_factory=lambda: CompartmentState("PD"))
    ph: CompartmentState = field(default_factory=lambda: CompartmentState("PH"))
    extinction_year: dict[str, int] = field(default_factory=dict)


def initial_state(
    eo: str,
    individuals: list[tuple[str, ...]],
    census: dict,
    area_ha: float,
    empirical: EmpiricalParams | None,
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

    # K: area-bounded (the empirically-supported constraint)
    if area_ha > 0.001:
        K_area = area_ha * AREA_MAX_DENSITY_PER_HA
        K_demo = census_N * K_MULT_MEAN
        K_val = min(K_area, K_demo)
    else:
        log_K = rng.gauss(math.log(census_N * K_MULT_MEAN), K_MULT_SIGMA)
        K_val = max(census_N, math.exp(log_K))

    if census_N <= len(individuals):
        adults_init = rng.sample(individuals, census_N)
    else:
        adults_init = [rng.choice(individuals) for _ in range(census_N)]

    n_aaaa = sum(1 for g in adults_init if class_of(g) == "AAAA")
    prop_aaaa = n_aaaa / len(adults_init) if adults_init else 0.0
    L0 = min(L_CAP, prop_aaaa / 3.5)

    f_SC = SC_INDIVIDUALS.get(eo, 0) / max(census_N, 1)

    # Empirical environmental stochasticity (per-EO)
    if empirical is not None:
        sigma_env = empirical.sigma_env
        cat_rate = empirical.catastrophe_rate
    else:
        sigma_env = SIGMA_ENV_FALLBACK
        cat_rate = 0.0

    state = EOState(
        eo=eo, adults=adults_init, K0=K_val, fec_per_plant=fec,
        f_ND=f_ND_n, f_PD=f_PD_n, f_PH=f_PH_n, f_NV=f_NV, L=L0, f_SC=f_SC,
        sigma_env=sigma_env, catastrophe_rate=cat_rate,
    )

    init_bank = max(0.0, float(census.get("seeds_banked", 0.0)))
    if init_bank > 0 and adults_init:
        rep_seeds = _initial_bank_representatives(adults_init, REPRESENTATIVE_N, rng)
        for compartment, share in [
            (state.nd, f_ND), (state.pd, f_PD), (state.ph, f_PH),
        ]:
            n_in_pool = init_bank * share
            if n_in_pool < 0.5:
                continue
            compartment.cohorts.append(SeedCohort(
                age=1, n=n_in_pool, representatives=list(rep_seeds), f_sc=0.0,
            ))
    return state


def _initial_bank_representatives(
    individuals: list[tuple[str, ...]], n: int, rng: random.Random
) -> list[tuple[str, ...]]:
    out: list[tuple[str, ...]] = []
    if not individuals:
        return out
    for _ in range(n):
        if len(individuals) == 1:
            out.append(sample_selfed_offspring(individuals[0], rng))
            continue
        # Initial bank uses SI-respecting outcross for consistency with v1.6
        child = sample_outcross_offspring_si(individuals, rng)
        if child is None:
            out.append(sample_selfed_offspring(rng.choice(individuals), rng))
        else:
            out.append(child)
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
    K_t = state.K0 * math.exp(rng.gauss(0, state.sigma_env))
    K_t_int = max(0, int(round(K_t)))
    is_catastrophe = (
        state.catastrophe_rate > 0 and rng.random() < state.catastrophe_rate
    )

    if state.adults:
        freqs = gamete_freqs_from_adults(state.adults)
        p_compat = p_compat_tetraploid(freqs)
        N_adult = len(state.adults)
        sc_share = state.f_SC

        seeds_sc = N_adult * sc_share * state.fec_per_plant
        seeds_si_outcr = N_adult * (1 - sc_share) * state.fec_per_plant * p_compat * (1 - state.L)
        seeds_si_self = N_adult * (1 - sc_share) * state.fec_per_plant * state.L
        seeds_total = seeds_sc + seeds_si_outcr + seeds_si_self
        viable_seeds = seeds_total * (1 - state.f_NV)

        if viable_seeds >= 1 and N_adult >= 1:
            outcr_share = seeds_si_outcr / max(seeds_total, 1e-9)
            si_self_share = seeds_si_self / max(seeds_total, 1e-9)
            sc_repr_share = seeds_sc / max(seeds_total, 1e-9)

            n_outcr = int(round(REPRESENTATIVE_N * outcr_share))
            n_si_self = int(round(REPRESENTATIVE_N * si_self_share))
            n_sc = REPRESENTATIVE_N - n_outcr - n_si_self

            reps: list[tuple[str, ...]] = []
            if n_outcr > 0 and N_adult >= 2:
                succeeded = 0
                while succeeded < n_outcr:
                    child = sample_outcross_offspring_si(state.adults, rng)
                    if child is None:
                        break
                    reps.append(child)
                    succeeded += 1
                n_si_self += (n_outcr - succeeded)
            elif N_adult == 1:
                n_si_self += n_outcr

            for _ in range(n_si_self):
                parent = rng.choice(state.adults)
                reps.append(sample_selfed_offspring(parent, rng))
            for _ in range(n_sc):
                parent = rng.choice(state.adults)
                reps.append(sample_selfed_offspring(parent, rng))

            cohort_f_sc = (n_sc / REPRESENTATIVE_N) if REPRESENTATIVE_N else 0.0

            for comp, share in [
                (state.nd, state.f_ND),
                (state.pd, state.f_PD),
                (state.ph, state.f_PH),
            ]:
                n_in_pool = viable_seeds * share
                if n_in_pool >= 0.5 and share > 0:
                    comp.cohorts.append(SeedCohort(
                        age=0, n=n_in_pool, representatives=list(reps),
                        f_sc=cohort_f_sc * ID_SEVERITY,
                    ))

    for comp in (state.nd, state.pd, state.ph):
        comp.age_one_year()

    juveniles_with_mult: list[tuple[tuple[str, ...], float]] = []
    for comp in (state.nd, state.pd, state.ph):
        juveniles_with_mult.extend(comp.germinate(rng))

    survivors: list[tuple[str, ...]] = []
    base_surv = JUV_SURVIVAL
    if is_catastrophe:
        base_surv *= CATASTROPHE_RECRUIT_MULT
    for ind, mult in juveniles_with_mult:
        if rng.random() < base_surv * mult:
            survivors.append(ind)

    if len(survivors) > K_t_int:
        survivors = rng.sample(survivors, K_t_int)
    state.adults = survivors

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
    area_ha: float,
    empirical: EmpiricalParams | None,
    rng: random.Random,
) -> dict:
    state = initial_state(eo, base_individuals, census, area_ha, empirical, rng)
    ever_present: set[str] = set()
    for g in state.adults:
        ever_present.update(g)
    for comp in (state.nd, state.pd, state.ph):
        ever_present.update(comp.alleles_present())

    snapshots: dict[int, dict] = {}
    collapse_year: int | None = None
    min_N_to_date = len(state.adults)
    for year in range(1, HORIZON + 1):
        step_year(state, year, rng)
        N_now = len(state.adults)
        if N_now < min_N_to_date:
            min_N_to_date = N_now
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
                "N_adult": N_now,
                "k_alleles": len(present_now),
                "P_compat": p_compat,
                "collapsed": collapse_year is not None and collapse_year <= year,
                "L": state.L,
                "min_N_to_date": min_N_to_date,
            }
    for cp in CHECKPOINTS:
        snapshots.setdefault(cp, {
            "N_adult": 0, "k_alleles": 0, "P_compat": 0.0,
            "collapsed": True, "L": 0.0, "min_N_to_date": 0,
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


def run_eo(args) -> tuple[list[dict], list[dict]]:
    (eo, base_individuals, census, area_ha, empirical, n_reps, seed) = args
    rng = random.Random(seed)
    reps = [run_replicate(eo, base_individuals, census, area_ha, empirical, rng)
            for _ in range(n_reps)]

    summary_rows = []
    for cp in CHECKPOINTS:
        snaps = [r["snapshots"][cp] for r in reps]
        N_vals = [s["N_adult"] for s in snaps]
        k_vals = [s["k_alleles"] for s in snaps]
        pc_vals = [s["P_compat"] for s in snaps]
        L_vals = [s["L"] for s in snaps]
        min_N_vals = [s["min_N_to_date"] for s in snaps]
        n = max(1, len(snaps))
        p_collapsed = sum(1 for s in snaps if s["collapsed"]) / n
        # Demographic fragility: chronic (current N) vs lifetime (min ever)
        p_N_lt50 = sum(1 for v in N_vals if v < 50) / n
        p_N_lt10 = sum(1 for v in N_vals if v < 10) / n
        p_minN_lt10 = sum(1 for v in min_N_vals if v < 10) / n
        p_extinct_ever = sum(1 for v in min_N_vals if v == 0) / n
        row = {
            "EO": f"EO{eo}",
            "checkpoint_year": cp,
            "P_collapsed": f"{p_collapsed:.3f}",
            "P_N_lt50": f"{p_N_lt50:.3f}",
            "P_N_lt10": f"{p_N_lt10:.3f}",
            "P_minN_lt10_to_date": f"{p_minN_lt10:.3f}",
            "P_adults_extinct_ever": f"{p_extinct_ever:.3f}",
        }
        for q in PCTILES:
            row[f"N_adult_p{q}"] = f"{percentile(N_vals, q):.0f}"
            row[f"k_alleles_p{q}"] = f"{percentile(k_vals, q):.1f}"
            row[f"P_compat_p{q}"] = f"{percentile(pc_vals, q):.3f}"
            row[f"L_p{q}"] = f"{percentile(L_vals, q):.3f}"
        row["min_N_to_date_p10"] = f"{percentile(min_N_vals, 10):.0f}"
        row["min_N_to_date_p50"] = f"{percentile(min_N_vals, 50):.0f}"
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
    global N_REPS
    if len(sys.argv) > 1 and sys.argv[1].startswith("--n-reps="):
        N_REPS = int(sys.argv[1].split("=", 1)[1])
        print(f"Override N_REPS = {N_REPS}")

    print("Loading inputs ...")
    meta = load_metadata()
    individuals_by_eo = load_per_eo_individuals(meta)
    census_by_eo = load_eo_census_and_fecundity()
    areas_by_eo = load_eo_areas()
    empirical_by_eo = load_empirical_inputs()

    if empirical_by_eo:
        print("  Per-EO empirical parameters (from sensitive PVA inputs TSV):")
        print(f"    {'EO':<5}{'σ_env':>8}{'cat_rate':>10}{'ext':>6}")
        for eo in FOCAL_EOS:
            e = empirical_by_eo.get(eo)
            if e:
                print(f"    EO{eo:<3}{e.sigma_env:>8.2f}{e.catastrophe_rate:>10.2f}"
                      f"{('Y' if e.historical_extinction else '.'): >6}")
    else:
        print("  No empirical inputs; using fallback σ_env = 0.30, "
              "catastrophe_rate = 0")

    jobs = []
    for eo in FOCAL_EOS:
        inds = individuals_by_eo[eo]
        if not inds:
            continue
        census = census_by_eo.get(eo, {"census_N": len(inds),
                                        "fec_per_plant": 300.0,
                                        "seeds_banked": 0.0})
        area_ha = areas_by_eo.get(eo, 0.0)
        emp = empirical_by_eo.get(eo)
        seed = SEED + abs(hash(eo)) % 1_000_000
        jobs.append((eo, inds, census, area_ha, emp, N_REPS, seed))

    n_workers = min(len(jobs), max(1, mp.cpu_count() - 2))
    print(f"\nSubmitting {len(jobs)} EO jobs to pool ({n_workers} workers)\n")

    all_summary = []
    all_alleles = []
    with mp.Pool(n_workers) as pool:
        for i, (summ, alls) in enumerate(pool.imap_unordered(run_eo, jobs), 1):
            all_summary.extend(summ)
            all_alleles.extend(alls)
            print(f"  [{i}/{len(jobs)}] done — {summ[0]['EO']}")
            sys.stdout.flush()

    all_summary.sort(key=lambda r: (r["EO"], int(r["checkpoint_year"])))
    all_alleles.sort(key=lambda r: (r["EO"], r["allele"]))

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
    print(f"\nPrediction summary (medians):\n{line}\n{'-' * len(line)}")
    for r in all_summary:
        print("  ".join(str(r[k]).ljust(col_w[k]) for k in cols))


if __name__ == "__main__":
    main()
