#!/usr/bin/env python3
"""
SRK_TP1_pva_sensitivity.py — Population Viability Analysis v1.5
sensitivity sweep with parallel-across-EO execution.

Spec: SRK_TP1_pva_design.md §6b (SC + ID) and §8b (7 scenarios).

Extends the v1 baseline simulator (SRK_TP1_pva.py) with six alternative
scenarios and a "worst case" combination. Adds:
  * area-bounded K (area_capped)
  * catastrophic year events (catastrophe)
  * K-decline trend (K_decline_1pct)
  * per-mating SI rejection — proper NFDS (si_proper)
  * SC + inbreeding depression (sc_inbreeding) — EO76 only
  * combined worst-case envelope

Parallelism: multiprocessing.Pool, one worker per (EO × scenario) job —
each replicate is sequential within the worker. All 7 scenarios × 6 EOs
= 42 independent jobs running across up to 6 cores at once.

Outputs
-------
Tables/SRK_TP1_pva_sensitivity_summary.tsv  -- per (EO, scenario, checkpoint)
Tables/SRK_TP1_pva_sensitivity_alleles.tsv  -- per (EO, scenario, allele)
"""

from __future__ import annotations

import csv
import math
import multiprocessing as mp
import os
import random
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from pathlib import Path

REPO = Path(__file__).resolve().parent
TABLES = REPO / "Tables"

METRICS_TSV = TABLES / "SRK_EO_allele_richness.tsv"
GENOTYPES_TSV = TABLES / "SRK_individual_allele_genotypes.tsv"
METADATA_CSV = TABLES / "sampling_metadata.csv"
AREAS_CSV = TABLES / "EO_group_BL_summary.csv"
FECUNDITY_CSV = Path(
    "/Users/sven/Documents/Current_projects/NSF24-543_Self-incompatibility"
    "/Brainstorm_NSF/Data/Data_aggregation_by_location_fecundity.csv"
)

OUT_SUMMARY = TABLES / "SRK_TP1_pva_sensitivity_summary.tsv"
OUT_ALLELES = TABLES / "SRK_TP1_pva_sensitivity_alleles.tsv"

# ----- Core constants (kept in sync with SRK_TP1_pva.py) -------------------

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
SIGMA_ENV = 0.30
L_CAP = 0.50

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

# SC individuals per EO (5 in EO76 per Step 12c SI-escape candidates)
SC_INDIVIDUALS = {eo: 0 for eo in FOCAL_EOS}
SC_INDIVIDUALS["76"] = 5

# Per-mating SI rejection: max attempts per requested offspring before giving up
SI_MAX_ATTEMPTS = 100

# ----- Scenario configuration ----------------------------------------------


@dataclass(frozen=True)
class Scenario:
    name: str
    area_capped: bool = False
    max_density_per_ha: float = 5000.0   # used iff area_capped=True
    catastrophe_rate: float = 0.0        # P(year is catastrophe)
    catastrophe_recruit_mult: float = 0.05
    K_decline_rate: float = 0.0          # annual K decay (fraction)
    si_proper: bool = False              # per-mating SI rejection
    sc_inbreeding: bool = False          # SC + ID active
    id_severity: float = 0.5             # SC offspring survival multiplier reduction


SCENARIOS: dict[str, Scenario] = {
    "baseline": Scenario("baseline"),
    "area_capped": Scenario("area_capped", area_capped=True),
    "catastrophe": Scenario("catastrophe",
                              catastrophe_rate=0.05,
                              catastrophe_recruit_mult=0.05),
    "K_decline_1pct": Scenario("K_decline_1pct", K_decline_rate=0.01),
    "si_proper": Scenario("si_proper", si_proper=True),
    "sc_inbreeding": Scenario("sc_inbreeding", sc_inbreeding=True),
    "worst_case": Scenario("worst_case",
                             area_capped=True,
                             catastrophe_rate=0.05,
                             catastrophe_recruit_mult=0.05,
                             K_decline_rate=0.01,
                             sc_inbreeding=True),
}


# ============================================================================
# Data loading (same as v1)
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
    """Aggregate Area_ha per focal EO from EO_group_BL_summary.csv.
    Handles shared-polygon rows (e.g., 'EO118; EO76') by splitting area
    equally across listed EOs."""
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
    p1, p2 = rng.sample(parents, 2)
    g1 = rng.sample(p1, 2)
    g2 = rng.sample(p2, 2)
    return tuple(g1 + g2)


def sample_outcross_offspring_si(
    parents: list[tuple[str, ...]], rng: random.Random
) -> tuple[str, ...] | None:
    """Per-mating SI rejection: redraw parent pair until pollen gamete
    shares no alleles with stigma. Returns None if SI_MAX_ATTEMPTS
    consecutive attempts all fail (population is essentially incompatible)."""
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
    f_sc: float = 0.0  # fraction of cohort from SC parents (for ID survival adjust)

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
        """Returns list of (individual_genotype, survival_multiplier) where
        the multiplier captures cohort-average SC-derived inbreeding
        depression. The standard JUV_SURVIVAL is applied downstream."""
        rates = COMPARTMENT_RATES[self.name]
        gr = rates["germ"]
        recruits: list[tuple[tuple[str, ...], float]] = []
        for coh in self.cohorts:
            n_germ = int(round(coh.n * gr))
            if n_germ <= 0 or not coh.representatives:
                continue
            # Per-juvenile survival multiplier (accounts for cohort-average ID)
            picks = rng.choices(coh.representatives, k=n_germ)
            recruits.extend((p, 1.0 - coh.f_sc) for p in picks)
            # f_sc fraction of cohort suffers full id_severity reduction;
            # callers multiply JUV_SURVIVAL × multiplier with the cohort's f_sc
            # baked into the multiplier upstream.
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
    K0: float                    # initial K (used for K_decline)
    fec_per_plant: float
    f_ND: float
    f_PD: float
    f_PH: float
    f_NV: float
    L: float
    f_SC: float                  # fraction of adults that are SC (for sc_inbreeding)
    nd: CompartmentState = field(default_factory=lambda: CompartmentState("ND"))
    pd: CompartmentState = field(default_factory=lambda: CompartmentState("PD"))
    ph: CompartmentState = field(default_factory=lambda: CompartmentState("PH"))
    extinction_year: dict[str, int] = field(default_factory=dict)


def initial_state(
    eo: str,
    individuals: list[tuple[str, ...]],
    census: dict,
    area_ha: float,
    scenario: Scenario,
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

    # K calculation depends on scenario.
    # Baseline: K ~ LogNormal(census × 1.5, σ=0.30), bounded below by census.
    # area_capped: K = min(area × max_density, census × 1.5) — area is an
    # upper bound on the demographic ceiling, not a replacement. Falls back
    # to baseline when area_ha is zero (data gap).
    if scenario.area_capped and area_ha > 0.001:
        K_area = area_ha * scenario.max_density_per_ha
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

    f_SC = 0.0
    if scenario.sc_inbreeding:
        f_SC = SC_INDIVIDUALS.get(eo, 0) / max(census_N, 1)

    state = EOState(
        eo=eo, adults=adults_init, K0=K_val, fec_per_plant=fec,
        f_ND=f_ND_n, f_PD=f_PD_n, f_PH=f_PH_n, f_NV=f_NV, L=L0, f_SC=f_SC,
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


def step_year(
    state: EOState, year: int, scenario: Scenario, rng: random.Random,
) -> None:
    # Compute K_t for this year (K_decline + env noise)
    K_after_trend = state.K0 * ((1.0 - scenario.K_decline_rate) ** year)
    K_t = K_after_trend * math.exp(rng.gauss(0, SIGMA_ENV))
    K_t_int = max(0, int(round(K_t)))

    # Catastrophe roll (affects this year's recruitment)
    is_catastrophe = (
        scenario.catastrophe_rate > 0
        and rng.random() < scenario.catastrophe_rate
    )

    # ---- Reproduction ----
    if state.adults:
        freqs = gamete_freqs_from_adults(state.adults)
        p_compat = p_compat_tetraploid(freqs)
        N_adult = len(state.adults)
        sc_share = state.f_SC

        # Seeds from SC parents (always self) — separate cohort tag
        seeds_sc = N_adult * sc_share * state.fec_per_plant
        # Seeds from SI parents
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

            # Outcross representatives
            if n_outcr > 0 and N_adult >= 2:
                if scenario.si_proper:
                    succeeded = 0
                    while succeeded < n_outcr:
                        child = sample_outcross_offspring_si(state.adults, rng)
                        if child is None:
                            break  # population effectively incompatible
                        reps.append(child)
                        succeeded += 1
                    # If we couldn't fill outcross slots, redirect to selfing
                    n_si_self += (n_outcr - succeeded)
                else:
                    for _ in range(n_outcr):
                        reps.append(sample_outcross_offspring(state.adults, rng))
            elif N_adult == 1:
                # No outcross possible; everything becomes selfing
                n_si_self += n_outcr
                n_outcr = 0

            # SI selfed (leaky-self) representatives
            for _ in range(n_si_self):
                parent = rng.choice(state.adults)
                reps.append(sample_selfed_offspring(parent, rng))
            # SC representatives (always selfed)
            for _ in range(n_sc):
                parent = rng.choice(state.adults)
                reps.append(sample_selfed_offspring(parent, rng))

            # f_sc for the cohort: fraction of representatives from SC parents.
            # SC reps survive at reduced rate, captured via f_sc field.
            cohort_f_sc = (n_sc / REPRESENTATIVE_N) if REPRESENTATIVE_N else 0.0
            # ID is applied later as survival multiplier scaled by f_sc
            cohort_id_mult = 1.0 - cohort_f_sc * scenario.id_severity

            for comp, share in [
                (state.nd, state.f_ND),
                (state.pd, state.f_PD),
                (state.ph, state.f_PH),
            ]:
                n_in_pool = viable_seeds * share
                if n_in_pool >= 0.5 and share > 0:
                    comp.cohorts.append(SeedCohort(
                        age=0, n=n_in_pool, representatives=list(reps),
                        # store the cohort's average ID multiplier in f_sc field
                        # (re-purposed: multiplier on survival applied at germination)
                        f_sc=cohort_f_sc * scenario.id_severity,
                    ))

    # ---- Age cohorts ----
    for comp in (state.nd, state.pd, state.ph):
        comp.age_one_year()

    # ---- Germinate ----
    juveniles_with_mult: list[tuple[tuple[str, ...], float]] = []
    for comp in (state.nd, state.pd, state.ph):
        juveniles_with_mult.extend(comp.germinate(rng))

    # ---- Juvenile survival (with per-individual mult) ----
    survivors: list[tuple[str, ...]] = []
    base_surv = JUV_SURVIVAL
    if is_catastrophe:
        base_surv *= scenario.catastrophe_recruit_mult
    for ind, mult in juveniles_with_mult:
        if rng.random() < base_surv * mult:
            survivors.append(ind)

    # ---- Cap to K ----
    if len(survivors) > K_t_int:
        survivors = rng.sample(survivors, K_t_int)
    state.adults = survivors

    # ---- Update L ----
    if state.adults:
        n_aaaa = sum(1 for g in state.adults if class_of(g) == "AAAA")
        prop_aaaa = n_aaaa / len(state.adults)
        state.L = min(L_CAP, prop_aaaa / 3.5)
    # f_SC unchanged through generations in v1.5 (proportion is structural,
    # not mendelian-inherited; a future iteration could let it drift).


# ============================================================================
# Replicate runner
# ============================================================================


def run_replicate(
    eo: str,
    base_individuals: list[tuple[str, ...]],
    census: dict,
    area_ha: float,
    scenario: Scenario,
    rng: random.Random,
) -> dict:
    state = initial_state(eo, base_individuals, census, area_ha, scenario, rng)
    ever_present: set[str] = set()
    for g in state.adults:
        ever_present.update(g)
    for comp in (state.nd, state.pd, state.ph):
        ever_present.update(comp.alleles_present())

    snapshots: dict[int, dict] = {}
    collapse_year: int | None = None
    for year in range(1, HORIZON + 1):
        step_year(state, year, scenario, rng)
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


def run_eo_scenario(args) -> tuple[list[dict], list[dict]]:
    """Top-level worker for multiprocessing.Pool."""
    (eo, scenario_name, base_individuals, census, area_ha, n_reps, seed) = args
    scenario = SCENARIOS[scenario_name]
    rng = random.Random(seed)
    reps = [run_replicate(eo, base_individuals, census, area_ha, scenario, rng)
            for _ in range(n_reps)]

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
            "scenario": scenario_name,
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
            "scenario": scenario_name,
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
# Sensitivity orchestrator
# ============================================================================


def main() -> None:
    global N_REPS
    # Allow override from CLI for smoke testing
    if len(sys.argv) > 1 and sys.argv[1].startswith("--n-reps="):
        N_REPS = int(sys.argv[1].split("=", 1)[1])
        print(f"Override N_REPS = {N_REPS}")
    print("Loading inputs ...")
    meta = load_metadata()
    individuals_by_eo = load_per_eo_individuals(meta)
    census_by_eo = load_eo_census_and_fecundity()
    areas_by_eo = load_eo_areas()
    print(f"  areas: " +
          ", ".join(f"EO{eo}={areas_by_eo[eo]:.3f}ha" for eo in FOCAL_EOS))
    print(f"  SC individuals: " +
          ", ".join(f"EO{eo}={SC_INDIVIDUALS[eo]}" for eo in FOCAL_EOS))
    print(f"  scenarios: {list(SCENARIOS.keys())}\n")

    # Build job list
    jobs = []
    for eo in FOCAL_EOS:
        inds = individuals_by_eo[eo]
        if not inds:
            continue
        census = census_by_eo.get(eo, {"census_N": len(inds),
                                        "fec_per_plant": 300.0,
                                        "seeds_banked": 0.0})
        area_ha = areas_by_eo.get(eo, 0.0)
        for sname in SCENARIOS:
            # Per-job seed: master + hash of (eo, scenario)
            seed = SEED + abs(hash((eo, sname))) % 1_000_000
            jobs.append((eo, sname, inds, census, area_ha, N_REPS, seed))
    print(f"Submitting {len(jobs)} (EO × scenario) jobs to pool ...")

    n_workers = min(6, max(1, mp.cpu_count() - 2))
    print(f"  workers: {n_workers}\n")

    all_summary = []
    all_alleles = []
    with mp.Pool(n_workers) as pool:
        for i, (summ, alls) in enumerate(pool.imap_unordered(run_eo_scenario, jobs), 1):
            all_summary.extend(summ)
            all_alleles.extend(alls)
            print(f"  [{i}/{len(jobs)}] done — {summ[0]['EO']}/{summ[0]['scenario']}")
            sys.stdout.flush()

    # Add delta_P_collapsed_vs_baseline
    base_lookup: dict[tuple[str, int], float] = {}
    for r in all_summary:
        if r["scenario"] == "baseline":
            base_lookup[(r["EO"], int(r["checkpoint_year"]))] = float(r["P_collapsed"])
    for r in all_summary:
        base = base_lookup.get((r["EO"], int(r["checkpoint_year"])), 0.0)
        r["delta_P_collapsed_vs_baseline"] = f"{float(r['P_collapsed']) - base:+.3f}"

    # Sort for readability: EO, scenario (baseline first), checkpoint
    scenario_order = list(SCENARIOS.keys())
    all_summary.sort(key=lambda r: (r["EO"], scenario_order.index(r["scenario"]),
                                      int(r["checkpoint_year"])))
    all_alleles.sort(key=lambda r: (r["EO"], scenario_order.index(r["scenario"]),
                                      r["allele"]))

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

    # Compact at-a-glance summary table at y50
    cols = ["EO", "scenario", "P_collapsed", "N_adult_p50",
            "k_alleles_p50", "P_compat_p50",
            "delta_P_collapsed_vs_baseline"]
    y50_rows = [r for r in all_summary if int(r["checkpoint_year"]) == 50]
    if y50_rows:
        col_w = {k: max(len(k), max(len(str(r[k])) for r in y50_rows))
                 for k in cols}
        line = "  ".join(k.ljust(col_w[k]) for k in cols)
        print(f"\nSensitivity at y50 (medians):\n{line}\n{'-' * len(line)}")
        for r in y50_rows:
            print("  ".join(str(r[k]).ljust(col_w[k]) for k in cols))


if __name__ == "__main__":
    main()
