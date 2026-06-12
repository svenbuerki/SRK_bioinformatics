#!/usr/bin/env python3
"""
SRK_TP1_pva_v17_maternal_dormancy.py — PVA v1.7 with maternal-dormancy
inheritance.

Extends v1.6 (empirical per-EO parameters) by tagging each adult with a
dormancy class drawn at initialization from the EO's ND / PD / PH split (NV
excluded — it's non-viable, not a maternal phenotype). Seeds inherit their
mother's class deterministically; cohorts are routed into compartments by
maternal class rather than by fixed EO-level shares. Lineages can therefore
go locally extinct: if the PD-class mothers all die in a single bad year,
the PD compartment stops being refilled and drains over ~20 years.

This addresses the v1.6 simplification flagged on 2026-05-27: dormancy was
treated as a fixed environmental EO parameter rather than as a heritable
maternal trait. Most impactful for EO27 (60 % dormant-lineage mothers) and
EO67 (small enough that losing a few PD/PH mothers in one year empties the
bank).

New outputs vs v1.6:
  * P_no_PD_mothers_at_cp, P_no_PH_mothers_at_cp per EO × checkpoint
  * adult_class_share_p50 percentiles per class

Inputs
------
EO_data_sensitive_data/extracted_pva_inputs.tsv  ← required (same as v1.6)

Outputs
-------
Tables/SRK_TP1_pva_v17_summary.tsv             (per EO × checkpoint)
Tables/SRK_TP1_pva_v17_alleles.tsv             (per EO × allele)

Sensitive-data handling: same as v1.6 — reads sensitive census-derived
parameters, writes only aggregated PVA metrics.
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

OUT_SUMMARY = TABLES / "SRK_TP1_pva_v17_summary.tsv"
OUT_ALLELES = TABLES / "SRK_TP1_pva_v17_alleles.tsv"

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

# v1.7 refinement (2026-05-28): assign initial dormancy classes by the
# largest-remainder method instead of multinomial draws. This removes
# replicate-to-replicate stochasticity in starting class counts (which
# matters most for small EOs like EO67, N=20). When True, every replicate
# starts with exactly the same class proportions; the only difference
# across replicates is which individuals happen to be in each class.
DETERMINISTIC_INIT_CLASSES = True

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
    """Per-mating SI rejection (kept for initial-bank build only)."""
    for _ in range(SI_MAX_ATTEMPTS):
        p1, p2 = rng.sample(parents, 2)
        g2 = rng.sample(p2, 2)
        stigma_set = set(p1)
        if g2[0] not in stigma_set and g2[1] not in stigma_set:
            g1 = rng.sample(p1, 2)
            return tuple(g1 + g2)
    return None


def sample_outcross_from_mother(
    mother: tuple[str, ...],
    all_adults: list[tuple[str, ...]],
    rng: random.Random,
) -> tuple[str, ...] | None:
    """SI-respecting outcross where the mother is fixed; pollen is drawn
    from any other adult. Mother is the stigma carrier — offspring inherits
    her dormancy class."""
    n = len(all_adults)
    if n < 2:
        return None
    stigma = set(mother)
    for _ in range(SI_MAX_ATTEMPTS):
        i = rng.randrange(n)
        father = all_adults[i]
        if father is mother:
            continue
        g_father = rng.sample(father, 2)
        if g_father[0] not in stigma and g_father[1] not in stigma:
            g_mother = rng.sample(mother, 2)
            return tuple(g_mother + g_father)
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


DORMANCY_CLASSES = ("ND", "PD", "PH")


@dataclass
class SeedCohort:
    age: int
    n: float
    representatives: list[tuple[str, ...]]
    maternal_class: str        # ND / PD / PH — every seed in this cohort
                               # inherits the same class from its mother
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

    def germinate(self, rng: random.Random) -> list[tuple[tuple[str, ...], float, str]]:
        rates = COMPARTMENT_RATES[self.name]
        gr = rates["germ"]
        recruits: list[tuple[tuple[str, ...], float, str]] = []
        for coh in self.cohorts:
            n_germ = int(round(coh.n * gr))
            if n_germ <= 0 or not coh.representatives:
                continue
            picks = rng.choices(coh.representatives, k=n_germ)
            id_mult = 1.0 - coh.f_sc        # f_sc field holds the ID-derived survival reduction
            mclass = coh.maternal_class
            recruits.extend((p, id_mult, mclass) for p in picks)
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
    adult_dormancy: list[str]   # parallel to adults, index-aligned;
                                # values ∈ DORMANCY_CLASSES
    K0: float
    fec_per_plant: float
    f_ND: float                 # EO baseline (used only at init)
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

    # K: demographic-only (census × K_MULT_MEAN).
    #
    # 2026-05-28: dropped the area-based cap (K_area = area_ha × 5000).
    # Observed plant-per-ha densities in the dataset span 79-29552 (370x),
    # so a single AREA_MAX_DENSITY_PER_HA constant could not represent both
    # large diffuse populations (EO27 @ 79/ha) and dense compact ones
    # (EO70 @ 29552/ha). The census number already integrates whatever
    # density the local habitat supports, so a separate area cap was both
    # uncalibrated and redundant.  Area data is still used in the
    # spatial-fragmentation analyses (LEPA_EO_spatial_clustering), just
    # not as a per-year K constraint here.  See the discussion on
    # 2026-05-28 for the full reasoning.
    K_val = census_N * K_MULT_MEAN
    _ = area_ha  # kept in the call signature for compatibility; unused here

    if census_N <= len(individuals):
        adults_init = rng.sample(individuals, census_N)
    else:
        adults_init = [rng.choice(individuals) for _ in range(census_N)]

    # Assign each starter adult a dormancy class from the EO's ND/PD/PH split
    # (NV is non-viable, not a maternal phenotype, so it's excluded here).
    class_weights = (f_ND_n, f_PD_n, f_PH_n)
    N_init = len(adults_init)
    if sum(class_weights) > 0 and N_init > 0:
        if DETERMINISTIC_INIT_CLASSES:
            # Largest-remainder method: exact class proportions every replicate.
            raw = [w * N_init for w in class_weights]
            base_counts = [int(x) for x in raw]
            slack = N_init - sum(base_counts)
            residuals = sorted(
                ((raw[i] - base_counts[i], i) for i in range(3)),
                reverse=True,
            )
            for k in range(slack):
                base_counts[residuals[k][1]] += 1
            dormancy_init = []
            for cls, count in zip(DORMANCY_CLASSES, base_counts):
                dormancy_init.extend([cls] * count)
            rng.shuffle(dormancy_init)
        else:
            dormancy_init = rng.choices(
                DORMANCY_CLASSES, weights=class_weights, k=N_init
            )
    else:
        dormancy_init = ["ND"] * N_init

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
        eo=eo, adults=adults_init, adult_dormancy=dormancy_init,
        K0=K_val, fec_per_plant=fec,
        f_ND=f_ND_n, f_PD=f_PD_n, f_PH=f_PH_n, f_NV=f_NV, L=L0, f_SC=f_SC,
        sigma_env=sigma_env, catastrophe_rate=cat_rate,
    )

    # Starting seed bank: split by EO baseline shares (initial bank is a
    # historical snapshot, not generated by current adults). Each starting
    # cohort is tagged with its own compartment's class so it propagates
    # correctly if any seed survives to germinate.
    init_bank = max(0.0, float(census.get("seeds_banked", 0.0)))
    if init_bank > 0 and adults_init:
        rep_seeds = _initial_bank_representatives(adults_init, REPRESENTATIVE_N, rng)
        for compartment, share, mclass in [
            (state.nd, f_ND, "ND"),
            (state.pd, f_PD, "PD"),
            (state.ph, f_PH, "PH"),
        ]:
            n_in_pool = init_bank * share
            if n_in_pool < 0.5:
                continue
            compartment.cohorts.append(SeedCohort(
                age=1, n=n_in_pool, representatives=list(rep_seeds),
                maternal_class=mclass, f_sc=0.0,
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

        # Per-mother per-year seed counts (apply uniformly across classes).
        per_mother_total = state.fec_per_plant
        outcr_share = (1 - sc_share) * p_compat * (1 - state.L)
        si_self_share = (1 - sc_share) * state.L
        sc_repr_share = sc_share
        route_total = outcr_share + si_self_share + sc_repr_share
        if route_total <= 0:
            route_total = 1e-9

        # Mothers grouped by their inherited dormancy class.
        mothers_by_class: dict[str, list[tuple[str, ...]]] = {
            c: [] for c in DORMANCY_CLASSES
        }
        for ind, dclass in zip(state.adults, state.adult_dormancy):
            mothers_by_class[dclass].append(ind)

        for mclass in DORMANCY_CLASSES:
            mothers_C = mothers_by_class[mclass]
            n_C = len(mothers_C)
            if n_C == 0:
                continue
            f_C_t = n_C / N_adult
            seeds_C_total = n_C * per_mother_total * route_total
            viable_C = seeds_C_total * (1 - state.f_NV)
            if viable_C < 0.5:
                continue

            n_reps_C = max(1, int(round(REPRESENTATIVE_N * f_C_t)))
            n_outcr = int(round(n_reps_C * outcr_share / route_total))
            n_si_self = int(round(n_reps_C * si_self_share / route_total))
            n_sc = n_reps_C - n_outcr - n_si_self
            if n_sc < 0:
                # rounding can over-allocate; clip
                n_si_self += n_sc
                n_sc = 0

            reps: list[tuple[str, ...]] = []
            if n_outcr > 0 and N_adult >= 2:
                succeeded = 0
                while succeeded < n_outcr:
                    mother = rng.choice(mothers_C)
                    child = sample_outcross_from_mother(mother, state.adults, rng)
                    if child is None:
                        # bail out on persistent SI rejection — remaining slots
                        # become selfed (same maternal class)
                        break
                    reps.append(child)
                    succeeded += 1
                n_si_self += (n_outcr - succeeded)
            elif N_adult == 1:
                n_si_self += n_outcr

            for _ in range(n_si_self):
                mother = rng.choice(mothers_C)
                reps.append(sample_selfed_offspring(mother, rng))
            for _ in range(n_sc):
                mother = rng.choice(mothers_C)
                reps.append(sample_selfed_offspring(mother, rng))

            cohort_f_sc = (n_sc / max(1, n_reps_C))

            compartment = {"ND": state.nd, "PD": state.pd, "PH": state.ph}[mclass]
            compartment.cohorts.append(SeedCohort(
                age=0, n=viable_C, representatives=list(reps),
                maternal_class=mclass, f_sc=cohort_f_sc * ID_SEVERITY,
            ))

    for comp in (state.nd, state.pd, state.ph):
        comp.age_one_year()

    juveniles_with_mult: list[tuple[tuple[str, ...], float, str]] = []
    for comp in (state.nd, state.pd, state.ph):
        juveniles_with_mult.extend(comp.germinate(rng))

    survivors: list[tuple[str, ...]] = []
    survivors_dorm: list[str] = []
    base_surv = JUV_SURVIVAL
    if is_catastrophe:
        base_surv *= CATASTROPHE_RECRUIT_MULT
    for ind, mult, dclass in juveniles_with_mult:
        if rng.random() < base_surv * mult:
            survivors.append(ind)
            survivors_dorm.append(dclass)

    if len(survivors) > K_t_int:
        keep_idx = rng.sample(range(len(survivors)), K_t_int)
        survivors = [survivors[i] for i in keep_idx]
        survivors_dorm = [survivors_dorm[i] for i in keep_idx]
    state.adults = survivors
    state.adult_dormancy = survivors_dorm

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
            class_counts = {c: 0 for c in DORMANCY_CLASSES}
            for c in state.adult_dormancy:
                class_counts[c] = class_counts.get(c, 0) + 1
            snapshots[year] = {
                "N_adult": N_now,
                "k_alleles": len(present_now),
                "P_compat": p_compat,
                "collapsed": collapse_year is not None and collapse_year <= year,
                "L": state.L,
                "min_N_to_date": min_N_to_date,
                "n_ND_mothers": class_counts["ND"],
                "n_PD_mothers": class_counts["PD"],
                "n_PH_mothers": class_counts["PH"],
            }
    for cp in CHECKPOINTS:
        snapshots.setdefault(cp, {
            "N_adult": 0, "k_alleles": 0, "P_compat": 0.0,
            "collapsed": True, "L": 0.0, "min_N_to_date": 0,
            "n_ND_mothers": 0, "n_PD_mothers": 0, "n_PH_mothers": 0,
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
        # Lineage extinction (v1.7): probability that a given maternal class
        # has zero mothers at the checkpoint. Counts replicates where the
        # class was present at init but absent at the checkpoint.
        n_ND_vals = [s["n_ND_mothers"] for s in snaps]
        n_PD_vals = [s["n_PD_mothers"] for s in snaps]
        n_PH_vals = [s["n_PH_mothers"] for s in snaps]
        p_no_ND = sum(1 for v in n_ND_vals if v == 0) / n
        p_no_PD = sum(1 for v in n_PD_vals if v == 0) / n
        p_no_PH = sum(1 for v in n_PH_vals if v == 0) / n
        row = {
            "EO": f"EO{eo}",
            "checkpoint_year": cp,
            "P_collapsed": f"{p_collapsed:.3f}",
            "P_N_lt50": f"{p_N_lt50:.3f}",
            "P_N_lt10": f"{p_N_lt10:.3f}",
            "P_minN_lt10_to_date": f"{p_minN_lt10:.3f}",
            "P_adults_extinct_ever": f"{p_extinct_ever:.3f}",
            "P_no_ND_mothers": f"{p_no_ND:.3f}",
            "P_no_PD_mothers": f"{p_no_PD:.3f}",
            "P_no_PH_mothers": f"{p_no_PH:.3f}",
        }
        for q in PCTILES:
            row[f"N_adult_p{q}"] = f"{percentile(N_vals, q):.0f}"
            row[f"k_alleles_p{q}"] = f"{percentile(k_vals, q):.1f}"
            row[f"P_compat_p{q}"] = f"{percentile(pc_vals, q):.3f}"
            row[f"L_p{q}"] = f"{percentile(L_vals, q):.3f}"
        row["min_N_to_date_p10"] = f"{percentile(min_N_vals, 10):.0f}"
        row["min_N_to_date_p50"] = f"{percentile(min_N_vals, 50):.0f}"
        row["n_ND_mothers_p50"] = f"{percentile(n_ND_vals, 50):.0f}"
        row["n_PD_mothers_p50"] = f"{percentile(n_PD_vals, 50):.0f}"
        row["n_PH_mothers_p50"] = f"{percentile(n_PH_vals, 50):.0f}"
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
