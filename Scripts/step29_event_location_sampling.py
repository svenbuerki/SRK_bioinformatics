#!/usr/bin/env python3
"""Step 29 — Event- and location-level SRK sampling design.

Extends Step 28 one level up: given the per-mother recipe already in place,
Step 29 answers *how many mothers to sample per event*, and *how many events
to sample per location*, so that the location-level maternal SRK allele pool
is characterised at a stated coverage.

Applies the same random-mating coupon-collector logic used in Step 28, one
level higher: instead of drawing pollen alleles from an event's pollen pool,
we draw maternal alleles by sampling mothers across the location's events.

Inputs (read-only)
    /Users/sven/Documents/Current_projects/LEPA_fieldwork_protocol/SQL_DB/LEPA_SQL.db

Outputs (SAMPLING family — no prediction machinery here)
    Tables/Phase5/step29_sampling_per_event.tsv
        One row per event with recommended mothers to sample.
    Tables/Phase5/step29_sampling_per_location.tsv
        One row per location with total sampling design.
    Tables/Phase5/step29_location_coverage_curves.tsv
        Analytical + simulated coverage curves binned by total location size.
    figures/Phase5/step29_location_coverage_curves.png/pdf
        Expected maternal-allele coverage vs mothers sampled per location.
"""
from __future__ import annotations

import math
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
DEFAULT_TABLES = Path("Tables/Phase5")
DEFAULT_FIGURES = Path("figures/Phase5")
# Empirical prior over Fgs (same file used by Step 30).
DEFAULT_PRIOR_TSV = DEFAULT_TABLES / "step26i_L1_carrier_inventory.tsv"
# Per-mother sampling design from Step 28 — carries the ACTUAL per-mother
# seed budgets (n_achievable_exp / n_achievable_miss are already capped at
# each mother's germplasmQuantityEstimate).
DEFAULT_STEP28_TSV = DEFAULT_TABLES / "step28_seed_sampling_per_mother.tsv"

# LEPA is tetraploid (2n = 4x): every plant carries PLOIDY = 4 SRK allele
# copies, and each seed carries PATERNAL_ALLELES_PER_SEED = 2 paternal
# alleles. Sourced from Step 28 so both scripts stay in lock-step.
from step28_seed_sampling_per_mother import (
    PLOIDY,
    PATERNAL_ALLELES_PER_SEED,
    n_for_miss_probability,
)

EXPECTED_COVERAGE = 0.90
DEFAULT_SEEDS_PER_MOTHER = n_for_miss_probability()   # tetraploid Rule 2 = 15

LOCATION_SIZE_BINS = [
    ("<20",     0,   20),
    ("20-50",  21,   50),
    ("51-100", 51,  100),
    ("101-250", 101, 250),
    (">250",   251, 10_000_000),
]
N_SIM = 1_000
MOTHERS_GRID = np.arange(1, 121)


# ---------------------------------------------------------------------------
# Coverage math (uniform maternal-allele draws over K = PLOIDY * N_fertile)
# Each mother contributes PLOIDY maternal alleles under tetraploidy.
# ---------------------------------------------------------------------------
def n_mothers_for_coverage(K: int, target: float = EXPECTED_COVERAGE) -> int:
    """Smallest M s.t. E[coverage] >= target, drawing PLOIDY·M maternal
    alleles uniformly from K. Under uniform p_j = 1/K,
        E[coverage] = 1 - (1 - 1/K)^(PLOIDY·M)
    Solve for M: M = ceil( log(1-target) / (PLOIDY · log(1 - 1/K)) ).
    """
    if K <= 1:
        return 1
    return int(math.ceil(math.log(1 - target)
                         / (PLOIDY * math.log(1 - 1.0 / K))))


def expected_coverage(K: int, M: int) -> float:
    if K <= 1:
        return 1.0
    return 1.0 - (1.0 - 1.0 / K) ** (PLOIDY * M)


def simulate_coverage(K: int, M: int, n_sim: int = N_SIM,
                       rng: np.random.Generator | None = None) -> np.ndarray:
    if K <= 1:
        return np.ones(n_sim)
    rng = rng or np.random.default_rng(2029)
    draws = rng.integers(0, K, size=(n_sim, PLOIDY * M))
    counts = np.array([np.unique(r).size for r in draws])
    return counts / K


# ---------------------------------------------------------------------------
# Empirical prior helpers (P1 species-wide Fg frequencies)
# ---------------------------------------------------------------------------
def load_p1_prior(carrier_tsv: Path) -> np.ndarray:
    """Return the P1 species-wide Fg frequency vector as a 1-D numpy array
    (rows sum to 1)."""
    df = pd.read_csv(carrier_tsv, sep="\t", encoding="utf-8-sig")
    fg = df.groupby("Fg", as_index=False)["n_carriers"].sum()
    f = fg["n_carriers"].values.astype(float)
    return f / f.sum()


def expected_fg_coverage_empirical(f: np.ndarray, A: int) -> float:
    """E[fraction of distinct Fgs observed] after A allele draws under the
    empirical frequency vector f (of length K_fg). Non-uniform-aware.
    """
    if A <= 0:
        return 0.0
    K_fg = len(f)
    return float(np.sum(1.0 - (1.0 - f) ** A) / K_fg)


def alleles_for_fg_coverage_empirical(f: np.ndarray,
                                      target: float = EXPECTED_COVERAGE,
                                      max_A: int = 200_000) -> int:
    """Smallest A such that E[fraction of distinct Fgs observed] >= target
    under the empirical prior f. No closed form; solve by bracketed search.
    Returns max_A if the target is unreachable within that budget.
    """
    K_fg = len(f)
    # Coverage is monotone increasing in A → doubling search + bisection.
    hi = 1
    while (expected_fg_coverage_empirical(f, hi) < target) and (hi < max_A):
        hi *= 2
    if expected_fg_coverage_empirical(f, hi) < target:
        return max_A
    lo = hi // 2 if hi > 1 else 0
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if expected_fg_coverage_empirical(f, mid) >= target:
            hi = mid
        else:
            lo = mid
    return int(hi)


# ---------------------------------------------------------------------------
# 1 — pull event & location table from the LEPA SQL DB
# ---------------------------------------------------------------------------
def _extract_year(date_series: pd.Series) -> pd.Series:
    """LEPA `eventDate` is stored as MM-DD-YYYY. Extract the 4-digit year."""
    return (date_series.astype(str)
            .str.extract(r"(\d{4})$", expand=False)
            .astype("Int64"))


def load_events_by_location(db_path: Path,
                            year: int | None = None) -> pd.DataFrame:
    con = sqlite3.connect(str(db_path))
    try:
        df = pd.read_sql_query(
            """
            SELECT
                l.locationID   AS locationID,
                l.locationCode AS locationCode,
                e.eventID      AS eventID,
                e.eventDate    AS eventDate,
                e.organismQuantityFertile AS n_fertile
            FROM Locations l
            JOIN Events      e ON e.locationID = l.locationID
            JOIN Occurrences o ON o.eventID    = e.eventID
            JOIN Taxonomy    t ON o.taxonID    = t.taxonID
            WHERE t.genus='Lepidium' AND t.specificEpithet='papilliferum'
              AND e.organismQuantityFertile IS NOT NULL
              AND (o.provenance IS NULL OR o.provenance = 'in situ')
            GROUP BY l.locationID, e.eventID
            """,
            con,
        )
    finally:
        con.close()
    df["event_year"] = _extract_year(df["eventDate"])
    if year is not None:
        n_before = len(df)
        df = df[df["event_year"] == int(year)].copy()
        print(f"[step29] --year {year}: kept {len(df)} of {n_before} events.")
    # `organismQuantityFertile` is sometimes stored as free text like
    # ">200", "36 = 15 failed, 21 fruits...", or "98 - 17 failed, 81 fruited".
    # Take the FIRST integer token (that's the total count in every LEPA
    # convention I've seen so far), coerce to numeric, drop unparseable rows.
    n = (
        df["n_fertile"].astype(str)
        .str.extract(r"(\d+)", expand=False)
        .replace("", np.nan)
    )
    df["n_fertile"] = pd.to_numeric(n, errors="coerce")
    df = df.dropna(subset=["n_fertile"])
    df["n_fertile"] = df["n_fertile"].astype(int)
    return df


# ---------------------------------------------------------------------------
# 2 — per-event and per-location sampling design
# ---------------------------------------------------------------------------
def build_designs(events: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Per-location aggregates
    loc = events.groupby(["locationID", "locationCode"], as_index=False).agg(
        n_events=("eventID", "nunique"),
        total_n_fertile=("n_fertile", "sum"),
        min_n_fertile=("n_fertile", "min"),
        max_n_fertile=("n_fertile", "max"),
    )
    K_loc = (2 * loc["total_n_fertile"]).astype(int)
    loc["K_maternal_pool"] = K_loc
    loc["M_recommended_location"] = [
        n_mothers_for_coverage(int(k)) for k in K_loc
    ]
    loc["expected_coverage_at_M"] = [
        expected_coverage(int(k), int(m))
        for k, m in zip(K_loc, loc["M_recommended_location"])
    ]

    # Per-event: proportional allocation of M_recommended_location, capped at
    # that event's N_fertile.
    per_event = events.merge(
        loc[["locationID", "total_n_fertile", "M_recommended_location"]],
        on="locationID", how="left",
    )
    raw_alloc = (
        per_event["n_fertile"] / per_event["total_n_fertile"]
        * per_event["M_recommended_location"]
    )
    # ceil so we don't under-allocate rounding-wise
    per_event["M_alloc_raw"] = np.ceil(raw_alloc).astype(int)
    per_event["M_alloc_capped"] = np.minimum(
        per_event["M_alloc_raw"], per_event["n_fertile"]
    ).astype(int)

    # Locations: report the achievable M after per-event capping
    ach = per_event.groupby("locationID", as_index=False).agg(
        M_achievable_location=("M_alloc_capped", "sum"),
    )
    loc = loc.merge(ach, on="locationID", how="left")
    loc["expected_coverage_achieved"] = [
        expected_coverage(int(k), int(m))
        for k, m in zip(loc["K_maternal_pool"], loc["M_achievable_location"])
    ]

    per_event = per_event[[
        "locationID", "locationCode", "eventID", "n_fertile",
        "M_alloc_raw", "M_alloc_capped",
    ]].sort_values(["locationID", "eventID"]).reset_index(drop=True)

    loc = loc.sort_values("total_n_fertile").reset_index(drop=True)
    return per_event, loc


def augment_with_empirical_prior(loc: pd.DataFrame,
                                 prior_f: np.ndarray,
                                 seeds_per_mother: int = DEFAULT_SEEDS_PER_MOTHER
                                 ) -> pd.DataFrame:
    """Add empirical-prior (P1) columns to the per-location table.

    New columns
    -----------
    A_target_90pct_P1
        Total allele draws needed at the location to reach 90 % expected
        coverage of the 32 Fgs under the empirical species-wide prior.
    A_delivered_maternal_only
        PLOIDY x M_achievable — allele draws from mother genotypes alone
        (4 per mother under tetraploid).
    A_delivered_with_seeds
        M_achievable x (PLOIDY + PATERNAL_ALLELES_PER_SEED * seeds_per_mother)
        — mother genotypes plus paternal allele draws per genotyped seed.
        Default seeds_per_mother is `DEFAULT_SEEDS_PER_MOTHER` (Step 28
        Rule 2 floor; tetraploid = 15).
    exp_Fg_cov_maternal_only_P1
        E[fraction of distinct Fgs observed] at A_delivered_maternal_only.
    exp_Fg_cov_with_seeds_P1
        E[fraction of distinct Fgs observed] at A_delivered_with_seeds.
    seeds_per_mother_for_90pct_P1
        Smallest integer n such that
        M_achievable x (PLOIDY + PATERNAL_ALLELES_PER_SEED * n)
        >= A_target_90pct_P1. NaN if the target is unreachable given
        M_achievable.
    """
    A_target = alleles_for_fg_coverage_empirical(prior_f, EXPECTED_COVERAGE)
    loc = loc.copy()
    loc["A_target_90pct_P1"] = A_target

    M_ach = loc["M_achievable_location"].astype(int)
    A_mat_only = (PLOIDY * M_ach).astype(int)
    A_with_seeds = (M_ach * (PLOIDY
                             + PATERNAL_ALLELES_PER_SEED * seeds_per_mother)
                    ).astype(int)

    loc["A_delivered_maternal_only"] = A_mat_only
    loc["A_delivered_with_seeds"]    = A_with_seeds
    loc["exp_Fg_cov_maternal_only_P1"] = [
        expected_fg_coverage_empirical(prior_f, int(a)) for a in A_mat_only
    ]
    loc["exp_Fg_cov_with_seeds_P1"] = [
        expected_fg_coverage_empirical(prior_f, int(a)) for a in A_with_seeds
    ]

    # Smallest seeds/mother to hit A_target given M_achievable, under
    # tetraploid: m * (PLOIDY + PATERNAL_ALLELES_PER_SEED * n) >= target
    def seeds_for_target(m: int, target_A: int) -> int | float:
        if m <= 0:
            return np.nan
        need_draws = target_A / m - PLOIDY
        need_seeds = math.ceil(need_draws / PATERNAL_ALLELES_PER_SEED)
        return max(need_seeds, 0)

    loc["seeds_per_mother_for_90pct_P1"] = [
        seeds_for_target(int(m), int(A_target)) for m in M_ach
    ]
    return loc


def _spatial_realised_columns(step28: pd.DataFrame,
                               events_df: pd.DataFrame,
                               prior_f: np.ndarray) -> pd.DataFrame | None:
    """Per-location rollup of the SPATIAL-K version of Rule 1 already
    written to Step 28. For each radius R present in step28's columns
    (`n_achievable_exp_spatial_{R}m`), aggregate to location and compute
    the realised A_delivered and expected Fg coverage under P1.
    """
    import re
    r_cols = [c for c in step28.columns
              if re.match(r"n_achievable_exp_spatial_\d+m$", c)]
    if not r_cols:
        return None
    if "locationID" in step28.columns:
        joined = step28.copy()
    else:
        joined = step28.merge(
            events_df[["eventID", "locationID"]].drop_duplicates(),
            on="eventID", how="left",
        )
    out = None
    for col in r_cols:
        R = re.match(r"n_achievable_exp_spatial_(\d+)m$", col).group(1)
        roll = joined.groupby("locationID", as_index=False).agg(
            **{
                f"total_n_seeds_realised_exp_spatial_{R}m": (col, "sum"),
                "M_actual_in_step28_spatial":               ("occurrenceID", "count"),
            }
        )
        roll[f"A_delivered_realised_exp_spatial_{R}m"] = (
            2 * roll["M_actual_in_step28_spatial"]
            + roll[f"total_n_seeds_realised_exp_spatial_{R}m"]
        )
        roll[f"exp_Fg_cov_realised_exp_spatial_{R}m_P1"] = [
            expected_fg_coverage_empirical(prior_f, int(a))
            for a in roll[f"A_delivered_realised_exp_spatial_{R}m"]
        ]
        roll = roll.drop(columns=["M_actual_in_step28_spatial"])
        out = roll if out is None else out.merge(roll, on="locationID", how="outer")
    return out


def build_field_team_recipe(step28_tsv: Path,
                             events_df: pd.DataFrame) -> pd.DataFrame:
    """Distil Step 28's analyst-facing per-mother table into a compact
    field-team recipe: one row per germplasmID, joined with locationCode,
    with a single `n_seeds_recommended` number and a plain-language
    rationale telling the field team which rule applied.

    Decision rule per mother:
        n_reco = min(seeds_est, min(n_expected_cov_90, n_miss_prob_10pct_95))

    Interpretation:
      - If n_expected_cov_90 <= 29 (Rule 2 floor): Rule 1 is reachable and
        efficient — genotype n_expected_cov_90 seeds.
      - If n_expected_cov_90 > 29: Rule 1 is impractical — use the Rule 2
        floor of 29 seeds.
      - If either exceeds the mother's real seed budget, use everything
        available and flag as budget-limited.
    """
    if not step28_tsv.exists():
        return pd.DataFrame()
    s28 = pd.read_csv(step28_tsv, sep="\t", encoding="utf-8-sig")

    # Attach locationID + locationCode via eventID
    loc_key = events_df[["eventID", "locationID", "locationCode"]].drop_duplicates()
    df = s28.merge(loc_key, on="eventID", how="left")

    n_exp  = df["n_expected_cov_90"].astype(int)
    n_miss = df["n_miss_prob_10pct_95"].astype(int)
    seeds  = df["seeds_est"].astype(int)

    n_target = np.minimum(n_exp, n_miss)          # smaller of the two rules
    n_reco   = np.minimum(seeds, n_target)        # cap at real budget

    rationale = np.where(
        seeds < n_target,
        "budget-limited (all available seeds)",
        np.where(
            n_exp <= n_miss,
            "Rule 1 (90 % expected coverage — small event, cheap to characterise)",
            "Rule 2 (29-seed miss-probability floor — large event, Rule 1 impractical)",
        ),
    )

    out = pd.DataFrame({
        "germplasmID":         df["germplasmID"].astype(int),
        "occurrenceID":        df["occurrenceID"].astype(int),
        "eventID":             df["eventID"].astype(int),
        "locationID":          df["locationID"],
        "locationCode":        df["locationCode"],
        "n_fertile":           df["n_fertile"].astype(int),
        "seeds_available":     seeds,
        "n_seeds_recommended": n_reco.astype(int),
        "rationale":           rationale,
    }).sort_values(["locationCode", "eventID", "germplasmID"]).reset_index(drop=True)

    return out


def augment_with_realised_step28(loc: pd.DataFrame,
                                  events_df: pd.DataFrame,
                                  step28_tsv: Path,
                                  prior_f: np.ndarray) -> pd.DataFrame:
    """Add per-location columns that fold in the ACTUAL per-mother seed
    budgets already computed in Step 28.

    Rationale: seeds per mother is NOT constant — it is the smaller of the
    Rule-1/Rule-2 recommendation and each mother's real seed budget
    (`germplasmQuantityEstimate`). Step 28 already stores this in
    `n_achievable_exp` and `n_achievable_miss`. Step 29's location rollup
    should sum those actual values rather than assume a uniform seeds/mother.

    New columns
    -----------
    M_actual_in_step28
        Number of real LEPA mothers of this location with both a seed budget
        and a known N_fertile (from Step 28).
    total_n_seeds_realised_exp
        Sum over the location's mothers of `n_achievable_exp` — total seeds
        we would genotype under Rule 1 (event-size-dependent target).
    total_n_seeds_realised_miss
        Sum over the location's mothers of `n_achievable_miss` — total seeds
        we would genotype under Rule 2 (29-seed floor, capped by budget).
    A_delivered_realised_exp
        2 x M_actual + total_n_seeds_realised_exp.  Maternal alleles + one
        paternal allele per genotyped seed.
    A_delivered_realised_miss
        2 x M_actual + total_n_seeds_realised_miss.
    exp_Fg_cov_realised_exp_P1
        E[fraction of 32 Fgs observed] at A_delivered_realised_exp under P1.
    exp_Fg_cov_realised_miss_P1
        E[fraction of 32 Fgs observed] at A_delivered_realised_miss under P1.
    """
    if not step28_tsv.exists():
        return loc
    step28 = pd.read_csv(step28_tsv, sep="\t", encoding="utf-8-sig")
    step28 = step28.merge(
        events_df[["eventID", "locationID"]].drop_duplicates(),
        on="eventID", how="left",
    )
    roll = step28.groupby("locationID", as_index=False).agg(
        M_actual_in_step28=("occurrenceID", "count"),
        total_n_seeds_realised_exp=("n_achievable_exp", "sum"),
        total_n_seeds_realised_miss=("n_achievable_miss", "sum"),
    )
    # A_delivered per mother = PLOIDY maternal + PATERNAL_ALLELES_PER_SEED
    # per genotyped seed. Under tetraploid, per mother = 4 + 2·n_seeds.
    roll["A_delivered_realised_exp"] = (
        PLOIDY * roll["M_actual_in_step28"]
        + PATERNAL_ALLELES_PER_SEED * roll["total_n_seeds_realised_exp"]
    )
    roll["A_delivered_realised_miss"] = (
        PLOIDY * roll["M_actual_in_step28"]
        + PATERNAL_ALLELES_PER_SEED * roll["total_n_seeds_realised_miss"]
    )
    roll["exp_Fg_cov_realised_exp_P1"] = [
        expected_fg_coverage_empirical(prior_f, int(a))
        for a in roll["A_delivered_realised_exp"]
    ]
    roll["exp_Fg_cov_realised_miss_P1"] = [
        expected_fg_coverage_empirical(prior_f, int(a))
        for a in roll["A_delivered_realised_miss"]
    ]
    # Recommended seeds/mother using the PERMIT-REALISTIC M (mothers in
    # the seed bank), not the theoretical M_achievable ceiling.
    A_target = alleles_for_fg_coverage_empirical(prior_f, EXPECTED_COVERAGE)
    def _seeds_needed(m: int) -> int | float:
        if m <= 0:
            return np.nan
        need_draws = A_target / m - PLOIDY
        return max(math.ceil(need_draws / PATERNAL_ALLELES_PER_SEED), 0)
    roll["seeds_per_mother_for_90pct_P1_actual"] = [
        _seeds_needed(int(m)) for m in roll["M_actual_in_step28"]
    ]
    # Optional: spatial-K version of Rule 1 (present only if Step 28's
    # spatial extension has been run and its columns are in step28.tsv).
    spatial_roll = _spatial_realised_columns(step28, events_df, prior_f)
    if spatial_roll is not None:
        roll = roll.merge(spatial_roll, on="locationID", how="outer")
    return loc.merge(roll, on="locationID", how="left")


# ---------------------------------------------------------------------------
# 3 — location-level coverage curves by size bin
# ---------------------------------------------------------------------------
def build_curves(rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for label, lo, hi in LOCATION_SIZE_BINS:
        # Representative K: midpoint of the bin (cap the top bin at 500 total
        # fertile plants for display purposes)
        rep = (lo + min(hi, 500)) / 2
        K = max(2 * int(round(rep)), 2)
        for M in MOTHERS_GRID:
            e = expected_coverage(K, int(M))
            sim = simulate_coverage(K, int(M), n_sim=N_SIM, rng=rng)
            rows.append({
                "bin":       label,
                "N_rep":     int(round(rep)),
                "K":         K,
                "M_mothers": int(M),
                "E_cov":     e,
                "sim_mean":  sim.mean(),
                "sim_lo":    np.quantile(sim, 0.025),
                "sim_hi":    np.quantile(sim, 0.975),
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 4 — figure
# ---------------------------------------------------------------------------
BIN_COLOURS = {
    "<20":     "#1b7837",
    "20-50":   "#7fbc41",
    "51-100":  "#dfc27d",
    "101-250": "#bf812d",
    ">250":    "#8c510a",
}


def plot_curves_empirical(prior_f: np.ndarray, per_location: pd.DataFrame,
                          out_png: Path, out_pdf: Path,
                          seeds_per_mother: int = DEFAULT_SEEDS_PER_MOTHER):
    """Two-panel location-level view under the P1 empirical prior.

    LEFT: theoretical curve — E[Fg coverage] vs total allele draws A.
          Reference lines at 90 %-Fgs target and at A_target under P1.

    RIGHT: per-location horizontal bars — each location's achieved Fg
          coverage given its M_achievable and the default seeds/mother from
          Step 28. Bars sorted by achieved coverage, coloured by
          N_fertile bin, with a vertical line at 0.90 marking the target.
    """
    A_grid = np.arange(1, 2000, 5)
    cov = np.array([expected_fg_coverage_empirical(prior_f, int(a))
                    for a in A_grid])
    A_target = alleles_for_fg_coverage_empirical(prior_f, EXPECTED_COVERAGE)

    fig, (ax_curve, ax_bars) = plt.subplots(
        1, 2, figsize=(14.5, 8.0),
        gridspec_kw={"width_ratios": [1.0, 1.3]},
    )

    # ------------------- LEFT: theoretical curve --------------------------
    ax_curve.plot(A_grid, cov, color="#1f78b4", lw=2.4,
                  label=f"E[Fg coverage] under P1 (K_fg = {len(prior_f)})")
    ax_curve.axhline(EXPECTED_COVERAGE, color="#333", ls="--", lw=1.0,
                     alpha=0.6)
    ax_curve.axvline(A_target, color="#b2182b", ls=":", lw=1.4)
    ax_curve.text(A_target + 15, 0.05,
                  f"A_target = {A_target} allele draws\n"
                  f"= 90 %-Fgs coverage",
                  fontsize=9, color="#b2182b", va="bottom")
    ax_curve.set_xlim(0, 1500); ax_curve.set_ylim(0, 1.02)
    ax_curve.set_xlabel("Total allele draws per location "
                        f"A = M × (PLOIDY + PATERNAL_PER_SEED × n_seeds) "
                        f"= M × ({PLOIDY} + {PATERNAL_ALLELES_PER_SEED}·n_seeds)")
    ax_curve.set_ylabel("Expected fraction of the 32 Fgs observed")
    ax_curve.set_title("Prediction curve under P1 prior", fontsize=11)
    ax_curve.legend(loc="lower right", fontsize=9, frameon=True)
    ax_curve.spines["top"].set_visible(False)
    ax_curve.spines["right"].set_visible(False)

    # ------------------- RIGHT: per-location bars -------------------------
    # sort by realised (per-mother-budget-aware) coverage under Rule 2 if
    # available; fall back to the uniform-seeds column otherwise.
    cov_col = ("exp_Fg_cov_realised_miss_P1"
               if "exp_Fg_cov_realised_miss_P1" in per_location.columns
               else "exp_Fg_cov_with_seeds_P1")
    df = (per_location
          .dropna(subset=[cov_col])
          .sort_values(cov_col, ascending=True)
          .reset_index(drop=True))

    def _bin_of(n_fert: int) -> str:
        for label, lo, hi in LOCATION_SIZE_BINS:
            if lo <= n_fert <= hi:
                return label
        return LOCATION_SIZE_BINS[-1][0]

    df["bin"] = df["total_n_fertile"].apply(_bin_of)
    colours = df["bin"].map(BIN_COLOURS)
    y = np.arange(len(df))
    ax_bars.barh(y, df[cov_col],
                 color=colours, edgecolor="white", height=0.75)
    ax_bars.axvline(EXPECTED_COVERAGE, color="#b2182b", ls="--", lw=1.4)
    ax_bars.text(EXPECTED_COVERAGE + 0.005, len(df) - 0.5,
                 " 90 %-Fgs target",
                 color="#b2182b", fontsize=9, va="top", ha="left")
    ax_bars.set_yticks(y)
    ax_bars.set_yticklabels(df["locationCode"], fontsize=6.5)
    ax_bars.set_xlim(0, 1.02); ax_bars.set_ylim(-0.7, len(df) - 0.3)
    cov_label = ("realised per-mother seed budgets (Rule 2)"
                 if cov_col == "exp_Fg_cov_realised_miss_P1"
                 else f"uniform {seeds_per_mother} seeds / mother")
    ax_bars.set_xlabel(f"Achieved Fg coverage under P1  —  {cov_label}")
    ax_bars.set_title(
        f"Per-location result — "
        f"{(df[cov_col] >= EXPECTED_COVERAGE).sum()}"
        f" / {len(df)} locations hit target",
        fontsize=11,
    )
    ax_bars.spines["top"].set_visible(False)
    ax_bars.spines["right"].set_visible(False)

    # bin-colour legend (with location counts per bin)
    from matplotlib.patches import Patch
    counts = df["bin"].value_counts().reindex(
        [b[0] for b in LOCATION_SIZE_BINS], fill_value=0)
    handles = [
        Patch(facecolor=BIN_COLOURS[label],
              label=f"N_fertile {label}  (n = {counts[label]})")
        for label, _, _ in LOCATION_SIZE_BINS if counts[label] > 0
    ]
    ax_bars.legend(handles=handles, loc="lower right",
                    fontsize=8, frameon=True)

    fig.suptitle("Step 29 — location-level SRK sampling under P1 empirical prior",
                 fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)


def plot_recommended_seeds_per_mother(per_location: pd.DataFrame,
                                      out_png: Path, out_pdf: Path):
    """Per-location bar chart of the recommended seeds/mother to hit the
    90 %-of-species-wide-SRK-alleles target, given each location's
    number of mothers with seed records in the DB (permit-realistic M).

    Bars are coloured by achievability tier:
      green   : ≤ 29 seeds/mother          (fits Step 28 Rule 2 floor)
      amber   : 30 – 100 seeds/mother      (achievable with effort)
      red     : > 100 seeds/mother         (unrealistic — location too small)

    Vertical reference lines: 29 (Rule 2 floor), 100 (practical ceiling).
    Bars > 300 are capped for display with an "off-scale" mark, and the
    real value is printed to the right.
    """
    # Prefer the permit-realistic column when Step 28's mother rollup has
    # been merged in; fall back to the theoretical ceiling only if it is
    # missing.
    col = ("seeds_per_mother_for_90pct_P1_actual"
           if "seeds_per_mother_for_90pct_P1_actual" in per_location.columns
           else "seeds_per_mother_for_90pct_P1")
    df = per_location.dropna(subset=[col]).copy()
    df = df.sort_values(col, ascending=True).reset_index(drop=True)
    seeds = df[col].astype(float).values

    def _tier_colour(n: float) -> str:
        if n <= 29:      return "#1b7837"     # green
        if n <= 100:     return "#e08214"     # amber
        return "#b2182b"                       # red

    colours = [_tier_colour(n) for n in seeds]

    CAP = 300
    display_seeds = np.minimum(seeds, CAP)
    off_scale = seeds > CAP

    fig, ax = plt.subplots(figsize=(9.5, max(4.5, 0.18 * len(df) + 2)))
    y = np.arange(len(df))
    ax.barh(y, display_seeds, color=colours, edgecolor="white", height=0.75)

    # off-scale annotations
    for i, (val, is_off) in enumerate(zip(seeds, off_scale)):
        if is_off:
            ax.text(CAP + 2, i, f"→ {int(val)}",
                    fontsize=8, color="#b2182b",
                    va="center", ha="left")

    RULE2_SEEDS = DEFAULT_SEEDS_PER_MOTHER
    ax.axvline(RULE2_SEEDS, color="#333", ls="--", lw=1.0, alpha=0.6)
    ax.axvline(100,         color="#333", ls=":",  lw=1.0, alpha=0.6)
    ax.text(RULE2_SEEDS, len(df) - 0.5,
            f" {RULE2_SEEDS} (Rule 2 floor, tetraploid)",
            fontsize=8, color="#333", va="top", ha="left")
    ax.text(100, len(df) - 0.5, " 100 (practical ceiling)",
            fontsize=8, color="#333", va="top", ha="left")

    # Y-tick labels include the sampled-mother count so readers see the
    # denominator behind each bar directly.
    if "M_actual_in_step28" in df.columns:
        labels = [f"{code}   (n = {int(m)})"
                  for code, m in zip(df["locationCode"],
                                     df["M_actual_in_step28"].fillna(0))]
    else:
        labels = df["locationCode"].tolist()
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlim(0, CAP + 40)
    ax.set_ylim(-0.7, len(df) - 0.3)
    ax.set_xlabel(
        "Recommended number of seeds to genotype per mother plant  "
        "(to detect 90 % of the 32 species-wide SRK alleles at the "
        "location)"
    )
    ax.set_title(
        "Per-location seed-genotyping recommendation across LEPA  "
        "— one bar per location, based on the mothers actually in "
        "the seed bank",
        fontsize=12,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor="#1b7837",
              label=f"≤ {RULE2_SEEDS} seeds / mother  "
                    f"(n = {(seeds <= RULE2_SEEDS).sum()})"),
        Patch(facecolor="#e08214",
              label=f"{RULE2_SEEDS + 1} – 100 seeds / mother  "
                    f"(n = {((seeds > RULE2_SEEDS) & (seeds <= 100)).sum()})"),
        Patch(facecolor="#b2182b",
              label=f"> 100 seeds / mother  (n = {(seeds > 100).sum()})"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=9, frameon=True)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)


def plot_curves(curves: pd.DataFrame, out_png: Path, out_pdf: Path):
    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    for label, _, _ in LOCATION_SIZE_BINS:
        sub = curves[curves["bin"] == label]
        colour = BIN_COLOURS[label]
        ax.fill_between(sub["M_mothers"], sub["sim_lo"], sub["sim_hi"],
                        color=colour, alpha=0.18, linewidth=0)
        ax.plot(sub["M_mothers"], sub["E_cov"],
                color=colour, lw=2.0,
                label=f"Location N_fertile {label}  (K={int(sub['K'].iloc[0])})")
    ax.axhline(EXPECTED_COVERAGE, color="#333333",
               ls="--", lw=1.0, alpha=0.6)
    ax.text(MOTHERS_GRID.max(), EXPECTED_COVERAGE + 0.01,
            f"target = {int(EXPECTED_COVERAGE*100)} %  expected coverage",
            fontsize=9, color="#333333", ha="right", va="bottom")
    ax.set_xlim(1, MOTHERS_GRID.max())
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Number of mothers sampled per location")
    ax.set_ylabel("Expected maternal-allele coverage")
    ax.set_title("Location-level SRK sampling — coverage vs mothers sampled",
                 fontsize=12)
    ax.legend(loc="lower right", fontsize=9, frameon=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 5 — main
# ---------------------------------------------------------------------------
def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--year", type=int, default=None,
                    help="Restrict to events of the given eventDate year "
                         "(e.g. 2025). Default: all years combined.")
    args = ap.parse_args()

    tables_dir  = Path(DEFAULT_TABLES); tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = Path(DEFAULT_FIGURES); figures_dir.mkdir(parents=True, exist_ok=True)

    events = load_events_by_location(DEFAULT_DB, year=args.year)
    scope = f"year {args.year}" if args.year is not None else "all years"
    print(f"[step29] Loaded {len(events)} events ({scope}) across "
          f"{events['locationID'].nunique()} LEPA locations "
          f"(N_fertile per event 1-{events['n_fertile'].max()}).")

    per_event, per_location = build_designs(events)

    # Empirical P1 prior — augments the per-location table with 90 %-of-Fgs
    # target under the real (skewed) LEPA Fg frequency vector.
    prior_f = None
    if DEFAULT_PRIOR_TSV.exists():
        prior_f = load_p1_prior(DEFAULT_PRIOR_TSV)
        per_location = augment_with_empirical_prior(
            per_location, prior_f,
            seeds_per_mother=DEFAULT_SEEDS_PER_MOTHER,
        )
        A_target = int(per_location["A_target_90pct_P1"].iloc[0])
        print(f"[step29] P1 empirical prior loaded ({len(prior_f)} Fgs). "
              f"90 %-of-Fgs target = {A_target} total allele draws / location.")
    else:
        print(f"[step29] WARNING: {DEFAULT_PRIOR_TSV} not found — empirical "
              "prior columns skipped.")

    # Realised (per-mother-budget-aware) coverage — reads the ACTUAL
    # `n_achievable_*` from Step 28 so seeds/mother varies with each mother's
    # `germplasmQuantityEstimate` instead of assuming a uniform value.
    if prior_f is not None and DEFAULT_STEP28_TSV.exists():
        per_location = augment_with_realised_step28(
            per_location, events, DEFAULT_STEP28_TSV, prior_f,
        )
        pct_hit_real = (
            per_location["exp_Fg_cov_realised_miss_P1"]
            >= EXPECTED_COVERAGE
        ).mean()
        print(f"[step29] Realised (per-mother-budget) coverage under Rule 2: "
              f"{pct_hit_real:.0%} of locations hit the 90 %-Fgs target.")
    elif prior_f is not None:
        print(f"[step29] WARNING: {DEFAULT_STEP28_TSV} not found — "
              "per-mother-budget-aware columns skipped. Run step28 first.")

    per_event_path = tables_dir / "step29_sampling_per_event.tsv"
    per_location_path = tables_dir / "step29_sampling_per_location.tsv"
    per_event.to_csv(per_event_path, sep="\t", index=False)
    per_location.to_csv(per_location_path, sep="\t", index=False)
    print(f"[step29] Wrote {per_event_path}")
    print(f"[step29] Wrote {per_location_path}")

    # Field-team recipe — Phase-A actionable per-germplasmID table.
    recipe = build_field_team_recipe(DEFAULT_STEP28_TSV, events)
    if len(recipe):
        recipe_path = tables_dir / "step29_field_team_sampling_recipe.tsv"
        recipe.to_csv(recipe_path, sep="\t", index=False)
        n_r1 = int(recipe["rationale"].str.startswith("Rule 1").sum())
        n_r2 = int(recipe["rationale"].str.startswith("Rule 2").sum())
        n_bl = int(recipe["rationale"].str.startswith("budget").sum())
        print(f"[step29] Wrote {recipe_path} "
              f"({len(recipe)} mothers: {n_r1} Rule 1, "
              f"{n_r2} Rule 2, {n_bl} budget-limited).")

    n_short = (per_location["M_achievable_location"]
               < per_location["M_recommended_location"]).sum()
    print(f"[step29] {n_short}/{len(per_location)} locations cannot reach the "
          f"uniform-prior recommended M via proportional allocation.")

    if "exp_Fg_cov_with_seeds_P1" in per_location.columns:
        pct_hit = (per_location["exp_Fg_cov_with_seeds_P1"] >= EXPECTED_COVERAGE).mean()
        print(f"[step29] Under empirical P1 prior + {DEFAULT_SEEDS_PER_MOTHER} "
              f"seeds/mother: {pct_hit:.0%} of locations hit the 90 %-Fgs "
              f"target with their achievable M.")

    rng = np.random.default_rng(2029)
    curves = build_curves(rng)
    curves_path = tables_dir / "step29_location_coverage_curves.tsv"
    curves.to_csv(curves_path, sep="\t", index=False)
    print(f"[step29] Wrote {curves_path}")
    # The coverage-curve figures and per-mother tier figure were dropped
    # in the 2026-09-30 declutter — superseded by
    # step29c_sampling_comparison (target vs already collected) and by
    # step30_A_diversity_unbiased_vs_sampling (finite-population truth
    # vs sampling recovery). Coverage curves TSV kept for audit only.


if __name__ == "__main__":
    main()
