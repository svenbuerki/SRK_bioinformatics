#!/usr/bin/env python3
"""Step 28 — Per-mother seed-sampling design for SRK diversity estimation.

Answers, per mother plant, how many of her seeds must be genotyped so that the
pollen SRK allele pool she was exposed to is characterised at a target
resolution. Uses the "reverse random-mating" logic: assume compatible pollen is
sampled uniformly from the (organismQuantityFertile - 1) other fertile plants
at her event; compute expected sire-allele coverage as a function of the number
of genotyped seeds.

Two decision rules are reported side-by-side per mother:

  * n_expected_cov  — seeds needed for 90 % expected sire-allele coverage,
    under uniform pollen weights with K = 2 * (N_fertile - 1) potential sire
    alleles. Closed-form: n = log(1 - 0.9) / log(1 - 1/K).

  * n_miss_prob     — seeds needed to be 95 % sure of detecting any pollen
    allele contributing at least 10 % of siring events. Closed-form (K-
    independent): n = log(0.05) / log(1 - 0.10) ≈ 29.

Both are compared against `germplasmQuantityEstimate` (with Low/Upr bracket),
flagging mothers whose seed budget is smaller than the recommended n.

Inputs (read-only)
    /Users/sven/Documents/Current_projects/LEPA_fieldwork_protocol/SQL_DB/LEPA_SQL.db
    Tables used:
        Germplasm     — germplasmID, occurrenceID, germplasmQuantityEstimate,
                        germplasmQuantityEstimateLow, germplasmQuantityEstimateUpr
        Occurrences   — occurrenceID, eventID, taxonID
        Events        — eventID, organismQuantityFertile
        Taxonomy      — taxonID, genus, specificEpithet  (filter to LEPA)

Outputs
    Tables/Phase5/step28_seed_sampling_per_mother.tsv
        Per-mother design table (785 LEPA mothers).
    Tables/Phase5/step28_coverage_curves_by_Nfertile.tsv
        Analytical + simulated coverage curves at fixed N_fertile bins.
    figures/Phase5/step28_coverage_curves.png/pdf
        Expected coverage vs #seeds, one line per N_fertile bin, with
        simulation-CI ribbons.
    figures/Phase5/step28_per_mother_budget.png/pdf
        Seed budget vs recommended n, coloured by N_fertile bin — flags
        mothers whose budget is insufficient.
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

# Rule 1 — expected-coverage target
EXPECTED_COVERAGE = 0.90
# Rule 2 — miss-probability guarantee at fixed relative frequency threshold
P_MIN = 0.10           # detect any allele siring at least 10 % of offspring
MISS_ALPHA = 0.05      # with at least 95 % probability

# Biological ceiling on distinct SRK alleles at any single mating event:
# no local pollen pool can carry more distinct Fgs than the species-wide
# functional-SRK inventory (32 Fgs from the Canu-amplicon preliminary
# study). Used to cap the per-bin K in the coverage-curves figure so
# large events don't visually promise more allele diversity than the
# species actually holds.
K_SPECIES_FG = 32

# Spatial mating-neighbourhood radii (metres). R = 10 is the primary
# assumption; R = 25 and R = 50 are sensitivity checks.
SPATIAL_RADII_M = (10.0, 25.0, 50.0)
PRIMARY_SPATIAL_RADIUS_M = 25.0

# For the aggregate curves figure
NFERTILE_BINS = [
    ("1-2",   1,  2),
    ("3-5",   3,  5),
    ("6-10",  6, 10),
    ("11-20", 11, 20),
    ("21-50", 21, 50),
    (">50",   51, 10_000),
]
N_SIM = 2_000                     # bootstrap draws per (K, n) point
SEEDS_GRID = np.arange(1, 201)    # x-axis for the curves


# ---------------------------------------------------------------------------
# Analytical helpers
# ---------------------------------------------------------------------------
def n_for_expected_coverage(K: int, target: float = EXPECTED_COVERAGE) -> int:
    """Smallest n such that E[distinct sire alleles]/K >= target, uniform p.

    Under uniform p_j = 1/K, E[D]/K = 1 - (1 - 1/K)^n. Solve for n:
        n = ceil( log(1 - target) / log(1 - 1/K) )
    K must be at least 1. For K == 1 the mother sees exactly one father, whose
    2 alleles are both observed once any seed is genotyped.
    """
    if K <= 1:
        return 1
    p = 1.0 / K
    return int(math.ceil(math.log(1 - target) / math.log(1 - p)))


def n_for_miss_probability(p_min: float = P_MIN,
                           alpha: float = MISS_ALPHA) -> int:
    """Smallest n such that any allele of frequency >= p_min is missed with
    probability <= alpha. Closed-form and K-independent:
        n = ceil( log(alpha) / log(1 - p_min) )
    """
    return int(math.ceil(math.log(alpha) / math.log(1 - p_min)))


def expected_coverage(K: int, n: int) -> float:
    """E[distinct alleles observed]/K at sample size n, uniform pollen."""
    if K <= 1:
        return 1.0
    return 1.0 - (1.0 - 1.0 / K) ** n


def simulate_coverage(K: int, n: int, n_sim: int = N_SIM,
                       rng: np.random.Generator | None = None) -> np.ndarray:
    """Bootstrap coverage: draw n multinomial pollen alleles from K uniform
    weights, count distinct hits, repeat n_sim times. Returns array of
    coverage fractions (each in [0, 1]).
    """
    if K <= 1:
        return np.ones(n_sim)
    rng = rng or np.random.default_rng(2026)
    draws = rng.integers(0, K, size=(n_sim, n))
    # count distinct values per row
    def n_unique(row):
        return np.unique(row).size
    counts = np.apply_along_axis(n_unique, 1, draws)
    return counts / K


# ---------------------------------------------------------------------------
# Spatial helpers — event mating-neighbourhood (used by Step 28 + Step 29)
# ---------------------------------------------------------------------------
EARTH_RADIUS_M = 6_371_008.8


def haversine_meters(lat1: np.ndarray, lon1: np.ndarray,
                     lat2: np.ndarray, lon2: np.ndarray) -> np.ndarray:
    """Great-circle distance in metres between arrays of (lat, lon) points.
    Inputs in decimal degrees, broadcastable. Returns metres."""
    lat1_r = np.deg2rad(lat1); lon1_r = np.deg2rad(lon1)
    lat2_r = np.deg2rad(lat2); lon2_r = np.deg2rad(lon2)
    dlat = lat2_r - lat1_r
    dlon = lon2_r - lon1_r
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1_r) * np.cos(lat2_r) \
        * np.sin(dlon / 2.0) ** 2
    return 2 * EARTH_RADIUS_M * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def load_all_events(db_path: Path, year: int | None = None) -> pd.DataFrame:
    """All LEPA events with coordinates and organismQuantityFertile (parsed
    leniently — the DB stores free-text notes in that field for some rows).
    Includes the event's own N_fertile so we can build the mating pool.

    If `year` is given, restrict to events of that eventDate year — the
    spatial neighbourhood should never mix survey years (§ 1.5 of the
    design doc).
    """
    con = sqlite3.connect(str(db_path))
    try:
        df = pd.read_sql_query(
            """
            SELECT DISTINCT
                e.eventID                 AS eventID,
                e.locationID              AS locationID,
                e.eventDate               AS eventDate,
                e.eventDecimalLatitude    AS lat,
                e.eventDecimalLongitude   AS lon,
                e.organismQuantityFertile AS n_fertile_raw
            FROM Events e
            JOIN Occurrences o ON o.eventID = e.eventID
            JOIN Taxonomy    t ON o.taxonID = t.taxonID
            WHERE t.genus='Lepidium' AND t.specificEpithet='papilliferum'
              AND e.eventDecimalLatitude  IS NOT NULL
              AND e.eventDecimalLongitude IS NOT NULL
              AND e.organismQuantityFertile IS NOT NULL
              AND (o.provenance IS NULL OR o.provenance = 'in situ')
            """,
            con,
        )
    finally:
        con.close()
    df["event_year"] = _extract_year(df["eventDate"])
    if year is not None:
        df = df[df["event_year"] == int(year)].copy()
    # Parse the leading integer of organismQuantityFertile (free-text safe).
    n = (df["n_fertile_raw"].astype(str)
         .str.extract(r"(\d+)", expand=False)
         .replace("", np.nan))
    df["n_fertile"] = pd.to_numeric(n, errors="coerce")
    df = df.dropna(subset=["n_fertile"]).copy()
    df["n_fertile"] = df["n_fertile"].astype(int)
    # A handful of legacy rows store lat/lon as free text like
    # "43, 39 10 28" — coerce leniently and drop anything unparseable.
    df["lat"] = pd.to_numeric(df["lat"], errors="coerce")
    df["lon"] = pd.to_numeric(df["lon"], errors="coerce")
    df = df.dropna(subset=["lat", "lon"]).copy()
    # LEPA (SW Idaho) sits near 43 N, -116 E; drop obvious dumps outside a
    # generous North-America bounding box so a mistyped row does not push
    # a neighbourhood off-planet.
    ok = (df["lat"].between(30, 55)) & (df["lon"].between(-125, -100))
    df = df[ok].copy()
    return df[["eventID", "locationID", "event_year",
                "lat", "lon", "n_fertile"]]


def compute_pairwise_distances(events: pd.DataFrame) -> np.ndarray:
    """Full haversine distance matrix (metres) over all events."""
    lat = events["lat"].values
    lon = events["lon"].values
    return haversine_meters(
        lat[:, None], lon[:, None],
        lat[None, :], lon[None, :],
    )


def spatial_neighborhood_stats(events: pd.DataFrame,
                                dist_m: np.ndarray,
                                radius_m: float) -> pd.DataFrame:
    """For each event, count how many *other* events lie within `radius_m`
    metres and sum their N_fertile. Returns a DataFrame indexed by eventID
    with columns:
        N_neighbor_events_<R>m
        N_fertile_in_neighbors_<R>m
        N_compatible_spatial_<R>m
            = (self N_fertile - 1) + neighbors' N_fertile
        K_spatial_<R>m = 2 * N_compatible_spatial_<R>m
    """
    within = (dist_m > 0) & (dist_m <= radius_m)
    n_neighbors = within.sum(axis=1)
    n_fert_neighbors = (within * events["n_fertile"].values[None, :]).sum(axis=1)
    n_comp_spatial = (events["n_fertile"].values - 1) + n_fert_neighbors
    K_spatial = 2 * np.clip(n_comp_spatial, 0, None)
    R = int(round(radius_m))
    return pd.DataFrame({
        "eventID": events["eventID"].values,
        f"N_neighbor_events_{R}m":       n_neighbors.astype(int),
        f"N_fertile_in_neighbors_{R}m":  n_fert_neighbors.astype(int),
        f"N_compatible_spatial_{R}m":    n_comp_spatial.astype(int),
        f"K_spatial_{R}m":               K_spatial.astype(int),
    })


def build_spatial_frame(db_path: Path,
                        radii_m: tuple[float, ...] = SPATIAL_RADII_M,
                        year: int | None = None,
                        ) -> pd.DataFrame:
    """Return one row per event with all spatial-neighbourhood columns for
    each radius in `radii_m`. Also carries lat / lon / n_fertile.
    Restricts to events of `year` if given.
    """
    events = load_all_events(db_path, year=year)
    dist = compute_pairwise_distances(events)
    frame = events[["eventID", "locationID", "event_year",
                    "lat", "lon", "n_fertile"]].copy()
    for R in radii_m:
        stats = spatial_neighborhood_stats(events, dist, R)
        frame = frame.merge(stats, on="eventID", how="left")
    return frame


# ---------------------------------------------------------------------------
# 1 — pull the per-mother table from the LEPA SQL DB
# ---------------------------------------------------------------------------
def _extract_year(date_series: pd.Series) -> pd.Series:
    """LEPA `eventDate` is stored as MM-DD-YYYY (free text). Extract the
    four-digit year. Unparseable dates return NaN."""
    return (date_series.astype(str)
            .str.extract(r"(\d{4})$", expand=False)
            .astype("Int64"))


def load_mother_table(db_path: Path, year: int | None = None) -> pd.DataFrame:
    con = sqlite3.connect(str(db_path))
    try:
        df = pd.read_sql_query(
            """
            SELECT
                o.occurrenceID                 AS occurrenceID,
                o.eventID                      AS eventID,
                e.eventDate                    AS eventDate,
                e.organismQuantityFertile      AS n_fertile,
                g.germplasmID                  AS germplasmID,
                g.germplasmQuantityEstimate    AS seeds_est,
                g.germplasmQuantityEstimateLow AS seeds_low,
                g.germplasmQuantityEstimateUpr AS seeds_upr
            FROM Germplasm g
            JOIN Occurrences o ON g.occurrenceID = o.occurrenceID
            JOIN Events      e ON o.eventID      = e.eventID
            JOIN Taxonomy    t ON o.taxonID      = t.taxonID
            WHERE g.germplasmQuantityEstimate IS NOT NULL
              AND e.organismQuantityFertile   IS NOT NULL
              AND t.genus = 'Lepidium'
              AND t.specificEpithet = 'papilliferum'
              AND g.biologicalStatus = 'Wild'
              AND (o.provenance IS NULL OR o.provenance = 'in situ')
            """,
            con,
        )
    finally:
        con.close()
    df["event_year"] = _extract_year(df["eventDate"])
    if year is not None:
        n_before = len(df)
        df = df[df["event_year"] == int(year)].copy()
        print(f"[step28] --year {year}: kept {len(df)} of {n_before} mothers.")
    return df


# ---------------------------------------------------------------------------
# 2 — per-mother design table
# ---------------------------------------------------------------------------
def build_per_mother_table(df: pd.DataFrame) -> pd.DataFrame:
    # No SI filter yet: N_compatible = organismQuantityFertile - 1
    n_comp = (df["n_fertile"].astype(float) - 1).clip(lower=0).astype(int)
    K = (2 * n_comp).clip(lower=1).astype(int)   # 2 SRK alleles per mate

    n_exp = np.array([n_for_expected_coverage(k) for k in K])
    n_miss = n_for_miss_probability()    # constant
    seeds_est = df["seeds_est"].astype(float)
    seeds_low = df["seeds_low"].astype(float)
    seeds_upr = df["seeds_upr"].astype(float)

    n_achievable_exp  = np.minimum(n_exp,  seeds_est).round().astype(int)
    n_achievable_miss = np.minimum(n_miss, seeds_est).round().astype(int)
    budget_ok_exp  = seeds_est >= n_exp
    budget_ok_miss = seeds_est >= n_miss

    # --- Seed-production-aware quantities --------------------------------
    # (1) Achieved coverage per mother, given how many of her seeds she
    #     actually has: 1 - (1 - 1/K)^n_use, where n_use = min(n_rec, S).
    #     Under uniform pollen weights.
    K_arr = K.astype(float).values
    n_use_exp  = n_achievable_exp.astype(float).values
    n_use_miss = n_achievable_miss.astype(float).values
    S = seeds_est.astype(float).values
    with np.errstate(divide="ignore", invalid="ignore"):
        one_minus_1_over_K = np.where(K_arr > 0, 1.0 - 1.0 / K_arr, 0.0)
    achieved_cov_exp  = 1.0 - one_minus_1_over_K ** n_use_exp
    achieved_cov_miss = 1.0 - one_minus_1_over_K ** n_use_miss

    # (2) Expected distinct paternal SRK alleles present in the mother's
    #     complete seed lot (the biological ceiling — no amount of
    #     genotyping can exceed this): K * (1 - (1 - 1/K)^S).
    expected_distinct_in_lot = K_arr * (1.0 - one_minus_1_over_K ** S)

    out = pd.DataFrame({
        "occurrenceID":          df["occurrenceID"],
        "eventID":               df["eventID"],
        "germplasmID":           df["germplasmID"],
        "n_fertile":             df["n_fertile"].astype(int),
        "n_compatible_event":    n_comp,
        "K_event":               K,
        "n_compatible":          n_comp,       # legacy alias — same as *_event
        "K_potential_sire_alleles": K,          # legacy alias
        "seeds_est":             seeds_est.round().astype(int),
        "seeds_low":             seeds_low.round().astype("Int64"),
        "seeds_upr":             seeds_upr.round().astype("Int64"),
        "n_expected_cov_90":     n_exp,
        "n_miss_prob_10pct_95":  n_miss,
        "n_achievable_exp":      n_achievable_exp,
        "n_achievable_miss":     n_achievable_miss,
        "budget_ok_exp":         budget_ok_exp,
        "budget_ok_miss":        budget_ok_miss,
        # NEW seed-production-aware columns
        "achieved_coverage_exp":   np.round(achieved_cov_exp,  4),
        "achieved_coverage_miss":  np.round(achieved_cov_miss, 4),
        "expected_distinct_alleles_in_seed_lot":
                                   np.round(expected_distinct_in_lot, 2),
    })
    return out.sort_values(["n_fertile", "occurrenceID"]).reset_index(drop=True)


def augment_per_mother_with_spatial(per_mother: pd.DataFrame,
                                     spatial_frame: pd.DataFrame
                                     ) -> pd.DataFrame:
    """Join the per-mother table (indexed by eventID) with the event-level
    spatial frame, then add — for each SPATIAL_RADII_M value — the mother's
    spatial K, the Rule 1 recommended n at that spatial K, the achievable n
    given her real seed budget, and the achieved coverage.

    Rule 2 does not depend on K, so n_miss_prob_10pct_95 is unchanged.
    """
    keep_cols = ["eventID", "lat", "lon"] + [
        c for c in spatial_frame.columns if c.startswith((
            "N_neighbor_events_",
            "N_fertile_in_neighbors_",
            "N_compatible_spatial_",
            "K_spatial_",
        ))
    ]
    df = per_mother.merge(spatial_frame[keep_cols], on="eventID", how="left")

    for R in SPATIAL_RADII_M:
        R_int = int(round(R))
        K_col = f"K_spatial_{R_int}m"
        K_arr = df[K_col].fillna(df["K_event"]).astype(int).values

        # Rule 1 recommended n at spatial K
        n_exp_col   = f"n_expected_cov_90_spatial_{R_int}m"
        n_ach_col   = f"n_achievable_exp_spatial_{R_int}m"
        cov_col     = f"achieved_coverage_exp_spatial_{R_int}m"
        exp_lot_col = f"expected_distinct_alleles_in_seed_lot_spatial_{R_int}m"
        budget_col  = f"budget_ok_exp_spatial_{R_int}m"

        n_exp_arr = np.array([n_for_expected_coverage(int(k)) for k in K_arr])
        S = df["seeds_est"].astype(float).values
        n_use = np.minimum(n_exp_arr, S).astype(int)
        with np.errstate(divide="ignore", invalid="ignore"):
            one_minus_1_over_K = np.where(K_arr > 0, 1.0 - 1.0 / K_arr, 0.0)
        achieved_cov = 1.0 - one_minus_1_over_K ** n_use
        exp_in_lot = K_arr * (1.0 - one_minus_1_over_K ** S)

        df[n_exp_col]   = n_exp_arr
        df[n_ach_col]   = n_use
        df[cov_col]     = np.round(achieved_cov, 4)
        df[exp_lot_col] = np.round(exp_in_lot, 2)
        df[budget_col]  = S >= n_exp_arr
    return df


# ---------------------------------------------------------------------------
# 3 — aggregate coverage curves by N_fertile bin
# ---------------------------------------------------------------------------
def build_curves_by_bin(rng: np.random.Generator) -> pd.DataFrame:
    """One row per (bin, n) with analytical E[coverage] + simulation CI.

    For each bin we use the midpoint of the bin's N_fertile range as the
    representative pool size. The pool of *distinct* pollen-donor SRK
    alleles is capped at the species-wide ceiling (K_SPECIES_FG = 32
    Fgs) because a local mating neighbourhood cannot carry more
    functional alleles than the species actually holds — a 100-plant
    event doesn't produce 200 distinct alleles, it produces at most 32.
    Both the fraction (E_cov, ∈ [0, 1]) and the absolute allele count
    (E_alleles = E_cov · K) are stored so the figure can plot on the
    honest absolute scale while the TSV still carries the fraction.
    """
    rows = []
    for label, lo, hi in NFERTILE_BINS:
        # representative K: average of 2*(N-1) over the bin range, rounded,
        # then capped at the species-wide Fg ceiling.
        n_rep = (lo + min(hi, 200)) / 2      # cap >50 bin at 200 for display
        K_uncapped = max(2 * (int(round(n_rep)) - 1), 1)
        K = min(K_uncapped, K_SPECIES_FG)
        for n in SEEDS_GRID:
            e = expected_coverage(K, int(n))
            sim = simulate_coverage(K, int(n), n_sim=N_SIM, rng=rng)
            rows.append({
                "bin":            label,
                "N_rep":          int(round(n_rep)),
                "K_uncapped":     K_uncapped,
                "K":              K,
                "n_seeds":        int(n),
                "E_cov":          e,
                "E_alleles":      e * K,
                "sim_mean":       sim.mean(),
                "sim_lo":         float(np.quantile(sim, 0.025)),
                "sim_hi":         float(np.quantile(sim, 0.975)),
                "sim_alleles_lo": float(np.quantile(sim, 0.025)) * K,
                "sim_alleles_hi": float(np.quantile(sim, 0.975)) * K,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 4 — figures
# ---------------------------------------------------------------------------
BIN_COLOURS = {
    "1-2":   "#1b7837",   # dark green — smallest
    "3-5":   "#7fbc41",
    "6-10":  "#f6e8c3",   # neutral
    "11-20": "#dfc27d",
    "21-50": "#bf812d",
    ">50":   "#8c510a",   # dark brown — largest
}


def _bin_of(n_fertile: int) -> str:
    for label, lo, hi in NFERTILE_BINS:
        if lo <= n_fertile <= hi:
            return label
    return ">50"


def plot_coverage_curves(curves: pd.DataFrame, out_png: Path, out_pdf: Path):
    """Theoretical coverage curves on the absolute-allele scale, with the
    operational 29-seed cap.

    Each curve shows the **expected number of distinct SRK alleles** a
    mother's seed lot would reveal as we increase the number of seeds
    genotyped, with the local pool of distinct alleles capped at the
    species-wide ceiling (32 Fgs). Plotting on the absolute scale keeps
    the curves honestly comparable across event sizes: small events
    saturate quickly because there are few alleles to find, not because
    they are easier to sample. The dashed grey horizontal line marks
    the 32-allele species ceiling; the vertical red line marks the
    29-seed operational cap (Rule 2); coloured dots on each curve show
    the number of distinct alleles each event-size bin **actually
    delivers per mother** at the 29-seed cap.
    """
    RULE2 = 29
    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    # Shade the "not-operational" region beyond the 29-seed cap.
    ax.axvspan(RULE2, SEEDS_GRID.max(),
               color="#f0f0f0", alpha=0.6, zorder=0)
    for label, _, _ in NFERTILE_BINS:
        sub = curves[curves["bin"] == label]
        colour = BIN_COLOURS[label]
        K_bin = int(sub["K"].iloc[0])
        ax.fill_between(sub["n_seeds"],
                        sub["sim_alleles_lo"], sub["sim_alleles_hi"],
                        color=colour, alpha=0.18, linewidth=0)
        ax.plot(sub["n_seeds"], sub["E_alleles"],
                color=colour, lw=2.0,
                label=f"{label} fertile plants at event  "
                      f"(local pool = {K_bin} distinct alleles)")
        # Mark the alleles delivered by the 29-seed cap — this is the
        # sampling-strategy dot the field team can point to.
        try:
            alleles_at_cap = float(
                sub.loc[sub["n_seeds"] == RULE2, "E_alleles"].iloc[0]
            )
        except IndexError:
            continue
        ax.scatter([RULE2], [alleles_at_cap], s=90, c=colour,
                   edgecolor="white", linewidth=1.2, zorder=5)
        ax.text(RULE2 + 3, alleles_at_cap,
                f" {alleles_at_cap:.1f} of {K_bin}",
                fontsize=9, color=colour, va="center", ha="left")

    # Species-wide SRK allele ceiling (32 Fgs from the preliminary study).
    ax.axhline(K_SPECIES_FG, color="#333333",
               ls="--", lw=1.0, alpha=0.6)
    ax.text(SEEDS_GRID.max() - 2, K_SPECIES_FG + 0.3,
            f"32 alleles — species-wide SRK ceiling",
            fontsize=9, color="#333333", ha="right", va="bottom")
    ax.axvline(RULE2, color="#b2182b", ls="-", lw=1.6, alpha=0.85)
    ax.text(RULE2 - 1.5, 0.5,
            f"29 seeds — operational cap (Rule 2)\n"
            f"the recipe never asks for more than this per mother",
            fontsize=9, color="#b2182b",
            ha="right", va="bottom")
    ax.set_xlim(1, SEEDS_GRID.max())
    ax.set_ylim(0, K_SPECIES_FG + 2)
    ax.set_xlabel("Number of seeds genotyped per mother")
    ax.set_ylabel("Expected number of distinct SRK alleles detected")
    ax.set_title(
        "Per-mother seed sampling — distinct SRK alleles revealed at each event size",
        fontsize=12,
    )
    ax.legend(loc="lower right", fontsize=9, frameon=True, title="Legend")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)


def plot_per_mother_budget(per_mother: pd.DataFrame,
                           out_png: Path, out_pdf: Path):
    df = per_mother.copy()
    df["bin"] = df["n_fertile"].apply(_bin_of)

    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    for label, _, _ in NFERTILE_BINS:
        sub = df[df["bin"] == label]
        if sub.empty:
            continue
        ax.scatter(
            sub["n_expected_cov_90"], sub["seeds_est"],
            s=18, c=BIN_COLOURS[label], alpha=0.75,
            edgecolor="white", linewidth=0.4,
            label=f"N_fertile {label}  (n={len(sub)})",
        )
    lim = max(1, int(df["seeds_est"].quantile(0.99)))
    ax.plot([0, lim], [0, lim], color="#333333", ls="--", lw=1.0, alpha=0.6)
    ax.text(lim, lim, " budget = required", fontsize=9,
            color="#333333", ha="left", va="bottom", rotation=45)
    ax.axvline(n_for_miss_probability(),
               color="#b2182b", ls=":", lw=1.2, alpha=0.7)
    ax.text(n_for_miss_probability() + 1, lim * 0.02,
            f"miss-prob rule n = {n_for_miss_probability()}",
            fontsize=9, color="#b2182b", ha="left", va="bottom")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, max(200, int(df["n_expected_cov_90"].max()) + 10))
    ax.set_ylim(1, lim * 1.2)
    ax.set_xlabel("Recommended seeds for 90 % expected coverage")
    ax.set_ylabel("Seed budget per mother (germplasmQuantityEstimate)")
    ax.set_title(
        "Per-mother seed budget vs required sample size — coloured by event size",
        fontsize=12,
    )
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
                    help="Restrict to mothers / events of the given "
                         "eventDate year (e.g. 2025). Also restricts the "
                         "spatial neighbourhood to same-year events. "
                         "Default: all years combined.")
    args = ap.parse_args()

    tables_dir  = Path(DEFAULT_TABLES); tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = Path(DEFAULT_FIGURES); figures_dir.mkdir(parents=True, exist_ok=True)

    raw = load_mother_table(DEFAULT_DB, year=args.year)
    scope = f"year {args.year}" if args.year is not None else "all years"
    print(f"[step28] Loaded {len(raw)} LEPA mothers ({scope}) with seed + "
          f"N_fertile across {raw['eventID'].nunique()} events "
          f"(N_fertile range {raw['n_fertile'].min():g}"
          f" – {raw['n_fertile'].max():g}).")

    per_mother = build_per_mother_table(raw)

    # Spatial extension: pull all LEPA events + coordinates once, build the
    # neighbourhood stats at each SPATIAL_RADII_M radius, and augment the
    # per-mother table so the field team can compare event-only vs
    # spatial-neighbourhood recommendations side by side.
    spatial_frame = build_spatial_frame(DEFAULT_DB, SPATIAL_RADII_M,
                                         year=args.year)
    per_mother = augment_per_mother_with_spatial(per_mother, spatial_frame)

    # Write the event-level spatial frame as its own reference table
    spatial_frame_path = tables_dir / "step28_events_spatial_neighborhood.tsv"
    spatial_frame.to_csv(spatial_frame_path, sep="\t", index=False)
    print(f"[step28] Wrote {spatial_frame_path}")

    per_mother_path = tables_dir / "step28_seed_sampling_per_mother.tsv"
    per_mother.to_csv(per_mother_path, sep="\t", index=False)
    print(f"[step28] Wrote {per_mother_path}")

    # Report how much the spatial extension changed K vs the event-only view.
    R = int(round(PRIMARY_SPATIAL_RADIUS_M))
    K_ev  = per_mother["K_event"].astype(int)
    K_sp  = per_mother[f"K_spatial_{R}m"].fillna(K_ev).astype(int)
    n_grew = int((K_sp > K_ev).sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = (K_sp / K_ev.replace(0, np.nan)).dropna()
    n_no_coord = int(per_mother[f"K_spatial_{R}m"].isna().sum())
    print(f"[step28] Spatial neighbourhood (R = {R} m): "
          f"{n_grew}/{len(per_mother)} mothers see a larger K than the "
          f"event-only view. Median K_spatial / K_event = "
          f"{frac.median():.2f}. "
          f"({n_no_coord} mothers had no usable event coords → K_spatial "
          f"defaults to K_event.)")

    n_miss = n_for_miss_probability()
    n_short_exp  = (~per_mother["budget_ok_exp"]).sum()
    n_short_miss = (~per_mother["budget_ok_miss"]).sum()
    print(f"[step28] miss-prob rule n = {n_miss} seeds "
          f"({n_short_miss}/{len(per_mother)} mothers below budget).")
    print(f"[step28] expected-coverage rule (90 %): "
          f"{n_short_exp}/{len(per_mother)} mothers below budget.")
    ac_exp  = per_mother["achieved_coverage_exp"]
    ac_miss = per_mother["achieved_coverage_miss"]
    print(f"[step28] achieved coverage under Rule 1: "
          f"mean = {ac_exp.mean():.2%}, "
          f"median = {ac_exp.median():.2%}, "
          f"min = {ac_exp.min():.2%}.")
    print(f"[step28] achieved coverage under Rule 2: "
          f"mean = {ac_miss.mean():.2%}, "
          f"median = {ac_miss.median():.2%}, "
          f"min = {ac_miss.min():.2%}.")

    rng = np.random.default_rng(2026)
    curves = build_curves_by_bin(rng)
    curves_path = tables_dir / "step28_coverage_curves_by_Nfertile.tsv"
    curves.to_csv(curves_path, sep="\t", index=False)
    print(f"[step28] Wrote {curves_path}")

    plot_coverage_curves(
        curves,
        out_png=figures_dir / "step28_coverage_curves.png",
        out_pdf=figures_dir / "step28_coverage_curves.pdf",
    )
    plot_per_mother_budget(
        per_mother,
        out_png=figures_dir / "step28_per_mother_budget.png",
        out_pdf=figures_dir / "step28_per_mother_budget.pdf",
    )
    print(f"[step28] Figures in {figures_dir}/")


if __name__ == "__main__":
    main()
