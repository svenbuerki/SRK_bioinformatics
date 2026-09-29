#!/usr/bin/env python3
"""Sensitivity of Phase 5 predictions + sampling to the pollinator radius.

Sweeps the primary pollen-flight radius across 10, 25, 50, 75, 100, 150,
200 m and reports at each radius:

  * **Landscape / connectivity** — per-location connected_share
    (fraction of adults in a multi-event pollen-flow component) and
    largest_component_share (fraction in the biggest component).
  * **§ B.4.2 fragmentation-aware sampling** — per-location
    M_frag_aware (mothers needed for 90 % chance to observe every
    predicted allele), plus total across the 39 locations.
  * **§ A.8 pollen compatibility** — per-location mean predicted
    P_compat under the sporophytic + empirical LEPA zygosity model.
  * **Traffic-light classification** — fraction of locations at
    P_compat mean ≥ sustainable threshold (2/3 of species mean).

Purpose
-------
Identify the radius at which the sampling cost stabilises AND the
fragmentation signal plateaus — that is the biologically sound
primary radius.

Outputs
-------
    tables/Phase5/step30_A_radius_sensitivity_per_location.tsv
    tables/Phase5/step30_A_radius_sensitivity_summary.tsv
    figures/Phase5/step30_A_radius_sensitivity.png/pdf
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from srk_bl_constants import BL_COLORS, BL_ORDER, locationCode_to_bl
from step28_seed_sampling_per_mother import (
    haversine_meters, K_pool, K_SPECIES_FG, PLOIDY,
    mothers_for_full_detection, n_for_miss_probability,
)
from srk_si_model import (
    load_class_map, build_class_i_mask,
    load_zygosity_dist, sample_genotypes_empirical,
    p_compat_sporophytic_empirical, species_mean_p_compat_empirical,
    traffic_light_bands,
)

RADII_M = [10, 25, 50, 75, 100, 150, 200]

TABLES = Path("Tables/Phase5")
FIGURES = Path("figures/Phase5")
EVENT_TSV = TABLES / "step28_events_spatial_neighborhood.tsv"
PRIOR_TSV = TABLES / "step26i_L1_carrier_inventory.tsv"
LOCATIONS_TSV = TABLES / "step29_sampling_per_location.tsv"

TARGET_PROB = 0.90            # § B.4.2 90 %-chance-to-see-every-allele target
N_DRAWS_MC = 120              # simulation replicates per location
N_FATHERS_MC = 200            # MC candidate fathers per replicate


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def connected_components(sub: pd.DataFrame, radius_m: float) -> list[list[int]]:
    """BFS on the radius-m adjacency graph over a location's events.
    Returns list-of-lists of POSITIONAL indices into `sub`."""
    n = len(sub)
    if n == 0:
        return []
    lat = sub["lat"].values
    lon = sub["lon"].values
    d = haversine_meters(lat[:, None], lon[:, None],
                         lat[None, :], lon[None, :])
    adj = (d > 0) & (d <= radius_m)
    seen = np.zeros(n, bool)
    comps: list[list[int]] = []
    for start in range(n):
        if seen[start]:
            continue
        stack = [start]
        comp: list[int] = []
        while stack:
            v = stack.pop()
            if seen[v]:
                continue
            seen[v] = True
            comp.append(v)
            for u in np.where(adj[v] & ~seen)[0]:
                stack.append(int(u))
        comps.append(comp)
    return comps


def location_pcompat_mean(N_fertile: int, M: int,
                           prior_f: np.ndarray, class_i_mask: np.ndarray,
                           zygosity_probs: np.ndarray,
                           rng: np.random.Generator) -> float:
    """Mean predicted per-mother P_compat across MC replicates for one
    location, using its radius-dependent effective N."""
    N_fertile = max(int(N_fertile), 1)
    M = max(int(M), 1)
    K_fg = len(prior_f)
    pool_size = PLOIDY * N_fertile
    replicate_means = []
    for _ in range(N_DRAWS_MC):
        local_alleles = rng.choice(K_fg, size=pool_size, p=prior_f)
        local_f = np.bincount(local_alleles, minlength=K_fg) / pool_size
        if not (local_f > 0).any():
            continue
        mothers = sample_genotypes_empirical(M, local_f, zygosity_probs, rng)
        pc = p_compat_sporophytic_empirical(
            mothers, local_f, class_i_mask, zygosity_probs,
            n_fathers=N_FATHERS_MC, rng=rng,
        )
        replicate_means.append(pc.mean())
    return float(np.mean(replicate_means)) if replicate_means else 0.0


def compute_per_location(events: pd.DataFrame, m_lookup: dict,
                          prior_f: np.ndarray, class_i_mask: np.ndarray,
                          zygosity_probs: np.ndarray,
                          radius_m: float,
                          rng: np.random.Generator) -> pd.DataFrame:
    """One row per location at the given radius."""
    rows = []
    for loc_id, sub in events.groupby("locationID"):
        sub = sub.reset_index(drop=True)
        comps = connected_components(sub, radius_m)
        n_events = len(sub)
        total_nf = int(sub["n_fertile"].astype(int).sum())
        comp_nf = np.array([int(sub.iloc[c]["n_fertile"].astype(int).sum())
                             for c in comps])
        comp_ne = np.array([len(c) for c in comps])
        largest_share = (float(comp_nf.max() / total_nf)
                         if total_nf > 0 else 0.0)
        connected_share = (float(comp_nf[comp_ne > 1].sum() / total_nf)
                           if total_nf > 0 else 0.0)
        F_location = 1.0 - connected_share
        n_fert_eff = int(round(total_nf * largest_share))

        # § B.4.2 fragmentation-aware M
        M_frag = 0
        for nf_c, ne_c in zip(comp_nf, comp_ne):
            K_c = min(K_pool(int(nf_c)), K_SPECIES_FG)
            M_pat = (1 if K_c <= 1
                     else mothers_for_full_detection(K_c, TARGET_PROB))
            M_c = max(M_pat, int(ne_c))
            M_frag += M_c

        # Predicted per-mother P_compat mean at this location + radius
        M_actual = int(m_lookup.get(int(loc_id), 1))
        pc_mean = location_pcompat_mean(
            n_fert_eff, M_actual,
            prior_f, class_i_mask, zygosity_probs, rng,
        )

        rows.append({
            "locationID":                 int(loc_id),
            "radius_m":                   int(radius_m),
            "n_events":                   n_events,
            "total_n_fertile":            total_nf,
            "n_components":               len(comps),
            "largest_component_share":    largest_share,
            "connected_share":            connected_share,
            "F_location":                 F_location,
            "N_fertile_effective":        n_fert_eff,
            "M_actual_in_step28":         M_actual,
            "M_frag_aware":               int(M_frag),
            "predicted_P_compat_mean":    pc_mean,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def plot_sensitivity(summary: pd.DataFrame, per_loc: pd.DataFrame,
                     bands: dict, out_png: Path, out_pdf: Path):
    """Four-panel sensitivity view: connectivity, sampling, compatibility,
    sustainable fraction, all as a function of pollinator radius."""
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 10.0))
    (axA, axB), (axC, axD) = axes

    x = summary["radius_m"].values

    # A. Connectivity — mean and median connected_share across locations.
    axA.plot(x, summary["connected_share_median"],
             "o-", color="#1b7837", lw=2.0, label="median across locations")
    axA.plot(x, summary["connected_share_mean"],
             "s--", color="#1b7837", lw=1.4, alpha=0.7,
             label="mean across locations")
    axA.axvline(50, color="#666", ls=":", lw=1.0)
    axA.text(50, 0.02, "  adopted primary\n  = 50 m",
             fontsize=8, color="#666", va="bottom")
    axA.set_xlabel("Pollinator radius (m)"); axA.set_ylabel(
        "Fraction of adults in a multi-event component")
    axA.set_title("A. Landscape connectivity", fontsize=12, loc="left",
                  weight="bold")
    axA.set_ylim(0, 1.02)
    axA.legend(loc="lower right", fontsize=9)
    axA.spines["top"].set_visible(False); axA.spines["right"].set_visible(False)

    # B. Fragmentation-aware sampling — total mothers across the 39 locations.
    axB.plot(x, summary["M_frag_total"], "o-", color="#c73030", lw=2.0,
             label="total M_frag_aware across 39 locations")
    axB.axvline(50, color="#666", ls=":", lw=1.0)
    axB.text(25, summary["M_frag_total"].max() * 0.02,
             "  50 m", fontsize=8, color="#666", va="bottom")
    axB.set_xlabel("Pollinator radius (m)")
    axB.set_ylabel("Total mothers required (M_frag_aware)")
    axB.set_title("B. § B.4.2 sampling cost", fontsize=12, loc="left",
                  weight="bold")
    axB.legend(loc="upper right", fontsize=9)
    axB.spines["top"].set_visible(False); axB.spines["right"].set_visible(False)

    # C. Species-mean pollen compatibility (over locations) at each radius.
    axC.plot(x, summary["P_compat_median"], "o-", color="#5A5A9F", lw=2.0,
             label="median location P_compat mean")
    axC.plot(x, summary["P_compat_mean"], "s--", color="#5A5A9F", lw=1.4,
             alpha=0.7, label="mean location P_compat mean")
    axC.axhline(bands["species_mean"], color="#1b7837", ls=":", lw=1.2,
                label=f"sporophytic species mean ({bands['species_mean']:.3f})")
    axC.axhline(bands["struggling_max"], color="#c73030", ls=":", lw=1.0,
                alpha=0.6, label=f"sustainable threshold "
                                  f"({bands['struggling_max']:.3f})")
    axC.axvline(50, color="#666", ls=":", lw=1.0)
    axC.set_xlabel("Pollinator radius (m)")
    axC.set_ylabel("Predicted per-mother P_compat mean")
    axC.set_title("C. Pollen-compatibility prediction",
                  fontsize=12, loc="left", weight="bold")
    axC.set_ylim(0, 1.0)
    axC.legend(loc="lower right", fontsize=8, frameon=True)
    axC.spines["top"].set_visible(False); axC.spines["right"].set_visible(False)

    # D. Fraction of locations in sustainable band.
    axD.plot(x, summary["frac_sustainable"], "o-", color="#009E73", lw=2.0)
    axD.axvline(50, color="#666", ls=":", lw=1.0)
    axD.set_xlabel("Pollinator radius (m)")
    axD.set_ylabel("Fraction of locations in 'sustainable' band")
    axD.set_title("D. Locations at or above sustainable",
                  fontsize=12, loc="left", weight="bold")
    axD.set_ylim(0, 1.02)
    axD.spines["top"].set_visible(False); axD.spines["right"].set_visible(False)

    fig.suptitle(
        "Sensitivity of Phase 5 predictions + sampling to the "
        "pollinator radius\n"
        "(fragmentation, § B.4.2 M_frag, § A.8 pollen compatibility, "
        "sustainable fraction — 39 LEPA locations, 2025 field)",
        fontsize=13, y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    fig.savefig(out_png, dpi=200); fig.savefig(out_pdf); plt.close(fig)


def main() -> None:
    if not EVENT_TSV.exists():
        raise SystemExit(f"Missing {EVENT_TSV} — run step28 first.")
    if not PRIOR_TSV.exists():
        raise SystemExit(f"Missing {PRIOR_TSV}.")

    events = pd.read_csv(EVENT_TSV, sep="\t", encoding="utf-8-sig")
    # Only the year-filtered events (step28 writes only 2025 when --year set)
    events = events.dropna(subset=["lat", "lon", "n_fertile"])
    events["n_fertile"] = events["n_fertile"].astype(int)

    # P1 prior
    p1 = pd.read_csv(PRIOR_TSV, sep="\t", encoding="utf-8-sig")
    fg_freq = (p1.groupby("Fg")["n_carriers"].sum().reset_index()
                 .sort_values("Fg").reset_index(drop=True))
    prior_f = (fg_freq["n_carriers"] / fg_freq["n_carriers"].sum()).values
    fg_labels = fg_freq["Fg"].astype(str).tolist()
    class_i_mask = build_class_i_mask(fg_labels, load_class_map())
    zygosity_probs = load_zygosity_dist()

    # Species-mean bands (independent of radius — used only as reference)
    rng0 = np.random.default_rng(2026)
    species_mean_pc = species_mean_p_compat_empirical(
        prior_f, class_i_mask, zygosity_probs,
        n_mothers=5_000, n_fathers=1_000, rng=rng0)
    bands = traffic_light_bands(species_mean_pc)
    print(f"[sensitivity] Sporophytic + empirical-zygosity species-mean = "
          f"{species_mean_pc:.4f}; bands failed < {bands['failed_max']:.3f}, "
          f"struggling < {bands['struggling_max']:.3f}, sustainable ≥ "
          f"{bands['struggling_max']:.3f}.")

    # Per-location M lookup from step29
    m_lookup = {}
    if LOCATIONS_TSV.exists():
        loc_df = pd.read_csv(LOCATIONS_TSV, sep="\t", encoding="utf-8-sig")
        col = ("M_actual_in_step28" if "M_actual_in_step28" in loc_df.columns
               else "M_mothers_in_db" if "M_mothers_in_db" in loc_df.columns
               else None)
        if col is not None:
            m_lookup = dict(zip(loc_df["locationID"].astype(int),
                                loc_df[col].fillna(1).astype(int)))
    if not m_lookup:
        # Fallback: use 5 mothers per location as a placeholder
        m_lookup = {int(k): 5 for k in events["locationID"].unique()}

    all_rows = []
    rng = np.random.default_rng(2029)
    for r in RADII_M:
        print(f"[sensitivity] radius {r} m …")
        rows = compute_per_location(events, m_lookup, prior_f, class_i_mask,
                                     zygosity_probs, r, rng)
        all_rows.append(rows)
    per_loc = pd.concat(all_rows, ignore_index=True)

    # Attach locationCode from step29 for readability
    if LOCATIONS_TSV.exists():
        loc_meta = pd.read_csv(LOCATIONS_TSV, sep="\t", encoding="utf-8-sig")
        per_loc = per_loc.merge(
            loc_meta[["locationID", "locationCode"]].drop_duplicates(),
            on="locationID", how="left")

    out_ploc = TABLES / "step30_A_radius_sensitivity_per_location.tsv"
    per_loc.to_csv(out_ploc, sep="\t", index=False)
    print(f"[sensitivity] Wrote {out_ploc}")

    # Summary — one row per radius across the 39 locations
    summary = (per_loc.groupby("radius_m")
               .agg(n_locations=("locationID", "nunique"),
                    connected_share_mean=("connected_share", "mean"),
                    connected_share_median=("connected_share", "median"),
                    F_location_mean=("F_location", "mean"),
                    n_components_median=("n_components", "median"),
                    N_fertile_effective_median=("N_fertile_effective", "median"),
                    M_frag_total=("M_frag_aware", "sum"),
                    M_frag_median=("M_frag_aware", "median"),
                    P_compat_mean=("predicted_P_compat_mean", "mean"),
                    P_compat_median=("predicted_P_compat_mean", "median"),
                    frac_sustainable=(
                        "predicted_P_compat_mean",
                        lambda s: float((s >= bands["struggling_max"]).mean())),
                    frac_struggling=(
                        "predicted_P_compat_mean",
                        lambda s: float(((s >= bands["failed_max"])
                                          & (s < bands["struggling_max"])).mean())),
                    frac_failed=(
                        "predicted_P_compat_mean",
                        lambda s: float((s < bands["failed_max"]).mean())),
                    )
               .reset_index())
    summary["species_mean_reference"] = bands["species_mean"]
    summary["sustainable_threshold"] = bands["struggling_max"]
    summary["failed_threshold"] = bands["failed_max"]

    out_sum = TABLES / "step30_A_radius_sensitivity_summary.tsv"
    summary.to_csv(out_sum, sep="\t", index=False)
    print(f"[sensitivity] Wrote {out_sum}")

    print("[sensitivity] Summary preview:")
    print(summary[["radius_m", "connected_share_median", "M_frag_total",
                   "P_compat_median", "frac_sustainable"]].round(3).to_string(index=False))

    FIGURES.mkdir(parents=True, exist_ok=True)
    plot_sensitivity(
        summary, per_loc, bands,
        out_png=FIGURES / "step30_A_radius_sensitivity.png",
        out_pdf=FIGURES / "step30_A_radius_sensitivity.pdf",
    )
    print(f"[sensitivity] Figure in {FIGURES}/")


if __name__ == "__main__":
    main()
