#!/usr/bin/env python3
"""Step 29a — Pollinator-radius sensitivity analysis.

Biological justification for the 50 m primary pollen-flight radius
adopted by every downstream Phase 5 script (§ 29b connectivity,
§ 29c fragmentation-aware sampling, § 30 SRK predictions, § 30b
fragmentation index).

Dependency order: run **Step 28 → Step 29 → Step 29a**.
This is a one-time validation — it does not need to be rerun on
every pipeline iteration once the 50 m primary radius is locked in.

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
    load_all_events, DEFAULT_DB,
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
def plot_sensitivity(summaries_by_year: dict,
                     bands: dict, out_png: Path, out_pdf: Path,
                     n_locations_by_year: dict | None = None):
    """Four-panel sensitivity view: connectivity, sampling, compatibility,
    sustainable fraction, all as a function of pollinator radius. 2025 and
    2026 are overlaid in each panel following the project convention —
    **2025 = open marker + dashed line**, **2026 = filled marker + solid
    line** (same convention as the merged Phase A prediction figures,
    Figures 8 and 10 of the compact doc)."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 10.5))
    (axA, axB), (axC, axD) = axes

    # Shared style per year — kept consistent across all panels.
    year_style = {
        2025: dict(marker="o", linestyle="--", markerfacecolor="white",
                   markeredgewidth=1.6, label_tag="2025 (open)"),
        2026: dict(marker="o", linestyle="-",  markerfacecolor=None,
                   markeredgewidth=0.6, label_tag="2026 (filled)"),
    }
    # Per-panel colour (same between years)
    PANEL_COL = {"A": "#1b7837", "B": "#c73030",
                 "C": "#5A5A9F", "D": "#009E73"}

    for yr, summary in summaries_by_year.items():
        if summary is None or summary.empty:
            continue
        st = year_style[yr]
        x = summary["radius_m"].values
        tag = st["label_tag"]
        n_loc = (n_locations_by_year or {}).get(yr, "?")

        mfc_A = ("white" if yr == 2025 else PANEL_COL["A"])
        mfc_B = ("white" if yr == 2025 else PANEL_COL["B"])
        mfc_C = ("white" if yr == 2025 else PANEL_COL["C"])
        mfc_D = ("white" if yr == 2025 else PANEL_COL["D"])

        axA.plot(x, summary["connected_share_median"],
                 marker=st["marker"], linestyle=st["linestyle"],
                 color=PANEL_COL["A"], markerfacecolor=mfc_A,
                 markeredgewidth=st["markeredgewidth"],
                 lw=1.8, label=f"{tag} · median ({n_loc} pops)")
        axB.plot(x, summary["M_frag_total"],
                 marker=st["marker"], linestyle=st["linestyle"],
                 color=PANEL_COL["B"], markerfacecolor=mfc_B,
                 markeredgewidth=st["markeredgewidth"],
                 lw=1.8, label=f"{tag} · total ({n_loc} pops)")
        axC.plot(x, summary["P_compat_median"],
                 marker=st["marker"], linestyle=st["linestyle"],
                 color=PANEL_COL["C"], markerfacecolor=mfc_C,
                 markeredgewidth=st["markeredgewidth"],
                 lw=1.8, label=f"{tag} · median")
        axD.plot(x, summary["frac_sustainable"],
                 marker=st["marker"], linestyle=st["linestyle"],
                 color=PANEL_COL["D"], markerfacecolor=mfc_D,
                 markeredgewidth=st["markeredgewidth"],
                 lw=1.8, label=f"{tag}")

    for ax in (axA, axB, axC, axD):
        ax.axvline(50, color="#666", ls=":", lw=1.0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axA.set_xlabel("Pollinator radius (m)")
    axA.set_ylabel("Fraction of adults in a multi-event component")
    axA.set_title("A. Landscape connectivity", fontsize=12, loc="left",
                  weight="bold")
    axA.set_ylim(0, 1.02)
    axA.legend(loc="lower right", fontsize=8)

    axB.set_xlabel("Pollinator radius (m)")
    axB.set_ylabel("Total mothers required (M_frag_aware)")
    axB.set_title("B. § B.4.2 sampling cost", fontsize=12, loc="left",
                  weight="bold")
    axB.legend(loc="upper right", fontsize=8)

    axC.axhline(bands["species_mean"], color="#1b7837", ls=":", lw=1.2,
                label=f"sporophytic species mean ({bands['species_mean']:.3f})")
    axC.axhline(bands["struggling_max"], color="#c73030", ls=":", lw=1.0,
                alpha=0.6,
                label=f"sustainable threshold ({bands['struggling_max']:.3f})")
    axC.set_xlabel("Pollinator radius (m)")
    axC.set_ylabel("Predicted per-mother P_compat mean")
    axC.set_title("C. Pollen-compatibility prediction",
                  fontsize=12, loc="left", weight="bold")
    axC.set_ylim(0, 1.0)
    axC.legend(loc="lower right", fontsize=8, frameon=True)

    axD.set_xlabel("Pollinator radius (m)")
    axD.set_ylabel("Fraction of locations in 'sustainable' band")
    axD.set_title("D. Locations at or above sustainable",
                  fontsize=12, loc="left", weight="bold")
    axD.set_ylim(0, 1.02)
    axD.legend(loc="lower right", fontsize=8)

    years_tag = " + ".join(str(y) for y in sorted(summaries_by_year.keys()))
    fig.suptitle(
        "Sensitivity of Phase 5 predictions + sampling to the "
        "pollinator radius\n"
        "(fragmentation, § B.4.2 M_frag, § A.8 pollen compatibility, "
        f"sustainable fraction — LEPA Snake River Plain, {years_tag} field)",
        fontsize=13, y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    fig.savefig(out_png, dpi=200); fig.savefig(out_pdf); plt.close(fig)


def _sweep_one_year(events: pd.DataFrame, m_lookup: dict,
                      prior_f: np.ndarray, class_i_mask: np.ndarray,
                      zygosity_probs: np.ndarray, bands: dict,
                      rng_seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the pollinator-radius sweep over `RADII_M` for one year's
    event set. Returns `(per_loc, summary)`."""
    rng = np.random.default_rng(rng_seed)
    rows_per_r = []
    for r in RADII_M:
        print(f"[sensitivity]   radius {r} m …")
        rows = compute_per_location(events, m_lookup, prior_f, class_i_mask,
                                     zygosity_probs, r, rng)
        rows_per_r.append(rows)
    per_loc = pd.concat(rows_per_r, ignore_index=True)
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
    return per_loc, summary


def main() -> None:
    if not PRIOR_TSV.exists():
        raise SystemExit(f"Missing {PRIOR_TSV}.")

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

    # Per-location M lookup (static — same for both years; falls back to
    # 5 if the step29 lookup is missing).
    m_lookup: dict[int, int] = {}
    if LOCATIONS_TSV.exists():
        loc_df = pd.read_csv(LOCATIONS_TSV, sep="\t", encoding="utf-8-sig")
        col = ("M_actual_in_step28" if "M_actual_in_step28" in loc_df.columns
               else "M_mothers_in_db" if "M_mothers_in_db" in loc_df.columns
               else None)
        if col is not None:
            m_lookup = dict(zip(loc_df["locationID"].astype(int),
                                loc_df[col].fillna(1).astype(int)))

    # locationCode lookup (for the per-location TSV)
    loc_meta = None
    if LOCATIONS_TSV.exists():
        loc_meta = (pd.read_csv(LOCATIONS_TSV, sep="\t", encoding="utf-8-sig")
                      [["locationID", "locationCode"]]
                      .drop_duplicates())

    # ---- Loop over 2025 and 2026, each year queried directly from the DB
    FIGURES.mkdir(parents=True, exist_ok=True)
    summaries_by_year: dict[int, pd.DataFrame] = {}
    n_locations_by_year: dict[int, int] = {}
    for yr in (2025, 2026):
        print(f"[sensitivity] =============================================="
              f"\n[sensitivity] Year {yr} — loading events from LEPA DB")
        events = load_all_events(DEFAULT_DB, year=yr)
        if events.empty:
            print(f"[sensitivity]   WARNING: no events for {yr}; skipping.")
            continue
        events = events.dropna(subset=["lat", "lon", "n_fertile"]).copy()
        events["n_fertile"] = events["n_fertile"].astype(int)
        n_loc = int(events["locationID"].nunique())
        n_evt = len(events)
        print(f"[sensitivity]   {n_evt} events across {n_loc} locationIDs")
        n_locations_by_year[yr] = n_loc

        # Fallback m_lookup for locationIDs missing from the static lookup
        local_m_lookup = dict(m_lookup)
        for lid in events["locationID"].unique():
            local_m_lookup.setdefault(int(lid), 5)

        per_loc, summary = _sweep_one_year(
            events, local_m_lookup, prior_f, class_i_mask, zygosity_probs,
            bands, rng_seed=2029 + yr)
        per_loc["year"] = yr
        summary["year"] = yr
        if loc_meta is not None:
            per_loc = per_loc.merge(loc_meta, on="locationID", how="left")

        out_ploc = TABLES / f"step30_A_radius_sensitivity_per_location_{yr}.tsv"
        per_loc.to_csv(out_ploc, sep="\t", index=False)
        out_sum  = TABLES / f"step30_A_radius_sensitivity_summary_{yr}.tsv"
        summary.to_csv(out_sum, sep="\t", index=False)
        print(f"[sensitivity]   wrote {out_ploc.name}, {out_sum.name}")
        print(summary[["radius_m", "connected_share_median", "M_frag_total",
                        "P_compat_median", "frac_sustainable"]]
              .round(3).to_string(index=False))
        summaries_by_year[yr] = summary

    # Combined summary TSV with year column (easy for downstream reporting)
    if summaries_by_year:
        combined = pd.concat(list(summaries_by_year.values()),
                              ignore_index=True)
        out_combined = TABLES / "step30_A_radius_sensitivity_summary.tsv"
        combined.to_csv(out_combined, sep="\t", index=False)
        print(f"[sensitivity] Combined summary → {out_combined.name}")

    # Overlay figure — 2025 open+dashed, 2026 filled+solid
    plot_sensitivity(
        summaries_by_year, bands,
        out_png=FIGURES / "step30_A_radius_sensitivity.png",
        out_pdf=FIGURES / "step30_A_radius_sensitivity.pdf",
        n_locations_by_year=n_locations_by_year,
    )
    print(f"[sensitivity] Figure in {FIGURES}/step30_A_radius_sensitivity.png")


if __name__ == "__main__":
    main()
