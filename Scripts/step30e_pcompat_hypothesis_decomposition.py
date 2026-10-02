"""Step 30 Part C — SRK diversity collapse → pollen-compatibility decomposition.

This diagnostic is **triggered by the per-location SRK diversity
discrepancy** surfaced in Panel A of the Part C anchor figure
(`step30_B_partC_clean_overlap.png`). When a location's observed
distinct-allele count sits far below the Phase 5 prediction, the
local pool has drifted and the next question is *how that drift
redistributes into the pollen compatibility prediction*. Two
locations with the same number of lost alleles can land at very
different observed pollen compatibilities depending on **which**
alleles were lost. This script quantifies that.

For each clean-overlap location (EO67, EO70, EO76) the observed
pollen compatibility can deviate from Phase 5's P1-based prediction
through three distinct biological mechanisms:

  1. **Class I / II mass imbalance** — the Class I share of the
     pool differs from species-wide. Because between-class crosses
     are always compatible, raising / lowering the Class I mass
     changes how much "between-class rescue" is available to Class II
     mothers.
  2. **Within-class allele spread (concentration)** — the TOTAL
     Class I and Class II masses can match species-wide, but all of
     Class I's mass may be stuck in a single allele (drift
     monoculture) rather than spread across several. Within-class
     monoculture kills within-class compatibility (a Class I × Class I
     cross is rejected whenever parents share an allele).
  3. **Zygosity composition** — the fraction of mothers that carry
     1 / 2 / 3 distinct SRK identities. Multi-identity mothers
     express a bigger set and face more compatible fathers.

This script decomposes each location's observed pollen compatibility
by swapping ONE factor at a time against the species-wide P1
reference and recomputing pollen compatibility against the SAME
observed mother genotypes:

    observed_full     — observed f_local + observed zygosity (truth)
    swap_within_class — class balance kept at observed, within-class
                        spread reshaped to the P1 shape
    swap_class_bal    — within-class spread kept at observed, class
                        balance (p_I) reshaped to the P1 value
    swap_zygosity     — f_local kept at observed, father zygosity
                        reshaped to species-wide
    pred_P1_full      — species-wide P1 f + species-wide zygosity
                        (= Phase 5 prediction with observed mothers)

If `swap_within_class` moves pollen compatibility close to
`pred_P1_full`, within-class concentration is the dominant driver.
If `swap_class_bal` does, class balance is the driver. If
`swap_zygosity` does, zygosity dominates. Several factors can
contribute additively.

Example (worked out in the compact doc §): EO67 observed = 0.70,
EO70 observed = 0.53, both with identical Class I mass (~0.37) and
similar zygosity (~0.70 one-identity). The decomposition confirms
that EO70's deficit is driven almost entirely by within-class
concentration — one Class I allele (FG024, 35 %) absorbs almost all
of Class I's mass, killing Class I × Class I compatibility.

Outputs
-------
Tables/Phase5/step30_B_partC_hypothesis_decomposition.tsv
    One row per clean-overlap location with observed, three single-
    swap counterfactuals, and the full P1 prediction, each with
    bootstrap 95 % CIs. The columns `driver_within_class`,
    `driver_class_balance`, `driver_zygosity` give the fraction of
    the observed→P1 gap each factor explains in isolation.
figures/Phase5/step30_B_partC_hypothesis_decomposition.png / .pdf
    Grouped bar chart per location.

CLI
---
    python step30e_pcompat_hypothesis_decomposition.py
        Runs on the three clean-overlap EOs by default.
    python step30e_pcompat_hypothesis_decomposition.py --eos EO67,EO76
        Restrict to a subset.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from srk_si_model import (
    load_class_map,
    build_class_i_mask,
    load_zygosity_dist,
    p_compat_sporophytic_empirical,
)
from step30c_srk_validation_at_eo_level import (
    load_allele_to_fg,
    load_individuals,
    compute_local_f,
    bootstrap_pcompat_mean,
)

DEFAULT_TABLES  = Path("Tables/Phase5")
DEFAULT_FIGURES = Path("figures/Phase5")
DEFAULT_PRIOR_TSV = DEFAULT_TABLES / "step30_A_prediction_prior_frequencies.tsv"
DEFAULT_ALLELE_FG_TSV = DEFAULT_TABLES / "step26i_L1_carrier_inventory.tsv"
DEFAULT_INDIV_TSV = Path(
    "Tables/Phase4/step23_individual_allele_genotypes_with_nulls.tsv")
DEFAULT_BANDS_TSV = DEFAULT_TABLES / "step30_A_traffic_light_bands.tsv"
# Link to the Phase 5 diversity prediction so the trigger metric (how
# many Fgs did the location *lose*?) can be reported alongside each
# decomposition.
DEFAULT_PRED_DIVERSITY = DEFAULT_TABLES / "step30_A_prediction_location_diversity.tsv"

CLEAN_OVERLAP_EOS = ["EO67", "EO70", "EO76"]

N_MC_FATHERS = 800
N_BOOTSTRAP  = 400
RNG_SEED     = 2034


# ---------------------------------------------------------------------------
# Counterfactual f_local builders
# ---------------------------------------------------------------------------
def swap_within_class_spread(f_obs: np.ndarray,
                              f_P1: np.ndarray,
                              class_i_mask: np.ndarray) -> np.ndarray:
    """Keep the observed Class I / Class II TOTAL masses, but reshape
    the within-class spread to match the P1 shape.

    Fgs that are present in P1 but absent in the observed pool get
    their P1 share within the class; nothing fancier.
    """
    f = np.zeros_like(f_obs)
    for mask in (class_i_mask, ~class_i_mask):
        obs_mass = float(f_obs[mask].sum())
        p1_mass = float(f_P1[mask].sum())
        if p1_mass <= 0.0:
            continue
        # P1 shape, scaled to the observed mass of this class
        f[mask] = f_P1[mask] * (obs_mass / p1_mass)
    s = f.sum()
    return f / s if s > 0 else f


def swap_class_balance(f_obs: np.ndarray,
                       f_P1: np.ndarray,
                       class_i_mask: np.ndarray) -> np.ndarray:
    """Keep the observed WITHIN-CLASS shape, but scale the two class
    totals so Class I mass equals the P1 Class I mass."""
    f = np.zeros_like(f_obs)
    target_mass = {True: float(f_P1[class_i_mask].sum()),
                   False: float(f_P1[~class_i_mask].sum())}
    for cls_val, mask in ((True, class_i_mask), (False, ~class_i_mask)):
        obs_mass = float(f_obs[mask].sum())
        if obs_mass <= 0.0:
            # No observed alleles of this class — fall back to P1
            if f_P1[mask].sum() > 0:
                f[mask] = f_P1[mask] / f_P1[mask].sum() * target_mass[cls_val]
            continue
        f[mask] = f_obs[mask] / obs_mass * target_mass[cls_val]
    s = f.sum()
    return f / s if s > 0 else f


# ---------------------------------------------------------------------------
# Decomposition per EO
# ---------------------------------------------------------------------------
def per_eo_decomposition(indiv_df: pd.DataFrame,
                          prior: pd.DataFrame,
                          class_i_mask: np.ndarray,
                          zygosity_species: np.ndarray,
                          rng: np.random.Generator,
                          eos: List[str],
                          pred_diversity: pd.DataFrame | None = None) -> pd.DataFrame:
    f_P1 = prior["f_mean"].values
    pred_div_by_loc = (pred_diversity.set_index("locationCode")
                        if pred_diversity is not None
                        else None)
    rows = []
    for eo in eos:
        sub_all = indiv_df[indiv_df["EO"] == eo]
        sub = sub_all[sub_all["n_functional_copies"] > 0]
        n_fun = len(sub)
        if n_fun == 0:
            continue
        genotypes = np.stack(sub["genotype"].values, axis=0)
        fg_copies_stack = np.stack(sub["fg_copies"].values, axis=0)
        f_obs = compute_local_f(fg_copies_stack)

        # Observed zygosity from the data (indices 0..4 → 0..4 distinct)
        zyg_counts = np.bincount(
            sub["n_distinct_functional"].values, minlength=5)[:5]
        zyg_obs = zyg_counts / max(zyg_counts.sum(), 1)
        # Match the length used by the model (typically length 4 or 5)
        zyg_obs = zyg_obs[:len(zygosity_species)]
        if zyg_obs.sum() > 0:
            zyg_obs = zyg_obs / zyg_obs.sum()

        # Counterfactual local-frequency vectors
        f_swap_within = swap_within_class_spread(f_obs, f_P1, class_i_mask)
        f_swap_class  = swap_class_balance(f_obs, f_P1, class_i_mask)

        # Bootstrap pollen compatibility under each scenario, same
        # mother genotypes throughout — so the comparison isolates the
        # father-drawing side of the model.
        scenarios = {
            "observed_full":     (f_obs, zyg_obs),
            "swap_within_class": (f_swap_within, zyg_obs),
            "swap_class_bal":    (f_swap_class, zyg_obs),
            "swap_zygosity":     (f_obs, zygosity_species),
            "pred_P1_full":      (f_P1, zygosity_species),
        }
        # Diversity-trigger metrics — this diagnostic is motivated
        # by the observed vs Phase 5 SRK diversity gap at this
        # location. Record it alongside each pollen compatibility
        # row so the chain is explicit.
        n_distinct_obs = int((f_obs > 0).sum())
        if pred_div_by_loc is not None and eo in pred_div_by_loc.index:
            div_row = pred_div_by_loc.loc[eo]
            pred_div = float(div_row["predicted_local_pool_size_mean"])
            pred_div_lo = float(div_row["predicted_local_pool_size_lo"])
            pred_div_hi = float(div_row["predicted_local_pool_size_hi"])
        else:
            pred_div = np.nan
            pred_div_lo = np.nan
            pred_div_hi = np.nan
        row = {
            "locationCode":          eo,
            "n_functional_carriers": n_fun,
            # SRK diversity trigger — the gap that motivated the
            # decomposition in the first place.
            "obs_distinct_SRK_alleles":           n_distinct_obs,
            "pred_phase5_distinct_SRK_alleles":   pred_div,
            "pred_phase5_distinct_SRK_alleles_lo": pred_div_lo,
            "pred_phase5_distinct_SRK_alleles_hi": pred_div_hi,
            "diversity_gap_obs_minus_pred":       (n_distinct_obs - pred_div
                                                    if not np.isnan(pred_div)
                                                    else np.nan),
            # Local-pool composition metrics
            "obs_p_ClassI":          float(f_obs[class_i_mask].sum()),
            "p1_p_ClassI":           float(f_P1[class_i_mask].sum()),
            "obs_zyg_1_distinct":    float(zyg_obs[0]) if len(zyg_obs) > 0 else np.nan,
            "obs_zyg_2_distinct":    float(zyg_obs[1]) if len(zyg_obs) > 1 else np.nan,
            "obs_zyg_3plus_distinct":float(zyg_obs[2:].sum()) if len(zyg_obs) > 2 else np.nan,
        }
        for scen_name, (f_use, zyg_use) in scenarios.items():
            m, lo, hi = bootstrap_pcompat_mean(
                genotypes, f_use, class_i_mask, zyg_use,
                N_MC_FATHERS, N_BOOTSTRAP, rng)
            row[f"pcompat_{scen_name}_mean"] = m
            row[f"pcompat_{scen_name}_lo95"] = lo
            row[f"pcompat_{scen_name}_hi95"] = hi

        # How much of the observed → P1 gap does each single-swap
        # explain? Positive = that factor moves pollen compatibility
        # toward the P1 prediction (so it was part of the deficit).
        total_gap = row["pcompat_pred_P1_full_mean"] - row["pcompat_observed_full_mean"]
        for factor in ("swap_within_class", "swap_class_bal", "swap_zygosity"):
            swap_m = row[f"pcompat_{factor}_mean"]
            if abs(total_gap) < 1e-6:
                row[f"driver_{factor}_fraction"] = 0.0
            else:
                row[f"driver_{factor}_fraction"] = float(
                    (swap_m - row["pcompat_observed_full_mean"]) / total_gap)
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def draw_figure(stats: pd.DataFrame,
                bands: dict,
                out_png: Path, out_pdf: Path) -> None:
    # Order: baseline (prediction) → ground truth (observed) → three
    # single-swap scenarios that explain the gap between them.
    scen_labels = [
        ("pred_P1_full",      "Predicted (species-wide)",                "#3c8f4c"),
        ("observed_full",     "Observed",                                "#8a8a8a"),
        ("swap_within_class", "Swap within-class spread → species-wide", "#2b6cb0"),
        ("swap_class_bal",    "Swap Class I/II balance → species-wide",  "#c94b4b"),
        ("swap_zygosity",     "Swap zygosity → species-wide",            "#e6b325"),
    ]
    n_loc = len(stats)
    n_scen = len(scen_labels)

    fig, ax = plt.subplots(figsize=(12.0, 5.5))

    # Traffic-light background bands (horizontal strips, pollen
    # compatibility is on the y-axis)
    ax.axhspan(0, bands["failed_max"],
                color="#c94b4b", alpha=0.10, zorder=0)
    ax.axhspan(bands["failed_max"], bands["struggling_max"],
                color="#e6b325", alpha=0.10, zorder=0)
    ax.axhspan(bands["struggling_max"], 1.0,
                color="#3c8f4c", alpha=0.10, zorder=0)
    ax.axhline(bands["species_mean"], color="#555555",
                lw=0.8, ls="--", alpha=0.5)

    width = 1.0
    # Horizontal offset for each scenario WITHIN a group.  Visual
    # separator between the two reference bars (Predicted + Observed)
    # and the three single-swap explanations lives in `GAP`.
    GAP = 0.9
    scen_offsets = [0.0, 1.0, 2.0 + GAP, 3.0 + GAP, 4.0 + GAP]
    # Group spacing: each group takes (max offset + width + inter-
    # group padding).
    group_pitch = max(scen_offsets) + width + 2.0
    group_centres = np.arange(n_loc) * group_pitch
    local_centre = (max(scen_offsets) + 0.0) / 2.0

    # Pre-compute, per location, which single-swap scenario moves
    # pollen compatibility furthest from Observed — the dominant
    # explanation of the gap between Predicted and Observed. Ties
    # break by the first-listed scenario.
    swap_scen_names = ["swap_within_class", "swap_class_bal",
                        "swap_zygosity"]
    best_swap_per_loc: list[str] = []
    for _, r in stats.iterrows():
        obs = r["pcompat_observed_full_mean"]
        best_scen = max(swap_scen_names,
                        key=lambda s: abs(r[f"pcompat_{s}_mean"] - obs))
        best_swap_per_loc.append(best_scen)

    group_tops = np.full(n_loc, -np.inf)

    for s_idx, (scen, label, colour) in enumerate(scen_labels):
        xs = group_centres + scen_offsets[s_idx] - local_centre
        means = stats[f"pcompat_{scen}_mean"].values
        lows  = stats[f"pcompat_{scen}_lo95"].values
        highs = stats[f"pcompat_{scen}_hi95"].values
        yerr = np.array([means - lows, highs - means])
        yerr = np.clip(yerr, 0.0, None)
        ax.bar(xs, means, width=width, color=colour,
               edgecolor="white", linewidth=0.5, label=label, zorder=2)
        ax.errorbar(xs, means, yerr=yerr, fmt="none",
                     ecolor="#333333", elinewidth=0.8, capsize=2,
                     alpha=0.8, zorder=3)
        group_tops = np.maximum(group_tops, highs)

    # Stars on the dominant explanation
    for loc_idx, best_scen in enumerate(best_swap_per_loc):
        s_idx = [name for name, _, _ in scen_labels].index(best_scen)
        x_star = group_centres[loc_idx] + scen_offsets[s_idx] - local_centre
        y_star = group_tops[loc_idx] + 0.03
        ax.scatter(x_star, y_star, marker="*", s=200,
                    color="#f6b800", edgecolor="#333333",
                    linewidth=0.8, zorder=5)
    # Legend entry for the star — reused across groups.
    star_handle = plt.Line2D(
        [], [], linestyle="None", marker="*", markersize=12,
        markerfacecolor="#f6b800", markeredgecolor="#333333",
        label="Dominant explanation of the gap",
    )

    ax.set_xticks(group_centres)
    ax.set_xticklabels(
        [f"{code}\n(n={int(n)})" for code, n in zip(
            stats["locationCode"], stats["n_functional_carriers"])],
        fontsize=10,
    )
    ax.set_ylim(0.0, max(0.9, 1.05 * bands["struggling_max"] + 0.3))
    ax.set_ylabel("Pollen compatibility under random mating", fontsize=11)
    ax.set_title(
        "Pollen compatibility — competing-hypothesis decomposition "
        "per location",
        fontsize=12, loc="left",
    )
    handles, labels = ax.get_legend_handles_labels()
    handles.append(star_handle)
    labels.append(star_handle.get_label())
    ax.legend(handles=handles, labels=labels,
              loc="lower center", bbox_to_anchor=(0.5, -0.30),
              ncol=3, fontsize=9, frameon=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.subplots_adjust(bottom=0.30, left=0.08, right=0.98, top=0.90)
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eos", type=str, default=",".join(CLEAN_OVERLAP_EOS))
    args = ap.parse_args()
    eos = [e.strip() for e in args.eos.split(",") if e.strip()]

    rng = np.random.default_rng(RNG_SEED)
    prior = pd.read_csv(DEFAULT_PRIOR_TSV, sep="\t", encoding="utf-8-sig")
    fg_labels = prior["Fg"].astype(str).tolist()
    allele_to_fg = load_allele_to_fg(DEFAULT_ALLELE_FG_TSV)
    indiv_df = load_individuals(DEFAULT_INDIV_TSV, allele_to_fg, fg_labels)
    class_map = load_class_map()
    class_i_mask = build_class_i_mask(fg_labels, class_map)
    zyg_species = load_zygosity_dist()
    bands_df = pd.read_csv(DEFAULT_BANDS_TSV, sep="\t", encoding="utf-8-sig")
    bands = {
        "species_mean":   float(bands_df["species_mean"].iloc[0]),
        "failed_max":     float(bands_df["failed_max"].iloc[0]),
        "struggling_max": float(bands_df["struggling_max"].iloc[0]),
    }
    pred_div = pd.read_csv(DEFAULT_PRED_DIVERSITY, sep="\t",
                             encoding="utf-8-sig")

    stats = per_eo_decomposition(indiv_df, prior, class_i_mask,
                                   zyg_species, rng, eos,
                                   pred_diversity=pred_div)
    DEFAULT_TABLES.mkdir(parents=True, exist_ok=True)
    DEFAULT_FIGURES.mkdir(parents=True, exist_ok=True)
    tsv_path = DEFAULT_TABLES / "step30_B_partC_hypothesis_decomposition.tsv"
    stats.to_csv(tsv_path, sep="\t", index=False)
    print(f"[step30e] Wrote {tsv_path} ({len(stats)} rows)")

    fig_png = DEFAULT_FIGURES / "step30_B_partC_hypothesis_decomposition.png"
    fig_pdf = DEFAULT_FIGURES / "step30_B_partC_hypothesis_decomposition.pdf"
    draw_figure(stats, bands, fig_png, fig_pdf)
    print(f"[step30e] Wrote {fig_png} + .pdf")

    # Console summary — lead with the diversity trigger, then decompose.
    print()
    print("[step30e] ===== Diversity-trigger diagnostic: SRK diversity → P_compat =====")
    for _, r in stats.iterrows():
        gap = r["pcompat_pred_P1_full_mean"] - r["pcompat_observed_full_mean"]
        print(f"[step30e] {r['locationCode']}  (n={int(r['n_functional_carriers'])} adults)")
        print(f"[step30e]   TRIGGER — SRK diversity:")
        print(f"[step30e]      observed = {int(r['obs_distinct_SRK_alleles'])}/32 alleles  |  "
              f"Phase 5 predicted = {r['pred_phase5_distinct_SRK_alleles']:.1f}  "
              f"[{r['pred_phase5_distinct_SRK_alleles_lo']:.1f}, "
              f"{r['pred_phase5_distinct_SRK_alleles_hi']:.1f}]  |  "
              f"gap = {r['diversity_gap_obs_minus_pred']:+.1f}")
        print(f"[step30e]   DOWNSTREAM — pollen compatibility:")
        print(f"[step30e]      observed = {r['pcompat_observed_full_mean']:.3f}  |  "
              f"Phase 5 pred (P1) = {r['pcompat_pred_P1_full_mean']:.3f}  |  "
              f"gap = {gap:+.3f}")
        print(f"[step30e]   DECOMPOSITION — which factor explains the gap?")
        print(f"[step30e]      within-class allele spread   {r['driver_swap_within_class_fraction'] * 100:+6.1f}%   "
              f"(swap →  {r['pcompat_swap_within_class_mean']:.3f})")
        print(f"[step30e]      Class I/II balance           {r['driver_swap_class_bal_fraction'] * 100:+6.1f}%   "
              f"(swap →  {r['pcompat_swap_class_bal_mean']:.3f})")
        print(f"[step30e]      zygosity composition         {r['driver_swap_zygosity_fraction'] * 100:+6.1f}%   "
              f"(swap →  {r['pcompat_swap_zygosity_mean']:.3f})")


if __name__ == "__main__":
    main()
