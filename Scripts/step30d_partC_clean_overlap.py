"""Step 30 Part C — Clean-overlap EOs → Phase 5 location-scale comparison.

Three EOs sit 1:1 with a Phase 5 locationCode AND carry n >= 10
functional-carrier adult genotypes in Phase 4 Step 23 (the clean
overlap set): EO76, EO70, EO67. For these three the Phase 4 adult SRK
data can feed Phase 5 Part C *directly*, at the Phase 5 location
scale, with no event-level remap required.

What this script does, per clean-overlap location:
  - Observed Fg set + frequencies from the adult genotypes
  - Observed per-mother P_compat under sporophytic Class I / II +
    empirical zygosity, with FATHERS drawn from OBSERVED local Fg
    frequencies
  - Observed SRK diversity = |observed Fg set|
  - Compared to the Phase 5 per-location prediction tables
    (`step30_A_prediction_location_pcompat.tsv` +
    `step30_A_prediction_location_diversity.tsv`)

Outputs (all in `Tables/Phase5/` with the Phase-B naming convention,
so they live alongside demo / future real seed-based outputs):
  step30_B_partC_clean_overlap_per_location.tsv
      One row per clean location with observed vs predicted P_compat
      and SRK diversity.
  step30_B_partC_clean_overlap_fg_frequencies.tsv
      Long form (location x Fg) with observed f, species-wide P1 f,
      and whether the Fg is present at the location.

And a figure in `figures/Phase5/`:
  step30_B_partC_clean_overlap.png / .pdf
      Three panels per location: observed vs predicted P_compat (dot
      + 95 % CI), observed vs predicted distinct Fg count, observed
      vs P1 Fg frequency spectrum.

Caveats
-------
  * Three locations only — this is an anchor for Part C, not the full
    Phase B regression; the mate-limitation regression needs seed set
    + per-mother P_compat and will come online once seed genotypes
    arrive at the other clean-overlap locations.
  * P1 was built from these same individuals, so the species-mean
    comparison is overlapping — the honest test is per-location
    frequency skew (same caveat as step30c).
  * The three split EOs (EO18, EO25, EO27) are deliberately excluded
    here; they need an Individual -> germplasmID -> eventID remap
    before their adults can be used at the Phase 5 location scale.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

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

# step30c already has the per-individual loader and the local-f /
# bootstrap helpers; reuse them to keep the model identical.
from step30c_srk_validation_at_eo_level import (
    load_allele_to_fg,
    load_individuals,
    compute_local_f,
    bootstrap_pcompat_mean,
)

DEFAULT_TABLES  = Path("Tables/Phase5")
DEFAULT_FIGURES = Path("figures/Phase5")
DEFAULT_PRIOR_TSV   = DEFAULT_TABLES / "step30_A_prediction_prior_frequencies.tsv"
DEFAULT_ALLELE_FG_TSV = DEFAULT_TABLES / "step26i_L1_carrier_inventory.tsv"
DEFAULT_INDIV_TSV   = Path(
    "Tables/Phase4/step23_individual_allele_genotypes_with_nulls.tsv")
DEFAULT_PRED_PCOMPAT = DEFAULT_TABLES / "step30_A_prediction_location_pcompat.tsv"
DEFAULT_PRED_DIVERSITY = DEFAULT_TABLES / "step30_A_prediction_location_diversity.tsv"
DEFAULT_BANDS_TSV   = DEFAULT_TABLES / "step30_A_traffic_light_bands.tsv"

# EO <-> Phase 5 locationCode is 1:1 for these three.
CLEAN_OVERLAP_EOS = ["EO67", "EO70", "EO76"]

N_MC_FATHERS = 800
N_BOOTSTRAP  = 400
RNG_SEED     = 2033

COLOUR_OBS   = "#333333"
COLOUR_PRED  = "#2b6cb0"
COLOUR_FAILED     = "#c94b4b"
COLOUR_STRUGGLING = "#e6b325"
COLOUR_SUSTAIN    = "#3c8f4c"


def per_location_stats(indiv_df: pd.DataFrame,
                       prior: pd.DataFrame,
                       class_i_mask: np.ndarray,
                       zygosity_probs: np.ndarray,
                       pred_pcompat: pd.DataFrame,
                       pred_diversity: pd.DataFrame,
                       rng: np.random.Generator,
                       eos: List[str]) -> pd.DataFrame:
    """One row per clean-overlap EO (= Phase 5 locationCode)."""
    prior_f = prior["f_mean"].values
    fg_labels = prior["Fg"].astype(str).tolist()
    rows = []
    for eo in eos:
        sub_all = indiv_df[indiv_df["EO"] == eo]
        sub = sub_all[sub_all["n_functional_copies"] > 0]
        n_total = len(sub_all)
        n_fun = len(sub)
        n_sc = int((sub_all["n_functional_copies"] == 0).sum())
        if n_fun == 0:
            continue
        genotypes = np.stack(sub["genotype"].values, axis=0)
        fg_copies_stack = np.stack(sub["fg_copies"].values, axis=0)
        local_f_obs = compute_local_f(fg_copies_stack)
        present_mask = local_f_obs > 0
        n_distinct_obs = int(present_mask.sum())

        # Observed P_compat: real mothers x observed local Fg frequencies
        obs_mean, obs_lo, obs_hi = bootstrap_pcompat_mean(
            genotypes, local_f_obs, class_i_mask, zygosity_probs,
            N_MC_FATHERS, N_BOOTSTRAP, rng)
        # Predicted-at-EO-scale (species-wide P1 fathers, same mothers)
        pred_adult_mean, pred_adult_lo, pred_adult_hi = bootstrap_pcompat_mean(
            genotypes, prior_f, class_i_mask, zygosity_probs,
            N_MC_FATHERS, N_BOOTSTRAP, rng)

        # Phase 5 per-location prediction (component-based)
        pc_row = pred_pcompat[pred_pcompat["locationCode"] == eo]
        div_row = pred_diversity[pred_diversity["locationCode"] == eo]
        if pc_row.empty or div_row.empty:
            print(f"[step30d] {eo} not found in Phase 5 prediction tables "
                  "— skipping.")
            continue
        pc_row = pc_row.iloc[0]
        div_row = div_row.iloc[0]

        # Sample-size-adjusted upper bound — what we would detect at
        # 4·n_fun adult allele draws from P1 if there were NO drift.
        # The observed distinct count sits BELOW this ceiling for two
        # reasons: finite-sample coupon-collecting (small) and genuine
        # drift collapse (large). At n ≥ 40 adults from P1 the
        # upper bound is close to the 32-Fg species ceiling, so any
        # large deficit is drift, not sampling.
        n_allele_draws_adults = 4 * n_fun
        exp_distinct_P1 = float(
            (1.0 - (1.0 - prior_f) ** n_allele_draws_adults).sum())

        rows.append({
            "locationCode":            eo,
            "locationID":              int(pc_row["locationID"]),
            "EO":                      eo,
            "n_components_50m":        int(pc_row.get("n_components_50m", 1)),
            "N_fertile_effective":     int(pc_row["N_fertile_effective"]),
            "n_adults_total":          n_total,
            "n_functional_carriers":   n_fun,
            "n_SC_fully_null":         n_sc,
            # Observed
            "obs_distinct_Fgs":        n_distinct_obs,
            "obs_p_ClassI":            float(local_f_obs[class_i_mask].sum()),
            # Sample-size-adjusted no-drift upper bound (species-wide
            # P1 at 4·n_fun adult allele draws). Deficit from this
            # upper bound to the observed count isolates DRIFT from
            # SAMPLING.
            "exp_distinct_adult_sample_from_P1_upper": exp_distinct_P1,
            "obs_pcompat_mean":        obs_mean,
            "obs_pcompat_lo95":        obs_lo,
            "obs_pcompat_hi95":        obs_hi,
            # Predicted at adult scale (same mothers, P1 fathers)
            "pred_adult_pcompat_mean": pred_adult_mean,
            "pred_adult_pcompat_lo95": pred_adult_lo,
            "pred_adult_pcompat_hi95": pred_adult_hi,
            # Phase 5 per-location prediction (component-based union
            # for diversity, size-weighted mean for P_compat)
            "pred_phase5_distinct_Fgs":        float(
                div_row["predicted_local_pool_size_mean"]),
            "pred_phase5_distinct_Fgs_lo":     float(
                div_row["predicted_local_pool_size_lo"]),
            "pred_phase5_distinct_Fgs_hi":     float(
                div_row["predicted_local_pool_size_hi"]),
            "pred_phase5_pcompat_mean":        float(
                pc_row["predicted_P_compat_mean"]),
            "pred_phase5_pcompat_lo95":        float(
                pc_row["predicted_P_compat_lo"]),
            "pred_phase5_pcompat_hi95":        float(
                pc_row["predicted_P_compat_hi"]),
        })
    return pd.DataFrame(rows)


def per_fg_frequencies(indiv_df: pd.DataFrame,
                       prior: pd.DataFrame,
                       eos: List[str]) -> pd.DataFrame:
    """Long table — one row per (location, Fg) with observed f and
    species-wide P1 f. Keep ALL 32 Fgs per location so missing-at-
    location entries get f_observed = 0 (drift signature)."""
    prior_f = prior["f_mean"].values
    fg_labels = prior["Fg"].astype(str).tolist()
    rows = []
    for eo in eos:
        sub = indiv_df[(indiv_df["EO"] == eo)
                       & (indiv_df["n_functional_copies"] > 0)]
        if len(sub) == 0:
            continue
        fg_stack = np.stack(sub["fg_copies"].values, axis=0)
        local_f = compute_local_f(fg_stack)
        for j, fg in enumerate(fg_labels):
            rows.append({
                "locationCode":    eo,
                "EO":              eo,
                "Fg":              fg,
                "f_observed":      float(local_f[j]),
                "f_P1_species":    float(prior_f[j]),
                "present_in_location": bool(local_f[j] > 0),
            })
    return pd.DataFrame(rows)


def draw_figure(stats: pd.DataFrame,
                freqs: pd.DataFrame,
                bands: dict,
                prior: pd.DataFrame,
                out_png: Path, out_pdf: Path) -> None:
    """Two-row figure.

    Top row (two square panels):
      A — predicted vs observed SRK diversity per location
      B — predicted vs observed pollen compatibility per location
    Bottom row (per-location drift residual, full width):
      C — per-Fg (observed − species-wide P1) frequency deviation,
          one lane per clean-overlap location. Positive bars = Fgs
          drift-enriched beyond P1; negative bars = depleted or
          absent at the location (direct visual of drift divergence).
    """
    stats = stats.sort_values("n_functional_carriers", ascending=False)
    fg_labels = prior["Fg"].astype(str).tolist()
    n_fg = len(fg_labels)
    prior_order = (prior.sort_values("f_mean", ascending=False)
                        ["Fg"].astype(str).tolist())
    fg_pos = {fg: i for i, fg in enumerate(prior_order)}

    fig = plt.figure(figsize=(13.5, 8.5))
    gs = fig.add_gridspec(
        2, 2,
        height_ratios=[1.0, 1.15],
        width_ratios=[1.0, 1.0],
        hspace=0.42, wspace=0.32,
    )
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC_gs = gs[1, :].subgridspec(len(stats), 1, hspace=0.18)
    axC_rows = [fig.add_subplot(axC_gs[i, 0]) for i in range(len(stats))]

    # ------------------------------------------------------------------
    # Panel A — SRK diversity: predicted vs observed scatter (clean)
    # ------------------------------------------------------------------
    hiA = 1.08 * max(stats["obs_distinct_Fgs"].max(),
                     stats["pred_phase5_distinct_Fgs_hi"].max(),
                     n_fg)
    axA.plot([0, hiA], [0, hiA], color="#888888", lw=1.0, ls="--",
             alpha=0.6, zorder=1)
    axA.axvline(n_fg, color="#aaaaaa", ls=":", lw=0.8, alpha=0.5)
    axA.axhline(n_fg, color="#aaaaaa", ls=":", lw=0.8, alpha=0.5)
    for _, r in stats.iterrows():
        axA.errorbar(r["pred_phase5_distinct_Fgs"], r["obs_distinct_Fgs"],
                     xerr=[[r["pred_phase5_distinct_Fgs"]
                            - r["pred_phase5_distinct_Fgs_lo"]],
                           [r["pred_phase5_distinct_Fgs_hi"]
                            - r["pred_phase5_distinct_Fgs"]]],
                     fmt="s", color=COLOUR_OBS, ecolor="#999999",
                     elinewidth=1.0, capsize=3, markersize=8, zorder=2)
        label = (f"{r['locationCode']}  (n={int(r['n_functional_carriers'])}, "
                 f"obs {int(r['obs_distinct_Fgs'])}/{n_fg})")
        if r["pred_phase5_distinct_Fgs"] > 0.65 * hiA:
            axA.annotate(
                label + "  ",
                xy=(r["pred_phase5_distinct_Fgs"], r["obs_distinct_Fgs"]),
                fontsize=9, va="center", ha="right")
        else:
            axA.annotate(
                "  " + label,
                xy=(r["pred_phase5_distinct_Fgs"], r["obs_distinct_Fgs"]),
                fontsize=9, va="center", ha="left")
    axA.set_xlim(0, hiA); axA.set_ylim(0, hiA)
    axA.set_xlabel("Phase 5 predicted (union across 50 m components)",
                   fontsize=10)
    axA.set_ylabel("Observed in adults (distinct Fgs)", fontsize=10)
    axA.set_title(f"A  SRK diversity per location (ceiling = {n_fg} Fgs)",
                   fontsize=11, loc="left")
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)

    # ------------------------------------------------------------------
    # Panel B — Pollen compatibility: predicted vs observed scatter
    # ------------------------------------------------------------------
    hiB = max(0.05,
             1.08 * max(stats["obs_pcompat_hi95"].max(),
                        stats["pred_phase5_pcompat_hi95"].max(),
                        bands["species_mean"]))
    axB.add_patch(Rectangle((0, 0), hiB, bands["failed_max"],
                             color=COLOUR_FAILED, alpha=0.10, zorder=0))
    axB.add_patch(Rectangle((0, bands["failed_max"]), hiB,
                             bands["struggling_max"] - bands["failed_max"],
                             color=COLOUR_STRUGGLING, alpha=0.10, zorder=0))
    axB.add_patch(Rectangle((0, bands["struggling_max"]), hiB,
                             hiB - bands["struggling_max"],
                             color=COLOUR_SUSTAIN, alpha=0.10, zorder=0))
    axB.plot([0, hiB], [0, hiB], color="#888888", lw=1.0, ls="--",
             alpha=0.6, zorder=1)
    for _, r in stats.iterrows():
        axB.errorbar(r["pred_phase5_pcompat_mean"], r["obs_pcompat_mean"],
                     xerr=[[r["pred_phase5_pcompat_mean"]
                            - r["pred_phase5_pcompat_lo95"]],
                           [r["pred_phase5_pcompat_hi95"]
                            - r["pred_phase5_pcompat_mean"]]],
                     yerr=[[r["obs_pcompat_mean"] - r["obs_pcompat_lo95"]],
                           [r["obs_pcompat_hi95"] - r["obs_pcompat_mean"]]],
                     fmt="o", color=COLOUR_OBS, ecolor="#999999",
                     elinewidth=1.0, capsize=3, markersize=8, zorder=2)
        # Flip label to the left of the dot when the dot sits in the
        # right third of the axis, otherwise it clips.
        label = r["locationCode"]
        if r["pred_phase5_pcompat_mean"] > 0.65 * hiB:
            axB.annotate(
                label + "  ",
                xy=(r["pred_phase5_pcompat_mean"], r["obs_pcompat_mean"]),
                fontsize=9, va="center", ha="right")
        else:
            axB.annotate(
                "  " + label,
                xy=(r["pred_phase5_pcompat_mean"], r["obs_pcompat_mean"]),
                fontsize=9, va="center", ha="left")
    axB.set_xlim(0, hiB); axB.set_ylim(0, hiB)
    axB.set_xlabel("Phase 5 predicted (component-weighted mean)",
                   fontsize=10)
    axB.set_ylabel("Observed in adults", fontsize=10)
    axB.set_title("B  Pollen compatibility per location",
                   fontsize=11, loc="left")
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)

    # ------------------------------------------------------------------
    # Panel C — Drift residual: f_observed − f_P1 per Fg, per location
    # ------------------------------------------------------------------
    # Fgs sorted by P1 descending (most-common on the left, matches
    # the reader's intuition that drift most dramatically alters the
    # tails of the distribution).
    p1_order = np.array([prior.set_index("Fg")["f_mean"]
                          .to_dict().get(fg, 0.0)
                          for fg in prior_order])
    x = np.arange(n_fg)
    palette = ["#2b6cb0", "#c94b4b", "#3c8f4c"]
    y_abs_max = 0.0
    locs = stats["locationCode"].tolist()
    for ax_i, loc in zip(axC_rows, locs):
        sub = freqs[freqs["locationCode"] == loc].set_index("Fg")
        f_loc = np.array([sub.loc[fg, "f_observed"] if fg in sub.index else 0.0
                           for fg in prior_order])
        residual = f_loc - p1_order
        y_abs_max = max(y_abs_max, float(np.abs(residual).max()))

        # Signed bars: green = enriched (observed > P1); red = depleted.
        colours = ["#2ca25f" if v > 0 else "#de2d26" for v in residual]
        ax_i.bar(x, residual, color=colours, edgecolor="white",
                  linewidth=0.4)
        ax_i.axhline(0.0, color="#333333", lw=0.6)
        # Mark Fgs that are ABSENT at the location with a short tick
        # at the x-axis so readers see which species-wide Fgs the
        # location has lost entirely.
        absent_idx = np.where(f_loc == 0.0)[0]
        ax_i.scatter(absent_idx,
                      np.full_like(absent_idx, -0.005, dtype=float),
                      marker="x", color="#555555", s=16, zorder=3,
                      label="Absent at location")
        ax_i.set_xlim(-0.6, n_fg - 0.4)
        ax_i.set_xticks(x)
        if ax_i is axC_rows[-1]:
            ax_i.set_xticklabels(prior_order, rotation=90, fontsize=7)
            ax_i.set_xlabel(
                "32 Fgs, sorted by species-wide (P1) frequency, "
                "most-common → rarest. × marker = Fg absent at this location.",
                fontsize=9)
        else:
            ax_i.set_xticklabels([])
        n_fun = int(stats.set_index("locationCode").loc[loc,
                                                         "n_functional_carriers"])
        obs_fgs = int(stats.set_index("locationCode").loc[loc,
                                                           "obs_distinct_Fgs"])
        ax_i.text(0.995, 0.90,
                   f"{loc}  (n={n_fun} adults, {obs_fgs}/{n_fg} Fgs present)",
                   transform=ax_i.transAxes,
                   fontsize=10, ha="right", va="top", fontweight="bold")
        ax_i.set_ylabel("f_obs − f_P1", fontsize=9)
        ax_i.spines["top"].set_visible(False)
        ax_i.spines["right"].set_visible(False)

    # Symmetric y-range so enrichment and depletion are directly
    # comparable across the three locations.
    y_lim = 1.15 * y_abs_max
    for ax_i in axC_rows:
        ax_i.set_ylim(-y_lim, y_lim)

    axC_rows[0].set_title(
        "C  Drift residual per Fg per location  "
        "(green = observed > P1; red = observed < P1; × = absent)",
        fontsize=11, loc="left")

    fig.suptitle(
        "Phase 5 Part C — Clean-overlap EOs (1:1 with Phase 5 location)  "
        f"|  Phase 4 adult SRK genotypes  |  species-mean P_compat = "
        f"{bands['species_mean']:.3f}",
        fontsize=12, y=0.995,
    )
    fig.subplots_adjust(left=0.07, right=0.98, top=0.94, bottom=0.08)
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eos", type=str, default=",".join(CLEAN_OVERLAP_EOS),
                    help="Comma-separated EOs to include (defaults to the "
                         "three clean-overlap EOs).")
    args = ap.parse_args()
    eos = [e.strip() for e in args.eos.split(",") if e.strip()]

    rng = np.random.default_rng(RNG_SEED)
    prior = pd.read_csv(DEFAULT_PRIOR_TSV, sep="\t", encoding="utf-8-sig")
    fg_labels = prior["Fg"].astype(str).tolist()
    allele_to_fg = load_allele_to_fg(DEFAULT_ALLELE_FG_TSV)
    indiv_df = load_individuals(DEFAULT_INDIV_TSV, allele_to_fg, fg_labels)
    print(f"[step30d] Loaded {len(indiv_df)} adults; targeting EOs: "
          f"{', '.join(eos)}")

    pred_pcompat = pd.read_csv(DEFAULT_PRED_PCOMPAT, sep="\t",
                                encoding="utf-8-sig")
    pred_diversity = pd.read_csv(DEFAULT_PRED_DIVERSITY, sep="\t",
                                  encoding="utf-8-sig")

    class_map = load_class_map()
    class_i_mask = build_class_i_mask(fg_labels, class_map)
    zygosity_probs = load_zygosity_dist()
    bands_df = pd.read_csv(DEFAULT_BANDS_TSV, sep="\t", encoding="utf-8-sig")
    bands = {
        "species_mean":   float(bands_df["species_mean"].iloc[0]),
        "failed_max":     float(bands_df["failed_max"].iloc[0]),
        "struggling_max": float(bands_df["struggling_max"].iloc[0]),
    }

    stats = per_location_stats(indiv_df, prior, class_i_mask,
                                zygosity_probs, pred_pcompat,
                                pred_diversity, rng, eos)
    freqs = per_fg_frequencies(indiv_df, prior, eos)
    DEFAULT_TABLES.mkdir(parents=True, exist_ok=True)
    DEFAULT_FIGURES.mkdir(parents=True, exist_ok=True)

    stats_path = DEFAULT_TABLES / "step30_B_partC_clean_overlap_per_location.tsv"
    freqs_path = DEFAULT_TABLES / "step30_B_partC_clean_overlap_fg_frequencies.tsv"
    fig_png = DEFAULT_FIGURES / "step30_B_partC_clean_overlap.png"
    fig_pdf = DEFAULT_FIGURES / "step30_B_partC_clean_overlap.pdf"

    stats.to_csv(stats_path, sep="\t", index=False)
    freqs.to_csv(freqs_path, sep="\t", index=False)
    print(f"[step30d] Wrote {stats_path} ({len(stats)} rows)")
    print(f"[step30d] Wrote {freqs_path} ({len(freqs)} rows)")

    draw_figure(stats, freqs, bands, prior, fig_png, fig_pdf)
    print(f"[step30d] Wrote {fig_png} + .pdf")

    # Short console summary — observed vs Phase 5 predicted at the
    # three clean locations.
    print()
    print("[step30d] ============== SUMMARY — observed vs Phase 5 predicted ==============")
    for _, r in stats.iterrows():
        print(f"[step30d] {r['locationCode']}  (adults = "
              f"{int(r['n_functional_carriers'])})")
        print(f"[step30d]   SRK diversity (distinct Fgs):")
        print(f"[step30d]      observed in adults          = "
              f"{int(r['obs_distinct_Fgs'])} of 32")
        print(f"[step30d]      no-drift upper bound at n   = "
              f"{r['exp_distinct_adult_sample_from_P1_upper']:.1f} of 32  "
              f"(if P1 held locally)")
        print(f"[step30d]      Phase 5 predicted (full location) = "
              f"{r['pred_phase5_distinct_Fgs']:.1f}  "
              f"[{r['pred_phase5_distinct_Fgs_lo']:.1f}, "
              f"{r['pred_phase5_distinct_Fgs_hi']:.1f}]")
        print(f"[step30d]   Pollen compatibility:")
        print(f"[step30d]      observed in adults          = "
              f"{r['obs_pcompat_mean']:.3f}  "
              f"[{r['obs_pcompat_lo95']:.3f}, {r['obs_pcompat_hi95']:.3f}]")
        print(f"[step30d]      Phase 5 predicted (component) = "
              f"{r['pred_phase5_pcompat_mean']:.3f}  "
              f"[{r['pred_phase5_pcompat_lo95']:.3f}, "
              f"{r['pred_phase5_pcompat_hi95']:.3f}]")


if __name__ == "__main__":
    main()
