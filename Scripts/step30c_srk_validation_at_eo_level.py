"""Step 30 Part C — Empirical validation at EO scale.

Uses Canu-amplicon adult SRK genotypes (Phase 4 Step 23) to compute
OBSERVED per-EO P_compat under the sporophytic + empirical-zygosity
Class I / II model (§ A.7 of the Phase 5 doc) and compare to what the
finite-population P1 model predicts at the same scale.

Data scope
----------
Per-individual allele calls come from
`Tables/Phase4/step23_individual_allele_genotypes_with_nulls.tsv` (368
individuals, 23 EOs, 49 functional alleles + Allele_NULL). The
allele → Fg mapping comes from `Tables/Phase5/step26i_L1_carrier_inventory.tsv`
(the same inventory that P1 was built from — 49 alleles across 32
Fgs). Inclusion: EOs with n ≥ 10 individuals having ≥ 1 functional
copy. Individuals with 0 functional copies (fully-null tetraploids)
are Self-Compatible by definition and reported separately.

The predicted P_compat uses P1 and empirical zygosity as father-drawing
priors, applied to the actual observed mothers at each EO. The
observed P_compat replaces P1 with the OBSERVED local Fg frequencies
at each EO. Deviation between the two isolates the effect of local
frequency drift on random-mating compatibility, holding the mothers
themselves fixed.

Caveat
------
P1 was built from these same individuals, so the species-mean is
guaranteed to line up on average. What this figure genuinely tests is
whether local (per-EO) frequency skew moves P_compat in the predicted
direction — i.e. whether the sporophytic Class I / II model is
robust to per-EO drift, not whether it is calibrated at the species
level. EO ≠ location — this is an intermediate scale of validation
until seed-derived per-location genotypes come in.
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
    p_compat_sporophytic_batch,
    p_compat_sporophytic_empirical,
    sample_genotypes_empirical,
    species_mean_p_compat_empirical,
    traffic_light_bands,
)

# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
DEFAULT_TABLES  = Path("Tables/Phase5")
DEFAULT_FIGURES = Path("Figures/Phase5")
DEFAULT_ALLELE_FG_TSV = DEFAULT_TABLES / "step26i_L1_carrier_inventory.tsv"
DEFAULT_INDIV_TSV = Path(
    "Tables/Phase4/step23_individual_allele_genotypes_with_nulls.tsv")
DEFAULT_PRIOR_TSV = DEFAULT_TABLES / "step30_A_prediction_prior_frequencies.tsv"

MIN_N_INDIV_PER_EO = 10
N_MC_FATHERS       = 800
N_BOOTSTRAP        = 400
RNG_SEED           = 2032

# Traffic-light colours (match § A.7 palette)
COLOUR_FAILED      = "#c94b4b"
COLOUR_STRUGGLING  = "#e6b325"
COLOUR_SUSTAIN     = "#3c8f4c"
COLOUR_PRED        = "#2b6cb0"
COLOUR_OBS         = "#333333"


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------
def load_allele_to_fg(tsv: Path) -> Dict[str, str]:
    df = pd.read_csv(tsv, sep="\t", encoding="utf-8-sig")
    return dict(zip(df["Allele"].astype(str), df["Fg"].astype(str)))


def load_prior(tsv: Path) -> pd.DataFrame:
    return pd.read_csv(tsv, sep="\t", encoding="utf-8-sig")


def load_individuals(tsv: Path,
                     allele_to_fg: Dict[str, str],
                     fg_labels: List[str]) -> pd.DataFrame:
    """Parse per-individual allele calls into per-individual Fg genotypes.

    Returns rows with:
      Individual, EO, BL, SI_status, fg_copies (K_fg,), genotype (4,),
      n_functional_copies, n_null_copies, n_distinct_functional
    Individuals with 0 functional copies are kept but excluded from
    P_compat computations (their genotype vector is all -1).
    """
    df = pd.read_csv(tsv, sep="\t", encoding="utf-8-sig")
    allele_cols = [c for c in df.columns
                   if c.startswith("Allele_") and c != "Allele_NULL"]
    fg_to_idx = {fg: i for i, fg in enumerate(fg_labels)}
    K_fg = len(fg_labels)
    records = []
    for _, row in df.iterrows():
        n_null = int(row.get("Allele_NULL", 0)) if pd.notna(
            row.get("Allele_NULL", 0)) else 0
        fg_copies = np.zeros(K_fg, dtype=int)
        for a in allele_cols:
            v = row[a]
            if pd.isna(v):
                continue
            c = int(v)
            if c <= 0:
                continue
            fg = allele_to_fg.get(a)
            if fg is None or fg not in fg_to_idx:
                continue
            fg_copies[fg_to_idx[fg]] += c
        n_functional = int(fg_copies.sum())
        distinct_present = np.where(fg_copies > 0)[0]
        raw: list[int] = []
        for j in distinct_present:
            raw.extend([int(j)] * int(fg_copies[j]))
        if n_functional == 0:
            padded = np.full(4, -1, dtype=int)
        else:
            while len(raw) < 4:
                raw.append(raw[0])
            padded = np.array(raw[:4], dtype=int)
        records.append({
            "Individual": row["Individual"],
            "EO":         str(row.get("EO_normalised", "Unassigned")),
            "BL":         str(row.get("BL_inferred", "Unassigned")),
            "SI_status":  row.get("SI_status", ""),
            "fg_copies":  fg_copies,
            "genotype":   padded,
            "n_functional_copies":    n_functional,
            "n_null_copies":          n_null,
            "n_distinct_functional":  int((fg_copies > 0).sum()),
        })
    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Per-EO computations
# ---------------------------------------------------------------------------
def compute_local_f(fg_copies_stack: np.ndarray) -> np.ndarray:
    """Aggregate per-individual Fg copy counts to a normalised local Fg
    frequency vector (functional copies only)."""
    totals = fg_copies_stack.sum(axis=0)
    s = totals.sum()
    if s <= 0:
        return None
    return totals / s


def bootstrap_pcompat_mean(genotypes: np.ndarray,
                            local_f: np.ndarray,
                            class_i_mask: np.ndarray,
                            zygosity_probs: np.ndarray,
                            n_fathers: int,
                            n_boot: int,
                            rng: np.random.Generator) -> tuple[float, float, float]:
    """Bootstrap 95 % CI on the mean P_compat over the observed mothers.
    Fathers are re-drawn each bootstrap replicate (father-pool + mother
    resampling both contribute to uncertainty)."""
    n = genotypes.shape[0]
    means = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        p = p_compat_sporophytic_empirical(
            genotypes[idx], local_f, class_i_mask,
            zygosity_probs, n_fathers=n_fathers, rng=rng)
        means[b] = float(p.mean())
    return float(means.mean()), float(np.quantile(means, 0.025)), float(
        np.quantile(means, 0.975))


def per_eo_stats(indiv_df: pd.DataFrame,
                  prior_f: np.ndarray,
                  class_i_mask: np.ndarray,
                  zygosity_probs: np.ndarray,
                  rng: np.random.Generator,
                  min_n: int = MIN_N_INDIV_PER_EO) -> pd.DataFrame:
    """One row per EO with n ≥ min_n functional-carrier individuals."""
    K_fg = prior_f.size
    rows = []
    eos = sorted(indiv_df["EO"].unique())
    for eo in eos:
        sub_all = indiv_df[indiv_df["EO"] == eo]
        sub = sub_all[sub_all["n_functional_copies"] > 0]
        n_total = len(sub_all)
        n_fun   = len(sub)
        n_sc    = int((sub_all["n_functional_copies"] == 0).sum())
        if n_fun < min_n:
            continue
        genotypes = np.stack(sub["genotype"].values, axis=0)
        fg_copies_stack = np.stack(sub["fg_copies"].values, axis=0)
        local_f_obs = compute_local_f(fg_copies_stack)

        # ---- Observed: real mothers × observed local Fg frequencies ----
        p_obs = p_compat_sporophytic_empirical(
            genotypes, local_f_obs, class_i_mask,
            zygosity_probs, n_fathers=N_MC_FATHERS, rng=rng)
        obs_mean, obs_lo, obs_hi = bootstrap_pcompat_mean(
            genotypes, local_f_obs, class_i_mask, zygosity_probs,
            N_MC_FATHERS, N_BOOTSTRAP, rng)

        # ---- Predicted: real mothers × P1 species-wide frequencies ----
        p_pred = p_compat_sporophytic_empirical(
            genotypes, prior_f, class_i_mask,
            zygosity_probs, n_fathers=N_MC_FATHERS, rng=rng)
        pred_mean, pred_lo, pred_hi = bootstrap_pcompat_mean(
            genotypes, prior_f, class_i_mask, zygosity_probs,
            N_MC_FATHERS, N_BOOTSTRAP, rng)

        # ---- Analytical observed (for reference) ----
        p_obs_anly = p_compat_sporophytic_batch(
            genotypes, local_f_obs, class_i_mask).mean()

        # ---- Observed zygosity distribution ----
        zyg_counts = np.bincount(
            sub["n_distinct_functional"].values, minlength=5)[:5]
        zyg_frac = zyg_counts / max(zyg_counts.sum(), 1)

        # ---- Observed Class I mass ----
        p_I_obs = float(local_f_obs[class_i_mask].sum())

        rows.append({
            "EO":                     eo,
            "n_total":                n_total,
            "n_functional_carriers":  n_fun,
            "n_SC_fully_null":        n_sc,
            "n_distinct_1_frac":      float(zyg_frac[1]),
            "n_distinct_2_frac":      float(zyg_frac[2]),
            "n_distinct_3_frac":      float(zyg_frac[3]),
            "n_distinct_4_frac":      float(zyg_frac[4]),
            "p_ClassI_obs":           p_I_obs,
            "obs_pcompat_mean":       obs_mean,
            "obs_pcompat_lo95":       obs_lo,
            "obs_pcompat_hi95":       obs_hi,
            "obs_pcompat_analytical": float(p_obs_anly),
            "pred_pcompat_mean":      pred_mean,
            "pred_pcompat_lo95":      pred_lo,
            "pred_pcompat_hi95":      pred_hi,
        })
    return pd.DataFrame(rows).sort_values(
        "n_functional_carriers", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def draw_figure(stats: pd.DataFrame,
                 species_mean: float,
                 bands: dict,
                 zygosity_species: np.ndarray,
                 out_path: Path) -> None:
    """Two-panel figure — (A) observed vs predicted P_compat scatter,
    (B) observed distinct-identity distribution per EO with species-wide
    reference lines."""
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.6))

    # ---- Panel A: observed vs predicted scatter ----
    lo = 0.0
    hi = max(0.05, 1.05 * max(stats["obs_pcompat_hi95"].max(),
                              stats["pred_pcompat_hi95"].max(),
                              bands["struggling_max"]))

    # Traffic-light bands (horizontal, on Y = observed axis)
    axA.axhspan(0, bands["failed_max"],
                color=COLOUR_FAILED, alpha=0.10, zorder=0,
                label=f"Failed (< {bands['failed_max']:.2f})")
    axA.axhspan(bands["failed_max"], bands["struggling_max"],
                color=COLOUR_STRUGGLING, alpha=0.10, zorder=0,
                label=f"Struggling ({bands['failed_max']:.2f}–{bands['struggling_max']:.2f})")
    axA.axhspan(bands["struggling_max"], 1.0,
                color=COLOUR_SUSTAIN, alpha=0.08, zorder=0,
                label=f"Sustainable (> {bands['struggling_max']:.2f})")

    # 1:1 diagonal
    axA.plot([lo, hi], [lo, hi], "--", color="#777", lw=1.0, zorder=1,
             label="1:1 (perfect match)")
    # Species-mean cross-hair
    axA.axvline(species_mean, color=COLOUR_PRED, lw=0.8, alpha=0.5, zorder=1)
    axA.axhline(species_mean, color=COLOUR_PRED, lw=0.8, alpha=0.5, zorder=1,
                label=f"Species mean = {species_mean:.2f}")

    # Per-EO dots with error bars in both axes
    for _, r in stats.iterrows():
        size = 40 + 6 * np.sqrt(r["n_functional_carriers"])
        axA.errorbar(r["pred_pcompat_mean"], r["obs_pcompat_mean"],
                     xerr=[[r["pred_pcompat_mean"] - r["pred_pcompat_lo95"]],
                           [r["pred_pcompat_hi95"] - r["pred_pcompat_mean"]]],
                     yerr=[[r["obs_pcompat_mean"] - r["obs_pcompat_lo95"]],
                           [r["obs_pcompat_hi95"] - r["obs_pcompat_mean"]]],
                     fmt="o", color=COLOUR_OBS, ecolor="#888",
                     markersize=np.sqrt(size), zorder=3, capsize=2)
        axA.annotate(r["EO"], xy=(r["pred_pcompat_mean"], r["obs_pcompat_mean"]),
                     xytext=(6, 4), textcoords="offset points",
                     fontsize=9, zorder=4)

    axA.set_xlim(lo, hi)
    axA.set_ylim(lo, hi)
    axA.set_xlabel(
        "Predicted mean pollen compatibility per EO\n"
        "(real mothers vs simulated fathers drawn from species-wide prior)")
    axA.set_ylabel(
        "Observed mean pollen compatibility per EO\n"
        "(real mothers vs simulated fathers drawn from observed local frequencies)")
    axA.set_title("A. EO-level compatibility — observed vs predicted")
    axA.legend(loc="lower right", fontsize=8, framealpha=0.85)
    axA.grid(alpha=0.25)

    # ---- Panel B: observed zygosity distribution per EO ----
    eo_names = stats["EO"].tolist()
    x = np.arange(len(eo_names))
    w = 0.22
    axB.bar(x - 1.5 * w, stats["n_distinct_1_frac"], w,
            color="#0e5c8f", label="1 identity (homozygote-like)")
    axB.bar(x - 0.5 * w, stats["n_distinct_2_frac"], w,
            color="#3aa7e0", label="2 identities")
    axB.bar(x + 0.5 * w, stats["n_distinct_3_frac"], w,
            color="#f6c76a", label="3 identities")
    axB.bar(x + 1.5 * w, stats["n_distinct_4_frac"], w,
            color="#c94b4b", label="4 identities")

    # Species-wide reference lines
    for k, y in enumerate(zygosity_species, start=1):
        axB.axhline(y, ls=":", lw=0.8, color="#555", alpha=0.6)
    axB.text(len(eo_names) - 0.5, zygosity_species[0] + 0.02,
             f"Species-wide (Step 23) — 1 id = {zygosity_species[0]:.0%}, "
             f"2 = {zygosity_species[1]:.0%}, 3 = {zygosity_species[2]:.0%}",
             fontsize=8, ha="right", color="#555")

    axB.set_xticks(x)
    axB.set_xticklabels([f"{eo}\n(n = {n})" for eo, n in zip(
        eo_names, stats["n_functional_carriers"])], fontsize=9)
    axB.set_ylim(0, 1.0)
    axB.set_ylabel("Fraction of individuals")
    axB.set_title("B. Distinct functional SRK identities per plant, by EO")
    axB.legend(loc="upper right", fontsize=8, framealpha=0.85)
    axB.grid(alpha=0.25, axis="y")

    fig.suptitle(
        "Step 30 Part C — First empirical validation of the sporophytic "
        "Class I / II + empirical-zygosity P_compat model at EO scale\n"
        "(Canu-amplicon adult SRK genotypes, Phase 4 Step 23)",
        fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_path, dpi=200)
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-n", type=int, default=MIN_N_INDIV_PER_EO,
                    help="Minimum functional-carrier individuals per EO "
                         "(default: %(default)s)")
    ap.add_argument("--indiv-tsv", type=Path, default=DEFAULT_INDIV_TSV)
    ap.add_argument("--allele-fg-tsv", type=Path, default=DEFAULT_ALLELE_FG_TSV)
    ap.add_argument("--prior-tsv", type=Path, default=DEFAULT_PRIOR_TSV)
    args = ap.parse_args()

    tables_dir  = DEFAULT_TABLES;  tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = DEFAULT_FIGURES; figures_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)

    # ---- Load prior + Fg labels + Class map + zygosity ----
    prior = load_prior(args.prior_tsv)
    fg_labels = prior["Fg"].astype(str).tolist()
    prior_f = prior["f_mean"].values
    class_map = load_class_map()
    class_i_mask = build_class_i_mask(fg_labels, class_map)
    zygosity_probs = load_zygosity_dist()

    # ---- Load individuals ----
    allele_to_fg = load_allele_to_fg(args.allele_fg_tsv)
    indiv_df = load_individuals(args.indiv_tsv, allele_to_fg, fg_labels)
    n_all  = len(indiv_df)
    n_fun  = int((indiv_df["n_functional_copies"] > 0).sum())
    n_sc   = int((indiv_df["n_functional_copies"] == 0).sum())
    print(f"[step30c] Loaded {n_all} individuals across "
          f"{indiv_df['EO'].nunique()} EOs "
          f"({n_fun} with ≥1 functional copy, {n_sc} fully null).")

    # ---- Species mean under P1 + empirical zygosity ----
    species_mean = species_mean_p_compat_empirical(
        prior_f, class_i_mask, zygosity_probs,
        n_mothers=8_000, n_fathers=N_MC_FATHERS, rng=rng)
    bands = traffic_light_bands(species_mean)
    print(f"[step30c] Species-mean P_compat under P1 + empirical zygosity: "
          f"{species_mean:.3f} (bands: "
          f"failed < {bands['failed_max']:.3f}, "
          f"struggling < {bands['struggling_max']:.3f}).")

    # ---- Per-EO stats ----
    stats = per_eo_stats(indiv_df, prior_f, class_i_mask,
                          zygosity_probs, rng, min_n=args.min_n)
    if stats.empty:
        raise SystemExit(f"[step30c] No EO meets the n ≥ {args.min_n} threshold.")
    print(f"[step30c] {len(stats)} EOs meet the n ≥ {args.min_n} threshold.")

    # Add species reference for downstream readability
    stats["species_mean_pcompat"] = species_mean
    stats["band_failed_max"]      = bands["failed_max"]
    stats["band_struggling_max"]  = bands["struggling_max"]

    out_tsv = tables_dir / "step30_C_pcompat_validation_at_eo.tsv"
    stats.to_csv(out_tsv, sep="\t", index=False,
                 float_format="%.4f")
    print(f"[step30c] Wrote {out_tsv}")

    # ---- Figure ----
    out_png = figures_dir / "step30_C_pcompat_observed_vs_predicted.png"
    draw_figure(stats, species_mean, bands, zygosity_probs, out_png)
    print(f"[step30c] Wrote {out_png}")

    # ---- Console summary ----
    print("\n[step30c] Per-EO summary:")
    with pd.option_context("display.max_columns", 20,
                            "display.width", 180,
                            "display.float_format", "{:.3f}".format):
        print(stats[["EO", "n_functional_carriers",
                     "obs_pcompat_mean", "obs_pcompat_lo95", "obs_pcompat_hi95",
                     "pred_pcompat_mean", "pred_pcompat_lo95", "pred_pcompat_hi95",
                     "p_ClassI_obs",
                     "n_distinct_1_frac", "n_distinct_2_frac"]])


if __name__ == "__main__":
    main()
