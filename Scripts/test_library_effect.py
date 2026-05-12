#!/usr/bin/env python3
"""
test_library_effect.py
======================

Tests whether **library identity** introduces systematic technical bias into
the SRK pipeline results, with particular focus on the most recent libraries
(009 and 010) which have different processing histories from libraries 001–008.

Why this matters
----------------
Library 009 was assembled before Step 4b (BLAST coverage filter) and Step 7b
(N-content filter) were introduced — its sequences passed through the legacy
pipeline (length filter only). Library 010 was the first to be processed
through the full new pipeline. A "library effect" — i.e., systematic
differences between libraries that go beyond random sampling — would
contaminate every downstream biological inference (genotype distributions,
SI-functional status, allele richness, the candidate SI-escape signal).

Test logic
----------
The cleanest test is **within-EO between-library** comparison. If multiple
libraries have sampled the SAME Element Occurrence (EO), their genotype and
SI-status distributions should be statistically indistinguishable (under the
null hypothesis of no library effect). Significant differences within an EO
indicate technical bias rather than real population signal.

Tests performed:

  Test 1 — Global Library × SI_functional_status independence
              (chi-square: are libraries enriched/depleted in Functional /
              Partial_loss / Complete_loss?)

  Test 2 — Global Library × Dominant_failure_mode independence
              (chi-square: are specific libraries enriched in
              premature_stop vs ambiguous_aa vs mixed failures?)

  Test 3 — Within-EO Library × Genotype (AAAA vs non-AAAA)
              (Fisher's exact test per EO between Library 009/010 and the
              pooled "other" libraries within that EO)

  Test 4 — Within-EO Library × SI_functional_status
              (Fisher's exact test per EO between Library 009/010 and pooled
              "other" libraries — most directly addresses whether the EO76
              SI-escape cluster is contaminated by a library bias)

Outputs
-------
  SRK_library_effect_tests.tsv   — full test results, one row per test
  SRK_library_effect_summary.pdf — multi-panel diagnostic figure
  stdout                          — human-readable summary with interpretation
"""

from __future__ import annotations

import os
import sys
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


AUDIT_TSV   = "SRK_sample_exclusion_audit.tsv"
ZYG_TSV     = "SRK_individual_zygosity.tsv"
GFS_TSV     = "SRK_individual_GFS.tsv"
BL_TSV      = "SRK_individual_BL_assignments.tsv"

OUT_TSV     = "SRK_library_effect_tests.tsv"
OUT_PDF     = "SRK_library_effect_summary.pdf"

ALPHA = 0.05    # significance threshold for printed callouts


# ─── Load data ────────────────────────────────────────────────────────────────
audit = pd.read_csv(AUDIT_TSV, sep="\t")
audit["Library"] = audit["Library"].astype(str)

# Merge in genotype (from zygosity) and GFS for individuals in the final dataset
zyg = pd.read_csv(ZYG_TSV, sep="\t")[["Individual", "Genotype"]]
gfs = pd.read_csv(GFS_TSV, sep="\t")[["Individual", "GFS"]]
df  = audit.merge(zyg, left_on="Sample_ID", right_on="Individual", how="left")
df  = df.merge(gfs, on="Individual", how="left")

print(f"Loaded {len(df)} samples ({df['Final_status'].eq('Included').sum()} included).\n")


# ─── Test 1: Library × SI_functional_status (global) ─────────────────────────
print("=" * 72)
print("Test 1 — Global Library × SI_functional_status independence")
print("=" * 72)

s9_df = df[df["SI_functional_status"].isin(["Functional", "Partial_loss", "Complete_loss"])]
ct = pd.crosstab(s9_df["Library"], s9_df["SI_functional_status"])
print("\nContingency table (Library × SI_functional_status):")
print(ct)

chi2, pval, dof, expected = stats.chi2_contingency(ct.values)
print(f"\nChi-square = {chi2:.2f}  df = {dof}  p = {pval:.4g}")
test1_sig = pval < ALPHA

# Identify libraries with largest standardized residual contributions
residuals = (ct.values - expected) / np.sqrt(expected)
res_df = pd.DataFrame(residuals, index=ct.index, columns=ct.columns)
print("\nStandardised residuals (|>2| = notable deviation):")
print(res_df.round(2))

test_rows = [{
    "Test": "1_global_SI_status",
    "Comparison": "Library × SI_functional_status (all EOs)",
    "Statistic": "chi2",
    "Value": chi2, "df": dof, "p_value": pval,
    "Significant": test1_sig,
    "Notes": "Standardised residuals saved in PDF figure",
}]


# ─── Test 2: Library × Dominant_failure_mode (global) ────────────────────────
print()
print("=" * 72)
print("Test 2 — Global Library × Dominant_failure_mode independence")
print("=" * 72)

mode_df = df[df["Dominant_failure_mode"].isin(["premature_stop", "ambiguous_aa", "mixed", "other"])]
ct2 = pd.crosstab(mode_df["Library"], mode_df["Dominant_failure_mode"])
print("\nContingency table (Library × Dominant_failure_mode):")
print(ct2)

chi2_2, pval_2, dof_2, expected_2 = stats.chi2_contingency(ct2.values)
print(f"\nChi-square = {chi2_2:.2f}  df = {dof_2}  p = {pval_2:.4g}")
test2_sig = pval_2 < ALPHA

residuals_2 = (ct2.values - expected_2) / np.sqrt(expected_2)
res_df_2 = pd.DataFrame(residuals_2, index=ct2.index, columns=ct2.columns)
print("\nStandardised residuals (|>2| = notable deviation):")
print(res_df_2.round(2))

test_rows.append({
    "Test": "2_global_failure_mode",
    "Comparison": "Library × Dominant_failure_mode (all EOs)",
    "Statistic": "chi2",
    "Value": chi2_2, "df": dof_2, "p_value": pval_2,
    "Significant": test2_sig,
    "Notes": "Expected: Library 009 enriched in ambiguous_aa (no Step 7b)",
})


# ─── Test 3 & 4: Within-EO library effect for AAAA proportion and SI status ──
print()
print("=" * 72)
print("Test 3 & 4 — Within-EO library effect (focus: Library 009 / 010)")
print("=" * 72)

# Test 3 (AAAA proportion) uses only INCLUDED individuals — Complete_loss samples
# are excluded by definition (they have no genotype call), so we test whether the
# AAAA/non-AAAA split differs by library within an EO using the surviving cohort.
inc = df[df["Final_status"] == "Included"].copy()
inc["is_AAAA"] = (inc["Genotype"] == "AAAA").astype(int)

# Test 4 (Complete_loss proportion) uses ALL samples that REACHED Step 9 — i.e.
# Functional + Partial_loss + Complete_loss. This includes excluded samples and
# is the right denominator for asking "does Library X produce more Complete_loss
# *signals* per sample-that-was-translatable than Library Y?"
reach9 = df[df["SI_functional_status"].isin(["Functional", "Partial_loss", "Complete_loss"])].copy()
reach9["is_CompleteLoss"] = (reach9["SI_functional_status"] == "Complete_loss").astype(int)

# For each focus library (009 and 010), test within each multi-library EO
within_results = []
for focus_lib in ["9", "10"]:
    # === Test 3: AAAA proportion among INCLUDED individuals ===
    for eo, grp in inc.groupby("EO"):
        libs_in_eo = grp["Library"].unique().tolist()
        if focus_lib not in libs_in_eo or len(libs_in_eo) < 2:
            continue
        focus_grp = grp[grp["Library"] == focus_lib]
        other_grp = grp[grp["Library"] != focus_lib]
        if len(focus_grp) < 3 or len(other_grp) < 3:
            continue

        a_focus = focus_grp["is_AAAA"].sum()
        n_focus = len(focus_grp)
        a_other = other_grp["is_AAAA"].sum()
        n_other = len(other_grp)
        _, p_aaaa = stats.fisher_exact([[a_focus, n_focus - a_focus],
                                        [a_other, n_other - a_other]])

        # Test 4: Complete_loss proportion among ALL samples reaching Step 9
        eo_reach9 = reach9[reach9["EO"] == eo]
        focus_r9 = eo_reach9[eo_reach9["Library"] == focus_lib]
        other_r9 = eo_reach9[eo_reach9["Library"] != focus_lib]
        if len(focus_r9) < 3 or len(other_r9) < 3:
            n_focus_r9 = len(focus_r9); n_other_r9 = len(other_r9)
            c_focus = focus_r9["is_CompleteLoss"].sum()
            c_other = other_r9["is_CompleteLoss"].sum()
            p_cl = np.nan
        else:
            n_focus_r9 = len(focus_r9); n_other_r9 = len(other_r9)
            c_focus = focus_r9["is_CompleteLoss"].sum()
            c_other = other_r9["is_CompleteLoss"].sum()
            _, p_cl = stats.fisher_exact([[c_focus, n_focus_r9 - c_focus],
                                          [c_other, n_other_r9 - c_other]])

        within_results.append({
            "EO":               eo,
            "Focus_library":    f"Library {focus_lib}",
            "n_focus_incl":     n_focus,
            "n_other_incl":     n_other,
            "prop_AAAA_focus":  a_focus / n_focus,
            "prop_AAAA_other":  a_other / n_other,
            "p_AAAA_diff":      p_aaaa,
            "n_focus_reach9":   n_focus_r9,
            "n_other_reach9":   n_other_r9,
            "prop_CL_focus":    (c_focus / n_focus_r9) if n_focus_r9 else np.nan,
            "prop_CL_other":    (c_other / n_other_r9) if n_other_r9 else np.nan,
            "p_CL_diff":        p_cl,
        })

within_df = pd.DataFrame(within_results)
print()
print("Within-EO comparisons (Library 009 / 010 vs other libraries in same EO):")
if not within_df.empty:
    print(within_df.round(3).to_string(index=False))

    for row in within_df.itertuples():
        if row.p_AAAA_diff < ALPHA:
            test_rows.append({
                "Test":        "3_within_EO_AAAA",
                "Comparison":  f"{row.Focus_library} vs others in EO{row.EO}, AAAA proportion",
                "Statistic":   "fisher_exact",
                "Value":       np.nan, "df": np.nan, "p_value": row.p_AAAA_diff,
                "Significant": True,
                "Notes":       f"{row.Focus_library}: {row.prop_AAAA_focus:.2f}; others: {row.prop_AAAA_other:.2f}",
            })
        if pd.notna(row.p_CL_diff) and row.p_CL_diff < ALPHA:
            test_rows.append({
                "Test":        "4_within_EO_CompleteLoss",
                "Comparison":  f"{row.Focus_library} vs others in EO{row.EO}, Complete_loss proportion (Step-9-reached cohort)",
                "Statistic":   "fisher_exact",
                "Value":       np.nan, "df": np.nan, "p_value": row.p_CL_diff,
                "Significant": True,
                "Notes":       f"{row.Focus_library}: {row.prop_CL_focus:.2f}; others: {row.prop_CL_other:.2f}",
            })
else:
    print("  (no multi-library EOs with sufficient sample size for testing)")

print()
print(f"  Tests run: {len(within_df) * 2}  Significant: "
      f"{(within_df['p_AAAA_diff'] < ALPHA).sum() + (within_df['p_CL_diff'] < ALPHA).sum()}")


# ─── Write tests TSV ──────────────────────────────────────────────────────────
out = pd.DataFrame(test_rows)
out.to_csv(OUT_TSV, sep="\t", index=False)
print()
print(f"Test results table written to: {OUT_TSV}")


# ─── Diagnostic figure ────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: Heatmap of standardised residuals — Library × SI_status
im = axes[0, 0].imshow(res_df.values, cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto")
axes[0, 0].set_xticks(range(len(res_df.columns)))
axes[0, 0].set_xticklabels(res_df.columns, rotation=30, ha="right")
axes[0, 0].set_yticks(range(len(res_df.index)))
axes[0, 0].set_yticklabels([f"Library {i}" for i in res_df.index])
axes[0, 0].set_title(f"Test 1 — Library × SI_functional_status\n"
                     f"χ² = {chi2:.1f}, df = {dof}, p = {pval:.3g}")
for i in range(len(res_df.index)):
    for j in range(len(res_df.columns)):
        v = res_df.values[i, j]
        if abs(v) > 1.5:
            axes[0, 0].text(j, i, f"{v:.1f}", ha="center", va="center",
                            color="white" if abs(v) > 2.5 else "black", fontsize=9)
plt.colorbar(im, ax=axes[0, 0], label="Std residual")

# Panel B: Heatmap of residuals — Library × Dominant_failure_mode
im2 = axes[0, 1].imshow(res_df_2.values, cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto")
axes[0, 1].set_xticks(range(len(res_df_2.columns)))
axes[0, 1].set_xticklabels(res_df_2.columns, rotation=30, ha="right")
axes[0, 1].set_yticks(range(len(res_df_2.index)))
axes[0, 1].set_yticklabels([f"Library {i}" for i in res_df_2.index])
axes[0, 1].set_title(f"Test 2 — Library × Dominant_failure_mode\n"
                     f"χ² = {chi2_2:.1f}, df = {dof_2}, p = {pval_2:.3g}")
for i in range(len(res_df_2.index)):
    for j in range(len(res_df_2.columns)):
        v = res_df_2.values[i, j]
        if abs(v) > 1.5:
            axes[0, 1].text(j, i, f"{v:.1f}", ha="center", va="center",
                            color="white" if abs(v) > 2.5 else "black", fontsize=9)
plt.colorbar(im2, ax=axes[0, 1], label="Std residual")

# Panel C: Within-EO AAAA proportion (Library 9/10 vs others)
if not within_df.empty:
    bars = []
    labels = []
    colors = []
    for row in within_df.itertuples():
        labels.append(f"EO{row.EO}\n({row.Focus_library})")
        bars.append((row.prop_AAAA_focus, row.prop_AAAA_other))
        colors.append("crimson" if row.p_AAAA_diff < ALPHA else "steelblue")

    x = np.arange(len(bars))
    w = 0.4
    focus_vals = [b[0] for b in bars]
    other_vals = [b[1] for b in bars]
    axes[1, 0].bar(x - w/2, focus_vals, width=w, label="Focus library", color=colors, edgecolor="black")
    axes[1, 0].bar(x + w/2, other_vals, width=w, label="Other libraries in same EO",
                   color="lightgrey", edgecolor="black")
    axes[1, 0].set_xticks(x)
    axes[1, 0].set_xticklabels(labels, fontsize=8, rotation=0)
    axes[1, 0].set_ylabel("Proportion AAAA")
    axes[1, 0].set_title("Test 3 — Within-EO AAAA proportion\n(red bar = significant difference)")
    axes[1, 0].legend(loc="best", fontsize=8)
    axes[1, 0].set_ylim(0, 1)
else:
    axes[1, 0].axis("off")

# Panel D: Within-EO Complete_loss proportion
if not within_df.empty:
    bars = []
    labels = []
    colors = []
    for row in within_df.itertuples():
        labels.append(f"EO{row.EO}\n({row.Focus_library})")
        bars.append((row.prop_CL_focus, row.prop_CL_other))
        colors.append("crimson" if row.p_CL_diff < ALPHA else "steelblue")

    x = np.arange(len(bars))
    w = 0.4
    focus_vals = [b[0] for b in bars]
    other_vals = [b[1] for b in bars]
    axes[1, 1].bar(x - w/2, focus_vals, width=w, label="Focus library", color=colors, edgecolor="black")
    axes[1, 1].bar(x + w/2, other_vals, width=w, label="Other libraries in same EO",
                   color="lightgrey", edgecolor="black")
    axes[1, 1].set_xticks(x)
    axes[1, 1].set_xticklabels(labels, fontsize=8, rotation=0)
    axes[1, 1].set_ylabel("Proportion Complete_loss")
    axes[1, 1].set_title("Test 4 — Within-EO Complete_loss proportion\n(red bar = significant difference)")
    axes[1, 1].legend(loc="best", fontsize=8)
    axes[1, 1].set_ylim(0, max(0.3, axes[1, 1].get_ylim()[1]))
else:
    axes[1, 1].axis("off")

fig.tight_layout()
fig.savefig(OUT_PDF, dpi=200)
print(f"Diagnostic figure written to: {OUT_PDF}")


# ─── Interpretation summary ───────────────────────────────────────────────────
print()
print("=" * 72)
print("INTERPRETATION")
print("=" * 72)

if test1_sig:
    print("• Test 1 (SI_functional_status): LIBRARIES DIFFER significantly. Examine "
          "the residuals table to identify which library(ies) deviate. A single "
          "library with a strong positive residual in 'Complete_loss' would be a "
          "red flag for technical artefact contaminating the SI-escape signal.")
else:
    print("• Test 1 (SI_functional_status): no significant global library effect.")

if test2_sig:
    print("• Test 2 (Dominant_failure_mode): LIBRARIES DIFFER. This is expected "
          "for Library 009 (no Step 7b processing → elevated ambiguous_aa); but a "
          "library showing elevated 'premature_stop' beyond expectation would "
          "indicate a true bias in the LoF signal.")
else:
    print("• Test 2 (Dominant_failure_mode): no significant library effect.")

if not within_df.empty:
    sig_aaaa = (within_df["p_AAAA_diff"] < ALPHA).sum()
    sig_cl   = (within_df["p_CL_diff"] < ALPHA).sum()
    print(f"• Tests 3+4 (within-EO Library 009/010 vs others): {sig_aaaa} significant "
          f"AAAA-proportion differences, {sig_cl} significant Complete_loss differences "
          f"out of {len(within_df)} comparisons. If 0 → no within-EO bias.")

    cl_focus_higher = within_df[within_df["prop_CL_focus"] > within_df["prop_CL_other"]]
    if not cl_focus_higher.empty:
        print(f"  Within-EO comparisons where focus library has HIGHER Complete_loss:")
        for row in cl_focus_higher.itertuples():
            flag = " ← SIGNIFICANT" if row.p_CL_diff < ALPHA else ""
            print(f"    EO{row.EO} ({row.Focus_library}): "
                  f"{row.prop_CL_focus:.2f} vs {row.prop_CL_other:.2f}  "
                  f"p = {row.p_CL_diff:.3f}{flag}")
