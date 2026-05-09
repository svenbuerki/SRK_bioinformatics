#!/usr/bin/env python3
"""
srk_wgroup_collapse_test.py   (Step 22c — diagnostic, optional)

Synonymy group collapse diagnostic for the cross-Brassicaceae variability landscape.

Hypothesis tested
-----------------
LEPA's per-column Shannon entropy is markedly lower than Brassica's or
Arabidopsis's despite a comparable HV-column count. One mechanism is that
many LEPA "alleles" are HV-identical (W-groups), so they contribute the same
residue at every HV column and inflate the column-variability count without
adding per-column diversity.

This script tests that hypothesis directly: it re-runs the per-species
Shannon-entropy scan on a deduplicated LEPA set in which each Synonymy group is
represented by a single allele (the one with the most AAAA individuals;
ties broken by lowest allele number for reproducibility) plus all isolated
alleles. If the LEPA entropy profile rises substantially toward the Brassica
profile after collapsing, Synonymy group redundancy is the dominant explanation.

Pipeline placement
------------------
This script must be run AFTER Steps 22a and 22b:

  python3 srk_variability_landscape.py     # Step 22a — produces SRK_LEPA_HV_positions.tsv
  python3 srk_allele_hypotheses.py         # Step 22b — produces SRK_synonymy_groups.csv
  python3 srk_wgroup_collapse_test.py      # Step 22c — this diagnostic

It reads only existing outputs and produces only new files; it never
overwrites Step 22a's authoritative landscape figure or HV-position TSV.

Inputs
------
  SRK_synonymy_groups.csv         from Step 22b (Synonymy group membership)
  SRK_combined_alignment.fasta    from Step 22a (LEPA + Brassica + Arabidopsis)
  SRK_LEPA_HV_positions.tsv       from Step 22a (canonical LEPA HV columns)

Outputs
-------
  SRK_LEPA_synonymy_group_representatives.tsv               which alleles were retained / dropped
  SRK_variability_landscape_wgroup_collapsed.pdf    side-by-side figure (PDF)
  figures/SRK_variability_landscape_wgroup_collapsed.png
  SRK_wgroup_collapse_entropy_summary.tsv           numeric before/after summary
"""

from __future__ import annotations
import os
import sys
import math
import csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
from Bio import SeqIO

# =============================================================================
# Settings — must match srk_variability_landscape.py for apples-to-apples test
# =============================================================================
COMBINED_ALN     = "SRK_combined_alignment.fasta"
SYN_GROUPS_CSV   = "SRK_synonymy_groups.csv"
LEPA_HV_TSV      = "SRK_LEPA_HV_positions.tsv"
DOMAIN_REGION    = (31, 430)
WINDOW_SIZE      = 20
PEAK_SD_FACTOR   = 1.0
MIN_HV_RUN       = 3

# Species classification (same logic as srk_variability_landscape.py)
SPECIES_KEYS = [
    ("LEPA",        "Allele_"),
    ("Brassica",    "Brassica_"),
    ("Brassica",    "BrSRK9"),
    ("Arabidopsis", "Arabidopsis_"),
]
SPECIES_COLOURS = {
    "LEPA":         "#1f77b4",
    "LEPA_collapsed": "#9467bd",
    "Brassica":     "#d62728",
    "Arabidopsis":  "#2ca02c",
}

OUT_REP_TSV   = "SRK_LEPA_synonymy_group_representatives.tsv"
OUT_PDF       = "SRK_variability_landscape_wgroup_collapsed.pdf"
OUT_PNG       = os.path.join("figures", "SRK_variability_landscape_wgroup_collapsed.png")
OUT_SUMMARY   = "SRK_wgroup_collapse_entropy_summary.tsv"

os.makedirs("figures", exist_ok=True)

# =============================================================================
# 1. Sanity-check inputs and choose Synonymy group representatives
# =============================================================================
for f in (COMBINED_ALN, SYN_GROUPS_CSV, LEPA_HV_TSV):
    if not os.path.exists(f):
        sys.exit(f"ERROR: required input missing: {f}\n"
                 f"Run Steps 22a and 22b first (see script docstring).")

print(f"Reading Synonymy group memberships from {SYN_GROUPS_CSV} ...")
syn = pd.read_csv(SYN_GROUPS_CSV)
print(f"  {len(syn)} alleles, {syn['Synonymy_group'].nunique()} groups "
      f"(incl. 'Isolated' singletons)")

# For each W-group, pick the allele with the most AAAA individuals
# (ties broken by lowest allele ID for reproducibility). Isolated alleles
# are kept as-is (each is its own representative).
def pick_representative(grp: pd.DataFrame) -> str:
    grp = grp.sort_values(["AAAA_count", "Allele"],
                          ascending=[False, True])
    return grp.iloc[0]["Allele"]

representatives = []
dropped_log = []
for gname, grp in syn.groupby("Synonymy_group", sort=False):
    if gname == "Isolated":
        # Each isolated allele is its own representative
        for a in grp["Allele"]:
            representatives.append({"Allele": a, "Group": "Isolated",
                                    "Group_size": 1, "Role": "isolated"})
        continue
    rep = pick_representative(grp)
    n = len(grp)
    for a in grp["Allele"]:
        if a == rep:
            representatives.append({"Allele": a, "Group": gname,
                                    "Group_size": n, "Role": "representative"})
        else:
            dropped_log.append({"Allele": a, "Group": gname,
                                "Group_size": n, "Role": "redundant",
                                "Representative": rep})

# Build a clean output table
rows = []
for r in representatives:
    rows.append({**r, "Representative": r["Allele"]})
for r in dropped_log:
    rows.append(r)
rep_df = pd.DataFrame(rows).sort_values("Allele").reset_index(drop=True)
rep_df.to_csv(OUT_REP_TSV, sep="\t", index=False)
n_kept = (rep_df["Role"].isin(["representative", "isolated"])).sum()
n_drop = (rep_df["Role"] == "redundant").sum()
print(f"  kept {n_kept} allele(s) ({n_kept - (rep_df['Role']=='isolated').sum()} "
      f"Synonymy group representatives + {(rep_df['Role']=='isolated').sum()} isolated)")
print(f"  dropped {n_drop} Synonymy group redundant allele(s)")
print(f"  Written {OUT_REP_TSV}")

kept_set = set(rep_df.loc[rep_df["Role"].isin(["representative", "isolated"]),
                          "Allele"])

# =============================================================================
# 2. Load combined alignment, build per-species sequence buckets
# =============================================================================
print(f"\nLoading combined alignment from {COMBINED_ALN} ...")
records = list(SeqIO.parse(COMBINED_ALN, "fasta"))
aln_len = len(records[0].seq)


def classify(rec_id: str) -> str | None:
    for species, key in SPECIES_KEYS:
        if key in rec_id:
            return species
    return None


lepa_full      = []
lepa_collapsed = []
brassica       = []
arabidopsis    = []
for r in records:
    sp = classify(r.id)
    s = str(r.seq).upper()
    if sp == "LEPA":
        lepa_full.append(s)
        if r.id in kept_set or r.id.split(".")[0] in kept_set:
            lepa_collapsed.append(s)
    elif sp == "Brassica":
        brassica.append(s)
    elif sp == "Arabidopsis":
        arabidopsis.append(s)

# Sanity-check — kept alleles should match
print(f"  LEPA full       : {len(lepa_full)} sequences")
print(f"  LEPA collapsed  : {len(lepa_collapsed)} sequences "
      f"(target {len(kept_set)})")
print(f"  Brassica        : {len(brassica)} sequences")
print(f"  Arabidopsis     : {len(arabidopsis)} sequences")
if len(lepa_collapsed) != len(kept_set):
    print(f"  WARNING: LEPA collapsed set size differs from kept allele list")

# =============================================================================
# 3. Identify LEPA-original alignment columns and S-domain restriction
# =============================================================================
lepa_orig_cols = [c for c in range(aln_len)
                  if any(seq[c] != "-" for seq in lepa_full)]
exp_to_lepa1 = {}
i = 0
for c in range(aln_len):
    if c in set(lepa_orig_cols):
        i += 1
        exp_to_lepa1[c] = i
n_lepa_orig = i

lo1, hi1 = DOMAIN_REGION
sd_cols_exp = [c for c in lepa_orig_cols
               if lo1 <= exp_to_lepa1[c] <= hi1]
plot_x = np.array([exp_to_lepa1[c] for c in sd_cols_exp])
print(f"\n  LEPA-original columns: {n_lepa_orig}")
print(f"  S-domain (cols {lo1}-{hi1}) covers {len(sd_cols_exp)} positions")

# =============================================================================
# 4. Per-species Shannon entropy
# =============================================================================
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def shannon(col: list[str]) -> float:
    counts: dict[str, int] = {}
    n = 0
    for r in col:
        if r in AA_ALPHABET:
            counts[r] = counts.get(r, 0) + 1
            n += 1
    if n == 0:
        return 0.0
    H = 0.0
    for k in counts.values():
        p = k / n
        H -= p * math.log2(p)
    return H


def entropy_profile(seqs: list[str], cols: list[int]) -> np.ndarray:
    H = np.zeros(len(cols))
    for i, c in enumerate(cols):
        H[i] = shannon([s[c] for s in seqs])
    return H


print("\nComputing per-species Shannon entropy ...")
profiles: dict[str, dict] = {}
for sp, seqs in [
    ("LEPA",           lepa_full),
    ("LEPA_collapsed", lepa_collapsed),
    ("Brassica",       brassica),
    ("Arabidopsis",    arabidopsis),
]:
    if len(seqs) < 2:
        continue
    H        = entropy_profile(seqs, sd_cols_exp)
    H_smooth = uniform_filter1d(H, size=WINDOW_SIZE)
    thr      = H_smooth.mean() + PEAK_SD_FACTOR * H_smooth.std()
    profiles[sp] = {"n": len(seqs), "H": H, "H_smooth": H_smooth, "thr": thr}
    print(f"  {sp:16s}: n={len(seqs):3d}  median H = {np.median(H):.3f} bits  "
          f"max H = {H.max():.3f}  threshold = {thr:.3f}")

# =============================================================================
# 5. Identify HV regions per species (same MIN_HV_RUN as Step 22a)
# =============================================================================
def find_runs(profile: np.ndarray, x_lepa: np.ndarray,
              thr: float, min_run: int):
    runs = []
    mask = profile > thr
    i, n = 0, len(mask)
    while i < n:
        if mask[i]:
            j = i
            while j + 1 < n and mask[j + 1]:
                j += 1
            if j - i + 1 >= min_run:
                runs.append((int(x_lepa[i]), int(x_lepa[j]), j - i + 1))
            i = j + 1
        else:
            i += 1
    return runs


hv_regions = {}
for sp, prof in profiles.items():
    runs = find_runs(prof["H_smooth"], plot_x, prof["thr"], MIN_HV_RUN)
    hv_regions[sp] = runs

# =============================================================================
# 6. Numeric summary table
# =============================================================================
summary_rows = []
for sp, prof in profiles.items():
    runs = hv_regions[sp]
    summary_rows.append({
        "Species":        sp,
        "n_sequences":    prof["n"],
        "median_H_bits":  f"{np.median(prof['H']):.4f}",
        "mean_H_bits":    f"{prof['H'].mean():.4f}",
        "max_H_bits":     f"{prof['H'].max():.4f}",
        "threshold_bits": f"{prof['thr']:.4f}",
        "n_HV_runs":      len(runs),
        "n_HV_cols":      sum(r[2] for r in runs),
    })
summ = pd.DataFrame(summary_rows)
summ.to_csv(OUT_SUMMARY, sep="\t", index=False)

print(f"\n-- Numeric summary ({OUT_SUMMARY}) --")
print(summ.to_string(index=False))

# Headline: how much does collapsing change LEPA?
if "LEPA" in profiles and "LEPA_collapsed" in profiles:
    lepa_med  = np.median(profiles["LEPA"]["H"])
    lepac_med = np.median(profiles["LEPA_collapsed"]["H"])
    bra_med   = np.median(profiles["Brassica"]["H"]) if "Brassica" in profiles else float("nan")
    delta = lepac_med - lepa_med
    print(f"\n-- Headline --")
    print(f"  LEPA full       median H = {lepa_med:.3f} bits  "
          f"(n={len(lepa_full)})")
    print(f"  LEPA collapsed  median H = {lepac_med:.3f} bits  "
          f"(n={len(lepa_collapsed)})  Δ = {delta:+.3f} bits")
    print(f"  Brassica        median H = {bra_med:.3f} bits  "
          f"(n={len(brassica)})")
    if abs(delta) < 0.05:
        verdict = "negligible — Synonymy group redundancy does NOT explain LEPA's low entropy"
    elif delta > 0 and lepac_med >= 0.5 * bra_med:
        verdict = ("substantial — Synonymy group collapse moves LEPA toward Brassica-like "
                   "entropy → redundancy IS a dominant mechanism")
    elif delta > 0:
        verdict = ("partial — Synonymy group collapse raises LEPA entropy but a gap with "
                   "Brassica remains; redundancy contributes but other mechanisms "
                   "(shallow phylogeny, drift purging of rare residues) also active")
    else:
        verdict = "decreased — unexpected, requires further investigation"
    print(f"  Verdict: {verdict}")

# =============================================================================
# 7. Side-by-side variability-landscape figure
# =============================================================================
print(f"\nWriting figure to {OUT_PDF} ...")

fig, ax = plt.subplots(figsize=(17, 5.6))
plot_order = [("Arabidopsis", "Arabidopsis"),
              ("Brassica",    "Brassica"),
              ("LEPA",        "LEPA"),
              ("LEPA_collapsed", "LEPA (Synonymy group collapsed)")]
for sp, label in plot_order:
    if sp not in profiles:
        continue
    prof = profiles[sp]
    col = SPECIES_COLOURS[sp]
    ax.fill_between(plot_x, prof["H"], alpha=0.07, color=col)
    ax.plot(plot_x, prof["H_smooth"], color=col, linewidth=1.6,
            label=f"{label} (n={prof['n']}, threshold = {prof['thr']:.2f})")
    ax.axhline(prof["thr"], color=col, linestyle=":",
               linewidth=0.8, alpha=0.5)

ax.set_xlim(plot_x[0], plot_x[-1])
ax.set_xlabel("LEPA alignment position (aa, 1-based; S-domain)", fontsize=11)
ax.set_ylabel("Smoothed Shannon entropy (bits)", fontsize=11)
ax.set_title(
    "Synonymy group collapse diagnostic — does Synonymy group redundancy explain LEPA's "
    "low Shannon entropy?\n"
    f"LEPA full ({len(lepa_full)} alleles) vs LEPA collapsed "
    f"({len(lepa_collapsed)} Synonymy group representatives + isolated) vs Brassica vs Arabidopsis",
    fontsize=11.5,
)
ax.legend(fontsize=9, loc="upper right", framealpha=0.95, edgecolor="grey")

# Annotate headline
if "LEPA" in profiles and "LEPA_collapsed" in profiles:
    ax.text(
        0.005, 0.97,
        f"LEPA full        median H = {lepa_med:.3f} bits  (n={len(lepa_full)})\n"
        f"LEPA collapsed median H = {lepac_med:.3f} bits  "
        f"(n={len(lepa_collapsed)}; Δ = {delta:+.3f})\n"
        f"Brassica         median H = {bra_med:.3f} bits  (n={len(brassica)})",
        transform=ax.transAxes, ha="left", va="top", fontsize=9,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="grey", lw=0.6))

plt.tight_layout()
plt.savefig(OUT_PDF, format="pdf", bbox_inches="tight")
plt.savefig(OUT_PNG, format="png", dpi=200, bbox_inches="tight")
plt.close()
print(f"  Saved {OUT_PDF}")
print(f"  Saved {OUT_PNG}")

print("\nDone.")
