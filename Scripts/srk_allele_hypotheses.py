#!/usr/bin/env python3
"""
srk_allele_hypotheses.py  —  Step 21 (revised)

Data-driven allele hypothesis testing and crossing design for Objective 2.

Replaces the arbitrary N_SUPERGROUPS parameter in the previous approach with
a moving-window variability scan that identifies which positions in the S-domain
alignment actually discriminate alleles (hypervariable, HV positions), then
clusters allele representatives using only those positions.

Cross categories — biologically motivated:
  W — same original allele bin (Step 10 definition)
      Expected: incompatible (same SI specificity by sequence identity)
  N — different original bin, same HV functional group
      Expected: synonymy test — alleles identical in HV residues but split
      by the full-domain clustering; do they share recognition specificity?
  P — different HV functional group
      Expected: compatible (positive control — functionally distinct alleles)

Workflow
--------
  1  Load aligned allele representative sequences (Step 10 output)
  2  Moving-window per-column variability scan → identify HV positions
  3  Compute HV-only pairwise distance matrix
  4  UPGMA clustering in HV space → functional groups (auto threshold)
  5  Compare with original Step 10 bins → flag synonymy candidates
  6  Load AAAA individuals (Step 12 zygosity output)
  7  Generate W / N / P cross design matrix
  8  Output figures and TSVs
  9  (Optional) Cross result analysis when data available

Outputs
-------
  SRK_variability_landscape.pdf     per-column and smoothed variability profile
  SRK_HV_allele_distances.tsv       HV-only pairwise distance matrix
  SRK_functional_allele_groups.tsv  allele bin → functional group mapping
  SRK_synonymy_candidates.tsv       pairs in same HV group but different original bins
  SRK_AAAA_cross_design_HV.tsv      AAAA × AAAA cross design (W/N/P)
  SRK_HV_cluster_figure.pdf         UPGMA dendrogram (HV distances) + AAAA availability
  SRK_cross_result_analysis_HV.pdf  seed yield by category (requires cross data)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from scipy.ndimage import uniform_filter1d
from Bio import SeqIO

# ─────────────────────────────────────────────────────────────────────────────
# Settings — adjust before running
# ─────────────────────────────────────────────────────────────────────────────

# Input files
REPS_FASTA   = "SRK_protein_allele_representatives.fasta"
ZYGO_TSV     = "SRK_individual_zygosity.tsv"
ALLELE_TABLE = "SRK_individual_allele_table.tsv"
CROSS_TSV    = None   # set to cross results TSV to activate Part 4

# S-domain region (1-based, inclusive) — must match Step 10
DOMAIN_REGION = (31, 430)

# Moving-window variability scan
WINDOW_SIZE      = 20    # sliding window width (aa) for smoothing the profile
PEAK_SD_FACTOR   = 1.0   # positions above mean + k*SD of smoothed profile → HV
MIN_HV_POSITIONS = 10    # warn if fewer HV positions found; lower PEAK_SD_FACTOR

# Functional group clustering on HV-only distances (identifies phylogenetic classes)
# Priority: N_GROUPS > DISTANCE_THRESHOLD > auto (largest gap in linkage heights)
# The auto threshold is recommended on first run — it detects the Class I / Class II split.
N_GROUPS           = None   # fixed number of groups; set to integer to override auto
DISTANCE_THRESHOLD = None   # HV p-distance cutoff; used when N_GROUPS is None

# Within-class W / N / P thresholds (applied after class split is identified)
# W: HV distance = 0  (HV-identical alleles — strongest incompatibility prediction)
# N: 0 < HV distance < WITHIN_CLASS_THRESHOLD  (small HV differences — synonymy test)
# P_within: HV distance >= WITHIN_CLASS_THRESHOLD, same class  (distinct within-class alleles)
# P_between: different class  (between-class positive control — functionally guaranteed distinct)
WITHIN_CLASS_THRESHOLD = 0.04   # adjust based on within-class HV distance distribution

# Cross category labels
CAT_W        = "W"          # HV-identical — expected incompatible
CAT_N        = "N"          # small HV differences — synonymy test
CAT_P_WITHIN = "P_within"   # same class, substantial HV differences — predicted compatible
CAT_P_CROSS  = "P_cross"    # different class — guaranteed positive control

# Output files
OUT_VARIABILITY  = "SRK_variability_landscape.pdf"
OUT_HV_DISTS     = "SRK_HV_allele_distances.tsv"
OUT_FUNC_GROUPS  = "SRK_functional_allele_groups.tsv"
OUT_SYNONYMY     = "SRK_synonymy_candidates.tsv"
OUT_CROSS_DESIGN = "SRK_AAAA_cross_design_HV.tsv"
OUT_CLUSTER_FIG  = "SRK_HV_cluster_figure.pdf"
OUT_RESULTS_FIG  = "SRK_cross_result_analysis_HV.pdf"

os.makedirs("figures", exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# PART 1 — Variability landscape and HV position identification
# ─────────────────────────────────────────────────────────────────────────────

print("=" * 60)
print("PART 1 — Variability landscape and HV position identification")
print("=" * 60)

# Load aligned representative sequences
records     = list(SeqIO.parse(REPS_FASTA, "fasta"))
allele_names = [r.id.split()[0] for r in records]
seqs         = [str(r.seq).upper() for r in records]
n_alleles    = len(allele_names)
aln_len      = len(seqs[0])

print(f"\nLoaded {n_alleles} allele representative sequences ({aln_len} aligned positions)")

# Restrict to S-domain
start_0  = max(0, DOMAIN_REGION[0] - 1)
end_0    = min(aln_len, DOMAIN_REGION[1])
seqs_sd  = [s[start_0:end_0] for s in seqs]
sd_len   = end_0 - start_0
print(f"S-domain: columns {DOMAIN_REGION[0]}–{DOMAIN_REGION[1]} ({sd_len} positions)")

# Per-column mean pairwise distance
print("Computing per-column variability ...")
col_var = np.zeros(sd_len)
for col in range(sd_len):
    residues = [s[col] for s in seqs_sd]
    diffs = valid = 0
    for i in range(n_alleles):
        for j in range(i + 1, n_alleles):
            a, b = residues[i], residues[j]
            if a == "-" or b == "-":
                continue
            valid += 1
            if a != b:
                diffs += 1
    col_var[col] = diffs / valid if valid > 0 else 0.0

# Smooth with uniform moving window
col_var_smooth = uniform_filter1d(col_var, size=WINDOW_SIZE)

# Identify HV positions from the smoothed profile
hv_mean      = col_var_smooth.mean()
hv_sd        = col_var_smooth.std()
hv_threshold = hv_mean + PEAK_SD_FACTOR * hv_sd
hv_mask      = col_var_smooth > hv_threshold
hv_positions = np.where(hv_mask)[0]   # 0-based within S-domain

print(f"\nVariability profile (smoothed, window = {WINDOW_SIZE} aa):")
print(f"  Mean per-column distance  : {hv_mean:.4f}")
print(f"  SD                        : {hv_sd:.4f}")
print(f"  HV threshold (mean+{PEAK_SD_FACTOR}×SD): {hv_threshold:.4f}")
print(f"  HV positions identified   : {len(hv_positions)} / {sd_len}")

if len(hv_positions) < MIN_HV_POSITIONS:
    print(f"\nWARNING: Only {len(hv_positions)} HV positions found (minimum {MIN_HV_POSITIONS}).")
    print(f"Consider lowering PEAK_SD_FACTOR (currently {PEAK_SD_FACTOR}).")

# Summarise HV positions as consecutive runs (1-based alignment coords)
hv_runs = []
if len(hv_positions) > 0:
    run_start = hv_positions[0]
    prev = hv_positions[0]
    for p in hv_positions[1:]:
        if p > prev + 1:
            hv_runs.append((int(run_start + start_0 + 1), int(prev + start_0 + 1)))
            run_start = p
        prev = p
    hv_runs.append((int(run_start + start_0 + 1), int(prev + start_0 + 1)))

    print("\nHV regions (1-based alignment positions):")
    for r in hv_runs:
        width = r[1] - r[0] + 1
        print(f"  {r[0]}–{r[1]}  ({width} aa)")

# ── Variability landscape figure ─────────────────────────────────────────────

print(f"\nWriting variability landscape to {OUT_VARIABILITY} ...")

positions = np.arange(sd_len) + start_0 + 1   # 1-based

fig, ax = plt.subplots(figsize=(16, 4))
ax.fill_between(positions, col_var, alpha=0.2, color="steelblue",
                label="Per-column mean distance")
ax.plot(positions, col_var_smooth, color="steelblue", linewidth=1.2,
        label=f"Smoothed (window = {WINDOW_SIZE} aa)")
ax.axhline(hv_threshold, color="crimson", linestyle="--", linewidth=1.0,
           label=f"HV threshold (mean + {PEAK_SD_FACTOR}×SD = {hv_threshold:.3f})")
for r in hv_runs:
    ax.axvspan(r[0], r[1], alpha=0.12, color="crimson", zorder=0)
ax.set_xlabel("Alignment position (aa, 1-based)", fontsize=11)
ax.set_ylabel("Mean pairwise distance", fontsize=11)
ax.set_title(
    f"S-domain variability landscape — {len(hv_positions)} HV positions identified",
    fontsize=12)
ax.legend(fontsize=9, loc="upper right")
ax.set_xlim(positions[0], positions[-1])
plt.tight_layout()
plt.savefig(OUT_VARIABILITY, format="pdf", dpi=150, bbox_inches="tight")
plt.savefig(os.path.join("figures", OUT_VARIABILITY.replace(".pdf", ".png")),
            format="png", dpi=200, bbox_inches="tight")
plt.close()
print(f"Saved {OUT_VARIABILITY}")

# ─────────────────────────────────────────────────────────────────────────────
# PART 2 — HV-only distance matrix and UPGMA functional group clustering
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("PART 2 — HV-only distances and functional group clustering")
print("=" * 60)

if len(hv_positions) == 0:
    print("ERROR: No HV positions found. Lower PEAK_SD_FACTOR and re-run.")
    sys.exit(1)

# Extract HV columns from S-domain sequences
seqs_hv = ["".join(s[p] for p in hv_positions) for s in seqs_sd]


def p_distance(s1, s2):
    """Fraction of differing non-gap positions."""
    diff = valid = 0
    for a, b in zip(s1, s2):
        if a == "-" or b == "-":
            continue
        valid += 1
        if a != b:
            diff += 1
    return diff / valid if valid > 0 else np.nan


print(f"\nComputing pairwise distances on {len(hv_positions)} HV positions ...")
dist_hv = np.zeros((n_alleles, n_alleles))
for i in range(n_alleles):
    for j in range(i + 1, n_alleles):
        d = p_distance(seqs_hv[i], seqs_hv[j])
        dist_hv[i, j] = d
        dist_hv[j, i] = d

# Replace any NaN with column mean (safeguard)
col_means = np.nanmean(dist_hv, axis=0)
for i in range(n_alleles):
    for j in range(n_alleles):
        if np.isnan(dist_hv[i, j]):
            dist_hv[i, j] = (col_means[i] + col_means[j]) / 2

upper_hv = dist_hv[np.triu_indices(n_alleles, k=1)]
print(f"\nHV-only p-distance summary:")
print(f"  Min    : {upper_hv.min():.4f}")
print(f"  Median : {np.median(upper_hv):.4f}")
print(f"  Mean   : {upper_hv.mean():.4f}")
print(f"  Max    : {upper_hv.max():.4f}")
print(f"  Pairs with dist = 0 (HV-identical) : {(upper_hv == 0).sum()}")

# Write HV distance matrix
pd.DataFrame(dist_hv, index=allele_names, columns=allele_names).to_csv(
    OUT_HV_DISTS, sep="\t")
print(f"\nHV distance matrix written to {OUT_HV_DISTS}")

# ── UPGMA clustering ──────────────────────────────────────────────────────────

condensed = squareform(dist_hv)
Z = linkage(condensed, method="average")   # UPGMA
heights = Z[:, 2]

if N_GROUPS is not None:
    labels   = fcluster(Z, N_GROUPS, criterion="maxclust")
    n_groups = len(set(labels))
    mode_str = f"N_GROUPS = {N_GROUPS} → {n_groups} groups"

elif DISTANCE_THRESHOLD is not None:
    labels   = fcluster(Z, DISTANCE_THRESHOLD, criterion="distance")
    n_groups = len(set(labels))
    mode_str = f"DISTANCE_THRESHOLD = {DISTANCE_THRESHOLD:.4f} → {n_groups} groups"

else:
    # Auto: cut at the midpoint of the largest gap in merge heights
    sorted_h  = np.sort(heights)
    gaps      = np.diff(sorted_h)
    gap_idx   = np.argmax(gaps)
    auto_thr  = (sorted_h[gap_idx] + sorted_h[gap_idx + 1]) / 2
    labels    = fcluster(Z, auto_thr, criterion="distance")
    n_groups  = len(set(labels))
    mode_str  = (f"Auto threshold = {auto_thr:.4f} "
                 f"(gap: {sorted_h[gap_idx]:.4f} → {sorted_h[gap_idx+1]:.4f}) "
                 f"→ {n_groups} functional groups")
    print(f"\nAuto-detected HV distance threshold: {auto_thr:.4f}")
    print("  Largest gap in UPGMA merge heights used as natural cut point.")
    print(f"  Top-5 merge height gaps:")
    top5_idx = np.argsort(gaps)[::-1][:5]
    for idx in top5_idx:
        print(f"    {sorted_h[idx]:.4f} → {sorted_h[idx+1]:.4f}  (gap {gaps[idx]:.4f})")

print(f"\nFunctional group clustering: {mode_str}")

allele_to_fg = {allele_names[i]: f"FG{labels[i]:02d}" for i in range(n_alleles)}

# UPGMA leaf order (used by heatmap and cluster figure)
from scipy.cluster.hierarchy import leaves_list
leaf_order_idx   = leaves_list(Z)
leaf_order_names = [allele_names[i] for i in leaf_order_idx]

# ─────────────────────────────────────────────────────────────────────────────
# PART 2b — Allele similarity heatmap
# ─────────────────────────────────────────────────────────────────────────────

OUT_HEATMAP = "SRK_allele_similarity_heatmap.pdf"
print(f"\nWriting allele similarity heatmap to {OUT_HEATMAP} ...")

# Reorder distance matrix by UPGMA leaf order
dist_ordered = dist_hv[np.ix_(leaf_order_idx, leaf_order_idx)]
sim_ordered  = 1.0 - dist_ordered   # convert to similarity

# Functional group and class colour strips
fg_ids_ord  = [allele_to_fg[n] for n in leaf_order_names]
fg_ids_uniq = sorted(set(allele_to_fg.values()))
cmap_fg     = matplotlib.colormaps.get_cmap("tab10").resampled(max(len(fg_ids_uniq), 2))
fg_col_map  = {fg: cmap_fg(i) for i, fg in enumerate(fg_ids_uniq)}
strip_cols  = np.array([fg_col_map[fg] for fg in fg_ids_ord])   # (n, 4) RGBA

fig = plt.figure(figsize=(11, 10))

# Layout: left colour strip | heatmap | colourbar
# (right strip removed — left strip + colourbar are sufficient)
gs = fig.add_gridspec(
    2, 3,
    width_ratios=[0.025, 1, 0.045],
    height_ratios=[1, 0.025],
    hspace=0.025, wspace=0.025,
)

ax_hm     = fig.add_subplot(gs[0, 1])   # main heatmap
ax_lstrip = fig.add_subplot(gs[0, 0])  # left class strip
ax_bstrip = fig.add_subplot(gs[1, 1])  # bottom class strip
ax_cb     = fig.add_subplot(gs[0, 2])  # colourbar

# Colour scale: stretch across within-class similarity range so W/N/P_within
# variation is visible. Between-class values (Allele_061) fall below vmin and
# saturate to the darkest colour — correctly showing them as maximally different.
within_class_sims = sim_ordered[sim_ordered < 1.0]   # exclude diagonal
within_class_sims = within_class_sims[within_class_sims > 0.5]  # exclude between-class
vmin_hm = max(0.0, within_class_sims.min() - 0.02)   # small buffer below within-class min
vmax_hm = 1.0

im = ax_hm.imshow(sim_ordered, aspect="auto", interpolation="nearest",
                  cmap="RdYlGn", vmin=vmin_hm, vmax=vmax_hm, origin="upper")

sim_thr = 1.0 - WITHIN_CLASS_THRESHOLD
# Find positions in the ordered matrix where similarity transitions
for pos in range(n_alleles - 1):
    if fg_ids_ord[pos] != fg_ids_ord[pos + 1]:
        ax_hm.axhline(pos + 0.5, color="white", linewidth=1.5, zorder=3)
        ax_hm.axvline(pos + 0.5, color="white", linewidth=1.5, zorder=3)

# Tick labels
ax_hm.set_xticks(range(n_alleles))
ax_hm.set_yticks(range(n_alleles))
ax_hm.set_xticklabels(leaf_order_names, rotation=90, fontsize=5.5)
ax_hm.set_yticklabels(leaf_order_names, fontsize=5.5)
ax_hm.tick_params(length=0)

# Left class strip
strip_img_left = strip_cols[:, np.newaxis, :]   # (n, 1, 4)
ax_lstrip.imshow(strip_img_left, aspect="auto", interpolation="nearest", origin="upper")
ax_lstrip.set_xticks([])
ax_lstrip.set_yticks([])
for spine in ax_lstrip.spines.values():
    spine.set_visible(False)

# Bottom class strip
strip_img_bot = strip_cols[np.newaxis, :, :]    # (1, n, 4)
ax_bstrip.imshow(strip_img_bot, aspect="auto", interpolation="nearest", origin="upper")
ax_bstrip.set_xticks([])
ax_bstrip.set_yticks([])
for spine in ax_bstrip.spines.values():
    spine.set_visible(False)

# Colourbar
cb = plt.colorbar(im, cax=ax_cb)
cb.set_label("HV similarity  (1 − p-distance)", fontsize=9, labelpad=6)
cb.ax.tick_params(labelsize=8)

# Annotate W / N / P_within thresholds on colourbar
# Convert similarity values to normalised colourbar position [0, 1]
def cb_pos(sim_val):
    return (sim_val - vmin_hm) / (vmax_hm - vmin_hm)

if vmin_hm <= sim_thr <= vmax_hm:
    cb.ax.axhline(cb_pos(sim_thr), color="black", linewidth=1.2, linestyle="--")
    cb.ax.text(0.5, cb_pos(sim_thr) + 0.02, f"N/P_within  d={WITHIN_CLASS_THRESHOLD:.2f}",
               ha="center", va="bottom", fontsize=6.5, color="black",
               transform=cb.ax.transAxes)
cb.ax.axhline(cb_pos(1.0), color="#d62728", linewidth=1.8)
cb.ax.text(0.5, cb_pos(1.0) - 0.02, "W  d=0",
           ha="center", va="top", fontsize=6.5, color="#d62728",
           transform=cb.ax.transAxes)
cb.ax.text(0.5, 0.02, "Between-class\n(saturated)",
           ha="center", va="bottom", fontsize=6.5, color="grey",
           transform=cb.ax.transAxes)

# Class legend — inside the heatmap (lower-right corner)
legend_patches = [mpatches.Patch(color=fg_col_map[fg], label=fg) for fg in fg_ids_uniq]
ax_hm.legend(handles=legend_patches, loc="lower right", fontsize=8,
             title="Class", title_fontsize=8, framealpha=0.9)

ax_hm.set_title(
    f"Allele pairwise HV similarity  ({len(hv_positions)} HV positions)\n"
    f"Ordered by UPGMA · white lines = class boundaries · "
    f"colour scale spans within-class range ({vmin_hm:.3f}–1.0)",
    fontsize=10, pad=8)

plt.savefig(OUT_HEATMAP, format="pdf", dpi=150, bbox_inches="tight")
plt.savefig(os.path.join("figures", OUT_HEATMAP.replace(".pdf", ".png")),
            format="png", dpi=200, bbox_inches="tight")
plt.close()
print(f"Saved {OUT_HEATMAP}")

# ─────────────────────────────────────────────────────────────────────────────
# PART 3 — Functional groups, synonymy candidates, and cross design
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("PART 3 — Functional groups, synonymy candidates, and cross design")
print("=" * 60)

# Load AAAA individuals
df_zyg          = pd.read_csv(ZYGO_TSV, sep="\t", encoding="utf-8-sig")
aaaa_individuals = set(df_zyg.loc[df_zyg["Genotype"] == "AAAA", "Individual"])
print(f"\nAAA individuals in zygosity data: {len(aaaa_individuals)}")

df_allele_table = pd.read_csv(ALLELE_TABLE, sep="\t", encoding="utf-8-sig")

aaaa_allele = {}
for ind in aaaa_individuals:
    rows = df_allele_table[df_allele_table["Individual"] == ind]
    if rows.empty:
        continue
    bins = rows["Allele"].unique()
    if len(bins) == 1:
        aaaa_allele[ind] = bins[0]
    else:
        counts = rows.groupby("Allele")["Count"].sum()
        aaaa_allele[ind] = counts.idxmax()

print(f"AAAA individuals with allele assignment: {len(aaaa_allele)}")
unassigned = len(aaaa_individuals) - len(aaaa_allele)
if unassigned:
    print(f"  ({unassigned} AAAA individuals not found in allele table — excluded)")

aaaa_count = {}
for allele in aaaa_allele.values():
    aaaa_count[allele] = aaaa_count.get(allele, 0) + 1

# Number of distinct individuals carrying each allele (all genotypes)
obs_count = (df_allele_table.groupby("Allele")["Individual"]
             .nunique().to_dict())

# ── Functional group table ────────────────────────────────────────────────────

fg_rows = []
for allele in allele_names:
    fg     = allele_to_fg[allele]
    n_aaaa = aaaa_count.get(allele, 0)
    fg_rows.append({
        "Allele":          allele,
        "FunctionalGroup": fg,
        "N_AAAA":          n_aaaa,
        "CrossPower":      ("full"      if n_aaaa >= 2 else
                            "singleton" if n_aaaa == 1 else "none"),
    })

df_fg = pd.DataFrame(fg_rows).sort_values(["FunctionalGroup", "Allele"])
df_fg.to_csv(OUT_FUNC_GROUPS, sep="\t", index=False)
print(f"\nFunctional group table written to {OUT_FUNC_GROUPS}")

print("\n── Functional group summary ──")
for fg, grp in df_fg.groupby("FunctionalGroup"):
    n_aaaa_fg    = grp["N_AAAA"].sum()
    alleles_list = ", ".join(grp["Allele"].tolist())
    tag = " ← synonymy candidates" if len(grp) > 1 else ""
    print(f"  {fg}: {len(grp):2d} allele(s), {n_aaaa_fg:3d} AAAA individuals{tag}")
    print(f"       {alleles_list}")

print("\n── Cross power ──")
for power, grp in df_fg.groupby("CrossPower"):
    notes = {"full": "≥2 AAAA → full 3-tier design",
             "singleton": "1 AAAA → cross partner only",
             "none": "no AAAA → must use AAAB/AABB"}
    print(f"  {power:9s}: {len(grp):2d} allele bins  ({notes[power]})")

# ── Synonymy candidates ───────────────────────────────────────────────────────

syn_rows = []
fg_allele_map = df_fg.groupby("FunctionalGroup")["Allele"].apply(list).to_dict()
for fg, alleles in fg_allele_map.items():
    if len(alleles) < 2:
        continue
    for ii, a1 in enumerate(alleles):
        for a2 in alleles[ii + 1:]:
            idx1 = allele_names.index(a1)
            idx2 = allele_names.index(a2)
            n1   = aaaa_count.get(a1, 0)
            n2   = aaaa_count.get(a2, 0)
            hv_d = round(dist_hv[idx1, idx2], 4)
            if hv_d == 0.0:
                tier = CAT_W
                hyp  = "HV-identical — strong synonymy prediction (W cross)"
            elif hv_d < WITHIN_CLASS_THRESHOLD:
                tier = CAT_N
                hyp  = "Small HV difference — synonymy test (N cross)"
            else:
                tier = CAT_P_WITHIN
                hyp  = "Substantial within-class HV divergence — predict compatible (P_within)"
            syn_rows.append({
                "Allele_A":        a1,
                "Allele_B":        a2,
                "FunctionalGroup": fg,
                "HV_distance":     hv_d,
                "Cross_tier":      tier,
                "AAAA_A":          n1,
                "AAAA_B":          n2,
                "Testable":        "yes" if (n1 >= 1 and n2 >= 1) else "no",
                "Hypothesis":      hyp,
            })

df_syn = (pd.DataFrame(syn_rows)
          if syn_rows
          else pd.DataFrame(columns=["Allele_A", "Allele_B", "FunctionalGroup",
                                     "HV_distance", "AAAA_A", "AAAA_B",
                                     "Testable", "Hypothesis"]))
df_syn = df_syn.sort_values(["FunctionalGroup", "HV_distance"])
df_syn.to_csv(OUT_SYNONYMY, sep="\t", index=False)

n_testable = int((df_syn["Testable"] == "yes").sum()) if len(df_syn) else 0
print(f"\nSynonymy candidates: {len(df_syn)} pairs written to {OUT_SYNONYMY}")
print(f"  Testable with available AAAA individuals: {n_testable}")

# ── W / N / P cross design matrix ────────────────────────────────────────────

aaaa_list  = sorted(aaaa_allele.keys())
cross_rows = []

for i, ind_a in enumerate(aaaa_list):
    for ind_b in aaaa_list[i + 1:]:
        allele_a = aaaa_allele[ind_a]
        allele_b = aaaa_allele[ind_b]
        fg_a     = allele_to_fg.get(allele_a, "?")
        fg_b     = allele_to_fg.get(allele_b, "?")

        if fg_a != fg_b:
            category = CAT_P_CROSS
            expected = "Compatible — between-class positive control (guaranteed distinct)"
        elif d == 0.0:
            category = CAT_W
            expected = "Incompatible — HV-identical alleles (strongest synonymy prediction)"
        elif d < WITHIN_CLASS_THRESHOLD:
            category = CAT_N
            expected = "Synonymy test — small HV differences; does this substitution matter?"
        else:
            category = CAT_P_WITHIN
            expected = "Compatible — within-class but substantial HV divergence"

        idx_a = allele_names.index(allele_a) if allele_a in allele_names else -1
        idx_b = allele_names.index(allele_b) if allele_b in allele_names else -1
        d = dist_hv[idx_a, idx_b] if idx_a >= 0 and idx_b >= 0 else np.nan

        cross_rows.append({
            "Mother":        ind_a,
            "Father":        ind_b,
            "Mother_allele": allele_a,
            "Father_allele": allele_b,
            "Mother_FG":     fg_a,
            "Father_FG":     fg_b,
            "HV_dist":       round(d, 4),
            "Category":      category,
            "Expected":      expected,
        })

df_cross = pd.DataFrame(cross_rows)
cat_order = {CAT_W: 0, CAT_N: 1, CAT_P_WITHIN: 2, CAT_P_CROSS: 3}
df_cross["_sort"] = df_cross["Category"].map(cat_order)
df_cross = df_cross.sort_values(["_sort", "HV_dist"]).drop(columns="_sort")
df_cross.to_csv(OUT_CROSS_DESIGN, sep="\t", index=False)

n_w  = (df_cross["Category"] == CAT_W).sum()
n_n  = (df_cross["Category"] == CAT_N).sum()
n_pw = (df_cross["Category"] == CAT_P_WITHIN).sum()
n_pc = (df_cross["Category"] == CAT_P_CROSS).sum()

print(f"\nCross design matrix written to {OUT_CROSS_DESIGN}")
print(f"  W        (HV-identical, expected incompatible)     : {n_w:5d} pairs")
print(f"  N        (small HV diff, synonymy test)            : {n_n:5d} pairs")
print(f"  P_within (within-class, substantial HV divergence) : {n_pw:5d} pairs")
print(f"  P_cross  (between-class, guaranteed compatible)    : {n_pc:5d} pairs")
print(f"  Total                                              : {len(df_cross):5d} pairs")

# ─────────────────────────────────────────────────────────────────────────────
# PART 3b — Cross design summary figure
# ─────────────────────────────────────────────────────────────────────────────

OUT_SUMMARY_FIG = "SRK_cross_design_summary.pdf"
print(f"\nWriting cross design summary figure to {OUT_SUMMARY_FIG} ...")

from matplotlib.patches import FancyBboxPatch

# Within-class pairwise HV distances (majority class only)
fg_majority      = df_fg.groupby("FunctionalGroup").size().idxmax()
majority_alleles = df_fg[df_fg["FunctionalGroup"] == fg_majority]["Allele"].tolist()

within_dists = []
for ii, a1 in enumerate(majority_alleles):
    for a2 in majority_alleles[ii + 1:]:
        within_dists.append(dist_hv[allele_names.index(a1), allele_names.index(a2)])
within_dists = np.array(within_dists)

n_w_d  = (within_dists == 0).sum()
n_n_d  = ((within_dists > 0) & (within_dists < WITHIN_CLASS_THRESHOLD)).sum()
n_pw_d = (within_dists >= WITHIN_CLASS_THRESHOLD).sum()
d_max  = within_dists.max()

fig, axes = plt.subplots(1, 3, figsize=(18, 6.5))
fig.subplots_adjust(wspace=0.38)

# ── Panel A: HV distance histogram with category zones ───────────────────────
ax_a = axes[0]

bins_fine = np.linspace(0, d_max * 1.06, 55)
ax_a.hist(within_dists[within_dists >= WITHIN_CLASS_THRESHOLD], bins=bins_fine,
          color="#1f77b4", alpha=0.85,
          label=f"P_within   d ≥ {WITHIN_CLASS_THRESHOLD:.2f}   ({n_pw_d:,} pairs)")
ax_a.hist(within_dists[(within_dists > 0) & (within_dists < WITHIN_CLASS_THRESHOLD)],
          bins=bins_fine, color="#ff7f0e", alpha=0.85,
          label=f"N   0 < d < {WITHIN_CLASS_THRESHOLD:.2f}   ({n_n_d:,} pairs)")
ax_a.bar(0, n_w_d, width=d_max * 0.018, color="#d62728", alpha=0.85,
         label=f"W   d = 0   ({n_w_d:,} pairs)")
ax_a.axvline(WITHIN_CLASS_THRESHOLD, color="black", linestyle="--", linewidth=1.3)
ax_a.text(WITHIN_CLASS_THRESHOLD + d_max * 0.01,
          ax_a.get_ylim()[1] if ax_a.get_ylim()[1] > 0 else 1,
          f"threshold\n{WITHIN_CLASS_THRESHOLD:.2f}", va="top", fontsize=8, color="black")
ax_a.set_xlabel("HV p-distance (within-class pairs)", fontsize=11)
ax_a.set_ylabel("Number of allele pairs", fontsize=11)
ax_a.set_title("A   HV distance distribution", fontsize=11, fontweight="bold")
ax_a.legend(fontsize=8.5, loc="upper right", framealpha=0.9)
ax_a.set_xlim(-d_max * 0.02, d_max * 1.08)

# ── Panel B: Cross category schematic ────────────────────────────────────────
ax_b = axes[1]
ax_b.set_xlim(0, 10)
ax_b.set_ylim(0, 10)
ax_b.axis("off")
ax_b.set_title("B   Cross category logic", fontsize=11, fontweight="bold")

def fbox(ax, cx, cy, w, h, fc, ec, alpha=0.18, radius=0.25):
    p = FancyBboxPatch((cx - w/2, cy - h/2), w, h,
                       boxstyle=f"round,pad={radius}",
                       facecolor=fc, edgecolor=ec, linewidth=1.5,
                       alpha=alpha, zorder=2)
    ax.add_patch(p)

def fbox_solid(ax, cx, cy, w, h, fc, ec, alpha=0.9, radius=0.2):
    p = FancyBboxPatch((cx - w/2, cy - h/2), w, h,
                       boxstyle=f"round,pad={radius}",
                       facecolor=fc, edgecolor=ec, linewidth=1.5,
                       alpha=alpha, zorder=4)
    ax.add_patch(p)

def arc_arrow(ax, x1, y1, x2, y2, color, rad=0.35, lw=2.2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                connectionstyle=f"arc3,rad={rad}",
                                mutation_scale=15), zorder=5)

# Class backgrounds
fbox(ax_b, 3.8, 5.5, 6.8, 8.2, "#1f77b4", "#1f77b4", alpha=0.07)
ax_b.text(3.8, 9.55, "Class I  (majority)", ha="center", fontsize=9,
          color="#1f77b4", fontweight="bold")
fbox(ax_b, 8.5, 5.5, 2.2, 3.5, "#9467bd", "#9467bd", alpha=0.1)
ax_b.text(8.5, 7.45, "Class II", ha="center", fontsize=9,
          color="#9467bd", fontweight="bold")

# Allele bins — Class I
bin_pos = {"A": (2.0, 8.2), "B": (5.5, 8.2), "C": (2.0, 5.5), "D": (5.5, 5.5)}
for lbl, (bx, by) in bin_pos.items():
    fbox_solid(ax_b, bx, by, 1.6, 0.75, "#1f77b4", "white", alpha=0.85)
    ax_b.text(bx, by, f"Bin {lbl}", ha="center", va="center",
              fontsize=9, color="white", fontweight="bold", zorder=5)

# Class II bin
fbox_solid(ax_b, 8.5, 5.5, 1.6, 0.75, "#9467bd", "white", alpha=0.85)
ax_b.text(8.5, 5.5, "Bin E", ha="center", va="center",
          fontsize=9, color="white", fontweight="bold", zorder=5)

# W: Bin A → Bin A (self-loop, same bin, d=0)
ax_b.annotate("", xy=(2.55, 8.55), xytext=(1.45, 8.55),
              arrowprops=dict(arrowstyle="-|>", color="#d62728", lw=2.2,
                              connectionstyle="arc3,rad=-1.2", mutation_scale=14), zorder=5)
ax_b.text(2.0, 9.35, "W  d = 0\nHV-identical → Incompatible",
          ha="center", fontsize=7.5, color="#d62728", fontweight="bold",
          bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#d62728", alpha=0.92))

# N: Bin A → Bin B (different bin, small HV diff)
arc_arrow(ax_b, 2.8, 8.2, 4.7, 8.2, "#ff7f0e", rad=-0.4)
ax_b.text(3.75, 9.2, "N  0 < d < 0.04\nSmall HV diff → Synonymy test",
          ha="center", fontsize=7.5, color="#ff7f0e", fontweight="bold",
          bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#ff7f0e", alpha=0.92))

# P_within: Bin A → Bin C (same class, large HV diff)
arc_arrow(ax_b, 2.0, 7.83, 2.0, 5.88, "#1f77b4", rad=0.5)
ax_b.text(0.55, 6.8, "P_within\nd ≥ 0.04\n→ Compatible",
          ha="center", fontsize=7.5, color="#1f77b4", fontweight="bold",
          bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#1f77b4", alpha=0.92))

# P_cross: Bin B → Bin E (different class)
arc_arrow(ax_b, 6.3, 8.2, 7.7, 5.88, "#2ca02c", rad=-0.3)
ax_b.text(7.6, 7.5, "P_cross\nDiff. class\n→ Compatible\n(guaranteed)",
          ha="center", fontsize=7.5, color="#2ca02c", fontweight="bold",
          bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#2ca02c", alpha=0.92))

# HV distance label
ax_b.text(3.75, 1.0,
          "HV distance reflects divergence at 75 hypervariable S-domain positions",
          ha="center", fontsize=8, color="grey", style="italic")

# ── Panel C: N-cross outcome interpretation ───────────────────────────────────
ax_c = axes[2]
ax_c.set_xlim(0, 10)
ax_c.set_ylim(0, 10)
ax_c.axis("off")
ax_c.set_title("C   Interpreting N-cross outcomes", fontsize=11, fontweight="bold")

def cbox(ax, cx, cy, w, h, fc, ec, txt, fs=8.5):
    fbox_solid(ax, cx, cy, w, h, fc, ec, alpha=0.15, radius=0.2)
    ax.text(cx, cy, txt, ha="center", va="center", fontsize=fs,
            color=ec, fontweight="bold", zorder=5, multialignment="center")

def varrow(ax, x, y1, y2, color="grey", label=""):
    ax.annotate("", xy=(x, y2), xytext=(x, y1),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.6,
                                mutation_scale=12), zorder=5)
    if label:
        ax.text(x + 0.2, (y1 + y2) / 2, label, va="center",
                fontsize=8, color=color)

def harrow(ax, x1, x2, y, color="grey"):
    ax.annotate("", xy=(x2, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.6,
                                mutation_scale=12), zorder=5)

# Top: N cross
cbox(ax_c, 5, 9.1, 7.5, 1.0, "#ff7f0e", "#ff7f0e",
     "N cross:  Bin A × Bin B  (0 < d < 0.04)", fs=9)

varrow(ax_c, 5, 8.55, 7.85, "grey")
ax_c.text(5, 7.65, "Seed yield?", ha="center", fontsize=9.5,
          color="grey", style="italic")

# Branch arrows
ax_c.annotate("", xy=(2.5, 6.55), xytext=(5, 7.4),
              arrowprops=dict(arrowstyle="-|>", color="#d62728", lw=1.6,
                              mutation_scale=12), zorder=5)
ax_c.text(2.5, 7.1, "Low / zero\n(incompatible)", ha="center",
          fontsize=8, color="#d62728", fontweight="bold")

ax_c.annotate("", xy=(7.5, 6.55), xytext=(5, 7.4),
              arrowprops=dict(arrowstyle="-|>", color="#2ca02c", lw=1.6,
                              mutation_scale=12), zorder=5)
ax_c.text(7.5, 7.1, "High\n(compatible)", ha="center",
          fontsize=8, color="#2ca02c", fontweight="bold")

# Left branch
cbox(ax_c, 2.5, 5.9, 4.2, 1.0, "#d62728", "#d62728",
     "Bins A & B share\nSI specificity", fs=8.5)
varrow(ax_c, 2.5, 5.35, 4.45, "#d62728")
cbox(ax_c, 2.5, 4.0, 4.2, 0.85, "#d62728", "#d62728",
     "Merge bins A + B\n→ Synonymous alleles", fs=8.5)

# Right branch
cbox(ax_c, 7.5, 5.9, 4.2, 1.0, "#2ca02c", "#2ca02c",
     "Bins A & B have\ndistinct specificities", fs=8.5)
varrow(ax_c, 7.5, 5.35, 4.45, "#2ca02c")
cbox(ax_c, 7.5, 4.0, 4.2, 0.85, "#2ca02c", "#2ca02c",
     "Retain separate bins\n→ Allele bins confirmed", fs=8.5)

# Converge at bottom
ax_c.annotate("", xy=(5, 2.65), xytext=(2.5, 3.55),
              arrowprops=dict(arrowstyle="-|>", color="grey", lw=1.4,
                              mutation_scale=11), zorder=5)
ax_c.annotate("", xy=(5, 2.65), xytext=(7.5, 3.55),
              arrowprops=dict(arrowstyle="-|>", color="grey", lw=1.4,
                              mutation_scale=11), zorder=5)
cbox(ax_c, 5, 2.1, 7.0, 1.0, "grey", "grey",
     "Refined allele bin definitions\n→ Updated crossing design", fs=8.5)

ax_c.text(5, 0.55,
          "W crosses validate bins · P crosses confirm HV threshold",
          ha="center", fontsize=8, color="grey", style="italic")

fig.suptitle("SRK allele hypothesis testing — cross design framework",
             fontsize=13, fontweight="bold", y=1.01)

plt.tight_layout()
plt.savefig(OUT_SUMMARY_FIG, format="pdf", dpi=150, bbox_inches="tight")
plt.savefig(os.path.join("figures", OUT_SUMMARY_FIG.replace(".pdf", ".png")),
            format="png", dpi=200, bbox_inches="tight")
plt.close()
print(f"Saved {OUT_SUMMARY_FIG}")

# ─────────────────────────────────────────────────────────────────────────────
# PART 3c — Synonymy network (two-panel: W-only  vs  W+N)
#
# Panel A — W-only edges (d=0, HV-identical): tight synonymy groups; layout
#           separates them so each cluster is visually distinct.
# Panel B — same node positions; N-edges (dashed) added as overlay to show
#           how small HV differences chain W-groups into larger components.
#
# Node colour encodes W-only component membership → colour is constant across
# both panels, making the chaining effect immediately legible.
# ─────────────────────────────────────────────────────────────────────────────

try:
    import networkx as nx
    _HAS_NX = True
except ImportError:
    print("\nWARNING: networkx not installed — skipping synonymy network figure.")
    print("  Install with: pip install networkx")
    _HAS_NX = False

OUT_NET_W        = "SRK_synonymy_network_W.pdf"
OUT_NET_N        = "SRK_synonymy_network_N.pdf"
OUT_SYN_GROUPS   = "SRK_synonymy_groups.csv"

if _HAS_NX:
    print(f"\nBuilding synonymy network ...")

    # ── Build two graphs ──────────────────────────────────────────────────────
    G_w  = nx.Graph()   # W edges only  (d = 0)
    G_wn = nx.Graph()   # W + N edges   (d < WITHIN_CLASS_THRESHOLD)

    for allele in allele_names:
        for G_ in (G_w, G_wn):
            G_.add_node(allele,
                        fg=allele_to_fg[allele],
                        n_aaaa=aaaa_count.get(allele, 0))

    for i, a1 in enumerate(allele_names):
        for j in range(i + 1, n_alleles):
            a2   = allele_names[j]
            d_ij = dist_hv[i, j]
            if allele_to_fg[a1] != allele_to_fg[a2]:
                continue                              # between-class: skip
            if d_ij >= WITHIN_CLASS_THRESHOLD:
                continue                              # P_within: skip
            cat = CAT_W if d_ij == 0.0 else CAT_N
            G_wn.add_edge(a1, a2, distance=round(d_ij, 4), category=cat)
            if d_ij == 0.0:
                G_w.add_edge(a1, a2, distance=0.0, category=CAT_W)

    # ── Connected components ──────────────────────────────────────────────────
    w_comps  = sorted(nx.connected_components(G_w),  key=len, reverse=True)
    wn_comps = sorted(nx.connected_components(G_wn), key=len, reverse=True)
    w_multi  = [c for c in w_comps  if len(c) > 1]
    w_iso    = [list(c)[0] for c in w_comps if len(c) == 1]
    wn_multi = [c for c in wn_comps if len(c) > 1]
    wn_iso   = [list(c)[0] for c in wn_comps if len(c) == 1]

    n_w_edges   = G_w.number_of_edges()
    n_wn_edges  = G_wn.number_of_edges()
    n_n_edges   = n_wn_edges - n_w_edges

    print(f"  W-only graph   : {n_w_edges} edges  →  {len(w_multi)} groups, {len(w_iso)} isolated")
    print(f"  W+N graph      : {n_wn_edges} edges  →  {len(wn_multi)} groups, {len(wn_iso)} isolated")
    print(f"    W (d = 0)    : {n_w_edges}")
    print(f"    N (0 < d)    : {n_n_edges}")

    # ── Node colours: one colour per W-only multi-node component; grey = isolated ─
    n_w_multi  = len(w_multi)
    cmap_comp  = matplotlib.colormaps.get_cmap("tab20").resampled(max(n_w_multi, 2))
    w_node_col = {}
    for ci, comp in enumerate(w_multi):
        col = cmap_comp(ci % 20)
        for allele in comp:
            w_node_col[allele] = col
    for allele in w_iso:
        w_node_col[allele] = (0.75, 0.75, 0.75, 1.0)

    # ── Layout: per-W-component tiling ───────────────────────────────────────
    # Each W-only multi-node component is laid out independently with
    # kamada_kawai_layout, then tiled in a grid.  Isolated nodes go in a
    # dense row below.  This keeps the W-groups visually separated so
    # Panel A is clean, and N-edges in Panel B visibly bridge them.
    CELL  = 2.5     # grid cell size (data units) per W-component
    NCOLS = max(1, int(np.ceil(np.sqrt(max(n_w_multi, 1)))))

    pos = {}
    w_comp_layouts = []   # list of (comp_set, sub_pos_dict)

    for comp in w_multi:
        sub = G_w.subgraph(comp)
        if len(comp) == 2:
            nodes_c = list(comp)
            sub_pos = {nodes_c[0]: np.array([0.0, 0.0]),
                       nodes_c[1]: np.array([1.0, 0.0])}
        else:
            sub_pos = nx.kamada_kawai_layout(sub, weight=None)
        w_comp_layouts.append((comp, sub_pos))

    for ci, (comp, sub_pos) in enumerate(w_comp_layouts):
        col_ci = ci % NCOLS
        row_ci = ci // NCOLS
        offset = np.array([col_ci * CELL, -row_ci * CELL])
        for node, xy in sub_pos.items():
            pos[node] = np.array(xy) + offset

    n_rows_used = (n_w_multi - 1) // NCOLS + 1 if w_multi else 0
    iso_y_start = -(n_rows_used * CELL + 1.5)
    ISO_COLS    = max(1, min(12, len(w_iso)))
    for ii, node in enumerate(w_iso):
        ic = ii % ISO_COLS
        ir = ii // ISO_COLS
        pos[node] = np.array([ic * 1.2, iso_y_start - ir * 1.2])

    # ── Canvas size from pos extents ─────────────────────────────────────────
    all_xy = np.array(list(pos.values()))
    x_min, x_max = all_xy[:, 0].min(), all_xy[:, 0].max()
    y_min, y_max = all_xy[:, 1].min(), all_xy[:, 1].max()
    pad = 1.0
    ax_xlim = (x_min - pad, x_max + pad)
    ax_ylim = (y_min - pad, y_max + pad)

    node_list   = list(G_w.nodes())   # consistent node order for both panels
    node_colors = [w_node_col[n] for n in node_list]
    node_sizes  = [max(100, obs_count.get(n, 0) * 55 + 120) for n in node_list]

    # ── Synonymy-collapsed allele count (W-only, HV-identical) ──────────────
    n_merged_W        = len(w_multi) + len(w_iso)
    n_collapsed       = sum(len(c) for c in w_multi) - len(w_multi)
    print(f"\n── Synonymy-collapsed allele count (W-only, HV-identical) ──")
    print(f"  Original allele count : {n_alleles}")
    print(f"  After merging W-groups: {n_merged_W}  "
          f"({len(w_multi)} groups + {len(w_iso)} singletons)")
    print(f"  Alleles collapsed     : {n_collapsed} bins removed "
          f"({sum(len(c) for c in w_multi)} alleles → {len(w_multi)} groups)")

    # ── Figure A: two panels — W-groups (left) | isolated alleles (right) ────
    w_iso_set     = set(w_iso)
    grp_node_list = [n for n in node_list if n not in w_iso_set]
    grp_colors    = [w_node_col[n] for n in grp_node_list]
    grp_sizes     = [max(100, obs_count.get(n, 0) * 55 + 120) for n in grp_node_list]

    # Main-panel limits: W-group nodes only (removes the iso row from extent)
    grp_xy   = np.array([pos[n] for n in grp_node_list]) if grp_node_list else np.array([[0, 0]])
    xg_min, xg_max = grp_xy[:, 0].min(), grp_xy[:, 0].max()
    yg_min, yg_max = grp_xy[:, 1].min(), grp_xy[:, 1].max()
    ax_main_xlim = (xg_min - pad, xg_max + pad)
    ax_main_ylim = (yg_min - pad, yg_max + pad)
    fig_h = max(8, (yg_max - yg_min + 2 * pad) * 1.8)

    # Isolated alleles: grid positions sorted by obs_count descending
    w_iso_sorted  = sorted(w_iso, key=lambda a: obs_count.get(a, 0), reverse=True)
    ISO_GRID_COLS = 3
    iso_grid_pos  = {}
    for ii, node in enumerate(w_iso_sorted):
        ic = ii % ISO_GRID_COLS
        ir = ii // ISO_GRID_COLS
        iso_grid_pos[node] = np.array([ic * 1.8, -ir * 1.8], dtype=float)
    iso_sizes = [max(100, obs_count.get(n, 0) * 55 + 120) for n in w_iso_sorted]

    fig_a, (ax_main, ax_iso) = plt.subplots(
        1, 2, figsize=(fig_h * 1.5, fig_h),
        gridspec_kw={"width_ratios": [3, 1]})
    for ax in (ax_main, ax_iso):
        ax.set_aspect("equal")
        ax.axis("off")

    ax_main.set_xlim(*ax_main_xlim)
    ax_main.set_ylim(*ax_main_ylim)

    # W-group nodes and edges in ax_main
    nx.draw_networkx_edges(G_w, pos, edgelist=list(G_w.edges()),
                           edge_color="#d62728", width=3.0, alpha=0.80, ax=ax_main)
    nx.draw_networkx_nodes(G_w, pos, nodelist=grp_node_list,
                           node_color=grp_colors, node_size=grp_sizes,
                           alpha=0.92, linewidths=0.8, edgecolors="white", ax=ax_main)
    nx.draw_networkx_labels(G_w, pos,
                            labels={n: n for n in grp_node_list},
                            font_size=5.5, font_color="white",
                            font_weight="bold", ax=ax_main)

    for comp, _ in w_comp_layouts:
        gxs = [pos[n][0] for n in comp]
        gys = [pos[n][1] for n in comp]
        ax_main.text(np.mean(gxs), max(gys) + 0.28,
                     f"W-group  ({len(comp)} alleles)",
                     ha="center", va="bottom", fontsize=5.5,
                     color="dimgrey", style="italic")

    ax_main.set_title(
        f"HV-identical synonymy groups  (d = 0)\n"
        f"{n_w_edges} W-edges  ·  {len(w_multi)} synonymy groups",
        fontsize=11, fontweight="bold", pad=8)

    # Isolated alleles in ax_iso: grey nodes, dark labels (allele ID + obs count)
    nx.draw_networkx_nodes(G_w, iso_grid_pos, nodelist=w_iso_sorted,
                           node_color=[(0.75, 0.75, 0.75, 1.0)] * len(w_iso_sorted),
                           node_size=iso_sizes,
                           alpha=0.92, linewidths=0.8, edgecolors="white", ax=ax_iso)
    for node in w_iso_sorted:
        x, y = iso_grid_pos[node]
        ax_iso.text(x, y, f"{node.replace('Allele_', 'A')}\n({obs_count.get(node, 0)})",
                    ha="center", va="center", fontsize=5.5,
                    color="#111111", fontweight="bold")

    if iso_grid_pos:
        iso_xy = np.array(list(iso_grid_pos.values()))
        ax_iso.set_xlim(iso_xy[:, 0].min() - 1.2, iso_xy[:, 0].max() + 1.2)
        ax_iso.set_ylim(iso_xy[:, 1].min() - 1.2, iso_xy[:, 1].max() + 1.2)
    ax_iso.set_title(
        f"Isolated  ({len(w_iso)} alleles)\nno HV-identical partner",
        fontsize=11, fontweight="bold", pad=8)

    # Synonymy count annotation
    fig_a.text(
        0.5, 0.01,
        f"Synonymy-collapsed allele count:  {n_alleles} observed alleles  →  "
        f"{n_merged_W} functional groups  (merging {n_collapsed} HV-identical bins)",
        ha="center", va="bottom", fontsize=9,
        style="italic", color="dimgrey",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85, ec="lightgrey"))

    legend_elems_a = [
        plt.Line2D([0], [0], color="#d62728", linewidth=3,
                   label="W-edge  d = 0  (HV-identical)"),
        mpatches.Patch(facecolor=(0.75, 0.75, 0.75), edgecolor="white",
                       label="Isolated allele (no W-edges)"),
        mpatches.Patch(facecolor="none", edgecolor="none",
                       label="Node colour = W-only group  ·  node size ∝ individuals observed  ·  label = allele ID (obs. count)"),
    ]
    fig_a.legend(handles=legend_elems_a, loc="lower center", ncol=3,
                 fontsize=9, framealpha=0.9, bbox_to_anchor=(0.5, -0.04))

    fig_a.suptitle(
        "SRK allele synonymy network — HV-identical W-groups",
        fontsize=12, fontweight="bold", y=1.01)

    plt.tight_layout()
    plt.savefig(OUT_NET_W, format="pdf", dpi=150, bbox_inches="tight")
    plt.savefig(os.path.join("figures", OUT_NET_W.replace(".pdf", ".png")),
                format="png", dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved {OUT_NET_W}")

    # ── Figure B: Condensed super-node graph — N-connectivity between W-groups ─
    # Each node = one W-only group (or one isolated allele).
    # An edge exists when ≥1 N-edge (0 < d < threshold) bridges the two groups.
    # Edge width ∝ log(number of N-bridges).

    allele_to_cid = {}
    cid_alleles   = {}
    cid_color     = {}
    cid_n_alleles = {}
    cid_obs       = {}

    for ci, comp in enumerate(w_multi):
        cid = f"Group {ci + 1}"
        for a in comp:
            allele_to_cid[a] = cid
        cid_alleles[cid]   = sorted(comp)
        cid_color[cid]     = cmap_comp(ci % 20)
        cid_n_alleles[cid] = len(comp)
        cid_obs[cid]       = sum(obs_count.get(a, 0) for a in comp)

    for a in w_iso:
        allele_to_cid[a] = a
        cid_alleles[a]   = [a]
        cid_color[a]     = (0.75, 0.75, 0.75, 1.0)
        cid_n_alleles[a] = 1
        cid_obs[a]       = obs_count.get(a, 0)

    G_super = nx.Graph()
    for cid in cid_alleles:
        G_super.add_node(cid)

    n_bridges_dict = {}
    for u, v, d in G_wn.edges(data=True):
        if d["category"] != CAT_N:
            continue
        cu, cv = allele_to_cid[u], allele_to_cid[v]
        if cu == cv:
            continue
        key = tuple(sorted([cu, cv]))
        n_bridges_dict[key] = n_bridges_dict.get(key, 0) + 1

    for (cu, cv), cnt in n_bridges_dict.items():
        G_super.add_edge(cu, cv, n_bridges=cnt)

    super_pos = nx.kamada_kawai_layout(G_super)

    super_nodes   = list(G_super.nodes())
    super_colors  = [cid_color[n] for n in super_nodes]
    super_sizes   = [max(250, cid_obs[n] * 55 + 120) for n in super_nodes]

    s_edge_list   = list(G_super.edges(data=True))
    s_edge_widths = [max(0.8, np.log1p(d["n_bridges"]) * 1.2)
                     for _, _, d in s_edge_list]

    super_comps   = list(nx.connected_components(G_super))
    n_super_multi = len([c for c in super_comps if len(c) > 1])
    n_super_iso   = len([c for c in super_comps if len(c) == 1])

    fig_b, ax_b = plt.subplots(1, 1, figsize=(12, 12))
    ax_b.set_aspect("equal")
    ax_b.axis("off")

    nx.draw_networkx_edges(G_super, super_pos,
                           edgelist=[(u, v) for u, v, _ in s_edge_list],
                           width=s_edge_widths,
                           edge_color="#ff7f0e", alpha=0.75, ax=ax_b)
    nx.draw_networkx_nodes(G_super, super_pos, nodelist=super_nodes,
                           node_color=super_colors, node_size=super_sizes,
                           alpha=0.92, linewidths=0.8, edgecolors="white", ax=ax_b)

    grp_nodes  = [n for n in super_nodes if cid_n_alleles[n] > 1]
    iso_snodes = [n for n in super_nodes if cid_n_alleles[n] == 1]

    def _slabel(cid):
        if cid_n_alleles[cid] > 1:
            return f"{cid}\n{cid_n_alleles[cid]} alleles\n{cid_obs[cid]} indiv."
        return cid.replace("Allele_", "A")

    if grp_nodes:
        nx.draw_networkx_labels(G_super, super_pos,
                                labels={n: _slabel(n) for n in grp_nodes},
                                font_size=6.5, font_color="white",
                                font_weight="bold", ax=ax_b)
    if iso_snodes:
        nx.draw_networkx_labels(G_super, super_pos,
                                labels={n: _slabel(n) for n in iso_snodes},
                                font_size=5.0, font_color="#333333",
                                font_weight="normal", ax=ax_b)

    ax_b.set_title(
        f"N-connectivity between W-groups (condensed)\n"
        f"Each node = one W-group or isolated allele  ·  "
        f"{G_super.number_of_edges()} inter-group N-bridges  ·  "
        f"{n_super_multi} merged component(s)  ·  {n_super_iso} isolated",
        fontsize=11, fontweight="bold", pad=8)

    legend_elems_b = [
        plt.Line2D([0], [0], color="#ff7f0e", linewidth=2, alpha=0.75,
                   label=f"N-bridge  0 < d < {WITHIN_CLASS_THRESHOLD}  (synonymy test)"),
        mpatches.Patch(facecolor=(0.75, 0.75, 0.75), edgecolor="white",
                       label="Isolated allele (no W-edges)"),
        mpatches.Patch(facecolor="none", edgecolor="none",
                       label="Node colour = W-only group  ·  node size ∝ individuals observed  ·  edge width ∝ log(N-bridges)"),
    ]
    fig_b.legend(handles=legend_elems_b, loc="lower center", ncol=3,
                 fontsize=9, framealpha=0.9, bbox_to_anchor=(0.5, -0.04))

    fig_b.suptitle(
        "SRK allele synonymy network — N-connectivity between W-groups",
        fontsize=12, fontweight="bold", y=1.01)

    plt.tight_layout()
    plt.savefig(OUT_NET_N, format="pdf", dpi=150, bbox_inches="tight")
    plt.savefig(os.path.join("figures", OUT_NET_N.replace(".pdf", ".png")),
                format="png", dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved {OUT_NET_N}")

    # ── Synonymy groups CSV ───────────────────────────────────────────────────
    syn_grp_rows = []
    for ci, comp in enumerate(w_multi, 1):
        grp_label   = f"W-group {ci}"
        grp_size    = len(comp)
        grp_obs_tot = sum(obs_count.get(a, 0) for a in comp)
        for a in sorted(comp):
            syn_grp_rows.append({
                "Allele":           a,
                "Synonymy_group":   grp_label,
                "Group_size":       grp_size,
                "Obs_count":        obs_count.get(a, 0),
                "AAAA_count":       aaaa_count.get(a, 0),
                "Group_obs_total":  grp_obs_tot,
            })
    for a in sorted(w_iso):
        syn_grp_rows.append({
            "Allele":           a,
            "Synonymy_group":   "Isolated",
            "Group_size":       1,
            "Obs_count":        obs_count.get(a, 0),
            "AAAA_count":       aaaa_count.get(a, 0),
            "Group_obs_total":  obs_count.get(a, 0),
        })
    df_syn_grp = (pd.DataFrame(syn_grp_rows)
                    .sort_values(["Synonymy_group", "Allele"])
                    .reset_index(drop=True))
    df_syn_grp.to_csv(OUT_SYN_GROUPS, index=False)
    print(f"Saved {OUT_SYN_GROUPS}  ({len(df_syn_grp)} alleles: "
          f"{len(w_multi)} W-groups + {len(w_iso)} isolated)")

    # ── Console summary ───────────────────────────────────────────────────────
    print("\n── W-only synonymy groups (HV-identical; strong merge-bin candidates) ──")
    for ci, comp in enumerate(w_multi, 1):
        alleles_sorted = sorted(comp)
        n_aaaa_grp = sum(aaaa_count.get(a, 0) for a in alleles_sorted)
        print(f"  Group {ci}: {len(comp)} alleles,  {n_aaaa_grp} AAAA total")
        print(f"    {', '.join(alleles_sorted)}")

    print("\n── W+N components (chaining via small HV differences) ──")
    for ci, comp in enumerate(wn_multi, 1):
        alleles_sorted = sorted(comp)
        sub = G_wn.subgraph(comp)
        w_cnt = sum(1 for *_, d in sub.edges(data=True) if d["category"] == CAT_W)
        n_cnt = sum(1 for *_, d in sub.edges(data=True) if d["category"] == CAT_N)
        n_aaaa_grp = sum(aaaa_count.get(a, 0) for a in alleles_sorted)
        print(f"  Component {ci}: {len(comp)} alleles  "
              f"({w_cnt} W-edges, {n_cnt} N-edges, {n_aaaa_grp} AAAA total)")
        print(f"    {', '.join(alleles_sorted)}")

# ─────────────────────────────────────────────────────────────────────────────
# PART 4 — Cluster figure: HV dendrogram + AAAA availability
# ─────────────────────────────────────────────────────────────────────────────

print(f"\nWriting HV cluster figure to {OUT_CLUSTER_FIG} ...")

fg_ids     = fg_ids_uniq          # already computed in Part 2b
fg_colours = fg_col_map           # reuse colours from heatmap for consistency

fig, axes = plt.subplots(1, 2, figsize=(16, max(8, n_alleles * 0.35)),
                         gridspec_kw={"width_ratios": [3, 1]})
ax_dendro = axes[0]
ax_bar    = axes[1]

dend = dendrogram(
    Z,
    labels=allele_names,
    orientation="left",
    ax=ax_dendro,
    color_threshold=0,
    above_threshold_color="black",
    leaf_font_size=8,
)

leaf_order = dend["ivl"]   # same order as leaf_order_names
for lbl in ax_dendro.get_yticklabels():
    allele = lbl.get_text()
    fg     = allele_to_fg.get(allele, "?")
    lbl.set_color(fg_colours.get(fg, "black"))
    lbl.set_fontsize(7)

ax_dendro.set_xlabel("HV p-distance", fontsize=10)
ax_dendro.set_title(
    f"Allele functional group clustering\n(UPGMA, {len(hv_positions)} HV positions, "
    f"threshold = auto)",
    fontsize=11)

fg_patches = [mpatches.Patch(color=fg_colours[fg], label=fg) for fg in fg_ids]
ax_dendro.legend(handles=fg_patches, loc="lower right", fontsize=7,
                 title="Functional group", ncol=3)

counts_ordered = [aaaa_count.get(a, 0) for a in leaf_order]
colours_bar    = [fg_colours.get(allele_to_fg.get(a, "?"), "grey") for a in leaf_order]

ax_bar.barh(range(len(leaf_order)), counts_ordered,
            color=colours_bar, edgecolor="none", height=0.7)
ax_bar.axvline(2, color="black", linestyle="--", linewidth=0.8,
               label="n = 2 minimum")
ax_bar.set_yticks([])
ax_bar.set_xlabel("AAAA individuals", fontsize=9)
ax_bar.set_title("AAAA\navailability", fontsize=10)
ax_bar.legend(fontsize=7, loc="lower right")
ax_bar.set_xlim(0, max(counts_ordered) * 1.15 if counts_ordered else 5)

plt.tight_layout()
plt.savefig(OUT_CLUSTER_FIG, format="pdf", dpi=150, bbox_inches="tight")
plt.savefig(os.path.join("figures", OUT_CLUSTER_FIG.replace(".pdf", ".png")),
            format="png", dpi=200, bbox_inches="tight")
plt.close()
print(f"Saved {OUT_CLUSTER_FIG}")

# ─────────────────────────────────────────────────────────────────────────────
# PART 5 — Cross result analysis (when cross data available)
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("PART 5 — Cross result analysis")
print("=" * 60)

if not (CROSS_TSV and os.path.exists(CROSS_TSV)):
    print("\nCROSS_TSV is not set or file not found. Skipping Part 5.")
    print("Set CROSS_TSV = '<filename>' in settings when cross data are available.")
    print("\n── DONE ──")
    print(f"  {OUT_VARIABILITY}    — S-domain variability landscape")
    print(f"  {OUT_HV_DISTS}       — HV-only pairwise distances")
    print(f"  {OUT_FUNC_GROUPS}    — allele functional group assignments")
    print(f"  {OUT_SYNONYMY}       — synonymy candidate pairs")
    print(f"  {OUT_CROSS_DESIGN}   — AAAA × AAAA cross design (W/N/P)")
    print(f"  {OUT_CLUSTER_FIG}    — HV-distance dendrogram")
    if _HAS_NX:
        print(f"  {OUT_NET_W}   — synonymy network (W-only groups)")
        print(f"  {OUT_NET_N}   — synonymy network (N-connectivity condensed)")
        print(f"  {OUT_SYN_GROUPS}        — per-allele synonymy group membership")
    sys.exit(0)

df_cx = pd.read_csv(CROSS_TSV, sep="\t", encoding="utf-8-sig")
print(f"\nLoaded {len(df_cx)} cross records from {CROSS_TSV}")


def classify_cross(row):
    """Assign W/N/P based on HV functional groups."""
    mother   = str(row.get("Mother", "")).strip()
    father   = str(row.get("Father", "")).strip()
    allele_m = aaaa_allele.get(mother)
    allele_f = aaaa_allele.get(father)

    if allele_m is None:
        rows_m = df_allele_table[df_allele_table["Individual"] == mother]
        if not rows_m.empty:
            allele_m = rows_m.sort_values("Count", ascending=False).iloc[0]["Allele"]
    if allele_f is None:
        rows_f = df_allele_table[df_allele_table["Individual"] == father]
        if not rows_f.empty:
            allele_f = rows_f.sort_values("Count", ascending=False).iloc[0]["Allele"]

    if allele_m is None or allele_f is None:
        return "?", "?", "?"

    fg_m = allele_to_fg.get(allele_m, "?")
    fg_f = allele_to_fg.get(allele_f, "?")

    idx_m = allele_names.index(allele_m) if allele_m in allele_names else -1
    idx_f = allele_names.index(allele_f) if allele_f in allele_names else -1
    d = dist_hv[idx_m, idx_f] if idx_m >= 0 and idx_f >= 0 else np.nan

    if fg_m != fg_f:
        return CAT_P_CROSS, allele_m, allele_f
    elif np.isnan(d) or d == 0.0:
        return CAT_W, allele_m, allele_f
    elif d < WITHIN_CLASS_THRESHOLD:
        return CAT_N, allele_m, allele_f
    else:
        return CAT_P_WITHIN, allele_m, allele_f


cats, alleles_m, alleles_f = zip(*df_cx.apply(classify_cross, axis=1))
df_cx["Category"]    = cats
df_cx["Allele_Mother"] = alleles_m
df_cx["Allele_Father"] = alleles_f

# Identify seed yield column
seed_col = None
for col in ["Seeds", "seed_count", "germplasmQuantityEstimate", "seeds"]:
    if col in df_cx.columns:
        seed_col = col
        break
if seed_col is None:
    numeric_cols = df_cx.select_dtypes(include="number").columns.tolist()
    seed_col = numeric_cols[0] if numeric_cols else None

if seed_col is None:
    print("WARNING: No seed count column found. Cannot run result analysis.")
    sys.exit(0)

df_cx[seed_col] = pd.to_numeric(df_cx[seed_col], errors="coerce")
print(f"Seed yield column: '{seed_col}'")

df_known  = df_cx[df_cx["Category"] != "?"].copy()
n_unknown = (df_cx["Category"] == "?").sum()
if n_unknown:
    print(f"Note: {n_unknown} cross(es) excluded — individual not in allele table")

print("\n── Seed yield by cross category ──")
print(f"  {'Category':<8} {'N':>6} {'N seeds>0':>10} {'Mean':>8} {'Median':>8} {'% success':>10}")
cat_summary = []
all_cats = [CAT_W, CAT_N, CAT_P_WITHIN, CAT_P_CROSS]
for cat in all_cats:
    sub    = df_known[df_known["Category"] == cat][seed_col].dropna()
    n_with = int((sub > 0).sum())
    pct    = round(100 * n_with / len(sub), 1) if len(sub) > 0 else np.nan
    cat_summary.append({
        "Category": cat, "N_crosses": len(sub), "N_with_seeds": n_with,
        "Pct_successful": pct,
        "Mean_seeds": round(sub.mean(), 2) if len(sub) > 0 else np.nan,
        "Median_seeds": round(sub.median(), 2) if len(sub) > 0 else np.nan,
    })
    print(f"  {cat:<10} {len(sub):>6} {n_with:>10} "
          f"{sub.mean() if len(sub) else float('nan'):>8.2f} "
          f"{sub.median() if len(sub) else float('nan'):>8.2f} {pct:>9.1f}%")

# Statistical tests
groups = [df_known[df_known["Category"] == cat][seed_col].dropna().values
          for cat in all_cats]
groups_valid = [g for g in groups if len(g) >= 2]

if len(groups_valid) >= 2:
    from scipy.stats import kruskal, mannwhitneyu
    stat, pval = kruskal(*groups_valid)
    print(f"\nKruskal-Wallis across categories: H = {stat:.3f}, p = {pval:.4f}")
    print("Pairwise Mann-Whitney U (two-sided):")
    for ii in range(len(all_cats)):
        for jj in range(ii + 1, len(all_cats)):
            g1 = df_known[df_known["Category"] == all_cats[ii]][seed_col].dropna()
            g2 = df_known[df_known["Category"] == all_cats[jj]][seed_col].dropna()
            if len(g1) >= 2 and len(g2) >= 2:
                u, p = mannwhitneyu(g1, g2, alternative="two-sided")
                sig  = "**" if p < 0.01 else ("*" if p < 0.05 else "ns")
                print(f"  {all_cats[ii]} vs {all_cats[jj]}: "
                      f"U = {u:.1f}, p = {p:.4f} {sig}")

# Results figure
print(f"\nWriting results figure to {OUT_RESULTS_FIG} ...")
cat_colours = {
    CAT_W:        "#d62728",   # red
    CAT_N:        "#ff7f0e",   # orange
    CAT_P_WITHIN: "#1f77b4",   # blue
    CAT_P_CROSS:  "#2ca02c",   # green
}
cat_labels_full = {
    CAT_W:        "W — HV-identical\n(expected incompatible)",
    CAT_N:        "N — Small HV diff\n(synonymy test)",
    CAT_P_WITHIN: "P_within — Same class, distinct HV\n(predicted compatible)",
    CAT_P_CROSS:  "P_cross — Between class\n(guaranteed compatible)",
}

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
ax_box, ax_bar = axes

plot_data = [df_known[df_known["Category"] == cat][seed_col].dropna().values
             for cat in all_cats]
plot_cols = [cat_colours[c] for c in all_cats]

bp = ax_box.boxplot(plot_data, positions=range(len(all_cats)), patch_artist=True,
                    widths=0.4, showfliers=False,
                    medianprops=dict(color="white", linewidth=2))
for patch, col in zip(bp["boxes"], plot_cols):
    patch.set_facecolor(col)
    patch.set_alpha(0.6)
for w in bp["whiskers"] + bp["caps"]:
    w.set_color("grey")

rng = np.random.default_rng(42)
for i, (vals, col) in enumerate(zip(plot_data, plot_cols)):
    jitter = rng.uniform(-0.12, 0.12, size=len(vals))
    ax_box.scatter(np.full(len(vals), i) + jitter, vals,
                   color=col, alpha=0.7, s=25, zorder=3, edgecolors="none")

ax_box.set_xticks(range(len(all_cats)))
ax_box.set_xticklabels([f"{c}\n(n={len(d)})" for c, d in zip(all_cats, plot_data)],
                       fontsize=9)
ax_box.set_ylabel("Seed yield", fontsize=11)
ax_box.set_title("Seed yield by cross category", fontsize=12)

cats_bar = [r["Category"] for r in cat_summary if r["N_crosses"] > 0]
pcts_bar = [r["Pct_successful"] for r in cat_summary if r["N_crosses"] > 0]
ns_bar   = [r["N_crosses"] for r in cat_summary if r["N_crosses"] > 0]

ax_bar.bar(range(len(cats_bar)), pcts_bar,
           color=[cat_colours[c] for c in cats_bar], alpha=0.8, edgecolor="white")
ax_bar.set_xticks(range(len(cats_bar)))
ax_bar.set_xticklabels([f"{c}\n(n={n})" for c, n in zip(cats_bar, ns_bar)], fontsize=9)
ax_bar.set_ylabel("Successful crosses (%)", fontsize=11)
ax_bar.set_title("Cross success rate by category", fontsize=12)
ax_bar.set_ylim(0, 110)
for i, pct in enumerate(pcts_bar):
    if not np.isnan(pct):
        ax_bar.text(i, pct + 2, f"{pct:.0f}%", ha="center", va="bottom", fontsize=10)

legend_patches = [mpatches.Patch(color=cat_colours[c], alpha=0.8,
                                  label=cat_labels_full[c])
                  for c in all_cats]
fig.legend(handles=legend_patches, loc="lower center", ncol=3,
           fontsize=8, bbox_to_anchor=(0.5, -0.08))

plt.tight_layout()
plt.savefig(OUT_RESULTS_FIG, format="pdf", dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved {OUT_RESULTS_FIG}")

print("\n── DONE ──")
