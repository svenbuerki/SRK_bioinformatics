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
