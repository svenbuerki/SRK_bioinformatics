#!/usr/bin/env python3
"""
test_allele_definitions_from_crosses.py

Two-part module for Objective 2 — validating sequence-based S-allele bins
using controlled crosses between AAAA individuals.

PART 1 — Allele super-group clustering
  Groups the allele bins (defined in Obj1 by define_SRK_alleles_from_distance.py)
  into super-groups via a second-level UPGMA clustering on S-domain p-distances
  between allele representative sequences.  Each AAAA individual is assigned a
  bin + super-group label, and all pairwise AAAA combinations are ranked into
  three cross categories:
    W — within-bin    : same allele bin; expected incompatible (SI rejection)
    N — within-cluster: different bin, same super-group; hypothesis test
                        (are closely related bins truly distinct specificities?)
    P — between-cluster: different super-group; expected compatible (positive control)

PART 2 — Cross result analysis
  Reads completed cross records and tests whether the W / N / P category
  predicts seed yield.  Produces a summary table and figure.

Rationale
---------
Sequence-based allele bins are hypotheses about recognition specificity.  AAAA
plants are the cleanest test subjects: all four gene copies carry the same SRK
allele, so pollen identity is unambiguous.  Within-bin crosses (W) should fail;
between-cluster crosses (P) should succeed.  Within-cluster crosses (N) probe
the allele-bin boundaries — closely related bins could share recognition
specificity (synonymous alleles) or could be truly distinct.

Usage
-----
    cd /path/to/Canu_amplicon
    python test_allele_definitions_from_crosses.py

Outputs
-------
    SRK_allele_supergroups.tsv        allele bin → super-group + AAAA count
    SRK_AAAA_cross_design.tsv         all AAAA × AAAA pairs ranked by cross category
    SRK_allele_cluster_figure.pdf     dendrogram of allele bins + AAAA availability
    SRK_cross_result_analysis.pdf     seed yield by cross category (requires cross data)
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
from Bio import SeqIO

# ─────────────────────────────────────────────────────────────────────────────
# Settings — adjust here before running
# ─────────────────────────────────────────────────────────────────────────────

# Input files (paths relative to Canu_amplicon/)
REPS_FASTA   = "SRK_protein_allele_representatives.fasta"  # one seq per allele bin
GFS_TSV      = "SRK_individual_GFS.tsv"          # individual GFS and genotype pattern
ALLELE_TABLE = "SRK_individual_allele_table.tsv"  # individual → protein → allele → count
CROSS_TSV    = None   # set to "SRK_cross_compatibility.tsv" when ready to analyse results

# Super-group clustering
# Choose N_SUPERGROUPS (fixed count) OR SUPERGROUP_THRESHOLD (distance cutoff).
# N_SUPERGROUPS takes priority when set to an integer.
N_SUPERGROUPS        = 8     # target number of super-groups; set to None to use threshold
SUPERGROUP_THRESHOLD = 0.10  # p-distance used only when N_SUPERGROUPS is None

# S-domain region for distance calculation (1-based, inclusive).
# Must match the region used in define_SRK_alleles_from_distance.py.
DOMAIN_REGION = (31, 430)

# Cross category labels
CAT_W = "W"   # within-bin        — same allele, expected incompatible
CAT_N = "N"   # within-cluster    — different bin, same super-group, hypothesis test
CAT_P = "P"   # between-cluster   — different super-group, expected compatible

# Output files
OUT_SUPERGROUPS  = "SRK_allele_supergroups.tsv"
OUT_CROSS_DESIGN = "SRK_AAAA_cross_design.tsv"
OUT_CLUSTER_FIG  = "SRK_allele_cluster_figure.pdf"
OUT_RESULTS_FIG  = "SRK_cross_result_analysis.pdf"

# ─────────────────────────────────────────────────────────────────────────────
# PART 1 — Allele super-group clustering
# ─────────────────────────────────────────────────────────────────────────────

print("=" * 60)
print("PART 1 — Allele super-group clustering")
print("=" * 60)

# ── 1a. Load allele representative sequences ──────────────────────────────────

records = list(SeqIO.parse(REPS_FASTA, "fasta"))
# Representative headers: ">Allele_XXX  source=..."  — take first token as allele name
allele_names = [r.id.split()[0] for r in records]
seqs         = [str(r.seq).upper() for r in records]
n_alleles    = len(allele_names)
aln_len      = len(seqs[0])

print(f"\nLoaded {n_alleles} allele representative sequences from {REPS_FASTA}")

# ── 1b. Restrict to S-domain ──────────────────────────────────────────────────

start_0      = max(0, DOMAIN_REGION[0] - 1)   # convert to 0-based
end_0        = min(aln_len, DOMAIN_REGION[1])
seqs_sdomain = [s[start_0:end_0] for s in seqs]
print(f"Using S-domain columns {DOMAIN_REGION[0]}–{DOMAIN_REGION[1]} "
      f"({end_0 - start_0} positions)")

# ── 1c. Pairwise p-distance between allele representatives ───────────────────

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

print("Computing pairwise distances between allele representatives …")
dist = np.zeros((n_alleles, n_alleles))
for i in range(n_alleles):
    for j in range(i + 1, n_alleles):
        d = p_distance(seqs_sdomain[i], seqs_sdomain[j])
        dist[i, j] = d
        dist[j, i] = d

# Replace any NaN with the column mean (rare; safeguard only)
col_means = np.nanmean(dist, axis=0)
for i in range(n_alleles):
    for j in range(n_alleles):
        if np.isnan(dist[i, j]):
            dist[i, j] = (col_means[i] + col_means[j]) / 2

# Print distance summary to help choose N_SUPERGROUPS
print(f"\nS-domain p-distance summary across all {n_alleles}×{n_alleles} allele pairs:")
upper = dist[np.triu_indices(n_alleles, k=1)]
print(f"  Min    : {upper.min():.4f}")
print(f"  Median : {np.median(upper):.4f}")
print(f"  Mean   : {upper.mean():.4f}")
print(f"  Max    : {upper.max():.4f}")
print(f"  Pairs with dist < 0.05 : {(upper < 0.05).sum()}")
print(f"  Pairs with dist < 0.10 : {(upper < 0.10).sum()}")
print(f"  Pairs with dist < 0.20 : {(upper < 0.20).sum()}")

# ── 1d. Second-level hierarchical clustering into super-groups ────────────────

condensed = squareform(dist)
Z = linkage(condensed, method="average")   # UPGMA

if N_SUPERGROUPS is not None:
    labels   = fcluster(Z, N_SUPERGROUPS, criterion="maxclust")
    n_groups = len(set(labels))
    mode_label = f"N_SUPERGROUPS = {N_SUPERGROUPS} → {n_groups} groups obtained"
else:
    labels   = fcluster(Z, SUPERGROUP_THRESHOLD, criterion="distance")
    n_groups = len(set(labels))
    mode_label = f"SUPERGROUP_THRESHOLD = {SUPERGROUP_THRESHOLD} → {n_groups} groups"

print(f"\nSuper-group clustering: {mode_label}")

allele_to_supergroup = {allele_names[i]: f"SG{labels[i]:02d}"
                        for i in range(n_alleles)}

# ── 1e. Load AAAA individual data ─────────────────────────────────────────────

df_gfs = pd.read_csv(GFS_TSV, sep="\t")
aaaa_individuals = set(df_gfs.loc[df_gfs["Genotype_Pattern"] == "AAAA", "Individual"])
print(f"\nAAA individuals in GFS data: {len(aaaa_individuals)}")

# Map each AAAA individual to their single allele bin
df_allele_table = pd.read_csv(ALLELE_TABLE, sep="\t")
aaaa_allele = {}   # individual → allele bin
for ind in aaaa_individuals:
    rows = df_allele_table[df_allele_table["Individual"] == ind]
    if rows.empty:
        continue
    bins = rows["Allele"].unique()
    if len(bins) == 1:
        aaaa_allele[ind] = bins[0]
    else:
        # AAAA should carry one bin; if multiple appear take the most frequent
        counts = rows.groupby("Allele")["Count"].sum()
        aaaa_allele[ind] = counts.idxmax()

print(f"AAAA individuals with allele assignment: {len(aaaa_allele)}")

# ── 1f. Count AAAA individuals per allele bin ─────────────────────────────────

aaaa_count = {}
for ind, allele in aaaa_allele.items():
    aaaa_count[allele] = aaaa_count.get(allele, 0) + 1

# ── 1g. Build allele super-group table ───────────────────────────────────────

rows_sg = []
for allele in allele_names:
    sg     = allele_to_supergroup[allele]
    n_aaaa = aaaa_count.get(allele, 0)
    rows_sg.append({
        "Allele":     allele,
        "SuperGroup": sg,
        "N_AAAA":     n_aaaa,
        "CrossPower": "full" if n_aaaa >= 2 else ("singleton" if n_aaaa == 1 else "none"),
    })

df_sg = pd.DataFrame(rows_sg).sort_values(["SuperGroup", "Allele"])
df_sg.to_csv(OUT_SUPERGROUPS, sep="\t", index=False)
print(f"\nAllele super-group table written to {OUT_SUPERGROUPS}")

print("\n── Super-group summary ──")
for sg, grp in df_sg.groupby("SuperGroup"):
    n_aaaa_sg    = grp["N_AAAA"].sum()
    alleles_list = ", ".join(grp["Allele"].tolist())
    print(f"  {sg}: {len(grp):2d} allele(s), {n_aaaa_sg:3d} AAAA individuals")
    print(f"       {alleles_list}")

print(f"\n── AAAA cross power ──")
for power, grp in df_sg.groupby("CrossPower"):
    if power == "full":
        note = "≥2 AAAA → full 3-tier design"
    elif power == "singleton":
        note = "1 AAAA → cross partner only"
    else:
        note = "no AAAA → must use AAAB/AABB"
    print(f"  {power:9s}: {len(grp):2d} allele bins  ({note})")

# ── 1h. Generate AAAA cross design matrix ─────────────────────────────────────

print("\nGenerating AAAA × AAAA cross design matrix …")

aaaa_list  = sorted(aaaa_allele.keys())
cross_rows = []

for i, ind_a in enumerate(aaaa_list):
    for ind_b in aaaa_list[i + 1:]:
        allele_a = aaaa_allele[ind_a]
        allele_b = aaaa_allele[ind_b]
        sg_a     = allele_to_supergroup.get(allele_a, "?")
        sg_b     = allele_to_supergroup.get(allele_b, "?")

        if allele_a == allele_b:
            category = CAT_W
            expected = "Incompatible — same SI specificity"
        elif sg_a == sg_b:
            category = CAT_N
            expected = "Test — closely related alleles; specificity unknown"
        else:
            category = CAT_P
            expected = "Compatible — positive control"

        idx_a = allele_names.index(allele_a) if allele_a in allele_names else -1
        idx_b = allele_names.index(allele_b) if allele_b in allele_names else -1
        d = dist[idx_a, idx_b] if idx_a >= 0 and idx_b >= 0 else np.nan

        cross_rows.append({
            "Mother":        ind_a,
            "Father":        ind_b,
            "Mother_allele": allele_a,
            "Father_allele": allele_b,
            "Mother_SG":     sg_a,
            "Father_SG":     sg_b,
            "S_domain_dist": round(d, 4),
            "Category":      category,
            "Expected":      expected,
        })

df_cross = pd.DataFrame(cross_rows)

cat_order    = {CAT_W: 0, CAT_N: 1, CAT_P: 2}
df_cross["_sort"] = df_cross["Category"].map(cat_order)
df_cross = df_cross.sort_values(["_sort", "S_domain_dist"]).drop(columns="_sort")
df_cross.to_csv(OUT_CROSS_DESIGN, sep="\t", index=False)
print(f"Cross design matrix written to {OUT_CROSS_DESIGN}")

n_w = (df_cross["Category"] == CAT_W).sum()
n_n = (df_cross["Category"] == CAT_N).sum()
n_p = (df_cross["Category"] == CAT_P).sum()
print(f"\n  W (within-bin, expected incompatible) : {n_w:4d} pairs")
print(f"  N (within-cluster, hypothesis test)   : {n_n:4d} pairs")
print(f"  P (between-cluster, positive control) : {n_p:4d} pairs")
print(f"  Total                                 : {len(df_cross):4d} pairs")

# ── 1i. Figure: dendrogram + AAAA availability ────────────────────────────────

print(f"\nWriting cluster figure to {OUT_CLUSTER_FIG} …")

sg_ids     = sorted(set(allele_to_supergroup.values()))
cmap_sg    = matplotlib.colormaps.get_cmap("tab10").resampled(len(sg_ids))
sg_colours = {sg: cmap_sg(i) for i, sg in enumerate(sg_ids)}

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

leaf_order   = dend["ivl"]
ytick_labels = ax_dendro.get_yticklabels()
for lbl in ytick_labels:
    allele = lbl.get_text()
    sg     = allele_to_supergroup.get(allele, "?")
    lbl.set_color(sg_colours.get(sg, "black"))
    lbl.set_fontsize(7)
ax_dendro.set_yticklabels(ytick_labels)
ax_dendro.set_xlabel("S-domain p-distance", fontsize=10)
ax_dendro.set_title("Allele bin super-group clustering\n(UPGMA, S-domain p-distance)",
                    fontsize=11)

sg_patches = [mpatches.Patch(color=sg_colours[sg], label=sg) for sg in sg_ids]
ax_dendro.legend(handles=sg_patches, loc="lower right", fontsize=8,
                 title="Super-group", ncol=2)

counts_ordered = [aaaa_count.get(a, 0) for a in leaf_order]
colours_bar    = [sg_colours.get(allele_to_supergroup.get(a, "?"), "grey")
                  for a in leaf_order]

ax_bar.barh(range(len(leaf_order)), counts_ordered,
            color=colours_bar, edgecolor="none", height=0.7)
ax_bar.axvline(2, color="black", linestyle="--", linewidth=0.8,
               label="Minimum for\nwithin-bin test (n=2)")
ax_bar.set_yticks([])
ax_bar.set_xlabel("AAAA individuals", fontsize=9)
ax_bar.set_title("AAAA\navailability", fontsize=10)
ax_bar.legend(fontsize=7, loc="lower right")
ax_bar.set_xlim(0, max(counts_ordered) * 1.15 if counts_ordered else 5)

plt.tight_layout()
plt.savefig(OUT_CLUSTER_FIG, format="pdf", dpi=150, bbox_inches="tight")

# Save PNG to figures/ folder for use in reports
os.makedirs("figures", exist_ok=True)
png_path = os.path.join("figures", os.path.splitext(OUT_CLUSTER_FIG)[0] + ".png")
plt.savefig(png_path, format="png", dpi=200, bbox_inches="tight")
plt.close()
print(f"Cluster figure saved ({OUT_CLUSTER_FIG}, {png_path})")

# ─────────────────────────────────────────────────────────────────────────────
# PART 2 — Cross result analysis
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("PART 2 — Cross result analysis")
print("=" * 60)

if not (CROSS_TSV and os.path.exists(CROSS_TSV)):
    print(f"\nCROSS_TSV is not set or file not found. Skipping Part 2.")
    print("Set CROSS_TSV = '<filename>' in settings when ready to analyse results.")
    sys.exit(0)

df_cx = pd.read_csv(CROSS_TSV, sep="\t")
print(f"\nLoaded {len(df_cx)} cross records from {CROSS_TSV}")

# ── 2a. Assign W / N / P to each cross ───────────────────────────────────────

def classify_cross(row):
    """Assign W/N/P using allele bin and super-group; falls back to most
    frequent allele for non-AAAA parents."""
    mother = str(row.get("Mother", "")).strip()
    father = str(row.get("Father", "")).strip()

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

    sg_m = allele_to_supergroup.get(allele_m, "?")
    sg_f = allele_to_supergroup.get(allele_f, "?")

    if allele_m == allele_f:
        return CAT_W, allele_m, allele_f
    elif sg_m == sg_f:
        return CAT_N, allele_m, allele_f
    else:
        return CAT_P, allele_m, allele_f

cats, alleles_m, alleles_f = zip(*df_cx.apply(classify_cross, axis=1))
df_cx["Category"]      = cats
df_cx["Allele_Mother"] = alleles_m
df_cx["Allele_Father"] = alleles_f
df_cx["SG_Mother"]     = [allele_to_supergroup.get(a, "?") for a in alleles_m]
df_cx["SG_Father"]     = [allele_to_supergroup.get(a, "?") for a in alleles_f]

# ── 2b. Seed yield summary by category ───────────────────────────────────────

seed_col = None
for col in ["Seeds", "germplasmQuantityEstimate", "seed_count", "seeds"]:
    if col in df_cx.columns:
        seed_col = col
        break
if seed_col is None:
    numeric_cols = df_cx.select_dtypes(include="number").columns.tolist()
    seed_col = numeric_cols[0] if numeric_cols else None

if seed_col is None:
    print("WARNING: No seed count column found in cross data. Cannot run Part 2.")
    sys.exit(0)

df_cx[seed_col] = pd.to_numeric(df_cx[seed_col], errors="coerce")
print(f"Using '{seed_col}' as seed yield column")

df_known  = df_cx[df_cx["Category"] != "?"].copy()
n_unknown = (df_cx["Category"] == "?").sum()
if n_unknown:
    print(f"Note: {n_unknown} cross(es) excluded — individual not in allele table")

print("\n── Seed yield by cross category ──")
print(f"  {'Category':<8} {'N':>6} {'N seeds>0':>10} {'Mean':>8} {'Median':>8}")
cat_summary = []
for cat in [CAT_W, CAT_N, CAT_P]:
    sub    = df_known[df_known["Category"] == cat][seed_col].dropna()
    n_with = int((sub > 0).sum())
    cat_summary.append({
        "Category":       cat,
        "N_crosses":      len(sub),
        "N_with_seeds":   n_with,
        "Pct_successful": round(100 * n_with / len(sub), 1) if len(sub) > 0 else np.nan,
        "Mean_seeds":     round(sub.mean(), 2) if len(sub) > 0 else np.nan,
        "Median_seeds":   round(sub.median(), 2) if len(sub) > 0 else np.nan,
    })
    print(f"  {cat:<8} {len(sub):>6} {n_with:>10} "
          f"{sub.mean():>8.2f} {sub.median():>8.2f}")

# ── 2c. Statistical tests ─────────────────────────────────────────────────────

groups = [df_known[df_known["Category"] == cat][seed_col].dropna().values
          for cat in [CAT_W, CAT_N, CAT_P]]
groups = [g for g in groups if len(g) >= 2]

if len(groups) >= 2:
    from scipy.stats import kruskal, mannwhitneyu
    stat, pval = kruskal(*groups)
    print(f"\nKruskal-Wallis: H = {stat:.3f}, p = {pval:.4f}")
    cat_labels = [CAT_W, CAT_N, CAT_P]
    print("Pairwise Mann-Whitney U (two-sided):")
    for i in range(len(cat_labels)):
        for j in range(i + 1, len(cat_labels)):
            g1 = df_known[df_known["Category"] == cat_labels[i]][seed_col].dropna()
            g2 = df_known[df_known["Category"] == cat_labels[j]][seed_col].dropna()
            if len(g1) >= 2 and len(g2) >= 2:
                u, p = mannwhitneyu(g1, g2, alternative="two-sided")
                sig = "**" if p < 0.01 else ("*" if p < 0.05 else "ns")
                print(f"  {cat_labels[i]} vs {cat_labels[j]}: "
                      f"U = {u:.1f}, p = {p:.4f} {sig}")
else:
    print("Insufficient data for statistical tests (need ≥2 crosses per category).")

# ── 2d. Figure: seed yield by cross category ─────────────────────────────────

print(f"\nWriting results figure to {OUT_RESULTS_FIG} …")

cat_colours = {CAT_W: "#d62728", CAT_N: "#ff7f0e", CAT_P: "#2ca02c"}
cat_labels_full = {
    CAT_W: "W — Within-bin\n(same allele, expected incompatible)",
    CAT_N: "N — Within-cluster\n(different bin, same SG; hypothesis test)",
    CAT_P: "P — Between-cluster\n(positive control, expected compatible)",
}

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
ax_box = axes[0]
ax_bar = axes[1]

plot_data = [df_known[df_known["Category"] == cat][seed_col].dropna().values
             for cat in [CAT_W, CAT_N, CAT_P]]
plot_cols = [cat_colours[c] for c in [CAT_W, CAT_N, CAT_P]]

bp = ax_box.boxplot(
    plot_data, positions=range(3), patch_artist=True,
    widths=0.4, showfliers=False,
    medianprops=dict(color="white", linewidth=2),
)
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

ax_box.set_xticks([0, 1, 2])
ax_box.set_xticklabels(
    [f"{c}\n(n={len(d)})" for c, d in zip([CAT_W, CAT_N, CAT_P], plot_data)],
    fontsize=10)
ax_box.set_ylabel("Seed yield", fontsize=11)
ax_box.set_title("Seed yield by cross category", fontsize=12)

cats_bar = [r["Category"] for r in cat_summary if r["N_crosses"] > 0]
pcts_bar = [r["Pct_successful"] for r in cat_summary if r["N_crosses"] > 0]
ns_bar   = [r["N_crosses"] for r in cat_summary if r["N_crosses"] > 0]

ax_bar.bar(range(len(cats_bar)), pcts_bar,
           color=[cat_colours[c] for c in cats_bar], alpha=0.8, edgecolor="white")
ax_bar.set_xticks(range(len(cats_bar)))
ax_bar.set_xticklabels(
    [f"{c}\n(n={n})" for c, n in zip(cats_bar, ns_bar)], fontsize=10)
ax_bar.set_ylabel("Successful crosses (%)", fontsize=11)
ax_bar.set_title("Cross success rate by category", fontsize=12)
ax_bar.set_ylim(0, 110)
for i, pct in enumerate(pcts_bar):
    if not np.isnan(pct):
        ax_bar.text(i, pct + 2, f"{pct:.0f}%", ha="center", va="bottom", fontsize=10)

legend_patches = [mpatches.Patch(color=cat_colours[c], alpha=0.8,
                                  label=cat_labels_full[c])
                  for c in [CAT_W, CAT_N, CAT_P]]
fig.legend(handles=legend_patches, loc="lower center", ncol=3,
           fontsize=8, bbox_to_anchor=(0.5, -0.08))

plt.tight_layout()
plt.savefig(OUT_RESULTS_FIG, format="pdf", dpi=150, bbox_inches="tight")
plt.close()
print("Results figure saved.")

print("\n── DONE ──")
print(f"  {OUT_SUPERGROUPS}   — allele bin super-group assignments")
print(f"  {OUT_CROSS_DESIGN}  — AAAA × AAAA cross design matrix")
print(f"  {OUT_CLUSTER_FIG}   — dendrogram + AAAA availability figure")
print(f"  {OUT_RESULTS_FIG}   — seed yield by cross category")
