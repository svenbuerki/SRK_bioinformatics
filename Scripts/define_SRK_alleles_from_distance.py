#!/usr/bin/env python3
"""
define_SRK_alleles_from_distance.py

Groups SRK functional protein sequences into allele bins using pairwise amino
acid p-distance and average-linkage hierarchical clustering.

Rationale: proteins that share >X% identity in the ectodomain (S-domain) are
likely to present the same recognition specificity in crossing experiments.
Computing distance on the S-domain only avoids dilution by the conserved kinase.

Outputs
-------
  SRK_protein_allele_assignments.tsv   — protein → allele mapping
  SRK_protein_allele_representatives.fasta — one sequence per allele
  SRK_protein_distance_analysis.pdf   — sensitivity curve + distance heatmap

Usage
-----
    python define_SRK_alleles_from_distance.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from Bio import SeqIO
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

# ─────────────────────────────────────────────────────────────────────────────
# User settings
# ─────────────────────────────────────────────────────────────────────────────

INPUT_FASTA    = "SRK_functional_proteins_aligned.fasta"

# ── Clustering mode: choose ONE of the two options below ──────────────────
#
# Option A — fix the NUMBER OF ALLELES (recommended when you can see an elbow
#             on the sensitivity curve but cannot read the exact distance).
#             Set N_ALLELES to the desired count; DIST_THRESHOLD is ignored.
N_ALLELES = 58   # e.g. 50  — set to None to use Option B instead
#
# Option B — fix the DISTANCE THRESHOLD (p-distance = fraction of differing AA).
#             Used only when N_ALLELES is None.
#             Lower = stricter (more alleles); higher = more lumping (fewer alleles).
DIST_THRESHOLD = 0.01   # e.g. 0.01 = 1 %
# ──────────────────────────────────────────────────────────────────────────

# Alignment columns to use for distance (1-based, inclusive).
# Restrict to S-domain ectodomain to focus on the specificity-determining region.
# Approximate positions for Brassicaceae SRK (~858 aa protein):
#   Signal peptide : 1-30
#   S-domain       : 31-430   ← key for allele specificity
#   Transmembrane  : 431-455
#   Kinase domain  : 456-858
# Set to None to use the full alignment.
DOMAIN_REGION  = (31, 430)   # (start, end) 1-based inclusive

OUT_TSV        = "SRK_protein_allele_assignments.tsv"
OUT_FASTA      = "SRK_protein_allele_representatives.fasta"
OUT_PDF        = "SRK_protein_distance_analysis.pdf"

# Sensitivity analysis range.
# Start at 0 so the curve shows the full picture from exact duplicates upward.
THRESH_MIN, THRESH_MAX, THRESH_STEP = 0.0, 0.05, 0.0005

# ─────────────────────────────────────────────────────────────────────────────
# Load alignment
# ─────────────────────────────────────────────────────────────────────────────

records = list(SeqIO.parse(INPUT_FASTA, "fasta"))
names   = [r.id for r in records]
seqs    = [str(r.seq).upper() for r in records]
n       = len(seqs)
aln_len = len(seqs[0])

print(f"Loaded {n} sequences, alignment length {aln_len}")

# ─────────────────────────────────────────────────────────────────────────────
# Optionally restrict to S-domain columns
# ─────────────────────────────────────────────────────────────────────────────

if DOMAIN_REGION is not None:
    start, end = DOMAIN_REGION
    start = max(0, start - 1)          # convert to 0-based
    end   = min(aln_len, end)
    seqs_for_dist = [s[start:end] for s in seqs]
    print(f"Computing distances on alignment columns {start+1}-{end} "
          f"({end - start} positions, S-domain)")
else:
    seqs_for_dist = seqs
    print("Computing distances on full alignment")

# ─────────────────────────────────────────────────────────────────────────────
# Pairwise p-distance (gaps excluded from denominator)
# ─────────────────────────────────────────────────────────────────────────────

def p_distance(s1, s2):
    """Fraction of differing positions, gaps excluded from both numerator and
    denominator. Returns NaN if no valid columns."""
    diff = valid = 0
    for a, b in zip(s1, s2):
        if a == "-" or b == "-":
            continue
        valid += 1
        if a != b:
            diff += 1
    return diff / valid if valid > 0 else np.nan

print("Computing pairwise distance matrix …")
dist = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        d = p_distance(seqs_for_dist[i], seqs_for_dist[j])
        dist[i, j] = d
        dist[j, i] = d

# ─────────────────────────────────────────────────────────────────────────────
# Hierarchical clustering at chosen threshold
# ─────────────────────────────────────────────────────────────────────────────

condensed = squareform(dist)
Z = linkage(condensed, method="average")   # UPGMA

if N_ALLELES is not None:
    clusters  = fcluster(Z, N_ALLELES, criterion="maxclust")
    n_alleles = len(set(clusters))
    print(f"\nN_ALLELES = {N_ALLELES}  →  {n_alleles} alleles obtained "
          f"(from {n} proteins)")
else:
    clusters  = fcluster(Z, DIST_THRESHOLD, criterion="distance")
    n_alleles = len(set(clusters))
    print(f"\nDIST_THRESHOLD = {DIST_THRESHOLD:.4f}  →  {n_alleles} alleles "
          f"(from {n} proteins)")

# ─────────────────────────────────────────────────────────────────────────────
# Sensitivity analysis — allele count across a range of thresholds
# ─────────────────────────────────────────────────────────────────────────────

thresholds = np.arange(THRESH_MIN, THRESH_MAX + THRESH_STEP, THRESH_STEP)
allele_counts = []
for t in thresholds:
    # fcluster requires t > 0; at t=0 use a tiny positive value to get the
    # same result (only exact duplicates merged)
    t_eff = t if t > 0 else 1e-12
    c = fcluster(Z, t_eff, criterion="distance")
    allele_counts.append(len(set(c)))

# Count of proteins that are S-domain identical to at least one other protein
# (these are already merged at threshold=0, which is why the curve starts below n)
n_unique_sdomain = allele_counts[0]   # allele count at threshold=0
n_sdomain_dupes  = n - n_unique_sdomain

# Derive the distance that corresponds to the chosen allele count.
# This is the smallest threshold at which allele_count <= n_alleles.
# Used for plot annotation regardless of which clustering mode was chosen.
allele_counts_arr = np.array(allele_counts)
idx_match = np.where(allele_counts_arr <= n_alleles)[0]
implied_threshold = float(thresholds[idx_match[0]]) if len(idx_match) > 0 else None

# ─────────────────────────────────────────────────────────────────────────────
# Allele names and cluster statistics
# ─────────────────────────────────────────────────────────────────────────────

allele_map = {names[i]: f"Allele_{clusters[i]:03d}" for i in range(n)}

# Representative = sequence with most central position (min mean dist to cluster)
cluster_ids = sorted(set(clusters))
representative = {}
for cid in cluster_ids:
    members = [i for i, c in enumerate(clusters) if c == cid]
    if len(members) == 1:
        representative[cid] = members[0]
    else:
        sub = dist[np.ix_(members, members)]
        mean_dists = sub.mean(axis=1)
        representative[cid] = members[int(np.argmin(mean_dists))]

# Cluster size summary
sizes = pd.Series(clusters).value_counts().sort_index()

# ─────────────────────────────────────────────────────────────────────────────
# Write TSV
# ─────────────────────────────────────────────────────────────────────────────

rows = []
for i, name in enumerate(names):
    cid     = clusters[i]
    allele  = allele_map[name]
    is_repr = (representative[cid] == i)
    rows.append({"Protein": name, "Allele": allele,
                 "Cluster_id": int(cid), "Is_representative": is_repr})

df_out = pd.DataFrame(rows)
df_out.to_csv(OUT_TSV, sep="\t", index=False)
print(f"Allele assignments written to {OUT_TSV}")

# ─────────────────────────────────────────────────────────────────────────────
# Write representative FASTA
# ─────────────────────────────────────────────────────────────────────────────

with open(OUT_FASTA, "w") as fh:
    for cid in cluster_ids:
        idx    = representative[cid]
        allele = f"Allele_{cid:03d}"
        seq_no_gaps = seqs[idx].replace("-", "")
        fh.write(f">{allele}  source={names[idx]}\n{seq_no_gaps}\n")
print(f"Representative sequences written to {OUT_FASTA}")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1: sensitivity curve
# ─────────────────────────────────────────────────────────────────────────────

# Order proteins by cluster for a clean heatmap
order = np.argsort(clusters)
dist_ordered = dist[np.ix_(order, order)]
labels_ordered = [names[i] for i in order]

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# — sensitivity analysis ——————————————————————
ax = axes[0]
ax.plot(thresholds, allele_counts, color="steelblue", linewidth=2)

# Reference line: total proteins in input
ax.axhline(n, color="grey", linestyle=":", linewidth=1.2,
           label=f"Total proteins in input (n={n})")

# Reference line: unique S-domain sequences (curve starting point)
ax.axhline(n_unique_sdomain, color="#888888", linestyle="--", linewidth=1,
           label=f"Unique S-domain sequences at threshold=0 (n={n_unique_sdomain})\n"
                 f"  → {n_sdomain_dupes} proteins share their S-domain with another protein\n"
                 f"     (differ only in kinase / other regions)")

# Chosen clustering — annotate differently depending on mode
if N_ALLELES is not None:
    # Horizontal line at the target allele count
    ax.axhline(n_alleles, color="tomato", linestyle="--", linewidth=1.5,
               label=f"N_ALLELES = {N_ALLELES}  →  {n_alleles} alleles obtained")
    # Vertical line at the implied distance threshold
    if implied_threshold is not None:
        ax.axvline(implied_threshold, color="tomato", linestyle=":",
                   linewidth=1.2,
                   label=f"Implied distance threshold ≈ {implied_threshold:.4f}")
        ax.text(implied_threshold, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else n,
                f" {implied_threshold:.4f}", color="tomato",
                fontsize=7, va="top", ha="left")
else:
    # Vertical line at the chosen distance threshold
    ax.axvline(DIST_THRESHOLD, color="tomato", linestyle="--", linewidth=1.5,
               label=f"Current DIST_THRESHOLD = {DIST_THRESHOLD:.4f}  →  {n_alleles} alleles\n"
                     f"  (look for the elbow in this curve to guide your choice)")

ax.set_xlabel("Distance threshold (p-distance)", fontsize=12)
ax.set_ylabel("Number of alleles", fontsize=12)
ax.set_title("Sensitivity: allele count vs. clustering threshold\n"
             "(S-domain p-distance, average linkage)", fontsize=12)
ax.legend(fontsize=8, loc="upper right")
ax.set_xlim(0, THRESH_MAX)

# — pairwise distance heatmap ——————————————————
ax = axes[1]

# Cluster boundary lines
boundary_pos = []
prev = clusters[order[0]]
for k, idx in enumerate(order):
    if clusters[idx] != prev:
        boundary_pos.append(k - 0.5)
        prev = clusters[idx]

im = ax.imshow(dist_ordered, aspect="auto", cmap="YlOrRd",
               vmin=0, vmax=dist.max())
plt.colorbar(im, ax=ax, label="p-distance")

for bp in boundary_pos:
    ax.axhline(bp, color="white", linewidth=0.5)
    ax.axvline(bp, color="white", linewidth=0.5)

ax.set_xticks([])
ax.set_yticks([])
mode_label = (f"N_ALLELES={N_ALLELES}" if N_ALLELES is not None
              else f"threshold={DIST_THRESHOLD:.4f}")
ax.set_title(f"Pairwise AA p-distance (S-domain)\nordered by cluster "
             f"({mode_label}, {n_alleles} alleles)",
             fontsize=11)
ax.set_xlabel(f"{n} proteins ordered by allele cluster", fontsize=10)
ax.set_ylabel(f"{n} proteins ordered by allele cluster", fontsize=10)

plt.tight_layout()
plt.savefig(OUT_PDF, format="pdf", dpi=150, bbox_inches="tight")
plt.close()
print(f"Distance analysis plot saved to {OUT_PDF}")

# ─────────────────────────────────────────────────────────────────────────────
# Summary to stdout
# ─────────────────────────────────────────────────────────────────────────────

print("\n── S-domain uniqueness ──")
print(f"  Total proteins       : {n}")
print(f"  Unique S-domains     : {n_unique_sdomain}  "
      f"({n_sdomain_dupes} proteins share their S-domain with ≥1 other protein)")

print("\n── Allele size distribution ──")
singleton = (sizes == 1).sum()
print(f"  Total alleles        : {n_alleles}")
print(f"  Singleton alleles    : {singleton}  "
      f"({100*singleton/n_alleles:.0f} % of alleles)")
print(f"  Max cluster size     : {sizes.max()}")
if N_ALLELES is not None:
    print(f"\n  Mode                 : N_ALLELES = {N_ALLELES}")
    if implied_threshold is not None:
        print(f"  Implied distance     : ≈ {implied_threshold:.4f}")
else:
    print(f"\n  Mode                 : DIST_THRESHOLD = {DIST_THRESHOLD}")
print(f"  Domain region        : {DOMAIN_REGION}")
