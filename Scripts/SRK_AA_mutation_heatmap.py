#!/usr/bin/env python3
"""
SRK_AA_mutation_heatmap.py

Produces two complementary heatmaps for SRK functional proteins, saved as
separate pages in a single PDF.

Figure 1 — Protein × position heatmap
  Rows    = proteins (ordered by allele cluster if assignment file present)
  Columns = variable alignment positions
  Color   = amino acid physicochemical class
  Header  = domain annotation band

Figure 2 — Amino acid frequency heatmap ("mutation rate" view)
  Rows    = 20 standard amino acids (grouped by physicochemical class)
  Columns = same variable alignment positions
  Color   = frequency of that AA at that position (0–1, white → dark red)
  Markers = black dot on the consensus (most frequent) AA per position
  Bottom  = Shannon entropy bar (position-level diversity)
  Header  = domain annotation band (same as Figure 1)

Outputs
-------
  SRK_AA_mutation_heatmap.pdf       — Figure 1: protein × position (AA class)
  SRK_AA_frequency_heatmap.pdf      — Figure 2: AA × position (frequency)
  SRK_AA_variable_positions.tsv     — variable positions with entropy & composition

Usage
-----
    python SRK_AA_mutation_heatmap.py

Dependencies: biopython, numpy, pandas, matplotlib
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm
from collections import Counter
from Bio import SeqIO
import os

# ─────────────────────────────────────────────────────────────────────────────
# User settings
# ─────────────────────────────────────────────────────────────────────────────

INPUT_FASTA      = "SRK_functional_proteins_aligned.fasta"

# Allele assignment file produced by define_SRK_alleles_from_distance.py.
# If present, rows are sorted by allele cluster. Set to None to sort by name.
ALLELE_TSV       = "SRK_protein_allele_assignments.tsv"

# Minimum gap frequency allowed at a position to be included (0-1).
# Positions where more than this fraction of sequences have a gap are excluded.
MAX_GAP_FREQ     = 0.20

# Only show variable positions (entropy > 0). Set ENTROPY_MIN > 0 to further
# restrict to positions above a minimum diversity level (0-2.77 bits for 20 AA).
ENTROPY_MIN      = 0.0

# Domain boundaries (1-based alignment positions, approximate).
# Adjust to match your species/amplicon if known from UniProt / literature.
# Format: list of (name, start, end, color)
DOMAIN_REGIONS = [
    ("Signal peptide",  1,   30,  "#b0c4de"),   # light steel blue
    ("S-domain",        31,  430, "#90ee90"),    # light green
    ("TM domain",       431, 455, "#ffa07a"),    # light salmon
    ("Kinase domain",   456, 848, "#d3a4d3"),    # light purple
]

OUT_PDF          = "SRK_AA_mutation_heatmap.pdf"
OUT_FREQ_PDF     = "SRK_AA_frequency_heatmap.pdf"
OUT_VAR_TSV      = "SRK_AA_variable_positions.tsv"

# Maximum number of variable positions to display (most variable first by entropy).
# Increase if you want to see more; a very wide figure may become hard to read.
MAX_POSITIONS    = 300

# Figure height per protein (inches). Reduce for many sequences.
ROW_HEIGHT       = 0.18
COL_WIDTH        = 0.10   # per position

# ─────────────────────────────────────────────────────────────────────────────
# Amino acid color scheme (physicochemical classes)
# ─────────────────────────────────────────────────────────────────────────────

AA_CLASS = {
    # Aliphatic / hydrophobic
    "G": "aliphatic", "A": "aliphatic", "V": "aliphatic",
    "L": "aliphatic", "I": "aliphatic", "M": "aliphatic",
    # Aromatic
    "F": "aromatic",  "Y": "aromatic",  "W": "aromatic",
    # Positively charged
    "K": "positive",  "R": "positive",  "H": "positive",
    # Negatively charged
    "D": "negative",  "E": "negative",
    # Polar uncharged
    "S": "polar",     "T": "polar",     "N": "polar",   "Q": "polar",
    # Special
    "C": "cysteine",
    "P": "proline",
    # Gap / unknown
    "-": "gap",       "X": "gap",
}

CLASS_COLORS = {
    "aliphatic": "#d3d3d3",   # light grey
    "aromatic":  "#ffd700",   # gold
    "positive":  "#6495ed",   # cornflower blue
    "negative":  "#ff6347",   # tomato
    "polar":     "#90ee90",   # light green
    "cysteine":  "#ffa500",   # orange
    "proline":   "#da70d6",   # orchid
    "gap":       "#ffffff",   # white
}

CLASS_ORDER = ["aliphatic", "aromatic", "positive", "negative",
               "polar", "cysteine", "proline", "gap"]
COLOR_LIST  = [CLASS_COLORS[c] for c in CLASS_ORDER]
CLASS_INDEX = {c: i for i, c in enumerate(CLASS_ORDER)}

# ─────────────────────────────────────────────────────────────────────────────
# Amino acid order for Figure 2 (grouped by physicochemical class)
# ─────────────────────────────────────────────────────────────────────────────

# 20 standard AAs in class groups; order matches CLASS_COLORS for consistency
AA_ORDER = list("GAVLIM" "FYW" "KRH" "DE" "STNQ" "C" "P")

# Class membership for each row (used to draw group dividers and row colors)
AA_ROW_CLASS = {aa: AA_CLASS.get(aa, "gap") for aa in AA_ORDER}

# Boundaries between physicochemical groups in AA_ORDER (indices after which
# a divider line is drawn: after aliphatic, aromatic, positive, negative, polar)
AA_GROUP_BREAKS = [5, 8, 11, 13, 17]   # 0-based last index of each group

# ─────────────────────────────────────────────────────────────────────────────
# Helper: draw the domain annotation track onto an axis
# ─────────────────────────────────────────────────────────────────────────────

def draw_domain_track(ax, pos_indices):
    """Draw domain color bands on a pre-configured axis (xlim already set)."""
    ax.set_ylim(0, 1)
    ax.axis("off")
    for dom_name, dom_start, dom_end, dom_color in DOMAIN_REGIONS:
        in_dom = [j for j, pos in enumerate(pos_indices)
                  if dom_start <= pos + 1 <= dom_end]
        if not in_dom:
            continue
        x0 = in_dom[0]  - 0.5
        x1 = in_dom[-1] + 0.5
        ax.add_patch(mpatches.FancyBboxPatch(
            (x0, 0.1), x1 - x0, 0.8,
            boxstyle="round,pad=0.01",
            facecolor=dom_color, edgecolor="white", linewidth=0.8
        ))
        if x1 - x0 > 5:
            ax.text((x0 + x1) / 2, 0.5, dom_name,
                    ha="center", va="center", fontsize=7,
                    fontweight="bold", color="black")

# ─────────────────────────────────────────────────────────────────────────────
# Load alignment
# ─────────────────────────────────────────────────────────────────────────────

records  = list(SeqIO.parse(INPUT_FASTA, "fasta"))
prot_ids = [r.id for r in records]
seqs     = [str(r.seq).upper() for r in records]
n        = len(seqs)
aln_len  = len(seqs[0])

print(f"Loaded {n} sequences, alignment length {aln_len}")

# ─────────────────────────────────────────────────────────────────────────────
# Row order: by allele cluster if available, otherwise by protein name
# ─────────────────────────────────────────────────────────────────────────────

if ALLELE_TSV and os.path.isfile(ALLELE_TSV):
    allele_df = pd.read_csv(ALLELE_TSV, sep="\t", encoding="utf-8-sig")
    allele_map = dict(zip(allele_df["Protein"], allele_df["Allele"]))
    # Sort by allele name (cluster), then protein name within cluster
    order = sorted(range(n), key=lambda i: (allele_map.get(prot_ids[i], "Z"),
                                             prot_ids[i]))
    row_labels  = [f"{allele_map.get(prot_ids[i], '?')} | {prot_ids[i]}"
                   for i in order]
    print(f"Row order: by allele cluster from {ALLELE_TSV}")
else:
    order      = list(range(n))
    row_labels = prot_ids
    print("Row order: by protein name (allele TSV not found)")

seqs_ordered = [seqs[i] for i in order]
ids_ordered  = [prot_ids[i] for i in order]

# ─────────────────────────────────────────────────────────────────────────────
# Identify variable positions
# ─────────────────────────────────────────────────────────────────────────────

def shannon_entropy(col):
    """Shannon entropy (bits) of AA composition at one alignment column."""
    counts = {}
    for aa in col:
        if aa not in ("-", "X"):
            counts[aa] = counts.get(aa, 0) + 1
    total = sum(counts.values())
    if total == 0:
        return 0.0
    H = 0.0
    for cnt in counts.values():
        p = cnt / total
        H -= p * np.log2(p)
    return H

variable_positions = []
for pos in range(aln_len):
    col      = [s[pos] for s in seqs]
    gap_freq = col.count("-") / n
    if gap_freq > MAX_GAP_FREQ:
        continue
    H = shannon_entropy(col)
    if H > ENTROPY_MIN:
        variable_positions.append((pos, H))

variable_positions.sort(key=lambda x: -x[1])  # sort by entropy descending
print(f"Variable positions (entropy>{ENTROPY_MIN}, gap<{MAX_GAP_FREQ}): "
      f"{len(variable_positions)}")

# Limit to top MAX_POSITIONS
if len(variable_positions) > MAX_POSITIONS:
    print(f"Showing top {MAX_POSITIONS} most variable positions")
    variable_positions = variable_positions[:MAX_POSITIONS]

# Re-sort by alignment position for left-to-right display
variable_positions.sort(key=lambda x: x[0])
pos_indices = [p for p, _ in variable_positions]
entropies   = [h for _, h in variable_positions]

# ─────────────────────────────────────────────────────────────────────────────
# Save variable positions table
# ─────────────────────────────────────────────────────────────────────────────

var_rows = []
for pos, H in zip(pos_indices, entropies):
    col = [s[pos] for s in seqs]
    aas = [aa for aa in col if aa not in ("-", "X")]
    aa_counts = Counter(aas)
    var_rows.append({
        "alignment_pos": pos + 1,
        "entropy_bits": round(H, 3),
        "n_variants": len(aa_counts),
        "aa_composition": ";".join(f"{aa}:{cnt}" for aa, cnt in
                                   sorted(aa_counts.items()))
    })

pd.DataFrame(var_rows).to_csv(OUT_VAR_TSV, sep="\t", index=False)
print(f"Variable positions table saved to {OUT_VAR_TSV}")

# ─────────────────────────────────────────────────────────────────────────────
# Build numeric matrix for heatmap
# ─────────────────────────────────────────────────────────────────────────────

n_pos = len(pos_indices)
matrix = np.zeros((n, n_pos), dtype=int)

for row_i, seq in enumerate(seqs_ordered):
    for col_j, pos in enumerate(pos_indices):
        aa = seq[pos]
        cls = AA_CLASS.get(aa, "gap")
        matrix[row_i, col_j] = CLASS_INDEX[cls]

# ─────────────────────────────────────────────────────────────────────────────
# Figure layout
# ─────────────────────────────────────────────────────────────────────────────

fig_height = max(6, n * ROW_HEIGHT + 2.5)   # +2.5 for domain track + colorbars
fig_width  = max(10, n_pos * COL_WIDTH + 3)
fig_width  = min(fig_width, 40)               # cap width

fig = plt.figure(figsize=(fig_width, fig_height))

# GridSpec: row 0 = domain annotation, row 1 = heatmap
gs = fig.add_gridspec(2, 1, height_ratios=[0.06, 1], hspace=0.02)
ax_dom  = fig.add_subplot(gs[0])
ax_heat = fig.add_subplot(gs[1])

# ── domain annotation track ──────────────────────────────────────────────────
ax_dom.set_xlim(0, n_pos)
draw_domain_track(ax_dom, pos_indices)

# ── heatmap ──────────────────────────────────────────────────────────────────
cmap = ListedColormap(COLOR_LIST)
norm = BoundaryNorm(boundaries=np.arange(-0.5, len(CLASS_ORDER) + 0.5, 1),
                    ncolors=len(CLASS_ORDER))

ax_heat.imshow(matrix, aspect="auto", cmap=cmap, norm=norm,
               interpolation="nearest")

# x-axis: alignment position labels (every ~25 ticks to avoid crowding)
tick_every = max(1, n_pos // 20)
xticks = list(range(0, n_pos, tick_every))
ax_heat.set_xticks(xticks)
ax_heat.set_xticklabels([str(pos_indices[j] + 1) for j in xticks],
                        fontsize=6, rotation=90)
ax_heat.set_xlabel("Alignment position (variable only)", fontsize=10)

# y-axis: protein labels
ax_heat.set_yticks(range(n))
ax_heat.set_yticklabels(row_labels, fontsize=5)
ax_heat.set_ylabel("Protein (allele | id)", fontsize=10)

# entropy strip along x-axis (small bar below heatmap within same axis space)
# — skip to keep figure clean; entropy saved to TSV instead

ax_heat.set_title(
    f"SRK amino acid variation heatmap  "
    f"({n} proteins × {n_pos} variable positions)",
    fontsize=11, pad=8
)

# ── AA class legend ───────────────────────────────────────────────────────────
legend_handles = [
    mpatches.Patch(facecolor=CLASS_COLORS[c], edgecolor="grey",
                   linewidth=0.5, label=c.capitalize())
    for c in CLASS_ORDER
]
ax_heat.legend(handles=legend_handles, title="AA class",
               loc="lower right", bbox_to_anchor=(1.0, -0.12),
               ncol=len(CLASS_ORDER), fontsize=7, title_fontsize=8,
               framealpha=0.9)

plt.savefig(OUT_PDF, format="pdf", dpi=150, bbox_inches="tight")
plt.close()
print(f"Heatmap (Figure 1) saved to {OUT_PDF}")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2: per-AA frequency heatmap
# Rows = 20 standard AAs (grouped by class), columns = variable positions,
# values = frequency of that AA at that position among non-gap sequences.
# The most frequent (consensus) AA per position is marked with a black dot.
# A Shannon entropy bar at the bottom shows per-position diversity.
# ─────────────────────────────────────────────────────────────────────────────

n_aa = len(AA_ORDER)

# Build frequency matrix (n_aa × n_pos) using all sequences (not just ordered)
freq_matrix = np.zeros((n_aa, n_pos))
for col_j, pos in enumerate(pos_indices):
    col   = [s[pos] for s in seqs]          # use original seqs, order irrelevant
    valid = [aa for aa in col if aa not in ("-", "X")]
    total = len(valid)
    if total == 0:
        continue
    counts = Counter(valid)
    for row_i, aa in enumerate(AA_ORDER):
        freq_matrix[row_i, col_j] = counts.get(aa, 0) / total

# Consensus row index per position (highest frequency AA)
consensus_row = np.argmax(freq_matrix, axis=0)

# ── figure layout: domain track / heatmap / entropy bar ──────────────────────
fig2_height = max(6, n_aa * 0.35 + 3.0)
fig2 = plt.figure(figsize=(fig_width, fig2_height))

gs2 = fig2.add_gridspec(
    3, 1,
    height_ratios=[0.05, 1, 0.08],
    hspace=0.04
)
ax2_dom    = fig2.add_subplot(gs2[0])
ax2_freq   = fig2.add_subplot(gs2[1])
ax2_entrop = fig2.add_subplot(gs2[2])

# ── domain annotation track ───────────────────────────────────────────────────
ax2_dom.set_xlim(-0.5, n_pos - 0.5)
draw_domain_track(ax2_dom, pos_indices)

# ── frequency heatmap ─────────────────────────────────────────────────────────
im2 = ax2_freq.imshow(
    freq_matrix,
    aspect="auto",
    cmap="YlOrRd",
    vmin=0, vmax=1,
    interpolation="nearest"
)

# Mark consensus AA per position with a filled black circle
for col_j in range(n_pos):
    row_i = consensus_row[col_j]
    if freq_matrix[row_i, col_j] > 0:
        ax2_freq.plot(col_j, row_i, "k.", markersize=4, alpha=0.7)

# Horizontal dividers between AA physicochemical groups
for brk in AA_GROUP_BREAKS:
    ax2_freq.axhline(brk + 0.5, color="white", linewidth=1.2)

# Left-side class color strip: thin colored rectangles outside the main axis
for row_i, aa in enumerate(AA_ORDER):
    cls   = AA_ROW_CLASS[aa]
    color = CLASS_COLORS[cls]
    ax2_freq.add_patch(mpatches.Rectangle(
        (-1.6, row_i - 0.5), 0.9, 1.0,
        transform=ax2_freq.transData,
        facecolor=color, edgecolor="none", clip_on=False
    ))

# y-axis: single-letter AA labels
ax2_freq.set_yticks(range(n_aa))
ax2_freq.set_yticklabels(list(AA_ORDER), fontsize=8, fontfamily="monospace")
ax2_freq.set_ylabel("Amino acid", fontsize=10)
ax2_freq.set_xlim(-0.5, n_pos - 0.5)

# x-axis: position labels (suppress — shared with entropy bar below)
ax2_freq.set_xticks([])

ax2_freq.set_title(
    "SRK amino acid frequency at variable positions",
    fontsize=11, pad=6
)

# Colorbar
cbar = plt.colorbar(im2, ax=ax2_freq, fraction=0.015, pad=0.01)
cbar.set_label("Frequency", fontsize=9)
cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])

# ── Shannon entropy bar ───────────────────────────────────────────────────────
ax2_entrop.bar(range(n_pos), entropies, color="steelblue",
               width=0.9, linewidth=0)
ax2_entrop.set_xlim(-0.5, n_pos - 0.5)
ax2_entrop.set_ylim(0, max(entropies) * 1.15)
ax2_entrop.set_xticks(xticks)
ax2_entrop.set_xticklabels([str(pos_indices[j] + 1) for j in xticks],
                             fontsize=6, rotation=90)
ax2_entrop.set_xlabel("Alignment position (variable only)", fontsize=10)
ax2_entrop.set_ylabel("H (bits)", fontsize=7)
ax2_entrop.tick_params(axis="y", labelsize=6)

plt.savefig(OUT_FREQ_PDF, format="pdf", dpi=150, bbox_inches="tight")
plt.close()
print(f"Frequency heatmap (Figure 2) saved to {OUT_FREQ_PDF}")

print(f"\nDone. Figure 1: {n} proteins × {n_pos} positions. "
      f"Figure 2: {n_aa} AAs × {n_pos} positions.")
