#!/usr/bin/env python3
"""
SRK_AAAA_cascade_hypothesis.py

Conceptual figure showing the five-stage cascade hypothesis explaining
AAAA predominance (56%) in LEPA populations.

Output
------
    figures/SRK_AAAA_cascade_hypothesis.png  (300 dpi)

Usage
-----
    python SRK_AAAA_cascade_hypothesis.py
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch, Polygon
from matplotlib.lines import Line2D
import numpy as np
import os

np.random.seed(42)

# ─── colour constants ──────────────────────────────────────────────────────────
A_COL = '#D62728'   # dominant A haplotype — red
B_COL = '#1F77B4'   # blue
C_COL = '#2CA02C'   # green
D_COL = '#FF7F0E'   # orange
E_COL = '#9467BD'   # purple
F_COL = '#8C564B'   # brown
G_COL = '#E377C2'   # pink
H_COL = '#BCBD22'   # olive

ALL_COLS = [A_COL, B_COL, C_COL, D_COL, E_COL, F_COL, G_COL, H_COL]

# Stage colour triplets: (background, border, badge-fill)
STAGE_PALETTE = [
    ('#FFF8E1', '#F9A825', '#E65100'),   # amber   — Stage 1
    ('#E8F5E9', '#2E7D32', '#1B5E20'),   # green   — Stage 2
    ('#E3F2FD', '#1565C0', '#0D47A1'),   # blue    — Stage 3
    ('#FCE4EC', '#C62828', '#B71C1C'),   # red     — Stage 4
    ('#EDE7F6', '#6A1B9A', '#4A148C'),   # purple  — Stage 5
]

OUTCOME_BG   = '#FFEBEE'
OUTCOME_EDGE = '#C62828'
ARROW_COL    = '#37474F'

# ─── stage text content ────────────────────────────────────────────────────────
TITLES = [
    'Ancestral Bottleneck',
    'Balancing Selection Breaks Down',
    'Mate Limitation and Reproductive Skew',
    'Genetic Drift Overwhelms Balancing Selection',
    'Polyploid-Specific Self-Incompatibility Breakdown',
]

BODY = [
    ("A severe demographic bottleneck strips S-allele richness faster than neutral\n"
     "models predict. Each S-allele is individually rare under balancing selection and\n"
     "is easily lost when founder group size is small. LEPA retains only 9–23% of the\n"
     "estimated species optimum of 65 alleles; 60–89% of the allele pool has been\n"
     "permanently lost per population (Wright 1939; Schierup et al. 1997)."),

    ("Under balancing selection, rare S-alleles are protected because their carriers\n"
     "have more compatible mates. Below a diversity threshold, balancing selection\n"
     "breaks down: S-allele A becomes numerically dominant, producing proportionally\n"
     "more offspring, rare S-alleles cannot recover, and the collapse accelerates\n"
     "rather than self-correcting (Castric & Vekemans 2004; Vekemans & Slatkin 1994)."),

    ("As homozygous individuals (AAAA) become predominant, those carrying rare S-alleles\n"
     "find progressively fewer compatible mates — the self-incompatibility system excludes them.\n"
     "Mate availability collapses sharply in small self-incompatible populations (Byers & Meagher 1992).\n"
     "Lineages carrying rare S-alleles fail to set seed and are removed from the breeding\n"
     "pool, irrespective of selection (Willi et al. 2005; Elam et al. 2007)."),

    ("At N = 25–40 individuals per population, genetic drift overwhelms balancing selection.\n"
     "The effective S-allele number advantage under balancing selection collapses when\n"
     "diversity is already low (Ne/N = 0.41–0.55 in LEPA). Rare S-alleles are\n"
     "stochastically lost and cannot be recovered without external gene flow\n"
     "(Willi et al. 2005; Aguilar et al. 2006)."),

    ("In diploids, self-incompatibility (SI) blocks same-S-allele crosses. In tetraploids, AABB produces\n"
     "AB pollen (4 of 6) carrying both A-SCR and B-SCR proteins. Both proteins\n"
     "compete to bind their matching pistil S-receptors, producing\n"
     "a signal too weak to trigger rejection — SI fails, self-fertilisation\n"
     "occurs, and AAAB offspring accumulate toward full homozygosity (AAAA; Mable et al. 2004)."),
]

# ─── figure geometry (data units = inches) ────────────────────────────────────
FW, FH  = 15.0, 15.2   # figure dimensions in inches
BOX_H   = 2.31          # stage box height (was 3.6)
GAP     = 0.35          # gap between boxes / arrow space (was 0.42)
MX      = 0.30          # left/right margin
BOX_W   = FW - 2 * MX  # stage box width

# x split within box: badge | schematic | text
BADGE_X  = MX + 0.52           # badge circle centre x
SCH_X0   = MX + 0.95           # schematic left edge
SCH_W    = BOX_W * 0.40        # schematic width
TEXT_X0  = SCH_X0 + SCH_W + 0.25  # text left edge

# y positions of top edge of each stage box (top-down layout)
TITLE_Y  = FH - 0.40
STAGE_Y  = [TITLE_Y - 1.20 - i * (BOX_H + GAP) for i in range(5)]
OUTCOME_Y = STAGE_Y[4] - BOX_H - GAP - 0.05

# bottom panel: outcome + conservation side by side
HALF_W   = (BOX_W - 0.20) / 2   # each panel ~7.1 inches wide
CI_X     = MX + HALF_W + 0.20   # left edge of conservation panel

# ─── create figure ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(FW, FH), facecolor='white')

# Single master axes covering the whole figure
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, FW)
ax.set_ylim(0, FH)
ax.axis('off')
ax.set_facecolor('white')


# ─── helpers ──────────────────────────────────────────────────────────────────
def rounded_box(x, y, w, h, fc, ec, lw=2.2, radius=0.18, zorder=2):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f'round,pad={radius}',
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder
    ))


def badge(cx, cy, n, fc):
    ax.add_patch(Circle((cx, cy), 0.36, fc=fc, ec='white', lw=2.2, zorder=5))
    ax.text(cx, cy, str(n), ha='center', va='center',
            fontsize=16, fontweight='bold', color='white', zorder=6)


def down_arrow(x, y_top, y_bot):
    ax.annotate('', xy=(x, y_bot + 0.05), xytext=(x, y_top - 0.05),
                arrowprops=dict(
                    arrowstyle='->', color=ARROW_COL,
                    lw=2.5, mutation_scale=22
                ), zorder=10)


def inset(x0_data, y0_data, w_data, h_data):
    """Add an inset axes using data coordinates (converted to figure fraction)."""
    x0f = x0_data / FW
    y0f = y0_data / FH
    wf  = w_data  / FW
    hf  = h_data  / FH
    a = fig.add_axes([x0f, y0f, wf, hf])
    a.set_facecolor('none')
    for sp in a.spines.values():
        sp.set_visible(False)
    a.set_xticks([])
    a.set_yticks([])
    return a


# ─── TITLE ────────────────────────────────────────────────────────────────────
ax.text(FW / 2, TITLE_Y,
        'The Compatibility Collapse Cascade (C3) Hypothesis:\n'
        'A five-stage mechanism linking habitat fragmentation to reproductive failure in self-incompatible plants',
        ha='center', va='top', fontsize=13, fontweight='bold', color='#1A1A1A')

ax.text(FW / 2, TITLE_Y - 0.72,
        'Each stage is initiated by the stage above and amplifies the next, '
        'creating a self-reinforcing cascade from habitat fragmentation to reproductive failure.',
        ha='center', va='top', fontsize=9.5, color='#444444', style='italic')


# ─── STAGE BOXES ──────────────────────────────────────────────────────────────
def draw_stage(i):
    bg, edge, bfc = STAGE_PALETTE[i]
    y_top = STAGE_Y[i]
    y_bot = y_top - BOX_H

    # Background rectangle
    rounded_box(MX, y_bot, BOX_W, BOX_H, bg, edge)

    # Badge
    badge(BADGE_X, y_bot + BOX_H / 2, i + 1, bfc)

    # Stage title — left-aligned, flush with schematic left edge (beside badge)
    ax.text(SCH_X0, y_top - 0.24, TITLES[i],
            ha='left', va='top', fontsize=12, fontweight='bold',
            color='#1A1A1A', zorder=4)

    # Body text — dropped below title line to prevent overlap
    ax.text(TEXT_X0, y_top - 0.52, BODY[i],
            ha='left', va='top', fontsize=11, color='#333333',
            linespacing=1.50, zorder=4)

    # Return inset bounds for the schematic (leave 0.3 padding top/bottom)
    pad = 0.30
    return inset(SCH_X0, y_bot + pad, SCH_W, BOX_H - 2 * pad - 0.55)


# ─── SCHEMATIC 1: Bottleneck ──────────────────────────────────────────────────
ax1 = draw_stage(0)
ax1.set_xlim(0, 10)
ax1.set_ylim(0, 4)

# Pre-bottleneck cloud (8 alleles, 2 copies each)
pre_cols = ALL_COLS
pre_xy = [
    (0.6, 3.2), (1.3, 3.5), (0.5, 2.6), (1.5, 2.8),
    (0.8, 2.0), (1.6, 2.2), (0.4, 1.4), (1.4, 1.5),
]
for (px, py), col in zip(pre_xy, pre_cols):
    ax1.add_patch(Circle((px, py), 0.28, fc=col, ec='white', lw=1.2, zorder=3))
    ax1.text(px, py, list('ABCDEFGH')[pre_cols.index(col)],
             ha='center', va='center', fontsize=7, fontweight='bold',
             color='white', zorder=4)

# Funnel
funnel = Polygon(
    [(2.1, 3.8), (3.9, 3.8), (3.2, 2.1), (2.8, 2.1)],
    closed=True, fc='#B0BEC5', ec='#607D8B', lw=1.5, zorder=2, alpha=0.7
)
ax1.add_patch(funnel)
ax1.text(3.0, 2.9, 'BOTTLE-\nNECK', ha='center', va='center',
         fontsize=6.5, color='#263238', fontweight='bold')
ax1.annotate('', xy=(3.0, 1.8), xytext=(3.0, 2.0),
             arrowprops=dict(arrowstyle='->', color='#607D8B', lw=1.5))

# Post-bottleneck (3 circles: mostly A)
post = [(A_COL, 'A'), (A_COL, 'A'), (B_COL, 'B')]
post_xy = [(4.5, 3.0), (5.4, 2.4), (4.8, 1.7)]
for (px, py), (col, lbl) in zip(post_xy, post):
    ax1.add_patch(Circle((px, py), 0.32, fc=col, ec='white', lw=1.2, zorder=3))
    ax1.text(px, py, lbl, ha='center', va='center', fontsize=8,
             fontweight='bold', color='white', zorder=4)

# labels
ax1.text(1.0, 0.3, '8 alleles', ha='center', fontsize=8, color='#555')
ax1.text(5.0, 0.3, '2 alleles', ha='center', fontsize=8, color='#555')

# ─── SCHEMATIC 2: Balancing Selection Breakdown ───────────────────────────────
ax2 = draw_stage(1)
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 4.5)

alleles_eq  = list('ABCDE')
freqs_eq    = [1.0, 0.9, 1.1, 0.95, 1.0]
alleles_col = [A_COL, B_COL, C_COL, D_COL, E_COL]

freqs_col   = [3.8, 0.3, 0.2, 0.15, 0.1]

# "Equilibrium" bar chart (left)
x0_eq = 0.3
bw = 0.55
for j, (f, col) in enumerate(zip(freqs_eq, alleles_col)):
    ax2.bar(x0_eq + j * (bw + 0.15), f, width=bw, color=col,
            ec='white', linewidth=0.8, bottom=0)
ax2.text(1.6, -0.5, 'Bal. sel. equilibrium', ha='center', fontsize=7.5,
         color='#333', style='italic')
ax2.text(1.6, 4.2, '← Equal frequencies\n   (bal. sel. protecting all)', ha='center',
         fontsize=6.8, color='#1B5E20')

# Arrow
ax2.annotate('', xy=(6.8, 2.5), xytext=(5.0, 2.5),
             arrowprops=dict(arrowstyle='->', color='#B71C1C', lw=2.0,
                             mutation_scale=18))
ax2.text(5.9, 2.8, 'diversity\ncollapse', ha='center', fontsize=7,
         color='#B71C1C', style='italic')

# "Collapsed" bar chart (right)
x0_col = 7.0
for j, (f, col) in enumerate(zip(freqs_col, alleles_col)):
    ax2.bar(x0_col + j * (bw + 0.1), f, width=bw, color=col,
            ec='white', linewidth=0.8, bottom=0)
ax2.text(8.3, -0.5, 'Post-collapse', ha='center', fontsize=7.5,
         color='#333', style='italic')
ax2.text(8.3, 4.2, '← A dominates;\n   rare alleles cannot\n   recover',
         ha='center', fontsize=6.8, color='#B71C1C')

# ─── SCHEMATIC 3: Mate Limitation ────────────────────────────────────────────
ax3 = draw_stage(2)
ax3.set_xlim(0, 10)
ax3.set_ylim(0, 4.5)

# Population of individuals: 7 AAAA (red) + 3 other (mixed)
ind_pos = [
    (1.2, 3.5, A_COL, 'AAAA'), (2.4, 3.7, A_COL, 'AAAA'),
    (3.6, 3.4, A_COL, 'AAAA'), (1.5, 2.5, A_COL, 'AAAA'),
    (3.0, 2.6, A_COL, 'AAAA'), (2.0, 1.6, A_COL, 'AAAA'),
    (3.8, 1.8, A_COL, 'AAAA'),
    (5.4, 3.6, B_COL, 'AABB'), (6.0, 2.5, C_COL, 'AABC'),
    (5.6, 1.4, B_COL, 'AABB'),
]

# Draw compatible pair lines first (between non-AAAA)
compatible_pairs = [(7, 8), (7, 9), (8, 9)]
for pi, pj in compatible_pairs:
    xi, yi = ind_pos[pi][0], ind_pos[pi][1]
    xj, yj = ind_pos[pj][0], ind_pos[pj][1]
    ax3.plot([xi, xj], [yi, yj], color='#43A047', lw=1.8,
             alpha=0.7, zorder=1)

# Draw incompatible "no-connection" X for AAAA pairs (just one example)
ax3.annotate('', xy=(2.9, 3.55), xytext=(1.4, 3.5),
             arrowprops=dict(arrowstyle='-', color='#EF9A9A',
                             lw=1.2, linestyle='dashed'))
ax3.text(2.15, 3.75, '✗', ha='center', fontsize=9, color='#C62828')

# Draw the individuals
for (ix, iy, col, lbl) in ind_pos:
    ax3.add_patch(Circle((ix, iy), 0.30, fc=col, ec='white', lw=1.2, zorder=3))
    ax3.text(ix, iy - 0.62, lbl, ha='center', fontsize=5.8,
             color='#333', zorder=4)

# Legend
ax3.plot([7.8, 8.5], [1.2, 1.2], color='#43A047', lw=1.8)
ax3.text(8.6, 1.2, 'compatible', va='center', fontsize=7, color='#43A047')
ax3.plot([7.8, 8.5], [0.6, 0.6], color='#EF9A9A', lw=1.2, linestyle='dashed')
ax3.text(8.6, 0.6, 'incompatible', va='center', fontsize=7, color='#C62828')

ax3.text(2.5, 0.3, '7 × homozygous (AAAA) — isolated', ha='center', fontsize=7.5,
         color='#B71C1C', fontweight='bold')
ax3.text(5.7, 0.3, '3 × diverse — connected', ha='center', fontsize=7.5,
         color='#1565C0', fontweight='bold')

# ─── SCHEMATIC 4: Genetic Drift ──────────────────────────────────────────────
ax4 = draw_stage(3)
ax4.set_xlim(0, 10)
ax4.set_ylim(0, 4.5)

# Show two generations side by side with allele frequency bars + small N label
gen1_alleles = [A_COL] * 4 + [B_COL] * 2 + [C_COL] * 2 + [D_COL] * 1
gen2_alleles = [A_COL] * 5 + [B_COL] * 2 + [C_COL] * 1   # D lost by drift

def mini_pop(ax_in, x0, allele_list, gen_label, n_label, lost=None):
    cols_unique = list(dict.fromkeys(allele_list))
    counts = {c: allele_list.count(c) for c in cols_unique}
    total = len(allele_list)
    bar_w = 1.1
    freqs = [counts[c] / total for c in cols_unique]
    bottoms = np.cumsum([0] + freqs[:-1])
    for fc, f, bot in zip(cols_unique, freqs, bottoms):
        ax_in.bar(x0, f, width=bar_w, bottom=bot, color=fc,
                  ec='white', lw=0.8)
    ax_in.text(x0, -0.45, gen_label, ha='center', fontsize=7.8, color='#333')
    ax_in.text(x0, -0.85, n_label, ha='center', fontsize=7.2, color='#777')
    if lost:
        ax_in.text(x0, 1.08, f'"{lost}" lost\nby drift', ha='center',
                   fontsize=7, color='#B71C1C', style='italic')

mini_pop(ax4, 2.0, gen1_alleles, 'Generation t', 'N = 9 allele copies',)
mini_pop(ax4, 5.5, gen2_alleles, 'Generation t+1', 'N = 8 allele copies', lost='D')

# Arrow between generations
ax4.annotate('', xy=(4.7, 0.6), xytext=(3.0, 0.6),
             arrowprops=dict(arrowstyle='->', color=ARROW_COL,
                             lw=1.8, mutation_scale=16))
ax4.text(3.85, 0.82, 'stochastic\nsampling', ha='center',
         fontsize=7, color='#555', style='italic')

# Small N emphasis
ax4.text(7.0, 2.2, 'N = 25–40\nper population', ha='center', fontsize=10,
         fontweight='bold', color='#B71C1C',
         bbox=dict(boxstyle='round,pad=0.4', fc='#FFCDD2', ec='#C62828', lw=1.5))
ax4.text(7.0, 0.3, 'drift\noverwhelms\nbal. selection', ha='center', fontsize=8,
         color='#555', style='italic')

ax4.set_ylim(-1.0, 4.5)

# ─── SCHEMATIC 5: Polyploid SI Breakdown — clean 3-zone layout ────────────────
# Use a custom inset with minimal padding so the schematic is taller
draw_stage(4)
_y5_bot = STAGE_Y[4] - BOX_H
ax5 = inset(SCH_X0, _y5_bot + 0.12, SCH_W, BOX_H - 0.57)
ax5.set_xlim(0, 10)
ax5.set_ylim(0.0, 3.5)

CM5 = {'A': A_COL, 'B': B_COL, 'C': C_COL, 'D': D_COL}

def chrom_stack(ax, cx, cy, alleles, label=None, sw=0.54, sh=0.33):
    """Tetraploid genotype as stacked horizontal coloured bars."""
    n  = len(alleles)
    y0 = cy - n * sh / 2
    for k, al in enumerate(alleles):
        ax.add_patch(mpatches.FancyBboxPatch(
            (cx - sw/2, y0 + k*sh), sw, sh*0.86,
            boxstyle='round,pad=0.02', fc=CM5[al], ec='white', lw=1.0))
        ax.text(cx, y0 + k*sh + sh*0.43, al,
                ha='center', va='center', fontsize=9,
                fontweight='bold', color='white')
    if label:
        ax.text(cx, y0 - 0.28, label, ha='center', va='top',
                fontsize=8.5, color='#222', fontweight='bold')
    return y0

def pollen_grain(ax, cx, cy, a1, a2, r=0.21, alpha=1.0):
    """Two side-by-side circles = one diploid pollen grain."""
    gap = r * 1.15
    for k, al in enumerate([a1, a2]):
        ax.add_patch(Circle((cx - gap/2 + k*gap, cy), r,
                             fc=CM5[al], ec='white', lw=0.8, alpha=alpha))
        ax.text(cx - gap/2 + k*gap, cy, al,
                ha='center', va='center', fontsize=8,
                fontweight='bold', color='white', alpha=alpha)

# ── Zone 1: AABB parent ───────────────────────────────────────────────────────
chrom_stack(ax5, 1.4, 2.41, ['A','A','B','B'], 'AABB\nparent', sw=0.52, sh=0.30)

ax5.annotate('', xy=(2.85, 2.41), xytext=(2.05, 2.41),
             arrowprops=dict(arrowstyle='->', color='#555',
                             lw=1.8, mutation_scale=15))
ax5.text(2.45, 2.63, 'makes\npollen', ha='center', fontsize=8,
         color='#555', style='italic')

# ── Zone 2: AB pollen in dashed box ──────────────────────────────────────────
ax5.add_patch(mpatches.FancyBboxPatch(
    (2.90, 1.94), 2.50, 0.90,
    boxstyle='round,pad=0.08',
    fc='#FFEBEE', ec='#B71C1C', lw=2.2, linestyle='dashed', zorder=1))

pollen_grain(ax5, 4.15, 2.38, 'A', 'B', r=0.27)

ax5.text(4.15, 2.73, 'AB pollen  (4 of 6)', ha='center', va='center',
         fontsize=9, fontweight='bold', color='#B71C1C')

ax5.text(4.15, 1.55, 'A & B compete → signal too weak',
         ha='center', va='center', fontsize=7.5, color='#B71C1C',
         fontweight='bold')
ax5.text(4.15, 1.18, 'SI fails → self-fertilises',
         ha='center', va='center', fontsize=7.5, color='#555', style='italic')

ax5.annotate('', xy=(5.85, 2.38), xytext=(5.42, 2.38),
             arrowprops=dict(arrowstyle='->', color='#B71C1C',
                             lw=2.0, mutation_scale=16))

# ── Zone 3: offspring and fixation ───────────────────────────────────────────
chrom_stack(ax5, 7.6, 2.84, ['A','A','A','B'], 'AAAB', sw=0.52, sh=0.30)

ax5.text(8.55, 2.84, 'generation 1', ha='left', va='center',
         fontsize=7.5, color='#555', style='italic')

ax5.annotate('', xy=(7.6, 1.62), xytext=(7.6, 1.95),
             arrowprops=dict(arrowstyle='->', color='#B71C1C',
                             lw=2.0, mutation_scale=15))

ax5.text(8.20, 1.78, 'AA pollen\nnow 3 of 6\n→ drift to AAAA',
         ha='left', va='center', fontsize=7.5, color='#B71C1C',
         linespacing=1.35)

bot5 = chrom_stack(ax5, 7.6, 0.88, ['A','A','A','A'], None, sw=0.52, sh=0.28)

ax5.add_patch(mpatches.FancyBboxPatch(
    (7.14, bot5 - 0.10), 0.98, 4*0.28 + 0.20,
    boxstyle='round,pad=0.06', fc='none', ec='#B71C1C', lw=2.5))

ax5.text(7.6, bot5 - 0.22, 'AAAA (homozygous)  ·  GFS = 0',
         ha='center', va='top', fontsize=8,
         fontweight='bold', color='#B71C1C')

# ─── CONNECTING ARROWS ────────────────────────────────────────────────────────
for i in range(4):
    y_from = STAGE_Y[i] - BOX_H
    y_to   = STAGE_Y[i + 1]
    down_arrow(FW / 2, y_from, y_to)


# ─── SAVE ─────────────────────────────────────────────────────────────────────
os.makedirs('figures', exist_ok=True)
outpath = 'figures/SRK_AAAA_cascade_hypothesis.png'
fig.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
print(f'Saved: {outpath}')
plt.close(fig)
