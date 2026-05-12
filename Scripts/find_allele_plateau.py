#!/usr/bin/env python3
"""
Identify the plateau in the SRK allele-sensitivity curve quantitatively.

Re-uses the alignment loading, S-domain distance calculation, and UPGMA
linkage from `define_SRK_alleles_from_distance.py` and then applies two
plateau-detection methods to the (threshold, n_alleles) curve:

  1. Slope-magnitude minimum (sliding window) — locates the threshold
     interval where dN/dt is smallest in absolute terms. The plateau N
     is reported as the median allele count over that interval.

  2. Kneedle (Satopaa et al. 2011) — locates the point of maximum
     perpendicular distance from the chord connecting the endpoints of
     the curve. This is the canonical "knee/elbow" definition.

The N_ALLELES at the plateau is the recommended setting for
`define_SRK_alleles_from_distance.py` (Option A).
"""

import sys
import numpy as np
from Bio import SeqIO
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster

# ─── inputs / config ──────────────────────────────────────────────────────────
INPUT_FASTA = "SRK_functional_proteins_aligned.fasta"
DOMAIN_REGION = (31, 430)
THRESH_MIN, THRESH_MAX, THRESH_STEP = 0.0, 0.05, 0.0005

# Constrain plateau search to a biologically plausible allele-count range so
# that the trivial plateaus at the curve's tails (N near n_proteins, or N→1)
# don't dominate. Adjust these if the LEPA dataset shifts substantially.
N_MIN_REASONABLE = 30
N_MAX_REASONABLE = 150

# Window size for slope-magnitude minimum (in number of threshold steps).
# Wider window = more robust to noise; narrower = sharper localisation.
SLOPE_WINDOW = 8   # 8 × 0.0005 = 0.004 threshold-units

# ─── load alignment ───────────────────────────────────────────────────────────
records = list(SeqIO.parse(INPUT_FASTA, "fasta"))
seqs = [str(r.seq).upper() for r in records]
n = len(seqs)
aln_len = len(seqs[0])

start, end = DOMAIN_REGION
seqs_for_dist = [s[start - 1 : end] for s in seqs]

# ─── pairwise p-distance (gaps excluded) ─────────────────────────────────────
def p_distance(s1, s2):
    diffs = comp = 0
    for a, b in zip(s1, s2):
        if a == "-" or b == "-":
            continue
        comp += 1
        if a != b:
            diffs += 1
    return diffs / comp if comp else 0.0

dist = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        d = p_distance(seqs_for_dist[i], seqs_for_dist[j])
        dist[i, j] = dist[j, i] = d

# ─── UPGMA linkage + sensitivity scan ────────────────────────────────────────
Z = linkage(squareform(dist), method="average")

thresholds = np.arange(THRESH_MIN, THRESH_MAX + THRESH_STEP, THRESH_STEP)
allele_counts = np.array(
    [len(set(fcluster(Z, max(t, 1e-12), criterion="distance"))) for t in thresholds]
)

# ─── Method 1 — slope-magnitude minimum within the reasonable N range ────────
# Use a centred sliding-window absolute mean slope.
dN_dt = np.gradient(allele_counts, thresholds)
abs_slope = np.abs(dN_dt)

# rolling mean of |slope|
def rolling_mean(arr, window):
    pad = window // 2
    out = np.full_like(arr, np.nan, dtype=float)
    for i in range(pad, len(arr) - pad):
        out[i] = arr[i - pad : i + pad + 1].mean()
    return out

smooth_slope = rolling_mean(abs_slope, SLOPE_WINDOW)

# Mask: keep only points where N is in the reasonable range
mask = (allele_counts >= N_MIN_REASONABLE) & (allele_counts <= N_MAX_REASONABLE)
search_slope = np.where(mask, smooth_slope, np.inf)

# argmin of smoothed slope within the mask
i_plateau = int(np.nanargmin(search_slope))
t_plateau = thresholds[i_plateau]
N_plateau = allele_counts[i_plateau]

# Report the plateau interval — contiguous neighbours within 1.2× the min slope
near_min = smooth_slope <= smooth_slope[i_plateau] * 1.2
# expand left/right contiguously from i_plateau
L = i_plateau
while L > 0 and near_min[L - 1] and mask[L - 1]:
    L -= 1
R = i_plateau
while R < len(near_min) - 1 and near_min[R + 1] and mask[R + 1]:
    R += 1

t_lo, t_hi = thresholds[L], thresholds[R]
N_lo, N_hi = allele_counts[R], allele_counts[L]  # decreasing in t
N_median = int(np.median(allele_counts[L : R + 1]))

# ─── Method 2 — Kneedle (max perpendicular distance from chord) ──────────────
# Normalise to [0,1] over the reasonable-N region.
if mask.any():
    x = thresholds[mask]
    y = allele_counts[mask].astype(float)
    xn = (x - x.min()) / (x.max() - x.min())
    yn = (y - y.min()) / (y.max() - y.min())
    # Convex-decreasing curve → invert y so the knee is a max-distance point
    # from the chord (0,1) → (1,0).
    dist_to_chord = np.abs(xn + yn - 1) / np.sqrt(2)
    i_knee = int(np.argmax(dist_to_chord))
    t_knee = x[i_knee]
    N_knee = int(y[i_knee])
else:
    t_knee = N_knee = None

# ─── Method 3 — Longest persistent plateau near a target N ───────────────────
# For each candidate centre N_c, scan the threshold axis and find the longest
# contiguous run where allele_counts ∈ [N_c - TOL, N_c + TOL]. The candidate N
# with the widest run is the "most persistent" plateau.
TOL = 3        # tolerance in alleles around plateau centre
N_RANGE = range(N_MIN_REASONABLE, N_MAX_REASONABLE + 1)

best = {"N": None, "run_steps": 0, "t_lo": None, "t_hi": None, "median_N": None}
for N_c in N_RANGE:
    in_band = (allele_counts >= N_c - TOL) & (allele_counts <= N_c + TOL)
    # longest run of True in in_band
    run = 0; best_run = 0; start_idx = end_idx = None; cur_start = None
    for i, v in enumerate(in_band):
        if v:
            if run == 0:
                cur_start = i
            run += 1
            if run > best_run:
                best_run = run
                start_idx = cur_start
                end_idx = i
        else:
            run = 0
    if best_run > best["run_steps"]:
        best.update(
            N=N_c, run_steps=best_run,
            t_lo=thresholds[start_idx], t_hi=thresholds[end_idx],
            median_N=int(np.median(allele_counts[start_idx : end_idx + 1])),
        )

# ─── print result ─────────────────────────────────────────────────────────────
print(f"Sensitivity scan : {len(thresholds)} thresholds in "
      f"[{THRESH_MIN}, {THRESH_MAX}] step {THRESH_STEP}")
print(f"Proteins         : {n}")
print(f"Range searched   : {N_MIN_REASONABLE} ≤ N ≤ {N_MAX_REASONABLE}\n")

print("METHOD 1 — Slope-magnitude minimum (smoothed |dN/dt|)")
print(f"  Plateau centre : threshold ≈ {t_plateau:.4f}  →  N = {N_plateau} alleles")
print(f"  Plateau range  : threshold [{t_lo:.4f}, {t_hi:.4f}]"
      f"  →  N ∈ [{N_lo}, {N_hi}]  (median {N_median})")
print(f"  |dN/dt| at centre: {smooth_slope[i_plateau]:.1f} alleles per unit threshold\n")

print("METHOD 2 — Kneedle (max distance from chord)")
print(f"  Elbow          : threshold ≈ {t_knee:.4f}  →  N = {N_knee} alleles\n")

print(f"METHOD 3 — Longest persistent plateau (N ± {TOL})")
print(f"  Plateau centre : N = {best['N']}  (median over run = {best['median_N']})")
print(f"  Persistence    : {best['run_steps']} consecutive threshold steps "
      f"({best['run_steps'] * THRESH_STEP:.4f} threshold-units wide)")
print(f"  Threshold range: [{best['t_lo']:.4f}, {best['t_hi']:.4f}]\n")

# ─── per-threshold table around the Kneedle elbow ─────────────────────────────
print("Per-threshold curve near the elbow:")
print(f"  {'threshold':>10}  {'N_alleles':>10}  {'|dN/dt|':>10}")
for i, t in enumerate(thresholds):
    if 0.003 <= t <= 0.012:
        print(f"  {t:>10.4f}  {allele_counts[i]:>10d}  {abs_slope[i]:>10.0f}")
print()

print("Recommendation:")
candidates = sorted({N_median, N_knee, best["median_N"]})
print(f"  Three methods returned: slope-min={N_median}, Kneedle={N_knee}, "
      f"persistence={best['median_N']}")
if max(candidates) - min(candidates) <= 5:
    print(f"  Methods agree (Δ ≤ 5). Use N_ALLELES = {int(np.median(candidates))}.")
else:
    print(f"  Methods diverge. The most defensible choice for 'a plateau around N=60' is:")
    print(f"    N_ALLELES = {best['N']}  (longest persistent run near that target),")
    print(f"  or N_ALLELES = {N_knee}  (canonical Kneedle elbow).")
