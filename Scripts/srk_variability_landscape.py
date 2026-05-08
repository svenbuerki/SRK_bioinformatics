#!/usr/bin/env python3
"""
srk_variability_landscape.py   (Step 22 — Part 1, three-species version)

Three-species cross-Brassicaceae HV-region comparison for the SRK S-domain.

Pipeline:
  1. Load the combined mafft-add alignment (LEPA representatives + Brassica +
     Arabidopsis SRK references, all in one alignment).
  2. Identify LEPA-original alignment columns (every column where ≥1 LEPA
     sequence has a non-gap residue) — defines the coordinate system.
  3. Compute per-column Shannon entropy AND Wu-Kabat variability separately
     for each species group (LEPA, Brassica, Arabidopsis), restricted to the
     LEPA-original columns.
  4. Identify HV regions per species (sliding-window smoothing, then
     mean + k×SD threshold on the smoothed Shannon entropy profile).
  5. Overlay the three variability profiles + per-species HV bands +
     Brassica SCR9-contact residues from Ma et al. 2016 (PDB 5GYY).
  6. Permutation test: is LEPA HV ↔ Brassica HV overlap greater than expected
     by chance, given the number and total length of HV regions in each?

INPUTS
  SRK_combined_alignment.fasta            — mafft --add output containing
                                             LEPA + Brassica + Arabidopsis SRK
  SRK_brassica_hv_mapping.tsv             — Brassica-residue → LEPA-aln-col
                                             mapping (from srk_brassica_hv_mapping.py)

OUTPUTS (root)
  SRK_variability_landscape.pdf           — three-layer comparison plot
  SRK_HV_regions_per_species.tsv          — HV region spans per species
  SRK_HV_overlap_permutation.tsv          — LEPA vs Brassica overlap stats

OUTPUTS (figures/)
  SRK_variability_landscape.png
"""

from __future__ import annotations
import os
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
# Settings
# =============================================================================
COMBINED_ALN  = "SRK_combined_alignment.fasta"
BRASSICA_TSV  = "SRK_brassica_hv_mapping.tsv"
DOMAIN_REGION = (31, 430)        # LEPA-original 1-based column range to plot
WINDOW_SIZE   = 20               # sliding-window width for smoothing
PEAK_SD_FACTOR = 1.0             # HV threshold = mean + k×SD
MIN_HV_RUN     = 3               # minimum consecutive HV cols to call a region
PERM_N         = 10000           # number of permutations for overlap test

# Species classification by record-id substring
SPECIES_KEYS = [
    ("LEPA",        "Allele_"),                # LEPA representatives
    ("Brassica",    "Brassica_"),              # B. rapa, B. oleracea
    ("Brassica",    "BrSRK9"),                 # explicit Ma 2016 reference
    ("Arabidopsis", "Arabidopsis_"),           # A. lyrata + A. halleri
]

# Plot palette
SPECIES_COLOURS = {
    "LEPA":        "#1f77b4",   # blue (primary signal)
    "Brassica":    "#d62728",   # red
    "Arabidopsis": "#2ca02c",   # green
}
SPECIES_BAND_ALPHA = {"LEPA": 0.18, "Brassica": 0.16, "Arabidopsis": 0.14}

# Single uniform colour for the Ma 2016 SCR9-contact residues. We deliberately
# do NOT split by hvI/hvII/hvIII because Ma et al. 2016 does not define those
# regions numerically in the text we could verify; the HV classification per
# residue came from an AI summarisation of the paper and could not be confirmed.
# Presenting the 12 contacts as one set is the conservative scientific choice.
MA2016_CONTACT_COLOUR = "#444444"

OUT_PDF = "SRK_variability_landscape.pdf"
OUT_PNG = os.path.join("figures", "SRK_variability_landscape.png")
OUT_HV_TSV  = "SRK_HV_regions_per_species.tsv"
OUT_PERM_TSV = "SRK_HV_overlap_permutation.tsv"
OUT_LEPA_HV_COLS = "SRK_LEPA_HV_positions.tsv"   # canonical input to Step 22 Parts 2-5

os.makedirs("figures", exist_ok=True)

# =============================================================================
# 1. Load combined alignment and tag rows by species
# =============================================================================
print("Loading combined alignment ...")
records = list(SeqIO.parse(COMBINED_ALN, "fasta"))
aln_len = len(records[0].seq)
print(f"  {len(records)} sequences × {aln_len} alignment columns")


def classify(rec_id: str) -> str | None:
    for species, key in SPECIES_KEYS:
        if key in rec_id:
            return species
    return None


species_seqs: dict[str, list[str]] = {"LEPA": [], "Brassica": [], "Arabidopsis": []}
for r in records:
    sp = classify(r.id)
    if sp:
        species_seqs[sp].append(str(r.seq).upper())
    else:
        print(f"  WARNING: unclassified record {r.id}")

for sp, seqs in species_seqs.items():
    print(f"  {sp:12s}: {len(seqs)} sequences")

# =============================================================================
# 2. Identify LEPA-original alignment columns
# =============================================================================
# A column is "LEPA-original" iff at least one LEPA sequence has a non-gap
# residue there. All other columns are reference-only insertions.
lepa_seqs = species_seqs["LEPA"]
lepa_orig_cols = [c for c in range(aln_len)
                  if any(seq[c] != "-" for seq in lepa_seqs)]
print(f"\n  LEPA-original columns: {len(lepa_orig_cols)} "
      f"(out of {aln_len} combined columns)")

# Build a 1-based LEPA-coordinate axis (column index in the original LEPA
# representatives alignment, which had 843 columns).
exp_to_lepa1 = {}
lepa_idx = 0
for c in range(aln_len):
    if c in set(lepa_orig_cols):
        lepa_idx += 1
        exp_to_lepa1[c] = lepa_idx
n_lepa_orig = lepa_idx

# Restrict to S-domain LEPA columns
lo1, hi1 = DOMAIN_REGION
sd_cols_exp = [c for c in lepa_orig_cols
               if lo1 <= exp_to_lepa1[c] <= hi1]
print(f"  S-domain ({lo1}-{hi1}) covers {len(sd_cols_exp)} LEPA columns")

# Per-LEPA-position (1-based) axis for plotting
plot_x = np.array([exp_to_lepa1[c] for c in sd_cols_exp])

# =============================================================================
# 3. Per-species variability metrics
# =============================================================================
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def shannon_entropy(column: list[str]) -> float:
    """Per-column Shannon entropy (bits) over standard 20 AA."""
    counts = {}
    n = 0
    for r in column:
        if r == "-" or r not in AA_ALPHABET:
            continue
        counts[r] = counts.get(r, 0) + 1
        n += 1
    if n == 0:
        return 0.0
    H = 0.0
    for k in counts.values():
        p = k / n
        H -= p * math.log2(p)
    return H


def wu_kabat(column: list[str]) -> float:
    """
    Wu-Kabat variability: V = N · k / n_max
      N = number of sequences (with valid AA at this column)
      k = number of distinct AAs observed
      n_max = count of the most common AA
    Yields V = 1 for fully invariant; V increases with variability.
    """
    counts = {}
    N = 0
    for r in column:
        if r == "-" or r not in AA_ALPHABET:
            continue
        counts[r] = counts.get(r, 0) + 1
        N += 1
    if N == 0 or not counts:
        return 0.0
    k = len(counts)
    n_max = max(counts.values())
    return N * k / n_max


def col_metric(seqs: list[str], cols: list[int], func) -> np.ndarray:
    out = np.zeros(len(cols))
    for i, c in enumerate(cols):
        column = [s[c] for s in seqs]
        out[i] = func(column)
    return out


print("\nComputing per-species variability profiles ...")
profiles = {}
for sp in ("LEPA", "Brassica", "Arabidopsis"):
    seqs = species_seqs[sp]
    if len(seqs) < 2:
        print(f"  {sp:12s}: only {len(seqs)} seq — skipping")
        continue
    H  = col_metric(seqs, sd_cols_exp, shannon_entropy)
    WK = col_metric(seqs, sd_cols_exp, wu_kabat)
    H_smooth  = uniform_filter1d(H,  size=WINDOW_SIZE)
    WK_smooth = uniform_filter1d(WK, size=WINDOW_SIZE)
    profiles[sp] = {"n": len(seqs), "H": H, "WK": WK,
                    "H_smooth": H_smooth, "WK_smooth": WK_smooth}
    print(f"  {sp:12s}: n={len(seqs):2d}  median H = {np.median(H):.3f} bits  "
          f"max H = {H.max():.3f}  median WK = {np.median(WK):.2f}")

# =============================================================================
# 4. Per-species HV-region detection
# =============================================================================
print(f"\nIdentifying HV regions (Shannon entropy > mean + "
      f"{PEAK_SD_FACTOR}×SD; min run = {MIN_HV_RUN}) ...")


def find_hv_runs(profile: np.ndarray, x_lepa: np.ndarray,
                 thr: float, min_run: int) -> list[tuple[int, int, int]]:
    """
    Identify HV regions. Returns list of (start_lepa1, end_lepa1, length)
    where both bounds are 1-based LEPA alignment coordinates and length
    counts S-domain positions in the run. Only consecutive S-domain positions
    are considered (i.e., consecutive in plot_x order).
    """
    mask = profile > thr
    runs = []
    i = 0
    n = len(mask)
    while i < n:
        if mask[i]:
            j = i
            while j + 1 < n and mask[j + 1]:
                j += 1
            run_len = j - i + 1
            if run_len >= min_run:
                runs.append((int(x_lepa[i]), int(x_lepa[j]), run_len))
            i = j + 1
        else:
            i += 1
    return runs


hv_regions = {}
for sp, prof in profiles.items():
    H = prof["H_smooth"]
    thr = H.mean() + PEAK_SD_FACTOR * H.std()
    runs = find_hv_runs(H, plot_x, thr, MIN_HV_RUN)
    hv_regions[sp] = {"threshold": thr, "runs": runs}
    print(f"  {sp:12s}: threshold = {thr:.3f} bits  →  {len(runs)} HV runs  "
          f"({sum(r[2] for r in runs)} HV cols)")
    for r in runs:
        print(f"      LEPA cols {r[0]}-{r[1]}  ({r[2]} aa)")

# Wu-Kabat sanity check: compute parallel HV regions and report concordance
print("\nSanity check — Wu-Kabat HV calls vs Shannon entropy HV calls:")
for sp, prof in profiles.items():
    WK = prof["WK_smooth"]
    thr_wk = WK.mean() + PEAK_SD_FACTOR * WK.std()
    runs_wk = find_hv_runs(WK, plot_x, thr_wk, MIN_HV_RUN)
    cols_h  = set()
    for r in hv_regions[sp]["runs"]:
        cols_h.update(range(r[0], r[1] + 1))
    cols_wk = set()
    for r in runs_wk:
        cols_wk.update(range(r[0], r[1] + 1))
    if cols_h or cols_wk:
        jaccard = len(cols_h & cols_wk) / max(len(cols_h | cols_wk), 1)
    else:
        jaccard = float("nan")
    print(f"  {sp:12s}: Shannon HV cols = {len(cols_h):3d}, "
          f"Wu-Kabat HV cols = {len(cols_wk):3d}, "
          f"Jaccard = {jaccard:.2f}")

# Persist HV region table (per-species, run-level summary)
with open(OUT_HV_TSV, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["Species", "Run_index", "Start_lepa1", "End_lepa1",
                "Length_aa", "Threshold_H_bits"])
    for sp, info in hv_regions.items():
        for k, r in enumerate(info["runs"]):
            w.writerow([sp, k + 1, r[0], r[1], r[2], f"{info['threshold']:.4f}"])
print(f"  Written {OUT_HV_TSV}")

# Persist LEPA HV columns (one row per HV column) — canonical input to
# Step 22 Parts 2-5 (HV distance matrix, UPGMA classes, cross design,
# synonymy network).
lepa_runs = hv_regions.get("LEPA", {}).get("runs", [])
lepa_hv_cols_sorted = []
for r in lepa_runs:
    lepa_hv_cols_sorted.extend(range(r[0], r[1] + 1))
lepa_hv_cols_sorted.sort()
lepa_threshold = hv_regions.get("LEPA", {}).get("threshold", float("nan"))

with open(OUT_LEPA_HV_COLS, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["LEPA_aln_col_1based"])
    for c in lepa_hv_cols_sorted:
        w.writerow([c])
print(f"  Written {OUT_LEPA_HV_COLS} ({len(lepa_hv_cols_sorted)} HV columns; "
      f"threshold = {lepa_threshold:.4f} bits)")

# =============================================================================
# 5. Brassica SCR9-contact residues (Ma et al. 2016) — for figure annotation
# =============================================================================
br_contacts = []
br_hv_spans: dict[str, tuple[int, int]] = {}
if os.path.exists(BRASSICA_TSV):
    print(f"\nLoading Ma 2016 SCR9-contact residues from {BRASSICA_TSV} ...")
    df = pd.read_csv(BRASSICA_TSV, sep="\t")
    for _, row in df.iterrows():
        col = row["LEPA_aln_col"]
        if str(col).startswith("NA"):
            continue
        col_int = int(col)
        if lo1 <= col_int <= hi1:
            br_contacts.append({
                "residue": int(row["Brassica_residue"]),
                "aa": row["Brassica_aa"],
                "hv": row["HV_region"],
                "col": col_int,
            })
    for c in br_contacts:
        prev = br_hv_spans.get(c["hv"])
        if prev is None:
            br_hv_spans[c["hv"]] = (c["col"], c["col"])
        else:
            br_hv_spans[c["hv"]] = (min(prev[0], c["col"]), max(prev[1], c["col"]))
    print(f"  {len(br_contacts)} Ma 2016 contact residues mapped to LEPA coords")
else:
    print(f"\n  ({BRASSICA_TSV} not found — skipping Ma 2016 markers)")

# =============================================================================
# 6. Permutation tests: HV-region overlap significance, all pairwise
# =============================================================================
print(f"\nPermutation tests ({PERM_N} permutations) for HV-region overlap ...")

species_hv_cols: dict[str, set[int]] = {}
for sp in ("LEPA", "Brassica", "Arabidopsis"):
    cols = set()
    for r in hv_regions.get(sp, {}).get("runs", []):
        cols.update(range(r[0], r[1] + 1))
    species_hv_cols[sp] = cols

domain_cols = list(plot_x)
domain_n = len(domain_cols)


def perm_overlap(set_a: set[int], set_b: set[int],
                 universe: list[int], n: int, seed: int = 42):
    if not set_a or not set_b:
        return None, None, None, None
    rng = np.random.default_rng(seed)
    nA = len(set_a)
    nB = len(set_b)
    null = np.zeros(n, dtype=int)
    for k in range(n):
        a = rng.choice(universe, size=nA, replace=False)
        b = rng.choice(universe, size=nB, replace=False)
        null[k] = len(set(a.tolist()) & set(b.tolist()))
    obs = len(set_a & set_b)
    p_one_sided = float((null >= obs).sum() + 1) / (n + 1)
    return obs, float(null.mean()), float(null.std()), p_one_sided


pairs = [("LEPA", "Brassica"), ("LEPA", "Arabidopsis"),
         ("Brassica", "Arabidopsis")]
perm_rows = []
for a, b in pairs:
    obs, mu, sd, p_val = perm_overlap(species_hv_cols[a], species_hv_cols[b],
                                       domain_cols, PERM_N)
    if obs is None:
        print(f"  {a} vs {b}: skipped (no HV runs)")
        continue
    print(f"  {a} ↔ {b}: obs={obs} cols  null={mu:.1f}±{sd:.1f}  "
          f"p = {p_val:.4g}")
    perm_rows.append({"pair": f"{a}__{b}",
                      "obs": obs, "null_mean": mu, "null_sd": sd,
                      "p_one_sided": p_val,
                      "n_a": len(species_hv_cols[a]),
                      "n_b": len(species_hv_cols[b])})

if perm_rows:
    with open(OUT_PERM_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=list(perm_rows[0].keys()))
        w.writeheader()
        for row in perm_rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v)
                        for k, v in row.items()})
    print(f"  Written {OUT_PERM_TSV}")

# Pull out the LEPA-Brassica result for the plot annotation
lepa_brassica = next((r for r in perm_rows if r["pair"] == "LEPA__Brassica"),
                    None)
if lepa_brassica:
    observed   = lepa_brassica["obs"]
    null_mean  = lepa_brassica["null_mean"]
    null_sd    = lepa_brassica["null_sd"]
    p_one_sided = lepa_brassica["p_one_sided"]
else:
    observed = null_mean = null_sd = p_one_sided = None

# =============================================================================
# 7. Multi-panel figure: variability profiles + HV-region tracks + Ma 2016
# =============================================================================
print(f"\nWriting figure to {OUT_PDF} ...")

# Layout: top main panel for variability lines, then 3 narrow species
# track strips, then Ma 2016 contact-residue strip at the bottom.
fig = plt.figure(figsize=(17, 7.0))
gs = fig.add_gridspec(
    nrows=5, ncols=1,
    height_ratios=[8, 0.6, 0.6, 0.6, 1.4],
    hspace=0.18,
)
ax_main = fig.add_subplot(gs[0])
ax_lepa = fig.add_subplot(gs[1], sharex=ax_main)
ax_bra  = fig.add_subplot(gs[2], sharex=ax_main)
ax_ara  = fig.add_subplot(gs[3], sharex=ax_main)
ax_ma   = fig.add_subplot(gs[4], sharex=ax_main)

species_axes = {"LEPA": ax_lepa, "Brassica": ax_bra, "Arabidopsis": ax_ara}

# ── Main panel: per-species smoothed Shannon entropy + threshold lines ──────
for sp in ("Arabidopsis", "Brassica", "LEPA"):
    if sp not in profiles:
        continue
    prof = profiles[sp]
    col = SPECIES_COLOURS[sp]
    ax_main.fill_between(plot_x, prof["H"], alpha=0.10, color=col, zorder=1)
    ax_main.plot(plot_x, prof["H_smooth"], color=col, linewidth=1.7,
                 label=f"{sp} (n={prof['n']}, threshold = "
                       f"{hv_regions[sp]['threshold']:.2f} bits)",
                 zorder=2)
    ax_main.axhline(hv_regions[sp]["threshold"],
                    color=col, linestyle=":", linewidth=0.8, alpha=0.6,
                    zorder=1)

ax_main.set_xlim(plot_x[0], plot_x[-1])
ax_main.set_ylabel("Smoothed Shannon entropy (bits)", fontsize=11)
ax_main.set_title(
    "S-domain variability landscape across Brassicaceae — "
    "LEPA vs Brassica vs Arabidopsis SRK alleles, "
    f"sliding-window Shannon entropy (window = {WINDOW_SIZE} aa)",
    fontsize=11.5)
ax_main.legend(fontsize=9, loc="upper right",
               framealpha=0.95, edgecolor="grey")
if p_one_sided is not None:
    ax_main.text(0.005, 0.97,
                 f"LEPA ↔ Brassica HV overlap = {observed} cols\n"
                 f"null mean = {null_mean:.1f} ± {null_sd:.1f}\n"
                 f"p = {p_one_sided:.3g}  ({PERM_N} permutations)",
                 transform=ax_main.transAxes, ha="left", va="top", fontsize=9,
                 bbox=dict(boxstyle="round,pad=0.35",
                           fc="white", ec="grey", lw=0.6))
ax_main.tick_params(axis="x", labelbottom=False)

# ── Per-species HV tracks ───────────────────────────────────────────────────
for sp in ("LEPA", "Brassica", "Arabidopsis"):
    ax_t = species_axes[sp]
    col = SPECIES_COLOURS[sp]
    n_runs = 0
    n_cols = 0
    for r in hv_regions.get(sp, {}).get("runs", []):
        ax_t.axvspan(r[0], r[1], color=col, alpha=0.85, zorder=2)
        n_runs += 1
        n_cols += r[2]
    ax_t.set_yticks([])
    ax_t.set_ylabel(f"{sp}\n{n_runs} HV / {n_cols} aa",
                    rotation=0, ha="right", va="center", fontsize=9,
                    labelpad=6)
    ax_t.tick_params(axis="x", labelbottom=False)
    for spine in ("top", "right", "left"):
        ax_t.spines[spine].set_visible(False)
    ax_t.set_facecolor("#fafafa")

# ── Ma 2016 contact residues (triangles + labels) ──────────────────────────
ax_ma.set_yticks([])
ax_ma.set_ylim(0, 1)
ax_ma.set_facecolor("#fafafa")
ax_ma.set_ylabel("Ma 2016\nSCR9-contact",
                 rotation=0, ha="right", va="center", fontsize=9, labelpad=6)
for spine in ("top", "right", "left"):
    ax_ma.spines[spine].set_visible(False)

if br_contacts:
    for c in br_contacts:
        ax_ma.plot(c["col"], 0.78, marker="^",
                   color=MA2016_CONTACT_COLOUR, markersize=8,
                   markeredgecolor="black", markeredgewidth=0.4,
                   clip_on=False, zorder=5)
        ax_ma.text(c["col"], 0.62, f"{c['aa']}{c['residue']}",
                   ha="center", va="top", fontsize=7,
                   color=MA2016_CONTACT_COLOUR, rotation=90, zorder=5)
    # Add contact-residue marker proxy to the main legend
    ax_main.plot([], [], marker="^", linestyle="None",
                 color=MA2016_CONTACT_COLOUR, markersize=8,
                 markeredgecolor="black", markeredgewidth=0.4,
                 label="Brassica eSRK9 SCR9-contact residues (Ma 2016)")
    # rebuild legend to include the proxy
    ax_main.legend(fontsize=9, loc="upper right",
                   framealpha=0.95, edgecolor="grey")

ax_ma.set_xlabel("LEPA alignment position (aa, 1-based; S-domain)",
                 fontsize=11)

plt.savefig(OUT_PDF, format="pdf", bbox_inches="tight")
plt.savefig(OUT_PNG, format="png", dpi=200, bbox_inches="tight")
plt.close()
print(f"  Saved {OUT_PDF}")
print(f"  Saved {OUT_PNG}")

print("\nDone.")
