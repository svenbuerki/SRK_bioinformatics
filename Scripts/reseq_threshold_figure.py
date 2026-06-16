#!/usr/bin/env python3
"""
reseq_threshold_figure.py — focused figure supporting the
~20 000 raw-reads-per-barcode re-sequencing target.

Reads:  Tables/Phase4/step22c_reseq_calibration_per_sample.tsv
Writes: figures/Phase4/step22c_reseq_threshold.{pdf,png}

Restricts to the COVERAGE-DRIVEN cohort (excludes pipeline_downstream_filter,
chimera_dominated, DNA_contamination, partial_recovery, not_in_SI_status_table)
so the curve reflects only samples where more reads is the right intervention.
"""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")          # non-interactive backend — no dock icon
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
IN_TSV  = HERE / "Tables" / "Phase4" / "step22c_reseq_calibration_per_sample.tsv"
FIG_DIR = HERE / "figures" / "Phase4"
OUT_PDF = FIG_DIR / "step22c_reseq_threshold.pdf"
OUT_PNG = FIG_DIR / "step22c_reseq_threshold.png"
FIG_DIR.mkdir(parents=True, exist_ok=True)

EXCLUDE = {"pipeline_downstream_filter", "chimera_dominated", "DNA_contamination",
           "partial_recovery", "not_in_SI_status_table"}
TARGET_READS = 20_000          # the recommended threshold

df = pd.read_csv(IN_TSV, sep="\t", encoding="utf-8-sig")
cov = df[~df["failure_mode"].isin(EXCLUDE)].dropna(subset=["raw_reads"]).copy()
print(f"Coverage-driven cohort: n = {len(cov)} of {len(df)} samples "
      f"(overall pass rate {cov['pass_si_threshold'].mean():.1%})")

# Bin samples for the dose-response bar chart
edges  = [0, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000, 25000, 30000, 1e9]
labels = ["<2.5k", "2.5-5k", "5-7.5k", "7.5-10k", "10-12.5k", "12.5-15k",
          "15-17.5k", "17.5-20k", "20-25k", "25-30k", "30k+"]
cov["bin"] = pd.cut(cov["raw_reads"], bins=edges, labels=labels, include_lowest=True)
summary = cov.groupby("bin", observed=True).agg(
    n=("Sample_ID", "size"),
    pass_rate=("pass_si_threshold", "mean"),
).reset_index()
summary["pct_pass"] = (summary["pass_rate"] * 100).round(1)

# Cumulative thresholds — smallest raw_reads at which P(pass | reads ≥ x) ≥ target
def cumulative_threshold(d: pd.DataFrame, target: float) -> int | None:
    d = d.sort_values("raw_reads")
    for r in sorted(d["raw_reads"].unique()):
        s = d[d["raw_reads"] >= r]
        if len(s) >= 20 and s["pass_si_threshold"].mean() >= target:
            return int(r)
    return None

T80 = cumulative_threshold(cov, 0.80)
T90 = cumulative_threshold(cov, 0.90)
T95 = cumulative_threshold(cov, 0.95)
print(f"Cumulative thresholds — 80 %: {T80:,}  90 %: {T90:,}  95 %: {T95:,}")

# ─── Build the figure ────────────────────────────────────────────────────────
fig, (axL, axR) = plt.subplots(1, 2, figsize=(15, 5.8),
                                gridspec_kw={"width_ratios": [1.0, 1.0]})

# Left — bin-level pass rate
ax = axL
xs = np.arange(len(summary))
# Green when the bin's pass rate clears the 80 % target; blue otherwise.
bar_colors = ["#4DAF4A" if g else "#377EB8" for g in summary["pct_pass"] >= 80]
ax.bar(xs, summary["pct_pass"], color=bar_colors, edgecolor="black", linewidth=0.5)
for x, (n, pct) in enumerate(zip(summary["n"], summary["pct_pass"])):
    ax.text(x, pct + 1.5, f"n={n}", ha="center", fontsize=8)
    ax.text(x, max(pct - 6, 4), f"{pct:.0f}%", ha="center", fontsize=8,
            color="white" if pct >= 25 else "black", fontweight="bold")

ax.axhline(80, color="#FF7F00", linestyle="--", linewidth=1.0, label="80 % target")
ax.axhline(90, color="#E41A1C", linestyle="--", linewidth=1.0, label="90 % target")
ax.axhline(95, color="#984EA3", linestyle="--", linewidth=1.0, label="95 % target")

# Highlight the recommended 20 k bin
target_bin_idx = labels.index("17.5-20k")
ax.axvspan(target_bin_idx - 0.5, target_bin_idx + 0.5,
           color="#4DAF4A", alpha=0.12, zorder=0)
ax.text(target_bin_idx, 102, "↑ recommended\noperational target",
        ha="center", fontsize=9, color="#1B5E20", fontweight="bold")

ax.set_xticks(xs)
ax.set_xticklabels(summary["bin"], rotation=45, ha="right", fontsize=9)
ax.set_ylabel("% of samples reaching an SI call (n_haps_total ≥ 4)", fontsize=10)
ax.set_xlabel("Raw reads per barcode (FASTQ count)", fontsize=10)
ax.set_title("A. Dose-response — pass rate vs raw read budget\n"
             "(coverage-driven cohort only — pipeline-failure / chimera / "
             "DNA-contamination samples excluded)",
             fontsize=10)
ax.set_ylim(0, 112)
ax.legend(loc="lower right", fontsize=8, framealpha=0.95)
ax.grid(axis="y", alpha=0.25, linestyle=":")

# Right — cumulative pass rate vs threshold (the "set the bar here" view)
ax = axR
sorted_reads = np.sort(cov["raw_reads"].values)
cum_pass = []
cum_x    = []
for r in np.unique(sorted_reads):
    s = cov[cov["raw_reads"] >= r]
    if len(s) >= 20:
        cum_pass.append(s["pass_si_threshold"].mean() * 100)
        cum_x.append(r)
cum_x = np.array(cum_x); cum_pass = np.array(cum_pass)
ax.plot(cum_x, cum_pass, color="#377EB8", linewidth=2)
ax.fill_between(cum_x, 0, cum_pass, color="#377EB8", alpha=0.10)

# Mark target thresholds
for thr, color, label in [(T80, "#FF7F00", "80 %"),
                           (T90, "#E41A1C", "90 %"),
                           (T95, "#984EA3", "95 %")]:
    if thr is None: continue
    ax.axhline(int(label.split()[0]), color=color, linestyle="--", linewidth=1)
    ax.axvline(thr, color=color, linestyle=":", linewidth=1)
    ax.annotate(f"{label} reached at\n{thr:,} reads",
                xy=(thr, int(label.split()[0])), xytext=(thr + 1500, int(label.split()[0]) - 12),
                fontsize=8.5, color=color, fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=color, lw=0.8))

# Recommended operational target shaded band
ax.axvspan(TARGET_READS - 2500, TARGET_READS + 2500, color="#4DAF4A", alpha=0.15)
ax.axvline(TARGET_READS, color="#2E7D32", linewidth=1.5)
ax.text(TARGET_READS, 16,
        f"Recommended\noperational target\n= {TARGET_READS:,} raw reads",
        ha="center", fontsize=9, color="#1B5E20", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor="#2E7D32", linewidth=0.8))

ax.set_xlim(0, max(35_000, cum_x.max()))
ax.set_ylim(0, 105)
ax.set_xlabel("Minimum raw reads per barcode (threshold)", fontsize=10)
ax.set_ylabel("Pass rate among samples meeting the threshold (%)", fontsize=10)
ax.set_title("B. Cumulative pass rate — choose a threshold, read off the success rate",
             fontsize=10)
ax.grid(alpha=0.25, linestyle=":")

fig.suptitle("Re-sequencing read-count threshold — empirical support for "
             f"the {TARGET_READS:,}-read-per-barcode operational target",
             fontsize=12, fontweight="bold", y=1.02)
fig.text(0.5, -0.02,
         f"Source: Tables/Phase4/step22c_reseq_calibration_per_sample.tsv  ·  "
         f"n = {len(cov)} coverage-driven samples across 11 libraries  ·  "
         "thresholds = smallest raw_reads at which P(pass | reads ≥ threshold) "
         "first crosses target (n ≥ 20)",
         ha="center", fontsize=8, color="grey")

fig.tight_layout()
fig.savefig(OUT_PDF, bbox_inches="tight")
fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
print(f"Saved {OUT_PDF}")
print(f"Saved {OUT_PNG}")
