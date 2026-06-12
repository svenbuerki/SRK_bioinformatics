#!/usr/bin/env python3
"""
chimera_coverage_filter.py — post-Canu chimera detection via read coverage

Reads a per-base depth file from `samtools depth -a coverage.bam` and computes
per-contig coverage statistics that detect chimeric Canu contigs:

  mean_cov         mean read depth across the contig
  median_cov       median read depth
  min_window_cov   minimum mean depth across overlapping sliding windows of
                   size WINDOW (default 200 bp) — chimeric contigs have
                   localised low-coverage junctions that depress this value
  uniformity       min_window_cov / median_cov

A contig is flagged DROP if:
    mean_cov < MIN_MEAN_COV   OR   uniformity < MIN_UNIFORMITY

In `report` mode every contig is written to the output FASTA unchanged
(stats TSV only — used for threshold calibration on Library 005).
In `filter` mode only KEEP contigs are written to the output FASTA.

Usage
─────
    chimera_coverage_filter.py \
        --depth depth.tsv \
        --in   combined_unique_contigs.fasta \
        --out  combined_unique_contigs_chimera_filtered.fasta \
        --stats chimera_filter_stats.tsv \
        --mode {report|filter} \
        [--window 200] [--min-mean-cov 20] [--min-uniformity 0.2]

Stats TSV columns:
    contig  length  mean_cov  median_cov  min_window_cov  uniformity  decision
"""
from __future__ import annotations

import argparse
import statistics
from collections import defaultdict
from pathlib import Path


def parse_depth(depth_path: str) -> dict[str, list[int]]:
    """Return {contig: dense_depth_array} from a `samtools depth -a` file.

    samtools depth output is 1-indexed (pos, depth) per row, one row per
    base. With `-a`, positions with zero coverage are included.
    """
    per_contig_pos_depth: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open(depth_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            contig, pos, depth = parts[0], int(parts[1]), int(parts[2])
            per_contig_pos_depth[contig].append((pos, depth))
    dense: dict[str, list[int]] = {}
    for contig, items in per_contig_pos_depth.items():
        items.sort()
        max_pos = items[-1][0]
        arr = [0] * max_pos
        for pos, depth in items:
            arr[pos - 1] = depth
        dense[contig] = arr
    return dense


def sliding_window_min_mean(arr: list[int], window: int) -> float:
    """Minimum mean depth over any window-size window. O(n) running sum."""
    if not arr:
        return 0.0
    if len(arr) <= window:
        return sum(arr) / len(arr)
    window_sum = sum(arr[:window])
    min_mean = window_sum / window
    for i in range(window, len(arr)):
        window_sum += arr[i] - arr[i - window]
        mean = window_sum / window
        if mean < min_mean:
            min_mean = mean
    return min_mean


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--depth", required=True,
                    help="Per-base depth TSV from `samtools depth -a`")
    ap.add_argument("--in", dest="fasta_in", required=True,
                    help="Input Canu FASTA")
    ap.add_argument("--out", dest="fasta_out", required=True,
                    help="Output filtered FASTA")
    ap.add_argument("--stats", required=True,
                    help="Output per-contig stats TSV")
    ap.add_argument("--mode", choices=["report", "filter"], default="report",
                    help="report = keep all contigs (calibration); "
                         "filter = drop DROP-flagged contigs")
    ap.add_argument("--window", type=int, default=200,
                    help="Sliding window size for min_window_cov (bp)")
    ap.add_argument("--min-mean-cov", type=float, default=20.0,
                    help="Drop contigs with mean_cov below this")
    ap.add_argument("--min-uniformity", type=float, default=0.2,
                    help="Drop contigs with min_window_cov / median_cov below this")
    args = ap.parse_args()

    depth_arrays = parse_depth(args.depth)

    stats: dict[str, dict] = {}
    for contig, arr in depth_arrays.items():
        length = len(arr)
        mean_cov = sum(arr) / length if length else 0.0
        median_cov = statistics.median(arr) if length else 0.0
        min_window = sliding_window_min_mean(arr, args.window)
        uniformity = (min_window / median_cov) if median_cov > 0 else 0.0
        drop = (mean_cov < args.min_mean_cov) or (uniformity < args.min_uniformity)
        stats[contig] = {
            "length": length,
            "mean_cov": mean_cov,
            "median_cov": median_cov,
            "min_window_cov": min_window,
            "uniformity": uniformity,
            "decision": "DROP" if drop else "KEEP",
        }

    with open(args.stats, "w") as fh:
        fh.write("contig\tlength\tmean_cov\tmedian_cov\tmin_window_cov\t"
                 "uniformity\tdecision\n")
        for contig, s in sorted(stats.items()):
            fh.write(f"{contig}\t{s['length']}\t"
                     f"{s['mean_cov']:.2f}\t{s['median_cov']:.0f}\t"
                     f"{s['min_window_cov']:.2f}\t{s['uniformity']:.3f}\t"
                     f"{s['decision']}\n")

    if args.mode == "filter":
        keep_set = {c for c, s in stats.items() if s["decision"] == "KEEP"}
    else:
        keep_set = set(stats.keys())   # report: pass-through

    with open(args.fasta_in) as fin, open(args.fasta_out, "w") as fout:
        skip = False
        for line in fin:
            if line.startswith(">"):
                contig = line[1:].split()[0]
                skip = (args.mode == "filter" and contig not in keep_set)
            if not skip:
                fout.write(line)

    n_kept = sum(1 for s in stats.values() if s["decision"] == "KEEP")
    n_drop = sum(1 for s in stats.values() if s["decision"] == "DROP")
    print(f"chimera_coverage_filter [{args.mode}]: "
          f"{n_kept} KEEP / {n_drop} DROP "
          f"(thresholds: mean_cov ≥ {args.min_mean_cov}, "
          f"uniformity ≥ {args.min_uniformity}, window = {args.window} bp)")


if __name__ == "__main__":
    main()
