#!/usr/bin/env python3
"""
SRK_null_allele_assignment.py — Step 28

Recover the AA sequence of every REMOVED (premature-stop) haplotype dropped by
Step 7 and assign each to its nearest *functional* SRK allele by AA p-distance.
Adds a nucleotide-identity check at the first stop codon (from the aligned
backfilled DNA) so the downstream IBD scorer can distinguish identity-by-descent
(same allele + same stop position + same DNA mutation) from identity-by-state
(same allele + independent stop mutations).

Why this step exists
────────────────────
The genotyping pipeline labels every broken copy as `Allele_NULL` (Step 26).
That label collapses 49 different broken-allele identities into one bin and
discards the IBD signal that distinguishes inherited broken alleles from
independently-arising LoF events. For breeding purposes we want to avoid
crossing two carriers of the SAME broken allele inherited from a recent
common ancestor — they risk producing offspring with two recombinantly-
trapped copies of that broken haplotype.

Pipeline
────────
  1. Parse every `all_Library*_frame1_stopcodon_log.tsv`; capture each
     `Status == REMOVED` record (Sequence_ID, stop_positions, n_stops).
  2. Extract the AA sequence of each REMOVED hap from the corresponding
     `_frame1_AA_raw.fasta`; replace `*` characters with `X`.
  3. Extract the DNA codon at the first stop position from the corresponding
     `_backfilled.fasta` (ungapped, frame 1).
  4. Filter `SRK_protein_allele_representatives.fasta` to the 49 canonical
     alleles present in `SRK_individual_allele_genotypes_with_nulls.tsv`.
  5. Align the 49 reps with MAFFT, then add the REMOVED haps via
     `mafft --add --keeplength` so the reference coordinates are preserved.
  6. For each REMOVED hap, compute p-distance to every functional rep
     (vectorised, ignoring `X` and `-`); assign to argmin distance.
  7. Confidence flag from (AA_distance, n_stops):
       high   : AA_distance ≤ 0.005 AND n_stops ≤ 3
       medium : AA_distance ≤ 0.05  AND n_stops ≤ 10
       low    : otherwise (likely chimeric Canu assembly)
  8. Write `Tables/SRK_null_allele_assignments.tsv`.

Outputs
───────
  Tables/SRK_REMOVED_haplotypes_AA.fasta              (combined query set)
  Tables/SRK_functional_allele_reps_canonical.fasta   (49 reps, post-filter)
  Tables/SRK_functional_allele_reps_canonical.aln     (MAFFT alignment of reps)
  Tables/SRK_null_allele_combined.aln                 (reps + REMOVED, post-add)
  Tables/SRK_null_allele_assignments.tsv              (the deliverable)
"""
from __future__ import annotations

import csv
import glob
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

REPO        = Path(__file__).resolve().parent
TABLES      = REPO / "Tables"
TABLES.mkdir(exist_ok=True)

STOP_LOG_GLOB = "all_Library*_frame1_stopcodon_log.tsv"
AA_RAW_GLOB   = "all_Library*_frame1_AA_raw.fasta"
BACKFILL_GLOB = "all_Library*_aligned_exons_backfilled.fasta"

GENO_TSV   = TABLES / "SRK_individual_allele_genotypes_with_nulls.tsv"
REP_FASTA  = REPO / "SRK_protein_allele_representatives.fasta"

OUT_QUERY_FA = TABLES / "SRK_REMOVED_haplotypes_AA.fasta"
OUT_REPS_FA  = TABLES / "SRK_functional_allele_reps_canonical.fasta"
OUT_REPS_ALN = TABLES / "SRK_functional_allele_reps_canonical.aln"
OUT_COMBO    = TABLES / "SRK_null_allele_combined.aln"
OUT_ASSIGN   = TABLES / "SRK_null_allele_assignments.tsv"

# Confidence thresholds (Sven, 2026-06-11)
THRESH_HIGH_DIST = 0.005   # Step 10 Kneedle elbow (per allele cluster)
THRESH_HIGH_STOP = 3
THRESH_MED_DIST  = 0.05
THRESH_MED_STOP  = 10

OUTGROUP_TAG = "SRK_BEA"


# ─── Helpers ──────────────────────────────────────────────────────────────────
def parse_library(path: str) -> str:
    m = re.search(r"all_(Library\d+)_", path)
    return m.group(1) if m else "Unknown"


def individual_from_seq_id(seq_id: str) -> str:
    """Sequence_ID is e.g. 'Library005_barcode06_hap1_rep3_tig00000003'."""
    parts = seq_id.split("_")
    return "_".join(parts[:2])


def first_stop_position(stop_positions_str: str) -> int | None:
    if not stop_positions_str:
        return None
    try:
        return int(stop_positions_str.split(",")[0])
    except ValueError:
        return None


# ─── Stage 1: collect REMOVED records ────────────────────────────────────────
print("Stage 1 — Parsing stop-codon logs")
removed_records: dict[str, dict] = {}
for log_path in sorted(glob.glob(str(REPO / STOP_LOG_GLOB))):
    library = parse_library(log_path)
    with open(log_path, "r", encoding="utf-8-sig") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            if row["Status"] != "REMOVED":
                continue
            sid = row["Sequence_ID"]
            if sid.startswith(OUTGROUP_TAG):
                continue
            removed_records[sid] = {
                "Sequence_ID": sid,
                "Library": library,
                "Individual": individual_from_seq_id(sid),
                "AA_length": int(row["AA_length"]),
                "n_stops": int(row["Stop_count"]),
                "stop_positions": row["Stop_positions"],
                "first_stop_position": first_stop_position(row["Stop_positions"]),
            }

print(f"  REMOVED haps (excluding outgroup): {len(removed_records)}")


# ─── Stage 2: extract AA sequences ────────────────────────────────────────────
print("Stage 2 — Extracting AA sequences from _frame1_AA_raw.fasta")
aa_seqs: dict[str, str] = {}
target_sids = set(removed_records)
for fa_path in sorted(glob.glob(str(REPO / AA_RAW_GLOB))):
    for rec in SeqIO.parse(fa_path, "fasta"):
        if rec.id in target_sids:
            # Replace stops with X (ambiguous AA — MAFFT-friendly)
            aa_seqs[rec.id] = str(rec.seq).replace("*", "X").upper()

missing_aa = target_sids - set(aa_seqs)
if missing_aa:
    print(f"  WARNING: {len(missing_aa)} REMOVED haps not found in AA_raw fastas")
print(f"  Recovered {len(aa_seqs)} AA sequences")


# ─── Stage 3: extract DNA codon at first stop ────────────────────────────────
print("Stage 3 — Extracting DNA codons at first stop position")
dna_codons: dict[str, str] = {}
for fa_path in sorted(glob.glob(str(REPO / BACKFILL_GLOB))):
    for rec in SeqIO.parse(fa_path, "fasta"):
        if rec.id in target_sids:
            stop_aa = removed_records[rec.id]["first_stop_position"]
            if stop_aa is None or stop_aa < 1:
                continue
            ungapped = str(rec.seq).replace("-", "").upper()
            codon_start = (stop_aa - 1) * 3
            codon = ungapped[codon_start:codon_start + 3]
            dna_codons[rec.id] = codon

print(f"  DNA codons retrieved for {len(dna_codons)} REMOVED haps")


# ─── Stage 4: filter functional reps to canonical 49 ─────────────────────────
print("Stage 4 — Filtering functional allele representatives to canonical set")
geno = pd.read_csv(GENO_TSV, sep="\t", encoding="utf-8-sig", nrows=1)
canonical_alleles = sorted(
    [c for c in geno.columns if c.startswith("Allele_") and c != "Allele_NULL"],
    key=lambda s: int(s.split("_")[1])
)
print(f"  Canonical alleles in genotype matrix: {len(canonical_alleles)}")

with OUT_REPS_FA.open("w") as out:
    n_written = 0
    for rec in SeqIO.parse(REP_FASTA, "fasta"):
        allele_id = rec.id.split()[0]
        if allele_id in canonical_alleles:
            out.write(f">{allele_id}\n{str(rec.seq).upper()}\n")
            n_written += 1
print(f"  Wrote {n_written} reps to {OUT_REPS_FA.name}")
if n_written != len(canonical_alleles):
    print(f"  WARNING: {len(canonical_alleles) - n_written} canonical alleles "
          "missing from representative FASTA")


# ─── Stage 5: align reps then add REMOVED haps via MAFFT --add ───────────────
print("Stage 5a — Aligning canonical reps with MAFFT (--auto)")
with OUT_REPS_ALN.open("w") as out_fh:
    subprocess.run(["mafft", "--auto", "--quiet", str(OUT_REPS_FA)],
                   check=True, stdout=out_fh)
print(f"  Wrote {OUT_REPS_ALN.name}")

print("Stage 5b — Writing combined query FASTA")
with OUT_QUERY_FA.open("w") as out:
    for sid, seq in aa_seqs.items():
        out.write(f">{sid}\n{seq}\n")
print(f"  Wrote {len(aa_seqs)} REMOVED hap sequences to {OUT_QUERY_FA.name}")

print("Stage 5c — Adding REMOVED haps via `mafft --add --keeplength`")
with OUT_COMBO.open("w") as out_fh:
    subprocess.run([
        "mafft", "--add", str(OUT_QUERY_FA), "--keeplength",
        "--quiet", str(OUT_REPS_ALN)
    ], check=True, stdout=out_fh)
print(f"  Wrote {OUT_COMBO.name}")


# ─── Stage 6: load combined alignment and compute p-distances ────────────────
print("Stage 6 — Loading combined alignment + computing p-distances")
ref_seqs: dict[str, np.ndarray] = {}
query_seqs: dict[str, np.ndarray] = {}
for rec in SeqIO.parse(OUT_COMBO, "fasta"):
    arr = np.array(list(str(rec.seq).upper()), dtype="<U1")
    if rec.id in canonical_alleles:
        ref_seqs[rec.id] = arr
    else:
        query_seqs[rec.id] = arr

print(f"  References loaded: {len(ref_seqs)}")
print(f"  Queries loaded:    {len(query_seqs)}")

aln_len = next(iter(ref_seqs.values())).shape[0]
ref_ids = sorted(ref_seqs)
ref_mat = np.vstack([ref_seqs[a] for a in ref_ids])         # (49, L)
ref_valid = (ref_mat != "X") & (ref_mat != "-")             # (49, L)

assignments_rows = []
for qid, qarr in query_seqs.items():
    q_valid = (qarr != "X") & (qarr != "-")                 # (L,)
    pair_valid = ref_valid & q_valid                        # (49, L)
    eq = (ref_mat == qarr)                                  # (49, L)
    n_compared = pair_valid.sum(axis=1)                     # (49,)
    n_diff     = (pair_valid & ~eq).sum(axis=1)             # (49,)
    safe = n_compared > 0
    dists = np.full(len(ref_ids), np.nan, dtype=float)
    dists[safe] = n_diff[safe] / n_compared[safe]
    if not safe.any():
        continue
    order = np.argsort(np.where(safe, dists, np.inf))
    best_idx = int(order[0])
    second_idx = int(order[1]) if len(order) > 1 else best_idx
    rec = removed_records[qid]
    best_d = dists[best_idx]
    second_d = dists[second_idx]
    if best_d <= THRESH_HIGH_DIST and rec["n_stops"] <= THRESH_HIGH_STOP:
        confidence = "high"
    elif best_d <= THRESH_MED_DIST and rec["n_stops"] <= THRESH_MED_STOP:
        confidence = "medium"
    else:
        confidence = "low"
    assignments_rows.append({
        "Sequence_ID":           qid,
        "Library":               rec["Library"],
        "Individual":            rec["Individual"],
        "n_stops":               rec["n_stops"],
        "first_stop_position":   rec["first_stop_position"],
        "stop_codon_DNA":        dna_codons.get(qid, ""),
        "assigned_allele":       ref_ids[best_idx],
        "AA_distance":           round(float(best_d), 5),
        "second_best_allele":    ref_ids[second_idx],
        "second_best_distance":  round(float(second_d), 5),
        "margin":                round(float(second_d - best_d), 5),
        "n_compared_columns":    int(n_compared[best_idx]),
        "confidence":            confidence,
    })

print(f"  Assigned {len(assignments_rows)} REMOVED haps")


# ─── Stage 7: write output table ──────────────────────────────────────────────
df = pd.DataFrame(assignments_rows)
df = df.sort_values(["Individual", "AA_distance"]).reset_index(drop=True)
df.to_csv(OUT_ASSIGN, sep="\t", index=False)
print(f"\nWrote {len(df)} rows → {OUT_ASSIGN}")

# Console summary
print("\nConfidence distribution:")
print(df["confidence"].value_counts().to_string())
print("\nAA distance summary by confidence:")
print(df.groupby("confidence")["AA_distance"]
        .describe()[["count","mean","50%","max"]].to_string())
print("\nAssigned-allele counts (top 10):")
print(df["assigned_allele"].value_counts().head(10).to_string())
print("\nDistinct REMOVED haps per individual (top 10):")
print(df.groupby("Individual").size().sort_values(ascending=False).head(10).to_string())
