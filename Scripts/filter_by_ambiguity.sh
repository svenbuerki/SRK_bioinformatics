#!/bin/bash
set -euo pipefail

# ==========================================
# Step 7b — Ambiguity (N-content) filter
# ==========================================
# Run AFTER Step 7 (backfill) and BEFORE Step 8 (AA translation).
#
# Step 7 backfill pads contigs that are shorter than the canonical
# reference with N's at the 3' end. Those N's translate to X residues
# in Step 8, which MAFFT then accommodates by inflating the AA alignment
# width and inserting gaps in the canonical references.
#
# Library 010 (2026-05-11): 128 sequences had >=10 N's at the 3' end
# (concentrated in a region equivalent to AA cols 812-835). Dropping
# them removes the trailing-gap artefact in BEA canonical refs.
# ==========================================

read -p "Enter FASTA file to filter (Step 7 backfilled output): " INPUT
read -p "Maximum N count per sequence [default 9]: " MAX_N
read -p "Maximum N fraction [0-1, default 0.005]: " MAX_FRAC

MAX_N=${MAX_N:-9}
MAX_FRAC=${MAX_FRAC:-0.005}

[[ -f "$INPUT" ]] || { echo "ERROR: file not found: $INPUT" >&2; exit 1; }

OUT="${INPUT%.fasta}_Nfilt.fasta"
LOG="${INPUT%.fasta}_Nfilt_log.tsv"

awk -v MAX_N="$MAX_N" -v MAX_FRAC="$MAX_FRAC" -v LOG="$LOG" -v OUT="$OUT" '
BEGIN {
  printf "query_id\tlength\tN_count\tN_fraction\tdecision\n" > LOG;
}
function flush(    n, L, frac, keep, decision) {
  if (hdr == "") return;
  n = gsub(/[Nn]/, "&", seq);
  L = length(seq);
  frac = (L > 0) ? n / L : 0;
  keep = (n <= MAX_N && frac <= MAX_FRAC) ? 1 : 0;
  decision = keep ? "KEEP" : "DROP";
  printf "%s\t%d\t%d\t%.4f\t%s\n", hdr, L, n, frac, decision >> LOG;
  if (keep) {
    print ">" full_hdr > OUT;
    print seq          > OUT;
  }
}
/^>/ {
  flush();
  full_hdr = substr($0, 2);
  hdr = full_hdr; sub(/[[:space:]].*$/, "", hdr);
  seq = "";
  next;
}
{ seq = seq $0 }
END { flush(); }
' "$INPUT"

n_in=$(grep -c '^>' "$INPUT")
n_out=$(grep -c '^>' "$OUT" || echo 0)
n_dropped=$((n_in - n_out))

cat <<SUMMARY

Ambiguity filter complete.
  Input sequences:    $n_in
  Kept:               $n_out
  Dropped (N-rich):   $n_dropped
  Thresholds:         max_N=$MAX_N  max_N_fraction=$MAX_FRAC

Outputs:
  Filtered FASTA:     $OUT
  Per-query log:      $LOG
SUMMARY
