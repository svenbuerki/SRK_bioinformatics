#!/bin/bash
set -euo pipefail

# ==========================================
# Summarize haplotypes by barcode from FASTA
# Header format:
# hapX_Library_XXXX_barcodeYY_tigZZZZZZZZ
# ==========================================

# ----------------------------
# User input
# ----------------------------
read -p "Enter FASTA file: " FASTA

if [[ ! -f "$FASTA" ]]; then
    echo "ERROR: FASTA file not found: $FASTA" >&2
    exit 1
fi

OUT="${FASTA%.fasta}_haplotype_summary.tsv"

# ----------------------------
# Parse FASTA headers
# ----------------------------
awk '
  /^>/ {
    sub(/^>/,"",$0)
    split($0, a, "_")

    hap = a[1]
    barcode = a[4]
    tig = a[5]

    seqs[barcode]++
    haps[barcode, hap] = 1
    tigs[barcode, tig] = 1
  }
  END {
    print "Barcode\tSequences\tHaplotypes\tTigs"
    for (b in seqs) {
      nh=0; nt=0
      for (k in haps)
        if (k ~ "^" b SUBSEP) nh++
      for (k in tigs)
        if (k ~ "^" b SUBSEP) nt++
      print b "\t" seqs[b] "\t" nh "\t" nt
    }
  }
' "$FASTA" | sort -V > "$OUT"

# ----------------------------
# Report
# ----------------------------
echo
echo "======================================="
echo "Haplotype summary written to:"
echo "→ $OUT"
echo "======================================="
echo
echo "Preview:"
column -t "$OUT"

