#!/bin/bash
set -euo pipefail

# ==========================================
# Filter FASTA sequences by minimum and maximum length
# and add canonical references
# ==========================================

# ----------------------------
# User input
# ----------------------------
read -p "Enter merged haplotype FASTA file: " MERGED_FASTA
read -p "Enter canonical reference FASTA file: " REF_FASTA
read -p "Enter minimum sequence length (bp): " MINLEN
read -p "Enter maximum sequence length (bp): " MAXLEN

# Check if input files exist
for f in "$MERGED_FASTA" "$REF_FASTA"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: File not found: $f" >&2
        exit 1
    fi
done

# Determine output file name
BASENAME="${MERGED_FASTA%.fasta}"
OUT="${BASENAME}_filtered_min${MINLEN}_max${MAXLEN}.fasta"

# ----------------------------
# Step 1: Combine canonical references with merged FASTA
# ----------------------------
TMP_COMBINED=$(mktemp)
cat "$REF_FASTA" "$MERGED_FASTA" > "$TMP_COMBINED"

# ----------------------------
# Step 2: Filter sequences by minimum and maximum length using seqkit
# ----------------------------
seqkit seq -m "$MINLEN" -M "$MAXLEN" "$TMP_COMBINED" > "$OUT"

# Cleanup
rm -f "$TMP_COMBINED"

echo "Filtered FASTA written to: $OUT"

