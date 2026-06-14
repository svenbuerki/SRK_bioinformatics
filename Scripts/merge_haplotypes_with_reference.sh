#!/bin/bash
set -euo pipefail

# ==========================================
# Merge all oriented haplotypes and canonical references
# ==========================================

# ----------------------------
# User input
# ----------------------------
read -p "Enter Library name (e.g., Library_0001): " Library
read -p "Enter canonical reference FASTA: " REF_FASTA

# Resolve absolute path to reference FASTA
REF_FASTA=$(realpath "$REF_FASTA")

if [[ ! -f "$REF_FASTA" ]]; then
    echo "ERROR: Reference FASTA not found: $REF_FASTA" >&2
    exit 1
fi

FINAL_OUT="all_${Library}_Phased_haplotypes_with_ref.fasta"

# Create / truncate output file at project root
> "$FINAL_OUT"

echo
echo "Merging phased haplotypes from Library: $Library"
echo "Adding canonical references: $REF_FASTA"
echo

# ----------------------------
# Step 1: Append canonical references first
# ----------------------------
cat "$REF_FASTA" >> "$FINAL_OUT"

# ----------------------------
# Step 2: Append all oriented haplotypes from all barcodes
# ----------------------------
find "$Library" -type f -name "*_Phased_haplotypes.oriented.fasta" | sort | while read -r HAP_FILE; do
    echo "→ Adding $HAP_FILE"
    cat "$HAP_FILE" >> "$FINAL_OUT"
done

echo
echo "======================================="
echo "Merged FASTA (haplotypes + references) written to:"
echo "→ $FINAL_OUT"
echo "======================================="

