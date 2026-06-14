#!/bin/bash
set -euo pipefail

# ==========================================
# Merge AA FASTAs, keep one copy of canonical
# sequences (by ID), then realign
# ==========================================

# ----------------------------
# User input
# ----------------------------
read -p "Enter directory to search (e.g., .): " SEARCH_DIR
read -p "Enter canonical FASTA (DNA or AA, IDs only are used): " CANONICAL_FASTA
read -p "Enter output prefix (e.g., SRK_all_AA): " PREFIX

# Resolve paths
SEARCH_DIR=$(realpath "$SEARCH_DIR")
CANONICAL_FASTA=$(realpath "$CANONICAL_FASTA")

if [[ ! -d "$SEARCH_DIR" ]]; then
    echo "ERROR: Directory not found: $SEARCH_DIR" >&2
    exit 1
fi

if [[ ! -f "$CANONICAL_FASTA" ]]; then
    echo "ERROR: Canonical FASTA not found: $CANONICAL_FASTA" >&2
    exit 1
fi

MERGED_RAW="${PREFIX}_merged_raw.fasta"
MERGED_DEDUP="${PREFIX}_merged_dedup.fasta"
ALIGNED_FINAL="${PREFIX}_merged_aligned.fasta"

# ----------------------------
# Step 1: Extract canonical IDs
# ----------------------------
echo "Extracting canonical IDs..."
grep "^>" "$CANONICAL_FASTA" | sed 's/^>//' > canonical_ids.txt

# ----------------------------
# Step 2: Merge AA FASTAs (NON-RECURSIVE)
# ----------------------------
echo "Merging *_AA_filtered_aligned.fasta files (current directory only)..."
> "$MERGED_RAW"

find "$SEARCH_DIR" -maxdepth 1 -type f -name "*_AA_filtered_aligned.fasta" | sort | while read -r FAA; do
    echo "→ Adding $FAA"
    cat "$FAA" >> "$MERGED_RAW"
done

# Sanity check
if [[ ! -s "$MERGED_RAW" ]]; then
    echo "ERROR: No *_AA_filtered_aligned.fasta files found in $SEARCH_DIR" >&2
    exit 1
fi

# ----------------------------
# Step 3: Deduplicate canonical sequences
#   - keep first occurrence of each canonical ID
#   - keep all non-canonical sequences
# ----------------------------
echo "Deduplicating canonical sequences..."

awk '
BEGIN {
    while ((getline < "canonical_ids.txt") > 0) {
        canon[$1] = 1
    }
}
{
    if ($0 ~ /^>/) {
        id = substr($0, 2)
        keep = !(id in canon_seen)
        if (id in canon) {
            canon_seen[id] = 1
        }
    }
    if (keep) {
        print
    }
}
' "$MERGED_RAW" > "$MERGED_DEDUP"

# ----------------------------
# Step 4: Realign merged AA FASTA
# ----------------------------
echo "Realigning merged AA FASTA with MAFFT..."
mafft --auto "$MERGED_DEDUP" > "$ALIGNED_FINAL"

# ----------------------------
# Cleanup
# ----------------------------
rm -f canonical_ids.txt "$MERGED_RAW"

echo
echo "======================================="
echo "Final aligned AA FASTA:"
echo "→ $ALIGNED_FINAL"
echo "======================================="

