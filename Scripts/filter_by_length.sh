#!/bin/bash
set -euo pipefail

# Ask user for inputs
read -p "Enter input FASTA file: " FASTA
read -p "Enter minimum length (bp): " MINLEN
read -p "Enter maximum length (bp): " MAXLEN

OUT="${FASTA%.fasta}_min${MINLEN}_max${MAXLEN}.fasta"

seqkit seq -m "$MINLEN" -M "$MAXLEN" "$FASTA" > "$OUT"

echo "Filtered FASTA written to: $OUT"
