#!/bin/bash
set -euo pipefail

# ==========================================
# Orient haplotypes for ALL barcodes
# ==========================================


read -p "Enter Library name: " Library
read -p "Enter reference FASTA: " REF_FASTA


REF_FASTA=$(realpath "$REF_FASTA")

FINAL_OUT="all_${Library}_Phased_haplotypes.fasta"

> "$FINAL_OUT"


echo
echo "Processing $Library"
echo


for BARCODE_DIR in "$Library"/barcode*/; do


BARCODE=$(basename "$BARCODE_DIR")

FULL_FASTA="${BARCODE_DIR}/${Library}_${BARCODE}_Phased_haplotypes.fasta"

FILE=$(basename "$FULL_FASTA")


if [[ ! -f "$FULL_FASTA" ]]; then

echo "WARNING missing $FULL_FASTA"

continue

fi


echo "Processing $FULL_FASTA"


cd "$BARCODE_DIR"


ORIENTED="${Library}_${BARCODE}_Phased_haplotypes.oriented.fasta"



# Align

minimap2 -a "$REF_FASTA" "$FILE" > hap_vs_ref.sam



# Get sequences with a primary alignment to the reference
# Flags excluded: 4 (unmapped), 256 (secondary), 2048 (supplementary)

samtools view -F 2308 hap_vs_ref.sam \
| cut -f1 \
| sort -u \
> mapped.txt

# Report filtered contigs
total=$(grep -c "^>" "$FILE" || true)
mapped_count=$(wc -l < mapped.txt)
filtered=$((total - mapped_count))
if [[ $filtered -gt 0 ]]; then
    echo "  WARNING: $filtered of $total contigs had no primary alignment to reference and were removed"
fi

# Subset to on-target sequences only
seqkit grep -f mapped.txt "$FILE" > mapped.fasta



# Detect reverse among mapped sequences

samtools view -f 16 -F 2308 hap_vs_ref.sam \
| cut -f1 \
| sort -u \
> to_reverse.txt



# Reverse if needed

if [[ -s to_reverse.txt ]]; then


seqkit grep -f to_reverse.txt mapped.fasta \
| seqkit seq -r -p \
> reversed.fasta


seqkit grep -v -f to_reverse.txt mapped.fasta \
> forward.fasta


cat forward.fasta reversed.fasta \
> "$ORIENTED"


else


cp mapped.fasta "$ORIENTED"


fi


rm -f hap_vs_ref.sam forward.fasta reversed.fasta mapped.fasta mapped.txt


cat "$ORIENTED" >> "$OLDPWD/$FINAL_OUT"


cd "$OLDPWD"


done



echo
echo "DONE"
echo "Output: $FINAL_OUT"
