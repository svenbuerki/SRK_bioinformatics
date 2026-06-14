#!/bin/bash
set -euo pipefail

# ==========================================
# Nanopore Amplicon Assembly & Phasing Script
# ==========================================

# Prompt user for library and sample
read -p "Enter library name (e.g., Library_0001): " library
read -p "Enter sample name (e.g., barcode03): " sample

# Navigate to sample folder
cd "$library/$sample"

# ----------------------------
# Merge Nanopore FASTQ files
# ----------------------------
cat *.fastq.gz > "${library}_${sample}.fastq.gz"

# ----------------------------
# De novo assembly with CANU
# ----------------------------
canu \
  -p "${library}_${sample}" \
  -d "${library}_${sample}_CANU" \
  genomeSize=5k \
  -nanopore "${library}_${sample}.fastq.gz" \
  corOutCoverage=10000 \
  contigFilter="2 0 1.0 0.5 0" \
  correctedErrorRate=0.15 \
  batOptions="-dg 3 -db 3 -dr 1 -ca 500 -cp 50" \
  useGrid=false

# ----------------------------
# Edit contig headers
# ----------------------------
cd "${library}_${sample}_CANU"

for fasta in *.contigs.fasta; do
    tmp="${fasta}.tmp"
    awk -v prefix="${library}_${sample}_" '
        /^>/ { print ">" prefix substr($0,2); next }
        { print }
    ' "$fasta" > "$tmp" && mv "$tmp" "$fasta"
done

# ----------------------------
# Iterative RACON polishing
# ----------------------------
assembly="${library}_${sample}.contigs.fasta"
reads=../"${library}_${sample}.fastq.gz"

for round in 1 2 3; do
    polished="polished_r${round}.fasta"
    sam="alignments_r${round}.sam"
    bam="alignments_r${round}.bam"
    input_assembly=$([ $round -eq 1 ] && echo "$assembly" || echo "polished_r$((round-1)).fasta")

    minimap2 -ax map-ont "$input_assembly" "$reads" > "$sam"
    samtools view -bS "$sam" | samtools sort -o "$bam" -
    samtools index "$bam"
    racon "$reads" "$sam" "$input_assembly" > "$polished"
done

# Last polished assembly
polished_final="polished_r3.fasta"

# ----------------------------
# Align reads to polished assembly
# ----------------------------
minimap2 -ax map-ont -t 8 "$polished_final" "$reads" \
    | samtools view -bS - \
    | samtools sort -o aligned.bam - \
    && samtools index aligned.bam

samtools stats aligned.bam > stats.txt
samtools depth aligned.bam > aligned.coverage

# ----------------------------
# Variant calling (tetraploid)
# ----------------------------
freebayes -f "$polished_final" -p 4 aligned.bam > variants.vcf
bgzip -c variants.vcf > variants.vcf.gz
tabix -p vcf variants.vcf.gz

# ----------------------------
# Variant phasing with WhatsHap
# ----------------------------
whatshap polyphase \
    --ploidy 4 \
    --reference "$polished_final" \
    --ignore-read-groups \
    -o phased.vcf \
    variants.vcf.gz aligned.bam

bgzip -c phased.vcf > phased.vcf.gz
tabix -p vcf phased.vcf.gz

# ----------------------------
# Generate phased haplotype sequences
# ----------------------------
VCF="phased.vcf.gz"
REF="$polished_final"

for i in $(seq 1 4); do
    OUT="${library}_${sample}_Phased_hap${i}.fasta"
    bcftools consensus -H $i -f "$REF" "$VCF" > temp.fasta
    sed "s/^>/>hap${i}_/" temp.fasta > "$OUT"
    rm temp.fasta
done

# Concatenate all haplotypes into a single FASTA
cat *_Phased_hap*.fasta > "${library}_${sample}_Phased_haplotypes.fasta"

# Optional: remove individual haplotype files
# rm *_Phased_hap*.fasta

