#!/bin/bash
set -euo pipefail

# ==========================================
# Nanopore Amplicon Assembly & Phasing Script
# Multi-CANU Rare Haplotype Optimized Version
# WITH FASTQ CONCATENATION
# ==========================================

orig_dir=$(pwd)

read -p "Enter library name (e.g., Library001): " library
read -p "Enter starting barcode number (e.g., 01): " start_barcode
read -p "Enter ending barcode number (e.g., 10): " end_barcode

start_num=$((10#$start_barcode))
end_num=$((10#$end_barcode))

# Number of CANU replicates
REPS=4

for n in $(seq "$start_num" "$end_num"); do

    sample=$(printf "barcode%02d" "$n")

    echo ""
    echo "========================================="
    echo "Processing $library / $sample"
    echo "========================================="

    cd "$orig_dir/$library/$sample"

    reads="${library}_${sample}.fastq.gz"

    # ==========================================
    # FASTQ CONCATENATION LOGIC
    # ==========================================
    
    if [[ ! -f "$reads" ]]; then
        echo "Concatenated FASTQ not found. Creating $reads..."
        
        # Check if there are individual .fastq.gz files to concatenate
        if ls *.fastq.gz 1> /dev/null 2>&1; then
            # Count existing .fastq.gz files (excluding the target file if it exists)
            gz_count=$(ls *.fastq.gz 2>/dev/null | grep -v "^${reads}$" | wc -l || echo "0")
            
            if [[ $gz_count -gt 0 ]]; then
                echo "Found $gz_count .fastq.gz files to concatenate"
                
                # Concatenate all .fastq.gz files except the target file
                cat $(ls *.fastq.gz | grep -v "^${reads}$") > "$reads"
                
                echo "Created $reads from individual files"
            else
                echo "ERROR: No .fastq.gz files found to concatenate in $(pwd)"
                echo "Available files:"
                ls -la
                exit 1
            fi
        else
            echo "ERROR: No .fastq.gz files found in $(pwd)"
            echo "Available files:"
            ls -la
            exit 1
        fi
    else
        echo "Using existing $reads"
        
        # Optional: Check if the existing file is newer than individual files
        if ls *.fastq.gz 1> /dev/null 2>&1; then
            individual_files=$(ls *.fastq.gz | grep -v "^${reads}$" || true)
            if [[ -n "$individual_files" ]]; then
                # Check if any individual file is newer than the concatenated file
                newer_files=$(find . -name "*.fastq.gz" -not -name "$reads" -newer "$reads" 2>/dev/null || true)
                if [[ -n "$newer_files" ]]; then
                    echo "WARNING: Some individual .fastq.gz files are newer than $reads"
                    echo "Consider regenerating the concatenated file"
                fi
            fi
        fi
    fi

    # Final check that the reads file exists and is not empty
    if [[ ! -f "$reads" ]]; then
        echo "ERROR: $reads still not found after concatenation attempt"
        exit 1
    fi
    
    if [[ ! -s "$reads" ]]; then
        echo "ERROR: $reads exists but is empty"
        exit 1
    fi
    
    echo "Using reads file: $reads ($(du -h "$reads" | cut -f1))"

    # ==========================================
    # Run multiple CANU assemblies
    # ==========================================

    for rep in $(seq 1 $REPS); do

        echo ""
        echo "Running CANU replicate $rep"

        canu \
        -p "${library}_${sample}_rep${rep}" \
        -d "CANU_rep${rep}" \
        genomeSize=5k \
        -nanopore "$reads" \
        useGrid=false \
        correctedErrorRate=0.15 \
        minInputCoverage=0 \
        stopOnLowCoverage=0 \
        corOutCoverage=10000 \
        batOptions="-dg 3 -db 3 -dr 1 -ca 500 -cp 50"

    done

    # ==========================================
    # Combine contigs with proper renaming
    # ==========================================

    echo ""
    echo "Combining contigs"

    > combined_contigs.fasta

    rep=1

    for dir in CANU_rep*; do

        fasta=$(ls $dir/*.contigs.fasta)

        awk -v prefix="rep${rep}_" '
        /^>/ {print ">" prefix substr($0,2); next}
        {print}
        ' "$fasta" >> combined_contigs.fasta

        ((rep++))

    done

    echo "Removing duplicate sequences"

    seqkit rmdup -s combined_contigs.fasta \
    > combined_unique_contigs.fasta

    assembly="combined_unique_contigs.fasta"

    # ==========================================
    # RACON polishing
    # ==========================================

    echo ""
    echo "Running RACON polishing"

    for round in 1 2 3; do

        input=$([ $round -eq 1 ] && echo "$assembly" || echo "polished_r$((round-1)).fasta")

        minimap2 -ax map-ont "$input" "$reads" > align.sam

        racon "$reads" align.sam "$input" > polished_r${round}.fasta

    done

    polished_final="polished_r3.fasta"

    # ==========================================
    # Align reads
    # ==========================================

    echo ""
    echo "Aligning reads"

    minimap2 -ax map-ont -t 8 "$polished_final" "$reads" | \
    samtools view -bS - | \
    samtools sort -o aligned.bam -

    samtools index aligned.bam

    # ==========================================
    # Variant Calling
    # ==========================================

    echo ""
    echo "Calling variants"

    freebayes \
    -f "$polished_final" \
    -p 4 \
    aligned.bam > variants.vcf

    bgzip -f variants.vcf

    tabix -p vcf variants.vcf.gz

    # ==========================================
    # WhatsHap Phasing
    # ==========================================

    echo ""
    echo "Phasing variants"

    whatshap polyphase \
    --ploidy 4 \
    --reference "$polished_final" \
    --ignore-read-groups \
    -o phased.vcf \
    variants.vcf.gz aligned.bam

    bgzip -f phased.vcf

    tabix -p vcf phased.vcf.gz

    # ==========================================
    # Export haplotypes
    # ==========================================

    echo ""
    echo "Exporting haplotypes"

    for i in 1 2 3 4; do

        bcftools consensus \
        -H $i \
        -f "$polished_final" \
        phased.vcf.gz \
        > temp.fasta

        sed "s/^>/>${library}_${sample}_hap${i}_/" temp.fasta \
        > "${library}_${sample}_hap${i}.fasta"

    done

    cat ${library}_${sample}_hap*.fasta \
    > ${library}_${sample}_Phased_haplotypes.fasta

    rm temp.fasta

    echo ""
    echo "Finished $library $sample"

    cd "$orig_dir"

done

echo ""
echo "====================================="
echo "ALL SAMPLES COMPLETE"
echo "====================================="

