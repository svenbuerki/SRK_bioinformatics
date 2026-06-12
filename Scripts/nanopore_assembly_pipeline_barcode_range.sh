#!/bin/bash
# Note: deliberately NO `set -e` — per-sample failures must not abort the loop.
# `set -u` catches typos in variable names; `pipefail` catches mid-pipeline failures.
set -uo pipefail

# ==========================================
# Nanopore Amplicon Assembly & Phasing Script
# Multi-CANU Rare Haplotype Optimized Version
# WITH FASTQ CONCATENATION
# Failure-tolerant: a per-step failure aborts the affected sample only,
# the loop continues with the next sample, and a summary log records
# which samples succeeded and which step (if any) caused each failure.
# ==========================================

orig_dir=$(pwd)

read -p "Enter library name (e.g., Library001): " library
read -p "Enter starting barcode number (e.g., 01): " start_barcode
read -p "Enter ending barcode number (e.g., 10): " end_barcode

start_num=$((10#$start_barcode))
end_num=$((10#$end_barcode))

# Number of CANU replicates
REPS=4

# Embedded coverage-based chimera filter — see chimera_coverage_filter.py.
# Default mode is "report" so a calibration run on Library 005 produces
# the per-contig stats TSV without removing anything. Override via env vars:
#   CHIMERA_FILTER_MODE=filter MIN_MEAN_COV=20 MIN_UNIFORMITY=0.2 \
#   CHIMERA_WINDOW=200 ./nanopore_assembly_pipeline_barcode_range.sh
CHIMERA_FILTER_MODE="${CHIMERA_FILTER_MODE:-report}"
MIN_MEAN_COV="${MIN_MEAN_COV:-20}"
MIN_UNIFORMITY="${MIN_UNIFORMITY:-0.2}"
CHIMERA_WINDOW="${CHIMERA_WINDOW:-200}"
CHIMERA_FILTER_SCRIPT="$orig_dir/chimera_coverage_filter.py"

# Auto-detect a working Python interpreter (override via `export PYTHON_BIN=...`).
if [[ -z "${PYTHON_BIN:-}" ]]; then
    if command -v python3 &>/dev/null; then
        PYTHON_BIN="$(command -v python3)"
    elif command -v python &>/dev/null; then
        PYTHON_BIN="$(command -v python)"
    elif [[ -x /Users/sven/anaconda3/bin/python ]]; then
        PYTHON_BIN="/Users/sven/anaconda3/bin/python"
    else
        echo "ERROR: no python3 or python found on PATH; set PYTHON_BIN explicitly" >&2
        exit 1
    fi
fi

# Canu resume logic: by default, if CANU_rep<N>/*.contigs.fasta exists and is
# non-empty for a given rep, that rep is skipped (re-running with new chimera
# thresholds or restarting from an interrupted run is then much faster).
# Override with FORCE_CANU=1 to wipe and re-run every Canu rep from scratch.
FORCE_CANU="${FORCE_CANU:-0}"

# ==========================================
# Logging setup
# ==========================================
run_id="$(date +%Y%m%d_%H%M%S)"
log_dir="$orig_dir/logs"
mkdir -p "$log_dir"

summary_log="$log_dir/${library}_${start_barcode}-${end_barcode}_${run_id}.summary.tsv"
printf "sample\tstatus\tfailed_step\tlog_file\n" > "$summary_log"

echo ""
echo "Per-sample logs    : $log_dir/${library}_${run_id}_<barcode>.log"
echo "Run summary table  : $summary_log"
echo "Chimera filter     : mode=$CHIMERA_FILTER_MODE  "\
"min_mean_cov=$MIN_MEAN_COV  min_uniformity=$MIN_UNIFORMITY  window=${CHIMERA_WINDOW}bp"
echo "Canu resume        : FORCE_CANU=$FORCE_CANU"\
" (skips reps with existing CANU_rep<N>/*.contigs.fasta unless forced)"
echo "Python interpreter : $PYTHON_BIN"
echo ""

# Run one pipeline step, capturing stdout+stderr to the per-sample log.
# On success: returns 0.
# On failure: prints the step name to stdout (for capture by the caller) and returns the original exit code.
try_step() {
    local step_name="$1"; shift
    local sample_log="$1"; shift

    echo "" >> "$sample_log"
    echo "[$(date +%H:%M:%S)] STEP: $step_name" >> "$sample_log"

    "$@" >> "$sample_log" 2>&1
    local rc=$?

    if [[ $rc -ne 0 ]]; then
        echo "[$(date +%H:%M:%S)] FAIL: $step_name (exit $rc)" >> "$sample_log"
        return $rc
    fi
    return 0
}

# ==========================================
# Main loop
# ==========================================
for n in $(seq "$start_num" "$end_num"); do

    sample=$(printf "barcode%02d" "$n")

    echo ""
    echo "========================================="
    echo "Processing $library / $sample"
    echo "========================================="

    sample_log="$log_dir/${library}_${run_id}_${sample}.log"
    : > "$sample_log"
    failed_step=""

    # --- enter sample dir ---
    if ! cd "$orig_dir/$library/$sample" 2>>"$sample_log"; then
        echo "Directory $library/$sample not found — skipping" | tee -a "$sample_log"
        printf "%s\tFAIL\tenter_directory\t%s\n" "$sample" "$sample_log" >> "$summary_log"
        cd "$orig_dir"
        continue
    fi

    # --- skip sample if final output already exists (override with FORCE_SAMPLE=1) ---
    final_out="${library}_${sample}_Phased_haplotypes.fasta"
    if [[ "${FORCE_SAMPLE:-0}" != "1" && -s "$final_out" ]]; then
        echo "[SKIP] $library / $sample — $final_out already exists" | tee -a "$sample_log"
        printf "%s\tSKIP\texisting_output\t%s\n" "$sample" "$sample_log" >> "$summary_log"
        cd "$orig_dir"
        continue
    fi

    reads="${library}_${sample}.fastq.gz"

    # ==========================================
    # FASTQ CONCATENATION LOGIC
    # ==========================================
    if [[ ! -f "$reads" ]]; then
        echo "Concatenated FASTQ not found. Creating $reads..." | tee -a "$sample_log"

        if ls *.fastq.gz 1> /dev/null 2>&1; then
            gz_count=$(ls *.fastq.gz 2>/dev/null | grep -v "^${reads}$" | wc -l || echo "0")

            if [[ $gz_count -gt 0 ]]; then
                echo "Found $gz_count .fastq.gz files to concatenate" >> "$sample_log"
                if ! cat $(ls *.fastq.gz | grep -v "^${reads}$") > "$reads" 2>>"$sample_log"; then
                    echo "ERROR: FASTQ concatenation failed" >> "$sample_log"
                    failed_step="fastq_concat"
                fi
            else
                echo "ERROR: No .fastq.gz files found to concatenate in $(pwd)" >> "$sample_log"
                ls -la >> "$sample_log" 2>&1
                failed_step="no_input_fastq"
            fi
        else
            echo "ERROR: No .fastq.gz files found in $(pwd)" >> "$sample_log"
            ls -la >> "$sample_log" 2>&1
            failed_step="no_input_fastq"
        fi
    else
        echo "Using existing $reads" >> "$sample_log"
    fi

    # Final check that the reads file exists and is not empty
    if [[ -z "$failed_step" ]] && [[ ! -f "$reads" ]]; then
        echo "ERROR: $reads still not found after concatenation attempt" >> "$sample_log"
        failed_step="reads_missing"
    fi
    if [[ -z "$failed_step" ]] && [[ ! -s "$reads" ]]; then
        echo "ERROR: $reads exists but is empty" >> "$sample_log"
        failed_step="reads_empty"
    fi

    if [[ -z "$failed_step" ]]; then
        echo "Using reads file: $reads ($(du -h "$reads" | cut -f1))" | tee -a "$sample_log"
    fi

    # ==========================================
    # Run multiple CANU assemblies
    # ==========================================
    if [[ -z "$failed_step" ]]; then
        for rep in $(seq 1 $REPS); do
            # Resume check — skip this Canu rep if a non-empty contigs.fasta
            # already exists in CANU_rep<rep>/, unless FORCE_CANU=1.
            existing_contigs=$(ls "CANU_rep${rep}"/*.contigs.fasta 2>/dev/null | head -1 || true)
            if [[ "$FORCE_CANU" != "1" && -n "$existing_contigs" && -s "$existing_contigs" ]]; then
                echo "[SKIP] canu_rep${rep} — found existing $existing_contigs" \
                    | tee -a "$sample_log"
                continue
            fi
            echo "Running CANU replicate $rep" | tee -a "$sample_log"
            if ! try_step "canu_rep${rep}" "$sample_log" canu \
                -p "${library}_${sample}_rep${rep}" \
                -d "CANU_rep${rep}" \
                genomeSize=5k \
                -nanopore "$reads" \
                useGrid=false \
                correctedErrorRate=0.15 \
                minInputCoverage=0 \
                stopOnLowCoverage=0 \
                corOutCoverage=10000 \
                batOptions="-dg 3 -db 3 -dr 1 -ca 500 -cp 50"; then
                failed_step="canu_rep${rep}"
                break
            fi
        done
    fi

    # ==========================================
    # Combine contigs with proper renaming
    # ==========================================
    if [[ -z "$failed_step" ]]; then
        echo "Combining contigs" | tee -a "$sample_log"
        : > combined_contigs.fasta
        rep=1
        combine_ok=1
        for dir in CANU_rep*; do
            fasta=$(ls "$dir"/*.contigs.fasta 2>/dev/null || true)
            if [[ -z "$fasta" ]]; then
                echo "ERROR: no contigs.fasta in $dir" >> "$sample_log"
                combine_ok=0
                break
            fi
            if ! awk -v prefix="rep${rep}_" '
                /^>/ {print ">" prefix substr($0,2); next}
                {print}
                ' "$fasta" >> combined_contigs.fasta 2>>"$sample_log"; then
                combine_ok=0
                break
            fi
            ((rep++))
        done
        if [[ $combine_ok -eq 0 ]]; then
            failed_step="combine_contigs"
        fi
    fi

    if [[ -z "$failed_step" ]]; then
        if ! try_step "seqkit_rmdup" "$sample_log" \
            bash -c 'seqkit rmdup -s combined_contigs.fasta > combined_unique_contigs.fasta'; then
            failed_step="seqkit_rmdup"
        fi
    fi

    assembly="combined_unique_contigs.fasta"

    # ==========================================
    # Embedded coverage-based chimera filter (see chimera_coverage_filter.py)
    # ==========================================
    # Align reads to the Canu contigs, compute per-contig coverage stats, then
    # (report mode) write the stats TSV / (filter mode) drop chimeric contigs.
    # The stats TSV per sample lands at
    #   $log_dir/${library}_${run_id}_${sample}.chimera_stats.tsv
    # so a calibration pass on Library 005 produces one TSV per barcode you can
    # inspect before flipping to filter mode.
    if [[ -z "$failed_step" ]]; then
        if ! try_step "chimera_minimap2" "$sample_log" \
            bash -c "minimap2 -ax map-ont -t 8 $assembly \"$reads\" \
                     | samtools sort -o chimera_coverage.bam -"; then
            failed_step="chimera_minimap2"
        fi
    fi
    if [[ -z "$failed_step" ]]; then
        if ! try_step "chimera_samtools_index" "$sample_log" \
            samtools index chimera_coverage.bam; then
            failed_step="chimera_samtools_index"
        fi
    fi
    if [[ -z "$failed_step" ]]; then
        if ! try_step "chimera_samtools_depth" "$sample_log" \
            bash -c "samtools depth -a chimera_coverage.bam > chimera_depth.tsv"; then
            failed_step="chimera_samtools_depth"
        fi
    fi
    if [[ -z "$failed_step" ]]; then
        chimera_stats_tsv="$log_dir/${library}_${run_id}_${sample}.chimera_stats.tsv"
        if ! try_step "chimera_coverage_filter" "$sample_log" \
            "$PYTHON_BIN" "$CHIMERA_FILTER_SCRIPT" \
                --depth chimera_depth.tsv \
                --in combined_unique_contigs.fasta \
                --out combined_unique_contigs_chimera_filtered.fasta \
                --stats "$chimera_stats_tsv" \
                --mode "$CHIMERA_FILTER_MODE" \
                --window "$CHIMERA_WINDOW" \
                --min-mean-cov "$MIN_MEAN_COV" \
                --min-uniformity "$MIN_UNIFORMITY"; then
            failed_step="chimera_coverage_filter"
        else
            # In filter mode, hand the filtered FASTA to RACON; in report mode
            # keep the unfiltered FASTA so the rest of the pipeline is unchanged.
            if [[ "$CHIMERA_FILTER_MODE" == "filter" ]]; then
                assembly="combined_unique_contigs_chimera_filtered.fasta"
                if [[ ! -s "$assembly" ]]; then
                    echo "ERROR: chimera filter dropped every contig" >> "$sample_log"
                    failed_step="chimera_filter_empty_output"
                fi
            fi
        fi
    fi

    # ==========================================
    # RACON polishing
    # ==========================================
    if [[ -z "$failed_step" ]]; then
        for round in 1 2 3; do
            input=$([ $round -eq 1 ] && echo "$assembly" || echo "polished_r$((round-1)).fasta")

            if ! try_step "minimap2_polish_r${round}" "$sample_log" \
                bash -c "minimap2 -ax map-ont \"$input\" \"$reads\" > align.sam"; then
                failed_step="minimap2_polish_r${round}"
                break
            fi
            if ! try_step "racon_r${round}" "$sample_log" \
                bash -c "racon \"$reads\" align.sam \"$input\" > polished_r${round}.fasta"; then
                failed_step="racon_r${round}"
                break
            fi
        done
    fi

    polished_final="polished_r3.fasta"

    # ==========================================
    # Align reads
    # ==========================================
    if [[ -z "$failed_step" ]]; then
        if ! try_step "align_reads" "$sample_log" \
            bash -c "minimap2 -ax map-ont -t 8 \"$polished_final\" \"$reads\" | samtools view -bS - | samtools sort -o aligned.bam -"; then
            failed_step="align_reads"
        fi
    fi
    if [[ -z "$failed_step" ]]; then
        if ! try_step "samtools_index" "$sample_log" samtools index aligned.bam; then
            failed_step="samtools_index"
        fi
    fi

    # ==========================================
    # Variant Calling
    # ==========================================
    if [[ -z "$failed_step" ]]; then
        if ! try_step "freebayes" "$sample_log" \
            bash -c "freebayes -f \"$polished_final\" -p 4 aligned.bam > variants.vcf"; then
            failed_step="freebayes"
        fi
    fi
    if [[ -z "$failed_step" ]]; then
        if ! try_step "bgzip_variants" "$sample_log" bgzip -f variants.vcf; then
            failed_step="bgzip_variants"
        fi
    fi
    if [[ -z "$failed_step" ]]; then
        if ! try_step "tabix_variants" "$sample_log" tabix -p vcf variants.vcf.gz; then
            failed_step="tabix_variants"
        fi
    fi

    # ==========================================
    # WhatsHap Phasing
    # ==========================================
    if [[ -z "$failed_step" ]]; then
        if ! try_step "whatshap_polyphase" "$sample_log" whatshap polyphase \
            --ploidy 4 \
            --reference "$polished_final" \
            --ignore-read-groups \
            -o phased.vcf \
            variants.vcf.gz aligned.bam; then
            failed_step="whatshap_polyphase"
        fi
    fi
    if [[ -z "$failed_step" ]]; then
        if ! try_step "bgzip_phased" "$sample_log" bgzip -f phased.vcf; then
            failed_step="bgzip_phased"
        fi
    fi
    if [[ -z "$failed_step" ]]; then
        if ! try_step "tabix_phased" "$sample_log" tabix -p vcf phased.vcf.gz; then
            failed_step="tabix_phased"
        fi
    fi

    # ==========================================
    # Export haplotypes
    # ==========================================
    if [[ -z "$failed_step" ]]; then
        for i in 1 2 3 4; do
            if ! try_step "bcftools_consensus_hap${i}" "$sample_log" \
                bash -c "bcftools consensus -H $i -f \"$polished_final\" phased.vcf.gz > temp.fasta"; then
                failed_step="bcftools_consensus_hap${i}"
                break
            fi
            if ! sed "s/^>/>${library}_${sample}_hap${i}_/" temp.fasta \
                > "${library}_${sample}_hap${i}.fasta" 2>>"$sample_log"; then
                failed_step="rename_hap${i}"
                break
            fi
        done
    fi

    if [[ -z "$failed_step" ]]; then
        if ! cat ${library}_${sample}_hap*.fasta \
            > "${library}_${sample}_Phased_haplotypes.fasta" 2>>"$sample_log"; then
            failed_step="merge_haplotypes"
        else
            rm -f temp.fasta
        fi
    fi

    # ==========================================
    # Record per-sample result
    # ==========================================
    if [[ -z "$failed_step" ]]; then
        echo "Finished $library $sample" | tee -a "$sample_log"
        printf "%s\tOK\t-\t%s\n" "$sample" "$sample_log" >> "$summary_log"
    else
        echo "FAILED: $library $sample at step '$failed_step' — see $sample_log" | tee -a "$sample_log"
        printf "%s\tFAIL\t%s\t%s\n" "$sample" "$failed_step" "$sample_log" >> "$summary_log"
    fi

    cd "$orig_dir"

done

# ==========================================
# Run summary
# ==========================================
echo ""
echo "====================================="
echo "RUN COMPLETE — SUMMARY"
echo "====================================="
column -t -s $'\t' "$summary_log" 2>/dev/null || cat "$summary_log"
echo ""
n_ok=$(awk -F'\t' 'NR>1 && $2 == "OK"   {c++} END {print c+0}' "$summary_log")
n_fail=$(awk -F'\t' 'NR>1 && $2 == "FAIL" {c++} END {print c+0}' "$summary_log")
echo "Successful samples : $n_ok"
echo "Failed samples     : $n_fail"
echo ""
echo "Detailed per-sample logs : $log_dir/"
echo "Summary table            : $summary_log"
