#!/usr/bin/env bash
# ==========================================
# SRK within-library pipeline wrapper
# Runs Steps 2 → 3 → 4 → 4b → 5 → 6 → 7 → 7b → 8 for each library listed
# in a config TSV. Failure-tolerant: a failure in any step aborts the
# affected library only; the loop continues with the next library.
#
# Usage:
#   ./run_within_library_suite.sh [config.tsv]
#   (defaults to ./libraries.tsv)
#
# Config TSV format (tab-separated, no header, # for comments):
#   Library001    01    22
#   Library002    23    59
#   ...
#
# Env vars (optional):
#   CHIMERA_FILTER_MODE   report (default) | filter   (passed to Step 2)
#   MIN_MEAN_COV          default 20                  (passed to Step 2)
#   MIN_UNIFORMITY        default 0.2                 (passed to Step 2)
#   CHIMERA_WINDOW        default 200                 (passed to Step 2)
#   FORCE_CANU            default 0                   (passed to Step 2)
#   FORCE_SAMPLE          default 0  — re-run samples even if final FASTA exists (Step 2)
#   FORCE                 default 0  — re-run every step even if output exists
#   START_AT_STEP         resume from step N (default 2)
#   STOP_AT_STEP          stop after step N (default 8)
#   PYTHON_BIN            override Python interpreter
#
# Resume behaviour:
#   - If `libraries.tsv` does not exist, it is auto-generated from Library*/barcode??
#   - Per-library: each step's expected output is checked first; if present and
#     non-empty the step is SKIPPED. Set FORCE=1 to re-run anyway.
#   - Per-sample (within Step 2): if ${LIB}/${BC}/${LIB}_${BC}_Phased_haplotypes.fasta
#     exists, that barcode is skipped. Set FORCE_SAMPLE=1 to re-run.
# ==========================================

set -uo pipefail

orig_dir=$(pwd)
config_tsv="${1:-libraries.tsv}"

# --- Canonical resources (absolute paths) ---
REF_FASTA="$orig_dir/Canonical_sequences/SRK_canonical_haplotype_sequences_revcomp_ATGfixed.fasta"
AUGUSTUS_CSV="$orig_dir/Canonical_sequences/augustus_SRK_canonical_haplotype_sequences_revcomp_ATGfixed_arabidopsis.csv"
CANONICAL_REF_ID='SRK_BEA_hybrid_bp_hap1_p_ctg_fa|amp_1|h1tg000019l|1700582|3407'

# --- Standard pipeline thresholds (from Bioinformatics_pipeline.md) ---
MIN_LEN=3250
MAX_LEN=4000
MIN_COV=0.90
MIN_IDENT=95
BLAST_THREADS=4
MAX_N=9
MAX_N_FRAC=0.005
TRANS_FRAME=1
EXON_STRAND="+"
MIN_EXON_LEN=0

START_AT_STEP="${START_AT_STEP:-2}"
STOP_AT_STEP="${STOP_AT_STEP:-8}"

# --- Python interpreter auto-detect ---
if [[ -z "${PYTHON_BIN:-}" ]]; then
    if command -v python3 &>/dev/null; then PYTHON_BIN="$(command -v python3)"
    elif command -v python &>/dev/null; then PYTHON_BIN="$(command -v python)"
    elif [[ -x /Users/sven/anaconda3/bin/python ]]; then PYTHON_BIN="/Users/sven/anaconda3/bin/python"
    else echo "ERROR: no python3 or python on PATH; set PYTHON_BIN explicitly" >&2; exit 1; fi
fi

# --- Sanity checks on canonical resources ---
[[ -f "$REF_FASTA" ]]    || { echo "ERROR: canonical ref not found: $REF_FASTA" >&2; exit 1; }
[[ -f "$AUGUSTUS_CSV" ]] || { echo "ERROR: AUGUSTUS CSV not found: $AUGUSTUS_CSV" >&2; exit 1; }

# --- Auto-generate config TSV if missing ---
if [[ ! -f "$config_tsv" ]]; then
    echo "[INFO] $config_tsv not found — auto-generating from Library*/barcode?? layout..."
    {
        echo "# Auto-generated $(date '+%Y-%m-%d %H:%M:%S') from $orig_dir"
        echo "# Format: library<TAB>start_barcode<TAB>end_barcode"
        echo "# Edit / comment out lines (#) to control which libraries run."
        for lib_dir in Library*/; do
            [[ -d "$lib_dir" ]] || continue
            lib="${lib_dir%/}"
            bcs=$(ls -d "$lib"/barcode?? 2>/dev/null | sed 's/.*barcode//' | sort -n)
            [[ -z "$bcs" ]] && continue
            first=$(echo "$bcs" | head -1)
            last=$(echo "$bcs"  | tail -1)
            printf '%s\t%s\t%s\n' "$lib" "$first" "$last"
        done
    } > "$config_tsv"
    n_libs=$(grep -vc '^#' "$config_tsv" || true)
    echo "[INFO] Wrote $config_tsv with $n_libs libraries:"
    grep -v '^#' "$config_tsv" | sed 's/^/    /'
    echo ""
fi
[[ -s "$config_tsv" ]] || { echo "ERROR: $config_tsv empty / no Library*/barcode?? directories found" >&2; exit 1; }

# --- Logging setup ---
run_id="$(date +%Y%m%d_%H%M%S)"
master_log_dir="$orig_dir/logs/wrapper_${run_id}"
mkdir -p "$master_log_dir"
master_summary="$master_log_dir/master_summary.tsv"
echo -e "library\tstatus\tlast_completed_step\tfailed_step\tlog_file" > "$master_summary"

echo "==========================================="
echo "SRK within-library wrapper — run $run_id"
echo "Config: $config_tsv"
echo "Logs:   $master_log_dir"
echo "Steps:  $START_AT_STEP → $STOP_AT_STEP"
echo "Python: $PYTHON_BIN"
echo "==========================================="

# ----------------------------------------------------------------------
# try_step: run a pipeline step, log output, return its exit code.
#   $1 = step label (e.g. "4b")
#   $2 = per-library log path
#   $3+ = command + args to run
# ----------------------------------------------------------------------
try_step() {
    local step="$1"; shift
    local lib_log="$1"; shift
    {
        echo ""
        echo "----- [STEP $step] $(date '+%Y-%m-%d %H:%M:%S') -----"
        echo "CMD: $*"
    } >> "$lib_log"
    "$@" >> "$lib_log" 2>&1
    local rc=$?
    echo "[STEP $step] exit=$rc" >> "$lib_log"
    return $rc
}

# ----------------------------------------------------------------------
# check_output: verify a step produced its expected output file.
# ----------------------------------------------------------------------
check_output() {
    local step="$1"; shift
    local path="$1"
    if [[ ! -s "$path" ]]; then
        echo "[STEP $step] MISSING EXPECTED OUTPUT: $path" >&2
        return 1
    fi
    return 0
}

# ----------------------------------------------------------------------
# already_done: returns 0 if expected output exists and non-empty (and
# FORCE != 1), else returns 1. Used to skip already-completed steps.
#   $1 = step label, $2 = expected output path, $3 = per-library log
# ----------------------------------------------------------------------
already_done() {
    local step="$1" path="$2" lib_log="$3"
    if [[ "${FORCE:-0}" != "1" && -s "$path" ]]; then
        echo "[SKIP step $step] output exists: $path" | tee -a "$lib_log"
        return 0
    fi
    return 1
}

# ======================================================================
# Main loop over libraries
# ======================================================================
while IFS=$'\t' read -r library start_bc end_bc; do
    # Skip blanks and comments
    [[ -z "${library:-}" || "$library" =~ ^[[:space:]]*# ]] && continue

    cd "$orig_dir"

    lib_log="$master_log_dir/${library}.log"
    {
        echo "========================================="
        echo "Library: $library"
        echo "Barcodes: ${start_bc}-${end_bc}"
        echo "Start:    $(date)"
        echo "========================================="
    } | tee -a "$lib_log"

    last_step="none"
    failed_step=""

    # Filenames evolve through the pipeline — compute them up front.
    merged="all_${library}_Phased_haplotypes.fasta"
    lenfilt="all_${library}_Phased_haplotypes_filtered_min${MIN_LEN}_max${MAX_LEN}.fasta"
    blastfilt="${lenfilt%.fasta}_blastfilt.fasta"
    aligned="${blastfilt%.fasta}_aligned.fasta"
    exons="${aligned%.fasta}_exons.fasta"
    backfilled="${exons%.fasta}_backfilled.fasta"
    nfilt="${backfilled%.fasta}_Nfilt.fasta"

    # --- Step 2: assembly + chimera filter ---
    if [[ -z "$failed_step" && $START_AT_STEP -le 2 && $STOP_AT_STEP -ge 2 ]]; then
        echo "[INFO] Step 2 (assembly) — $library $start_bc-$end_bc" | tee -a "$lib_log"
        if printf '%s\n%s\n%s\n' "$library" "$start_bc" "$end_bc" \
              | try_step "2" "$lib_log" bash "$orig_dir/nanopore_assembly_pipeline_barcode_range.sh"; then
            last_step="2"
        else
            failed_step="2_assembly"
        fi
    fi

    # --- Step 3: orient & consolidate ---
    if [[ -z "$failed_step" && $START_AT_STEP -le 3 && $STOP_AT_STEP -ge 3 ]]; then
        if already_done "3" "$merged" "$lib_log"; then
            last_step="3"
        else
            echo "[INFO] Step 3 (orient) — $library" | tee -a "$lib_log"
            if printf '%s\n%s\n' "$library" "$REF_FASTA" \
                  | try_step "3" "$lib_log" bash "$orig_dir/orient_haplotypes_library.sh" \
                  && check_output "3" "$merged"; then
                last_step="3"
            else
                failed_step="3_orient"
            fi
        fi
    fi

    # --- Step 4: length filter ---
    if [[ -z "$failed_step" && $START_AT_STEP -le 4 && $STOP_AT_STEP -ge 4 ]]; then
        if already_done "4" "$lenfilt" "$lib_log"; then
            last_step="4"
        else
            echo "[INFO] Step 4 (length filter) — $library" | tee -a "$lib_log"
            if printf '%s\n%s\n%s\n%s\n' "$merged" "$REF_FASTA" "$MIN_LEN" "$MAX_LEN" \
                  | try_step "4" "$lib_log" bash "$orig_dir/filter_and_add_references.sh" \
                  && check_output "4" "$lenfilt"; then
                last_step="4"
            else
                failed_step="4_length"
            fi
        fi
    fi

    # --- Step 4b: BLAST coverage filter ---
    if [[ -z "$failed_step" && $START_AT_STEP -le 4 && $STOP_AT_STEP -ge 4 ]]; then
        if already_done "4b" "$blastfilt" "$lib_log"; then
            last_step="4b"
        else
            echo "[INFO] Step 4b (BLAST filter) — $library" | tee -a "$lib_log"
            if printf '%s\n%s\n%s\n%s\n%s\n' "$lenfilt" "$REF_FASTA" "$MIN_COV" "$MIN_IDENT" "$BLAST_THREADS" \
                  | try_step "4b" "$lib_log" bash "$orig_dir/filter_by_reference_similarity.sh" \
                  && check_output "4b" "$blastfilt"; then
                last_step="4b"
            else
                failed_step="4b_blast"
            fi
        fi
    fi

    # --- Step 5: MAFFT alignment ---
    if [[ -z "$failed_step" && $START_AT_STEP -le 5 && $STOP_AT_STEP -ge 5 ]]; then
        if already_done "5" "$aligned" "$lib_log"; then
            last_step="5"
        else
            echo "[INFO] Step 5 (MAFFT) — $library" | tee -a "$lib_log"
            if try_step "5" "$lib_log" bash -c "mafft --auto --adjustdirection '$blastfilt' > '$aligned'" \
                  && check_output "5" "$aligned"; then
                last_step="5"
            else
                failed_step="5_mafft"
            fi
        fi
    fi

    # --- Step 6: exon extraction ---
    if [[ -z "$failed_step" && $START_AT_STEP -le 6 && $STOP_AT_STEP -ge 6 ]]; then
        if already_done "6" "$exons" "$lib_log"; then
            last_step="6"
        else
            echo "[INFO] Step 6 (exon extraction) — $library" | tee -a "$lib_log"
            if printf '%s\n%s\n%s\n%s\n%s\n' "$aligned" "$CANONICAL_REF_ID" "$AUGUSTUS_CSV" "$EXON_STRAND" "$MIN_EXON_LEN" \
                  | try_step "6" "$lib_log" "$PYTHON_BIN" "$orig_dir/extract_exons_with_annotation.py" \
                  && check_output "6" "$exons"; then
                last_step="6"
            else
                failed_step="6_exons"
            fi
        fi
    fi

    # --- Step 7: gap backfill ---
    if [[ -z "$failed_step" && $START_AT_STEP -le 7 && $STOP_AT_STEP -ge 7 ]]; then
        if already_done "7" "$backfilled" "$lib_log"; then
            last_step="7"
        else
            echo "[INFO] Step 7 (backfill) — $library" | tee -a "$lib_log"
            if printf '%s\n%s\n' "$exons" "$CANONICAL_REF_ID" \
                  | try_step "7" "$lib_log" "$PYTHON_BIN" "$orig_dir/backfill_alignment_ends.py" \
                  && check_output "7" "$backfilled"; then
                last_step="7"
            else
                failed_step="7_backfill"
            fi
        fi
    fi

    # --- Step 7b: N-content filter ---
    if [[ -z "$failed_step" && $START_AT_STEP -le 7 && $STOP_AT_STEP -ge 7 ]]; then
        if already_done "7b" "$nfilt" "$lib_log"; then
            last_step="7b"
        else
            echo "[INFO] Step 7b (N-filter) — $library" | tee -a "$lib_log"
            if printf '%s\n%s\n%s\n' "$backfilled" "$MAX_N" "$MAX_N_FRAC" \
                  | try_step "7b" "$lib_log" bash "$orig_dir/filter_by_ambiguity.sh" \
                  && check_output "7b" "$nfilt"; then
                last_step="7b"
            else
                failed_step="7b_Nfilt"
            fi
        fi
    fi

    # --- Step 8: AA translation (optional but standard) ---
    if [[ -z "$failed_step" && $START_AT_STEP -le 8 && $STOP_AT_STEP -ge 8 ]]; then
        aa_out="${nfilt%.fasta}_frame${TRANS_FRAME}_AA_filtered_aligned.fasta"
        if already_done "8" "$aa_out" "$lib_log"; then
            last_step="8"
        else
            echo "[INFO] Step 8 (AA translation) — $library" | tee -a "$lib_log"
            if printf '%s\n%s\n' "$nfilt" "$TRANS_FRAME" \
                  | try_step "8" "$lib_log" "$PYTHON_BIN" "$orig_dir/translate_filter_align_AA.py" \
                  && check_output "8" "$aa_out"; then
                last_step="8"
            else
                failed_step="8_AA"
            fi
        fi
    fi

    # --- Record library result ---
    if [[ -z "$failed_step" ]]; then
        status="OK"
        echo "[DONE] $library — completed through Step $last_step" | tee -a "$lib_log"
    else
        status="FAIL"
        echo "[FAIL] $library — failed at $failed_step (last OK = Step $last_step)" | tee -a "$lib_log"
    fi
    echo -e "${library}\t${status}\t${last_step}\t${failed_step:-}\t${lib_log}" >> "$master_summary"
    echo "End: $(date)" | tee -a "$lib_log"

done < "$config_tsv"

# ======================================================================
# Final report
# ======================================================================
echo ""
echo "==========================================="
echo "MASTER SUMMARY ($master_summary)"
echo "==========================================="
column -t -s $'\t' "$master_summary"

n_ok=$(awk -F'\t'   'NR>1 && $2=="OK"   {c++} END {print c+0}' "$master_summary")
n_fail=$(awk -F'\t' 'NR>1 && $2=="FAIL" {c++} END {print c+0}' "$master_summary")
echo ""
echo "Total: $n_ok OK / $n_fail FAIL"
echo "Per-library logs: $master_log_dir/<library>.log"
