#!/bin/bash
set -euo pipefail

# ==========================================
# Step 4b — Reference-similarity filter
# ==========================================
# Removes paralogs and chimeric assemblies that pass the Step 4 length
# filter but are not bona-fide SRK alleles. Discrimination is by BLAST
# query-coverage against the LEPA canonical SRK reference set.
#
# Rationale (Library 010, 2026-05-11): length-filtered FASTA contained a
# bimodal cluster — true SRK (~3300 bp) plus a paralog cluster (~3900 bp)
# carrying a ~500 bp insertion. Both clusters are 98-100% identical
# where they align, so identity does not discriminate. Coverage does:
# true SRK aligns over >=90% of query, paralogs over <=85%.
# ==========================================

# ----------------------------
# Inputs
# ----------------------------
read -p "Enter length-filtered FASTA (Step 4 output): " INPUT
read -p "Enter canonical SRK reference FASTA: " REF
read -p "Minimum query coverage [0-1, default 0.90]: " MINCOV
read -p "Minimum percent identity [0-100, default 95]: " MINID
read -p "BLAST threads [default 4]: " THREADS

MINCOV=${MINCOV:-0.90}
MINID=${MINID:-95}
THREADS=${THREADS:-4}

for f in "$INPUT" "$REF"; do
    [[ -f "$f" ]] || { echo "ERROR: file not found: $f" >&2; exit 1; }
done

OUT="${INPUT%.fasta}_blastfilt.fasta"
LOG="${INPUT%.fasta}_blastfilt_log.tsv"
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

# ----------------------------
# 1. Build BLAST db from reference
# ----------------------------
makeblastdb -in "$REF" -dbtype nucl -out "${TMPDIR}/refdb" >/dev/null

# ----------------------------
# 2. Run blastn (best hit per query)
# ----------------------------
blastn -query "$INPUT" -db "${TMPDIR}/refdb" \
       -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' \
       -max_target_seqs 1 -max_hsps 1 \
       -num_threads "$THREADS" \
       -out "${TMPDIR}/blast.tsv"

# ----------------------------
# 3. Per-query decision (best bitscore wins)
# ----------------------------
awk -F'\t' -v MINCOV="$MINCOV" -v MINID="$MINID" '
{
  qid=$1; pident=$3; alen=$4; qlen=$5; bits=$12;
  cov=alen/qlen;
  if (!(qid in bb) || bits > bb[qid]) {
    bb[qid]=bits; pp[qid]=pident; cc[qid]=cov; qq[qid]=qlen;
  }
}
END {
  for (q in bb) {
    decision = (cc[q] >= MINCOV && pp[q] >= MINID) ? "KEEP" : "DROP";
    printf "%s\t%d\t%.2f\t%.4f\t%s\n", q, qq[q], pp[q], cc[q], decision;
  }
}' "${TMPDIR}/blast.tsv" > "${TMPDIR}/decisions.tsv"

# Header + decisions (sorted by decision then qid)
{
  echo -e "query_id\tquery_length\tpct_identity\tcoverage\tdecision"
  sort -k5,5 -k1,1 "${TMPDIR}/decisions.tsv"
} > "$LOG"

# ----------------------------
# 4. Build keep-list (queries that passed + canonical refs)
# ----------------------------
awk -F'\t' '$5=="KEEP" {print $1}' "${TMPDIR}/decisions.tsv"  > "${TMPDIR}/keep.txt"
grep '^>' "$REF" | sed 's/^>//; s/[[:space:]].*//'           >> "${TMPDIR}/keep.txt"

# ----------------------------
# 5. Filter input FASTA
# ----------------------------
awk -v keepfile="${TMPDIR}/keep.txt" '
BEGIN {
  while ((getline line < keepfile) > 0) {
    sub(/[[:space:]].*$/, "", line);
    keep[line]=1;
  }
  close(keepfile);
  emit=0;
}
/^>/ {
  hdr=substr($0, 2); sub(/[[:space:]].*$/, "", hdr);
  emit = (hdr in keep) ? 1 : 0;
}
emit { print }
' "$INPUT" > "$OUT"

# ----------------------------
# 6. Summary
# ----------------------------
n_in=$(grep -c '^>' "$INPUT")
n_out=$(grep -c '^>' "$OUT")
n_dropped=$((n_in - n_out))

cat <<SUMMARY

Reference-similarity filter complete.
  Input sequences:    $n_in
  Kept (BLAST pass):  $n_out
  Dropped:            $n_dropped
  Thresholds:         min_coverage=$MINCOV  min_identity=$MINID

Outputs:
  Filtered FASTA:     $OUT
  Per-query log:      $LOG
SUMMARY
