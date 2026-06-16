#!/usr/bin/env python3
"""
reseq_calibration.py — Nanopore re-sequencing calibration + per-sample failure
diagnostics + recommended-action assignment.

For every Library*/barcode??/ on disk, pulls:
  - Canu rep1 input read count (post seqkit-rmdup) from CANU_rep1/*.report.
  - Raw FASTQ stats (count, total bases, median / p10 / p90 length, fraction
    of reads ≥ 3000 bp) from the merged Library###_barcode##.fastq.gz.
  - Per-stage pipeline diagnostics: per-rep Canu contig counts,
    combined_contigs.fasta count, combined_unique_contigs.fasta count
    (post seqkit rmdup), chimera_stats.tsv KEEP vs DROP counts +
    uniformity summary (if logs/*chimera_stats.tsv exists).

Then joins to:
  - Tables/Phase4/step22b_individual_SI_status.tsv  (haplotype counts + SI call)
  - Tables/Phase2/step12c_samples_redo.csv          (Lab_action, if flagged)

Classifies every sample into a `failure_mode` and assigns a
`recommended_action`. Both columns are written into the per-sample TSV so the
lab + bioinformatics team can sort and act directly.

Writes:
  Tables/Phase4/step22c_reseq_calibration_per_sample.tsv — joined per-sample table
  Tables/Phase4/step22c_reseq_calibration_summary.tsv    — bin-level dose-response
  figures/Phase4/step22c_reseq_calibration.pdf/.png      — 3-panel figure
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

AMPLICON_BP = 3000           # reads >= this length can span an SRK amplicon (~3.5 kb)
STEP4_MIN_BP = 3250          # Step 4 length filter — lower bound
STEP4_MAX_BP = 4000          # Step 4 length filter — upper bound

HERE = Path(__file__).resolve().parent
TABLES = HERE / "Tables"
FIGS   = HERE / "figures" / "Phase4"

IN_SI    = TABLES / "Phase4" / "step22b_individual_SI_status.tsv"
IN_REDO  = TABLES / "Phase2" / "step12c_samples_redo.csv"
OUT_PER  = TABLES / "Phase4" / "step22c_reseq_calibration_per_sample.tsv"
OUT_SUMM = TABLES / "Phase4" / "step22c_reseq_calibration_summary.tsv"
OUT_PDF  = FIGS / "step22c_reseq_calibration.pdf"
OUT_PNG  = FIGS / "step22c_reseq_calibration.png"

# Create every output directory if it doesn't exist yet (lab machine may not
# have the Phase-prefixed tree set up).
for d in [TABLES, TABLES / "Phase4", FIGS]:
    d.mkdir(parents=True, exist_ok=True)

# Verify required input tables exist with a friendly error if they don't —
# the lab machine may be on the older pre-2026-06-14 layout.
for label, path in [("Step 22b SI status", IN_SI),
                    ("Step 12c samples_redo", IN_REDO)]:
    if not path.exists():
        raise SystemExit(
            f"ERROR: required input not found:\n  {label}: {path}\n"
            f"Hint: copy this file from the dev workstation into the same path "
            f"on this machine, or update the IN_SI / IN_REDO constants at the top "
            f"of this script.")

print(f"Working directory : {HERE}")
print(f"Inputs:")
print(f"  {IN_SI}")
print(f"  {IN_REDO}")
print(f"Outputs:")
print(f"  {OUT_PER}")
print(f"  {OUT_SUMM}")
print(f"  {OUT_PDF}")
print(f"  {OUT_PNG}")
print()

# ── Pull input read count from each Canu rep1 .report ────────────────────────
READ_LINE_RE = re.compile(r"Found\s+([\d,]+)\s+reads", re.IGNORECASE)

def reads_from_report(path: Path) -> int | None:
    try:
        with open(path) as fh:
            for line in fh:
                m = READ_LINE_RE.search(line)
                if m:
                    return int(m.group(1).replace(",", ""))
    except FileNotFoundError:
        return None
    return None

# ── Raw FASTQ stats — count + length distribution from the merged .fastq.gz ──
def _n50(lengths_desc: np.ndarray) -> int:
    """N50 = the length such that 50 % of total bases sit in reads ≥ that length.
    Assumes lengths_desc is already sorted descending."""
    total = int(lengths_desc.sum())
    if total == 0:
        return 0
    half = total / 2
    cum = 0
    for L in lengths_desc:
        cum += int(L)
        if cum >= half:
            return int(L)
    return int(lengths_desc[-1])

def fastq_stats(path: Path) -> dict | None:
    """One pass over the gzipped FASTQ. Returns count + length quantiles + N50,
    fraction of reads >= AMPLICON_BP, and read-level quality stats (mean Phred
    across all bases, mean of per-read mean Phred).
    Truncated gzip files are tolerated — returns whatever was read before EOF."""
    if not path.exists():
        return None
    lengths: list[int] = []
    # incremental Phred accumulators (over ALL bases of ALL reads in this file)
    phred_sum_all_bases = 0
    phred_n_all_bases   = 0
    # per-read mean Phred — mean of those means is what most QC tools report
    per_read_mean_phred_sum = 0.0
    per_read_count          = 0
    try:
        with gzip.open(path, "rb") as fh:
            i = 0
            try:
                for line in fh:
                    mod = i % 4
                    if mod == 1:                       # sequence line
                        lengths.append(len(line.rstrip(b"\n")))
                    elif mod == 3:                     # quality line
                        q = line.rstrip(b"\n")
                        if q:
                            arr = np.frombuffer(q, dtype=np.uint8)
                            # Phred = ASCII − 33 (Sanger / standard Nanopore offset)
                            s = int(arr.sum()) - 33 * len(arr)
                            phred_sum_all_bases += s
                            phred_n_all_bases   += len(arr)
                            per_read_mean_phred_sum += s / len(arr)
                            per_read_count          += 1
                    i += 1
            except (EOFError, OSError):
                pass            # truncated file — keep partial stats
    except OSError:
        return None
    if not lengths:
        return None
    arr = np.array(lengths, dtype=np.int32)
    arr_desc = np.sort(arr)[::-1]
    return {
        "raw_reads":          int(arr.size),
        "raw_total_bases":    int(arr.sum()),
        "raw_p10_len":        int(np.percentile(arr, 10)),
        "raw_median_len":     int(np.median(arr)),
        "raw_p90_len":        int(np.percentile(arr, 90)),
        "raw_n50_len":        _n50(arr_desc),
        "raw_n_above_3kb":    int((arr >= AMPLICON_BP).sum()),
        "raw_pct_above_3kb":  float((arr >= AMPLICON_BP).mean() * 100),
        "raw_mean_phred_all_bases": (round(phred_sum_all_bases / phred_n_all_bases, 2)
                                     if phred_n_all_bases else None),
        "raw_mean_phred_per_read":  (round(per_read_mean_phred_sum / per_read_count, 2)
                                     if per_read_count else None),
    }

# ── Per-stage pipeline diagnostics — pinpoint where reads die ────────────────
LOGS_DIR = HERE / "logs"

def count_fasta_records(path: Path) -> int | None:
    """Count '>' header lines in a FASTA. Returns None if file missing."""
    if not path.exists():
        return None
    n = 0
    try:
        with open(path, "rb") as fh:
            for line in fh:
                if line.startswith(b">"):
                    n += 1
    except OSError:
        return None
    return n

def fasta_length_stats(path: Path) -> dict:
    """Return per-record length statistics from a FASTA. Useful for diagnosing
    which downstream filter (Step 4 length, 4b BLAST coverage, 7b N-content,
    Step 9 abundance) is killing contigs."""
    out = {
        "contig_min_len": None,
        "contig_median_len": None,
        "contig_max_len": None,
        "contig_n50": None,
        "contig_n_in_step4_range": None,
        "contig_n_below_step4_min": None,
        "contig_n_above_step4_max": None,
        "contig_pct_in_step4_range": None,
    }
    if not path.exists():
        return out
    lengths: list[int] = []
    cur = 0
    try:
        with open(path, "rb") as fh:
            for line in fh:
                if line.startswith(b">"):
                    if cur:
                        lengths.append(cur)
                    cur = 0
                else:
                    cur += len(line.rstrip(b"\n"))
            if cur:
                lengths.append(cur)
    except OSError:
        return out
    if not lengths:
        return out
    arr = np.array(lengths, dtype=np.int32)
    arr_desc = np.sort(arr)[::-1]
    n_in = int(((arr >= STEP4_MIN_BP) & (arr <= STEP4_MAX_BP)).sum())
    out.update({
        "contig_min_len":          int(arr.min()),
        "contig_median_len":       int(np.median(arr)),
        "contig_max_len":          int(arr.max()),
        "contig_n50":              _n50(arr_desc),
        "contig_n_in_step4_range": n_in,
        "contig_n_below_step4_min": int((arr < STEP4_MIN_BP).sum()),
        "contig_n_above_step4_max": int((arr > STEP4_MAX_BP).sum()),
        "contig_pct_in_step4_range": round(100.0 * n_in / len(arr), 1),
    })
    return out

def pipeline_diagnostics(bc_dir: Path, lib: str, bc: str) -> dict:
    """For one barcode, count contigs at each pipeline stage and pull the
    Step-2 chimera-filter summary. Missing files become None — downstream
    classification handles that gracefully."""
    out: dict = {}
    # Per-rep Canu contig counts
    rep_counts = []
    for rep in (1, 2, 3, 4):
        path = bc_dir / f"CANU_rep{rep}" / f"{lib}_{bc}_rep{rep}.contigs.fasta"
        c = count_fasta_records(path)
        out[f"canu_n_contigs_rep{rep}"] = c
        if c is not None:
            rep_counts.append(c)
    out["canu_n_contigs_total"] = sum(rep_counts) if rep_counts else None

    # Combined + deduplicated contig counts
    out["n_combined_contigs"] = count_fasta_records(bc_dir / "combined_contigs.fasta")
    out["n_unique_contigs"]   = count_fasta_records(bc_dir / "combined_unique_contigs.fasta")

    # Per-contig length stats (from combined_unique_contigs.fasta — the input
    # to the Step 2 chimera filter and the downstream Step 4 length filter)
    out.update(fasta_length_stats(bc_dir / "combined_unique_contigs.fasta"))

    # Step 2 chimera filter summary (logs/Library###_*_barcode##.chimera_stats.tsv)
    out["chim_n_keep"] = None
    out["chim_n_drop"] = None
    out["chim_median_uniformity"] = None
    out["chim_pct_dropped"] = None
    if LOGS_DIR.exists():
        chim_glob = list(LOGS_DIR.glob(f"{lib}_*_{bc}.chimera_stats.tsv"))
        if chim_glob:
            try:
                cs = pd.read_csv(chim_glob[0], sep="\t", encoding="utf-8-sig")
                if "decision" in cs.columns:
                    out["chim_n_keep"] = int((cs["decision"] == "KEEP").sum())
                    out["chim_n_drop"] = int((cs["decision"] == "DROP").sum())
                    n_total = len(cs)
                    if n_total:
                        out["chim_pct_dropped"] = round(
                            100.0 * out["chim_n_drop"] / n_total, 1)
                if "uniformity" in cs.columns and len(cs):
                    out["chim_median_uniformity"] = round(float(cs["uniformity"].median()), 3)
            except (OSError, pd.errors.EmptyDataError):
                pass
    return out

rows = []
print("Scanning Library*/barcode??/ — pulling FASTQ stats + pipeline-stage diagnostics...")
for lib_dir in sorted(HERE.glob("Library*")):
    if not lib_dir.is_dir():
        continue
    lib = lib_dir.name
    n_lib_samples = 0
    for bc_dir in sorted(lib_dir.glob("barcode??")):
        bc = bc_dir.name
        sample_id = f"{lib}_{bc}"
        report = bc_dir / "CANU_rep1" / f"{sample_id}_rep1.report"
        n_reads = reads_from_report(report)
        if n_reads is None:
            continue
        fq = bc_dir / f"{sample_id}.fastq.gz"
        fq_stats = fastq_stats(fq)
        diag = pipeline_diagnostics(bc_dir, lib, bc)
        row = {"Sample_ID": sample_id, "Library": lib, "Barcode": bc,
               "input_reads": n_reads}
        if fq_stats:
            row.update(fq_stats)
        row.update(diag)
        rows.append(row)
        n_lib_samples += 1
    print(f"  {lib}: {n_lib_samples} samples")

reads_df = pd.DataFrame(rows)
print(f"Pulled stats for {len(reads_df)} on-disk samples "
      f"across {reads_df['Library'].nunique()} libraries")
n_with_fastq = reads_df["raw_reads"].notna().sum() if "raw_reads" in reads_df.columns else 0
print(f"  ({n_with_fastq} have raw FASTQ stats; rest are Canu-report-only)")

# ── Join SI status + redo + BL/EO ────────────────────────────────────────────
si = pd.read_csv(IN_SI, sep="\t", encoding="utf-8-sig")
si = si[["Sample_ID", "EO_normalised", "BL_inferred",
         "n_haps_OK", "n_haps_REMOVED", "n_haps_chimeric", "n_haps_total",
         "SI_status"]].copy()

redo = pd.read_csv(IN_REDO, encoding="utf-8-sig")
redo = redo[["Sample_ID", "Lab_action", "Exclusion_stage"]].copy()

df = reads_df.merge(si, on="Sample_ID", how="left").merge(redo, on="Sample_ID", how="left")
df["Lab_action"] = df["Lab_action"].fillna("Functional / not flagged")
df["SI_status"] = df["SI_status"].fillna("Not_in_SI_table")
df["pass_si_threshold"] = (df["n_haps_total"].fillna(0) >= 4) & \
                          (df["SI_status"].isin(["SI", "pSI", "SC"]))

# ── Failure mode + recommended action ────────────────────────────────────────
# Decision tree (applied per sample). The ordering matters — the first matching
# rule wins. Each branch sets BOTH failure_mode and recommended_action.
LOW_RAW_READS = 10_000           # below this, the sample is genuinely under-sequenced
CHIM_DROP_PCT_AGGRESSIVE = 80    # chimera filter dropped >= 80 % of contigs → likely too aggressive

def classify_sample(r: pd.Series) -> tuple[str, str]:
    # 1. Samples that already produced a defensible SI call — no action.
    if r["pass_si_threshold"]:
        return ("none — SI call achieved",
                "Use directly (sample is functional in current dataset).")

    # 2. Samples not in the SI status table at all — outgroup / excluded individuals.
    if r["SI_status"] == "Not_in_SI_table":
        return ("not_in_SI_status_table",
                "Not part of the ingroup SI analysis — no action required.")

    # 3. DNA contamination (Step 9 too_many_alleles) — re-sequencing won't fix.
    if r["Lab_action"] == "Re-DNA-extraction":
        return ("DNA_contamination",
                "Re-extract DNA from fresh tissue, then re-PCR + re-sequence.")

    # 4. Re-PCR cases (Stage 4b paralogs, abundance failures, etc.) — flagged
    #    upstream of Phase 4; the SI failure is downstream of those problems.
    if r["Lab_action"] == "Re-PCR":
        return ("PCR_amplicon_failure",
                "Re-amplify the SRK amplicon from existing DNA; existing reads cannot be salvaged.")

    # 5. Chimera-dominated failure: lots of chimeric haplotypes, few clean.
    clean = (r["n_haps_OK"] if pd.notna(r["n_haps_OK"]) else 0) + \
            (r["n_haps_REMOVED"] if pd.notna(r["n_haps_REMOVED"]) else 0)
    chim  = r["n_haps_chimeric"] if pd.notna(r["n_haps_chimeric"]) else 0
    if chim >= 2 * max(clean, 1) and chim >= 4:
        return ("chimera_dominated",
                "Template-difficult amplicon — re-sequencing won't shift the chimera:clean ratio. "
                "Consider primer redesign, different polymerase, or longer-fragment library prep.")

    # 6. Partial recovery: between 1 and 3 clean haplotypes — close to the threshold.
    if 1 <= clean < 4:
        return ("partial_recovery",
                "Re-sequence on existing DNA + library at moderately higher coverage; "
                "a 2× read budget should push this sample over the 4-haplotype SI threshold.")

    # 7. From here on, clean == 0 — Canu output is empty for SI purposes.
    #    Decision splits by where the reads die.

    n_canu  = r.get("canu_n_contigs_total")
    n_uniq  = r.get("n_unique_contigs")
    chim_pct = r.get("chim_pct_dropped")

    # 7a. Step 2 chimera filter wiped everything.
    if pd.notna(chim_pct) and chim_pct >= CHIM_DROP_PCT_AGGRESSIVE \
            and pd.notna(n_uniq) and n_uniq > 0:
        return ("pipeline_chimera_filter_too_aggressive",
                f"Step 2 chimera filter dropped {chim_pct:.0f} % of contigs. "
                "Re-run Step 2 with a more permissive uniformity threshold OR investigate "
                "why the sample's coverage profile is so non-uniform.")

    # 7b. Canu produced contigs, seqkit rmdup or downstream filter wiped them.
    if pd.notna(n_canu) and n_canu > 0:
        if pd.notna(n_uniq) and n_uniq == 0:
            return ("pipeline_dedup_collapsed_all_contigs",
                    "Canu produced contigs but seqkit rmdup collapsed them to zero. "
                    "Inspect duplication / similarity of Canu reps.")
        # Now look at per-contig length. Step 4 keeps only contigs in
        # [STEP4_MIN_BP, STEP4_MAX_BP]. If unique contigs exist but none fall
        # in that range, Step 4 is the killer.
        n_in_step4 = r.get("contig_n_in_step4_range")
        n_below    = r.get("contig_n_below_step4_min")
        n_above    = r.get("contig_n_above_step4_max")
        if pd.notna(n_uniq) and n_uniq > 0 and pd.notna(n_in_step4) and n_in_step4 == 0:
            if pd.notna(n_below) and n_below > 0 and (pd.isna(n_above) or n_above == 0):
                return ("pipeline_step4_length_filter_too_short",
                        f"Canu produced {int(n_uniq)} unique contigs, all shorter than "
                        f"{STEP4_MIN_BP} bp (median = {int(r.get('contig_median_len', 0))} bp). "
                        "Step 4 length filter rejects every contig. Root cause is upstream of "
                        "Step 4: PCR / library prep produces fragmented amplicons; investigate "
                        "PCR conditions and library prep parameters.")
            if pd.notna(n_above) and n_above > 0 and (pd.isna(n_below) or n_below == 0):
                return ("pipeline_step4_length_filter_too_long",
                        f"Canu produced {int(n_uniq)} unique contigs, all longer than "
                        f"{STEP4_MAX_BP} bp (median = {int(r.get('contig_median_len', 0))} bp). "
                        "Step 4 length filter rejects every contig. Likely paralog amplification "
                        "with insertion (see Step 4b). Investigate primer specificity.")
            return ("pipeline_step4_length_filter_mixed",
                    f"Canu produced {int(n_uniq)} unique contigs; none fall in the Step 4 acceptance "
                    f"range ({STEP4_MIN_BP}–{STEP4_MAX_BP} bp). Mix of too-short and too-long contigs. "
                    "Investigate amplicon prep + paralog interference.")
        return ("pipeline_downstream_filter",
                "Canu + dedup + chimera filter retained contigs in the Step 4 length range, but "
                "the SI status pool is empty. Trace through Steps 4b / 7 / 7b / 9 to identify "
                "which filter dropped them.")

    # 7c. Canu itself failed.
    if pd.notna(n_canu) and n_canu == 0:
        if r.get("raw_reads") is not None and r["raw_reads"] < LOW_RAW_READS:
            return ("coverage_limited",
                    "Canu produced zero contigs and raw reads are low. "
                    "Re-sequence on existing DNA + library at higher coverage (target ≥ 25 000 raw reads).")
        # Adequate raw read count but Canu still empty — disambiguate with read
        # quality + length signals. The most actionable distinctions for the lab:
        #   * Low mean Phred → basecaller version / pore quality issue → re-basecall
        #     with a newer model before resequencing.
        #   * Short N50 → input DNA fragmented during library prep → fresh prep
        #     with attention to fragment-length control.
        #   * Very few reads ≥ 3 kb → reads exist but none span the amplicon →
        #     PCR / library-prep issue with amplicon recovery.
        mean_phred = r.get("raw_mean_phred_per_read")
        n50 = r.get("raw_n50_len")
        long_reads = r.get("raw_n_above_3kb")
        if pd.notna(mean_phred) and mean_phred < 10:
            return ("canu_assembly_failure_low_quality",
                    f"Canu produced zero contigs; raw read mean Phred = {mean_phred:.1f} (low). "
                    "Re-basecall with a newer Nanopore model before re-sequencing; the existing reads "
                    "may not be salvageable with current basecaller.")
        if pd.notna(n50) and n50 < 800:
            return ("canu_assembly_failure_fragmented",
                    f"Canu produced zero contigs; read N50 = {n50:,} bp (fragmented). "
                    "Reads are too short to span the SRK amplicon. Re-prep library with attention "
                    "to fragment-length retention before re-sequencing.")
        if pd.notna(long_reads) and long_reads < 50:
            return ("canu_assembly_failure_no_long_reads",
                    f"Canu produced zero contigs; only {int(long_reads)} reads ≥ 3 kb. "
                    "PCR amplicon recovery is poor — investigate PCR conditions before re-sequencing.")
        return ("canu_assembly_failure",
                "Canu produced zero contigs despite adequate raw reads, length, and quality. "
                "Investigate Canu parameters or sample-specific assembly behaviour before re-sequencing.")

    # 7d. No diagnostic data — fall back to crude raw-read rule.
    if r.get("raw_reads") is not None and r["raw_reads"] < LOW_RAW_READS:
        return ("coverage_limited",
                "Raw reads are below the under-sequenced threshold. "
                "Re-sequence at higher coverage on existing DNA + library.")
    return ("unresolved",
            "Diagnostic data missing — manual inspection required.")

df[["failure_mode", "recommended_action"]] = df.apply(
    lambda r: pd.Series(classify_sample(r)), axis=1)

df.to_csv(OUT_PER, sep="\t", index=False)
print(f"\nWrote {OUT_PER} ({len(df)} rows)")

print("\n── Failure mode distribution (all 500 samples) ──")
print(df["failure_mode"].value_counts().to_string())

print("\n── Failure mode × Lab_action ──")
print(pd.crosstab(df["failure_mode"], df["Lab_action"], margins=True).to_string())

# Subset focus: the 135 Re-sequence (SI status uncertain) cohort
reseq_subset = df[df["Lab_action"] == "Re-sequence (SI status uncertain)"]
if len(reseq_subset):
    print(f"\n── Re-sequence (SI status uncertain) cohort: failure mode × EO ──")
    print(pd.crosstab(reseq_subset["failure_mode"],
                      reseq_subset["EO_normalised"], margins=True).to_string())

# ── Dose-response: P(pass) vs input_reads, in deciles ────────────────────────
# Restrict to samples that AREN'T blocked by DNA contamination
#   (Re-DNA-extraction = more reads won't help; exclude from the curve)
curve_df = df[df["Lab_action"] != "Re-DNA-extraction"].copy()
curve_df = curve_df[curve_df["input_reads"] > 0].copy()

# Bin by reads (log-spaced)
edges = [0, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 1e9]
labels = ["<50", "50-100", "100-200", "200-400", "400-800", "800-1600",
          "1600-3200", "3200-6400", "6400-12800", "12800-25600", ">25600"]
curve_df["read_bin"] = pd.cut(curve_df["input_reads"], bins=edges, labels=labels,
                              include_lowest=True)

summ = curve_df.groupby("read_bin", observed=True).agg(
    n_samples=("Sample_ID", "size"),
    n_pass=("pass_si_threshold", "sum"),
    median_input_reads=("input_reads", "median"),
    median_n_haps_total=("n_haps_total", "median"),
    median_n_haps_chimeric=("n_haps_chimeric", "median"),
).reset_index()
summ["pct_pass"] = (summ["n_pass"] / summ["n_samples"] * 100).round(1)
summ.to_csv(OUT_SUMM, sep="\t", index=False)
print(f"\nDose-response (pass rate vs input reads):")
print(summ.to_string(index=False))

# ── Per-EO breakdown ─────────────────────────────────────────────────────────
print(f"\nPer-EO pass rate (samples with input_reads >= median):")
median_reads = curve_df["input_reads"].median()
hi = curve_df[curve_df["input_reads"] >= median_reads]
per_eo = hi.groupby("EO_normalised").agg(
    n=("Sample_ID", "size"),
    pass_rate=("pass_si_threshold", "mean"),
    median_raw_reads=("raw_reads", "median"),
    median_raw_len=("raw_median_len", "median"),
    median_pct_above_3kb=("raw_pct_above_3kb", "median"),
    median_chimeric_frac=("n_haps_chimeric", lambda s:
        (s / hi.loc[s.index, "n_haps_total"].clip(lower=1)).median()),
).round(3).sort_values("pass_rate")
print(per_eo.to_string())

# ── Read-length × chimera linkage ────────────────────────────────────────────
# The hypothesis from earlier: shorter reads → Canu has to bridge a longer span
# from fragmented evidence → more chimeric haplotypes. Test directly.
print(f"\n── Read-length vs chimera fraction (Spearman) ──")
has_fq = curve_df.dropna(subset=["raw_median_len", "n_haps_chimeric", "n_haps_total"]).copy()
has_fq["chim_frac"] = has_fq["n_haps_chimeric"] / has_fq["n_haps_total"].clip(lower=1)
if len(has_fq) >= 30:
    from scipy.stats import spearmanr
    rho_len, p_len = spearmanr(has_fq["raw_median_len"], has_fq["chim_frac"])
    rho_pct, p_pct = spearmanr(has_fq["raw_pct_above_3kb"], has_fq["chim_frac"])
    rho_rds, p_rds = spearmanr(has_fq["raw_reads"], has_fq["chim_frac"])
    print(f"  n = {len(has_fq)}")
    print(f"  rho(median read length,   chimera fraction) = {rho_len:+.3f}  p = {p_len:.3g}")
    print(f"  rho(% reads >= {AMPLICON_BP} bp, chimera fraction) = {rho_pct:+.3f}  p = {p_pct:.3g}")
    print(f"  rho(raw read count,       chimera fraction) = {rho_rds:+.3f}  p = {p_rds:.3g}")
else:
    print(f"  too few samples with raw FASTQ stats (n = {len(has_fq)}) to compute")

# ── Read-count target by quantile rule ───────────────────────────────────────
print(f"\n── Read-count target (90 % rule) ──")
# For each input_reads value, compute fraction passing among samples with
# >= that many reads. The smallest read-count where pass rate >= 90 % is the target.
sorted_reads = sorted(curve_df["input_reads"].unique())
targets = []
for r in sorted_reads:
    sub = curve_df[curve_df["input_reads"] >= r]
    if len(sub) < 10:
        continue
    rate = sub["pass_si_threshold"].mean()
    targets.append((r, len(sub), rate))
hits = [(r, n, rate) for r, n, rate in targets if rate >= 0.90]
if hits:
    target_reads, target_n, target_rate = hits[0]
    print(f"  smallest input_reads at which P(pass) >= 90 %: {target_reads:,} "
          f"(n = {target_n}, observed pass rate {target_rate:.1%})")
else:
    print(f"  no input_reads cutoff reaches 90 % pass rate in the on-disk data — "
          f"max observed pass rate = {max(r[2] for r in targets):.1%}")

hits_80 = [(r, n, rate) for r, n, rate in targets if rate >= 0.80]
if hits_80:
    t80, n80, rate80 = hits_80[0]
    print(f"  smallest input_reads at which P(pass) >= 80 %: {t80:,} "
          f"(n = {n80}, observed pass rate {rate80:.1%})")

# ── Failure mode at high coverage: probably-irrecoverable samples ────────────
hi_failures = curve_df[(curve_df["input_reads"] >= 5000) &
                       (~curve_df["pass_si_threshold"])]
print(f"\n── Probable irrecoverable samples (>= 5,000 reads, still failed) ──")
print(f"  n = {len(hi_failures)}")
if len(hi_failures) > 0:
    print(hi_failures[["Sample_ID", "EO_normalised", "BL_inferred", "input_reads",
                       "n_haps_OK", "n_haps_REMOVED", "n_haps_chimeric",
                       "SI_status", "Lab_action"]].to_string(index=False))

# ── Figure: 3 panels — read count, per-sample yield, read-length × chimera ──
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

# Panel A: bin-level pass rate
ax = axes[0]
xs = np.arange(len(summ))
ax.bar(xs, summ["pct_pass"], color="#377EB8", edgecolor="black", linewidth=0.4)
for i, (n, pct) in enumerate(zip(summ["n_samples"], summ["pct_pass"])):
    ax.text(i, pct + 1.5, f"n={n}", ha="center", fontsize=8)
ax.axhline(90, color="#E41A1C", linestyle="--", linewidth=1, label="90 % target")
ax.axhline(80, color="#FF7F00", linestyle="--", linewidth=1, label="80 % target")
ax.set_xticks(xs)
ax.set_xticklabels(summ["read_bin"], rotation=45, ha="right", fontsize=9)
ax.set_xlabel("Input reads per barcode (Canu, post seqkit-rmdup)")
ax.set_ylabel("% of samples reaching n_haps_total ≥ 4 (SI call)")
ax.set_title("A. Dose-response: read count vs SI-call success\n"
             "(Re-DNA-extraction samples excluded)",
             fontsize=10)
ax.set_ylim(0, 105)
ax.legend(loc="lower right", fontsize=8)
ax.grid(axis="y", alpha=0.25)

# Panel B: scatter input_reads vs n_haps_total, coloured by SI_status
ax = axes[1]
status_colour = {"SI": "#4DAF4A", "pSI": "#FF7F00", "SC": "#E41A1C",
                 "Insufficient_data": "#984EA3", "Not_in_SI_table": "grey"}
for status, color in status_colour.items():
    sub = curve_df[curve_df["SI_status"] == status]
    ax.scatter(sub["input_reads"], sub["n_haps_total"], s=14, alpha=0.55,
               color=color, edgecolor="none", label=f"{status} (n={len(sub)})")
ax.axhline(4, color="black", linestyle="--", linewidth=1, label="n_haps_total = 4 (SI call threshold)")
ax.set_xscale("log")
ax.set_xlabel("Input reads per barcode (log scale)")
ax.set_ylabel("Total clean haplotypes (OK + REMOVED_real)")
ax.set_title("B. Per-sample yield — what does extra coverage buy you?", fontsize=10)
ax.legend(loc="upper left", fontsize=8)
ax.grid(alpha=0.25)

# Panel C: read-length × chimera fraction. If reads are shorter than the
# amplicon, Canu has to bridge fragmented evidence → more chimeras.
ax = axes[2]
if len(has_fq) >= 30:
    for status, color in status_colour.items():
        sub = has_fq[has_fq["SI_status"] == status]
        if len(sub) == 0:
            continue
        ax.scatter(sub["raw_median_len"], sub["chim_frac"], s=14, alpha=0.55,
                   color=color, edgecolor="none", label=f"{status} (n={len(sub)})")
    ax.axvline(AMPLICON_BP, color="black", linestyle=":", linewidth=1,
               label=f"SRK amplicon length ≈ {AMPLICON_BP} bp")
    ax.set_xlabel("Raw FASTQ median read length (bp)")
    ax.set_ylabel("Chimeric-haplotype fraction (n_haps_chimeric / n_haps_total)")
    ax.set_title("C. Are short reads driving the chimeras?\n"
                 f"Spearman rho(median read length, chimera frac) = {rho_len:+.3f}  p = {p_len:.3g}",
                 fontsize=10)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.25)
else:
    ax.text(0.5, 0.5, f"Too few samples with FASTQ stats (n={len(has_fq)})",
            ha="center", va="center", transform=ax.transAxes)
    ax.set_axis_off()

fig.suptitle("Re-sequencing calibration — read count + read length vs SI-call success",
             fontsize=12, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(OUT_PDF, bbox_inches="tight")
fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
print(f"\nSaved {OUT_PDF}")
print(f"Saved {OUT_PNG}")
