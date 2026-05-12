#!/usr/bin/env python3
"""
audit_sample_exclusions.py
==========================

Per-sample audit of the SRK pipeline — for every barcode present in
`sampling_metadata.csv`, traces the sample through each filter stage and
reports (i) whether the sample is in the final genotyped dataset, and
(ii) if not, the EARLIEST pipeline stage at which it was lost and the
specific reason.

Conservation use: identifies which collected samples need to be re-sequenced
(stage 3 / step 4 / step 4b failures) versus which are biologically informative
even though excluded from the SI analyses (step 11 outgroup, step 13
BL-unassigned germplasm sub-codes).

Cascade (earliest exclusion wins):
    Stage 3 — no assembly             (sample in metadata, no raw haplotypes)
    Step 4 — length filter             (raw haps exist, 0 survive 3250 ≤ L ≤ 4000)
    Step 4b — BLAST coverage           (Library 010+; all haps dropped at 0.90 / 95 %)
    Step 7b — N-content                (Library 010+; all haps dropped at max_N = 9)
    Step 9 — translation/abundance     (no_functional_proteins / dropped_abundance_filter /
                                        too_many_alleles classifications)
    Step 11 — outgroup filter          (Ingroup = 0)
    NOT excluded                       (in final zygosity TSV)

BL_Unassigned is reported as an auxiliary flag (still Included; just not BL-stratified).

Two additional biological columns alongside the QC cascade:

    SI_functional_status:
        Functional                       Functional_rate ≥ 0.5 — SI system intact
                                         at the molecular level
        Partial_translation_failure      0 < rate < 0.5 (Step 9 low_functional_rate
                                         classification). DOES NOT IMPLY BIOLOGICAL
                                         SI LOSS. The most likely explanation is a
                                         technical artefact: Canu de novo assembly
                                         produces chimeric haplotypes whose spliced
                                         junctions introduce premature stops, which
                                         inflate the failure denominator without
                                         reflecting real biology. Step 9's abundance
                                         filter (min_count = 5) excludes chimeric
                                         proteins from the numerator but not from
                                         the denominator. The current dataset
                                         cannot distinguish technical artefact from
                                         real partial SI degradation; treat as
                                         uncertain.
        Complete_loss                    Total_sequences > 0 but 0 functional
                                         proteins — candidate SI breakdown /
                                         self-compatibility escape signal (when
                                         combined with premature_stop dominant
                                         failure mode; see Dominant_failure_mode).
        NA                               sample did not reach Step 9 (excluded
                                         earlier in pipeline).

    Dominant_failure_mode (parsed from Step 9 Top_failure_reasons):
        premature_stop    failures dominated by internal stop codons → real LoF
                          mutations, strong SI breakdown signal
        ambiguous_aa      failures dominated by X residues → likely N-rich data
                          quality issue (not biological)
        mixed             both modes contribute comparably
        other             different failure pattern
        (blank)           no failure data (no Step 9 entry)

Conservation interpretation: Complete_loss + premature_stop-dominated samples
are the strongest candidates for ongoing self-compatibility evolution and
deserve follow-up phenotyping (controlled selfing tests).

Inputs (auto-detected in CWD):
    sampling_metadata.csv
    all_LibraryXXX_Phased_haplotypes.fasta                (raw Canu output)
    all_LibraryXXX_*_min3250*.fasta                       (Step 4 output)
    all_LibraryXXX_*_blastfilt_log.tsv                    (Step 4b decisions; optional)
    all_LibraryXXX_*_backfilled_Nfilt_log.tsv             (Step 7b decisions; optional)
    SRK_individual_status_report.tsv                      (Step 9 classifications)
    SRK_individual_BL_assignments.tsv                     (Step 13 BL status)
    SRK_individual_zygosity.tsv                           (final included list)

Output:
    SRK_sample_exclusion_audit.tsv                        (one row per metadata sample)
"""

from __future__ import annotations

import csv
import glob
import os
import re
import sys
from collections import defaultdict


METADATA          = "sampling_metadata.csv"
STATUS_REPORT     = "SRK_individual_status_report.tsv"
BL_ASSIGNMENTS    = "SRK_individual_BL_assignments.tsv"
ZYGOSITY          = "SRK_individual_zygosity.tsv"
OUT_TSV           = "SRK_sample_exclusion_audit.tsv"


# ─── Helpers ──────────────────────────────────────────────────────────────────
def parse_individual_from_haplotype_id(hap_id: str) -> str | None:
    """Extract `LibraryXXX_barcodeYY` from a haplotype/protein FASTA header."""
    m = re.match(r"(Library\d+_barcode\d+)", hap_id)
    return m.group(1) if m else None


def count_haps_per_individual(fasta_path: str) -> dict[str, int]:
    """Count FASTA records per `LibraryXXX_barcodeYY` individual."""
    counts: dict[str, int] = defaultdict(int)
    if not os.path.exists(fasta_path):
        return counts
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                ind = parse_individual_from_haplotype_id(line[1:].strip())
                if ind:
                    counts[ind] += 1
    return counts


def count_keeps_per_individual(tsv_path: str, decision_col: str = "decision",
                               query_id_col: str = "query_id") -> dict[str, int]:
    """Count rows with decision == 'KEEP' per individual in a per-haplotype log TSV."""
    counts: dict[str, int] = defaultdict(int)
    if not os.path.exists(tsv_path):
        return counts
    with open(tsv_path, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get(decision_col, "").strip() == "KEEP":
                ind = parse_individual_from_haplotype_id(row.get(query_id_col, ""))
                if ind:
                    counts[ind] += 1
    return counts


# ─── Load metadata ────────────────────────────────────────────────────────────
samples: dict[str, dict] = {}
with open(METADATA, encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    for row in reader:
        sid = row.get("SampleID", "").strip()
        if not sid:
            continue
        samples[sid] = {
            "Library":       row.get("Library", "").strip(),
            "Barcode":       row.get("barcode", "").strip(),
            "EO":            row.get("EO_w_sub", row.get("Pop", "")).strip(),
            "Ingroup":       row.get("Ingroup", "1").strip(),
        }

print(f"Metadata: {len(samples)} samples loaded")


# ─── Load final zygosity (= INCLUDED individuals) ─────────────────────────────
included: set[str] = set()
if os.path.exists(ZYGOSITY):
    with open(ZYGOSITY, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ind = (row.get("Individual") or row.get("SampleID") or "").strip()
            if ind:
                included.add(ind)
print(f"Final zygosity (Included): {len(included)} individuals")


# ─── Load Step 9 status report ────────────────────────────────────────────────
step9: dict[str, dict] = {}
if os.path.exists(STATUS_REPORT):
    with open(STATUS_REPORT, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ind = row.get("Individual", "").strip()
            if ind:
                step9[ind] = {
                    "Total_sequences":      row.get("Total_sequences", ""),
                    "Functional_proteins":  row.get("Functional_proteins", ""),
                    "Functional_rate":      row.get("Functional_rate", ""),
                    "Proteins_in_final":    row.get("Proteins_in_final_data", ""),
                    "Classification":       row.get("Classification", ""),
                    "Top_failure":          row.get("Top_failure_reasons", ""),
                }
print(f"Step 9 status report: {len(step9)} individuals")


# ─── Load BL assignments ──────────────────────────────────────────────────────
bl: dict[str, str] = {}
if os.path.exists(BL_ASSIGNMENTS):
    with open(BL_ASSIGNMENTS, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ind = row.get("Individual", "").strip()
            if ind:
                bl[ind] = row.get("BL_status", "").strip()
print(f"BL assignments: {len(bl)} individuals "
      f"({sum(1 for v in bl.values() if v == 'Unassigned')} Unassigned)")


# ─── Aggregate per-library haplotype counts ───────────────────────────────────
raw_haps:    dict[str, int] = defaultdict(int)
step4_haps:  dict[str, int] = defaultdict(int)
step4b_haps: dict[str, int] = defaultdict(int)
step7b_haps: dict[str, int] = defaultdict(int)

for fasta in sorted(glob.glob("all_Library*_Phased_haplotypes.fasta")):
    for ind, n in count_haps_per_individual(fasta).items():
        raw_haps[ind] += n

# Step 4 outputs: prefer the max4000 variant when both exist for the same library
seen_libraries_step4: set[str] = set()
for fasta in sorted(glob.glob("all_Library*_Phased_haplotypes_filtered_min3250_max4000.fasta")):
    lib = re.search(r"Library\d+", fasta).group(0)
    seen_libraries_step4.add(lib)
    for ind, n in count_haps_per_individual(fasta).items():
        step4_haps[ind] += n
# Then handle libraries that only have the legacy _min3250.fasta (no max4000)
for fasta in sorted(glob.glob("all_Library*_Phased_haplotypes_filtered_min3250.fasta")):
    lib = re.search(r"Library\d+", fasta).group(0)
    if lib in seen_libraries_step4:
        continue
    for ind, n in count_haps_per_individual(fasta).items():
        step4_haps[ind] += n

for tsv in sorted(glob.glob("all_Library*_blastfilt_log.tsv")):
    for ind, n in count_keeps_per_individual(tsv).items():
        step4b_haps[ind] += n

for tsv in sorted(glob.glob("all_Library*_backfilled_Nfilt_log.tsv")):
    for ind, n in count_keeps_per_individual(tsv).items():
        step7b_haps[ind] += n


# ─── Decide outcome per sample ────────────────────────────────────────────────
# Step 9 classifications that ARE exclusions (i.e. drop the individual from the
# final genotype data). low_functional_rate and normal are not exclusions.
STEP9_EXCLUSIONS = {
    "no_functional_proteins": "all sequences failed translation / validation",
    "dropped_abundance_filter": "functional proteins present but below min_count = 5 abundance filter",
    "too_many_alleles": ">4 distinct proteins per individual — likely contamination / assembly artefact",
    "no_sequences": "no DNA sequences reached define_functional_proteins.py",
}


_FAIL_RE = re.compile(r"([a-z_]+)\((\d+)\)")


def parse_failure_modes(top_failure: str) -> dict[str, int]:
    """Parse 'premature_stop(28);ambiguous_aa(20)' → {'premature_stop': 28, 'ambiguous_aa': 20}."""
    return {m.group(1): int(m.group(2)) for m in _FAIL_RE.finditer(top_failure or "")}


def dominant_failure(top_failure: str) -> str:
    """Return the dominant failure mode label given Step 9 Top_failure_reasons."""
    modes = parse_failure_modes(top_failure)
    if not modes:
        return ""
    sorted_modes = sorted(modes.items(), key=lambda kv: -kv[1])
    top_mode, top_n = sorted_modes[0]
    if len(sorted_modes) == 1:
        return top_mode
    second_mode, second_n = sorted_modes[1]
    # "mixed" if the second-most-common mode is at least 50 % of the top one
    if second_n >= 0.5 * top_n:
        return "mixed"
    return top_mode


def si_functional_status(s9_row: dict) -> str:
    """Categorise SI functional status from a Step 9 status report row.

    Returns one of:
        Functional                   rate ≥ 0.5
        Partial_translation_failure  0 < rate < 0.5 — NEUTRAL framing. Does not
                                     imply biological SI loss; most likely
                                     reflects chimeric Canu assemblies inflating
                                     the failure denominator. See module
                                     docstring for full caveat.
        Complete_loss                rate = 0 with non-zero sequence count
        NA                           sample did not reach Step 9
    """
    if not s9_row:
        return "NA"
    try:
        n_total = int(s9_row.get("Total_sequences") or 0)
        n_func  = int(s9_row.get("Functional_proteins") or 0)
        rate    = float(s9_row.get("Functional_rate") or 0)
    except ValueError:
        return "NA"
    if n_total == 0:
        return "NA"
    if n_func == 0:
        return "Complete_loss"
    if rate < 0.5:
        return "Partial_translation_failure"
    return "Functional"


def decide(ind: str) -> tuple[str, str, str]:
    """Return (final_status, exclusion_stage, exclusion_reason)."""
    has_raw      = raw_haps.get(ind, 0)
    has_step4    = step4_haps.get(ind, 0)
    has_step4b   = step4b_haps.get(ind, 0)
    has_step7b   = step7b_haps.get(ind, 0)
    s9           = step9.get(ind, {})
    classification = s9.get("Classification", "")
    ingroup      = samples[ind]["Ingroup"]
    in_zyg       = ind in included
    lib_has_4b   = bool(glob.glob(f"all_{ind.split('_')[0]}_*_blastfilt_log.tsv"))
    lib_has_7b   = bool(glob.glob(f"all_{ind.split('_')[0]}_*_backfilled_Nfilt_log.tsv"))

    if in_zyg:
        return ("Included", "", "")

    # Earliest exclusion wins
    if has_raw == 0:
        return ("Excluded", "Stage 3 — no assembly",
                "sample present in metadata but no haplotypes assembled by Canu — needs re-sequencing")
    if has_step4 == 0:
        return ("Excluded", "Step 4 — length filter",
                f"{has_raw} raw haps, 0 passed 3250–4000 bp length filter — likely fragmented or chimeric assembly; re-sequencing recommended")
    if lib_has_4b and has_step4b == 0:
        return ("Excluded", "Step 4b — BLAST coverage filter",
                f"{has_step4} haps after Step 4, 0 passed coverage ≥ 0.90 / identity ≥ 95 vs canonical SRK — likely SRK paralogs / off-target; do NOT re-sequence (issue is sequence identity, not coverage)")
    if lib_has_7b and has_step7b == 0:
        return ("Excluded", "Step 7b — N-content filter",
                f"{has_step4b} haps after Step 4b, 0 passed N ≤ 9 — short Canu contigs filled with N's by Step 7 backfill; re-sequencing may help if more reads recover longer contigs")
    if classification in STEP9_EXCLUSIONS:
        base_reason = STEP9_EXCLUSIONS[classification]
        # Enrich no_functional_proteins reason with SI-breakdown interpretation
        if classification == "no_functional_proteins":
            dom = dominant_failure(s9.get("Top_failure", ""))
            if dom == "premature_stop":
                base_reason = ("all sequences carry premature stop codons — strong candidate "
                               "for ongoing SI breakdown / self-compatibility escape; "
                               "recommend controlled selfing test to phenotype")
            elif dom == "ambiguous_aa":
                base_reason = ("all sequences carry X residues from N-rich data — quality "
                               "issue, NOT a biological SI escape; re-sequencing may help "
                               "(would have been caught by Step 7b if applied)")
            elif dom == "mixed":
                base_reason = ("mixed premature_stop + ambiguous_aa failures — both LoF "
                               "mutations AND data-quality issues; interpret cautiously")
        return ("Excluded", f"Step 9 — {classification}",
                f"{base_reason} (Top failures: {s9.get('Top_failure', 'n/a')})")
    if ingroup == "0":
        return ("Excluded", "Step 11 — outgroup filter",
                "Ingroup = 0 in metadata; intentionally excluded from SI analyses (NOT a quality issue)")
    return ("Excluded", "Unknown",
            "no exclusion stage identified — investigate manually")


# ─── Write audit TSV ──────────────────────────────────────────────────────────
header = [
    "Sample_ID", "Library", "Barcode", "EO", "Ingroup",
    "n_raw_haps", "n_after_step4", "n_after_step4b", "n_after_step7b",
    "Step9_classification", "Proteins_in_final_data",
    "SI_functional_status", "Dominant_failure_mode",
    "BL_status",
    "Final_status", "Exclusion_stage", "Exclusion_reason",
]

rows_out = []
for sid in sorted(samples):
    meta = samples[sid]
    s9   = step9.get(sid, {})
    final_status, ex_stage, ex_reason = decide(sid)
    rows_out.append({
        "Sample_ID":              sid,
        "Library":                meta["Library"],
        "Barcode":                meta["Barcode"],
        "EO":                     meta["EO"],
        "Ingroup":                meta["Ingroup"],
        "n_raw_haps":             raw_haps.get(sid, 0),
        "n_after_step4":          step4_haps.get(sid, 0),
        "n_after_step4b":         step4b_haps.get(sid, ""),
        "n_after_step7b":         step7b_haps.get(sid, ""),
        "Step9_classification":   s9.get("Classification", ""),
        "Proteins_in_final_data": s9.get("Proteins_in_final", ""),
        "SI_functional_status":   si_functional_status(s9),
        "Dominant_failure_mode":  dominant_failure(s9.get("Top_failure", "")),
        "BL_status":              bl.get(sid, ""),
        "Final_status":           final_status,
        "Exclusion_stage":        ex_stage,
        "Exclusion_reason":       ex_reason,
    })

with open(OUT_TSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows_out)


# ─── Summary on stdout ────────────────────────────────────────────────────────
n_total       = len(rows_out)
n_included    = sum(1 for r in rows_out if r["Final_status"] == "Included")
n_excluded    = n_total - n_included

stage_counts: dict[str, int] = defaultdict(int)
for r in rows_out:
    if r["Final_status"] == "Excluded":
        stage_counts[r["Exclusion_stage"]] += 1

print()
print(f"Audit written to: {OUT_TSV}")
print(f"  Total samples in metadata: {n_total}")
print(f"  Included in final zygosity: {n_included}")
print(f"  Excluded: {n_excluded}")
print()
print("Excluded by stage (earliest exclusion wins):")
for stage in sorted(stage_counts, key=lambda s: -stage_counts[s]):
    print(f"  {stage_counts[stage]:>4}  {stage}")

# ─── SI functional status breakdown ───────────────────────────────────────────
print()
print("SI functional status (samples that reached Step 9):")
si_counts: dict[str, int] = defaultdict(int)
for r in rows_out:
    s = r["SI_functional_status"]
    if s and s != "NA":
        si_counts[s] += 1
for status in ["Functional", "Partial_loss", "Complete_loss"]:
    print(f"  {si_counts.get(status, 0):>4}  {status}")

# ─── Highlight candidate SI-escape samples ────────────────────────────────────
escape_candidates = [
    r for r in rows_out
    if r["SI_functional_status"] == "Complete_loss"
    and r["Dominant_failure_mode"] == "premature_stop"
]
quality_complete_loss = [
    r for r in rows_out
    if r["SI_functional_status"] == "Complete_loss"
    and r["Dominant_failure_mode"] == "ambiguous_aa"
]
mixed_complete_loss = [
    r for r in rows_out
    if r["SI_functional_status"] == "Complete_loss"
    and r["Dominant_failure_mode"] == "mixed"
]

if escape_candidates:
    print()
    print(f"Candidate SI-escape samples ({len(escape_candidates)}): "
          f"sequencing successful but ALL proteins carry premature stop codons")
    print(f"  → Strong loss-of-function signal; recommend phenotyping via controlled selfing")
    for r in escape_candidates:
        print(f"    {r['Sample_ID']}  EO={r['EO'] or '(none)'}  "
              f"N_seqs={r['n_raw_haps']}  failures={r['Exclusion_reason'].split('Top failures: ')[-1].rstrip(')')}")

if quality_complete_loss:
    print()
    print(f"Complete_loss but ambiguous_aa-dominated ({len(quality_complete_loss)}): "
          f"likely data quality issue, NOT biological SI escape")
    for r in quality_complete_loss:
        print(f"    {r['Sample_ID']}  EO={r['EO'] or '(none)'}  N_seqs={r['n_raw_haps']}")

if mixed_complete_loss:
    print()
    print(f"Complete_loss with mixed failure modes ({len(mixed_complete_loss)}): "
          f"interpret cautiously — both LoF and quality issues contribute")
    for r in mixed_complete_loss:
        print(f"    {r['Sample_ID']}  EO={r['EO'] or '(none)'}  N_seqs={r['n_raw_haps']}")

partial_loss = [r for r in rows_out if r["SI_functional_status"] == "Partial_loss"]
partial_lof = [r for r in partial_loss if r["Dominant_failure_mode"] == "premature_stop"]
if partial_lof:
    print()
    print(f"Partial_loss (included in dataset) with premature_stop-dominated failures: "
          f"{len(partial_lof)} samples")
    print(f"  → These individuals have functional SRK but >50 % of their SRK sequences carry")
    print(f"    premature stops — a softer signal of ongoing SI degradation worth tracking.")
