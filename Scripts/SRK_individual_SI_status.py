#!/usr/bin/env python3
"""
SRK_individual_SI_status.py — Step 25

Reconstruct per-individual self-incompatibility (SI) status by re-aggregating
the per-haplotype OK / REMOVED calls from every
`all_Library*_frame1_stopcodon_log.tsv` produced by Step 7 (translate_filter_align_AA.py).

Rationale
─────────
Step 7 tags every Canu-assembled haplotype as `OK` (no internal stop) or
`REMOVED` (premature stop → non-functional SRK protein) before dropping the
REMOVED records ahead of the collapse + genotyping steps. The functional
information is preserved in the stop-codon logs, but it never reaches the
per-individual genotype calls — those count only functional copies.

For a tetraploid species with sporophytic SI at the SRK locus the biologically
meaningful question is: of the four expected SRK copies, how many produce a
functional protein? An individual carrying 4 non-functional copies is fully
self-compatible (SC); 4 functional copies is self-incompatible (SI); 1-3
non-functional copies is partially SI (pSI). This script reconstructs that
counter without invalidating the canonical 49-allele catalog or Phase 1-2
outputs.

Categorisation rule (per Sven, 2026-06-11)
──────────────────────────────────────────
    copies_nonfunctional = round(4 × n_REMOVED / (n_OK + n_REMOVED))
    copies_functional    = 4 - copies_nonfunctional

    SI_status:
        SI                    copies_nonfunctional == 0  OR  n_REMOVED < 2
                              (fully SI, or single stop call treated as
                              chimeric-assembly noise — see assembly caveat
                              below)
        pSI                   1 ≤ copies_nonfunctional ≤ 3 AND n_REMOVED ≥ 2
                              (partially SI — leaky)
        SC                    copies_nonfunctional == 4   (self-compatible
                              escape; robust because every haplotype must
                              independently carry a stop)
        Insufficient_data     n_OK + n_REMOVED < 4

    pSI_confidence (only set when SI_status == "pSI"):
        high                  frac_nonfunctional ≥ 0.25
        low                   frac_nonfunctional < 0.25

Assembly caveat
───────────────
Canu can produce N > 4 haplotypes per individual (chimeras, indel-rich
low-coverage contigs, allele-specific over-clustering). A single REMOVED
haplotype in an otherwise clean set is exactly what an assembly chimera
produces, so we treat n_REMOVED == 1 as noise and call those individuals SI.
The SC class is unchanged and robust because EVERY haplotype must
independently carry a stop (very hard to explain by chance).
The pSI class is reported with a confidence flag so downstream users can
filter on high-confidence pSI (n_REMOVED ≥ 2 AND frac_NF ≥ 0.25).

The `SRK_BEA` Brassica reference haplotypes are excluded before aggregation
(outgroup, not part of the LEPA SI question).

Inputs
──────
    all_Library*_frame1_stopcodon_log.tsv   one per library (Sequence_ID, Status)
    Tables/SRK_data_quality_categories.tsv  EO / BL / Ingroup metadata

Outputs
───────
    Tables/SRK_individual_SI_status.tsv     one row per ingroup individual
"""
from __future__ import annotations

import glob
import os
import sys
from pathlib import Path

import pandas as pd

# ─── Paths ────────────────────────────────────────────────────────────────────
HERE          = Path(__file__).resolve().parent
TABLES_DIR    = HERE / "Tables"
QC_TSV        = TABLES_DIR / "SRK_data_quality_categories.tsv"
OUT_TSV       = TABLES_DIR / "SRK_individual_SI_status.tsv"
STOPLOG_GLOB  = "all_Library*_frame1_stopcodon_log.tsv"
OUTGROUP_TAG  = "SRK_BEA"

# ─── Load metadata (EO / BL / Ingroup) ────────────────────────────────────────
if not QC_TSV.exists():
    sys.exit(f"ERROR: missing {QC_TSV} — run evaluate_data_quality.py first")

qc = pd.read_csv(QC_TSV, sep="\t", encoding="utf-8-sig", dtype=str).fillna("")
qc = qc[qc["Ingroup"] == "1"].copy()
meta = qc.set_index("Sample_ID")[["EO_normalised", "BL_inferred"]].to_dict("index")
print(f"Loaded metadata for {len(meta)} ingroup individuals from {QC_TSV.name}")

# ─── Aggregate stop-codon logs ────────────────────────────────────────────────
log_paths = sorted(glob.glob(str(HERE / STOPLOG_GLOB)))
if not log_paths:
    sys.exit(f"ERROR: no stop-codon logs found matching {STOPLOG_GLOB}")
print(f"Reading {len(log_paths)} stop-codon log(s)")

counts: dict[str, dict[str, int]] = {}

for path in log_paths:
    df = pd.read_csv(path, sep="\t", encoding="utf-8-sig")
    df = df[~df["Sequence_ID"].str.startswith(OUTGROUP_TAG)].copy()
    df["Individual"] = df["Sequence_ID"].str.split("_").str[:2].str.join("_")
    for ind, sub in df.groupby("Individual"):
        rec = counts.setdefault(ind, {"n_OK": 0, "n_REMOVED": 0})
        rec["n_OK"]      += int((sub["Status"] == "OK").sum())
        rec["n_REMOVED"] += int((sub["Status"] == "REMOVED").sum())

print(f"Found {len(counts)} individuals across all libraries (before metadata filter)")

# ─── Build per-individual table ───────────────────────────────────────────────
PSI_HIGH_FRAC = 0.25
NOISE_FILTER_MIN_REMOVED = 2  # n_REMOVED < this → call SI (assembly noise)


def classify(n_ok: int, n_rm: int):
    """Return (copies_functional, copies_nonfunctional, SI_status, pSI_confidence).

    copies_* are empty strings when n_OK + n_REMOVED < 4 (Insufficient_data).
    pSI_confidence is set only for SI_status == "pSI" (else empty string).
    """
    n_total = n_ok + n_rm
    if n_total < 4:
        return ("", "", "Insufficient_data", "")
    cp_nf = round(4 * n_rm / n_total)
    cp_f  = 4 - cp_nf
    frac  = n_rm / n_total
    if cp_nf == 4:
        return (cp_f, cp_nf, "SC", "")
    # Assembly-noise filter: single REMOVED hap → call SI, not pSI
    if cp_nf == 0 or n_rm < NOISE_FILTER_MIN_REMOVED:
        return (cp_f, cp_nf, "SI", "")
    # 1 ≤ cp_nf ≤ 3 AND n_REMOVED ≥ 2 → genuine pSI
    conf = "high" if frac >= PSI_HIGH_FRAC else "low"
    return (cp_f, cp_nf, "pSI", conf)


rows = []
for ind, m in meta.items():
    c = counts.get(ind, {"n_OK": 0, "n_REMOVED": 0})
    n_ok, n_rm = c["n_OK"], c["n_REMOVED"]
    cp_f, cp_nf, status, conf = classify(n_ok, n_rm)
    n_total = n_ok + n_rm
    rows.append({
        "Sample_ID":             ind,
        "EO_normalised":         m["EO_normalised"],
        "BL_inferred":           m["BL_inferred"],
        "n_haps_OK":             n_ok,
        "n_haps_REMOVED":        n_rm,
        "n_haps_total":          n_total,
        "frac_nonfunctional":    round(n_rm / n_total, 4) if n_total else "",
        "copies_functional":     cp_f,
        "copies_nonfunctional":  cp_nf,
        "SI_status":             status,
        "pSI_confidence":        conf,
    })

out = pd.DataFrame(rows).sort_values(["BL_inferred", "EO_normalised", "Sample_ID"])
TABLES_DIR.mkdir(exist_ok=True)
out.to_csv(OUT_TSV, sep="\t", index=False)

# ─── Console summary ──────────────────────────────────────────────────────────
print()
print(f"Wrote {len(out)} rows → {OUT_TSV}")
print()
print("SI_status distribution (species-level):")
print(out["SI_status"].value_counts().to_string())
print()
print("pSI confidence (high = frac_NF ≥ 0.25; low = frac_NF < 0.25):")
psi = out[out["SI_status"] == "pSI"]
print(psi["pSI_confidence"].value_counts().to_string() if len(psi) else "  (no pSI)")
print()
print("SI_status × BL:")
print(pd.crosstab(out["BL_inferred"], out["SI_status"]).to_string())
print()
print("SI_status × EO:")
print(pd.crosstab(out["EO_normalised"], out["SI_status"]).to_string())
