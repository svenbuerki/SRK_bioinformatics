#!/usr/bin/env python3
"""
SRK_individual_SI_status.py — Step 25b (per-individual SI categorisation)

Reconstruct per-individual self-incompatibility (SI) status from the
per-haplotype OK / REMOVED calls produced by Step 7
(`translate_filter_align_AA.py`) and the per-haplotype broken-allele
identity assignments produced by Step 25a (`SRK_null_allele_assignment.py`).

Method
──────
Step 7 tags every Canu-assembled haplotype as `OK` (functional protein) or
`REMOVED` (premature stop). Canu also produces many chimeric / indel-rich
extra contigs per individual that fail the stop filter without representing
real broken SRK copies — these are quantitatively distinguishable from
genuine broken alleles by their high AA-distance from every functional
allele in the catalogue (Step 25a). Step 25b therefore:

  1. Counts OK haplotypes per individual (n_haps_OK).
  2. Counts REMOVED haplotypes that Step 25a assigned to a functional allele
     with `high` or `medium` confidence (n_haps_REMOVED) — these are real
     broken alleles.
  3. Counts REMOVED haplotypes flagged as chimeric / unassignable in
     Step 25a (n_haps_chimeric) — recorded for data-quality transparency
     but excluded from the tetraploid-copy signal pool.
  4. Tetraploid copy estimate uses only the clean signal pool:

        n_total              = n_haps_OK + n_haps_REMOVED
        copies_nonfunctional = round(4 × n_haps_REMOVED / n_total)
        copies_functional    = 4 − copies_nonfunctional

  5. SI_status classification:
        SI                copies_nonfunctional == 0
                          OR n_haps_REMOVED < 2  (single real broken hap
                          treated as residual noise — see floor below)
        pSI               1 ≤ copies_nonfunctional ≤ 3 AND n_haps_REMOVED ≥ 2
                          pSI_confidence: high if frac_NF ≥ 0.25, low otherwise
        SC                copies_nonfunctional == 4
        Insufficient_data n_total < 4

The `SRK_BEA` Brassica reference haplotypes are excluded before aggregation
(outgroup).

Inputs
──────
    all_Library*_frame1_stopcodon_log.tsv   one per library (Step 7)
    Tables/SRK_null_allele_assignments.tsv  per-haplotype identity + confidence (Step 25a)
    Tables/SRK_data_quality_categories.tsv  EO / BL / Ingroup metadata (Step 12c)

Outputs
───────
    Tables/SRK_individual_SI_status.tsv     one row per ingroup individual
"""
from __future__ import annotations

import glob
import sys
from pathlib import Path

import pandas as pd

# ─── Paths ────────────────────────────────────────────────────────────────────
HERE          = Path(__file__).resolve().parent
TABLES_DIR    = HERE / "Tables"
QC_TSV        = TABLES_DIR / "Phase2" / "step12c_data_quality_categories.tsv"
ASSIGN_TSV    = TABLES_DIR / "Phase5" / "step25a_null_allele_assignments.tsv"
OUT_TSV       = TABLES_DIR / "Phase5" / "step25b_individual_SI_status.tsv"
STOPLOG_GLOB  = "all_Library*_frame1_stopcodon_log.tsv"
OUTGROUP_TAG  = "SRK_BEA"

# Classification thresholds
PSI_HIGH_FRAC            = 0.25
NOISE_FILTER_MIN_REMOVED = 2   # n_real_broken < this → call SI

# ─── Load metadata (EO / BL / Ingroup) ────────────────────────────────────────
if not QC_TSV.exists():
    sys.exit(f"ERROR: missing {QC_TSV} — run evaluate_data_quality.py first")
if not ASSIGN_TSV.exists():
    sys.exit(f"ERROR: missing {ASSIGN_TSV} — run "
             "SRK_null_allele_assignment.py (Step 25a) first")

qc = pd.read_csv(QC_TSV, sep="\t", encoding="utf-8-sig", dtype=str).fillna("")
qc = qc[qc["Ingroup"] == "1"].copy()
meta = qc.set_index("Sample_ID")[["EO_normalised", "BL_inferred"]].to_dict("index")
print(f"Loaded metadata for {len(meta)} ingroup individuals from {QC_TSV.name}")

# ─── Load Step 25a assignments — split REMOVED haps into real vs chimeric ────
assign = pd.read_csv(ASSIGN_TSV, sep="\t", encoding="utf-8-sig")
real_set = set(assign.loc[
    assign["confidence"].isin(["high", "medium"]), "Sequence_ID"
])
chimeric_set = set(assign.loc[
    assign["confidence"] == "low", "Sequence_ID"
])
print(f"Step 25a assignments: {len(real_set)} real broken + "
      f"{len(chimeric_set)} chimeric")

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
        rec = counts.setdefault(ind, {"n_OK": 0, "n_REMOVED": 0, "n_chimeric": 0})
        rec["n_OK"] += int((sub["Status"] == "OK").sum())
        rm_sub = sub[sub["Status"] == "REMOVED"]
        n_real = int(rm_sub["Sequence_ID"].isin(real_set).sum())
        n_chim = int(rm_sub["Sequence_ID"].isin(chimeric_set).sum())
        # REMOVED haps not in either set (unlikely — should be 0) are not counted
        rec["n_REMOVED"]  += n_real
        rec["n_chimeric"] += n_chim

print(f"Found {len(counts)} individuals across all libraries (pre-metadata)")


# ─── Classifier ───────────────────────────────────────────────────────────────
def classify(n_ok: int, n_rm: int):
    """Return (copies_functional, copies_nonfunctional, SI_status, pSI_confidence)."""
    n_total = n_ok + n_rm
    if n_total < 4:
        return ("", "", "Insufficient_data", "")
    cp_nf = round(4 * n_rm / n_total)
    cp_f  = 4 - cp_nf
    frac  = n_rm / n_total
    if cp_nf == 4:
        return (cp_f, cp_nf, "SC", "")
    if cp_nf == 0 or n_rm < NOISE_FILTER_MIN_REMOVED:
        return (cp_f, cp_nf, "SI", "")
    conf = "high" if frac >= PSI_HIGH_FRAC else "low"
    return (cp_f, cp_nf, "pSI", conf)


# ─── Build per-individual table ───────────────────────────────────────────────
rows = []
for ind, m in meta.items():
    c = counts.get(ind, {"n_OK": 0, "n_REMOVED": 0, "n_chimeric": 0})
    n_ok = c["n_OK"]
    n_rm = c["n_REMOVED"]
    n_ch = c["n_chimeric"]
    cp_f, cp_nf, status, conf = classify(n_ok, n_rm)
    n_total = n_ok + n_rm
    rows.append({
        "Sample_ID":             ind,
        "EO_normalised":         m["EO_normalised"],
        "BL_inferred":           m["BL_inferred"],
        "n_haps_OK":             n_ok,
        "n_haps_REMOVED":        n_rm,
        "n_haps_chimeric":       n_ch,
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
