#!/usr/bin/env python3
"""
SRK_genotype_null_alleles.py — Step 26

Build null-aware tetraploid genotype tables from the canonical functional-only
genotype tables plus the per-individual SI status from Step 25.

Rationale
─────────
The canonical 335-individual dataset (`SRK_individual_allele_genotypes.tsv` +
`SRK_individual_zygosity.tsv`) counts only FUNCTIONAL SRK copies and pads
under-recovered individuals to four slots by homozygosity assumption. This
inflates apparent homozygosity for partial-SI (pSI) individuals and excludes
SC individuals entirely. Step 26 produces a parallel set of tables in which:

  • pSI individuals carry explicit `Allele_NULL` copies (= copies_nonfunctional
    from Step 25). Functional slots are scaled down proportionally so the row
    still sums to 4.
  • SC individuals are added with `Allele_NULL == 4` (genotype 0000).
  • SI individuals are kept as-is (no nulls introduced).
  • pSI_low individuals are *promoted to SI* — their non-functional fraction
    is too dilute (frac_NF < 0.25) to confidently assign nulls; introducing
    them would propagate Canu-noise uncertainty into the downstream allele
    frequencies. Sven, 2026-06-11.
  • Insufficient_data individuals are excluded from the null-aware canonical
    set AND written to a re-sequencing redo list.

The original Phase-1–2 tables remain unchanged and frozen as the
functional-only reference; downstream scripts may load either.

Inputs
──────
  Tables/SRK_individual_allele_genotypes.tsv   335 × 49 functional-only matrix
  Tables/SRK_individual_zygosity.tsv           335-row zygosity + Genotype
  Tables/SRK_individual_SI_status.tsv          401-row SI status from Step 25

Outputs
───────
  Tables/SRK_individual_allele_genotypes_with_nulls.tsv
      One row per individual in the null-aware canonical set; 49 functional
      Allele_xxx columns + Allele_NULL + Genotype_class_flag. Rows sum to 4.
  Tables/SRK_individual_zygosity_with_nulls.tsv
      Per-individual N_distinct, N_total (= 4 by construction), Genotype
      (e.g. AABC, AB00, 0000), Allele_composition, SI_status, pSI_confidence,
      Genotype_class_flag.
  Tables/SRK_samples_for_redo.tsv
      Insufficient_data individuals with EO / BL / hap counts + recommended
      action. These are the samples to prioritise for re-sequencing.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE       = Path(__file__).resolve().parent
TAB        = HERE / "Tables"
GENO_TSV   = TAB / "Phase2" / "step11_individual_allele_genotypes.tsv"
ZYG_TSV    = TAB / "Phase2" / "step12_individual_zygosity.tsv"
SI_TSV     = TAB / "Phase4" / "step22b_individual_SI_status.tsv"
OUT_GENO   = TAB / "Phase4" / "step23_individual_allele_genotypes_with_nulls.tsv"
OUT_ZYG    = TAB / "Phase4" / "step23_individual_zygosity_with_nulls.tsv"
OUT_REDO   = TAB / "Phase4" / "step23_samples_for_redo.tsv"


def scale_to_target(counts: dict[str, int], target: int) -> dict[str, int]:
    """Scale integer copy counts proportionally to sum to `target`.

    Used to allocate `copies_functional` slots when the original genotype
    row (which sums to 4) needs to fit into fewer slots after nulls are
    introduced. Uses largest-remainder rounding to keep totals exact.
    """
    total = sum(counts.values())
    if total == 0 or target == 0:
        return {k: 0 for k in counts}
    if total == target:
        return dict(counts)
    # Real-valued targets
    scaled = {k: (v / total) * target for k, v in counts.items()}
    floored = {k: int(v) for k, v in scaled.items()}
    deficit = target - sum(floored.values())
    if deficit > 0:
        # Distribute remainder to keys with largest fractional part
        remainders = sorted(scaled.keys(), key=lambda k: scaled[k] - floored[k],
                            reverse=True)
        for k in remainders:
            if deficit == 0:
                break
            floored[k] += 1
            deficit -= 1
    elif deficit < 0:
        # Trim from keys with smallest fractional part (and >0 count)
        remainders = sorted([k for k in scaled if floored[k] > 0],
                            key=lambda k: scaled[k] - floored[k])
        for k in remainders:
            if deficit == 0:
                break
            floored[k] -= 1
            deficit += 1
    return floored


def encode_genotype(counts: dict[str, int]) -> tuple[str, int]:
    """Convert a count dict (sums to 4) into a 4-character genotype label.

    Functional alleles get letters A, B, C in descending count order;
    Allele_NULL maps to 0. Returns (genotype_label, n_distinct_functional).
    """
    # Sort functional alleles by count desc
    func = [(a, c) for a, c in counts.items() if a != "Allele_NULL" and c > 0]
    func.sort(key=lambda kv: (-kv[1], kv[0]))
    letters = ["A", "B", "C", "D"]
    label_parts = []
    for i, (_, c) in enumerate(func):
        label_parts.append(letters[i] * c)
    null_n = counts.get("Allele_NULL", 0)
    label_parts.append("0" * null_n)
    label = "".join(label_parts)
    n_distinct = len(func)
    return label, n_distinct


def allele_composition_str(counts: dict[str, int]) -> str:
    """Reconstruct the Allele_composition string from a count dict."""
    parts = []
    # functional alleles by descending count, then by allele name
    func = [(a, c) for a, c in counts.items() if a != "Allele_NULL" and c > 0]
    func.sort(key=lambda kv: (-kv[1], kv[0]))
    for a, c in func:
        parts.append(f"{a}({c})")
    null_n = counts.get("Allele_NULL", 0)
    if null_n > 0:
        parts.append(f"Allele_NULL({null_n})")
    return "+".join(parts)


# ─── Load inputs ──────────────────────────────────────────────────────────────
geno = pd.read_csv(GENO_TSV, sep="\t", encoding="utf-8-sig")
zyg  = pd.read_csv(ZYG_TSV, sep="\t", encoding="utf-8-sig")
si   = pd.read_csv(SI_TSV, sep="\t", encoding="utf-8-sig")

allele_cols = [c for c in geno.columns if c.startswith("Allele_")]
print(f"Loaded {len(geno)} individuals × {len(allele_cols)} alleles from {GENO_TSV.name}")
print(f"Loaded {len(si)} individuals from {SI_TSV.name}")

si_idx = si.set_index("Sample_ID")[
    ["EO_normalised", "BL_inferred", "SI_status", "pSI_confidence",
     "copies_functional", "copies_nonfunctional",
     "n_haps_OK", "n_haps_REMOVED", "n_haps_total"]
]

# ─── Process every individual in the canonical 335 ────────────────────────────
geno_rows = []
zyg_rows  = []

for _, r in geno.iterrows():
    ind = r["Individual"]
    counts = {a: int(r[a]) for a in allele_cols}
    total = sum(counts.values())  # should be 4

    if ind not in si_idx.index:
        si_status = "Unknown"
        psi_conf  = ""
        cp_nf     = 0
        cp_f      = 4
        eo        = ""
        bl        = ""
    else:
        rec = si_idx.loc[ind]
        si_status = rec["SI_status"]
        psi_conf  = rec["pSI_confidence"] if pd.notna(rec["pSI_confidence"]) else ""
        cp_nf     = int(rec["copies_nonfunctional"]) if str(rec["copies_nonfunctional"]).strip() not in {"", "nan"} else 0
        cp_f      = int(rec["copies_functional"])    if str(rec["copies_functional"]).strip()    not in {"", "nan"} else 4
        eo        = rec["EO_normalised"]
        bl        = rec["BL_inferred"]

    # Decide treatment. Note: canonical `counts` sums to N_total_proteins,
    # NOT 4 — the homozygosity extension to four slots is implicit in the
    # Genotype field of SRK_individual_zygosity.tsv. We make padding explicit
    # by scaling to `target_func` and topping up with `target_null` NULLs.
    if si_status == "SI":
        target_func, target_null = 4, 0
        flag = "SI_canonical"
    elif si_status == "pSI" and psi_conf == "low":
        target_func, target_null = 4, 0
        flag = "SI_promoted_from_pSI_low"
    elif si_status == "pSI" and psi_conf == "high" and cp_nf > 0:
        target_func, target_null = cp_f, cp_nf
        flag = "pSI_null_introduced"
    elif si_status == "SC":
        target_func, target_null = 0, 4
        flag = "SC_all_null_from_genotype_table"
    elif si_status == "Insufficient_data":
        # It survived genotyping despite raw_haps<4 — the call is data-thin
        # but not nonsense. Keep existing functional alleles, no nulls; flag
        # for redo. Pad to 4 functional copies the same way SI does.
        target_func, target_null = 4, 0
        flag = "Insufficient_kept_existing"
    else:
        target_func, target_null = 4, 0
        flag = "No_SI_status"

    if target_func == 0:
        new_counts = {a: 0 for a in allele_cols}
    else:
        new_counts = scale_to_target(counts, target_func)
    new_counts["Allele_NULL"] = target_null

    label, n_distinct = encode_genotype(new_counts)
    geno_rows.append({"Individual": ind, **new_counts,
                      "Genotype_class_flag": flag,
                      "EO_normalised": eo, "BL_inferred": bl,
                      "SI_status": si_status, "pSI_confidence": psi_conf})
    zyg_rows.append({"Individual": ind,
                     "N_distinct_functional_alleles": n_distinct,
                     "N_total_copies": sum(new_counts.values()),
                     "N_null_copies": new_counts["Allele_NULL"],
                     "Genotype": label,
                     "Allele_composition": allele_composition_str(new_counts),
                     "SI_status": si_status,
                     "pSI_confidence": psi_conf,
                     "Genotype_class_flag": flag,
                     "EO_normalised": eo, "BL_inferred": bl})

# ─── Add SC individuals not yet in the canonical set ─────────────────────────
existing = set(geno["Individual"])
sc_extras = si[(si["SI_status"] == "SC") & (~si["Sample_ID"].isin(existing))]
for _, r in sc_extras.iterrows():
    ind = r["Sample_ID"]
    new_counts = {a: 0 for a in allele_cols}; new_counts["Allele_NULL"] = 4
    geno_rows.append({"Individual": ind, **new_counts,
                      "Genotype_class_flag": "SC_added_from_SI_status",
                      "EO_normalised": r["EO_normalised"], "BL_inferred": r["BL_inferred"],
                      "SI_status": "SC", "pSI_confidence": ""})
    zyg_rows.append({"Individual": ind,
                     "N_distinct_functional_alleles": 0,
                     "N_total_copies": 4, "N_null_copies": 4,
                     "Genotype": "0000",
                     "Allele_composition": "Allele_NULL(4)",
                     "SI_status": "SC", "pSI_confidence": "",
                     "Genotype_class_flag": "SC_added_from_SI_status",
                     "EO_normalised": r["EO_normalised"], "BL_inferred": r["BL_inferred"]})

# ─── Write outputs ───────────────────────────────────────────────────────────
out_geno_cols = ["Individual"] + allele_cols + ["Allele_NULL",
                 "Genotype_class_flag", "EO_normalised", "BL_inferred",
                 "SI_status", "pSI_confidence"]
out_zyg_cols  = ["Individual", "N_distinct_functional_alleles", "N_total_copies",
                 "N_null_copies", "Genotype", "Allele_composition",
                 "SI_status", "pSI_confidence", "Genotype_class_flag",
                 "EO_normalised", "BL_inferred"]

pd.DataFrame(geno_rows).reindex(columns=out_geno_cols).to_csv(
    OUT_GENO, sep="\t", index=False)
pd.DataFrame(zyg_rows).reindex(columns=out_zyg_cols).to_csv(
    OUT_ZYG, sep="\t", index=False)

# ─── Redo list ────────────────────────────────────────────────────────────────
redo = si[si["SI_status"] == "Insufficient_data"][
    ["Sample_ID", "EO_normalised", "BL_inferred",
     "n_haps_OK", "n_haps_REMOVED", "n_haps_total"]
].copy()
redo["redo_priority"] = redo["n_haps_total"].apply(
    lambda n: "high" if n == 0 else "medium" if n < 2 else "low"
)
redo["redo_action"] = "Re-sequencing recommended (raw_haps < 4 → cannot resolve tetraploid SI status)"
redo.to_csv(OUT_REDO, sep="\t", index=False)

# ─── Summary ──────────────────────────────────────────────────────────────────
df_g = pd.read_csv(OUT_GENO, sep="\t")
print()
print(f"Wrote {len(df_g)} rows → {OUT_GENO.name}")
print(f"Wrote {len(df_g)} rows → {OUT_ZYG.name}")
print(f"Wrote {len(redo)} rows → {OUT_REDO.name}")
print()
print("Genotype_class_flag distribution:")
print(df_g["Genotype_class_flag"].value_counts().to_string())
print()
print("Redo priority distribution:")
print(redo["redo_priority"].value_counts().to_string())
print()
print("Sanity: every row should sum to 4 across all 50 allele columns.")
sums = df_g[allele_cols + ["Allele_NULL"]].sum(axis=1)
print(f"  Row sums — min: {sums.min()}, max: {sums.max()}, "
      f"all == 4: {bool((sums == 4).all())}")
