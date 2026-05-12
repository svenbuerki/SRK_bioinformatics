#!/usr/bin/env python3
"""
evaluate_data_quality.py — Step 12c — Data Quality Evaluation

Integrating script at the end of Phase 2 that asks three sequential questions:

  Q1.  Is there a library effect impacting data interpretation?
       → consumes SRK_library_effect_tests.tsv (from test_library_effect.py)

  Q2.  Which samples failed sequencing/assembly and need to be re-sequenced?
       → CSV deliverable for lab colleagues

  Q3.  Which samples have non-functional SRK proteins (full to partial) and
       may have escaped self-incompatibility?
       → CSV for controlled-selfing phenotyping

Visualises:
  - Stacked bar chart of outcome categories per bottleneck lineage (BL1–BL5)
  - Stacked bar chart of outcome categories per focus Element Occurrence (N≥5)

Outcome categories (mutually exclusive; every sample assigned to exactly one):

  Functional                     included; ≥50 % of SRK sequences translated
                                 cleanly — SI system intact at the molecular
                                 level.
  Partial_translation_failure    included; <50 % translation rate. NEUTRAL
                                 framing — does NOT imply biological SI loss.
                                 The most likely explanation is a technical
                                 artefact: Canu de novo assembly produces
                                 chimeric haplotypes whose spliced junctions
                                 introduce premature stops, which mechanically
                                 inflate the failure denominator without
                                 reflecting real biology. The current data
                                 cannot distinguish technical artefact from
                                 real partial SI degradation; we have NO
                                 tangible evidence supporting biological SI
                                 loss in this category. Treat as uncertain.
  SI_escape_candidate            excluded; Complete_loss + premature_stop
                                 dominated — strongest molecular signal of
                                 LoF mutations. (Defensible because every
                                 observed haplotype must carry a stop, which
                                 is hard to explain by chimeras alone.)
  Re_PCR                         excluded; existing DNA stock is fine, the PCR
                                 product is the problem. Lab action: re-amplify
                                 SRK from the existing DNA stock.
  Re_DNA_extraction              excluded; DNA stock itself contaminated
                                 (too_many_alleles). Lab action: re-isolate
                                 single-plant tissue and re-extract DNA before
                                 any further PCR.
  Outgroup                       intentional exclusion; not plotted

Inputs (auto-detected in CWD):
  SRK_sample_exclusion_audit.tsv      (from audit_sample_exclusions.py)
  SRK_individual_BL_assignments.tsv   (Step 13; provides EO→BL map)
  SRK_library_effect_tests.tsv        (from test_library_effect.py; optional)

Outputs:
  Tables/SRK_data_quality_categories.tsv   master per-sample category assignment
  Tables/SRK_samples_redo.csv              lab deliverable — merged Re-PCR +
                                           Re-DNA-extraction with a `Lab_action`
                                           column distinguishing the two
  Tables/SRK_SI_escape_candidates.csv      phenotyping deliverable
  SRK_data_quality_per_BL.{pdf,png}        per-BL stacked bar (figures/)
  SRK_data_quality_per_EO.{pdf,png}        per-EO stacked bar (figures/)
"""

from __future__ import annotations

import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


AUDIT_TSV   = "SRK_sample_exclusion_audit.tsv"
BL_TSV      = "SRK_individual_BL_assignments.tsv"
LIB_TSV     = "SRK_library_effect_tests.tsv"

TABLES_DIR    = "Tables"
os.makedirs(TABLES_DIR, exist_ok=True)

OUT_CAT_TSV   = "Tables/SRK_data_quality_categories.tsv"
OUT_REDO_CSV  = "Tables/SRK_samples_redo.csv"
OUT_SI_CSV    = "Tables/SRK_SI_escape_candidates.csv"

FIG_DIR     = "figures"
os.makedirs(FIG_DIR, exist_ok=True)

# BL palette (Set1, locked across the pipeline)
BL_COLORS = {
    "BL1": "#E41A1C",  # red
    "BL2": "#377EB8",  # blue
    "BL3": "#4DAF4A",  # green
    "BL4": "#984EA3",  # purple
    "BL5": "#FF7F00",  # orange
    "Unassigned": "#999999",
}

# Categories + plotting order (bottom to top of the stack)
CATEGORY_ORDER = [
    "Functional",
    "Partial_translation_failure",
    "SI_escape_candidate",
    "Re_PCR",
    "Re_DNA_extraction",
]
CAT_COLORS = {
    "Functional":                  "#2ca02c",  # green — SI intact at molecular level
    "Partial_translation_failure": "#ffbb78",  # light orange — uncertain (likely technical)
    "SI_escape_candidate":         "#d62728",  # red — biological SI escape signal
    "Re_PCR":                      "#1f77b4",  # blue — DNA fine, re-amplify
    "Re_DNA_extraction":           "#9467bd",  # purple — DNA stock contaminated
}

# Focus EOs from Phase 3 — those with N ≥ 5 in the final dataset (sorted by BL)
FOCUS_EOS_BY_BL = {
    "BL2": ["EO70"],
    "BL3": ["EO76"],
    "BL4": ["EO27", "EO67"],
    "BL5": ["EO18", "EO25"],
}
FOCUS_EOS = [eo for lst in FOCUS_EOS_BY_BL.values() for eo in lst]


# ─── Load inputs ──────────────────────────────────────────────────────────────
if not os.path.exists(AUDIT_TSV):
    sys.exit(f"ERROR: missing input — run audit_sample_exclusions.py first to produce {AUDIT_TSV}")

audit = pd.read_csv(AUDIT_TSV, sep="\t")

# Build EO → BL map from BL_assignments (uses the standardised EO labels)
bl_map = {}
if os.path.exists(BL_TSV):
    bl_df = pd.read_csv(BL_TSV, sep="\t")
    # Each EO maps to a single BL — use mode (most common) within EO
    bl_map = bl_df.groupby("EO")["BL"].agg(lambda x: x.mode().iloc[0] if len(x) > 0 else "").to_dict()

# The audit's EO column has raw metadata values (e.g. "27", "98-2"); normalise to "EO27"
def normalise_eo(raw: str) -> str:
    s = str(raw).strip()
    if not s or s.lower() == "nan":
        return ""
    if s.startswith("EO"):
        return s
    # Handle compound codes like "26-3" → "EO26"
    return f"EO{s.split('-')[0]}"

audit["EO_normalised"] = audit["EO"].apply(normalise_eo)
audit["BL_inferred"]   = audit["EO_normalised"].map(bl_map).fillna("Unassigned")


# ─── Categorise every sample ──────────────────────────────────────────────────
def categorise(row: pd.Series) -> str:
    final  = row["Final_status"]
    stage  = str(row["Exclusion_stage"] or "")
    si     = row["SI_functional_status"]
    mode   = row["Dominant_failure_mode"]
    ingroup = str(row["Ingroup"])

    if ingroup == "0":
        return "Outgroup"

    if final == "Included":
        if si == "Functional":
            return "Functional"
        if si == "Partial_translation_failure":
            return "Partial_translation_failure"
        return "Functional"  # safety fallback

    # Excluded samples — categorise by lab action required
    if si == "Complete_loss" and mode == "premature_stop":
        return "SI_escape_candidate"

    # Re-DNA-extraction: DNA stock itself contaminated. Re-PCR won't help.
    if stage == "Step 9 — too_many_alleles":
        return "Re_DNA_extraction"

    # Re-PCR: DNA stock fine, PCR product is the problem. Re-amplify.
    # This covers: no assembly, fragmented product, paralog amplification
    # (could shift with stochastic PCR competition), N-rich short product,
    # low yield (singleton functional proteins), and ambiguous-AA/mixed
    # translation failures that reflect short / dirty PCR product.
    if stage.startswith("Stage 3"):
        return "Re_PCR"
    if stage.startswith("Step 4 — length"):
        return "Re_PCR"
    if stage.startswith("Step 4b"):
        return "Re_PCR"
    if stage.startswith("Step 7b"):
        return "Re_PCR"
    if stage == "Step 9 — dropped_abundance_filter":
        return "Re_PCR"
    if stage == "Step 9 — no_functional_proteins" and mode in ("ambiguous_aa", "mixed"):
        return "Re_PCR"

    # Should not reach here for ingroup samples — flag for investigation
    return "Re_PCR"


audit["Outcome_category"] = audit.apply(categorise, axis=1)


# ─── Write master per-sample categorisation TSV ───────────────────────────────
master_cols = [
    "Sample_ID", "Library", "Barcode", "EO", "EO_normalised", "BL_inferred",
    "Ingroup",
    "n_raw_haps", "n_after_step4", "n_after_step4b", "n_after_step7b",
    "Step9_classification", "Proteins_in_final_data",
    "SI_functional_status", "Dominant_failure_mode",
    "BL_status", "Final_status", "Exclusion_stage", "Exclusion_reason",
    "Outcome_category",
]
audit[master_cols].to_csv(OUT_CAT_TSV, sep="\t", index=False)


# ─── Q2: Merged Samples-to-Redo CSV (Re-PCR + Re-DNA-extraction combined) ─────
# Lab colleagues get ONE master file with a discriminator column.
re_pcr_df = audit[audit["Outcome_category"] == "Re_PCR"].copy()
re_dna_df = audit[audit["Outcome_category"] == "Re_DNA_extraction"].copy()

def re_pcr_action(row: pd.Series) -> str:
    stage = row["Exclusion_stage"]
    if "Stage 3" in stage:
        return "Re-PCR (no Canu assembly produced — try fresh primers / longer extension)"
    if "Step 4 — length" in stage:
        return "Re-PCR (current product fragmented — increase extension time)"
    if "Step 4b" in stage:
        return ("Re-PCR (current amplicon is SRK paralog rather than canonical SRK; "
                "PCR competition may resolve stochastically — if same paralog signature "
                "persists across replicates, consider primer re-design)")
    if "Step 7b" in stage:
        return "Re-PCR (short PCR product yields N-rich Canu contigs — aim for full-length amplicon)"
    if "dropped_abundance_filter" in stage:
        return "Re-PCR (low yield — singleton functional proteins below abundance threshold; re-amplify for more product)"
    if "no_functional_proteins" in stage:
        return "Re-PCR (short / dirty product translates to X residues — re-amplify, possibly with proofreading polymerase)"
    return "Re-PCR (see audit reason)"

re_pcr_df["Lab_action"]         = "Re-PCR"
re_pcr_df["Recommended_action"] = re_pcr_df.apply(re_pcr_action, axis=1)

re_dna_df["Lab_action"]         = "Re-DNA-extraction"
re_dna_df["Recommended_action"] = (
    "Re-isolate tissue from a single plant + re-extract DNA. "
    ">4 distinct functional proteins per individual indicates mixed/contaminated "
    "DNA stock; re-PCR will preserve the contamination. Verify single-plant origin "
    "of the source tissue before re-extracting."
)

# Merge and order: Re-DNA-extraction first (longer pipeline, higher operational
# priority for scheduling), then Re-PCR; within each block, sort by EO for
# field-collection grouping.
redo = pd.concat([re_dna_df, re_pcr_df], ignore_index=True)
redo_cols = [
    "Lab_action", "Sample_ID", "Library", "Barcode", "EO", "EO_normalised",
    "BL_inferred", "n_raw_haps", "Proteins_in_final_data", "Exclusion_stage",
    "Recommended_action",
]
redo = redo[redo_cols].sort_values(["Lab_action", "EO_normalised", "Sample_ID"])
redo.to_csv(OUT_REDO_CSV, index=False)


# ─── Q3: SI-escape candidates CSV ─────────────────────────────────────────────
si_escape = audit[audit["Outcome_category"] == "SI_escape_candidate"].copy()
si_escape["Recommended_action"] = (
    "Phenotype via controlled selfing test (no-pollen-deposition vs self-pollen "
    "vs cross-pollen seed-set comparison) to validate self-compatibility"
)
si_escape_cols = [
    "Sample_ID", "Library", "Barcode", "EO", "EO_normalised", "BL_inferred",
    "n_raw_haps", "Step9_classification", "Dominant_failure_mode",
    "Exclusion_reason", "Recommended_action",
]
si_escape[si_escape_cols].to_csv(OUT_SI_CSV, index=False)


# ─── Build per-BL and per-EO category matrices ────────────────────────────────
plot_df = audit[audit["Outcome_category"] != "Outgroup"].copy()

# Per-BL counts (BL1..BL5; samples with Unassigned BL pooled separately)
bl_order = ["BL1", "BL2", "BL3", "BL4", "BL5"]
bl_data = (
    plot_df[plot_df["BL_inferred"].isin(bl_order)]
    .groupby(["BL_inferred", "Outcome_category"])
    .size()
    .unstack(fill_value=0)
    .reindex(index=bl_order)
    .reindex(columns=CATEGORY_ORDER, fill_value=0)
)
bl_props = bl_data.div(bl_data.sum(axis=1), axis=0)
bl_totals = bl_data.sum(axis=1)

# Per-EO counts (focus EOs only; sorted by parent BL)
eo_order = FOCUS_EOS
eo_data = (
    plot_df[plot_df["EO_normalised"].isin(eo_order)]
    .groupby(["EO_normalised", "Outcome_category"])
    .size()
    .unstack(fill_value=0)
    .reindex(index=eo_order, fill_value=0)
    .reindex(columns=CATEGORY_ORDER, fill_value=0)
)
eo_props = eo_data.div(eo_data.sum(axis=1), axis=0)
eo_totals = eo_data.sum(axis=1)


# ─── Stacked-bar plotting helper ──────────────────────────────────────────────
def stacked_bar(props: pd.DataFrame, totals: pd.Series, title: str,
                xlabel: str, png_path: str, pdf_path: str,
                category_counts: pd.DataFrame, bl_label_map: dict | None = None):
    fig, ax = plt.subplots(figsize=(10, 6))
    bottom = np.zeros(len(props))
    x = np.arange(len(props))
    for cat in CATEGORY_ORDER:
        if cat not in props.columns:
            continue
        vals = props[cat].values
        ax.bar(x, vals, bottom=bottom, color=CAT_COLORS[cat], edgecolor="white",
               linewidth=0.5, label=cat.replace("_", " "))
        # Annotate non-trivial slices with absolute counts
        for i, (p, n) in enumerate(zip(vals, category_counts[cat].values)):
            if p >= 0.04:
                ax.text(x[i], bottom[i] + p / 2, f"{int(n)}",
                        ha="center", va="center",
                        color="white" if cat in ("SI_escape_candidate", "Excluded_technical", "Re-sequence")
                              else "black",
                        fontsize=9, fontweight="bold")
        bottom += vals

    # x-tick labels with N
    labels = [f"{idx}\n(N = {totals.loc[idx]})" for idx in props.index]
    if bl_label_map:
        # Recolour x-tick labels by BL
        for i, idx in enumerate(props.index):
            bl_label = bl_label_map.get(idx, "")
            if bl_label in BL_COLORS:
                pass  # we'll set after setting xticks
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    if bl_label_map:
        for i, lbl in enumerate(ax.get_xticklabels()):
            bl = bl_label_map.get(props.index[i], "")
            if bl in BL_COLORS:
                lbl.set_color(BL_COLORS[bl])
                lbl.set_fontweight("bold")
    elif "BL" in xlabel:
        # Colour BL labels directly
        for lbl, bl in zip(ax.get_xticklabels(), props.index):
            if bl in BL_COLORS:
                lbl.set_color(BL_COLORS[bl])
                lbl.set_fontweight("bold")

    ax.set_ylabel("Proportion of metadata samples")
    ax.set_xlabel(xlabel)
    ax.set_ylim(0, 1.02)
    ax.set_title(title)
    handles = [Patch(facecolor=CAT_COLORS[c], label=c.replace("_", " ")) for c in CATEGORY_ORDER]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    fig.tight_layout()
    fig.savefig(png_path, dpi=200)
    fig.savefig(pdf_path)
    plt.close(fig)


# BL plot
stacked_bar(
    bl_props, bl_totals,
    title="Step 12c — Sample outcome categories per Bottleneck Lineage",
    xlabel="Bottleneck Lineage",
    png_path=os.path.join(FIG_DIR, "SRK_data_quality_per_BL.png"),
    pdf_path="SRK_data_quality_per_BL.pdf",
    category_counts=bl_data,
)

# EO plot — colour x-tick labels by parent BL
eo_to_bl_label = {}
for bl, eo_list in FOCUS_EOS_BY_BL.items():
    for eo in eo_list:
        eo_to_bl_label[eo] = bl

stacked_bar(
    eo_props, eo_totals,
    title="Step 12c — Sample outcome categories per focus Element Occurrence",
    xlabel="Focus EO (label colour = parent BL)",
    png_path=os.path.join(FIG_DIR, "SRK_data_quality_per_EO.png"),
    pdf_path="SRK_data_quality_per_EO.pdf",
    category_counts=eo_data,
    bl_label_map=eo_to_bl_label,
)


# ─── Stdout summary ───────────────────────────────────────────────────────────
print()
print("=" * 72)
print("Step 12c — Data Quality Evaluation")
print("=" * 72)

# Q1: library effect
print()
print("Q1 — Library effect on data interpretation?")
if os.path.exists(LIB_TSV):
    lib_tests = pd.read_csv(LIB_TSV, sep="\t")
    n_sig = lib_tests["Significant"].sum()
    print(f"  Library effect tests loaded ({len(lib_tests)} tests). Significant: {n_sig}")
    sig_rows = lib_tests[lib_tests["Significant"]]
    for _, r in sig_rows.iterrows():
        print(f"    {r['Test']}: {r['Comparison']}  (p = {r['p_value']:.3g})")
    if n_sig == 0:
        print("  → No significant library effect detected.")
else:
    print(f"  ({LIB_TSV} not found — run test_library_effect.py first to test for library effects.)")

# Q2: Merged redo CSV (Re-PCR + Re-DNA-extraction in one file)
print()
print(f"Q2 — Samples needing lab redo: {len(redo)} → {OUT_REDO_CSV}")
print("  Lab_action breakdown:")
for action, n in redo["Lab_action"].value_counts().items():
    print(f"    {n:>3}  {action}")
re_pcr_stages = redo[redo["Lab_action"]=="Re-PCR"]["Exclusion_stage"].value_counts()
if len(re_pcr_stages) > 0:
    print("  Re-PCR by stage:")
    for stage, n in re_pcr_stages.items():
        print(f"    {n:>3}  {stage}")
re_dna_eos = redo[redo["Lab_action"]=="Re-DNA-extraction"]["EO_normalised"].value_counts()
if len(re_dna_eos) > 0:
    print("  Re-DNA-extraction by EO:")
    for eo, n in re_dna_eos.items():
        bl = bl_map.get(eo, "?")
        print(f"    {n:>3}  {eo or '(no EO)'}  ({bl})")

# Q3: SI-escape candidates
print()
print(f"Q3 — Candidate SI-escape samples (Complete_loss + premature_stop): {len(si_escape)} → {OUT_SI_CSV}")
if len(si_escape) > 0:
    print("  EO breakdown:")
    for eo, n in si_escape["EO_normalised"].value_counts().items():
        bl = bl_map.get(eo, "?")
        print(f"    {n:>3}  {eo or '(no EO)'}  ({bl})")

# Outcome category totals
print()
print("Outcome category totals (excluding outgroups):")
overall = audit[audit["Outcome_category"] != "Outgroup"]["Outcome_category"].value_counts().reindex(CATEGORY_ORDER, fill_value=0)
for cat, n in overall.items():
    pct = 100 * n / overall.sum() if overall.sum() else 0
    print(f"  {n:>4}  ({pct:4.1f} %)  {cat}")

print()
print(f"Outputs:")
print(f"  {OUT_CAT_TSV}     — master per-sample categorisation")
print(f"  {OUT_REDO_CSV}                — Re-PCR + Re-DNA-extraction (lab deliverable)")
print(f"  {OUT_SI_CSV}      — SI-escape candidates (phenotyping deliverable)")
print(f"  {FIG_DIR}/SRK_data_quality_per_BL.png   — per-BL stacked bar")
print(f"  {FIG_DIR}/SRK_data_quality_per_EO.png   — per-EO stacked bar")
print(f"  SRK_data_quality_per_BL.pdf, SRK_data_quality_per_EO.pdf  — vector versions")
