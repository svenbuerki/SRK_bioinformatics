#!/usr/bin/env python3
"""Re-annotate the HV-only Step 26e cross plan under the FUNCTIONAL definition.

For each cross in the existing Step 26e output (H0/H1a/H1b/H2/H3 files), look up
the mother's and father's test alleles in the functional synonymy groups and
functional-site distance table, then determine what cross category the pair
would fall into under the functional definition (Incompatible / Synonymy_test /
Compatible_within / Compatible_cross). Crosses whose category changes between
the two definitions are flagged as HYPOTHESIS-DISCRIMINATING — these are the
laboratory crosses that most efficiently distinguish which allele-identity
definition is empirically correct.

Inputs:
    Tables/Phase5/step26b_functional_synonymy_groups.csv     (from srk_functional_definitions.py)
    Tables/Phase5/step26b_functional_allele_distances.tsv    (from srk_functional_definitions.py)
    Tables/Phase5/step26e_cross_plan_H*.tsv                  (HV-only, from srk_cross_plan.py)

Outputs:
    Tables/Phase5/step26e_functional_cross_plan_annotated.tsv
    Tables/Phase5/step26e_hypothesis_discriminating_crosses.tsv
    Tables/Phase5/step26e_functional_cross_plan_transitions.tsv

The threshold for the Synonymy_test vs Compatible_within boundary defaults to
0.04 — the same value the HV-only plan uses. In the current dataset this may
need recalibration once functional distances have been benchmarked, but 0.04 is
the conservative starting point that matches the existing pipeline convention.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

BASE = Path("Tables/Phase5")
IN_FG          = BASE / "step26b_functional_synonymy_groups.csv"
IN_DIST        = BASE / "step26b_functional_allele_distances.tsv"
PLAN_FILES = {
    "H0":  BASE / "step26e_cross_plan_H0_SI_validation.tsv",
    "H1a": BASE / "step26e_cross_plan_H1a_within_class_baseline.tsv",
    "H1b": BASE / "step26e_cross_plan_H1b_between_class_baseline.tsv",
    "H2":  BASE / "step26e_cross_plan_H2_synonymy_tests.tsv",
    "H3":  BASE / "step26e_cross_plan_H3_hidden_bin_tests.tsv",
}
OUT_ANNOT     = BASE / "step26e_functional_cross_plan_annotated.tsv"
OUT_DISCRIM   = BASE / "step26e_hypothesis_discriminating_crosses.tsv"
OUT_TRANS     = BASE / "step26e_functional_cross_plan_transitions.tsv"

WITHIN_CLASS_THRESHOLD = 0.04         # matches existing pipeline default


def parse_test_allele(cell: str, hv_hypothesis: str, role: str) -> str:
    """Extract the primary test allele from a cross-plan cell.
    'Allele_046'                    → 'Allele_046'
    'Allele_050(3)+Allele_055(1)'   → 'Allele_055'   (heterozygous → low-count = test/hidden)
    """
    if not isinstance(cell, str):
        return ""
    cell = cell.strip()
    if "+" in cell:
        parts = re.split(r"\+", cell)
        alleles = []
        for p in parts:
            m = re.match(r"(Allele_\d+)\((\d+)\)", p.strip())
            if m:
                alleles.append((m.group(1), int(m.group(2))))
        if alleles:
            # For heterozygous cells, test allele = lowest-copy (the hidden/rare one)
            alleles.sort(key=lambda x: x[1])
            return alleles[0][0]
        return cell
    return cell


def reclassify_category(m_fg: str, f_fg: str,
                        f_dist: float | None,
                        hv_category: str) -> str:
    """Assign a cross category under the functional definition.
    Class information (Compatible_cross) is inherited from the HV plan (Class I/II
    split is deep and doesn't shift under the functional refinement)."""
    if hv_category == "Compatible_cross":
        return "Compatible_cross"       # different Class — unchanged
    if m_fg == f_fg:                    # same functional group
        return "Incompatible"
    if f_dist is None or pd.isna(f_dist):
        return "unknown"
    if f_dist < WITHIN_CLASS_THRESHOLD:
        return "Synonymy_test"
    return "Compatible_within"


def infer_new_hypothesis(hv_h: str, hv_category: str, functional_category: str) -> str:
    """Map (HV hypothesis, functional category) → the hypothesis level the cross
    would sit at under the functional plan. For heterozygous-parent crosses
    (H1b, H3) we keep the H-level (the genotype constraint doesn't change with
    the definition); only the category label may shift."""
    if hv_h in ("H1b", "H3"):
        return hv_h                     # heterozygous crosses keep their H-level
    # H0 / H1a / H2 all use AAAA × AAAA — reassigned purely by category
    if functional_category == "Incompatible":
        return "H0"
    if functional_category == "Synonymy_test":
        return "H2"
    if functional_category == "Compatible_within":
        return "H1a"
    return "unknown"


def main() -> None:
    # --- Load lookups ---
    fg = pd.read_csv(IN_FG)
    allele_to_fg = dict(zip(fg["Allele"], fg["Functional_group"]))

    dist_df = pd.read_csv(IN_DIST, sep="\t")
    dist_lookup: dict[tuple[str, str], float] = {}
    for _, r in dist_df.iterrows():
        a, b, d = r["Allele_a"], r["Allele_b"], float(r["p_distance"])
        dist_lookup[(a, b)] = d
        dist_lookup[(b, a)] = d

    # --- Iterate the HV-only cross plan and re-annotate each cross ---
    rows: list[dict] = []
    for hv_level, path in PLAN_FILES.items():
        if not path.exists():
            print(f"  WARN: {path} not found — skipping")
            continue
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            ma = parse_test_allele(r["Mother_allele(s)"], hv_level, "mother")
            fa = parse_test_allele(r["Father_allele(s)"], hv_level, "father")
            m_fg = allele_to_fg.get(ma, "?")
            f_fg = allele_to_fg.get(fa, "?")
            f_dist = dist_lookup.get((ma, fa))
            hv_cat = r.get("Cross_category_predicted", "")
            new_cat = reclassify_category(m_fg, f_fg, f_dist, hv_cat)
            new_h = infer_new_hypothesis(hv_level, hv_cat, new_cat)
            rows.append({
                "HV_hypothesis": hv_level,
                "Cross_id": r["Cross_id"],
                "Mother": r["Mother"],
                "Father": r["Father"],
                "Mother_test_allele": ma,
                "Father_test_allele": fa,
                "Mother_HV_synonymy_group": r.get("Mother_synonymy_group", ""),
                "Father_HV_synonymy_group": r.get("Father_synonymy_group", ""),
                "Mother_functional_group": m_fg,
                "Father_functional_group": f_fg,
                "HV_distance": r.get("HV_distance", ""),
                "Functional_distance": f_dist,
                "HV_category": hv_cat,
                "Functional_category": new_cat,
                "Category_changed": hv_cat != new_cat,
                "Functional_hypothesis": new_h,
                "Hypothesis_changed": hv_level != new_h,
                "Discriminating": (hv_cat != new_cat) or (hv_level != new_h),
                "BL_pairing": r.get("BL_pairing", ""),
            })

    annot = pd.DataFrame(rows)
    annot.to_csv(OUT_ANNOT, sep="\t", index=False)
    discrim = annot[annot["Discriminating"]].reset_index(drop=True)
    discrim.to_csv(OUT_DISCRIM, sep="\t", index=False)

    # --- Transitions summary ---
    trans = (annot.groupby(["HV_hypothesis", "Functional_hypothesis"], dropna=False)
                  .size().rename("n_crosses").reset_index())
    trans.to_csv(OUT_TRANS, sep="\t", index=False)

    # --- Console summary ---
    n_total = len(annot)
    n_discrim = len(discrim)
    n_unchanged = n_total - n_discrim
    print(f"\n=== Cross plan re-annotation ===")
    print(f"Total crosses annotated:     {n_total}")
    print(f"Unchanged under functional:  {n_unchanged}  ({100 * n_unchanged / n_total:.1f} %)")
    print(f"HYPOTHESIS-DISCRIMINATING:   {n_discrim}  ({100 * n_discrim / n_total:.1f} %)")

    print(f"\nTransition matrix (HV → Functional):")
    pivot = trans.pivot(index="HV_hypothesis",
                        columns="Functional_hypothesis",
                        values="n_crosses").fillna(0).astype(int)
    print(pivot)

    print(f"\nOutputs:")
    print(f"  {OUT_ANNOT}")
    print(f"  {OUT_DISCRIM}   ({n_discrim} rows)")
    print(f"  {OUT_TRANS}")


if __name__ == "__main__":
    main()
