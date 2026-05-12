#!/usr/bin/env python3
"""
srk_cross_plan.py   (Step 22e — hypothesis-testing cross plan generator)

Translates the Step 22b cross categorisation into a phased experimental
protocol that tests S-allele specificity hypotheses with explicit
genotype constraints. Generates one TSV per hypothesis level (H0, H1a,
H1b, H2, H3) listing concrete maternal/paternal individual IDs, plus a
summary figure showing phase counts and the decision tree.

Hypothesis framework
====================
Each sequence-defined allele bin is a hypothesis: "the proteins clustered
into this bin share one SI specificity." Crosses test these hypotheses
by comparing predicted outcome (Incompatible / Synonymy_test /
Compatible_within / Compatible_cross) with observed seed yield.

| Level | Question | Maternal | Paternal | Constraint | Why |
|-------|----------|----------|----------|------------|-----|
| H0 | Does SI rejection work? | AAAA | AAAA | both in SAME synonymy group | Unambiguous specificity both sides; 0 seeds is the prediction |
| H1a | Within-class compatible baseline | AAAA | AAAA | different synonymy groups, NO Synonymy_test edge (HV >= 0.04) | Distinct specificities, predicted compatible |
| H1b | Between-class compatible (max yield) | AAAA Class I | Heterozygous carrier (AAAB/AABB/AABC) of the Class II allele | mother allele != Class I main allele(s) of father | Class II has no AAAA carriers in the current dataset — heterozygous father gives at least 50% pollen carrying the Class II specificity |
| H2 | Synonymy bin boundaries | AAAA | AAAA | different synonymy groups, WITH Synonymy_test edge (0 < HV < 0.04) | The core synonymy test — outcome decides whether two bins represent one or two specificities |
| H3 | Hidden bins (no AAAA representative) | AAAA carrier of allele M | AAAB carrier of "hidden" allele B (and main allele A) | M != A AND M != B (so AA pollen is compatible — gives baseline), plus paired control AAAA × AAAA(A only) | 29 bins exist only in heterozygotes; AAAB father gives 50% AA + 50% AB pollen, requires paired baseline cross |

Inputs
------
  SRK_AAAA_cross_design_HV.tsv      Step 22b — pre-computed cross design for AAAA × AAAA pairs
  SRK_synonymy_groups.csv            Step 22b — synonymy-group membership per allele
  SRK_HV_allele_distances.tsv        Step 22b — HV distance matrix (for AAAB pairings)
  SRK_functional_allele_groups.tsv   Step 22b — allele -> Class I / Class II
  SRK_individual_BL_assignments.tsv  Step 13 — BL per individual
  SRK_individual_zygosity.tsv        Step 12 — genotype patterns
  SRK_individual_allele_table.tsv    Step 11 — allele copies per individual

Outputs
-------
  SRK_cross_plan_H0_SI_validation.tsv               H0 (Incompatible negative controls)
  SRK_cross_plan_H1a_within_class_baseline.tsv      H1a (Compatible_within positive controls)
  SRK_cross_plan_H1b_between_class_baseline.tsv     H1b (Compatible_cross via heterozygous Class II carrier)
  SRK_cross_plan_H2_synonymy_tests.tsv              H2 (Synonymy_test, one per inter-group bridge)
  SRK_cross_plan_H3_hidden_bin_tests.tsv            H3 (AAAB-mediated tests for the 29 hidden bins)
  SRK_cross_plan_summary.tsv                         per-phase counts + cumulative cross attempts
  SRK_cross_plan_summary.{pdf,png}                   phase timeline + decision tree figure
"""

from __future__ import annotations
import os
import re
import sys
import math
import itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from collections import defaultdict


_COMP_RE = re.compile(r"(Allele_\d+)\((\d+)\)")


def parse_allele_composition(comp_str: str) -> list[tuple[str, int]]:
    """Parse 'Allele_051(1)+Allele_049(2)' → [('Allele_051', 1), ('Allele_049', 2)].
    The (N) values are protein counts in the original genotyping output —
    NOT copy counts. For AAAB, the allele with the smaller protein count is
    the B (1-copy) allele; the other is the A (3-copy) main allele.
    """
    if not isinstance(comp_str, str):
        return []
    return [(m.group(1), int(m.group(2)))
            for m in _COMP_RE.finditer(comp_str)]

# =============================================================================
# Settings
# =============================================================================
CROSS_DESIGN_TSV  = "SRK_AAAA_cross_design_HV.tsv"
SYN_GROUPS_CSV    = "SRK_synonymy_groups.csv"
HV_DIST_TSV       = "SRK_HV_allele_distances.tsv"
FUNC_GROUPS_TSV   = "SRK_functional_allele_groups.tsv"
BL_TSV            = "SRK_individual_BL_assignments.tsv"
ZYGO_TSV          = "SRK_individual_zygosity.tsv"
ALLELE_TSV        = "SRK_individual_allele_table.tsv"

WITHIN_CLASS_THRESHOLD = 0.04   # must match Step 22b

# Per-phase replicate-count recommendations
REPLICATES = {
    "H0":  9,   # 3 mothers × 3 flowers
    "H1a": 9,
    "H1b": 9,
    "H2":  15,  # 5 mothers × 3 flowers (uncertain outcome → more reps)
    "H3":  15,
}

# Number of crosses to schedule per category (caps to keep workload tractable)
N_CROSSES = {
    "H0_within_allele":   5,    # same-allele AAAA × AAAA in Synonymy group 1
    "H0_within_group":    8,    # different-alleles-same-group in Synonymy group 1
    "H0_other_groups":   12,    # 2 per remaining synonymy group (8 groups)
    "H1a_far":           10,    # max-HV-distance pairs, between groups, no Synonymy_test edge
    "H1b":                5,    # heterozygous Class II carrier × AAAA Class I
    "H2_per_pair":        1,    # one cross per synonymy-group-pair Synonymy_test bridge
    "H3_per_bin":         1,    # one cross per "hidden" bin, plus a paired control
}

# Set1 BL palette (consistent with rest of pipeline)
BL_PALETTE = {
    "BL1": "#E41A1C", "BL2": "#377EB8", "BL3": "#4DAF4A",
    "BL4": "#984EA3", "BL5": "#FF7F00",
}

OUT_TSVS = {
    "H0":  "SRK_cross_plan_H0_SI_validation.tsv",
    "H1a": "SRK_cross_plan_H1a_within_class_baseline.tsv",
    "H1b": "SRK_cross_plan_H1b_between_class_baseline.tsv",
    "H2":  "SRK_cross_plan_H2_synonymy_tests.tsv",
    "H3":  "SRK_cross_plan_H3_hidden_bin_tests.tsv",
}
OUT_SUMMARY_TSV = "SRK_cross_plan_summary.tsv"
OUT_FIG_PDF     = "SRK_cross_plan_summary.pdf"
OUT_FIG_PNG     = os.path.join("figures", "SRK_cross_plan_summary.png")

CROSS_PLAN_COLS = [
    "Hypothesis", "Cross_id", "Priority",
    "Mother", "Mother_genotype", "Mother_allele(s)", "Mother_BL",
    "Mother_synonymy_group",
    "Father", "Father_genotype", "Father_allele(s)", "Father_BL",
    "Father_synonymy_group",
    "BL_pairing", "Cross_category_predicted", "HV_distance",
    "Replicates_recommended", "Expected_outcome", "Decision_rule",
    "Paired_control", "Notes",
]

os.makedirs("figures", exist_ok=True)

# =============================================================================
# 1. Load inputs and build per-individual lookups
# =============================================================================
print("Loading inputs ...")

for f in (CROSS_DESIGN_TSV, SYN_GROUPS_CSV, HV_DIST_TSV, FUNC_GROUPS_TSV,
          BL_TSV, ZYGO_TSV, ALLELE_TSV):
    if not os.path.exists(f):
        sys.exit(f"ERROR: missing input {f}")

cross_design = pd.read_csv(CROSS_DESIGN_TSV, sep="\t")
syn_groups   = pd.read_csv(SYN_GROUPS_CSV)
func_groups  = pd.read_csv(FUNC_GROUPS_TSV, sep="\t")
bl_df        = pd.read_csv(BL_TSV, sep="\t", encoding="utf-8-sig")
zygo_df      = pd.read_csv(ZYGO_TSV, sep="\t", encoding="utf-8-sig")
allele_df    = pd.read_csv(ALLELE_TSV, sep="\t", encoding="utf-8-sig")
hv_dist_df   = pd.read_csv(HV_DIST_TSV, sep="\t", index_col=0)

print(f"  cross_design:  {len(cross_design)} AAAA×AAAA pairs")
print(f"  syn_groups:    {len(syn_groups)} alleles, "
      f"{syn_groups['Synonymy_group'].nunique()} groups")
print(f"  individuals:   {len(zygo_df)} with zygosity, "
      f"{len(bl_df)} with BL")

# allele → synonymy_group, AAAA_count
allele_to_group = dict(zip(syn_groups["Allele"], syn_groups["Synonymy_group"]))
allele_aaaa_n   = dict(zip(syn_groups["Allele"], syn_groups["AAAA_count"]))
# allele → Class
allele_to_class = dict(zip(func_groups["Allele"],
                            func_groups["FunctionalGroup"]))

# individual → BL
ind_to_bl = {row["Individual"]: row["BL"] for _, row in bl_df.iterrows()
             if row.get("BL_status") != "Unassigned"}

# individual → genotype
ind_to_geno = dict(zip(zygo_df["Individual"], zygo_df["Genotype"]))

# individual → list of (allele, copy_count) tuples (sorted by copy desc)
ind_to_alleles = defaultdict(list)
for _, r in allele_df.iterrows():
    ind_to_alleles[r["Individual"]].append((r["Allele"], int(r["Count"])))
# sort each individual's alleles by copy count desc
for k in ind_to_alleles:
    ind_to_alleles[k].sort(key=lambda x: -x[1])

# Build helper dicts:
# AAAA individuals (one allele × 4): allele -> list of individuals carrying it
allele_to_aaaa_inds = defaultdict(list)
for ind, geno in ind_to_geno.items():
    if geno != "AAAA":
        continue
    if ind not in ind_to_alleles or not ind_to_alleles[ind]:
        continue
    # AAAA: pick the single allele
    allele = ind_to_alleles[ind][0][0]
    allele_to_aaaa_inds[allele].append(ind)

# Build per-individual parsed allele composition (uses zygosity TSV's authoritative
# Allele_composition string rather than the per-protein allele table)
ind_to_comp = {}
for _, row in zygo_df.iterrows():
    ind = row["Individual"]
    comp = parse_allele_composition(row.get("Allele_composition", ""))
    ind_to_comp[ind] = comp


def carriers_of(target_allele: str,
                allowed_genotypes: tuple[str, ...] = ("AAAB", "AABB", "AABC", "ABCD"),
                require_minor: bool = False):
    """Return list of (ind, target_n, other_alleles, BL, genotype) for
    individuals carrying target_allele in any of `allowed_genotypes`.

    target_n        — the (N) protein count of the target allele in
                      Allele_composition (informative; NOT a copy count)
    other_alleles   — list of (allele, N) for the rest of the composition
                      (this is what the AA/AB pollen will carry)
    require_minor   — if True, require target_allele to be the *minor* allele
                      (smallest protein count, i.e. B in AAAB → 1 copy of B);
                      this is what we want for H3 hidden-bin tests.
    """
    out = []
    for _, row in zygo_df.iterrows():
        ind = row["Individual"]
        geno = row.get("Genotype")
        if geno not in allowed_genotypes:
            continue
        comp = ind_to_comp.get(ind, [])
        if not comp:
            continue
        if not any(a == target_allele for a, _ in comp):
            continue
        target_n = next(n for a, n in comp if a == target_allele)
        if require_minor:
            min_n = min(n for _, n in comp)
            if target_n != min_n:
                continue
        others = [(a, n) for a, n in comp if a != target_allele]
        out.append((ind, target_n, others,
                    ind_to_bl.get(ind, "Unassigned"), geno))
    return out


def diff_bl_score(m_bl: str | None, f_bl: str | None) -> int:
    """Return 1 if mother and father are in different known BLs, 0 otherwise.
    Used as a tiebreaker to prefer between-BL pairings."""
    if m_bl is None or f_bl is None:
        return 0
    if m_bl == "Unassigned" or f_bl == "Unassigned":
        return 0
    return int(m_bl != f_bl)


def make_row(hyp: str, cross_id: str, priority: int,
             mother: str, father: str,
             mother_alleles: list[str], father_alleles: list[str],
             cross_cat: str, hv_d,
             expected: str, decision: str,
             paired: str = "", notes: str = "") -> dict:
    m_bl = ind_to_bl.get(mother, "Unassigned")
    f_bl = ind_to_bl.get(father, "Unassigned")
    pair = ("between-BL" if (m_bl != f_bl and m_bl != "Unassigned"
                              and f_bl != "Unassigned") else
            "within-BL"  if m_bl == f_bl and m_bl != "Unassigned" else "any")
    return {
        "Hypothesis":               hyp,
        "Cross_id":                 cross_id,
        "Priority":                 priority,
        "Mother":                   mother,
        "Mother_genotype":          ind_to_geno.get(mother, "?"),
        "Mother_allele(s)":         "+".join(mother_alleles),
        "Mother_BL":                m_bl,
        "Mother_synonymy_group":    "+".join(allele_to_group.get(a, "?")
                                             for a in mother_alleles),
        "Father":                   father,
        "Father_genotype":          ind_to_geno.get(father, "?"),
        "Father_allele(s)":         "+".join(father_alleles),
        "Father_BL":                f_bl,
        "Father_synonymy_group":    "+".join(allele_to_group.get(a, "?")
                                             for a in father_alleles),
        "BL_pairing":               pair,
        "Cross_category_predicted": cross_cat,
        "HV_distance":              hv_d if hv_d is not None else "",
        "Replicates_recommended":   REPLICATES[hyp],
        "Expected_outcome":         expected,
        "Decision_rule":            decision,
        "Paired_control":           paired,
        "Notes":                    notes,
    }


# =============================================================================
# 2. H0 — SI validation (Incompatible negative controls)
# =============================================================================
print("\nGenerating H0 (SI validation, predicted Incompatible) ...")
h0_rows: list[dict] = []
priority = 0

# Get the synonymy groups, sorted by AAAA count
group_to_alleles = defaultdict(list)
for _, r in syn_groups.iterrows():
    if r["Synonymy_group"] != "Isolated":
        group_to_alleles[r["Synonymy_group"]].append(r["Allele"])
sorted_groups = sorted(group_to_alleles.keys(),
                       key=lambda g: -sum(allele_aaaa_n[a]
                                          for a in group_to_alleles[g]))


def pick_aaaa_pair(allele_a: str, allele_b: str, used_pairs: set,
                   prefer_between_bl: bool = True):
    """Pick (mother, father) AAAA carriers of (allele_a, allele_b) preferring
    between-BL. Return None if no candidates."""
    moms = allele_to_aaaa_inds.get(allele_a, [])
    dads = allele_to_aaaa_inds.get(allele_b, [])
    if not moms or not dads:
        return None
    candidates = []
    for m in moms:
        for d in dads:
            if m == d:
                continue
            pkey = tuple(sorted([m, d]))
            if pkey in used_pairs:
                continue
            score = diff_bl_score(ind_to_bl.get(m), ind_to_bl.get(d))
            candidates.append((score, m, d))
    if not candidates:
        return None
    candidates.sort(key=lambda x: -x[0])
    return candidates[0][1], candidates[0][2]


used = set()

# 5 same-allele AAAA × AAAA crosses within Synonymy group 1
group1_alleles = group_to_alleles.get("Synonymy group 1", [])
group1_alleles_by_aaaa = sorted(group1_alleles,
                                 key=lambda a: -allele_aaaa_n.get(a, 0))
for allele in group1_alleles_by_aaaa[:8]:
    if len([r for r in h0_rows if "same-allele" in r["Notes"]]) >= N_CROSSES["H0_within_allele"]:
        break
    pair = pick_aaaa_pair(allele, allele, used)
    if pair is None:
        continue
    used.add(tuple(sorted(pair)))
    priority += 1
    h0_rows.append(make_row(
        "H0", f"H0-{priority:03d}", priority,
        pair[0], pair[1], [allele], [allele],
        "Incompatible", 0.0,
        "0 seeds (SI rejection — same allele on both sides)",
        "0 seeds → SI works; >0 seeds → SI broken at this allele",
        notes=f"same-allele AAAA × AAAA  ({allele})",
    ))

# 8 different-alleles-within-Synonymy-group-1 crosses
g1_pairs_seen = set()
for a, b in itertools.combinations(group1_alleles_by_aaaa, 2):
    if len([r for r in h0_rows if "different-alleles" in r["Notes"]
            and "Synonymy group 1" in r["Mother_synonymy_group"]]) >= N_CROSSES["H0_within_group"]:
        break
    pkey = tuple(sorted([a, b]))
    if pkey in g1_pairs_seen:
        continue
    g1_pairs_seen.add(pkey)
    pair = pick_aaaa_pair(a, b, used)
    if pair is None:
        continue
    used.add(tuple(sorted(pair)))
    priority += 1
    h0_rows.append(make_row(
        "H0", f"H0-{priority:03d}", priority,
        pair[0], pair[1], [a], [b],
        "Incompatible", 0.0,
        "0 seeds (HV-identical alleles in same Synonymy group)",
        "0 seeds → both alleles share specificity; >0 seeds → bin boundary breach",
        notes=f"different-alleles-same-group  (Synonymy group 1: {a} × {b})",
    ))

# 12 crosses across remaining synonymy groups (~2 per group, 6 groups)
n_other = 0
for g in sorted_groups:
    if g == "Synonymy group 1":
        continue
    if n_other >= N_CROSSES["H0_other_groups"]:
        break
    alleles_g = sorted(group_to_alleles[g],
                       key=lambda a: -allele_aaaa_n.get(a, 0))
    # within-group pairs (same allele, then different alleles)
    picks_in_group = 0
    for allele in alleles_g[:2]:
        if allele_aaaa_n.get(allele, 0) >= 2:
            pair = pick_aaaa_pair(allele, allele, used)
            if pair is not None:
                used.add(tuple(sorted(pair)))
                priority += 1
                n_other += 1
                picks_in_group += 1
                h0_rows.append(make_row(
                    "H0", f"H0-{priority:03d}", priority,
                    pair[0], pair[1], [allele], [allele],
                    "Incompatible", 0.0,
                    "0 seeds (same allele)",
                    "0 seeds → SI rejection on this specificity confirmed",
                    notes=f"same-allele AAAA × AAAA  ({g}: {allele})",
                ))
                if picks_in_group >= 2 or n_other >= N_CROSSES["H0_other_groups"]:
                    break
    if picks_in_group < 2 and len(alleles_g) >= 2 and n_other < N_CROSSES["H0_other_groups"]:
        a, b = alleles_g[0], alleles_g[1]
        pair = pick_aaaa_pair(a, b, used)
        if pair is not None:
            used.add(tuple(sorted(pair)))
            priority += 1
            n_other += 1
            h0_rows.append(make_row(
                "H0", f"H0-{priority:03d}", priority,
                pair[0], pair[1], [a], [b],
                "Incompatible", 0.0,
                "0 seeds (HV-identical, same Synonymy group)",
                "0 seeds → both alleles share specificity",
                notes=f"different-alleles-same-group  ({g}: {a} × {b})",
            ))

print(f"  H0 crosses: {len(h0_rows)}")

# =============================================================================
# 3. H1a — within-Class compatible baseline
# =============================================================================
print("\nGenerating H1a (within-Class compatible baseline) ...")
h1a_rows: list[dict] = []
priority = 0

# Filter cross_design for within-class, no Synonymy_test edge (HV >= 0.04)
h1a_pool = cross_design[(cross_design["Category"] == "Compatible_within")].copy()
# rank by HV distance descending — pick maximally divergent within-class pairs
h1a_pool = h1a_pool.sort_values("HV_dist", ascending=False)

# Track used (mother, father) pairs to avoid re-using same individuals
h1a_used_pairs: set = set()
h1a_used_inds: set = set()
for _, r in h1a_pool.iterrows():
    if len(h1a_rows) >= N_CROSSES["H1a_far"]:
        break
    m, f = r["Mother"], r["Father"]
    if m == f:
        continue
    if m in h1a_used_inds or f in h1a_used_inds:
        continue   # avoid reusing individuals
    pkey = tuple(sorted([m, f]))
    if pkey in h1a_used_pairs:
        continue
    h1a_used_pairs.add(pkey)
    h1a_used_inds.add(m)
    h1a_used_inds.add(f)
    priority += 1
    h1a_rows.append(make_row(
        "H1a", f"H1a-{priority:03d}", priority,
        m, f, [r["Mother_allele"]], [r["Father_allele"]],
        "Compatible_within", float(r["HV_dist"]),
        "high seed yield (different specificities, same Class)",
        "Anchors the Compatible_within baseline; below-baseline yield = "
        "unexpected — flag for repeat",
        notes=f"max-HV-divergent within-Class pair",
    ))

print(f"  H1a crosses: {len(h1a_rows)}")

# =============================================================================
# 4. H1b — between-Class compatible (max yield) via heterozygous Class II carrier
# =============================================================================
# Class II is identified as the minority functional group from Step 22b
# (whichever FG does NOT contain the majority of alleles).
fg_majority = max(set(allele_to_class.values()), key=lambda c: sum(1 for v in allele_to_class.values() if v == c))
fg_class2 = next((c for c in set(allele_to_class.values()) if c != fg_majority), "FG02")
class2_alleles = sorted(a for a, c in allele_to_class.items() if c == fg_class2)
class2_label = ", ".join(class2_alleles) if class2_alleles else "(none detected)"

print(f"\nGenerating H1b (between-Class compatible via heterozygous Class II carrier — "
      f"Class II allele(s): {class2_label}) ...")
h1b_rows: list[dict] = []
priority = 0

# Find heterozygous carriers of any Class II allele (AAAB/AABB/AABC/ABCD).
# The Class II allele typically has no AAAA carriers in LEPA datasets, so we
# rely on a heterozygous mother (AABB pollen distribution: 17% AA + 67% AB +
# 17% BB → 84% of pollen carries the target Class II specificity, AAAB: 50%).
class2_carriers = []
for c2 in class2_alleles:
    class2_carriers.extend(carriers_of(c2))
print(f"  Class II alleles: {class2_alleles}")
print(f"  Class II carriers found (any heterozygous genotype): "
      f"{len(class2_carriers)}")
for c in class2_carriers[:5]:
    print(f"    {c[0]}  geno={c[4]}  others={c[2]}")

if not class2_carriers:
    print("  WARNING: no Class II carriers found — H1b cannot be scheduled.")
else:
    # Mothers: AAAA Class I, allele ≠ Class II carrier's other Class I alleles
    class1_aaaa_inds = []
    for allele, inds in allele_to_aaaa_inds.items():
        if allele_to_class.get(allele) != fg_class2:
            class1_aaaa_inds.extend([(ind, allele) for ind in inds])

    used_pairs = set()
    for father_ind, target_n, others, f_bl, f_geno in class2_carriers:
        if len(h1b_rows) >= N_CROSSES["H1b"]:
            break
        # The father's "other" alleles are Class I alleles that AA pollen carries.
        # The mother allele must differ from any of them (so AA pollen is
        # compatible at the Class I level too).
        father_other_alleles = [a for a, _ in others]
        candidate_mothers = [
            (m, m_allele) for m, m_allele in class1_aaaa_inds
            if m_allele not in father_other_alleles
        ]
        if not candidate_mothers:
            continue
        candidate_mothers.sort(
            key=lambda mt: -diff_bl_score(ind_to_bl.get(mt[0]), f_bl))
        for m, m_allele in candidate_mothers:
            if (m, father_ind) in used_pairs:
                continue
            used_pairs.add((m, father_ind))
            priority += 1
            father_alleles = father_other_alleles + class2_alleles
            h1b_rows.append(make_row(
                "H1b", f"H1b-{priority:03d}", priority,
                m, father_ind, [m_allele], father_alleles,
                "Compatible_cross",
                "0.5–1.0 (between Class I and Class II)",
                "high seed yield — pollen carrying the Class II allele is "
                "between-class compatible with the Class I mother; pollen "
                "carrying only the father's other Class I alleles is also "
                "compatible because mother allele differs from them.",
                "Anchors the upper compatibility baseline (between-Class)",
                notes=f"{f_geno} Class II carrier; father other alleles = "
                      f"{father_other_alleles}; pollen distribution "
                      f"depends on father genotype",
            ))
            break

print(f"  H1b crosses: {len(h1b_rows)}")

# =============================================================================
# 5. H2 — synonymy tests (one cross per synonymy-group-pair Synonymy_test bridge)
# =============================================================================
print("\nGenerating H2 (synonymy tests) ...")
h2_rows: list[dict] = []
priority = 0

# Find inter-synonymy-group Synonymy_test pairs in cross_design.
syn_pool = cross_design[cross_design["Category"] == "Synonymy_test"].copy()
# Add synonymy_group columns
syn_pool["Mother_group"] = syn_pool["Mother_allele"].map(allele_to_group)
syn_pool["Father_group"] = syn_pool["Father_allele"].map(allele_to_group)
# Inter-group only (different groups, neither Isolated to the same allele)
syn_pool = syn_pool[syn_pool["Mother_group"] != syn_pool["Father_group"]]
# Build group-pair key (sorted)
syn_pool["GroupPair"] = syn_pool.apply(
    lambda r: tuple(sorted([str(r["Mother_group"]),
                            str(r["Father_group"])])), axis=1)

# For each unique group-pair, pick the cross with both individuals in
# different BLs (preferred), then highest combined AAAA availability
syn_pool["m_bl"] = syn_pool["Mother"].map(ind_to_bl)
syn_pool["f_bl"] = syn_pool["Father"].map(ind_to_bl)
syn_pool["bl_score"] = syn_pool.apply(
    lambda r: diff_bl_score(r["m_bl"], r["f_bl"]), axis=1)
syn_pool["m_aaaa"] = syn_pool["Mother_allele"].map(allele_aaaa_n)
syn_pool["f_aaaa"] = syn_pool["Father_allele"].map(allele_aaaa_n)
syn_pool["aaaa_total"] = syn_pool["m_aaaa"] + syn_pool["f_aaaa"]
syn_pool = syn_pool.sort_values(["bl_score", "aaaa_total"],
                                 ascending=[False, False])

picked_pairs = set()
seen_individuals = set()
for _, r in syn_pool.iterrows():
    if r["GroupPair"] in picked_pairs:
        continue
    if r["Mother"] in seen_individuals or r["Father"] in seen_individuals:
        # try to find a different individual pair for this group-pair
        continue
    picked_pairs.add(r["GroupPair"])
    seen_individuals.add(r["Mother"])
    seen_individuals.add(r["Father"])
    priority += 1
    h2_rows.append(make_row(
        "H2", f"H2-{priority:03d}", priority,
        r["Mother"], r["Father"], [r["Mother_allele"]], [r["Father_allele"]],
        "Synonymy_test", float(r["HV_dist"]),
        "uncertain — bimodal: 0 seeds (synonymy → merge groups) or "
        "near-baseline yield (boundary real → keep separate)",
        "0 seeds → MERGE the two synonymy groups (synonymous specificities). "
        ">0 seeds at H1a-baseline level → bin boundary confirmed.",
        notes=f"group-pair: {r['Mother_group']} × {r['Father_group']}",
    ))

# Allow second-pass picks for group-pairs that didn't get covered due to
# individual reuse — relax the seen-individuals constraint
for _, r in syn_pool.iterrows():
    if r["GroupPair"] in picked_pairs:
        continue
    picked_pairs.add(r["GroupPair"])
    priority += 1
    h2_rows.append(make_row(
        "H2", f"H2-{priority:03d}", priority,
        r["Mother"], r["Father"], [r["Mother_allele"]], [r["Father_allele"]],
        "Synonymy_test", float(r["HV_dist"]),
        "uncertain — bimodal",
        "0 seeds → merge; >0 seeds → keep separate",
        notes=f"group-pair: {r['Mother_group']} × {r['Father_group']} "
              f"(allows individual reuse)",
    ))

print(f"  H2 crosses: {len(h2_rows)}  "
      f"(group-pairs covered: {len(picked_pairs)})")

# =============================================================================
# 6. H3 — hidden bins (no AAAA representative)
# =============================================================================
print("\nGenerating H3 (hidden bins via AAAB donors) ...")
h3_rows: list[dict] = []
priority = 0

# Hidden bins = alleles with 0 AAAA carriers
hidden_alleles = [a for a in allele_to_class
                  if allele_aaaa_n.get(a, 0) == 0]
print(f"  hidden alleles (no AAAA): {len(hidden_alleles)}")

for hidden in sorted(hidden_alleles):
    # Prefer AAAB carriers where the hidden allele is the minor (1-copy) B;
    # fall back to AABB / AABC / ABCD if no AAAB carrier exists.
    carriers = carriers_of(hidden, allowed_genotypes=("AAAB",), require_minor=True)
    if not carriers:
        carriers = carriers_of(hidden,
                               allowed_genotypes=("AAAB", "AABB", "AABC", "ABCD"))
    if not carriers:
        continue

    father_ind, target_n, others, f_bl, f_geno = carriers[0]
    # main_a = the most common "other" allele in the father
    main_a = max(others, key=lambda x: x[1])[0] if others else None
    if main_a is None:
        continue
    # Find AAAA mothers carrying alleles M with M != main_a AND M != hidden
    # AND main_a vs M is Compatible_within (so the AA pollen baseline is compatible)
    candidate_mothers = []
    for ind, geno in ind_to_geno.items():
        if geno != "AAAA":
            continue
        if ind == father_ind:
            continue
        m_allele = ind_to_alleles[ind][0][0]
        if m_allele == main_a or m_allele == hidden:
            continue
        # check that M and main_a are Compatible_within
        try:
            d_main = float(hv_dist_df.loc[m_allele, main_a])
        except KeyError:
            continue
        if d_main < WITHIN_CLASS_THRESHOLD:
            continue   # not Compatible_within
        # both must be in same Class for Compatible_within
        if (allele_to_class.get(m_allele) != allele_to_class.get(main_a)):
            continue
        # measure mother vs hidden allele HV distance (the actual hypothesis)
        try:
            d_hidden = float(hv_dist_df.loc[m_allele, hidden])
        except KeyError:
            d_hidden = float("nan")
        candidate_mothers.append((ind, m_allele, d_hidden, d_main))
    if not candidate_mothers:
        continue
    # prefer between-BL with father
    candidate_mothers.sort(
        key=lambda c: -diff_bl_score(ind_to_bl.get(c[0]), f_bl))
    m_ind, m_allele, d_mhidden, d_main = candidate_mothers[0]

    # Paired control: same mother × an AAAA father carrying main_a only
    # (no rare allele) → measures the AA-only pollen baseline
    paired_father_ind = None
    if main_a in allele_to_aaaa_inds:
        # any AAAA carrier of main_a, prefer different individual from father_ind
        for cand in allele_to_aaaa_inds[main_a]:
            if cand != m_ind:
                paired_father_ind = cand
                break

    paired_id = ""
    if paired_father_ind:
        priority += 1
        paired_id = f"H3-{priority:03d}"
        h3_rows.append(make_row(
            "H3", paired_id, priority,
            m_ind, paired_father_ind, [m_allele], [main_a],
            "Compatible_within", float(d_main),
            "high seed yield (AA-only pollen baseline)",
            "Establishes the AA-pollen baseline yield for the paired H3 cross; "
            "subtract from the AAAB cross to isolate the hidden allele's effect",
            notes=f"PAIRED CONTROL for hidden allele {hidden}",
        ))

    # Main H3 cross
    priority += 1
    main_id = f"H3-{priority:03d}"
    cat = "Compatible_within" if (allele_to_class.get(m_allele) ==
                                   allele_to_class.get(hidden)) else "Compatible_cross"
    father_alleles_listed = [main_a, hidden]
    h3_rows.append(make_row(
        "H3", main_id, priority,
        m_ind, father_ind, [m_allele], father_alleles_listed,
        cat, d_mhidden,
        "uncertain — depends on whether the hidden allele is compatible "
        "with the mother. Interpret as: total yield ≈ paired-control yield "
        "→ hidden allele incompatible with M; total yield ≈ 2 × paired-control "
        "yield → hidden allele compatible with M",
        "Compare total yield to paired control. Higher than control = hidden "
        "allele's specificity is functional and compatible with M.",
        paired=paired_id if paired_id else "(no paired control found)",
        notes=f"hidden allele = {hidden}; father genotype = {f_geno}; "
              f"main_A = {main_a}; M = {m_allele}",
    ))

print(f"  H3 crosses: {len(h3_rows)} (hidden bins covered: "
      f"{len([r for r in h3_rows if 'PAIRED CONTROL' not in r['Notes']])})")

# =============================================================================
# 7. Write per-phase TSVs and summary
# =============================================================================
print("\nWriting per-phase TSVs ...")
phases = [
    ("H0",  h0_rows),
    ("H1a", h1a_rows),
    ("H1b", h1b_rows),
    ("H2",  h2_rows),
    ("H3",  h3_rows),
]

for hyp, rows in phases:
    if not rows:
        print(f"  {hyp}: empty (no rows)")
        continue
    df = pd.DataFrame(rows, columns=CROSS_PLAN_COLS)
    out = OUT_TSVS[hyp]
    df.to_csv(out, sep="\t", index=False)
    print(f"  {hyp}: {len(rows)} crosses → {out}")

# Phase summary
summ_rows = []
total_attempts = 0
for hyp, rows in phases:
    n_crosses = len(rows)
    reps      = REPLICATES[hyp]
    attempts  = n_crosses * reps
    total_attempts += attempts
    summ_rows.append({
        "Hypothesis":      hyp,
        "Description":     {
            "H0":  "SI validation (Incompatible negative controls)",
            "H1a": "Within-Class compatible baseline (Compatible_within positive)",
            "H1b": f"Between-Class compatible baseline (Compatible_cross via heterozygous carrier of Class II allele: {class2_label})",
            "H2":  "Synonymy bin-boundary tests (Synonymy_test)",
            "H3":  "Hidden-bin tests via AAAB donors (29 bins without AAAA)",
        }[hyp],
        "N_crosses":              n_crosses,
        "Replicates_per_cross":   reps,
        "Total_attempts":         attempts,
    })
summ_df = pd.DataFrame(summ_rows)
summ_df.loc[len(summ_df)] = {
    "Hypothesis": "TOTAL", "Description": "All hypothesis tests",
    "N_crosses": int(summ_df["N_crosses"].sum()),
    "Replicates_per_cross": "—",
    "Total_attempts": total_attempts,
}
summ_df.to_csv(OUT_SUMMARY_TSV, sep="\t", index=False)
print(f"\nSummary table written to {OUT_SUMMARY_TSV}")
print(summ_df.to_string(index=False))

# =============================================================================
# 8. Summary figure: phase counts + decision tree
# =============================================================================
print(f"\nWriting summary figure ...")

fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(15, 6.5),
                                  gridspec_kw={"width_ratios": [1, 1.3]})

# Left: bar chart of cross counts per phase
phase_labels = ["H0\nSI validation", "H1a\nC_within baseline",
                "H1b\nC_cross baseline", "H2\nSynonymy tests",
                "H3\nHidden bins"]
counts = [len(rows) for _, rows in phases]
attempts = [c * REPLICATES[h] for (h, _), c in zip(phases, counts)]
colors = ["#762a83", "#1b7837", "#5aae61", "#f1a340", "#3182bd"]
bars = ax_l.bar(phase_labels, counts, color=colors, edgecolor="black", linewidth=0.5)
for bar, n_c, n_a in zip(bars, counts, attempts):
    ax_l.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
              f"{n_c} crosses\n({n_a} attempts)", ha="center", va="bottom", fontsize=9)
ax_l.set_ylabel("Cross count", fontsize=11)
ax_l.set_title("Phased cross plan — hypothesis tests for S-allele specificity",
                fontsize=11)
ax_l.set_ylim(0, max(counts) * 1.25 if counts else 1)
ax_l.grid(axis="y", linestyle=":", alpha=0.4)
for spine in ("top", "right"):
    ax_l.spines[spine].set_visible(False)

# Right: decision tree (text-based hierarchical layout)
ax_r.set_xlim(0, 10)
ax_r.set_ylim(0, 10)
ax_r.axis("off")
ax_r.set_title("Decision logic per phase outcome", fontsize=11)

def box(x, y, w, h, text, fc, fontsize=8.5, weight="normal"):
    ax_r.add_patch(FancyBboxPatch((x, y), w, h,
                                   boxstyle="round,pad=0.18",
                                   fc=fc, ec="black", lw=0.6))
    ax_r.text(x + w/2, y + h/2, text, ha="center", va="center",
              fontsize=fontsize, fontweight=weight, wrap=True)

# H0 → if pass (0 seeds W) proceed; if fail (seeds in W) → SI broken
box(0.2, 8.6, 9.6, 1.0,
    "H0 — SI validation\n→ pass: ≥80% Incompatible crosses give 0 seeds → proceed.\n"
    "Fail: SI broken → revise framework before H1–H3.",
    "#e7d4e8", weight="bold")

box(0.2, 7.2, 4.5, 1.1,
    "H1a — Compatible_within\nMeasures within-Class compatible yield (baseline)",
    "#d9f0d3")
box(5.3, 7.2, 4.5, 1.1,
    "H1b — Compatible_cross\nMeasures between-Class compatible yield (max)",
    "#d9f0d3")

box(0.2, 5.6, 9.6, 1.2,
    "H2 — Synonymy test outcomes\n"
    "• 0 seeds → MERGE the two synonymy groups (synonymous specificities)\n"
    "• Yield ≥ H1a baseline → KEEP SEPARATE (bin boundary functionally real)\n"
    "• Intermediate → repeat with different individuals",
    "#fee08b", weight="normal")

box(0.2, 3.8, 9.6, 1.4,
    "H3 — Hidden-bin tests (AAAB donors)\n"
    "Compare main cross yield to paired-control yield:\n"
    "• Total ≈ paired control → hidden allele incompatible with mother\n"
    "• Total ≈ 2× paired control → hidden allele compatible (additive)\n"
    "• Total ≈ 1.5× → partial compatibility (mixed pollen accepted)",
    "#cfe1f2", weight="normal")

box(0.2, 1.6, 9.6, 1.6,
    "End of Phase 4 → Validated Functional S-allele Table\n"
    "  • Synonymy groups confirmed synonymous: merged into single specificities\n"
    "  • Bin boundaries confirmed: kept as separate specificities\n"
    "  • Hidden-bin alleles classified compatible / incompatible with major specificities\n"
    "→ Input to Phase 5 seed-orchard design (operational)",
    "#fff7bc", weight="bold")

# Connecting arrows
for x_start, y_start, x_end, y_end in [
    (5.0, 8.6, 5.0, 8.3),
    (2.5, 7.2, 2.5, 6.8),
    (7.5, 7.2, 7.5, 6.8),
    (5.0, 5.6, 5.0, 5.2),
    (5.0, 3.8, 5.0, 3.4),
]:
    ax_r.annotate("", xy=(x_end, y_end), xytext=(x_start, y_start),
                  arrowprops=dict(arrowstyle="->", color="#555555", lw=1.0))

plt.tight_layout()
plt.savefig(OUT_FIG_PDF, format="pdf", bbox_inches="tight")
plt.savefig(OUT_FIG_PNG, format="png", dpi=200, bbox_inches="tight")
plt.close()
print(f"  Saved {OUT_FIG_PDF}")
print(f"  Saved {OUT_FIG_PNG}")

print("\nDone.")
