#!/usr/bin/env python3
"""
SRK_BL_integration.py — Bridge script

Joins SRK individual-level data (Canu_amplicon pipeline) to bottleneck-lineage
(BL) assignments derived from the LEPA spatial connectivity analysis
(LEPA_EO_spatial_clustering project).

Why this exists
---------------
Until now the SRK pipeline stratifies analyses by Element Occurrence (EO).
The spatial-clustering project shows that EOs are nested within five
independent bottleneck lineages (BL1-BL5; Ward's D2, silhouette k=5),
which represent distinct demographic units that have evolved without gene
flow between them. Re-framing analyses around BLs lets us treat each
lineage as an independent replicate of the bottleneck-and-drift process
and aggregate the small EO/locality samples that EO-level statistics
cannot use.

Inputs
------
1. sampling_metadata.csv               — raw sample registry (Canu_amplicon)
2. SRK_individual_zygosity.tsv         — canonical post-QC list (272 ind.)
3. EO_group_BL_summary.csv             — LEPA_EO_spatial_clustering/data/

Output
------
SRK_individual_BL_assignments.tsv      — one row per post-QC individual:
    Individual, Pop, EO, Group, BL, Drift_index, BL_status

`BL_status` is one of:
    Assigned       — Pop matched directly to an EO present in the BL summary
    Inferred       — Pop matched via prefix split or manual override
    Unassigned     — no EO match; user must extend the override dict
"""

from __future__ import annotations

import csv
import os
import sys
from collections import Counter, defaultdict

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLING_META = os.path.join(HERE, "sampling_metadata.csv")
ZYGOSITY_TSV = os.path.join(HERE, "SRK_individual_zygosity.tsv")
BL_SUMMARY = os.path.join(
    os.path.dirname(os.path.dirname(HERE)),
    "LEPA_EO_spatial_clustering",
    "data",
    "EO_group_BL_summary.csv",
)
OUT_TSV = os.path.join(HERE, "SRK_individual_BL_assignments.tsv")

# ---------------------------------------------------------------------------
# Manual overrides for Pop codes that do not match an EO label directly.
#
# Resolution rules:
#   - Compound codes like "26-3" or "24-14" are auto-resolved by splitting
#     on "-" and matching the prefix to a known EO.
#   - "15" -> EO18 (locationID 15 is in EO18 group 12, BL5; no EO15 exists).
# Add new entries below as you cross-reference Peggy_EOs_Germplasm_*.csv
# from the LEPA_EO_spatial_clustering project.
# ---------------------------------------------------------------------------
POP_TO_EO_OVERRIDE: dict[str, str] = {
    "15": "EO18",
    # Add more as resolved:
    # "405":     "EO??",
    # "704":     "EO??",
    # "709":     "EO??",
    # "702-3":   "EO??",
    # "702-24":  "EO??",
    # "712-5":   "EO??",
    # "715-1":   "EO??",
    # "96-16":   "EO??",
    # "84-2":    "EO??",
    # "98-2":    "EO??",
    # "97-6":    "EO??",
}


def normalise_pop_to_eo(
    pop: str, known_eos: set[str]
) -> tuple[str | None, str]:
    """Map a raw `Pop` value from sampling_metadata.csv to an EO label.

    A label is only returned if it exists in `known_eos` (the set of EOs
    listed in EO_group_BL_summary.csv). Returns
    (EO_label_or_None, status) where status is "Assigned", "Inferred",
    or "Unassigned".
    """
    pop = (pop or "").strip()
    if not pop:
        return None, "Unassigned"
    if pop in POP_TO_EO_OVERRIDE:
        eo = POP_TO_EO_OVERRIDE[pop]
        return (eo, "Inferred") if eo in known_eos else (None, "Unassigned")

    def numeric_to_eo(s: str) -> str | None:
        if not s.isdigit():
            return None
        candidate = f"EO{int(s):02d}"
        return candidate if candidate in known_eos else None

    # Direct numeric match (e.g. "27" -> EO27).
    eo = numeric_to_eo(pop)
    if eo:
        return eo, "Assigned"

    # Compound codes like "26-3" or "24-14": split on "-" and match the
    # prefix as an EO number.
    if "-" in pop:
        prefix = pop.split("-", 1)[0]
        eo = numeric_to_eo(prefix)
        if eo:
            return eo, "Inferred"

    return None, "Unassigned"


def load_eo_to_bl(path: str) -> tuple[dict[str, dict], list[dict]]:
    """Read EO_group_BL_summary.csv -> per-EO BL mapping + raw rows.

    A single EO can span multiple groups within the same BL (e.g. EO27 in
    BL4 has groups 8/9/10/28). We collapse to one record per EO; the
    `Group` field becomes a comma-separated list and `Drift_index` is the
    mean across that EO's groups.
    """
    rows: list[dict] = []
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            rows.append(r)

    by_eo: dict[str, dict] = {}
    grouped: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        grouped[r["EO"]].append(r)

    for eo, group_rows in grouped.items():
        bls = sorted({g["BL"] for g in group_rows})
        if len(bls) > 1:
            sys.stderr.write(
                f"WARNING: EO {eo} spans multiple BLs {bls}; using first.\n"
            )
        groups = sorted({g["Group"] for g in group_rows}, key=lambda x: int(x))
        dis = [float(g["Drift_index"]) for g in group_rows if g["Drift_index"]]
        by_eo[eo] = {
            "BL": bls[0],
            "Group": ", ".join(groups),
            "Drift_index": round(sum(dis) / len(dis), 4) if dis else "",
        }
    return by_eo, rows


def load_individuals(path: str) -> list[dict]:
    """Read sampling_metadata.csv and return ingroup individuals only."""
    out: list[dict] = []
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            if (r.get("Ingroup") or "").strip() != "1":
                continue
            out.append(r)
    return out


def load_post_qc_ids(path: str) -> set[str]:
    """Return the canonical post-QC individual IDs (from zygosity TSV)."""
    ids: set[str] = set()
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            ids.add(r["Individual"])
    return ids


def main() -> int:
    for label, path in [
        ("BL summary", BL_SUMMARY),
        ("sampling metadata", SAMPLING_META),
        ("zygosity TSV", ZYGOSITY_TSV),
    ]:
        if not os.path.isfile(path):
            sys.stderr.write(f"ERROR: {label} not found at {path}\n")
            return 1

    eo_to_bl, _ = load_eo_to_bl(BL_SUMMARY)
    known_eos = set(eo_to_bl)
    individuals = load_individuals(SAMPLING_META)
    post_qc_ids = load_post_qc_ids(ZYGOSITY_TSV)

    individuals = [r for r in individuals if r["SampleID"] in post_qc_ids]
    missing = post_qc_ids - {r["SampleID"] for r in individuals}
    if missing:
        sys.stderr.write(
            f"WARNING: {len(missing)} post-QC individuals are missing from "
            f"sampling_metadata.csv: {sorted(missing)[:5]}...\n"
        )

    out_rows: list[dict] = []
    bl_counter: Counter[str] = Counter()
    bl_eo_counter: dict[str, Counter[str]] = defaultdict(Counter)
    unassigned_pops: Counter[str] = Counter()

    for r in individuals:
        sample_id = r["SampleID"]
        pop = (r.get("Pop") or "").strip()
        eo, status = normalise_pop_to_eo(pop, known_eos)

        bl = group = drift = ""
        if eo and eo in eo_to_bl:
            info = eo_to_bl[eo]
            bl = info["BL"]
            group = info["Group"]
            drift = info["Drift_index"]
        else:
            unassigned_pops[pop or "<empty>"] += 1

        out_rows.append(
            {
                "Individual": sample_id,
                "Pop": pop,
                "EO": eo or "",
                "Group": group,
                "BL": bl or "Unassigned",
                "Drift_index": drift,
                "BL_status": status,
            }
        )

        bl_counter[bl or "Unassigned"] += 1
        if eo:
            bl_eo_counter[bl or "Unassigned"][eo] += 1

    # Write output
    fieldnames = [
        "Individual",
        "Pop",
        "EO",
        "Group",
        "BL",
        "Drift_index",
        "BL_status",
    ]
    with open(OUT_TSV, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    # Console report
    print(f"Wrote {len(out_rows)} rows to {OUT_TSV}\n")

    print("BL × individual counts")
    print("-" * 50)
    for bl in sorted(bl_counter):
        n = bl_counter[bl]
        eos = bl_eo_counter.get(bl, {})
        eo_str = ", ".join(f"{eo}={c}" for eo, c in sorted(eos.items()))
        print(f"  {bl:<12} n={n:<4} EOs: {eo_str}")
    print()

    if unassigned_pops:
        print("Unassigned Pop codes (extend POP_TO_EO_OVERRIDE to resolve):")
        print("-" * 50)
        for p, c in unassigned_pops.most_common():
            print(f"  Pop={p!r:<10} n={c}")
        print()
    else:
        print("All individuals assigned to a BL.\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
