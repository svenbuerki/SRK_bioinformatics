#!/usr/bin/env python3
"""
SRK_injection_donor_ranking.py — Step 19 (TP1 companion)

For each recipient that the Step 17/18 ranking flags for allele injection,
rank candidate donor SOURCES (other EO / BL aggregates) by the number of
novel S-alleles they would contribute.

Recipients
----------
Every group in Tables/SRK_EO_allele_richness.tsv EXCEPT EO27. EO27 is the
only focal group that sits in the "INFORMED BREEDING (frequency skew)"
quadrant on the Depletion x P_compat figure (DI_predicted = 0.475, on the
mild side of 0.50), so it is the only one that does not require injection
at this stage.

Donors
------
Every OTHER group is a candidate, with two exclusions:
  - the recipient itself;
  - any aggregate that subsumes the recipient (the recipient EO's parent
    BL aggregate, or a BL recipient's child EOs) — these contain the
    recipient's own individuals and would just be self-donation.

Scoring
-------
  novel(D->R)   = |alleles(D) \\ alleles(R)|      # novel S-alleles brought
  shared(D->R)  = |alleles(D) & alleles(R)|       # already present in R
  k_after       = k(R) + novel(D->R)
  pct_species   = novel / (K_SPECIES - k(R)) * 100
                  (fraction of the species-level gap that this donor closes)
  pct_observed  = novel / (k_pool - k(R)) * 100
                  (fraction of the OBSERVED gap; k_pool = alleles seen across
                   all groups in this dataset)

Ranking tiers
-------------
  1. within-BL   donor BL == recipient BL  (local adaptation; minimal
                                              drift-load mismatch)
  2. cross-BL    donor BL != recipient BL
Within each tier, donors sort by: novel desc, shared desc, donor N desc.

Outputs
-------
  Tables/SRK_injection_donor_ranking.tsv   one row per recipient x donor
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parent
TABLES = REPO / "Tables"

GENOTYPES = TABLES / "Phase2" / "step11_individual_allele_genotypes.tsv"
METADATA = TABLES / "sampling_metadata.csv"
BL_FILE = TABLES / "Phase3" / "step13_individual_BL_assignments.tsv"
OUT = TABLES / "Phase4" / "step25_injection_donor_ranking.tsv"

K_SPECIES = 59
MIN_N = 15
EXCLUDED_RECIPIENTS = {("EO", "27")}
SUBCODE_POOL = {"27pv": "27", "27rt": "27"}


def load_metadata() -> dict[str, str]:
    out: dict[str, str] = {}
    with METADATA.open(encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            if (row.get("Ingroup") or "").strip() != "1":
                continue
            sid = row["SampleID"].strip()
            eo = (row.get("EO_w_sub") or "").strip()
            if not eo:
                continue
            out[sid] = SUBCODE_POOL.get(eo, eo)
    return out


def load_bl() -> dict[str, str]:
    out: dict[str, str] = {}
    with BL_FILE.open(encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["Individual"].strip()] = (row.get("BL") or "").strip()
    return out


def infer_tetraploid_set(raw: dict[str, int]) -> set[str]:
    """Identity of the distinct alleles after tetraploid copy inference.

    Mirrors SRK_zygosity_from_genotype.R: 1->AAAA, 2->AAAB/AABB, 3->AABC,
    4->ABCD. Identity of distinct alleles is unchanged by copy-filling,
    so we just return the set of observed allele IDs (capped at 4).
    """
    present = sorted([(a, c) for a, c in raw.items() if c > 0],
                     key=lambda x: -x[1])
    return {a for a, _ in present[:4]}


def load_allele_sets(
    meta: dict[str, str], bl: dict[str, str]
) -> tuple[set[str], dict[str, set[str]], dict[str, set[str]],
           Counter, Counter, dict[str, str]]:
    alleles_total: set[str] = set()
    indiv_alleles: dict[str, set[str]] = {}
    with GENOTYPES.open(encoding="utf-8-sig", newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        all_alleles = header[1:]
        for row in rdr:
            sid = row[0]
            if sid not in meta:
                continue
            raw = {a: int(v) for a, v in zip(all_alleles, row[1:])}
            s = infer_tetraploid_set(raw)
            if not s:
                continue
            indiv_alleles[sid] = s
            alleles_total |= s

    by_eo: dict[str, set[str]] = {}
    eo_n: Counter = Counter()
    by_bl: dict[str, set[str]] = {}
    bl_n: Counter = Counter()
    eo_bl_counter: dict[str, Counter] = {}
    for sid, alls in indiv_alleles.items():
        eo = meta[sid]
        by_eo.setdefault(eo, set()).update(alls)
        eo_n[eo] += 1
        b = bl.get(sid, "")
        if b and b != "Unassigned":
            by_bl.setdefault(b, set()).update(alls)
            bl_n[b] += 1
            eo_bl_counter.setdefault(eo, Counter())[b] += 1

    modal_bl = {eo: (c.most_common(1)[0][0] if c else "")
                for eo, c in eo_bl_counter.items()}
    return alleles_total, by_eo, by_bl, eo_n, bl_n, modal_bl


def main() -> None:
    meta = load_metadata()
    bl = load_bl()
    alleles_total, by_eo, by_bl, eo_n, bl_n, modal_bl = load_allele_sets(meta, bl)
    k_pool = len(alleles_total)

    groups: list[dict] = []
    for b in sorted(by_bl):
        groups.append({"level": "BL", "group": b, "bl": b,
                       "n": bl_n[b], "alleles": by_bl[b]})
    focal_eos = [eo for eo in by_eo if eo_n[eo] >= MIN_N]
    for eo in sorted(focal_eos, key=lambda x: -eo_n[x]):
        groups.append({"level": "EO", "group": eo, "bl": modal_bl.get(eo, ""),
                       "n": eo_n[eo], "alleles": by_eo[eo]})

    # Recipients: BL aggregates first (worst k first), then EOs (worst k first)
    bl_rec = [g for g in groups if g["level"] == "BL"
              and (g["level"], g["group"]) not in EXCLUDED_RECIPIENTS]
    eo_rec = [g for g in groups if g["level"] == "EO"
              and (g["level"], g["group"]) not in EXCLUDED_RECIPIENTS]
    bl_rec.sort(key=lambda g: len(g["alleles"]))
    eo_rec.sort(key=lambda g: len(g["alleles"]))
    recipients = bl_rec + eo_rec

    rows: list[dict] = []
    for recipient in recipients:
        rkey = (recipient["level"], recipient["group"])
        r_alleles = recipient["alleles"]
        k_r = len(r_alleles)

        candidates = []
        for donor in groups:
            dkey = (donor["level"], donor["group"])
            if dkey == rkey:
                continue
            # Recipient EO's parent BL aggregate -> contains recipient
            if recipient["level"] == "EO" and donor["level"] == "BL" \
               and donor["group"] == recipient["bl"]:
                continue
            # BL recipient's child EOs -> contained in recipient
            if recipient["level"] == "BL" and donor["level"] == "EO" \
               and donor["bl"] == recipient["group"]:
                continue
            novel = donor["alleles"] - r_alleles
            shared = donor["alleles"] & r_alleles
            same_bl = (donor["bl"] != "" and donor["bl"] == recipient["bl"])
            candidates.append({"donor": donor, "novel": novel,
                                "shared": shared, "same_bl": same_bl})

        for tier_same_bl in (True, False):
            tier = [c for c in candidates if c["same_bl"] == tier_same_bl]
            tier.sort(key=lambda c: (-len(c["novel"]), -len(c["shared"]),
                                       -c["donor"]["n"]))
            for rank, c in enumerate(tier, 1):
                novel_n = len(c["novel"])
                shared_n = len(c["shared"])
                k_after = k_r + novel_n
                pct_species = (novel_n / (K_SPECIES - k_r) * 100
                                if k_r < K_SPECIES else 0.0)
                pct_obs = (novel_n / (k_pool - k_r) * 100
                            if k_r < k_pool else 0.0)
                rows.append({
                    "recipient_level": recipient["level"],
                    "recipient_group": recipient["group"],
                    "recipient_BL": recipient["bl"],
                    "recipient_N": recipient["n"],
                    "k_recipient": k_r,
                    "donor_level": c["donor"]["level"],
                    "donor_group": c["donor"]["group"],
                    "donor_BL": c["donor"]["bl"],
                    "donor_N": c["donor"]["n"],
                    "k_donor": len(c["donor"]["alleles"]),
                    "tier": "within-BL" if tier_same_bl else "cross-BL",
                    "rank_in_tier": rank,
                    "novel_alleles": novel_n,
                    "shared_alleles": shared_n,
                    "k_after_injection": k_after,
                    "pct_species_gap_closed": f"{pct_species:.1f}",
                    "pct_observed_gap_closed": f"{pct_obs:.1f}",
                    "novel_allele_ids": (",".join(sorted(c["novel"]))
                                          if c["novel"] else "-"),
                })

    with OUT.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"k_pool (alleles observed across all groups) = {k_pool}")
    print(f"K_SPECIES (MM consensus)                    = {K_SPECIES}")
    print(f"Excluded recipients: {sorted(EXCLUDED_RECIPIENTS)}")
    print(f"Recipients ({len(recipients)}): "
          + ", ".join(f"{r['level']}{r['group']}" for r in recipients))
    print(f"\nWrote {OUT.relative_to(REPO)} ({len(rows)} rows)\n")

    # Top donor per recipient per tier
    top = [r for r in rows if r["rank_in_tier"] == 1]
    cols = ["recipient_level", "recipient_group", "recipient_BL",
            "k_recipient", "tier", "donor_level", "donor_group",
            "donor_BL", "novel_alleles", "shared_alleles",
            "k_after_injection", "pct_species_gap_closed"]
    col_w = {k: max(len(k), max(len(str(r[k])) for r in top)) for k in cols}
    line = "  ".join(k.ljust(col_w[k]) for k in cols)
    print("Top donor per recipient per tier:")
    print(line)
    print("-" * len(line))
    for r in top:
        print("  ".join(str(r[k]).ljust(col_w[k]) for k in cols))


if __name__ == "__main__":
    main()
