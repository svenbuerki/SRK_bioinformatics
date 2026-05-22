#!/usr/bin/env python3
"""
Per-EO S-allele richness, evenness, and tetraploid compatibility diagnostics.

Output: Tables/SRK_EO_allele_richness.tsv

Scope:
  Six EOs with N >= 15 genotyped ingroup individuals
  (18, 25, 27 [pooling 27pv + 27rt], 67, 70, 76).
  Rarefied to N = 30 over 1000 permutations.

Ploidy handling:
  LEPA is tetraploid. Amplicon sequencing under-recovers copies, so
  per-individual allele copies are inferred to 4 using the same rule as
  SRK_zygosity_from_genotype.R:
      1 distinct allele observed                  -> AAAA (4 copies of A)
      2 distinct, dominant >= 2x minor            -> AAAB (3+1)
      2 distinct, roughly equal                   -> AABB (2+2)
      3 distinct                                  -> AABC (2+1+1)
      4 distinct                                  -> ABCD (1+1+1+1)
  Every individual therefore contributes exactly 4 inferred copies.
  Denominator for allele frequencies in an EO is C = 4N.

Compatibility model:
  Sporophytic SI, tetraploid stigma (up to 4 distinct alleles), diploid
  pollen (up to 2 distinct alleles), co-dominance. Pollen is rejected if
  any of its 2 alleles matches any of the stigma's 4. Exact strict-SI
  formula under multinomial allele sampling:
      P_compat_strict = sum_{a,b} p_a p_b (1 - p_a - p_b*I[a!=b])^4
  Equal-frequency benchmark: 1 - 8/k.

SI leakage (C3 Stage 5):
  At leakage rate L in [0,1], a fraction L of would-be-incompatible
  matings nevertheless succeed (selfing or SI-leaky outcrossing):
      P_compat(L) = P_compat_strict + L * (1 - P_compat_strict)
  Reported at L in {0, 0.10, 0.25, 0.50}.

  Empirical leakage estimator from observed AAAA proportion:
      L_hat = prop_AAAA / 3.5
  Rationale: AAAA cannot arise under strict outcrossing SI; observed
  AAAA reflects historical selfing/sib-crossing events. The 3.5 divisor
  is a population-averaged tetrasomic correction (an AABB parent selfed
  yields AAAA in ~1/6 of offspring, AAAB in ~1/2; averaged across the
  AABB/AAAB mix typical of the data, the effective correction is ~3-4).
  L_hat is an UPPER BOUND because some apparent AAAA arise from amplicon
  under-recovery rather than true homozygosity.

Inheritance-mode assumption (IMPORTANT):
  The exact P_compat formula above assumes tetrasomic inheritance with
  random pairing across all four chromosomes (autotetraploid-like
  behaviour). If LEPA's SRK locus shows DISOMIC inheritance (i.e. the
  two homoeologous subgenomes from an allotetraploid origin segregate
  independently), the system behaves more like two superimposed diploid
  SI systems, which is generally MORE permissive than the unified
  tetraploid model used here. In that case P_compat as computed is a
  conservative lower bound. The formula should be revisited if cytology
  or segregation data establish disomic inheritance at SRK.

Output rows are mixed EO-level and BL-level:
  - EO level: one row per parent EO with N >= 15 (six EOs).
  - BL level: one row per assigned BL (BL1-BL5), pooling all genotyped
    individuals in that BL regardless of EO sample size.
  The `level` column flags which is which.

Columns:
  level                    "EO" or "BL"
  group                    EO id (e.g. "27") or BL label (e.g. "BL3")
  BL                       parent BL (same as group for BL rows;
                            modal BL among genotyped individuals for EO rows)
  N                        genotyped individuals in the group
  geno_mix                 inferred genotype tally (AAAA:AAAB:AABB:AABC:ABCD)
  prop_AAAA                fraction of individuals classed AAAA
  L_hat_from_AAAA          empirical leakage estimate = prop_AAAA / 3.5
  raw_copy_recovery        mean raw amplicon copies recovered per individual
  k_observed               distinct alleles present after inference
  k_rarefied30_mean        rarefied k (subsample 30 individuals, 1000 perms)
  k_rarefied30_sd
  shannon_H                Shannon entropy of inferred allele frequencies
  evenness_J               H / ln(k_observed); 1.0 = perfectly even
  chi2_uniform_p           p-value, chi-square vs uniform 1/k expectation
  P_compat_uniform_4x      1 - 8/k_observed (tetraploid equal-freq benchmark)
  P_compat_L0              strict SI: exact formula, L = 0
  P_compat_L0.10           with mild leakage, L = 0.10
  P_compat_L0.25           with substantial leakage, L = 0.25
  P_compat_L0.50           with heavy leakage, L = 0.50
"""

from __future__ import annotations

import csv
import math
import random
from collections import Counter
from pathlib import Path

from scipy.stats import chisquare

REPO = Path(__file__).resolve().parent
TABLES = REPO / "Tables"

GENOTYPES = TABLES / "SRK_individual_allele_genotypes.tsv"
METADATA = TABLES / "sampling_metadata.csv"
BL_FILE = REPO / "SRK_individual_BL_assignments.tsv"
OUT = TABLES / "SRK_EO_allele_richness.tsv"

MIN_N = 15
RAREFY_N = 30
N_PERMS = 1000
SEED = 20260522

SUBCODE_POOL = {"27pv": "27", "27rt": "27"}

GENOTYPE_CLASSES = ["AAAA", "AAAB", "AABB", "AABC", "ABCD"]

LEAK_LADDER = [0.0, 0.10, 0.25, 0.50]
AAAA_CORRECTION = 3.5  # divisor for L_hat = prop_AAAA / AAAA_CORRECTION


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


def infer_tetraploid(raw: dict[str, int]) -> tuple[dict[str, int], str]:
    """Mirror SRK_zygosity_from_genotype.R: fill to 4 copies."""
    present = sorted([(a, c) for a, c in raw.items() if c > 0],
                     key=lambda x: -x[1])
    n = len(present)
    if n == 0:
        return {}, "Unknown"
    if n == 1:
        a, _ = present[0]
        return {a: 4}, "AAAA"
    if n == 2:
        (a, ca), (b, cb) = present
        if ca >= 2 * cb:
            return {a: 3, b: 1}, "AAAB"
        return {a: 2, b: 2}, "AABB"
    if n == 3:
        (a, _), (b, _), (c, _) = present
        return {a: 2, b: 1, c: 1}, "AABC"
    # n >= 4 -> take top 4
    (a, _), (b, _), (c, _), (d, _) = present[:4]
    return {a: 1, b: 1, c: 1, d: 1}, "ABCD"


def load_inferred_genotypes(meta: dict[str, str]) -> tuple[
    list[str],                              # allele ids in raw order
    dict[str, dict[str, int]],              # SampleID -> inferred 4-copy dict
    dict[str, str],                         # SampleID -> genotype class
    dict[str, int],                         # SampleID -> raw copies recovered
]:
    with GENOTYPES.open(encoding="utf-8-sig", newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        alleles = header[1:]
        inferred: dict[str, dict[str, int]] = {}
        klass: dict[str, str] = {}
        raw_copies: dict[str, int] = {}
        for row in rdr:
            sid = row[0]
            if sid not in meta:
                continue
            raw = {a: int(v) for a, v in zip(alleles, row[1:])}
            inf, cls = infer_tetraploid(raw)
            if not inf:
                continue
            inferred[sid] = inf
            klass[sid] = cls
            raw_copies[sid] = sum(raw.values())
    return alleles, inferred, klass, raw_copies


def copy_vector(samples: list[str],
                inferred: dict[str, dict[str, int]]) -> Counter:
    c: Counter = Counter()
    for sid in samples:
        for a, n in inferred[sid].items():
            c[a] += n
    return c


def rarefied_k(samples: list[str],
                inferred: dict[str, dict[str, int]],
                n: int,
                perms: int,
                rng: random.Random) -> tuple[float, float]:
    if len(samples) <= n:
        return float(len(copy_vector(samples, inferred))), 0.0
    ks = []
    for _ in range(perms):
        sub = rng.sample(samples, n)
        ks.append(len(copy_vector(sub, inferred)))
    mean = sum(ks) / len(ks)
    var = sum((x - mean) ** 2 for x in ks) / (len(ks) - 1)
    return mean, math.sqrt(var)


def shannon(freqs: list[float]) -> float:
    return -sum(p * math.log(p) for p in freqs if p > 0)


def p_compat_tetraploid(freqs: list[float]) -> float:
    """Exact tetraploid sporophytic SI compatibility under multinomial sampling.

    Sums over ordered pollen allele pairs (a,b); stigma must contain neither.
    """
    total = 0.0
    n = len(freqs)
    for i in range(n):
        pi = freqs[i]
        # Pollen draws (i, i)
        total += pi * pi * (1 - pi) ** 4
        for j in range(n):
            if i == j:
                continue
            pj = freqs[j]
            # Pollen draws (i, j), i != j
            rem = 1 - pi - pj
            if rem < 0:
                rem = 0.0
            total += pi * pj * (rem ** 4)
    return total


def compute_row(level: str,
                 group: str,
                 bl_label: str,
                 sids: list[str],
                 inferred: dict[str, dict[str, int]],
                 klass: dict[str, str],
                 raw_copies: dict[str, int],
                 rng: random.Random) -> dict:
    n = len(sids)
    copies = copy_vector(sids, inferred)
    total = sum(copies.values())  # = 4n
    k_obs = len(copies)
    freqs = [copies[a] / total for a in copies]

    H = shannon(freqs)
    J = H / math.log(k_obs) if k_obs > 1 else float("nan")

    expected = total / k_obs
    chi2_p = chisquare(list(copies.values()),
                        f_exp=[expected] * k_obs).pvalue if k_obs > 1 else float("nan")

    k_mean, k_sd = rarefied_k(sids, inferred, RAREFY_N, N_PERMS, rng)

    p_uni = max(0.0, 1 - 8 / k_obs) if k_obs >= 2 else 0.0
    p_strict = p_compat_tetraploid(freqs)

    mix_counter = Counter(klass[s] for s in sids)
    mix_str = ":".join(str(mix_counter.get(g, 0)) for g in GENOTYPE_CLASSES)
    prop_aaaa = mix_counter.get("AAAA", 0) / n
    l_hat = prop_aaaa / AAAA_CORRECTION

    mean_raw = sum(raw_copies[s] for s in sids) / n

    row = {
        "level": level,
        "group": group,
        "BL": bl_label,
        "N": n,
        "geno_mix": mix_str,
        "prop_AAAA": f"{prop_aaaa:.3f}",
        "L_hat_from_AAAA": f"{l_hat:.3f}",
        "raw_copy_recovery": f"{mean_raw:.2f}",
        "k_observed": k_obs,
        "k_rarefied30_mean": f"{k_mean:.2f}",
        "k_rarefied30_sd": f"{k_sd:.2f}",
        "shannon_H": f"{H:.3f}",
        "evenness_J": f"{J:.3f}",
        "chi2_uniform_p": f"{chi2_p:.3g}",
        "P_compat_uniform_4x": f"{p_uni:.3f}",
    }
    for L in LEAK_LADDER:
        tag = "L0" if L == 0.0 else f"L{L:.2f}"
        row[f"P_compat_{tag}"] = f"{p_strict + L * (1 - p_strict):.3f}"
    return row


def main() -> None:
    rng = random.Random(SEED)
    meta = load_metadata()
    bl = load_bl()
    _, inferred, klass, raw_copies = load_inferred_genotypes(meta)

    # EO-level grouping (parent EO, sub-codes already pooled in meta)
    eo_to_samples: dict[str, list[str]] = {}
    for sid, eo in meta.items():
        if sid in inferred:
            eo_to_samples.setdefault(eo, []).append(sid)

    focal_eo = {eo: sids for eo, sids in eo_to_samples.items() if len(sids) >= MIN_N}
    if not focal_eo:
        raise SystemExit(f"No EO meets N >= {MIN_N}")

    # BL-level grouping (every genotyped individual with an assigned BL)
    bl_to_samples: dict[str, list[str]] = {}
    for sid in inferred:
        b = bl.get(sid, "")
        if b and b != "Unassigned":
            bl_to_samples.setdefault(b, []).append(sid)

    rows = []

    # BL rows first (top of table); sorted alphanumerically for now
    for b in sorted(bl_to_samples):
        sids = bl_to_samples[b]
        rows.append(compute_row("BL", b, b, sids, inferred, klass, raw_copies, rng))

    # EO rows; sorted by N desc within EO block
    for eo in sorted(focal_eo, key=lambda x: -len(focal_eo[x])):
        sids = focal_eo[eo]
        bls = [bl.get(s, "") for s in sids if bl.get(s, "") and bl.get(s) != "Unassigned"]
        bl_modal = Counter(bls).most_common(1)[0][0] if bls else ""
        rows.append(compute_row("EO", eo, bl_modal, sids, inferred, klass, raw_copies, rng))

    with OUT.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    col_w = {k: max(len(k), max(len(str(r[k])) for r in rows)) for k in rows[0]}
    line = "  ".join(k.ljust(col_w[k]) for k in rows[0])
    print(line)
    print("-" * len(line))
    for r in rows:
        print("  ".join(str(r[k]).ljust(col_w[k]) for k in rows[0]))
    n_bl = sum(1 for r in rows if r["level"] == "BL")
    n_eo = sum(1 for r in rows if r["level"] == "EO")
    print(f"\ngeno_mix order: {':'.join(GENOTYPE_CLASSES)}")
    print(f"Wrote {OUT.relative_to(REPO)} ({n_bl} BL rows + {n_eo} EO rows, "
          f"rarefied to N={RAREFY_N}, {N_PERMS} perms)")


if __name__ == "__main__":
    main()
