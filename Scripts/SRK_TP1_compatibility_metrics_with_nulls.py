#!/usr/bin/env python3
"""
SRK_TP1_compatibility_metrics_with_nulls.py — Step 17b

Null-aware variant of SRK_TP1_compatibility_metrics.py. Reads the augmented
`SRK_individual_allele_genotypes_with_nulls.tsv` from Step 26 (rows already
sum to 4 with explicit Allele_NULL copies) and recomputes every TP1 metric
per group.

Modelling decision (Sven, 2026-06-11)
─────────────────────────────────────
Allele_NULL is treated as a regular allele in the frequency vector
(p_NULL > 0 contributes to allele-frequency denominators per the
Hardy-Weinberg decision in this thread). The exact tetraploid P_compat
formula is applied unchanged. Biologically this is a heuristic — a true
null-aware compatibility model would treat NULL pollen as automatically
compatible and NULL stigma copies as unable to reject — but the heuristic
preserves comparability with the functional-only baseline, and the
direction of bias is conservative: the functional-only P_compat
*overestimates* the mating-pool size for populations with high p_NULL,
so the null-aware P_compat is the more pessimistic (more honest) number.

Side-by-side outputs allow direct comparison with the original baseline.

Inputs
──────
  Tables/SRK_individual_allele_genotypes_with_nulls.tsv  (Step 26 output)
  Tables/sampling_metadata.csv                            (ingroup + EO)
  SRK_individual_BL_assignments.tsv                       (BL per individual)
  SRK_allele_accumulation_stats.tsv                       (MM estimates)

Outputs
───────
  Tables/SRK_EO_allele_richness_with_nulls.tsv            (parallel to canonical)
"""
from __future__ import annotations

import csv
import math
import random
from collections import Counter
from pathlib import Path

from scipy.stats import chisquare

REPO    = Path(__file__).resolve().parent
TABLES  = REPO / "Tables"

GENO_TSV = TABLES / "Phase4" / "step26_individual_allele_genotypes_with_nulls.tsv"
META     = TABLES / "sampling_metadata.csv"
BL_FILE  = TABLES / "Phase3" / "step13_individual_BL_assignments.tsv"
ACC_STATS = TABLES / "Phase3" / "step15_allele_accumulation_stats.tsv"
OUT      = TABLES / "Phase4" / "step17b_EO_allele_richness_with_nulls.tsv"

K_SPECIES = 59
MIN_N     = 15
RAREFY_N  = 30
N_PERMS   = 1000
N_BOOT    = 1000
CI_LO     = 0.025
CI_HI     = 0.975
SEED      = 20260611

SUBCODE_POOL = {"27pv": "27", "27rt": "27"}
LEAK_LADDER  = [0.0, 0.10, 0.25, 0.50]
AAAA_CORRECTION = 3.5  # kept for compatibility with the existing schema


# ─── Helpers ──────────────────────────────────────────────────────────────────
def load_metadata() -> dict[str, str]:
    out: dict[str, str] = {}
    with META.open(encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            if (row.get("Ingroup") or "").strip() != "1":
                continue
            sid = row["SampleID"].strip()
            eo  = (row.get("EO_w_sub") or "").strip()
            if not eo:
                continue
            out[sid] = SUBCODE_POOL.get(eo, eo)
    return out


def load_bl_aug(meta: dict[str, str]) -> dict[str, str]:
    """BL lookup augmented with Step 26 EO/BL columns for SC additions."""
    out: dict[str, str] = {}
    with BL_FILE.open(encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["Individual"].strip()] = (row.get("BL") or "").strip()
    # Augment from the with-nulls genotype table (carries BL for SC additions)
    with GENO_TSV.open(encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sid = row["Individual"].strip()
            bl  = (row.get("BL_inferred") or "").strip()
            if sid not in out and bl:
                out[sid] = bl
    return out


def load_mm_estimates() -> dict[tuple[str, str], float]:
    out: dict[tuple[str, str], float] = {}
    if not ACC_STATS.exists():
        return out
    with ACC_STATS.open(encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                mm_val = float((row.get("MM_estimate") or "").strip())
            except ValueError:
                continue
            level = row["Level"].strip()
            pop   = row["Population"].strip()
            if level == "EO" and pop.startswith("EO"):
                pop = pop[2:]
            out[(level, pop)] = mm_val
    return out


def load_null_aware_genotypes(meta: dict[str, str]) -> tuple[
    list[str], dict[str, dict[str, int]], dict[str, str], dict[str, int],
    dict[str, str], dict[str, str], dict[str, str]
]:
    """Return (allele_cols, copies_by_ind, klass, raw_copies, eo, bl_aug,
    si_status).

    Rows already sum to 4 → no homozygosity inference. SC individuals added
    by Step 26 are kept (their meta-EO entry is patched from the genotype
    table's EO_normalised column when missing from sampling_metadata).
    """
    with GENO_TSV.open(encoding="utf-8-sig", newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        all_cols = rdr.fieldnames
        meta_cols = {"Individual", "Genotype_class_flag",
                     "EO_normalised", "BL_inferred",
                     "SI_status", "pSI_confidence"}
        allele_cols = [c for c in all_cols if c not in meta_cols]
        copies: dict[str, dict[str, int]] = {}
        klass:  dict[str, str] = {}
        raw_cp: dict[str, int] = {}
        eo_aug: dict[str, str] = {}
        si_aug: dict[str, str] = {}
        meta_local = dict(meta)
        for row in rdr:
            sid = row["Individual"].strip()
            eo  = row.get("EO_normalised", "").strip().lstrip("EO") or None
            if sid not in meta_local and eo:
                meta_local[sid] = SUBCODE_POOL.get(eo, eo)
            if sid not in meta_local:
                continue
            inf = {a: int(row[a]) for a in allele_cols if int(row[a]) > 0}
            if not inf or sum(inf.values()) != 4:
                continue
            # Genotype class taking nulls into account
            n_null  = inf.get("Allele_NULL", 0)
            n_func  = 4 - n_null
            func_counts = sorted(
                [c for a, c in inf.items() if a != "Allele_NULL"],
                reverse=True
            )
            if n_func == 0:
                cls = "0000"
            elif n_func == 4:
                if   func_counts == [4]:          cls = "AAAA"
                elif func_counts == [3, 1]:       cls = "AAAB"
                elif func_counts == [2, 2]:       cls = "AABB"
                elif func_counts == [2, 1, 1]:    cls = "AABC"
                else:                              cls = "ABCD"
            else:
                # nf nulls + func_counts functional copies
                letters_func = "ABCD"[:len(func_counts)]
                func_str = "".join(letters_func[i] * c
                                   for i, c in enumerate(func_counts))
                cls = func_str + "0" * n_null
            copies[sid] = inf
            klass[sid]  = cls
            raw_cp[sid] = 4
            eo_aug[sid] = meta_local[sid]
            si_aug[sid] = row.get("SI_status", "").strip()
    return allele_cols, copies, klass, raw_cp, eo_aug, load_bl_aug(eo_aug), si_aug


def shannon(freqs: list[float]) -> float:
    return -sum(p * math.log(p) for p in freqs if p > 0)


def p_compat_tetraploid(freqs: list[float]) -> float:
    total = 0.0
    n = len(freqs)
    for i in range(n):
        pi = freqs[i]
        total += pi * pi * (1 - pi) ** 4
        for j in range(n):
            if i == j:
                continue
            pj  = freqs[j]
            rem = max(0.0, 1 - pi - pj)
            total += pi * pj * (rem ** 4)
    return total


def copy_vector(samples: list[str],
                copies: dict[str, dict[str, int]]) -> Counter:
    c: Counter = Counter()
    for sid in samples:
        for a, n in copies[sid].items():
            c[a] += n
    return c


def rarefied_k(samples: list[str], copies, n: int, perms: int,
                rng: random.Random) -> tuple[float, float]:
    if len(samples) < n:
        n = len(samples)
    if n == 0:
        return 0.0, 0.0
    ks = []
    for _ in range(perms):
        pick = rng.sample(samples, n)
        ks.append(len(copy_vector(pick, copies)))
    mean = sum(ks) / len(ks)
    var = sum((k - mean) ** 2 for k in ks) / len(ks)
    return mean, math.sqrt(var)


def bootstrap_p_compat(samples, copies, leaks, n_boot, rng):
    out = {L: [] for L in leaks}
    if len(samples) == 0:
        return {L: (float("nan"), float("nan")) for L in leaks}
    for _ in range(n_boot):
        pick = [rng.choice(samples) for _ in range(len(samples))]
        cv = copy_vector(pick, copies)
        tot = sum(cv.values())
        freqs = [c / tot for c in cv.values()]
        p_s = p_compat_tetraploid(freqs)
        for L in leaks:
            out[L].append(p_s + L * (1 - p_s))
    res = {}
    for L in leaks:
        s = sorted(out[L])
        res[L] = (s[int(CI_LO * len(s))], s[int(CI_HI * len(s))])
    return res


def compute_row(level: str, group: str, bl_label: str, sids: list[str],
                copies, klass, raw_cp, mm_estimates, si_status, rng) -> dict:
    n = len(sids)
    cv = copy_vector(sids, copies)
    total = sum(cv.values())  # = 4n
    k_obs = len(cv)
    freqs = [cv[a] / total for a in cv]

    H = shannon(freqs)
    J = H / math.log(k_obs) if k_obs > 1 else float("nan")
    expected = total / k_obs
    chi2_p = chisquare(list(cv.values()),
                        f_exp=[expected] * k_obs).pvalue if k_obs > 1 else float("nan")

    k_mean, k_sd = rarefied_k(sids, copies, RAREFY_N, N_PERMS, rng)
    p_uni    = max(0.0, 1 - 8 / k_obs) if k_obs >= 2 else 0.0
    p_strict = p_compat_tetraploid(freqs)
    boot_ci  = bootstrap_p_compat(sids, copies, LEAK_LADDER, N_BOOT, rng)

    n_func_alleles = sum(1 for a in cv if a != "Allele_NULL" and cv[a] > 0)
    n_null_copies  = cv.get("Allele_NULL", 0)
    frac_null      = n_null_copies / total

    n_sc  = sum(1 for s in sids if si_status.get(s, "") == "SC")
    n_psi = sum(1 for s in sids if si_status.get(s, "") == "pSI")
    n_si  = sum(1 for s in sids if si_status.get(s, "") == "SI")

    mm_k = mm_estimates.get((level, group), float("nan"))
    di_obs  = 1 - (k_mean / K_SPECIES) if K_SPECIES > 0 else float("nan")
    di_pred = 1 - (mm_k / K_SPECIES) if (K_SPECIES > 0 and not math.isnan(mm_k)) else float("nan")

    row = {
        "level":             level,
        "group":             group,
        "BL":                bl_label,
        "N":                 n,
        "N_SI":              n_si,
        "N_pSI":             n_psi,
        "N_SC":              n_sc,
        "k_species":         K_SPECIES,
        "k_predicted_MM":    f"{mm_k:.1f}" if not math.isnan(mm_k) else "NA",
        "DI_observed":       f"{di_obs:.3f}",
        "DI_predicted":      f"{di_pred:.3f}" if not math.isnan(di_pred) else "NA",
        "k_observed":        k_obs,
        "k_functional":      n_func_alleles,
        "n_null_copies":     n_null_copies,
        "frac_null":         f"{frac_null:.4f}",
        "k_rarefied30_mean": f"{k_mean:.2f}",
        "k_rarefied30_sd":   f"{k_sd:.2f}",
        "shannon_H":         f"{H:.3f}",
        "evenness_J":        f"{J:.3f}",
        "chi2_uniform_p":    f"{chi2_p:.3g}",
        "P_compat_uniform_4x": f"{p_uni:.3f}",
    }
    for L in LEAK_LADDER:
        tag = "L0" if L == 0.0 else f"L{L:.2f}"
        row[f"P_compat_{tag}"]    = f"{p_strict + L * (1 - p_strict):.3f}"
        lo, hi = boot_ci[L]
        row[f"P_compat_{tag}_lo"] = f"{lo:.3f}"
        row[f"P_compat_{tag}_hi"] = f"{hi:.3f}"
    return row


def main() -> None:
    rng = random.Random(SEED)
    meta = load_metadata()
    mm_estimates = load_mm_estimates()
    _, copies, klass, raw_cp, eo_aug, bl_lookup, si_status = \
        load_null_aware_genotypes(meta)

    if not mm_estimates:
        print(f"WARNING: {ACC_STATS.name} not found; DI_predicted = NA.")

    # EO-level grouping
    eo_to_samples: dict[str, list[str]] = {}
    for sid, eo in eo_aug.items():
        if sid in copies:
            eo_to_samples.setdefault(eo, []).append(sid)
    focal_eo = {eo: sids for eo, sids in eo_to_samples.items()
                if len(sids) >= MIN_N}
    if not focal_eo:
        raise SystemExit(f"No EO meets N >= {MIN_N}")

    # BL-level grouping
    bl_to_samples: dict[str, list[str]] = {}
    for sid in copies:
        b = bl_lookup.get(sid, "")
        if b and b != "Unassigned":
            bl_to_samples.setdefault(b, []).append(sid)

    rows = []
    for b in sorted(bl_to_samples):
        sids = bl_to_samples[b]
        rows.append(compute_row("BL", b, b, sids, copies, klass, raw_cp,
                                mm_estimates, si_status, rng))
    for eo in sorted(focal_eo, key=lambda x: -len(focal_eo[x])):
        sids = focal_eo[eo]
        bls = [bl_lookup.get(s, "") for s in sids
               if bl_lookup.get(s, "") and bl_lookup.get(s) != "Unassigned"]
        bl_modal = Counter(bls).most_common(1)[0][0] if bls else ""
        rows.append(compute_row("EO", eo, bl_modal, sids, copies, klass,
                                raw_cp, mm_estimates, si_status, rng))

    with OUT.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n_bl = sum(1 for r in rows if r["level"] == "BL")
    n_eo = sum(1 for r in rows if r["level"] == "EO")
    print(f"\nWrote {OUT.relative_to(REPO)} ({n_bl} BL + {n_eo} EO rows)")
    print("\nFirst 5 BL rows (compact):")
    for r in rows[:5]:
        print(f"  {r['group']:8s} N={r['N']:>3}  k_func={r['k_functional']:>3}  "
              f"frac_null={r['frac_null']:>6}  P_compat_L0={r['P_compat_L0']:>6}  "
              f"P_compat_L0.25={r['P_compat_L0.25']:>6}")


if __name__ == "__main__":
    main()
