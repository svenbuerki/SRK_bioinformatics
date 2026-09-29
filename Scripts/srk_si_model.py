"""LEPA sporophytic SI model with Class I / Class II dominance (Phase 5 Part 2).

Used by Step 30 to compute the per-mother random-mating pollen
compatibility `P_compat` under the biologically correct sporophytic
model, replacing the diploid per-pollen-allele approximation used in
Phase 5 Part 1's § A.6.

Biological model
----------------
LEPA is tetraploid (see § A.3 of the Phase 5 doc). Each somatic plant
carries 4 SRK allele copies. Sporophytic self-incompatibility means
recognition is decided in the sporophyte generation (2n = 4x in LEPA)
— the pollen coat carries proteins deposited by the pollen parent's
sporophyte tissue, so the stigma decides against the whole parent's
expressed genotype, NOT against the individual haploid gamete's
allele. This module implements the classical Brassicaceae Class I /
Class II dominance rule:

  * **Within a plant**: Class I strictly dominant over Class II.
      - Plant carries ≥ 1 Class I allele → only Class I alleles are
        expressed on pollen and stigma.
      - Plant carries only Class II alleles → all its Class II alleles
        are expressed co-dominantly.
  * **Between plants** (SI check):
      - Pollen parent × stigma parent are COMPATIBLE iff their expressed
        allele sets are DISJOINT.
      - Between-class crosses (Class I × Class II) are therefore ALWAYS
        compatible (no shared Class), matching the H1b baseline in
        `step26e_cross_plan_H1b_between_class_baseline.tsv`.
      - Within-class crosses are compatible iff no shared expressed
        allele.

Configuration
-------------
The Fg → Class mapping is a first-class TSV at
`tables/Phase5/srk_fg_class.tsv` — editable by hand as the biology is
refined. `load_class_map()` reads it into a dict.

The provisional data-driven default (see the rationale column of the
TSV) puts Fgs with ≥ 2 sub-alleles in Class I (FG001–FG006, together
~65 % of P1 frequency) and single-allele Fgs in Class II (~35 %).
FG024 is flagged as REVIEW because it is single-allele but 18 % of the
species — the assignment may need manual reclassification once the SI
literature check is done.

Analytical formulas
-------------------
Under a local Fg-frequency vector `f` and Class I total mass
p_I = Σ_{j ∈ Class I} f_j, the sporophytic P_compat for a mother whose
expressed set M has total mass p(M) = Σ_{j ∈ M} f_j is:

  Case A — mother has ≥ 1 Class I allele (M ⊆ Class I):
      P_compat = (1 − p(M))^4

  Case B — mother has ONLY Class II alleles (M = mother's 4 alleles):
      P_compat = 1 − (1 − p_I)^4  +  (1 − p_I − p(M))^4

Derivation: the father draws 4 alleles i.i.d. from f; his expressed
set is his Class I alleles if any, else all four. Case A: any father
allele in M is a Class I allele in M → excluded from compatibility;
father with 0 Class I alleles has an all-Class-II expressed set → no
overlap with a Class I mother expressed set → compatible. Case B:
any father with ≥ 1 Class I allele expresses only Class I → no overlap
with Class II mother → compatible; a father with 0 Class I alleles is
compatible only if none of his 4 alleles are in the mother's expressed
Class II set.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_CLASS_TSV = Path("Tables/Phase5/srk_fg_class.tsv")
DEFAULT_ZYGOSITY_TSV = Path("Tables/Phase5/srk_zygosity_empirical.tsv")


def load_class_map(tsv_path: Path = DEFAULT_CLASS_TSV) -> dict[str, str]:
    """Return dict[Fg → 'I' | 'II']."""
    if not tsv_path.exists():
        raise FileNotFoundError(f"Missing Fg→Class TSV at {tsv_path}")
    df = pd.read_csv(tsv_path, sep="\t", encoding="utf-8-sig")
    return dict(zip(df["Fg"].astype(str), df["class"].astype(str)))


def build_class_i_mask(fg_labels: list[str],
                        class_map: dict[str, str]) -> np.ndarray:
    """Return a boolean numpy array of length K_fg: True at positions
    where the corresponding Fg is Class I."""
    return np.array(
        [class_map.get(str(fg), "II") == "I" for fg in fg_labels],
        dtype=bool,
    )


def p_compat_sporophytic(mother_fg_idx: np.ndarray,
                          local_f: np.ndarray,
                          class_i_mask: np.ndarray) -> float:
    """Analytical sporophytic P_compat for one mother given a local Fg
    frequency vector and the Class I mask.

    Parameters
    ----------
    mother_fg_idx : (4,) array of int
        The Fg indices of the mother's 4 SRK allele copies (may include
        duplicates — the model uses her distinct expressed set).
    local_f : (K_fg,) array of float
        Local Fg frequency vector (sums to 1).
    class_i_mask : (K_fg,) array of bool
        True where the Fg is Class I.
    """
    mother_fg_idx = np.asarray(mother_fg_idx).astype(int)
    is_i = class_i_mask[mother_fg_idx]
    p_I = float(local_f[class_i_mask].sum())

    if is_i.any():
        # Case A: mother expresses only her Class I alleles
        m_expr = np.unique(mother_fg_idx[is_i])
        p_M = float(local_f[m_expr].sum())
        return float(max(0.0, min(1.0, (1.0 - p_M) ** 4)))

    # Case B: mother has 0 Class I → expresses all her Class II alleles
    m_expr = np.unique(mother_fg_idx)
    p_M = float(local_f[m_expr].sum())
    # P(father has ≥1 Class I) + P(father has 0 Class I AND none in m_expr)
    return float(max(0.0, min(1.0, (1.0 - (1.0 - p_I) ** 4)
                              + max(0.0, 1.0 - p_I - p_M) ** 4)))


def p_compat_sporophytic_batch(mother_genotypes: np.ndarray,
                                local_f: np.ndarray,
                                class_i_mask: np.ndarray) -> np.ndarray:
    """Vectorised sporophytic P_compat over many mothers.

    mother_genotypes : (M, 4) array of Fg indices.
    Returns (M,) array of P_compat values.
    """
    M = mother_genotypes.shape[0]
    p_I = float(local_f[class_i_mask].sum())
    out = np.empty(M, dtype=float)
    # Broadcast approach: per-mother processing (M is small in practice).
    is_i_all = class_i_mask[mother_genotypes]           # (M, 4) bool
    has_i = is_i_all.any(axis=1)                        # (M,) bool

    # ---- Case A: mothers with any Class I ----
    if has_i.any():
        idx_a = np.where(has_i)[0]
        for i in idx_a:
            m_expr = np.unique(mother_genotypes[i][is_i_all[i]])
            p_M = float(local_f[m_expr].sum())
            out[i] = max(0.0, min(1.0, (1.0 - p_M) ** 4))

    # ---- Case B: mothers with 0 Class I ----
    if (~has_i).any():
        idx_b = np.where(~has_i)[0]
        for i in idx_b:
            m_expr = np.unique(mother_genotypes[i])
            p_M = float(local_f[m_expr].sum())
            out[i] = max(0.0, min(1.0,
                (1.0 - (1.0 - p_I) ** 4) + max(0.0, 1.0 - p_I - p_M) ** 4))
    return out


def species_mean_p_compat(prior_f: np.ndarray,
                           class_i_mask: np.ndarray,
                           n_mothers: int = 20_000,
                           rng: np.random.Generator | None = None) -> float:
    """Species-mean P_compat under the P1 prior: expected sporophytic
    P_compat for a random mother drawn from the species-wide
    distribution. Used to recalibrate the traffic-light bands.
    NOTE: this is the naive independent-draws model. Under empirical
    LEPA zygosity (~66 % of plants are single-identity homozygotes),
    the species mean is substantially higher — use
    `species_mean_p_compat_empirical` instead."""
    rng = rng or np.random.default_rng(2029)
    K_fg = len(prior_f)
    mothers = rng.choice(K_fg, size=(n_mothers, 4), p=prior_f)
    return float(p_compat_sporophytic_batch(mothers, prior_f, class_i_mask).mean())


# ---------------------------------------------------------------------------
# Empirical LEPA zygosity — Canu-amplicon Step 23 observation that
# ~66 % of individuals carry a single distinct functional SRK identity
# (AAAA-like), 32 % carry two, 2 % carry three. This is much higher than
# any independent-tetraploid-draw model would predict (~5 % under P1).
# ---------------------------------------------------------------------------


def load_zygosity_dist(tsv_path: Path = DEFAULT_ZYGOSITY_TSV) -> np.ndarray:
    """Return array [p_1, p_2, p_3, p_4] of empirical probabilities that
    a random LEPA plant carries 1, 2, 3, or 4 distinct functional SRK
    identities. Normalised to sum to 1."""
    if not tsv_path.exists():
        raise FileNotFoundError(f"Missing zygosity TSV at {tsv_path}")
    df = pd.read_csv(tsv_path, sep="\t", encoding="utf-8-sig")
    lookup = dict(zip(
        df["n_distinct_functional_alleles"].astype(int),
        df["fraction"].astype(float)))
    probs = np.array([lookup.get(k, 0.0) for k in (1, 2, 3, 4)],
                     dtype=float)
    s = probs.sum()
    if s <= 0:
        raise ValueError(f"Zygosity TSV at {tsv_path} sums to 0")
    return probs / s


def sample_genotypes_empirical(n: int,
                                local_f: np.ndarray,
                                zygosity_probs: np.ndarray,
                                rng: np.random.Generator) -> np.ndarray:
    """Sample n tetraploid genotypes from the local Fg frequency vector
    `local_f`, using the empirical LEPA zygosity distribution
    `zygosity_probs`.

    Each genotype is returned as a length-4 array of Fg indices, with
    duplicates padding out homozygous plants (a 1-distinct plant returns
    [A, A, A, A]; a 2-distinct plant returns [A, B, ?, ?] with the last
    two slots filled with the first allele to keep shape consistent).
    Downstream SI functions only inspect the DISTINCT set within each
    row, so the padding does not affect any calculation.

    Vectorised per n_distinct category; only the (rare) 2- and 3-distinct
    samples fall back to a per-row loop for without-replacement draws.
    """
    available_idx = np.where(local_f > 0)[0]
    if len(available_idx) == 0:
        raise ValueError("local_f has no non-zero entries")
    available_f = local_f[available_idx]
    available_f = available_f / available_f.sum()
    max_n = min(4, len(available_idx))
    probs_clipped = zygosity_probs[:max_n].copy()
    if probs_clipped.sum() <= 0:
        probs_clipped = np.ones(max_n)
    probs_clipped = probs_clipped / probs_clipped.sum()

    n_distinct_draws = rng.choice(
        np.arange(1, max_n + 1), size=n, p=probs_clipped)
    genotypes = np.empty((n, 4), dtype=int)
    for k in range(1, max_n + 1):
        mask = (n_distinct_draws == k)
        n_k = int(mask.sum())
        if n_k == 0:
            continue
        if k == 1:
            positions = rng.choice(
                len(available_idx), size=n_k, p=available_f)
            distinct = available_idx[positions][:, None]        # (n_k, 1)
        else:
            distinct_positions = np.empty((n_k, k), dtype=int)
            for i in range(n_k):
                distinct_positions[i] = rng.choice(
                    len(available_idx), size=k, replace=False, p=available_f)
            distinct = available_idx[distinct_positions]        # (n_k, k)
        if k < 4:
            pad_col = distinct[:, 0:1]
            padded = np.concatenate(
                [distinct, np.tile(pad_col, (1, 4 - k))], axis=1)
        else:
            padded = distinct
        genotypes[mask] = padded
    return genotypes


def p_compat_sporophytic_empirical(
        mother_genotypes: np.ndarray,
        local_f: np.ndarray,
        class_i_mask: np.ndarray,
        zygosity_probs: np.ndarray,
        n_fathers: int = 2_000,
        rng: np.random.Generator | None = None) -> np.ndarray:
    """Monte-Carlo P_compat for each mother against a population of
    candidate tetraploid fathers drawn from `local_f` under the
    empirical LEPA zygosity distribution.

    For each of the M mothers, sample `n_fathers` fathers with the same
    empirical zygosity structure that produced the mothers themselves;
    for each mother, return the fraction of fathers whose expressed
    SRK identity set does NOT overlap the mother's expressed set.

    This replaces the naive analytical formulas (which assumed
    independent-tetraploid-draw fathers with ~5 % homozygosity), giving
    a symmetric mother/father empirical-zygosity model."""
    rng = rng or np.random.default_rng(2029)
    M = mother_genotypes.shape[0]
    K_fg = len(local_f)

    # Sample fathers ONCE per call — used across all mothers for
    # efficiency; each mother sees the same simulated father sub-population.
    fathers = sample_genotypes_empirical(n_fathers, local_f, zygosity_probs, rng)
    father_class_i = class_i_mask[fathers]                    # (F, 4)
    father_has_i   = father_class_i.any(axis=1)               # (F,)
    # A father allele participates in SI iff:
    #   father has Class I AND allele is Class I, OR
    #   father has no Class I (Case B) — then ALL alleles participate.
    father_include = father_class_i | ~father_has_i[:, None]  # (F, 4)

    results = np.empty(M, dtype=float)
    for i in range(M):
        m = mother_genotypes[i]
        m_is_i = class_i_mask[m]
        if m_is_i.any():
            m_expr = np.unique(m[m_is_i])
        else:
            m_expr = np.unique(m)
        m_mask = np.zeros(K_fg, dtype=bool)
        m_mask[m_expr] = True
        # For each father, is any of his EXPRESSED alleles in M?
        father_in_M = m_mask[fathers]                          # (F, 4)
        overlap_positions = father_in_M & father_include       # (F, 4)
        compat = ~overlap_positions.any(axis=1)                # (F,)
        results[i] = float(compat.mean())
    return results


def species_mean_p_compat_empirical(
        prior_f: np.ndarray,
        class_i_mask: np.ndarray,
        zygosity_probs: np.ndarray | None = None,
        n_mothers: int = 20_000,
        n_fathers: int = 2_000,
        rng: np.random.Generator | None = None) -> float:
    """Species-mean P_compat under the P1 prior + empirical LEPA
    zygosity structure. This is the number to calibrate the § A.7
    traffic-light bands against."""
    rng = rng or np.random.default_rng(2029)
    if zygosity_probs is None:
        zygosity_probs = load_zygosity_dist()
    mothers = sample_genotypes_empirical(n_mothers, prior_f, zygosity_probs, rng)
    return float(p_compat_sporophytic_empirical(
        mothers, prior_f, class_i_mask, zygosity_probs, n_fathers, rng).mean())


def traffic_light_bands(species_mean: float,
                         failed_frac: float = 1.0 / 3.0,
                         struggling_frac: float = 2.0 / 3.0) -> dict[str, float]:
    """Recalibrate the failed / struggling / sustainable band boundaries
    proportional to the species-mean P_compat. Defaults are 1/3 and
    2/3 of the species mean — semantically clean and close to the
    diploid gametophytic bands (0.20 / 0.40, which were 0.32 and 0.63
    of the diploid species mean 0.63)."""
    return {
        "failed_max":      failed_frac * species_mean,
        "struggling_max":  struggling_frac * species_mean,
        "species_mean":    species_mean,
    }
