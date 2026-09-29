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
    distribution. Used to recalibrate the traffic-light bands."""
    rng = rng or np.random.default_rng(2029)
    K_fg = len(prior_f)
    mothers = rng.choice(K_fg, size=(n_mothers, 4), p=prior_f)
    return float(p_compat_sporophytic_batch(mothers, prior_f, class_i_mask).mean())


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
