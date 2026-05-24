"""Shared utility functions for the polyploid crossing model.

This module consolidates all functions previously duplicated across notebooks
01-06. Import with: `from polyploid_utils import *` (notebooks add ../src to
sys.path) or `from src.polyploid_utils import *` (scripts at project root).
"""

import itertools
from collections import Counter

import numpy as np
from scipy.optimize import minimize

__all__ = [
    # Core genotype helpers
    "canonical",
    "enumerate_genotypes",
    "allele_frequencies",
    # Crossing mechanics
    "form_gametes",
    "is_compatible",
    "gamete_acceptance_prob",
    "cross",
    "crossing_compatibility",
    "sample_offspring",
    "build_crossing_matrix",
    # Equilibrium analysis
    "target_frequencies",
    "distance_from_equilibrium",
    "evenness_J",
    "p_compat",
    # Genotype quality (GFS / TP2)
    "gfs",
    "genotype_class",
    "expected_offspring_gfs",
    "prop_AAAA",
    "mean_gfs",
    # Allele preservation
    "identify_rare_alleles",
    "get_mandatory_rare_crosses",
    "select_elites",
    # Simulation
    "simulate_generation",
    # Demography
    "logistic_n_offspring",
    "effective_population_size",
    "ne_harmonic_mean",
    # Optimization
    "enumerate_compatible_crosses",
    "compute_greedy_weights",
    "compute_optimal_weights",
]

# Mapping from genotype class label to ordinal tier rank, for thresholds
# expressed as "partner must be tier >= AABB" etc.
_TIER_RANK = {"AAAA": 0, "AAAB": 1, "AABB": 2, "AABC": 3, "ABCD": 4}

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

# Sven's empirical SI leakage rate, estimated from observed AAAA proportions
# (L_hat = prop_AAAA / 3.5 per group). This is a PAIR-LEVEL rate -- the
# fraction of strict-incompatible plant pairs that effectively produce seed.
# When passing to functions that use PER-GAMETE leakage (cross, simulate_*),
# convert via L_gamete ~ 1 - (1 - L_pair)^(1/n_gametes). For n=6 gametes
# (tetrasomic), L_pair = 0.18 corresponds to L_gamete ~ 0.033.
L_HAT_PAIR = 0.18      # Sven's empirical pair-level leakage
L_HAT_GAMETE = 0.033   # Approximate per-gamete equivalent for tetrasomic 6-gamete

# ---------------------------------------------------------------------------
# Core genotype helpers
# ---------------------------------------------------------------------------

def canonical(alleles):
    """Return the canonical (sorted tuple) form of a genotype."""
    return tuple(sorted(alleles))


def enumerate_genotypes(allele_pool, ploidy=4):
    """Enumerate all unique multiset genotypes for a given allele pool and ploidy.

    Uses combinations-with-replacement: C(n+k-1, k) genotypes total.
    """
    return list(itertools.combinations_with_replacement(sorted(allele_pool), ploidy))


def allele_frequencies(population, allele_pool=None):
    """Compute allele frequencies across a population.

    Counts every allele copy (4 per tetraploid) and normalizes to proportions.
    If allele_pool is provided, includes zero entries for absent alleles.
    """
    counts = Counter()
    total = 0
    for genotype in population:
        for allele in genotype:
            counts[allele] += 1
            total += 1
    freqs = {a: counts[a] / total for a in counts}
    if allele_pool is not None:
        for a in allele_pool:
            if a not in freqs:
                freqs[a] = 0.0
    return dict(sorted(freqs.items()))


# ---------------------------------------------------------------------------
# Crossing mechanics
# ---------------------------------------------------------------------------

def form_gametes(genotype, inheritance_mode="tetrasomic"):
    """Return diploid gametes from a tetraploid genotype with probabilities.

    Parameters
    ----------
    genotype : tuple
        4-tuple of allele IDs.
    inheritance_mode : {'tetrasomic', 'disomic-averaged'}
        - ``tetrasomic``: random pairing across all 4 chromosomes; all
          C(4,2)=6 ordered gametes equally probable (1/6 each). The
          default; matches the assumption in the upstream LEPA report.
        - ``disomic-averaged``: two homeologous subgenomes segregate
          independently (allopolyploid model). The gamete distribution
          is averaged uniformly across distinct subgenome partitions of
          the 4 alleles, since per-individual subgenome assignment is
          unobservable from the current genotyping data. See Q1 in
          ``Questions_for_Sven.md``.

    Returns
    -------
    list[tuple[tuple, float]]
        List of ``(gamete, probability)`` pairs. Probabilities sum to 1.
        For tetrasomic, returns 6 pairs each with probability 1/6.
        For disomic-averaged, the number of distinct gametes varies by
        genotype (e.g. AAAA -> 1; AAAB -> 2; AABB -> 3; AABC -> 4; ABCD -> 6),
        with probabilities reflecting the averaged disomic distribution.
    """
    if inheritance_mode == "tetrasomic":
        pairs = list(itertools.combinations(genotype, 2))
        per_pair = 1.0 / len(pairs)
        weights = Counter()
        for g in pairs:
            weights[tuple(sorted(g))] += per_pair
        return sorted(weights.items())
    if inheritance_mode == "disomic-averaged":
        return _disomic_averaged_gametes(genotype)
    raise ValueError(
        f"Unknown inheritance_mode {inheritance_mode!r}; "
        f"expected 'tetrasomic' or 'disomic-averaged'."
    )


def _disomic_averaged_gametes(genotype):
    """Disomic-averaged gametes for a tetraploid genotype.

    Enumerates all distinct partitions of the 4 alleles into two pairs
    (the homeologous subgenomes) and averages the gamete distribution
    uniformly across partitions. For each partition (S1, S2), gametes
    are formed by picking one allele from each subgenome -- 4 gametes
    of equal weight per partition.

    See ``Questions_for_Sven.md`` Q2 for why uniform-prior averaging is
    used (per-individual subgenome assignment is unobservable from the
    current genotyping pipeline).
    """
    g = list(genotype)
    if len(g) != 4:
        raise ValueError("Disomic-averaged inheritance requires a tetraploid genotype (4 alleles)")
    seen_partitions = set()
    partitions = []
    for indices in itertools.combinations(range(4), 2):
        s1 = tuple(sorted(g[i] for i in indices))
        s2 = tuple(sorted(g[i] for i in range(4) if i not in indices))
        # Canonical partition key: sort the two subgenome tuples so
        # (S1, S2) and (S2, S1) hash identically.
        key = tuple(sorted([s1, s2]))
        if key not in seen_partitions:
            seen_partitions.add(key)
            partitions.append(key)

    gamete_weights = Counter()
    per_partition_weight = 1.0 / len(partitions)
    for s1, s2 in partitions:
        for a1 in s1:
            for a2 in s2:
                gamete = tuple(sorted([a1, a2]))
                # 4 gametes per partition, equal probability within partition
                gamete_weights[gamete] += per_partition_weight * 0.25
    return sorted(gamete_weights.items())


def is_compatible(maternal_genotype, pollen_gamete):
    """Strict SI compatibility check: True if pollen shares no alleles with maternal plant.

    Deterministic boolean. For leaky-SI semantics, use
    :func:`gamete_acceptance_prob` (per-gamete) or pass ``leakage`` to
    :func:`cross` / :func:`crossing_compatibility`.
    """
    maternal_alleles = set(maternal_genotype)
    return not any(a in maternal_alleles for a in pollen_gamete)


def gamete_acceptance_prob(maternal_genotype, pollen_gamete, leakage=0.0):
    """Probability a paternal gamete is accepted by the maternal stigma.

    Parameters
    ----------
    maternal_genotype : tuple
        4-tuple of maternal allele IDs.
    pollen_gamete : tuple
        2-tuple of paternal allele IDs (one diploid gamete).
    leakage : float in [0, 1]
        Per-gamete SI escape probability for strict-rejected gametes.
        Aggregate of empirical mechanisms (Lewis 1947 competitive
        interaction, partial LoF, recognition errors). Default 0
        reproduces strict SI behaviour.

        NOTE: This is a PER-GAMETE rate, not Sven's pair-level L_hat.
        For pair-level matching at L_hat = 0.18 (tetrasomic, 6 gametes),
        use ``leakage = L_HAT_GAMETE`` (~ 0.033). See module docstring.

    Returns
    -------
    float
        1.0 if strictly compatible, ``leakage`` if strictly rejected.
    """
    if is_compatible(maternal_genotype, pollen_gamete):
        return 1.0
    return float(leakage)


def cross(parent_a, parent_b, inheritance_mode="tetrasomic", leakage=0.0):
    """Compute offspring genotype distribution for a directed cross.

    Returns a dict mapping offspring genotype -> probability (sums to 1).
    Returns an empty dict only when *no* paternal gamete has any
    acceptance probability (under strict SI, this happens when all
    gametes share alleles with the mother; under leakage > 0, every pair
    has at least some offspring probability).

    The offspring distribution is conditional on offspring being produced
    -- use :func:`crossing_compatibility` to obtain the unconditional
    acceptance rate for the pair.

    Parameters
    ----------
    parent_a : tuple
        Maternal genotype.
    parent_b : tuple
        Paternal genotype.
    inheritance_mode : {'tetrasomic', 'disomic-averaged'}
        See :func:`form_gametes`.
    leakage : float in [0, 1]
        Per-gamete SI escape probability. See :func:`gamete_acceptance_prob`.
    """
    maternal_gametes = form_gametes(parent_a, inheritance_mode)
    paternal_gametes = form_gametes(parent_b, inheritance_mode)

    offspring_weights = Counter()
    total_weight = 0.0
    mat_set = set(parent_a)
    for mg, mp in maternal_gametes:
        for pg, pp in paternal_gametes:
            # Strict-SI check inlined for speed; equivalent to
            # gamete_acceptance_prob(parent_a, pg, leakage).
            if any(a in mat_set for a in pg):
                accept = leakage
            else:
                accept = 1.0
            if accept <= 0:
                continue
            w = mp * pp * accept
            offspring_weights[canonical(mg + pg)] += w
            total_weight += w

    if total_weight <= 0:
        return {}
    return {g: w / total_weight for g, w in sorted(offspring_weights.items())}


def crossing_compatibility(parent_a, parent_b, inheritance_mode="tetrasomic", leakage=0.0):
    """Expected fraction of paternal gametes that yield viable offspring.

    Returns ``Σ_g pg_prob * acceptance(g, leakage)``. Under strict SI
    (``leakage=0``) this equals the probability of any paternal gamete
    being SI-compatible. Under ``leakage > 0`` it equals
    ``P_strict + leakage * (1 - P_strict)`` -- the per-gamete leakage
    interpretation of Sven's pair-level L_hat formula.

    Parameters
    ----------
    parent_a : tuple
        Maternal genotype.
    parent_b : tuple
        Paternal genotype.
    inheritance_mode : {'tetrasomic', 'disomic-averaged'}
    leakage : float in [0, 1]
        Per-gamete escape probability for strict-rejected gametes.
    """
    paternal_gametes = form_gametes(parent_b, inheritance_mode)
    mat_set = set(parent_a)
    total = 0.0
    for pg, pp in paternal_gametes:
        if any(a in mat_set for a in pg):
            total += pp * leakage
        else:
            total += pp
    return total


def sample_offspring(parent_a, parent_b, inheritance_mode="tetrasomic", leakage=0.0):
    """Sample a single offspring from a cross (None if no offspring possible).

    Convenience wrapper around :func:`cross` plus a stochastic draw.
    """
    offspring_dist = cross(parent_a, parent_b, inheritance_mode, leakage)
    if not offspring_dist:
        return None
    genotypes = list(offspring_dist.keys())
    probs = list(offspring_dist.values())
    idx = np.random.choice(len(genotypes), p=probs)
    return genotypes[idx]


def build_crossing_matrix(genotypes, inheritance_mode="tetrasomic", leakage=0.0):
    """Build compatibility and outcome matrices for all ordered genotype pairs.

    Returns ``(compatibility_dict, outcomes_dict)`` keyed by
    ``(genotype_a, genotype_b)``. Self-pairs (i == j) get compatibility
    0.0 and an empty outcomes dict by convention.
    """
    compatibility = {}
    outcomes = {}
    for i, ga in enumerate(genotypes):
        for j, gb in enumerate(genotypes):
            if i == j:
                compatibility[(ga, gb)] = 0.0
                outcomes[(ga, gb)] = {}
            else:
                compat = crossing_compatibility(ga, gb, inheritance_mode, leakage)
                compatibility[(ga, gb)] = compat
                outcomes[(ga, gb)] = cross(ga, gb, inheritance_mode, leakage)
    return compatibility, outcomes


# ---------------------------------------------------------------------------
# Equilibrium analysis
# ---------------------------------------------------------------------------

def target_frequencies(allele_pool):
    """Return the uniform (equilibrium) frequency distribution: 1/n for each allele."""
    n = len(allele_pool)
    return {a: 1.0 / n for a in sorted(allele_pool)}


def distance_from_equilibrium(population, allele_pool):
    """Compute distance metrics from equal-frequency equilibrium.

    Returns dict with: variance, chi_squared, kl_divergence, extinct_alleles,
    endangered_alleles.
    """
    freqs = allele_frequencies(population, allele_pool)
    n = len(allele_pool)
    target = 1.0 / n
    freq_values = np.array([freqs[a] for a in sorted(allele_pool)])
    variance = float(np.var(freq_values))
    chi_squared = float(np.sum((freq_values - target) ** 2 / target))
    eps = 1e-12
    freq_safe = np.maximum(freq_values, eps)
    kl_div = float(np.sum(freq_safe * np.log(freq_safe / target)))
    extinct = int(np.sum(freq_values == 0))
    endangered = int(np.sum((freq_values > 0) & (freq_values < 0.02)))
    return {
        "variance": variance,
        "chi_squared": chi_squared,
        "kl_divergence": kl_div,
        "extinct_alleles": extinct,
        "endangered_alleles": endangered,
    }


def evenness_J(population, allele_pool=None):
    """Shannon evenness J = H / ln(k) of allele frequencies in a population.

    Range: 0 (one allele dominates) -> 1 (all observed alleles equal frequency).
    J = 1 at NFDS equilibrium. This is the TP1 x-axis from the LEPA report Q4.

    Parameters
    ----------
    population : list[tuple]
        Population genotypes (each a 4-tuple of allele IDs).
    allele_pool : iterable or None
        Optional explicit allele pool. Unused for J (uses only observed
        alleles by definition), kept for API symmetry with allele_frequencies.

    Returns
    -------
    float
        Shannon evenness in [0, 1]. Returns 1.0 if k == 1 (single allele
        observed) by convention -- though k == 1 is degenerate biologically.

    Notes
    -----
    Matches the J metric in Sven's SRK_TP1_compatibility_metrics.py upstream.
    Uses observed alleles only (k = number of alleles with freq > 0), not
    the size of the allele pool. This is a property of the *observed*
    population, not of the (possibly empty) total pool.
    """
    freqs = allele_frequencies(population, allele_pool)
    p = np.array([f for f in freqs.values() if f > 0])
    k = len(p)
    if k <= 1:
        return 1.0
    shannon_h = float(-np.sum(p * np.log(p)))
    return shannon_h / float(np.log(k))


def p_compat(population, leakage=0.0, inheritance_mode="tetrasomic"):
    """Fraction of ordered (maternal, paternal) pairs that are SI-compatible.

    P_compat is the TP1 y-axis from the LEPA report Q4. A value of 0.40
    means ~40% of random plant pairs in the population can produce seed
    under tetraploid sporophytic SI. The strict-SI compatibility test
    matches :func:`is_compatible`: a pair is compatible if at least one
    paternal gamete shares no allele with the maternal genotype. Both
    ``i -> j`` and ``j -> i`` directed pairs are counted; self-pairs are
    excluded.

    Parameters
    ----------
    population : list[tuple]
        Population genotypes (each a 4-tuple of allele IDs).
    leakage : float in [0, 1]
        PAIR-LEVEL leakage (matches Sven's L_hat semantics): a fraction L
        of strict-incompatible pairs are treated as compatible. Default 0
        gives strict SI. Sven's L_hat ~ 0.18 reproduces empirical AAAA
        proportions.

        Transformation matches ``SRK_TP1_compatibility_metrics.py``::

            P_compat(L) = P_strict + L * (1 - P_strict)

        Note: This is *pair-level* leakage to match Sven's metric
        directly. :func:`cross` and :func:`crossing_compatibility` use
        *per-gamete* leakage instead -- a different semantic with the
        same parameter name. See module docstring for the conversion.
    inheritance_mode : {'tetrasomic', 'disomic-averaged'}
        Determines the gamete distribution used for the strict
        compatibility check. Default tetrasomic enumerates all C(4,2)=6
        gametes uniformly; disomic-averaged uses :func:`form_gametes`.

    Returns
    -------
    float
        Fraction in [0, 1]. Returns 0.0 if population has < 2 individuals.
    """
    n = len(population)
    if n < 2:
        return 0.0
    compat_strict = 0
    n_pairs = 0
    if inheritance_mode == "tetrasomic":
        # Fast path: uniform tetrasomic gamete enumeration
        for i, mat in enumerate(population):
            mat_set = set(mat)
            for j, pat in enumerate(population):
                if i == j:
                    continue
                n_pairs += 1
                for g in itertools.combinations(pat, 2):
                    if not any(a in mat_set for a in g):
                        compat_strict += 1
                        break
    else:
        for i, mat in enumerate(population):
            mat_set = set(mat)
            for j, pat in enumerate(population):
                if i == j:
                    continue
                n_pairs += 1
                for g, _ in form_gametes(pat, inheritance_mode):
                    if not any(a in mat_set for a in g):
                        compat_strict += 1
                        break
    p_strict = compat_strict / n_pairs if n_pairs > 0 else 0.0
    if leakage <= 0:
        return p_strict
    leakage = float(min(1.0, max(0.0, leakage)))
    return p_strict + leakage * (1.0 - p_strict)


# ---------------------------------------------------------------------------
# Genotype quality (GFS / TP2)
# ---------------------------------------------------------------------------

def gfs(genotype):
    """Genotypic Fitness Score: proportion of heterozygous diploid gametes.

    ``GFS = 1 - Σ_k n_k(n_k-1) / 12`` where n_k is the copy count of allele
    k in the individual. Equivalent to P(gamete is heterozygous) under
    tetrasomic inheritance. Matches the canonical formula in upstream
    ``Scripts/SRK_individual_GFS.R``.

    Tier values:

    +----------+--------------+-------+
    | Genotype | Copy counts  |  GFS  |
    +==========+==============+=======+
    | ABCD     | (1, 1, 1, 1) | 1.000 |
    +----------+--------------+-------+
    | AABC     | (2, 1, 1)    | 0.833 |
    +----------+--------------+-------+
    | AABB     | (2, 2)       | 0.667 |
    +----------+--------------+-------+
    | AAAB     | (3, 1)       | 0.500 |
    +----------+--------------+-------+
    | AAAA     | (4,)         | 0.000 |
    +----------+--------------+-------+

    Note: under disomic inheritance, GFS as defined above is not the
    correct fitness signal -- the proportion of heterozygous gametes
    depends on subgenome assignment. See ``Questions_for_Sven.md`` Q1.
    """
    counts = Counter(genotype).values()
    return 1.0 - sum(c * (c - 1) for c in counts) / 12.0


def genotype_class(genotype):
    """Return the tetraploid genotype tier label (AAAA / AAAB / AABB / AABC / ABCD).

    Raises ``KeyError`` for invalid (non-tetraploid) copy-count patterns.
    """
    counts = tuple(sorted(Counter(genotype).values(), reverse=True))
    return {
        (4,):           "AAAA",
        (3, 1):         "AAAB",
        (2, 2):         "AABB",
        (2, 1, 1):      "AABC",
        (1, 1, 1, 1):   "ABCD",
    }[counts]


def expected_offspring_gfs(parent_a, parent_b, inheritance_mode="tetrasomic", leakage=0.0):
    """Expected GFS of offspring from a directed cross.

    Computes the GFS of each possible offspring genotype weighted by its
    probability in :func:`cross`. Returns 0.0 if no offspring possible
    (e.g., strict-SI incompatible pair with leakage=0).

    Used by :func:`compute_optimal_weights` as the genotype-quality term
    in the GFS-aware objective, per Crossing Plan section 3 Phase 2 step 5.
    """
    dist = cross(parent_a, parent_b, inheritance_mode, leakage)
    if not dist:
        return 0.0
    return sum(prob * gfs(genotype) for genotype, prob in dist.items())


def prop_AAAA(population):
    """Fraction of individuals in the AAAA genotype tier (the C3 endpoint).

    Used in TP2: a group with prop_AAAA > 0.30 AND mean_GFS < 0.667 is
    flagged CRITICAL per the LEPA report Q5.
    """
    if not population:
        return 0.0
    aaaa = sum(1 for g in population if genotype_class(g) == "AAAA")
    return aaaa / len(population)


def mean_gfs(population):
    """Population mean GFS -- the TP2 x-axis from the LEPA report Q5."""
    if not population:
        return 0.0
    return sum(gfs(g) for g in population) / len(population)


# ---------------------------------------------------------------------------
# Allele preservation
# ---------------------------------------------------------------------------

def identify_rare_alleles(population, allele_pool, threshold=None, max_carriers=2):
    """Return set of alleles considered rare in this population.

    An allele is rare if EITHER condition holds:
      - It is observed in ``<= max_carriers`` individuals (default 2);
      - Its allele frequency is strictly below ``threshold`` (default None,
        meaning the frequency criterion is disabled).

    Default ``max_carriers=2, threshold=None`` is the recalibrated rule
    for the post-synonymy 27-biological-allele world: any allele in 0 < c <= 2
    carriers is endangered. The legacy ``threshold=0.05`` criterion is left
    available for backwards compatibility with the archived 94-protein-level
    notebooks -- but at 27 alleles, target_freq=1/27=0.037, so a 0.05
    threshold flags every allele below uniform as rare and breaks the
    optimizer. See ``Critical_Review.md`` finding B.6.
    """
    rare = set()
    if max_carriers is not None:
        carrier_counts = {a: 0 for a in allele_pool}
        for g in population:
            for a in set(g):
                if a in carrier_counts:
                    carrier_counts[a] += 1
        for a, c in carrier_counts.items():
            if 0 < c <= max_carriers:
                rare.add(a)
    if threshold is not None:
        freqs = allele_frequencies(population, allele_pool)
        for a, f in freqs.items():
            if 0 < f < threshold:
                rare.add(a)
    return rare


def get_mandatory_rare_crosses(population, allele_pool, threshold=0.05,
                               already_covered=None, min_partner_tier="AABB",
                               inheritance_mode="tetrasomic", leakage=0.0):
    """Find one compatible cross per endangered allele not already covered.

    An allele is endangered if 1 or 2 individuals carry it. For each
    endangered allele, search over (carrier, partner) pairs for the cross
    that maximises the probability of offspring carrying the target allele.

    Per Crossing Plan section 3 Phase 2 step 7, the partner is restricted
    to tier >= ``min_partner_tier`` (default "AABB"). AAAA-carrying rare
    alleles can still be USED -- but only as the maternal parent with an
    AABB+ donor -- which converts an AAAA dead-end lineage into a higher-GFS
    offspring lineage in one generation.

    Parameters
    ----------
    population : list[tuple]
        Population genotypes.
    allele_pool : iterable
        Full allele pool.
    threshold : float
        Legacy frequency threshold (currently unused; kept for backwards
        compatibility with older callers). The endangerment criterion is
        carrier-count <=2.
    already_covered : set or None
        Alleles already covered by elitism; skip them.
    min_partner_tier : {"AAAA", "AAAB", "AABB", "AABC", "ABCD"}
        Minimum genotype tier required of the partner (non-carrier side
        of the cross). Default "AABB" follows the Plan.
    inheritance_mode, leakage :
        Forwarded to :func:`cross` when scoring candidate pairs.

    Returns
    -------
    list[tuple[int, int]]
        ``(maternal_idx, paternal_idx)`` pairs.
    """
    carrier_counts = {a: 0 for a in allele_pool}
    for g in population:
        for a in set(g):
            if a in carrier_counts:
                carrier_counts[a] += 1
    endangered = {a for a, c in carrier_counts.items() if 0 < c <= 2}

    if already_covered:
        endangered -= already_covered

    if not endangered:
        return []

    min_rank = _TIER_RANK[min_partner_tier]

    mandatory = []
    covered = set()
    for allele in sorted(endangered):
        if allele in covered:
            continue
        carriers = [i for i, g in enumerate(population) if allele in g]
        best_cross, best_score = None, -1
        for ci in carriers:
            for j in range(len(population)):
                if j == ci:
                    continue
                # Per Plan section 3 step 7: partner must be tier >= min_partner_tier
                if _TIER_RANK[genotype_class(population[j])] < min_rank:
                    continue
                od = cross(population[ci], population[j], inheritance_mode, leakage)
                if not od:
                    continue
                score = sum(p for g, p in od.items() if allele in g)
                if score > best_score:
                    best_score, best_cross = score, (ci, j)
        if best_cross is not None:
            mandatory.append(best_cross)
            # Mark all endangered alleles this cross can produce
            od = cross(population[best_cross[0]], population[best_cross[1]],
                       inheritance_mode, leakage)
            for g, _ in od.items():
                for a in g:
                    if a in endangered:
                        covered.add(a)
    return mandatory


def select_elites(population, allele_pool, elite_frac=0.1, gfs_filter=True):
    """Select elite individuals carrying the rarest alleles, gated by GFS.

    Each individual's elite score is::

        score = quality_weight * sum(1 / freq(allele) for allele in genotype)

    where ``quality_weight`` is determined by the genotype tier per
    Crossing Plan section 3 Phase 2 step 6:

    +-----------+----------------+
    | Tier      | quality_weight |
    +===========+================+
    | AAAA      | 0.0  (excluded)|
    +-----------+----------------+
    | AAAB      | 0.5            |
    +-----------+----------------+
    | AABB/AABC/ABCD | 1.0       |
    +-----------+----------------+

    The discrete weighting (rather than raw GFS) follows the Plan literally.
    Rationale: an AAAA carrier of a rare allele scores ``4/freq`` -- the
    highest possible -- under the legacy rule, but is the *worst* possible
    carrier (only AA pollen, rejected wherever A is in the mother).
    Mandatory rare-allele crosses (see :func:`get_mandatory_rare_crosses`)
    cover any rare alleles that the GFS filter excludes from elitism.

    Set ``gfs_filter=False`` to recover legacy behaviour (no GFS gating).

    Number of elites is ``max(1, elite_frac * n)``.
    """
    freqs = allele_frequencies(population, allele_pool)
    n_elite = max(1, int(len(population) * elite_frac))
    scores = []
    for i, g in enumerate(population):
        rarity = sum(1.0 / max(freqs.get(a, 1e-12), 1e-12) for a in g)
        if gfs_filter:
            cls = genotype_class(g)
            quality = {"AAAA": 0.0, "AAAB": 0.5,
                       "AABB": 1.0, "AABC": 1.0, "ABCD": 1.0}[cls]
        else:
            quality = 1.0
        scores.append((rarity * quality, i))
    scores.sort(reverse=True)
    return [idx for s, idx in scores[:n_elite] if s > 0]


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------

def simulate_generation(population, n_offspring=None, crossing_plan=None,
                        allele_pool=None, preserve_rare=False,
                        rare_threshold=0.05, elite_frac=0.1,
                        inheritance_mode="tetrasomic", leakage=0.0,
                        safety_net=False,
                        gfs_filter=True, min_partner_tier="AABB"):
    """Simulate one generation with optional allele preservation.

    Parameters
    ----------
    population : list[tuple]
        Parental population genotypes.
    n_offspring : int or None
        Target offspring count. Defaults to len(population).
    crossing_plan : list[(maternal_idx, paternal_idx, weight)] or None
        Weighted crossing recommendations. If None, random mating fills
        the remaining slots.
    allele_pool : iterable or None
        Required when ``preserve_rare=True``; full allele pool used for
        elitism and mandatory rare-allele crosses.
    preserve_rare : bool
        Activates the preservation strategy: elitism + mandatory rare
        crosses.
    rare_threshold, elite_frac : float
        Parameters of the preservation strategy.
    inheritance_mode : {'tetrasomic', 'disomic-averaged'}
        Forwarded to :func:`sample_offspring`. See :func:`form_gametes`.
    leakage : float in [0, 1]
        Per-gamete SI escape probability. Forwarded to
        :func:`sample_offspring` for every cross attempted in this
        generation.
    safety_net : bool
        If True (and ``preserve_rare=True``), reinserts an allele carrier
        from the input population whenever an allele was lost during the
        generation. Default is False to honestly report optimizer
        performance; set True only for guaranteed-no-loss simulations.
        See ``Critical_Review.md`` finding D.5.
    """
    import random as _random

    if n_offspring is None:
        n_offspring = len(population)
    next_gen = []

    # Track which alleles were present in the input population
    present_alleles = None
    if preserve_rare and allele_pool is not None:
        present_alleles = set()
        for g in population:
            for a in g:
                if a in allele_pool:
                    present_alleles.add(a)

        # Step 1: Retain elites (10% of population, GFS-filtered by default)
        elite_indices = select_elites(population, allele_pool, elite_frac,
                                      gfs_filter=gfs_filter)
        for idx in elite_indices:
            next_gen.append(population[idx])

        # Step 2: Mandatory crosses for endangered alleles not covered by elites
        # (partner gated to tier >= min_partner_tier by default)
        elite_alleles = set()
        for idx in elite_indices:
            elite_alleles.update(population[idx])
        max_mandatory = max(1, int(n_offspring * 0.15))  # cap at 15% of population
        mandatory = get_mandatory_rare_crosses(
            population, allele_pool, rare_threshold,
            already_covered=elite_alleles,
            min_partner_tier=min_partner_tier,
            inheritance_mode=inheritance_mode,
            leakage=leakage,
        )
        for mi, pi in mandatory[:max_mandatory]:
            if len(next_gen) >= n_offspring:
                break
            child = sample_offspring(population[mi], population[pi], inheritance_mode, leakage)
            if child is not None:
                next_gen.append(child)

    # Step 3: Fill remaining slots
    if crossing_plan is not None and len(crossing_plan) > 0:
        indices = list(range(len(crossing_plan)))
        weights = np.array([w for _, _, w in crossing_plan])
        wsum = weights.sum()
        if wsum > 0:
            weights = weights / wsum
            attempts = 0
            while len(next_gen) < n_offspring and attempts < n_offspring * 20:
                idx = np.random.choice(indices, p=weights)
                mi, pi, _ = crossing_plan[idx]
                child = sample_offspring(population[mi], population[pi], inheritance_mode, leakage)
                if child is not None:
                    next_gen.append(child)
                attempts += 1
    elif crossing_plan is None:
        # Random mating
        n = len(population)
        if n >= 2:
            attempts = 0
            while len(next_gen) < n_offspring and attempts < n_offspring * 20:
                i, j = _random.sample(range(n), 2)
                child = sample_offspring(population[i], population[j], inheritance_mode, leakage)
                if child is not None:
                    next_gen.append(child)
                attempts += 1

    # Reproductive-collapse guard: if the strategy could not produce enough
    # viable offspring (e.g. EO70-style population where every pair is
    # SI-incompatible under strict SI), keep the parent population as a
    # static carryover so the simulation can continue recording metrics.
    # This is the model's empirical signal that within-population recovery
    # has failed -- the population persists without reproducing. The LEPA
    # report's BIOBANK + RESTORE designation corresponds exactly to this
    # regime; conservation requires inter-population gene flow.
    if len(next_gen) < 2:
        return list(population)

    # Step 4: Optional safety net -- verify no allele was lost.
    # Default off (safety_net=False) so strategy comparisons honestly
    # report optimizer performance; turn on for guaranteed-no-loss runs.
    if safety_net and preserve_rare and present_alleles is not None:
        next_gen_alleles = set()
        for g in next_gen:
            next_gen_alleles.update(g)
        missing = present_alleles - next_gen_alleles
        if missing:
            for allele in missing:
                carriers = [i for i, g in enumerate(population) if allele in g]
                if carriers:
                    chosen = carriers[_random.randint(0, len(carriers) - 1)]
                    next_gen.append(population[chosen])

    return next_gen


# ---------------------------------------------------------------------------
# Demography
# ---------------------------------------------------------------------------

def logistic_n_offspring(N, K, r=0.5, stochastic=True):
    """Compute next-generation population size using discrete logistic growth.

    Parameters
    ----------
    N : int
        Current census population size.
    K : int
        Carrying capacity (maximum sustainable population size).
    r : float
        Intrinsic growth rate (per-generation). Default 0.5.
    stochastic : bool
        If True, draw from Poisson(expected_N) to model demographic noise.

    Returns
    -------
    int
        Number of offspring for the next generation (minimum 2 to avoid
        extinction from stochasticity alone).
    """
    expected = N + N * r * (1 - N / K)
    expected = max(expected, 2.0)
    if stochastic:
        return max(2, int(np.random.poisson(expected)))
    return max(2, int(round(expected)))


def effective_population_size(population, allele_pool, inheritance_mode="tetrasomic", leakage=0.0):
    """Estimate effective population size (Ne) from SI compatibility structure.

    Computes Ne from the number of successful mating pairs possible under
    self-incompatibility.  Uses a heuristic Ne formula -- not a derived
    population-genetics quantity -- see ``Critical_Review.md`` finding B.7
    for caveats.

    Parameters
    ----------
    population : list of tuple
        Current population genotypes.
    allele_pool : set or list
        Complete set of allele IDs.
    inheritance_mode : {'tetrasomic', 'disomic-averaged'}
        Forwarded to :func:`crossing_compatibility`.
    leakage : float in [0, 1]
        Per-gamete SI escape probability. Under ``leakage > 0`` more
        pairs count as breeders so Ne increases.

    Returns
    -------
    dict
        Keys: 'Ne' (effective size), 'N' (census size),
        'ratio' (Ne/N), 'n_breeders' (individuals with >=1 mate),
        'mean_compatibility' (average pairwise compatibility),
        'total_compatible_pairs' (directed pair count).
    """
    N = len(population)
    if N < 2:
        return {"Ne": N, "N": N, "ratio": 1.0, "n_breeders": N,
                "mean_compatibility": 0.0, "total_compatible_pairs": 0}

    can_breed = set()
    total_compat = 0.0
    n_compatible_pairs = 0
    n_pairs = 0

    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            c = crossing_compatibility(population[i], population[j], inheritance_mode, leakage)
            n_pairs += 1
            total_compat += c
            if c > 0:
                can_breed.add(i)
                can_breed.add(j)
                n_compatible_pairs += 1

    n_breeders = len(can_breed)
    mean_compat = total_compat / n_pairs if n_pairs > 0 else 0.0

    # Ne from breeder count, discounted by mean compatibility
    # Rationale: SI restricts mate choice, reducing Ne below the breeder count
    Ne = n_breeders * mean_compat if n_breeders > 0 else 0.0
    # Floor at 1 if any breeders exist
    if n_breeders > 0:
        Ne = max(1.0, Ne)

    return {
        "Ne": Ne,
        "N": N,
        "ratio": Ne / N if N > 0 else 0.0,
        "n_breeders": n_breeders,
        "mean_compatibility": mean_compat,
        "total_compatible_pairs": n_compatible_pairs,
    }


def ne_harmonic_mean(ne_series):
    """Compute harmonic mean Ne across generations.

    The harmonic mean reflects that bottleneck generations dominate
    long-term effective size.  Zeros are excluded from the calculation.

    Parameters
    ----------
    ne_series : list or array of float
        Ne values for each generation.

    Returns
    -------
    float
        Harmonic mean of non-zero Ne values.
    """
    vals = [v for v in ne_series if v > 0]
    if not vals:
        return 0.0
    return len(vals) / sum(1.0 / v for v in vals)


# ---------------------------------------------------------------------------
# Optimization
# ---------------------------------------------------------------------------

def enumerate_compatible_crosses(pop, allele_pool, inheritance_mode="tetrasomic", leakage=0.0):
    """Enumerate all SI-compatible directed crosses and build effect arrays.

    Parameters
    ----------
    pop : list[tuple]
        Population genotypes.
    allele_pool : iterable
        Full allele pool (used to size the effect matrix consistently).
    inheritance_mode : {'tetrasomic', 'disomic-averaged'}
        See :func:`form_gametes`.
    leakage : float in [0, 1]
        Per-gamete SI escape probability. Under ``leakage > 0`` every
        ordered pair has nonzero compatibility, so all pairs are included;
        the ``compatibility`` value reflects the expected gamete-acceptance
        rate.

    Returns
    -------
    (list[tuple], np.ndarray, np.ndarray)
        ``compatible_crosses`` is a list of ``(maternal_idx, paternal_idx,
        compatibility)`` tuples. ``allele_effect_matrix`` has shape
        ``(n_crosses, n_alleles)`` with rows giving the expected per-allele
        frequency contribution of each cross. ``expected_gfs_per_cross``
        has shape ``(n_crosses,)`` and gives the expected GFS of offspring
        from each cross (used as the genotype-quality term in
        :func:`compute_optimal_weights`).

    Notes
    -----
    Return signature changed from 2-tuple to 3-tuple as of Task #7
    (GFS-aware optimizer). Callers that ignored GFS previously can simply
    drop the third element with a ``*_`` discard.
    """
    compatible_crosses = []
    cross_allele_effects = []
    cross_expected_gfs = []
    n_pop = len(pop)
    for i in range(n_pop):
        for j in range(n_pop):
            if i == j:
                continue
            compat = crossing_compatibility(pop[i], pop[j], inheritance_mode, leakage)
            if compat <= 0:
                continue
            offspring_dist = cross(pop[i], pop[j], inheritance_mode, leakage)
            expected_freqs = {a: 0.0 for a in allele_pool}
            expected_gfs = 0.0
            for genotype, prob in offspring_dist.items():
                expected_gfs += prob * gfs(genotype)
                for allele in genotype:
                    if allele in expected_freqs:
                        expected_freqs[allele] += prob / 4.0
            compatible_crosses.append((i, j, compat))
            cross_allele_effects.append(expected_freqs)
            cross_expected_gfs.append(expected_gfs)
    allele_effect_matrix = np.array([
        [effects[a] for a in sorted(allele_pool)]
        for effects in cross_allele_effects
    ])
    expected_gfs_array = np.array(cross_expected_gfs)
    return compatible_crosses, allele_effect_matrix, expected_gfs_array


def compute_greedy_weights(pop, allele_pool, compatible_crosses, allele_effect_matrix):
    """Compute greedy crossing weights that boost underrepresented alleles."""
    target_freq = 1.0 / len(allele_pool)
    current_freqs = allele_frequencies(pop, allele_pool)
    freq_array = np.array([current_freqs[a] for a in sorted(allele_pool)])
    n_crosses = len(compatible_crosses)
    greedy_scores = np.zeros(n_crosses)
    for k in range(n_crosses):
        deficit = target_freq - freq_array
        boost = allele_effect_matrix[k] - target_freq
        greedy_scores[k] = np.sum(np.maximum(deficit, 0) * np.maximum(boost, 0))
    if greedy_scores.sum() > 0:
        return greedy_scores / greedy_scores.sum()
    return np.ones(n_crosses) / n_crosses


def compute_optimal_weights(compatible_crosses, allele_effect_matrix, allele_pool,
                            maxiter=1000, rare_allele_indices=None,
                            preservation_weight=10.0,
                            expected_gfs_per_cross=None, gfs_weight=0.0):
    """Optimize crossing weights via L-BFGS-B with optional preservation and GFS penalties.

    Three objective terms (the latter two are optional):

    1. **Allele-frequency term** (always on): Pearson chi-squared distance
       from uniform target frequencies, ``Σ_a (expected_a - 1/n)^2 / (1/n)``.

    2. **Rare-allele preservation** (when ``rare_allele_indices`` provided):
       quadratic floor penalty per rare allele,
       ``preservation_weight * max(0, 1/n - expected_a)^2``.
       Pushes rare alleles up toward uniform (active only when below target).

    3. **Genotype-quality term** (when ``expected_gfs_per_cross`` and
       ``gfs_weight > 0`` are provided): penalty proportional to
       ``gfs_weight * mean(1 - expected_gfs)`` across the chosen crosses,
       weighted by ``w_norm``. Per Crossing Plan section 3 Phase 2 step 5.

    Parameters
    ----------
    compatible_crosses : list
        From :func:`enumerate_compatible_crosses`.
    allele_effect_matrix : np.ndarray, shape (n_crosses, n_alleles)
        From :func:`enumerate_compatible_crosses`.
    allele_pool : iterable
        Allele pool (used to compute target = 1/n).
    maxiter : int
        L-BFGS-B max iterations.
    rare_allele_indices : list[int] or None
        Column indices in ``allele_effect_matrix`` of alleles to preserve.
    preservation_weight : float
        Weight on the rare-allele preservation term. Default 10.0.
    expected_gfs_per_cross : np.ndarray or None
        From :func:`enumerate_compatible_crosses`. When provided alongside
        ``gfs_weight > 0``, activates the genotype-quality term.
    gfs_weight : float
        Weight on the GFS term. Crossing Plan suggests starting at 1.0
        (same order as ``preservation_weight``). Default 0.0 (disabled).

    Returns
    -------
    (np.ndarray, OptimizeResult)
        Normalised optimal weights (sum 1) and the raw scipy result.
    """
    target_freq = 1.0 / len(allele_pool)
    n_crosses = len(compatible_crosses)
    use_gfs = (expected_gfs_per_cross is not None and gfs_weight > 0)
    if use_gfs:
        expected_gfs_per_cross = np.asarray(expected_gfs_per_cross, dtype=float)

    def _fitness_and_grad(weights):
        w = np.abs(weights)
        w_sum = w.sum()
        if w_sum < 1e-12:
            return 1e6, np.zeros(n_crosses)
        w_norm = w / w_sum
        expected = w_norm @ allele_effect_matrix
        diff = expected - target_freq
        chi_sq = float(np.sum(diff ** 2 / target_freq))
        residuals = 2 * diff / target_freq
        g = allele_effect_matrix @ residuals
        grad = (g - np.dot(w_norm, g)) / w_sum

        if rare_allele_indices is not None and len(rare_allele_indices) > 0:
            for ai in rare_allele_indices:
                freq_i = expected[ai]
                penalty = preservation_weight * max(0, target_freq - freq_i) ** 2
                chi_sq += penalty
                if freq_i < target_freq:
                    pen_grad_i = -2 * preservation_weight * (target_freq - freq_i)
                    pen_g = allele_effect_matrix[:, ai] * pen_grad_i
                    grad += (pen_g - np.dot(w_norm, pen_g)) / w_sum

        if use_gfs:
            # mean(1 - expected_offspring_gfs) weighted by w_norm
            # = 1 - Σ_k w_norm[k] * expected_gfs[k]
            mean_gfs_loss = 1.0 - float(np.dot(w_norm, expected_gfs_per_cross))
            chi_sq += gfs_weight * mean_gfs_loss
            # ∂(mean_gfs_loss)/∂w_norm[k] = -expected_gfs[k]
            # Chain through normalisation: ∂/∂w_j = (-expected_gfs[j] - dot(w_norm, -expected_gfs)) / w_sum
            pen_g = -expected_gfs_per_cross
            grad += gfs_weight * (pen_g - np.dot(w_norm, pen_g)) / w_sum

        return chi_sq, grad

    w0 = np.ones(n_crosses) / n_crosses
    result = minimize(
        _fitness_and_grad, w0, method="L-BFGS-B",
        jac=True,
        bounds=[(0, 1)] * n_crosses,
        options={"maxiter": maxiter, "ftol": 1e-12},
    )
    optimal_weights = np.abs(result.x)
    optimal_weights = optimal_weights / optimal_weights.sum()
    return optimal_weights, result
