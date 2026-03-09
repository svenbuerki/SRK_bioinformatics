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
    "cross",
    "crossing_compatibility",
    "sample_offspring",
    "build_crossing_matrix",
    # Equilibrium analysis
    "target_frequencies",
    "distance_from_equilibrium",
    # Allele preservation
    "identify_rare_alleles",
    "get_mandatory_rare_crosses",
    "select_elites",
    # Simulation
    "simulate_generation",
    # Optimization
    "enumerate_compatible_crosses",
    "compute_greedy_weights",
    "compute_optimal_weights",
]

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

def form_gametes(genotype):
    """Return all C(4,2)=6 diploid gametes from a tetraploid genotype."""
    return list(itertools.combinations(genotype, 2))


def is_compatible(maternal_genotype, pollen_gamete):
    """Check SI compatibility: True if pollen shares no alleles with maternal plant."""
    maternal_alleles = set(maternal_genotype)
    return not any(a in maternal_alleles for a in pollen_gamete)


def cross(parent_a, parent_b):
    """Compute offspring genotype distribution for a directed cross (maternal x pollen).

    Returns dict mapping offspring genotype -> probability. Empty dict if incompatible.
    """
    maternal_gametes = form_gametes(parent_a)
    paternal_gametes = form_gametes(parent_b)
    compatible_paternal = [g for g in paternal_gametes if is_compatible(parent_a, g)]
    if not compatible_paternal:
        return {}
    offspring_counts = Counter()
    total_combinations = len(maternal_gametes) * len(compatible_paternal)
    for mg in maternal_gametes:
        for pg in compatible_paternal:
            offspring_genotype = canonical(mg + pg)
            offspring_counts[offspring_genotype] += 1
    return {g: count / total_combinations for g, count in sorted(offspring_counts.items())}


def crossing_compatibility(parent_a, parent_b):
    """Fraction of pollen gametes that pass SI (0.0 = incompatible, 1.0 = fully compatible)."""
    paternal_gametes = form_gametes(parent_b)
    compatible = sum(1 for g in paternal_gametes if is_compatible(parent_a, g))
    return compatible / len(paternal_gametes)


def sample_offspring(parent_a, parent_b):
    """Sample a single offspring from a cross (None if incompatible)."""
    offspring_dist = cross(parent_a, parent_b)
    if not offspring_dist:
        return None
    genotypes = list(offspring_dist.keys())
    probs = list(offspring_dist.values())
    idx = np.random.choice(len(genotypes), p=probs)
    return genotypes[idx]


def build_crossing_matrix(genotypes):
    """Build compatibility and outcome matrices for all ordered genotype pairs.

    Returns (compatibility_dict, outcomes_dict) where keys are (genotype_a, genotype_b).
    """
    compatibility = {}
    outcomes = {}
    for i, ga in enumerate(genotypes):
        for j, gb in enumerate(genotypes):
            if i == j:
                compatibility[(ga, gb)] = 0.0
                outcomes[(ga, gb)] = {}
            else:
                compat = crossing_compatibility(ga, gb)
                compatibility[(ga, gb)] = compat
                outcomes[(ga, gb)] = cross(ga, gb)
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


# ---------------------------------------------------------------------------
# Allele preservation
# ---------------------------------------------------------------------------

def identify_rare_alleles(population, allele_pool, threshold=0.05):
    """Return set of alleles with frequency below threshold."""
    freqs = allele_frequencies(population, allele_pool)
    return {a for a, f in freqs.items() if 0 < f < threshold}


def get_mandatory_rare_crosses(population, allele_pool, threshold=0.05,
                               already_covered=None):
    """Find one compatible cross per endangered allele not already covered.

    Only targets alleles with 1-2 carriers AND not already present in the
    elite set. Each cross is chosen to maximize probability of producing
    offspring carrying the target allele. Returns at most max_crosses entries.
    """
    carrier_counts = {a: 0 for a in allele_pool}
    for g in population:
        seen = set()
        for a in g:
            if a in carrier_counts and a not in seen:
                carrier_counts[a] += 1
                seen.add(a)
    endangered = {a for a, c in carrier_counts.items() if 0 < c <= 2}

    if already_covered:
        endangered -= already_covered

    if not endangered:
        return []

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
                od = cross(population[ci], population[j])
                if not od:
                    continue
                score = sum(p for g, p in od.items() if allele in g)
                if score > best_score:
                    best_score, best_cross = score, (ci, j)
        if best_cross is not None:
            mandatory.append(best_cross)
            # Mark all endangered alleles this cross can produce
            od = cross(population[best_cross[0]], population[best_cross[1]])
            for g, _ in od.items():
                for a in g:
                    if a in endangered:
                        covered.add(a)
    return mandatory


def select_elites(population, allele_pool, elite_frac=0.1):
    """Select elite individuals carrying the rarest alleles.

    Scores each individual by sum of inverse-frequency. Number of elites
    is max(1, elite_frac * n) — no inflation beyond that.
    """
    freqs = allele_frequencies(population, allele_pool)
    n_elite = max(1, int(len(population) * elite_frac))
    scores = [(sum(1.0 / max(freqs.get(a, 1e-12), 1e-12) for a in g), i)
              for i, g in enumerate(population)]
    scores.sort(reverse=True)
    return [idx for _, idx in scores[:n_elite]]


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------

def simulate_generation(population, n_offspring=None, crossing_plan=None,
                        allele_pool=None, preserve_rare=False,
                        rare_threshold=0.05, elite_frac=0.1):
    """Simulate one generation with optional allele preservation.

    When preserve_rare=True and allele_pool is provided:
    1. Elite individuals are retained — guarantees every present allele has a carrier.
    2. Mandatory crosses for rare alleles produce multiple offspring per rare allele.
    3. Remaining offspring slots filled by crossing plan or random mating.
    4. Final verification: any allele still missing gets its carrier directly inserted.
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

        # Step 1: Retain elites (10% of population)
        elite_indices = select_elites(population, allele_pool, elite_frac)
        for idx in elite_indices:
            next_gen.append(population[idx])

        # Step 2: Mandatory crosses for endangered alleles not covered by elites
        elite_alleles = set()
        for idx in elite_indices:
            elite_alleles.update(population[idx])
        max_mandatory = max(1, int(n_offspring * 0.15))  # cap at 15% of population
        mandatory = get_mandatory_rare_crosses(
            population, allele_pool, rare_threshold, already_covered=elite_alleles)
        for mi, pi in mandatory[:max_mandatory]:
            if len(next_gen) >= n_offspring:
                break
            child = sample_offspring(population[mi], population[pi])
            if child is not None:
                next_gen.append(child)

    # Step 3: Fill remaining slots
    if crossing_plan is not None:
        indices = list(range(len(crossing_plan)))
        weights = np.array([w for _, _, w in crossing_plan])
        weights = weights / weights.sum()
        attempts = 0
        while len(next_gen) < n_offspring and attempts < n_offspring * 20:
            idx = np.random.choice(indices, p=weights)
            mi, pi, _ = crossing_plan[idx]
            child = sample_offspring(population[mi], population[pi])
            if child is not None:
                next_gen.append(child)
            attempts += 1
    else:
        attempts = 0
        n = len(population)
        while len(next_gen) < n_offspring and attempts < n_offspring * 20:
            i, j = _random.sample(range(n), 2)
            child = sample_offspring(population[i], population[j])
            if child is not None:
                next_gen.append(child)
            attempts += 1

    # Step 4: Final safety net — verify no allele was lost
    if preserve_rare and present_alleles is not None:
        next_gen_alleles = set()
        for g in next_gen:
            next_gen_alleles.update(g)
        missing = present_alleles - next_gen_alleles
        if missing:
            # Insert carriers from the original population for any missing alleles
            for allele in missing:
                carriers = [i for i, g in enumerate(population) if allele in g]
                if carriers:
                    # Pick a random carrier and add it
                    chosen = carriers[_random.randint(0, len(carriers) - 1)]
                    next_gen.append(population[chosen])

    return next_gen


# ---------------------------------------------------------------------------
# Optimization
# ---------------------------------------------------------------------------

def enumerate_compatible_crosses(pop, allele_pool):
    """Enumerate all SI-compatible directed crosses and build allele effect matrix.

    Returns (compatible_crosses, allele_effect_matrix) where:
    - compatible_crosses: list of (maternal_idx, paternal_idx, compatibility) tuples
    - allele_effect_matrix: (n_crosses, n_alleles) array of expected allele frequencies
    """
    compatible_crosses = []
    cross_allele_effects = []
    n_pop = len(pop)
    for i in range(n_pop):
        for j in range(n_pop):
            if i == j:
                continue
            compat = crossing_compatibility(pop[i], pop[j])
            if compat > 0:
                offspring_dist = cross(pop[i], pop[j])
                expected_freqs = {a: 0.0 for a in allele_pool}
                for genotype, prob in offspring_dist.items():
                    for allele in genotype:
                        if allele in expected_freqs:
                            expected_freqs[allele] += prob / 4.0
                compatible_crosses.append((i, j, compat))
                cross_allele_effects.append(expected_freqs)
    allele_effect_matrix = np.array([
        [effects[a] for a in sorted(allele_pool)]
        for effects in cross_allele_effects
    ])
    return compatible_crosses, allele_effect_matrix


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
                            preservation_weight=10.0):
    """Optimize crossing weights via L-BFGS-B with optional rare-allele preservation.

    When rare_allele_indices is provided, adds a quadratic penalty for low expected
    frequency of those alleles, preventing the optimizer from sacrificing them.
    """
    target_freq = 1.0 / len(allele_pool)
    n_crosses = len(compatible_crosses)

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
