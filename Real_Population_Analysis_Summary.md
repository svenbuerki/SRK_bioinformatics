# Real Population Analysis: *Lepidium papilliferum* Crossing Strategy Optimization

## Overview

This analysis applies a full crossing, equilibrium, and optimization pipeline to empirical SRK genotype data from *Lepidium papilliferum* (slickspot peppergrass), a federally threatened tetraploid plant. The dataset includes **124 ingroup individuals** spanning **24 populations**, carrying **92 unique S-alleles** across the species. The central question is: *Can managed crossing strategies accelerate allele frequency convergence toward negative frequency-dependent selection (NFDS) equilibrium while preventing irreversible allele loss?*

---

## Biological Background

### Self-Incompatibility and NFDS

*L. papilliferum* is governed by a sporophytic self-incompatibility (SI) system encoded at the S-locus. Each tetraploid individual carries exactly **4 SRK allele copies** (diploid gametes carry 2). Pollen is rejected when it shares **any allele** with the maternal plant — preventing both selfing and near-kin matings.

This creates a powerful selection dynamic: rare alleles are compatible with nearly every individual in the population (their bearers have higher fitness), driving the system toward **NFDS equilibrium** — equal frequency of all alleles. The goal of managed crossing is to accelerate this natural process before rare alleles are lost to genetic drift.

### Genotype Structure

- **Ploidy**: Tetraploid (allopolyploid; 4 copies at S-locus)
- **Gamete formation**: Each parent produces C(4,2) = 6 equally probable diploid gamete types
- **Canonical representation**: Sorted tuple, e.g. `(1, 5, 23, 47)`
- **SI filtering**: Compatible cross requires paternal gamete shares **no alleles** with maternal plant

---

## Data

| Attribute | Value |
|---|---|
| Total ingroup individuals | 124 |
| Total S-alleles (ingroup) | 92 |
| Singletons (1 individual only) | 39 |
| Shared alleles (≥2 individuals) | 53 |
| Unique genotypes | 102 |
| Populations | 24 |

**Major populations** (≥5 individuals) selected for within-population analysis:

| Population | N | Alleles | Unique Genos | Baseline Chi² | Baseline Variance |
|---|---|---|---|---|---|
| Pop 25 | 31 | 26 | 21 | 1.181 | 0.001748 |
| Pop 27 | 15 | 17 | 14 | 0.133 | 0.000461 |
| Pop 67 | 22 | 20 | 21 | 0.271 | 0.000677 |
| Pop 76 | 27 | 23 | 26 | 0.412 | 0.000779 |

The remaining 20 minor populations (29 individuals total) are included in the cross-population analysis.

---

## Analysis Structure

The notebook runs two parallel pipelines:

1. **Within-population analysis** — Each major population treated independently. Optimization uses only that population's alleles and crossing pairs.
2. **Cross-population analysis** — All 124 individuals treated as one metapopulation. Enables inter-population crosses that can rescue singleton alleles.

---

## Algorithm

### Step 1: Enumerate Compatible Crosses

For each population (or the full metapopulation), all ordered maternal–paternal pairs are tested for SI compatibility. A cross `(i → j)` is compatible if every gamete that *j* can produce shares no allele with *i*.

```
Pop 25: 871 / 930 possible crosses compatible (93.7%)
Pop 27: 206 / 210 compatible (98.1%)
Pop 67: 448 / 462 compatible (97.0%)
Pop 76: 669 / 702 compatible (95.3%)
```

High compatibility rates reflect the diversity of the real allele pool — rare alleles open more mating doors.

### Step 2: Build Allele Effect Matrix

For each compatible cross, the expected allele frequency contribution to the next generation is computed analytically by enumerating all 36 gamete combinations (6 maternal × 6 paternal), filtering for SI compatibility, and summing the allele content of surviving offspring. This yields a matrix of shape `(n_crosses × n_alleles)`.

### Step 3: Optimize Crossing Weights — L-BFGS-B

Crossing weights `w` (a probability distribution over compatible crosses) are optimized using L-BFGS-B gradient descent to minimize chi-squared distance from uniform allele frequencies:

$$\text{minimize} \quad \chi^2 = \sum_k \frac{(f_k - f^*)^2}{f^*}$$

where $f_k$ is the expected frequency of allele $k$ under weights $w$, and $f^* = 1/n_{\text{alleles}}$ is the NFDS target.

**Convergence results:**

| Population | Baseline Chi² | Optimized Chi² | Converged |
|---|---|---|---|
| Pop 25 | 0.983 | ~0.000001 | Yes |
| Pop 27 | 0.123 | ~0.000000 | Yes |
| Pop 67 | 0.241 | ~0.000000 | Yes |
| Pop 76 | 0.355 | 0.022 | Yes |

Pop 76 retains a small residual due to SI structural constraints — some allele combinations cannot be equalized within this population alone.

---

## Four Crossing Strategies

All simulations run for **5 generations** (years), **10 stochastic trials** each, with identical random seeds per trial number across strategies (ensuring fair comparison). Crossing weights are **recomputed adaptively each generation**.

### 1. Random Mating (gray)
Baseline: uniform random pairing among SI-compatible pairs. No selection for allele frequency balance. Common alleles accumulate through sampling bias — SI filtering alone is insufficient to drive equilibrium.

### 2. Optimized Crossing — No Preservation (red)
L-BFGS-B weights applied each generation. Rapidly reduces variance but can drive rare alleles to zero: the optimizer may rationally sacrifice a rare allele if doing so reduces overall chi-squared. Allele extinctions accumulate over generations.

### 3. Optimized + Preservation (blue) *(recommended)*
Adds three protective mechanisms on top of the optimizer:

- **Elitism** (`elite_frac = 10%`): Individuals with the highest rarity scores (weighted inverse frequency of their alleles) are preserved unchanged into the next generation, guaranteeing rare-allele carriers survive.
- **Mandatory rare-allele crosses** (`rare_threshold = 5%`): For every allele below 5% frequency, the single highest-weight cross involving that allele is performed first, before any other crosses are selected.
- **Optimizer penalty** (`preservation_weight = 10.0`): A quadratic penalty term is added to the chi-squared objective for any rare allele expected to fall below its current frequency. The optimizer is discouraged from sacrificing rare alleles for global chi-squared gains:

$$\text{loss} = \chi^2 + \lambda \sum_{k \in \text{rare}} \max(0,\; f_k^{\text{target}} - f_k^{\text{expected}})^2$$

### 4. Optimized + Preservation + Demography (green)
Extends strategy 3 with realistic population dynamics:

- **Logistic growth**: $N_{t+1} = N_t + N_t \cdot r \cdot (1 - N_t / K)$, with Poisson noise
  - Carrying capacity $K = \max(3 N_0,\; 30)$; growth rate $r = 0.5$
- **Effective population size (Ne)** estimated each generation from the SI compatibility structure (individuals with many compatible mates contribute more to Ne)
- **Harmonic mean Ne** tracked as a long-term bottleneck indicator

Population growth provides more offspring per generation, giving the optimizer more degrees of freedom and accelerating convergence.

---

## Simulations and Results

### Within-Population: Variance Trajectory

**Year 5 results (10-trial mean):**

| Pop | Random var | Opt var | Pres var | Demo var | Pres reduction | Demo reduction |
|---|---|---|---|---|---|---|
| 25 | 0.001263 | 0.000328 | 0.000236 | 0.000100 | 81% | 92% |
| 27 | — | — | — | — | ~73% | ~84% |
| 67 | — | — | — | — | ~76% | ~87% |
| 76 | — | — | — | — | ~75% | ~85% |

**Year 1 summary for Pop 25:**
- Random: var = 0.001675 (4% improvement)
- Optimized: var = 0.000376 (78% improvement)
- Preservation: var = 0.000305 (82% improvement)
- Demography: var = 0.000270 (84% improvement); population grew 31 → 42

---

## Charts Explained

### Chart 1 — Allele Frequency Bar Charts (Set 1: Random vs. Preservation)

**What it shows**: Side-by-side allele frequencies at Years 0, 1, 2, and 5 for each major population.
- **Gray bars**: random mating
- **Blue bars**: optimized + preservation
- **Dashed horizontal line**: NFDS equilibrium target ($1/n_{\text{alleles}}$)

**How to read it**: Under random mating, common alleles remain dominant (tall bars persist). Under preservation strategy, bar heights converge toward the dashed target. Notably, even rare alleles (very short bars at Year 0) remain above zero through Year 5 — no allele is lost.

### Chart 2 — Allele Frequency Bar Charts (Set 2: Original Optimization vs. Preservation)

**What it shows**: Identical format but comparing the two optimization strategies head-to-head.
- **Red bars**: optimized, no preservation
- **Blue bars**: optimized + preservation

**How to read it**: In red, some bars disappear entirely by Year 2–5 (allele extinction). Blue bars show the same convergence trajectory but with all alleles maintained above zero. The difference is starkest in Pop 25 and 76, which have the most unequal initial frequency distributions.

### Chart 3 — Convergence Comparison Grid (Within-Population, 3 × 4)

A 3-row, 4-column figure (one column per major population):

- **Row 1 — Variance over time**: All 4 strategies plotted as lines with ±1 SD shading. The steep drop in blue/green vs. the shallow decline in gray illustrates how strongly optimization improves convergence speed.
- **Row 2 — Allele extinction count**: Red accumulates extinctions over time; blue and green hold at zero. This is the key allelic diversity cost of using the unpreserved optimizer.
- **Row 3 — Population size**: Green line only (demographic strategy), showing logistic approach to carrying capacity. Larger populations in later years drive the higher convergence rates seen in Row 1.

### Chart 4 — Cross-Population Allele Frequency Bar Charts

**What it shows**: All 92 alleles for the combined 124-individual metapopulation, formatted identically to Charts 1 and 2.

**Why it matters**: The metapopulation has 92 alleles vs. 17–26 per individual population. Many singletons (alleles present in only one individual species-wide) are visible at Year 0 as isolated tiny bars. Under preservation strategy, these singletons persist through Year 5; under random mating they are at high risk of loss.

Cross-population crosses can introduce rare alleles from one population into another — the optimizer identifies these inter-population crosses as high-value when a singleton allele carrier is paired with a compatible partner from another population.

### Chart 5 — Convergence Comparison Grid (Cross-Population, 3 × 1)

Same three-metric format as Chart 3, but for the full 124-individual metapopulation. The cross-population context achieves the highest absolute variance reduction because the optimizer has access to all 92 alleles simultaneously.

---

## Top Recommended Crosses

The notebook outputs the highest-weighted crosses from the optimizer for each population, ranked by crossing weight (probability of selection). These represent the crosses most likely to move allele frequencies toward equilibrium.

**Format**: `Individual A genotype → Individual B genotype | Weight (%)`

Inter-population crosses in the metapopulation analysis may link individuals from geographically separate populations — these are flagged as priorities if they involve carriers of singletons or near-extinct alleles.

---

## Summary of Key Findings

| Finding | Detail |
|---|---|
| Optimization dramatically outperforms random mating | 73–92% variance reduction vs. <10% for random by Year 5 |
| Unpreserved optimization causes allele loss | 1–2 allele extinctions per population over 5 years |
| Preservation strategy prevents all allele loss | 0 extinctions across all populations and trials |
| Cost of preservation is modest | <5–10% lower variance reduction vs. unpreserved optimizer |
| Demographic strategy achieves highest convergence | Population growth provides more offspring, accelerating equilibrium |
| Cross-population analysis rescues singleton alleles | 39 singletons species-wide; inter-population crosses are their only route to persistence |
| Pop 76 hardest to equilibrate within-population | Residual chi² = 0.022; some alleles require cross-population introductions |

---

## Validation

The notebook concludes with integrity checks:
- All genotypes verified as sorted length-4 tuples
- SI compatibility logic confirmed (no self-compatible crosses present)
- Allele frequency vectors sum to 1.0 per generation
- Population size consistency across trials

---

## Conclusions

Managed crossing under the **Optimized + Preservation** strategy achieves:
1. Rapid convergence toward NFDS equilibrium (allele frequency balance)
2. Complete prevention of allele extinction
3. Compatibility with realistic demographic constraints (logistic growth)

For practitioners: the top-ranked cross list from the optimizer provides a directly actionable management recommendation. Even partial implementation — prioritizing the top 10–20 crosses per season — provides substantial improvement over unmanaged random mating.

The **cross-population analysis** is especially critical for singleton alleles: without deliberate inter-population crossing, 39 unique alleles species-wide are at immediate extinction risk under any realistic demographic scenario.
