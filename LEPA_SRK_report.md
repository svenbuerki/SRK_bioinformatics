# SRK-Based Assessment of Self-Incompatibility in *LEPA*

### Background

Self-incompatibility (SI) in *LEPA* is controlled by the S-locus, where the extracellular S-domain of the S-receptor kinase (SRK) protein acts as the female determinant of pollen rejection. Each individual carries up to four allele copies (tetraploid), and compatible mating requires that pollen and pistil carry different SRK alleles. Loss of allele diversity through genetic drift therefore directly constrains reproductive success: individuals sharing the same allele are SI-incompatible and cannot produce seed together.

---

### 1. Allele Definition

We targeted the S-domain of the SRK protein — the functional "lock" in the lock-and-key recognition mechanism — to define alleles. We first identified all unique functional protein sequences within this domain across the dataset and visualised amino acid variation across positions (Figure 1: [SRK_AA_frequency_heatmap.pdf](figures/SRK_AA_frequency_heatmap.pdf)). We then applied a distance-based sensitivity analysis to cluster protein sequences into allele bins, selecting the clustering threshold that maximised biological resolution while minimising artefactual splitting (Figure 2: [SRK_protein_distance_analysis.pdf](figures/SRK_protein_distance_analysis.pdf)). Each resulting cluster represents an S-allele bin — an allele hypothesis that groups functionally equivalent proteins under a single identity.

---

### 2. Allelic Richness in *LEPA*

**Observed richness:** Across **189 individuals** sampled from Element Occurrences (EOs) spanning the species range, we identified **47 distinct S-allele bins** (Figure 3: [SRK_allele_accumulation_curves.pdf](figures/SRK_allele_accumulation_curves.pdf)).

**Predicted species richness:** The allele accumulation curve has not yet reached an asymptote, indicating that further sampling would discover additional alleles. Three complementary estimators were applied to the accumulation curve:

| Estimator | Predicted species richness |
|-----------|--------------------------|
| Michaelis-Menten (MM) | 65 alleles |
| Chao1 | 84 alleles |
| Consensus (MM + Chao1 mean) | **75 alleles** |

The consensus estimate of **75 alleles** is adopted as the species optimum — the allele richness expected in the absence of genetic drift. An additional **31 individuals** would need to be sampled to discover the next new allele, and reaching 80% of the MM estimate would require sampling approximately 123 more individuals, underscoring that the species SI repertoire remains substantially under-characterised.

---

### 3. Implications for Seed Production

The species optimum of 75 alleles is the critical baseline for seed production planning. Under negative frequency-dependent selection (NFDS), the SI system is self-stabilising when all alleles are present and approximately equally frequent: rare alleles are automatically favoured because they are compatible with more partners. This means:

- A population holding all 75 alleles at equal frequency maximises the proportion of compatible mating pairs and reproductive success.
- Any reduction in allele number directly reduces the fraction of compatible pairings and seed set.
- Cross-design for seed production must account for both the number of alleles present and their copy-count dosage (tetraploid genotype classes: AAAA, AAAB, AABB, AABC, ABCD).

---

### 4. Reproductive Status of Element Occurrences

**Allele richness deficit (Figure 3: [SRK_allele_accumulation_curves.pdf](figures/SRK_allele_accumulation_curves.pdf)):**

All five Element Occurrences (EOs) with sufficient sample sizes (≥25 individuals) fall far short of the species optimum of 75 alleles:

| Element Occurrence | N individuals | Observed alleles | Predicted (MM) | Predicted (Chao1) | % of species optimum (MM) |
|--------------------|:---:|:---:|:---:|:---:|:---:|
| EO25 | 32 | 11 | 15 | 12 | 23% |
| EO27 | 29 | 15 | 26 | 28 | 40% |
| EO67 | 32 | 12 | 16 | 20 | 25% |
| EO70 | 40 |  6 |  7 |  8 | 11% |
| EO76 | 25 |  9 | 12 | 12 | 18% |

EO27 retains the most allele diversity, yet even its asymptotic estimate of 26 alleles represents only 40% of the species optimum. EO70, despite being the most heavily sampled EO (40 individuals), harbours only 6 allele bins — the lowest richness of all five EOs and just 11% of the species optimum, suggesting severe historical bottlenecking or founder effects at this occurrence.

**Allele set composition and sharing (Figure: [SRK_allele_upset_EOs.pdf](figures/SRK_allele_upset_EOs.pdf), [SRK_allele_sharing_heatmap_EOs.pdf](figures/SRK_allele_sharing_heatmap_EOs.pdf)):**

S-allele sets are largely private to each Element Occurrence. Only **2 alleles — Allele_044 and Allele_048 — are shared across all five EOs**, and these are precisely the alleles with the highest copy counts species-wide, consistent with the severe frequency skew documented above. The remaining alleles are partitioned among EO subsets or are exclusive to a single EO:

| Element Occurrence | Private alleles (exclusive to this EO) |
|--------------------|:---:|
| EO27 | 10 |
| EO67 |  7 |
| EO76 |  6 |
| EO25 |  5 |
| EO70 |  3 |

EO27 holds the largest private allele set (10 alleles), reinforcing its status as the most allele-rich and irreplaceable contributor to the species SI repertoire. EO70, despite being the most depauperate EO, still retains 3 alleles found nowhere else. In pairwise comparisons, EO25 and EO27 share the most alleles (5), while EO70 shares only 2–3 alleles with any other EO — underscoring its compositional isolation and the importance of inter-EO crosses for redistributing allele diversity to this occurrence. No single Element Occurrence captures the full SI diversity of the species.

**Allele frequency imbalance (Figure 4: [SRK_chisq_species_population_frequency_plots.pdf](figures/SRK_chisq_species_population_frequency_plots.pdf)):**

Under NFDS, alleles are expected to be maintained at approximately equal frequencies. Chi-square goodness-of-fit tests against the uniform expectation reject this null at every level:

| Level | N individuals | N alleles | χ² | *p*-value |
|-------|:---:|:---:|:---:|:---:|
| Species | 189 | 47 | 2029.2 | < 10⁻⁴⁰⁰ |
| EO25 | 32 | 11 | 76.0 | 3.0 × 10⁻¹² |
| EO27 | 29 | 15 | 60.5 | 9.6 × 10⁻⁸ |
| EO67 | 32 | 12 | 97.8 | 4.9 × 10⁻¹⁶ |
| EO70 | 40 |  6 | 82.1 | 3.0 × 10⁻¹⁶ |
| EO76 | 25 |  9 | 47.1 | 1.5 × 10⁻⁷ |

A small number of alleles dominate in each Element Occurrence. Allele_044 and Allele_048 are the most prevalent across the species: Allele_044 contributes 54 copies in EO25 and 32 copies in EO27; Allele_048 dominates EO67 (46 copies), EO70 (74 copies), and EO76 (21 copies). EO70 is particularly extreme — Allele_048 alone accounts for 46% of all allele copies in that occurrence, while many alleles are absent entirely.

**Zygosity and self-compatibility risk:**

Genotype reconstruction from the tetraploid allele copy-count matrix reveals the following dosage classes across all 189 individuals:

| Genotype class | Description | N individuals | % |
|----------------|-------------|:---:|:---:|
| AAAA | Homozygous — one allele, four copies | 105 | 56% |
| AABB | Two alleles, two copies each | 40 | 21% |
| AAAB | Two alleles, 3:1 dosage | 39 | 21% |
| AABC | Three alleles, one doubled | 5 | 3% |

A majority of individuals (56%, 105/189) carry only a single allele bin (AAAA genotype), making them functionally equivalent to homozygotes with respect to SI. These individuals are candidates for self-compatibility and cannot contribute compatible pollen to any partner carrying the same allele. The proportion is consistent across EOs: EO25 (56%), EO27 (62%), EO67 (53%), EO70 (52%), EO76 (52%).

---

### 5. Crossing Strategy Simulations

To evaluate how different managed crossing strategies could restore allele frequency balance and prevent further allele loss, we simulated five years of crossing under four strategies across all five major EOs and for the combined species-wide metapopulation. All simulations used S-allele tetraploid genotypes with self-incompatibility constraints enforced, and were replicated across 10 stochastic trials.

**Strategies compared:**

1. **Random mating** — uniform random pairing among all SI-compatible pairs (baseline)
2. **Optimised** — crosses weighted by L-BFGS-B gradient descent to minimise χ² distance from equal allele frequencies, recomputed each generation
3. **Optimised + preservation** — adds rare-allele protection via elitism, mandatory rare-allele crosses, and an optimiser penalty
4. **Optimised + preservation + demography** — as above, with logistic population growth toward a carrying capacity *K* and Poisson demographic stochasticity

**Within-EO results (allele frequency variance reduction vs random mating):**

| EO | N | Allele bins | Initial χ² | Year 1 Pres. | Year 1 Demo. | Year 5 Pres. | Year 5 Demo. | Allele loss (rand. yr5) | Allele loss (pres. yr5) |
|----|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| EO25 | 32 | 11 | 1.53 | 88% | 91% | 77% | 91% | 0.8 | 0.0 |
| EO27 | 29 | 15 | 1.36 | 79% | 80% | 81% | 95% | 3.3 | 0.0 |
| EO67 | 32 | 12 | 1.31 | 77% | 80% | 89% | 93% | 1.9 | 0.0 |
| EO70 | 40 |  6 | 1.08 | 83% | 86% | 91% | 95% | 0.5 | 0.0 |
| EO76 | 25 |  9 | 0.58 | 68% | 74% | 84% | 90% | 1.1 | 0.0 |

Variance reduction percentages are relative to random mating at the same year. All optimised strategies achieved zero allele loss across all EOs and all years, whereas random mating caused progressive allele extinction in every EO except EO70 (which already has too few alleles to lose further). EO27 is the most at risk under random mating, losing on average 3.3 allele bins within five years.

**Species-wide (cross-population) results:**

Treating all 189 individuals as a single metapopulation, with 30,141 compatible crosses available (84.8% of all directed pairs), the optimised strategies achieved rapid convergence:

| Year | Random mating (var.) | Opt. + preserved (var.) | Variance reduction | Allele loss (random) | Allele loss (preserved) |
|:----:|:---:|:---:|:---:|:---:|:---:|
| 0 | 0.00242 | 0.00242 | — | 0.0 | 0.0 |
| 1 | 0.00182 | 0.00010 | 95% | 3.5 | 0.0 |
| 2 | 0.00130 | 0.00005 | 97% | 5.9 | 0.0 |
| 3 | 0.00088 | 0.00003 | 96% | 7.4 | 0.0 |

Random mating causes a mean loss of 7.4 allele bins by year 3 at the species level, whereas managed crosses with preservation eliminate allele loss entirely while driving allele frequency variance down by 95–97%. Cross-population crosses (inter-EO transfers) are critical: they are the primary mechanism for redistributing alleles that are private to single EOs and absent from the broader gene pool. Figure 5 ([sallele_cross_pop_allele_freq_random.png](modeling/outputs/figures/sallele_cross_pop_allele_freq_random.png) and [sallele_cross_pop_allele_freq_random_vs_preservation.png](modeling/outputs/figures/sallele_cross_pop_allele_freq_random_vs_preservation.png)) illustrates allele frequency distributions under each strategy over years 0–3.

Per-EO allele frequency trajectories and convergence plots are available in `modeling/outputs/figures/sallele_pop*`.

---

### 6. Conclusion

Genetic drift, promoted by habitat fragmentation, is severely eroding allele diversity and frequency balance in *LEPA* Element Occurrences. With 189 individuals sampled across five major EOs, only 47 of an estimated 75 species-level S-allele bins have been observed. Each Element Occurrence retains a small, largely private subset of the species SI repertoire (9–40% of the species optimum), allele frequencies are highly skewed from the equal distribution expected under NFDS (χ² *p* < 10⁻⁷ at every level), and 56% of individuals carry AAAA genotypes consistent with reduced SI function.

EO70, newly characterised in this analysis, is the most allele-depauperate of the five major EOs: with only 6 allele bins and extreme dominance by two alleles (Allele_048 and Allele_044), it faces the highest risk of SI collapse. At the same time, its 40 individuals represent the largest within-EO pool for managed crosses and it shows the fastest simulated convergence under optimisation (96% variance reduction by year 2), making it a tractable target for intervention.

Crossing simulations demonstrate that managed, optimisation-guided crosses can reduce allele frequency variance by 77–95% within a single generation compared to random mating, while completely preventing the allele extinction that random mating consistently produces. Cross-population transfers — particularly those moving rare alleles between EOs — are essential for approaching the species optimum, as no single EO holds more than 40% of estimated species allele richness. Seed production efforts should therefore prioritise inter-EO crosses targeting rare-allele carriers, guided by the optimal crossing weights computed in the simulation framework.
