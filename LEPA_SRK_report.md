# SRK-Based Assessment of Self-Incompatibility in *LEPA*

### Background

Self-incompatibility (SI) in *LEPA* is controlled by the S-locus, where the extracellular S-domain of the S-receptor kinase (SRK) protein acts as the female determinant of pollen rejection. Each individual carries up to four allele copies (tetraploid), and compatible mating requires that pollen and pistil carry different SRK alleles. Loss of allele diversity through genetic drift therefore directly constrains reproductive success: individuals sharing the same allele are SI-incompatible and cannot produce seed together.

---

### 1. Allele Definition

We targeted the S-domain of the SRK protein — the functional "lock" in the lock-and-key recognition mechanism — to define alleles. We first identified all unique functional protein sequences within this domain across the dataset and visualised amino acid variation across positions (Figure 1: [SRK_AA_frequency_heatmap.png](figures/SRK_AA_frequency_heatmap.png)). We then applied a distance-based sensitivity analysis to cluster protein sequences into allele bins, selecting the clustering threshold that maximised biological resolution while minimising artefactual splitting (Figure 2: [SRK_protein_distance_analysis.png](figures/SRK_protein_distance_analysis.png)). Each resulting cluster represents an S-allele bin — an allele hypothesis that groups functionally equivalent proteins under a single identity.

---

### 2. Allelic Richness in *LEPA*

**Observed richness:** Across **189 individuals** sampled from Element Occurrences (EOs) spanning the species range, we identified **47 distinct S-allele bins** (Figure 3: [SRK_allele_accumulation_species.png](figures/SRK_allele_accumulation_species.png)).

**Predicted species richness:** The allele accumulation curve has not yet reached an asymptote, indicating that further sampling would discover additional alleles. Three complementary estimators were applied to the accumulation curve:

| Estimator | Predicted species richness |
|-----------|--------------------------|
| Michaelis-Menten (MM) | 65 alleles |
| Chao1 | 84 alleles |
| Consensus (MM + Chao1 mean) | 75 alleles |

The **MM estimate of 65 alleles** is adopted as the species optimum — the allele richness expected under NFDS at evolutionary equilibrium, used here as the reference baseline against which population deficits are measured. An additional **31 individuals** would need to be sampled to discover the next new allele, and reaching 80% of the MM estimate would require sampling approximately 123 more individuals, underscoring that the species SI repertoire remains substantially under-characterised.

---

### 3. Implications for Seed Production

The species optimum of 65 alleles is the critical baseline for seed production planning. Under negative frequency-dependent selection (NFDS), the SI system is self-stabilising when all alleles are present and approximately equally frequent: rare alleles are automatically favoured because they are compatible with more partners. This means:

- A population holding all 65 alleles at equal frequency maximises the proportion of compatible mating pairs and reproductive success.
- Any reduction in allele number directly reduces the fraction of compatible pairings and seed set.
- Cross-design for seed production must account for both the number of alleles present and their copy-count dosage (tetraploid genotype classes: AAAA, AAAB, AABB, AABC, ABCD).

---

### 4. Reproductive Status of Element Occurrences

**Allele richness deficit (Figure 4: [SRK_allele_accumulation_combined.png](figures/SRK_allele_accumulation_combined.png)):**

All five Element Occurrences (EOs) with sufficient sample sizes (≥25 individuals) fall far short of the species optimum of 65 alleles:

| Element Occurrence | N individuals | Observed alleles | % of species optimum (MM) | Predicted (MM) | Predicted (Chao1) |
|--------------------|:---:|:---:|:---:|:---:|:---:|
| EO25 | 32 | 11 | 17% | 15 | 12 |
| EO27 | 29 | 15 | 23% | 26 | 28 |
| EO67 | 32 | 12 | 18% | 16 | 20 |
| EO70 | 40 |  6 |  9% |  7 |  8 |
| EO76 | 25 |  9 | 14% | 12 | 12 |

EO27 retains the most allele diversity, yet even its observed count of 15 alleles represents only 23% of the species optimum. EO70, despite being the most heavily sampled EO (40 individuals), harbours only 6 allele bins — the lowest richness of all five EOs and just 9% of the species optimum, suggesting severe historical bottlenecking or founder effects at this occurrence.

**Allele set composition and sharing (Figure 5: [SRK_allele_upset_EOs.png](figures/SRK_allele_upset_EOs.png), Figure 6: [SRK_allele_sharing_heatmap_EOs.png](figures/SRK_allele_sharing_heatmap_EOs.png)):**

S-allele sets are largely private to each Element Occurrence. Only **2 alleles — Allele_044 and Allele_048 — are shared across all five EOs**, and these are precisely the alleles with the highest copy counts species-wide, consistent with the severe frequency skew documented below. The remaining alleles are partitioned among EO subsets or are exclusive to a single EO:

| Element Occurrence | Private alleles (exclusive to this EO) |
|--------------------|:---:|
| EO27 | 10 |
| EO67 |  7 |
| EO76 |  6 |
| EO25 |  5 |
| EO70 |  3 |

EO27 holds the largest private allele set (10 alleles), reinforcing its status as the most allele-rich and irreplaceable contributor to the species SI repertoire. EO70, despite being the most depauperate EO, still retains 3 alleles found nowhere else. In pairwise comparisons, EO25 and EO27 share the most alleles (5), while EO70 shares only 2–3 alleles with any other EO — underscoring its compositional isolation and the importance of inter-EO crosses for redistributing allele diversity to this occurrence.

**Allele frequency imbalance (Figure 7: [SRK_chisq_species_population_frequency_plots.png](figures/SRK_chisq_species_population_frequency_plots.png)):**

Under NFDS, alleles are expected to be maintained at approximately equal frequencies. Chi-square goodness-of-fit tests against the uniform expectation reject this null at every level:

| Level | N individuals | N alleles | χ² | *p*-value |
|-------|:---:|:---:|:---:|:---:|
| Species | 189 | 47 | 2029.2 | < 10⁻⁴⁰⁰ |
| EO25 | 32 | 11 | 76.0 | 3.0 × 10⁻¹² |
| EO27 | 29 | 15 | 60.5 | 9.6 × 10⁻⁸ |
| EO67 | 32 | 12 | 97.8 | 4.9 × 10⁻¹⁶ |
| EO70 | 40 |  6 | 82.1 | 3.0 × 10⁻¹⁶ |
| EO76 | 25 |  9 | 47.1 | 1.5 × 10⁻⁷ |

A small number of alleles dominate in each Element Occurrence. Allele_044 and Allele_048 are the most prevalent across the species: Allele_044 contributes 25 copies in EO25 and 14 copies in EO27; Allele_048 dominates EO67 (26 copies), EO70 (36 copies), and EO76 (15 copies). EO70 is particularly extreme — Allele_048 alone accounts for a disproportionate share of all allele copies in that occurrence, while many alleles are absent entirely.

**Zygosity (Figure 8: [SRK_zygosity_distribution.png](figures/SRK_zygosity_distribution.png)):**

Genotype reconstruction from the tetraploid allele copy-count matrix reveals the following dosage classes across all 189 individuals:

| Genotype class | Description | N individuals | % |
|----------------|-------------|:---:|:---:|
| AAAA | Homozygous — one allele, four copies | 105 | 56% |
| AABB | Two alleles, two copies each | 40 | 21% |
| AAAB | Two alleles, 3:1 dosage | 39 | 21% |
| AABC | Three alleles, one doubled |  5 |  3% |

A majority of individuals (56%, 105/189) carry only a single allele bin (AAAA genotype), making them functionally equivalent to homozygotes with respect to SI. These individuals are candidates for self-compatibility and cannot contribute compatible pollen to any partner carrying the same allele. The proportion is consistent across EOs: EO25 (56%), EO27 (62%), EO67 (53%), EO70 (52%), EO76 (52%).

---

### 4a. Tipping Point 1 (TP1) — Allele Richness and Frequency Evenness

TP1 is breached when a population has lost so many S-alleles relative to the species optimum that inter-population allele transfers are required to restore richness. It is evaluated using two complementary metrics (Figure 9: [SRK_TP1_tipping_point.png](figures/SRK_TP1_tipping_point.png)):

- **Proportion of species optimum retained** (`prop_optimum` = N_alleles / 65): how much of the species-level SI repertoire each EO holds.
- **Frequency evenness** (`Ne / N_alleles`): how close allele frequencies are to the equal-frequency NFDS ideal. The effective allele number Ne = 1/Σpᵢ² answers *how many equally frequent alleles would produce the same diversity as observed*. A ratio of 1.0 means perfect evenness; drift pushes it downward as dominant alleles accumulate and rare alleles are marginalised.

An EO is flagged **CRITICAL** when both criteria are breached (< 50% of species optimum and Ne/N < 0.80), **AT RISK** when only one is breached, and **OK** when neither is.

| EO | N alleles | % of species optimum | Evenness (Ne/N) | TP1 status |
|----|:---:|:---:|:---:|:---:|
| EO70 |  6 |  9% | 0.51 | **CRITICAL** |
| EO76 |  9 | 14% | 0.55 | **CRITICAL** |
| EO25 | 11 | 17% | 0.47 | **CRITICAL** |
| EO67 | 12 | 18% | 0.41 | **CRITICAL** |
| EO27 | 15 | 23% | 0.46 | **CRITICAL** |

All five EOs are CRITICAL on both axes. No EO retains more than 23% of the species allele pool, and all have evenness values between 0.41 and 0.55 — meaning that even the alleles present are distributed so unevenly that fewer than half are effectively contributing to SI function. EO67 has the lowest evenness (0.41), indicating the most extreme dominance by a few alleles; EO76 has the highest (0.55), though this is still far from the NFDS ideal of 1.0. EO27 is least depauperate on richness (23%), while EO70 is the most allele-impoverished (9%).

---

### 4b. Individual Genotypic Fitness Score (GFS) and Tipping Point 2 (TP2)

The zygosity classification above groups individuals by the *number* of unique alleles, but for seed production the *dosage balance* of those alleles also matters. A tetraploid produces diploid gametes by randomly sampling 2 of its 4 allele copies — yielding C(4,2) = 6 equally probable gamete combinations. An individual with genotype AABB (two alleles, balanced 2+2) produces heterozygous gametes in 4 of 6 combinations (GFS = 0.667), whereas an AAAB individual (two alleles, unbalanced 3+1) produces heterozygous gametes in only 3 of 6 (GFS = 0.500). This distinction — invisible to zygosity analysis alone — is captured by the **Genotypic Fitness Score (GFS)**:

$$\text{GFS}_i = 1 - \frac{\sum_k n_k\,(n_k - 1)}{12}$$

where $n_k$ is the copy number of allele $k$ and the denominator normalises to the tetraploid gamete space.

| Genotype | GFS | Heterozygous gametes |
|----------|-----|----------------------|
| ABCD | 1.000 | 6 / 6 |
| AABC | 0.833 | 5 / 6 |
| AABB | 0.667 | 4 / 6 |
| AAAB | 0.500 | 3 / 6 |
| AAAA | 0.000 | 0 / 6 |

**Tipping Point 2 (TP2)** is breached when the distribution of alleles across individuals has degraded to the point that managed crossing within an EO cannot restore fitness without external allele introductions. Two criteria are evaluated per EO:

- **TP2-mean**: mean GFS < 0.667 (EO average below AABB level)
- **TP2-AAAA**: proportion AAAA > 30%

An EO breaching both simultaneously is flagged **CRITICAL**; one criterion is **AT RISK**; neither is **OK**.

**EO-level GFS results (Figure 10: [SRK_GFS_plots_p1_composition_proportional.png](figures/SRK_GFS_plots_p1_composition_proportional.png), Figure 11: [SRK_GFS_plots_p2_individual_jitter.png](figures/SRK_GFS_plots_p2_individual_jitter.png), Figure 12: [SRK_GFS_plots_p3_TP2_tipping_point.png](figures/SRK_GFS_plots_p3_TP2_tipping_point.png)):**

| EO | N | mean GFS | % AAAA | % AAAB | % AABB | % AABC | TP2 status |
|----|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| EO27 | 29 | 0.224 | 62% | 17% | 21% | 0% | **CRITICAL** |
| EO25 | 32 | 0.255 | 56% | 22% | 22% | 0% | **CRITICAL** |
| EO76 | 25 | 0.273 | 52% | 28% | 20% | 0% | **CRITICAL** |
| EO70 | 40 | 0.283 | 53% | 23% | 23% | 3% | **CRITICAL** |
| EO67 | 32 | 0.297 | 53% | 16% | 25% | 6% | **CRITICAL** |

All five EOs are CRITICAL. No EO approaches a mean GFS consistent with a functioning SI system. **EO67 is the least degraded**, with the highest mean GFS and the only AABC individuals (n = 2), alongside EO70 (n = 1). These five AABC individuals — producing 5/6 heterozygous gametes — are the highest-priority seed parents in the entire species.

**Seed production priority within each EO** (top-ranked individuals by GFS):

| EO | Individual | Genotype | GFS |
|----|-----------|----------|-----|
| EO67 | Library007_barcode82 | AABC | 0.833 |
| EO67 | Library007_barcode86 | AABC | 0.833 |
| EO70 | Library008_barcode22 | AABC | 0.833 |
| EO25 | Library002_barcode31 | AABB | 0.667 |
| EO27 | Library006_barcode38 | AABB | 0.667 |
| EO76 | Library001_barcode09 | AABB | 0.667 |

Full ranked lists per EO are available in `SRK_individual_GFS.tsv`; EO-level summaries and TP2 flags in `SRK_EO_GFS_summary.tsv`.

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

Random mating causes a mean loss of 7.4 allele bins by year 3 at the species level, whereas managed crosses with preservation eliminate allele loss entirely while driving allele frequency variance down by 95–97%. Cross-population crosses (inter-EO transfers) are critical: they are the primary mechanism for redistributing alleles that are private to single EOs and absent from the broader gene pool.

---

### 6. Conclusion

Genetic drift, promoted by habitat fragmentation, is severely eroding allele diversity and frequency balance in *LEPA* Element Occurrences at two distinct levels.

**Tipping Point 1 — Allele richness deficit (Figure 9: [SRK_TP1_tipping_point.png](figures/SRK_TP1_tipping_point.png)).** With 189 individuals sampled across five major EOs, only 47 of an estimated 65 species-level S-allele bins (MM estimate) have been observed. Each EO retains a small, largely private subset of the species SI repertoire (9–23% of the species optimum), allele frequencies are highly skewed from the equal distribution expected under NFDS (χ² *p* < 10⁻⁷ at every level), and frequency evenness (Ne/N) ranges from only 0.41 to 0.55 — meaning fewer than half of even the observed alleles are effectively contributing to SI function. Allele sets are almost entirely non-overlapping across EOs — only 2 alleles are shared across all five. All five EOs are flagged **CRITICAL** for TP1. EO70 is the most allele-depauperate (6 bins, 9% of optimum; evenness = 0.51); EO27 retains the most (15 bins, 23% of optimum) and the largest set of private alleles (n = 10).

**Tipping Point 2 — Genotypic fitness collapse (Figure 12: [SRK_GFS_plots_p3_TP2_tipping_point.png](figures/SRK_GFS_plots_p3_TP2_tipping_point.png)).** Beyond allele counts, the *distribution of alleles across individuals* has degraded severely. The Genotypic Fitness Score (GFS) — the proportion of heterozygous gametes a tetraploid individual produces — reveals that all five EOs are CRITICAL: mean GFS ranges from 0.22 to 0.30 (well below the AABB benchmark of 0.667), and 52–62% of individuals per EO are AAAA (producing zero heterozygous gametes). No ABCD individual exists in the dataset. Only five individuals across the entire species carry an AABC genotype (GFS = 0.833) — all in EO67 and EO70 — making them the highest-priority seed parents for near-term managed crossing.

Critically, these two tipping points interact: even if allele richness were restored through inter-EO transfers (TP1 intervention), the benefit would be limited if the incoming alleles are absorbed into AAAA or AAAB individuals. Effective restoration therefore requires simultaneously targeting allele richness (inter-EO transfers of rare alleles) and genotype quality (crosses designed to produce AABB, AABC, and ultimately ABCD offspring).

Crossing simulations demonstrate that managed, optimisation-guided crosses can reduce allele frequency variance by 77–95% within a single generation compared to random mating, with zero allele loss under preservation strategies. Seed production efforts should prioritise: (1) AABC and AABB individuals as seed parents, ranked by GFS within each EO; (2) inter-EO crosses targeting rare-allele carriers to address the richness deficit; and (3) crossing designs that pair AAAB individuals carrying complementary allele pairs, which can directly produce higher-GFS offspring in the next generation.
