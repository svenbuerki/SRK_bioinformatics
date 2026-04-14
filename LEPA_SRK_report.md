# SRK-Based Assessment of Self-Incompatibility in *LEPA*

### Background

Self-incompatibility (SI) in *LEPA* is controlled by the S-locus, where the extracellular S-domain of the S-receptor kinase (SRK) protein acts as the female determinant of pollen rejection. Each individual carries up to four allele copies (tetraploid), and compatible mating requires that pollen and pistil carry different SRK alleles. Under negative frequency-dependent selection (NFDS), all S-alleles are maintained at approximately equal frequencies, maximising the proportion of compatible mating pairs. In small or isolated populations, however, genetic drift counteracts NFDS — reducing allele richness, skewing allele frequencies, and degrading individual genotype quality — with direct consequences for reproductive success.

The analyses presented here address a single overarching conservation question: **what is the impact of genetic drift on reproductive success?** Three sequential questions structure the assessment:

1. **What is the S-allele richness of the species?** The species-level allele pool approximates the NFDS richness equilibrium and serves as the reference baseline against which population-level deficits are measured.
2. **How did genetic drift impact populations' S-allele diversity?** Tipping Point 1 (TP1) identifies Element Occurrences (EOs) where allele loss and frequency skew are so severe that inter-population allele transfers are required to restore the SI system.
3. **How did genetic drift impact the reproductive fitness of individuals within populations?** Tipping Point 2 (TP2) identifies EOs where allele copy-count imbalances at the individual level have degraded reproductive output to the point that within-population crosses alone cannot restore fitness.

---

### 1. Allele Definition

We targeted the S-domain of the SRK protein — the functional "lock" in the lock-and-key recognition mechanism — to define alleles. We first identified all unique functional protein sequences within this domain across the dataset and visualised amino acid variation across positions (Figure 1). We then applied a distance-based sensitivity analysis to cluster protein sequences into allele bins, selecting the clustering threshold that maximised biological resolution while minimising artefactual splitting (Figure 2). Each resulting cluster represents an S-allele bin — an allele hypothesis that groups functionally equivalent proteins under a single identity.

![Figure 1: SRK amino acid frequency heatmap](figures/SRK_AA_frequency_heatmap.png)

![Figure 2: SRK protein distance analysis](figures/SRK_protein_distance_analysis.png)

---

### 2. Allelic Richness in *LEPA*

**Observed richness:** Across **189 individuals** sampled from **26 population localities** spanning the species range, we identified **47 distinct S-allele bins** (Figure 3). The sample comprises five main EOs with sufficient individuals for population-level analysis (EO25, EO27, EO67, EO70, EO76; n = 32, 29, 32, 40, 25 respectively; total n = 158), plus 21 additional small localities (n = 1–3 each; total n = 31) contributing to species-level allele discovery.

![Figure 3: Species-level SRK allele accumulation curve](figures/SRK_allele_accumulation_species.png)

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

**Allele richness deficit (Figure 4):**

![Figure 4: SRK allele accumulation curves — EO comparison](figures/SRK_allele_accumulation_combined.png)

All five Element Occurrences (EOs) with sufficient sample sizes (≥25 individuals) fall far short of the species optimum of 65 alleles:

| Element Occurrence | N individuals | Observed alleles | % of species optimum (MM) | Predicted (MM) | Predicted (Chao1) |
|--------------------|:---:|:---:|:---:|:---:|:---:|
| EO25 | 32 | 11 | 17% | 15 | 12 |
| EO27 | 29 | 15 | 23% | 26 | 28 |
| EO67 | 32 | 12 | 18% | 16 | 20 |
| EO70 | 40 |  6 |  9% |  7 |  8 |
| EO76 | 25 |  9 | 14% | 12 | 12 |

EO27 retains the most allele diversity, yet even its observed count of 15 alleles represents only 23% of the species optimum. EO70, despite being the most heavily sampled EO (40 individuals), harbours only 6 allele bins — the lowest richness of all five EOs and just 9% of the species optimum, suggesting severe historical bottlenecking or founder effects at this occurrence.

**S-allele erosion by genetic drift (Figure 4b):**

![Figure 4b: S-allele erosion by genetic drift per EO](figures/SRK_allele_accumulation_drift_erosion.png)

The allele richness deficits documented above are not sampling artefacts — they reflect irreversible genetic erosion. Figure 4b decomposes each EO's deficit relative to the 65-allele species optimum into two components: alleles predicted to exist in the EO but not yet detected (light blue; EO MM estimate minus observed count), and alleles lost to genetic drift (red; species MM minus EO MM estimate). The red segment — alleles that are almost certainly absent from each EO — dominates catastrophically across all five EOs:

| Element Occurrence | Observed | Predicted undetected | Lost to genetic drift | % lost |
|--------------------|:---:|:---:|:---:|:---:|
| EO70 |  6 | ~1 | ~58 | **89%** |
| EO76 |  9 | ~3 | ~53 | **82%** |
| EO25 | 11 | ~4 | ~50 | **77%** |
| EO67 | 12 | ~4 | ~49 | **75%** |
| EO27 | 15 | ~11 | ~39 | **60%** |

Even in EO27 — the least affected occurrence — an estimated 39 of the 65 species-level S-allele bins (60%) have been permanently lost from the local gene pool. In EO70, the most severely eroded occurrence, approximately 58 of 65 alleles (89%) are gone. The predicted-undetected component is uniformly small (1–11 alleles), confirming that further sampling within these EOs will not close the deficit: the missing alleles are not hidden by insufficient effort — they no longer exist in these populations. Genetic drift has dismantled the SI system at its source, and inter-EO allele transfers are the only viable pathway to restoration.

**Allele set composition and sharing (Figures 5–6):**

![Figure 5: SRK allele upset plot — EO overlap](figures/SRK_allele_upset_EOs.png)

![Figure 6: SRK allele sharing heatmap — EOs](figures/SRK_allele_sharing_heatmap_EOs.png)

S-allele sets are largely private to each Element Occurrence. Only **2 alleles — Allele_044 and Allele_048 — are shared across all five EOs**, and these are precisely the alleles with the highest copy counts species-wide, consistent with the severe frequency skew documented below. The remaining alleles are partitioned among EO subsets or are exclusive to a single EO:

| Element Occurrence | Private alleles (exclusive to this EO) |
|--------------------|:---:|
| EO27 | 10 |
| EO67 |  7 |
| EO76 |  6 |
| EO25 |  5 |
| EO70 |  3 |

EO27 holds the largest private allele set (10 alleles), reinforcing its status as the most allele-rich and irreplaceable contributor to the species SI repertoire. EO70, despite being the most depauperate EO, still retains 3 alleles found nowhere else. In pairwise comparisons, EO25 and EO27 share the most alleles (5), while EO70 shares only 2–3 alleles with any other EO — underscoring its compositional isolation and the importance of inter-EO crosses for redistributing allele diversity to this occurrence.

**Allele frequency imbalance (Figure 7):**

![Figure 7: SRK allele frequency chi-square plots](figures/SRK_chisq_species_population_frequency_plots.png)

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

**Zygosity (Figure 8):**

![Figure 8: SRK zygosity distribution](figures/SRK_zygosity_distribution.png)

Genotype reconstruction from the tetraploid allele copy-count matrix reveals the following dosage classes across all 189 individuals:

| Genotype class | Description | N individuals | % |
|----------------|-------------|:---:|:---:|
| AAAA | Homozygous — one allele, four copies | 105 | 56% |
| AABB | Two alleles, two copies each | 40 | 21% |
| AAAB | Two alleles, 3:1 dosage | 39 | 21% |
| AABC | Three alleles, one doubled |  5 |  3% |

A majority of individuals (56%, 105/189) carry only a single allele bin (AAAA genotype), making them functionally equivalent to homozygotes with respect to SI. These individuals are candidates for self-compatibility and cannot contribute compatible pollen to any partner carrying the same allele. The proportion is consistent across EOs: EO25 (56%), EO27 (62%), EO67 (53%), EO70 (52%), EO76 (52%).

---

### 4a. Tipping Point 1 (TP1) — Health of the SI System

TP1 assesses the **health of the SI system** — the degree to which the SRK allele pool within a population is capable of sustaining compatible mating. It is breached when allele loss is so severe that inter-population allele transfers are required to restore SI function. Two complementary questions structure the assessment (Figure 9):

- **How many different alleles has a population retained?** (x-axis: `prop_optimum` = N_alleles / 65 — the proportion of the species-level SI repertoire still present in the EO. A population holding all alleles can offer every individual a large pool of compatible partners; as alleles are lost, compatible pairings become progressively rarer.)
- **How evenly are the remaining alleles distributed across individuals?** (y-axis: `Ne / N_alleles` — the ratio of effective to observed allele number. The effective allele number Ne = 1/Σpᵢ² answers *how many equally frequent alleles would produce the same level of diversity as observed*. A ratio of 1.0 means perfect evenness — the NFDS ideal in which every allele contributes equally to compatible crosses; drift and dominance push it downward as a few alleles monopolise copy numbers and rare alleles are marginalised.)

A population with all alleles present and evenly distributed maximises compatible mating pairs and SI system health. TP1 identifies where both dimensions have degraded beyond the point at which within-population crossing alone can restore SI function. An EO is flagged **CRITICAL** when both criteria are breached (< 50% of species optimum and Ne/N < 0.80), **AT RISK** when only one is breached, and **OK** when neither is.

![Figure 9: TP1 tipping point — health of the SI system](figures/SRK_TP1_tipping_point.png)

| EO | N alleles | % of species optimum | Evenness (Ne/N) | TP1 status |
|----|:---:|:---:|:---:|:---:|
| EO70 |  6 |  9% | 0.51 | **CRITICAL** |
| EO76 |  9 | 14% | 0.55 | **CRITICAL** |
| EO25 | 11 | 17% | 0.47 | **CRITICAL** |
| EO67 | 12 | 18% | 0.41 | **CRITICAL** |
| EO27 | 15 | 23% | 0.46 | **CRITICAL** |

All five EOs are CRITICAL on both axes. No EO retains more than 23% of the species allele pool, and all have evenness values between 0.41 and 0.55 — meaning that even the alleles present are distributed so unevenly that fewer than half are effectively contributing to SI function. EO67 has the lowest evenness (0.41), indicating the most extreme dominance by a few alleles; EO76 has the highest (0.55), though this is still far from the NFDS ideal of 1.0. EO27 is least depauperate on richness (23%), while EO70 is the most allele-impoverished (9%).

**Interpreting the two axes: bottleneck versus ongoing drift**

The richness and evenness axes of the TP1 plot capture different aspects of demographic history and together allow a partial — though not definitive — reading of the processes that generated the observed deficits. Low allele richness (prop_optimum) is the primary signature of a founder event or population bottleneck: the rapid, simultaneous loss of alleles when population size was severely reduced. Low frequency evenness (Ne/N) is more characteristic of ongoing genetic drift in small populations — the progressive accumulation of one or a few dominant alleles at the expense of rare survivors. An EO that experienced only a historical bottleneck with no subsequent drift would be expected to show low richness but moderate evenness (founders sampled stochastically, but the alleles they carried may be at similar frequencies). An EO experiencing only ongoing drift without a prior bottleneck would retain more alleles but show declining evenness as dominant frequencies grow.

LEPA EOs are severely deficient on **both** axes simultaneously, which is most consistent with a historical bottleneck that reduced richness sharply, followed by prolonged ongoing drift that further eroded evenness among the surviving alleles. This two-phase interpretation is reinforced by the allele composition analysis (Figures 5–6): if all EOs had been founded from the same bottlenecked ancestral pool, they would be missing the *same* alleles. The observation that each EO holds a largely private allele set — with only two alleles shared across all five EOs — points instead to **independent drift in isolated populations**, each losing different alleles stochastically, superimposed on a possible earlier species-level richness reduction. Formally resolving the relative contribution of each process requires genome-wide neutral marker data and demographic modelling (see Section 5).

---

### 4b. Individual Genotypic Fitness Score (GFS) and Tipping Point 2 (TP2)

The zygosity classification above groups individuals by the *number* of unique alleles, but for seed production the *dosage balance* of those alleles also matters. In self-incompatible plants, an individual's value as a breeding partner depends not only on which key/lock types it carries, but also on the number of distinct alleles present across its genome copies. Individuals with more key/lock types can participate in more compatible crosses, making them especially valuable in a managed breeding program. A tetraploid produces diploid gametes by randomly sampling 2 of its 4 allele copies — yielding C(4,2) = 6 equally probable gamete combinations. An individual with genotype AABB (two alleles, balanced 2+2) produces heterozygous gametes in 4 of 6 combinations (GFS = 0.667), whereas an AAAB individual (two alleles, unbalanced 3+1) produces heterozygous gametes in only 3 of 6 (GFS = 0.500). This distinction — invisible to zygosity analysis alone — is captured by the **Genotypic Fitness Score (GFS)**:

$$\text{GFS}_i = 1 - \frac{\sum_k n_k\,(n_k - 1)}{12}$$

where $n_k$ is the copy number of allele $k$ and the denominator normalises to the tetraploid gamete space.

| Genotype | GFS | Heterozygous gametes |
|----------|-----|----------------------|
| ABCD | 1.000 | 6 / 6 |
| AABC | 0.833 | 5 / 6 |
| AABB | 0.667 | 4 / 6 |
| AAAB | 0.500 | 3 / 6 |
| AAAA | 0.000 | 0 / 6 |

**Tipping Point 2 (TP2)** assesses **reproductive success** — the degree to which the individuals within a population can actually participate in compatible crosses and contribute allelic diversity to the next generation. It is breached when individual-level fitness has degraded to the point that managed crossing within an EO cannot restore reproductive output without external allele introductions. Two complementary questions structure the assessment:

- **What fraction of individuals are reproductive dead-ends?** (x-axis: `prop_AAAA` > 30% — the proportion of individuals producing only homotypic gametes. An AAAA individual cannot contribute to any compatible cross regardless of its mate, making it a complete loss to the breeding pool.)
- **How much reproductive capacity does the average individual still have?** (y-axis: `mean GFS` < 0.667 — the EO average falls below the AABB benchmark of producing heterozygous gametes in 4 of 6 combinations. When the mean GFS is low, even the individuals that are not dead-ends contribute relatively little diversity per cross.)

An EO breaching both simultaneously is flagged **CRITICAL**; one criterion is **AT RISK**; neither is **OK**.

**EO-level GFS results (Figures 10–12):**

![Figure 10: GFS genotype composition — proportional](figures/SRK_GFS_plots_p1_composition_proportional.png)

![Figure 11: GFS individual scores — jitter plot](figures/SRK_GFS_plots_p2_individual_jitter.png)

![Figure 12: TP2 tipping point — mean GFS vs proportion AAAA](figures/SRK_GFS_plots_p3_TP2_tipping_point.png)

| EO | N | mean GFS | % AAAA | % AAAB | % AABB | % AABC | TP2 status |
|----|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| EO27 | 29 | 0.224 | 62% | 17% | 21% | 0% | **CRITICAL** |
| EO25 | 32 | 0.255 | 56% | 22% | 22% | 0% | **CRITICAL** |
| EO76 | 25 | 0.273 | 52% | 28% | 20% | 0% | **CRITICAL** |
| EO70 | 40 | 0.283 | 53% | 23% | 23% | 3% | **CRITICAL** |
| EO67 | 32 | 0.297 | 53% | 16% | 25% | 6% | **CRITICAL** |

All five EOs are CRITICAL. No EO approaches a mean GFS consistent with a functioning SI system. **EO67 is the least degraded**, with the highest mean GFS and the only AABC individuals (n = 2), alongside EO70 (n = 1). These five AABC individuals — producing 5/6 heterozygous gametes — are the highest-priority seed parents in the entire species.

**Proportion of individuals supporting reproductive effort (Figure 13):**

![Figure 13: Proportion of individuals supporting reproductive effort per EO](figures/SRK_GFS_reproductive_effort.png)

Fewer than half of individuals in any EO carry more than one distinct SRK allele and can therefore contribute allelic diversity to compatible crosses. The proportion of "supporting" individuals (GFS > 0) ranges from 38% in EO27 to 48% in EO76, with mean GFS values uniformly well below the AABB benchmark (0.667). All five EOs breach the TP2 AAAA threshold: the red (AAAA) segment extends well past the 30% dashed line in every case, reaching 52–62% of each population. EO27 is the most degraded (only 38% supporting; mean GFS = 0.224); EO67 is least degraded (47% supporting; mean GFS = 0.297) and is the only EO with visible AABC and ABCD tier segments.

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

### 5. Hypothesized Processes Leading to AAAA Predominance

The finding that 56% of individuals across all five EOs carry an AAAA genotype — and that only five individuals in the entire species retain three or more distinct SRK alleles — requires mechanistic explanation. The pattern is unlikely to reflect a single event; instead, it is most consistent with a cascade of interacting demographic and genetic processes, each compounding the others over successive generations. The five stages of this cascade are summarised in [Figure 14](#figure-14) below.

**Stage 1 — Ancestral bottleneck: loss of SRK allele richness**

The most parsimonious starting point is a severe demographic bottleneck associated with historical habitat loss and fragmentation. Bottlenecks reduce allelic richness at the SRK locus faster than expected under neutral models, because S-alleles are individually rare even in healthy populations under NFDS and are easily lost when founder group size is small. The observed richness of 9–23% of the species optimum per EO, combined with the irreversible loss of 60–89% of the species allele pool per occurrence, is consistent with pronounced and prolonged founder effects (Wright 1939; Schierup et al. 1997).

**Bottleneck versus ongoing drift: what the pipeline can and cannot resolve**

The term "genetic drift" is used throughout this report in its broadest sense — encompassing both a discrete historical bottleneck (a sudden, severe reduction in population size) and the continuous stochastic erosion of allele frequencies that occurs in any finite population. These are conceptually distinct processes with different primary effects: a bottleneck principally reduces allele **richness** by eliminating rare alleles simultaneously at a single historical moment; ongoing drift in small populations subsequently erodes allele **frequency evenness**, causing surviving alleles to diverge progressively from the equal-frequency NFDS ideal.

The current pipeline measures the net outcome of both processes combined — the allele richness deficit relative to the MM species optimum — but cannot partition how much was lost in a historical bottleneck versus how much is being actively eroded today. This limitation is inherent to the SRK locus itself: it is under balancing selection (NFDS), which invalidates standard bottleneck-detection tools (e.g., BOTTLENECK software; Cornuet & Luikart 1996; Garza & Williamson 2001) that assume neutral allele dynamics. Resolving the relative contribution of each process would require genome-wide neutral marker data (e.g., SNPs from RADseq or whole-genome sequencing) analysed with demographic modelling approaches such as SMC++ or PSMC.

Nonetheless, the two TP1 axes — allele richness and frequency evenness — carry a weak but informative signal (see Section 4a). A bottleneck founder effect tends to reduce richness while leaving surviving alleles at roughly similar frequencies among founders; subsequent ongoing drift in small populations then additionally distorts those frequencies, further eroding evenness. LEPA shows both severely low richness (9–23% of species optimum) and very low evenness (Ne/N = 0.41–0.55), most consistent with a historical bottleneck followed by prolonged ongoing drift. The allele composition analysis (Figures 5–6) adds a second line of evidence: if losses were caused by a single shared ancestral bottleneck, EOs would be missing the *same* alleles. Instead, each EO holds a largely private allele set with minimal inter-EO sharing, consistent with independent drift in isolated EOs rather than one common founding event — though a prior reduction of the species-level pool cannot be excluded.

For conservation purposes, the distinction matters less than the intervention it implies: the alleles are absent regardless of when they were lost, and the urgency of inter-EO allele transfers is identical under either interpretation.

**Stage 2 — Paradox of NFDS under impoverished diversity**

Under normal conditions, NFDS protects rare S-alleles by giving individuals carrying them a reproductive advantage — they are compatible with more partners. However, this protective mechanism inverts once diversity collapses below a functional threshold. When the A haplotype becomes numerically dominant, the rare-allele advantage diminishes because too few rare alleles remain to recover. A-carrying individuals produce more offspring simply by numerical dominance, and their offspring are disproportionately likely to inherit one or more A copies. The result is a self-reinforcing collapse: drift reduces diversity, NFDS weakens, dominance of A increases, diversity collapses further (Castric & Vekemans 2004).

No universally agreed-upon threshold has been established for SSI systems, but simulation work specific to sporophytic SI indicates that the NFDS protective advantage begins to erode sharply once the dominant haplotype exceeds approximately 30–40% frequency in a population (Schierup et al. 1997). Below this point, the reproductive output gained through sheer numerical dominance outpaces the compatibility advantage conferred by rarity, and the system enters a state where drift and NFDS reinforce each other in the same direction rather than opposing each other. In LEPA, the dominant A haplotype is present in virtually every individual and AAAA individuals represent 56% of all genotypes — well past any threshold at which NFDS could act as a stabilising force. The 9–23% allele richness retained per EO places all five occurrences firmly in the regime where NFDS inversion is expected (Schierup et al. 1997; Castric & Vekemans 2004).

**Stage 3 — Mate limitation and reproductive skew**

In sporophytic SI systems, a cross is compatible only when the pistil does not recognise any haplotype expressed in the pollen donor. As AAAA frequency rises, individuals carrying rare haplotypes find progressively fewer compatible mates — not because the SI system has failed, but because it is working precisely as designed in a diversity-impoverished landscape. Byers & Meagher (1992) demonstrated analytically that in small SSI populations, mate availability collapses sharply as S-allele diversity decreases. In practice, individuals with rare haplotypes increasingly fail to set seed, removing those allele lineages from the next generation irrespective of selective pressure. Empirical confirmation of this Allee effect in SI species has been documented in *Ranunculus reptans* and invasive *Raphanus sativus* populations (Willi et al. 2005; Elam et al. 2007).

**Stage 4 — Genetic drift overwhelms NFDS at observed population sizes**

At census sizes of N = 25–40 individuals per EO, genetic drift is strong enough to overcome the frequency-dependent selection that would otherwise maintain S-allele diversity. The effective allele number (Ne) at the S-locus is theoretically elevated above neutral Ne under NFDS; however, this advantage collapses when diversity is already low. At the evenness values observed in LEPA (Ne/N = 0.41–0.55), a substantial fraction of individuals already carry duplicate allele copies, meaning effective recombination among distinct S-haplotypes is further constrained. Once a rare haplotype is lost by drift, it cannot be recovered without gene flow (Willi et al. 2005; Aguilar et al. 2006).

**Stage 5 — Polyploid-specific SI breakdown: the pathway to AAAA**

To understand this stage it helps to know the two molecular components involved. **SRK** (S-Receptor Kinase) is a protein on the surface of pistil cells that acts as the *receptor*: it detects pollen identity. **SCR** (S-locus Cysteine-Rich protein; also called SP11) is a small protein coating the surface of pollen grains that acts as the *ligand*: it identifies the pollen. When an SCR protein binds its matching SRK receptor on the pistil, the pistil recognises the pollen as self and rejects it; no binding means the pollen is accepted.

In diploid SSI plants, becoming homozygous at the S-locus is very difficult because SI prevents the same-haplotype crosses required to produce homozygous offspring. Tetraploidy fundamentally alters this constraint through the *competitive interaction model*. An AABB individual produces diploid pollen by randomly sampling two of its four allele copies. This generates three pollen types in predictable ratios: AA pollen (1 in 6), **AB pollen (4 in 6 — the most common type)**, and BB pollen (1 in 6). AB pollen, carrying both A-SCR and B-SCR proteins simultaneously, is the key to SI breakdown. When AB pollen lands on an AABB pistil — which expresses both SRK-A and SRK-B receptors — both SCR proteins attempt to bind their respective receptors at the same time. The competing signals interfere with each other, producing a recognition response too weak to trigger rejection. The pistil cannot tell the pollen is self, so it lets it through: self-fertilisation succeeds.

Self-fertilisation of an AABB plant via AB pollen produces predominantly **AAAB offspring** (e.g., AA egg × AB pollen, or AB egg × AA pollen). An AAAB individual now produces AA pollen in 3 of 6 combinations (versus 1 in 6 from AABB), which further increases the probability of SI-bypassing selfing events in the next generation. Over successive generations this runaway process drifts the genotype distribution toward AAAA. Mable (2004) and Mable et al. (2004) documented this competitive interaction mechanism in polyploid *Arabidopsis lyrata*, showing that neopolyploids frequently become partially self-compatible through dosage-mediated SCR competition. Busch & Schoen (2008) and Comai (2005) review the broader consequences of polyploidy for SI function.

**The five AABC individuals: three competing explanations**

The presence of five AABC individuals — all confined to EO67 and EO70 — is particularly informative. Three non-mutually exclusive explanations are consistent with the data:

- *Relict diversity.* These individuals are vestiges of a pre-bottleneck state, arising from crosses between individuals that still carried B, C, or D haplotypes before those haplotypes were lost from the broader population. Their AABC configuration requires two rare haplotypes from different parents — a cross that becomes increasingly improbable as the rare haplotype pool diminishes.

- *Inter-EO gene flow.* The co-occurrence of AABC genotypes exclusively in EO67 and EO70 may reflect historical gene flow between these two EOs, with the rarer haplotypes (B, C) derived from the alternate source population. Such inter-occurrence crosses would be fully compatible under SSI and would produce higher-fitness offspring, consistent with the AABC observation.

- *Mutational SI breakdown.* Rare loss-of-function mutations in SCR/SRK components can restore self-compatibility in otherwise SI species (Busch & Schoen 2008; Brennan et al. 2002). If a small number of LEPA individuals carry such mutations, they could participate in crosses normally rejected by the SI system, generating unusual allele combinations. However, AABC is more consistent with a compatible outcross than with selfing, making this the least parsimonious explanation.

**Conservation implication**

This mechanistic framework has a direct practical corollary: increasing census population size within an EO alone is insufficient to reverse the AAAA trajectory. If the prevailing AAAA frequency exceeds the reproductive dead-end threshold, within-population growth simply produces more AAAA offspring. Allele introduction via inter-EO crosses is the primary lever. Introducing B, C, and D haplotypes from the five AABC individuals into crosses with AAAA individuals simultaneously restores SI compatibility and re-engages NFDS — the self-reinforcing mechanism that, once functional, will favour the spread of introduced alleles through subsequent generations. The five AABC seed parents identified in Section 4b are therefore not merely high-priority for seed production in the current season; they represent the only endogenous genetic resource capable of reversing the collapse cascade described above.

<a name="figure-14"></a>

![Figure 14: Five-stage cascade hypothesis for AAAA predominance in LEPA](figures/SRK_AAAA_cascade_hypothesis.png)

---

### References

Aguilar, R., Ashworth, L., Galetto, L., & Aizen, M. A. (2006). Plant reproductive susceptibility to habitat fragmentation: Review and synthesis through a meta-analysis. *Ecology Letters*, 9(8), 968–980. https://doi.org/10.1111/j.1461-0248.2006.00927.x

Brennan, A. C., Harris, S. A., Tabah, D. A., & Hiscock, S. J. (2002). The population genetics of sporophytic self-incompatibility in *Senecio squalidus* L. (Asteraceae) I: S allele diversity in a natural population. *Heredity*, 89(6), 430–438. https://doi.org/10.1038/sj.hdy.6800154

Busch, J. W., & Schoen, D. J. (2008). The evolution of self-incompatibility when it is costly. *Trends in Plant Science*, 13(3), 128–135. https://doi.org/10.1016/j.tplants.2008.01.001

Byers, D. L., & Meagher, T. R. (1992). Mate availability in small populations of plant species with homomorphic sporophytic self-incompatibility. *Heredity*, 68(4), 353–359. https://doi.org/10.1038/hdy.1992.49

Castric, V., & Vekemans, X. (2004). Plant self-incompatibility in natural populations: A critical assessment of recent theoretical and empirical advances. *Molecular Ecology*, 13(10), 2873–2889. https://doi.org/10.1111/j.1365-294X.2004.02296.x

Comai, L. (2005). The advantages and disadvantages of being polyploid. *Nature Reviews Genetics*, 6(11), 836–846. https://doi.org/10.1038/nrg1711

Elam, D. R., Ridley, C. E., Goodell, K., & Ellstrand, N. C. (2007). Population size and relatedness affect fitness of a self-incompatible invasive plant. *Proceedings of the National Academy of Sciences*, 104(2), 549–552. https://doi.org/10.1073/pnas.0607259104

Mable, B. K. (2004). Polyploidy and self-compatibility: Is there an association? *New Phytologist*, 162(3), 803–811. https://doi.org/10.1111/j.1469-8137.2004.01055.x

Mable, B. K., Beland, J., & Di Berardo, C. (2004). Inheritance and dominance of self-incompatibility alleles in polyploid *Arabidopsis lyrata*. *Heredity*, 93(5), 476–486. https://doi.org/10.1038/sj.hdy.6800526

Schierup, M. H., Vekemans, X., & Christiansen, F. B. (1997). Evolutionary dynamics of sporophytic self-incompatibility alleles in plants. *Genetics*, 147(2), 835–846. https://doi.org/10.1093/genetics/147.2.835

Vekemans, X., & Slatkin, M. (1994). Gene and allelic genealogies at a gametophytic self-incompatibility locus. *Genetics*, 137(4), 1157–1165. https://doi.org/10.1093/genetics/137.4.1157

Willi, Y., Van Buskirk, J., & Fischer, M. (2005). A threefold genetic Allee effect: Population size affects cross-compatibility, inbreeding depression and drift load in the self-incompatible *Ranunculus reptans*. *Genetics*, 169(4), 2255–2265. https://doi.org/10.1534/genetics.104.034553

Wright, S. (1939). The distribution of self-sterility alleles in populations. *Genetics*, 24(4), 538–552.

---

### 6. Crossing Strategy Simulations

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

### 7. Conclusion

Genetic drift, promoted by habitat fragmentation, is severely eroding allele diversity and frequency balance in *LEPA* Element Occurrences at two distinct levels.

**Tipping Point 1 — Allele richness deficit (Figures 4b, 9).** With 189 individuals sampled across 26 population localities (five major EOs plus 21 additional small sites), only 47 of an estimated 65 species-level S-allele bins (MM estimate) have been observed. Each EO retains a small, largely private subset of the species SI repertoire (9–23% of the species optimum), allele frequencies are highly skewed from the equal distribution expected under NFDS (χ² *p* < 10⁻⁷ at every level), and frequency evenness (Ne/N) ranges from only 0.41 to 0.55 — meaning fewer than half of even the observed alleles are effectively contributing to SI function. Allele sets are almost entirely non-overlapping across EOs — only 2 alleles are shared across all five. All five EOs are flagged **CRITICAL** for TP1. EO70 is the most allele-depauperate (6 bins, 9% of optimum; evenness = 0.51); EO27 retains the most (15 bins, 23% of optimum) and the largest set of private alleles (n = 10). Critically, the allele erosion analysis (Figure 4b) confirms that these deficits are not sampling artefacts: between 60% (EO27) and 89% (EO70) of the species-level S-allele pool has been irreversibly lost from each EO through genetic drift. The predicted-undetected component is small (1–11 alleles per EO), meaning further sampling cannot close the gap — the missing alleles are genuinely absent from these populations.

**Tipping Point 2 — Genotypic fitness collapse (Figures 12–13).** Beyond allele counts, the *distribution of alleles across individuals* has degraded severely. The Genotypic Fitness Score (GFS) — the proportion of heterozygous gametes a tetraploid individual produces — reveals that all five EOs are CRITICAL: mean GFS ranges from 0.22 to 0.30 (well below the AABB benchmark of 0.667), and 52–62% of individuals per EO are AAAA (producing zero heterozygous gametes). Fewer than half of individuals in any EO carry more than one distinct SRK allele and can therefore support reproductive effort through compatible crosses (Figure 13): the proportion of supporting individuals (GFS > 0) ranges from only 38% in EO27 to 48% in EO76. No ABCD individual exists in the dataset. Only five individuals across the entire species carry an AABC genotype (GFS = 0.833) — all in EO67 and EO70 — making them the highest-priority seed parents for near-term managed crossing.

Critically, these two tipping points interact: even if allele richness were restored through inter-EO transfers (TP1 intervention), the benefit would be limited if the incoming alleles are absorbed into AAAA or AAAB individuals. Effective restoration therefore requires simultaneously targeting allele richness (inter-EO transfers of rare alleles) and genotype quality (crosses designed to produce AABB, AABC, and ultimately ABCD offspring).

Crossing simulations demonstrate that managed, optimisation-guided crosses can reduce allele frequency variance by 77–95% within a single generation compared to random mating, with zero allele loss under preservation strategies. Seed production efforts should prioritise: (1) AABC and AABB individuals as seed parents, ranked by GFS within each EO; (2) inter-EO crosses targeting rare-allele carriers to address the richness deficit; and (3) crossing designs that pair AAAB individuals carrying complementary allele pairs, which can directly produce higher-GFS offspring in the next generation.
