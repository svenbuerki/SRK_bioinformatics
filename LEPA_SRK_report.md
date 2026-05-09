# SRK-Based Assessment of Self-Incompatibility in *Lepidium papilliferum* (LEPA)

## Executive Summary

**The conservation challenge.** *Lepidium papilliferum* (LEPA) is a self-incompatible (SI) tetraploid plant restricted to fragmented populations occupying slickspots in southwestern Idaho (USA) (Buerki et al., 2026). This report uses an SRK haplotype dataset of **272 individuals across 19 Element Occurrences** to ask: *what is the impact of habitat fragmentation and genetic drift on the species' SI system, and how can we design a genetically informed breeding programme to recover reproductive fitness?* The analysis is structured around six questions, each with a defensible biological headline.

**Q1 — Habitat fragmentation and bottlenecks.** Spatial-connectivity analysis (sibling repository [LEPA_EO_spatial_clustering](https://github.com/svenbuerki/LEPA_EO_spatial_clustering)) identifies **five independent bottleneck lineages (BL1–BL5)** across the species range with **no between-EO pollinator connections anywhere in the dataset**, establishing the demographic strata used throughout the report.

**Q2 — Species S-allele richness.** 54 distinct S-alleles observed; Michaelis-Menten estimate = 69. Per-BL richness varies dramatically: **BL4 retains 31 alleles (57% — the diversity reservoir)** while BL1 and BL2 have lost 86–90% of the species pool. **33 of 54 alleles (61%) are private to a single BL** — direct empirical evidence of multiple independent founder events, not a single shared bottleneck.

**Q3 — Drift on populations' diversity (TP1).** All five BLs and 5 of 6 plotted EOs are flagged **CRITICAL** (richness < 50% of optimum AND frequency evenness Ne/N < 0.80). No lineage retains a balanced allele pool — restoration requires inter-BL allele transfers; within-lineage breeding alone cannot rebuild the SI system.

**Q4 — Drift on individual reproductive fitness (TP2).** All five BLs CRITICAL on TP2 (mean Genotypic Fitness Score < 0.667 AND > 30% AAAA homozygotes). Allele_050 and Allele_057 (Synonymy group 1) are pan-BL fixed in AAAA individuals — **convergent drift onto the same SI specificity despite independent bottlenecks**. **13 AABC individuals** species-wide carry the heterozygous-gamete potential needed for high-yield managed crossing — they are the immediate seed-parent priority.

**Q5 — Mechanism: convergent S-allele depletion.** A **five-stage Compatibility Collapse Cascade (C3)** explains how habitat fragmentation produces 60% reproductive dead-ends within a structurally intact SI system. Cross-Brassicaceae per-BL entropy decomposition reveals that the residue identities at LEPA hypervariable (HV) columns are **74% LEPA-specific** (different from Brassica AND Arabidopsis dominant residues), confirming the convergent allele depletion is **drift on a shared LEPA ancestral pool**, not pan-Brassicaceae selection convergence.

**Q6 — Cross-based hypothesis testing.** A **101-cross phased experimental plan (1 323 attempts)** validates the bioinformatic predictions. The plan is structured around five nested hypotheses (H0 SI validation; H1a/H1b compatibility baselines; H2 synonymy bin boundaries; H3 hidden bins via heterozygous donors with paired controls). Outcomes feed into a **validated functional S-allele table** that informs operational seed-orchard design.

**Action items.** (1) Cross all 13 AABC seed parents this season; (2) prioritise inter-BL allele transfers (BL4 → other BLs); (3) experimentally test Synonymy group 1 first (13 HV-identical alleles, 103 AAAA individuals — incompatible outcome would consolidate 40 % of all sequence bins into a single functional specificity).

---

## Background

Self-incompatibility (SI) in *L. papilliferum* is controlled by the S-locus, where the extracellular S-domain of the S-receptor kinase (SRK) protein acts as the female determinant of pollen rejection. Each individual carries up to four allele copies (tetraploid), and compatible mating requires that pollen and pistil carry different SRK alleles. Throughout this report, each functionally distinct SRK protein variant is referred to as an **S-allele**, and the tetraploid combination of S-alleles an individual carries is its **genotype** (e.g., `AAAA`, `AABB`, `AABC`). Under balancing selection — specifically, negative frequency-dependent selection (NFDS) — all S-alleles are maintained at approximately equal frequencies, maximising the proportion of compatible mating pairs. In small or isolated populations, however, genetic drift counteracts balancing selection, reducing allele richness, skewing allele frequencies, and degrading individual genotype quality, with direct consequences for reproductive success.

The analysis below addresses six sequential questions:

1. [**Q1 — How has habitat fragmentation produced demographic bottlenecks that shape the SRK self-incompatibility system in *L. papilliferum*?**](#q1-habitat-fragmentation-and-the-bottleneck-framework)
2. [**Q2 — What is the S-allele richness of *L. papilliferum*, and how is it distributed across bottleneck lineages?**](#q2-species-s-allele-richness-and-bl-distribution)
3. [**Q3 — How has genetic drift eroded S-allele diversity within populations and lineages? (Tipping Point 1)**](#q3-genetic-drift-on-populations-s-allele-diversity-tp1)
4. [**Q4 — How has genetic drift degraded the reproductive fitness of individuals within populations and lineages? (Tipping Point 2)**](#q4-genetic-drift-on-individual-reproductive-fitness-tp2)
5. [**Q5 — What mechanism explains the convergent S-allele depletion observed across all five bottleneck lineages of *L. papilliferum*?**](#q5-mechanism-convergent-s-allele-depletion)
6. [**Q6 — How can S-allele specificity hypotheses be tested with controlled crosses to support a genetically informed breeding programme?**](#q6-cross-based-hypothesis-testing-for-informed-breeding)

The bioinformatic pipeline producing the underlying data is documented in [Bioinformatics_pipeline.md](Bioinformatics_pipeline.md). A complete methodological reference, including reproducible commands and parameter tables, is in [index.Rmd](https://svenbuerki.github.io/SRK_bioinformatics/).

---

## Q1 — Habitat fragmentation and the bottleneck framework {#q1-habitat-fragmentation-and-the-bottleneck-framework}

### Why this question matters

Genetic drift acts within demographically isolated populations. Before evaluating how drift has impacted S-allele diversity (Q3) and individual reproductive fitness (Q4), we must establish how the species is spatially structured into demographically distinct units. This question defines the inferential strata used throughout the rest of the report.

### Method

The spatial-connectivity analysis was performed in the sibling GitHub repository [**LEPA_EO_spatial_clustering**](https://github.com/svenbuerki/LEPA_EO_spatial_clustering). It (1) constructed convex-hull polygons from 758 georeferenced collection events across 39 locations in 19 EOs (UTM Zone 11N); (2) partitioned the locations into geographic groups using the maximum pollinator dispersal distance for *L. papilliferum* (~500 m) as the connectivity threshold; (3) clustered the 32 resulting group centroids by Ward's D2 hierarchical clustering of geographic distances; and (4) selected k = 5 by silhouette optimisation (silhouette score = 0.73), yielding **five independent bottleneck lineages (BL1–BL5)** ([Figure 1](#figure-1)). A complementary drift index (DI) was derived from the union area of connected-group hull polygons (DI = 0: largest group, weakest drift; DI = 1: smallest, expected-strongest drift) as a relative spatial proxy for *Ne* and is visualised lineage-by-lineage in [Figure 2](#figure-2).

The BL × group × EO cross-reference is published as `tables/EO_group_BL_summary.csv` in the spatial-clustering repository and propagated to this dataset as the input to Step 13 of the SRK pipeline (`SRK_BL_integration.py`), which writes per-individual BL assignments consumed by every downstream analysis.

<a name="figure-1"></a>

![Figure 1: Ward's D2 hierarchical clustering of the 32 geographic groups (across 19 EOs of *L. papilliferum*) partitions the species into five independent bottleneck lineages (BL1–BL5; silhouette-optimal k = 5; silhouette score = 0.73). Lower strip: per-group drift index (red = strong drift / DI > 0.75; blue = weak drift / DI < 0.5). Bottom strip: census size N per group. Set1 palette matches the BL colour scheme used throughout this report. Source: `EO_clustering_dendrogram.png` from the [LEPA_EO_spatial_clustering](https://github.com/svenbuerki/LEPA_EO_spatial_clustering) repository.](figures/EO_clustering_dendrogram.png)

<a name="figure-2"></a>

![Figure 2: Predicted genetic drift intensity across the five independent bottleneck lineages of *L. papilliferum*. One panel per BL (BL1–BL5, in lineage order); each panel shows the geographic groups within that lineage as polygons coloured by drift index (red = strong drift / small habitat footprint; blue = weak drift / large footprint), sized by census population size, with within-BL pollinator connections shown as solid lines (≤ 500 m) and between-group separations as dashed lines (> 500 m). The lineage-level panels show how isolation and drift intensity are distributed within each BL — for example, BL4 contains the only large-area / low-DI group in the species (EO27 group 11 at 19.6 ha, DI = 0.000), while BL1 is dominated by extreme-DI singletons. Source: `EO_BL_drift_panel.png` from the [LEPA_EO_spatial_clustering](https://github.com/svenbuerki/LEPA_EO_spatial_clustering) repository.](figures/EO_BL_drift_panel.png)

### Bottleneck-lineage membership and SRK sample sizes

(262 of 272 BL-assigned ingroup individuals; 10 germplasm sub-codes remain Unassigned and contribute only to species-level baselines.)

| BL  | N geographic groups | N locations | N EOs | N individuals (SRK) | N alleles observed | % of 54-allele species pool retained |
|-----|--------------------:|------------:|------:|--------------------:|-------------------:|-------------------------------------:|
| BL1 | 10                  | 12          | 4     | 5                   | 4                  | 7%                                   |
| BL2 | 2                   | 2           | 2     | 51                  | 8                  | 15%                                  |
| BL3 | 4                   | 4           | 4     | 68                  | 16                 | 30%                                  |
| BL4 | 7                   | 8           | 4     | 78                  | **31**             | **57%**                              |
| BL5 | 9                   | 13          | 5     | 60                  | 18                 | 33%                                  |

### Key findings

- **Connectivity:** of 741 pairwise location comparisons, only 7 pairs (0.9 %) are connected within the 500 m pollinator dispersal threshold — and **all 7 are within the same EO**. There are *no between-EO pollinator connections anywhere in the dataset* under current conditions: every EO is a closed demographic unit at the pollination scale.
- **Habitat footprint:** ≈ 84 % of the 32 geographic groups occupy < 1 ha (DI > 0.95), placing the vast majority of locations in the extreme-drift regime. Only EO27 group 11 (19.6 ha, DI = 0.000), EO32 (14.6 ha), and EO18 (8.1 ha as a fully-connected 5-location chain) escape extreme spatial drift.
- **Five lineages, independent histories:** the BL framework provides three analytical advantages — it rescues 21 small localities (n = 1–3 each) that cannot support EO-level statistics on their own; it allows direct empirical testing of the independent-bottleneck hypothesis (confirmed in Q3 by allele-sharing analyses); and it preserves the full 272-individual species pool for any species-level baseline by stratifying only at the population/lineage level.

The locked Set1 colour palette (BL1 red `#E41A1C`, BL2 blue `#377EB8`, BL3 green `#4DAF4A`, BL4 purple `#984EA3`, BL5 orange `#FF7F00`) is used consistently across this report and the sibling spatial-clustering repository so that the dendrogram, drift panels, and every BL-stratified figure cross-reference visually.

---

## Q2 — Species S-allele richness and BL distribution {#q2-species-s-allele-richness-and-bl-distribution}

### Why this question matters

The species-level S-allele pool approximates the balancing-selection richness equilibrium and serves as the reference baseline against which all population-level deficits in Q3 and Q4 are measured. Per-BL stratification of the same data tests the independent-bottleneck hypothesis directly: if all five lineages had derived from a single shared ancestral bottleneck, each would have lost approximately the same alleles. If drift acted independently in each BL, each lineage's surviving allele set should be largely private.

### Method

S-alleles were defined by distance-based clustering of validated SRK protein sequences on the S-domain ectodomain (Step 10 of the bioinformatic pipeline), with a sensitivity-analysis-derived p-distance threshold of 0.005 (0.5 %). Each cluster represents an *S-allele bin* — a sequence-based hypothesis that proteins within the bin share recognition specificity. Allele accumulation curves were fit at species level (all 272 ingroup individuals) and at BL level (262 BL-assigned individuals) using rarefaction with Michaelis-Menten (MM), Chao1, and iNEXT asymptote estimators (Step 15).

<a name="figure-3"></a>

![Figure 3: Species-level S-allele accumulation curve. The curve has not yet reached an asymptote, indicating that further sampling would discover additional alleles. The Michaelis-Menten estimate of 69 alleles is adopted as the species optimum.](figures/SRK_allele_accumulation_species.png)

<a name="figure-4"></a>

![Figure 4: S-allele accumulation curves per bottleneck lineage. Each curve is coloured by parent BL (Set1 palette: BL1 red, BL2 blue, BL3 green, BL4 purple, BL5 orange — matching the LEPA_EO_spatial_clustering project). Species MM = 69 shown as a dashed grey reference line; BL-specific MM estimates as dashed coloured lines per curve. BL4 is the diversity reservoir (31 alleles observed, MM = 47); BL1 and BL2 have lost 86–90% of the species pool.](figures/SRK_allele_accumulation_BL_combined.png)

### Key findings

Across **272 individuals** sampled from 26 population localities, **54 distinct S-allele bins** were identified ([Figure 3](#figure-3)). Three asymptote estimators agree on the upper end:

| Estimator | Predicted species richness |
|-----------|---------------------------|
| Michaelis-Menten (MM) | **69 alleles** |
| Chao1 | 73 alleles |
| Consensus (mean of MM + Chao1) | 71 alleles |

The MM estimate of 69 alleles is adopted as the species optimum — the allele richness expected under balancing selection at evolutionary equilibrium, used as the reference baseline for Q3 and Q4. An additional ~49 individuals would need to be sampled to discover the next new allele, and reaching 80 % of the MM estimate would require sampling ~55 more individuals. The species SI repertoire remains substantially under-characterised.

**BL stratification reveals the diversity reservoir.** When the same individuals are partitioned into the five independent bottleneck lineages (Q1), S-allele richness varies dramatically across lineages despite comparable sampling effort ([Figure 4](#figure-4); see Q1 table for sample sizes):

- **BL4 acts as the diversity reservoir of the species** — its 78 individuals retain 31 of the 54 observed alleles (57 %), and its lineage-level MM asymptote (47) approaches the species pool (69). Further sampling within BL4 would still discover many additional alleles.
- The four other BLs have each collapsed to a fraction of their lineage-level potential: BL1 and BL2 have lost 86–90 % of the species pool to drift; BL3 and BL5 ~70 %.
- The > 2-fold gradient of richness loss across BLs can only arise if **the bottlenecks operating in each lineage have proceeded independently**: a single shared species-level bottleneck would predict comparable richness loss across BLs.

This independent-bottleneck signature is reinforced by the allele-sharing analyses in Q3.

---

## Q3 — Genetic drift on populations' S-allele diversity (TP1) {#q3-genetic-drift-on-populations-s-allele-diversity-tp1}

### Why this question matters

Tipping Point 1 (TP1) assesses the **health of the self-incompatibility system** at the population/lineage level — the degree to which the S-allele pool is capable of sustaining compatible mating. It is breached when allele loss is so severe that inter-population allele transfers are required to restore SI function.

### Method

Two complementary axes structure the assessment:

- **How many different alleles has a population retained?** (`prop_optimum` = N_alleles / 69 — the proportion of the species-level SI repertoire still present). A population holding all alleles can offer every individual a large pool of compatible partners; as alleles are lost, compatible pairings become progressively rarer.
- **How evenly are the remaining alleles distributed across individuals?** (`Ne / N_alleles` — the ratio of effective to observed allele number, where Ne = 1/Σpᵢ². A ratio of 1.0 means perfect evenness — the balancing-selection ideal in which every allele contributes equally to compatible crosses; drift and dominance push it downward.)

A population breaching both axes simultaneously (`prop_optimum < 0.50` AND `Ne/N < 0.80`) is flagged **CRITICAL**; one criterion is **AT RISK**; neither is **OK**. Both EO-level and BL-level results are computed in parallel.

A complementary χ² goodness-of-fit test against the equal-frequency NFDS expectation is computed at three levels (species, EO, BL).

### Key findings

**The raw drift signal — alleles lost per Element Occurrence ([Figure 5](#figure-5)) and per Bottleneck Lineage ([Figure 6](#figure-6)).** Before synthesising both axes into the TP1 diagnostic, it is worth visualising the raw drift signal directly. For each EO and each BL, we partition the deficit relative to the 69-allele species optimum (MM) into two components: alleles predicted to exist in the group but not yet detected (light blue; group-MM minus observed), and alleles **lost to genetic drift** (red; species-MM minus group-MM). The red component dominates catastrophically at both stratification levels.

<a name="figure-5"></a>

![Figure 5: S-allele erosion by genetic drift per Element Occurrence. For each EO, the bar height represents the species optimum (69 alleles); segments decompose this into observed alleles (dark blue), predicted-undetected alleles (light blue, group-MM minus observed), and alleles lost to genetic drift (red, species-MM minus group-MM). Between 45% (EO27) and 90% (EO70) of the species-level S-allele pool has been irreversibly lost from each EO. The predicted-undetected component is small (1–16 alleles), confirming that further sampling will not close the gap — the missing alleles are genuinely absent from these populations.](figures/SRK_allele_accumulation_drift_erosion.png)

<a name="figure-6"></a>

![Figure 6: S-allele erosion by genetic drift per Bottleneck Lineage. Same decomposition as Figure 5, aggregated by BL. The BL color strip below the bars matches the lineage palette used elsewhere in the report. BL4 (purple) retains the most alleles (31, 57% of species pool) but still shows substantial drift loss; BL1 and BL2 have lost 86–90% of the species pool to drift.](figures/SRK_allele_accumulation_BL_drift_erosion.png)

**Per-EO drift loss** ranges from 45% (EO27, the least eroded) to 90% (EO70, the most eroded). Even in EO27 — the most allele-rich EO — an estimated 31 of the 69 species-level S-allele bins have been permanently lost from the local gene pool. EO70 has lost ≈ 62 of 69 alleles. **These deficits are not sampling artefacts**: the predicted-undetected component is small (1–16 alleles per EO), confirming that further sampling within these EOs cannot close the gap.

**Per-BL drift loss** ranges from 32% (BL4, the diversity reservoir) to 86–90% (BL1, BL2). Even BL4's lineage-level MM asymptote (47) falls well short of the species pool (69), and BL1/BL2 have collapsed to a small fraction of their lineage-level potential. The > 2-fold gradient of richness loss across BLs is the **independent-bottleneck signature** documented in Q2: a single shared species-level bottleneck would predict comparable richness loss across BLs.

**Synthesis: the TP1 tipping point ([Figure 7](#figure-7)).** Combining the richness deficit with the frequency-evenness axis produces the TP1 diagnostic: a single scatter that classifies every EO and every BL as CRITICAL / AT RISK / OK based on whether each axis is breached.

<a name="figure-7"></a>

![Figure 7: TP1 tipping point — health of the SI system. Six EOs (circles) and 5 BL aggregates (triangles) plotted on the same scatter, both coloured by parent BL using the locked Set1 palette. All five BLs and five of six EOs fall in the CRITICAL zone (lower-left, both prop_optimum < 0.50 AND evenness < 0.80). EO18 (n = 5) escapes CRITICAL only because its small sample size inflates the evenness ratio; its richness deficit is among the most severe. Notably, BL4 — the species' diversity reservoir at 45% of optimum — is CRITICAL through the evenness axis (0.33), demonstrating that drift is operating along two independent dimensions even where allele richness is best preserved.](figures/SRK_TP1_tipping_point.png)

**All five BLs and 5 of 6 plotted EOs are CRITICAL on TP1** ([Figure 7](#figure-7)). Under the BL re-frame, TP1 reveals that **fragmentation has imposed a uniform CRITICAL status across all five bottleneck lineages**: BL1 and BL2 are CRITICAL through richness loss (≥ 86 % of the species pool absent); BL3, BL4 and BL5 are CRITICAL through frequency skew (Ne/N ≤ 0.58). The two axes capture different signatures of drift, but every BL fails on at least one — and most fail on both.

**This is the single most consequential finding of the population-genetic analysis: contemporary recovery requires inter-BL allele transfers because no lineage retains a balanced allele pool that within-lineage crossing alone could draw upon.**

**Allele-sharing patterns directly confirm the independent-bottleneck hypothesis.** If all five lineages had derived from a single shared species-level bottleneck, we would expect overlapping losses — every BL missing roughly the same alleles. Instead, the BL UpSet plot ([Figure 8](#figure-8)) and the pairwise heatmap ([Figure 9](#figure-9)) reveal the opposite pattern:

| Bottleneck lineage | N alleles observed | Private alleles | % private |
|---|---:|---:|---:|
| BL1 |  4 |  0 |  0% |
| BL2 |  8 |  2 | 25% |
| BL3 | 16 | 11 | 69% |
| BL4 | 31 | **14** | **45%** |
| BL5 | 18 |  6 | 33% |

**33 of 54 alleles (61 %) are present in only one BL.** Pairwise BL allele sharing ([Figure 9](#figure-9)) is highest between BL4 and BL5 (12 alleles), suggestive of limited historical contact between these two lineages, but otherwise BLs share between 2 and 6 alleles in each pairwise comparison — a near-disjoint partitioning that can only arise from independent demographic histories. The two alleles shared across all five BLs (Allele_050, Allele_057) likely represent the original species-wide pool that survived in every lineage by virtue of high ancestral frequency. The remaining 52 alleles either sit in private compartments per BL or are shared among at most two or three lineages.

<a name="figure-8"></a>

![Figure 8: S-allele sharing among Bottleneck Lineages (UpSet plot). Single-BL bars are coloured by BL using the locked Set1 palette (BL1 red, BL2 blue, BL3 green, BL4 purple, BL5 orange); multi-BL intersections are grey. The first three bars (private alleles in BL4, BL3, BL5) account for the bulk of the species' allele diversity and are absent from every other lineage.](figures/SRK_allele_upset_BLs.png)

<a name="figure-9"></a>

![Figure 9: Pairwise S-allele sharing heatmap between Bottleneck Lineages. Diagonal cells are total allele richness per BL (BL4 = 31, BL5 = 18, BL3 = 16, BL2 = 8, BL1 = 4); off-diagonal cells are shared-allele counts. The highest off-diagonal value (BL4 ↔ BL5 = 12) hints at limited historical contact between these two lineages.](figures/SRK_allele_sharing_heatmap_BLs.png)

**Frequency-distribution test of NFDS.** A χ² goodness-of-fit test of allele copy-count frequencies against the equal-frequency NFDS expectation confirms that drift has skewed allele frequencies away from the NFDS equilibrium at every analysis level: at the species level (χ² = 2438.84, df = 53, p ≈ 0), in every medium and large EO (all N ≥ 36 reject at p < 1 × 10⁻⁷), and in every BL with statistical power (BL2: χ² = 175, p < 1 × 10⁻³³; BL3: χ² = 202; BL4: χ² = 325 — the largest test statistic of any group despite holding the most alleles; BL5: χ² = 95). Notably **BL4's diversity is undermined by within-lineage frequency skew**, demonstrating that S-allele erosion can proceed via two complementary axes (loss of richness and frequency distortion of remaining alleles) and that even the diversity reservoir is not exempt from active drift.

---

## Q4 — Genetic drift on individual reproductive fitness (TP2) {#q4-genetic-drift-on-individual-reproductive-fitness-tp2}

### Why this question matters

Even where some allele diversity persists (Q3), drift causes allele copy-count imbalances at the *individual* level that directly reduce reproductive output. In a tetraploid, the proportion of compatible diploid gametes an individual produces — the **Genotypic Fitness Score (GFS)** — depends not only on which alleles are present but on their dosage balance across the four allele copies. Tipping Point 2 (TP2) marks the threshold at which this individual-level degradation is so advanced that within-population crosses alone cannot restore reproductive fitness.

### Method — Genotypic Fitness Score (GFS) and TP2

A tetraploid produces diploid gametes by sampling 2 of its 4 allele copies, yielding C(4,2) = 6 equally probable combinations. GFS is the proportion of those combinations that carry two distinct alleles:

$$\text{GFS}_i = 1 - \frac{\sum_k n_k\,(n_k - 1)}{12}$$

where $n_k$ is the copy number of allele $k$ and the denominator normalises to the tetraploid gamete space:

| Genotype | GFS | Heterozygous gametes |
|----------|-----|----------------------|
| ABCD | 1.000 | 6 / 6 |
| AABC | 0.833 | 5 / 6 |
| AABB | 0.667 | 4 / 6 |
| AAAB | 0.500 | 3 / 6 |
| AAAA | 0.000 | 0 / 6 |

TP2 is breached when (i) `mean GFS < 0.667` (the average individual has less reproductive capacity than an AABB genotype) AND (ii) `proportion AAAA > 0.30` (more than 30 % of individuals are reproductive dead-ends, producing only homotypic gametes). A group breaching both is flagged **CRITICAL**; one criterion is **AT RISK**; neither is **OK**.

### Key findings

**The raw reproductive-effort signal — what fraction of individuals can support compatible crosses, per EO ([Figure 10](#figure-10)) and per BL ([Figure 11](#figure-11)).** Before synthesising both axes into the TP2 diagnostic, it is worth visualising the underlying signal directly. For each EO and each BL, the proportional bar shows the GFS-tier composition of individuals: the red AAAA segment (GFS = 0) marks reproductive dead-ends; the orange/yellow/blue segments (AAAB through ABCD) mark individuals capable of contributing allelic diversity to compatible crosses.

<a name="figure-10"></a>

![Figure 10: Reproductive effort support per Element Occurrence (sorted by parent BL, BL-coloured y-axis labels). Each row decomposes the EO's individuals by GFS tier (red AAAA = dead-end through blue ABCD = maximally diverse). Dashed line marks the TP2 AAAA threshold (30%) — the red AAAA segment exceeds 30% in every plotted EO.](figures/SRK_GFS_reproductive_effort_EO.png)

<a name="figure-11"></a>

![Figure 11: Reproductive effort support per Bottleneck Lineage (sorted by mean GFS, worst at top). Y-axis labels coloured by BL (Set1 palette). Dashed line marks the TP2 AAAA threshold (30%) — the red AAAA segment exceeds 30% in every BL.](figures/SRK_GFS_reproductive_effort_BL.png)

**Fewer than half of individuals in any EO or BL carry more than one distinct SRK allele.** At the EO level ([Figure 10](#figure-10)), the proportion of "supporting" individuals (GFS > 0) ranges from 20 % in EO18 (n = 5) to 47 % in EO67 and EO25, with mean GFS values uniformly well below the AABB benchmark. At the BL level ([Figure 11](#figure-11)), the proportion ranges from 34 % in BL3 to 60 % in BL1 (small-N caveat). The remainder are reproductive dead-ends (AAAA, GFS = 0).

**Synthesis: the TP2 tipping point ([Figure 12](#figure-12)).** Combining mean GFS with the proportion of AAAA individuals produces the TP2 diagnostic — a single scatter that classifies every EO and every BL as CRITICAL / AT RISK / OK based on whether each axis is breached.

<a name="figure-12"></a>

![Figure 12: TP2 tipping point — mean GFS vs proportion AAAA. EOs as circles, BLs as triangles, all coloured by parent BL using the Set1 palette. All five BLs are CRITICAL (mean GFS < 0.667 AND > 30% AAAA); 5 of 6 plotted EOs are CRITICAL (EO67 is AT RISK, breaching only the AAAA threshold).](figures/SRK_GFS_plots_p3_TP2_tipping_point.png)

**All five BLs are CRITICAL on TP2** ([Figure 12](#figure-12)). The lineage-level pattern is robust: even when 262 BL-assigned individuals are pooled into independent bottleneck lineages, every lineage exceeds 30 % AAAA and falls well below the AABB-benchmark mean GFS of 0.667.

| BL | N | mean GFS | % AAAA | TP2 status |
|----|:---:|:---:|:---:|:---:|
| BL3 | 68 | 0.196 | 66% | **CRITICAL** |
| BL2 | 51 | 0.235 | 63% | **CRITICAL** |
| BL4 | 78 | 0.267 | 59% | **CRITICAL** |
| BL5 | 60 | 0.272 | 55% | **CRITICAL** |
| BL1 |  5 | 0.367 | 40% | **CRITICAL** |

EO-level results show 5 of 6 plotted EOs CRITICAL (EO67 is AT RISK, breaching only the AAAA threshold). **EO67 is the least degraded** EO with the highest mean GFS (0.319) and the highest proportion of AABC individuals (14 %).

**The AAAA majority is dominated by just two alleles species-wide — a pan-BL convergent fixation.** Allele_050 and Allele_057 together account for 24–69 % of AAAA individuals per BL ([Figure 13](#figure-13)) and are present as AAAA homozygotes in **every BL** (pan-BL, 5/5). This convergent fixation across five independently-bottlenecked lineages is the single most important biological result of the TP2 analysis — every lineage has independently fixed the same dominant allele family. If the synonymy hypothesis is confirmed by crossing (Q6), then the entire AAAA fraction species-wide represents a small number of effectively-identical functional locks, dramatically reducing the practical breeding pool.

<a name="figure-13"></a>

![Figure 13: Allele identity of AAAA individuals per Bottleneck Lineage. Allele_050 (orange) and Allele_057 (blue) — both members of Synonymy group 1 (Q6) — are pan-BL across all 5 BLs. The 16% pan-Brassicaceae conserved residues vs the 74% LEPA-specific residues at HV columns confirm this is shared-ancestry drift, not pan-genera selection (Q5).](figures/SRK_GFS_AAAA_allele_composition_BL.png)

**Top seed-parent priorities (13 AABC individuals, GFS = 0.833):**

| EO | Individual | Genotype | GFS |
|----|-----------|----------|-----|
| EO67 | Library001_barcode02 | AABC | 0.833 |
| EO67 | Library006_barcode59 | AABC | 0.833 |
| EO67 | Library007_barcode82 | AABC | 0.833 |
| EO70 | Library008_barcode01 | AABC | 0.833 |
| EO70 | Library008_barcode12 | AABC | 0.833 |
| EO27 | Library006_barcode37 | AABC | 0.833 |
| EO27 | Library006_barcode41 | AABC | 0.833 |
| EO27 | Library009_barcode82 | AABC | 0.833 |
| EO25 | Library009_barcode62 | AABC | 0.833 |
| EO76 | Library001_barcode09 | AABB | 0.667 |
| EO76 | Library001_barcode17 | AABB | 0.667 |
| EO76 | Library003_barcode66 | AABB | 0.667 |

Full ranked lists per EO are in `SRK_individual_GFS.tsv`; EO-level summaries and TP2 flags in `SRK_EO_GFS_summary.tsv`.

**TP1 and TP2 interact.** Even if allele richness were restored through inter-BL transfers (Q3 intervention), the benefit would be limited if incoming alleles were absorbed into AAAA or AAAB individuals. Effective restoration therefore requires simultaneously targeting allele richness (inter-BL transfers of rare alleles) AND genotype quality (crosses designed to produce AABB, AABC, and ABCD offspring).

---

## Q5 — Mechanism: convergent S-allele depletion {#q5-mechanism-convergent-s-allele-depletion}

### Why this question matters

Q1–Q4 document the *outcomes* of habitat fragmentation and drift on LEPA's SI system. The central question this section addresses is *how* — what mechanism explains both (i) the severe allele depletion observed in every lineage AND (ii) the convergent fixation of the same Synonymy group 1 alleles across five independently-bottlenecked lineages?

Two complementary lines of evidence are presented: a **mechanistic hypothesis** (the Compatibility Collapse Cascade, C3) explaining the demographic-genetic causal chain, and a **direct empirical test** (cross-Brassicaceae per-BL entropy decomposition) distinguishing drift on a shared LEPA ancestral pool from pan-Brassicaceae selection convergence.

### The Compatibility Collapse Cascade (C3) hypothesis

A critical baseline observation frames the entire mechanism: across all five EOs, functional SRK sequences were successfully recovered from all but five individuals in the dataset (< 4 % showing molecular SI failure). This rules out widespread loss-of-function mutation as the primary driver of the 60 % AAAA prevalence and confirms that the SI machinery itself remains structurally intact.

If the SI system is functional, how does a species accumulate 60 % reproductive dead-ends? The pattern is most consistent with a five-stage cascade of interacting demographic and genetic processes — each initiated by the stage above it and each amplifying the next ([Figure 14](#figure-14)):

<a name="figure-14"></a>

![Figure 14: The Compatibility Collapse Cascade (C3): a five-stage mechanism linking habitat fragmentation to reproductive failure in LEPA. Solid down-arrows = sequential causation; dashed feedback arc = Stage 5 → Stage 2 (AAAA accumulation worsens balancing-selection breakdown, making the cascade self-propelling once initiated).](figures/SRK_AAAA_cascade_hypothesis.png)

- **Stage 1 — Ancestral bottleneck: loss of S-allele richness.** Severe demographic bottlenecks associated with historical habitat loss reduce S-allele richness faster than expected under neutral models, because S-alleles are individually rare even in healthy populations under balancing selection and are easily lost when founder group size is small. The spatial-connectivity analysis (Q1) establishes the landscape context: 26 of 39 sampled locations (67 %) have no neighbour within the 500 m pollinator dispersal distance, and ≈ 84 % of geographic groups occupy < 1 ha. Each EO is an independent evolutionary unit in which S-allele erosion has proceeded in isolation. Five independent bottleneck lineages match the predicted independent founding events (Wright 1939; Schierup et al. 1997).

- **Stage 2 — Balancing selection breaks down.** Under normal conditions, balancing selection protects rare S-alleles by giving individuals carrying them a reproductive advantage. This protection erodes sharply once the dominant S-allele exceeds ≈ 30–40 % frequency (Schierup et al. 1997). In LEPA, AAAA individuals represent 60 % of all genotypes — well past any threshold at which balancing selection could act as a stabilising force.

- **Stage 3 — Mate limitation and reproductive skew.** As AAAA frequency rises, individuals carrying rare S-alleles find progressively fewer compatible mates — not because the SI system has failed, but because it is working precisely as designed in a diversity-impoverished landscape (Byers & Meagher 1992). Lineages carrying rare S-alleles increasingly fail to set seed, removing those allele lineages from the next generation (the Allee effect documented in *Ranunculus reptans* and *Raphanus sativus*; Willi et al. 2005; Elam et al. 2007).

- **Stage 4 — Genetic drift overwhelms balancing selection at observed population sizes.** At census sizes of N = 36–62 individuals per population, drift is strong enough to overcome the balancing selection that would otherwise maintain S-allele diversity. Once a rare S-allele is lost, it cannot be recovered without gene flow (Willi et al. 2005; Aguilar et al. 2006).

- **Stage 5 — Polyploid-specific SI breakdown — the pathway to AAAA.** In diploid SI plants, becoming homozygous at the S-locus is very difficult because SI prevents same-S-allele crosses. Tetraploidy alters this through the *competitive interaction model*. An AABB individual produces diploid pollen by sampling two of its four allele copies, generating three pollen types: AA (1/6), **AB (4/6 — most common)**, and BB (1/6). AB pollen carries both A-SCR and B-SCR proteins simultaneously; on an AABB pistil expressing both SRK-A and SRK-B receptors, the competing signals interfere with each other, producing a recognition response too weak to trigger rejection. Self-fertilisation succeeds, producing predominantly AAAB offspring. AAAB then produces AA pollen in 3/6 combinations (vs 1/6 from AABB), increasing the probability of SI-bypassing selfing in the next generation. Over successive generations this runaway process drifts the genotype distribution toward full homozygosity (AAAA). Mable (2004) and Mable et al. (2004) documented this competitive interaction mechanism in polyploid *Arabidopsis lyrata*.

**The self-reinforcing loop: Stage 5 → Stage 2.** The five stages are not strictly linear — Stage 5's output (AAAA accumulation) re-enters at Stage 2, further skewing S-allele frequencies and further eroding balancing selection. This feedback explains why the cascade is increasingly irreversible over time and why intervention must target the genetic system — not merely population size — to break the loop.

### Direct empirical test — drift on a shared LEPA ancestral pool, NOT pan-Brassicaceae selection

The C3 hypothesis predicts that LEPA's convergent allele depletion reflects independent drift in each BL acting on a shared ancestral pool — not selection convergence across genera. This prediction is directly testable by comparing LEPA's dominant residues at hypervariable (HV) columns to the corresponding residues in *Brassica* and *Arabidopsis* SRK alleles. If the convergence were driven by pan-Brassicaceae selection on a conserved SI-recognition surface, LEPA's dominant residues would match those of Brassica AND Arabidopsis at most HV columns. If it reflects drift on a shared LEPA ancestral pool, LEPA's dominant residues should be largely LEPA-specific.

The test is decisive ([Figure 15](#figure-15)). For every LEPA HV column, we asked: (i) within LEPA, do the 5 BLs share the same dominant residue among AAAA individuals? and (ii) does that LEPA-consensus residue match the dominant residue at the same alignment column in Brassica and Arabidopsis?

**The answer: 100 % within-LEPA concordance + 74 % LEPA-specific residues.** All 5 BLs share the same dominant residue at 73/73 HV columns (per-BL Shannon entropy near zero, 0.000–0.042 bits) — every BL has fixed the same dominant allele family (Synonymy group 1: Allele_050 + Allele_057 + 11 HV-identical relatives). But that LEPA-consensus residue matches the Brassica dominant residue at only **16/73 columns (22 %)**, the Arabidopsis dominant at 15/73 (21 %), and **both Brassica AND Arabidopsis at only 12/73 (16 %)**. The remaining **54/73 columns (74 %) are LEPA-specific** — the LEPA-consensus residue differs from the dominant residue in both reference genera.

<a name="figure-15"></a>

![Figure 15: Per-BL Shannon entropy decomposition with cross-genera dominant-residue comparison at LEPA HV columns. Top: heatmap of Shannon entropy at each HV column for each of the five BLs (Set1 palette); per-BL mean entropy 0.000–0.042 bits indicates each BL has fixed a single dominant residue at almost every HV column. Bottom: dominant-residue heatmap with all five BLs plus Brassica and Arabidopsis reference rows. The five LEPA BL rows are visually identical (100% within-LEPA concordance) — every BL has fixed the same dominant residue. Brassica and Arabidopsis rows show different colour patterns at most positions: 74% (54/73) of LEPA HV columns are LEPA-specific at the dominant-residue level — confirming drift on a shared LEPA ancestral pool, NOT pan-Brassicaceae selection convergence.](figures/SRK_perBL_entropy_figure.png)

**Interpretation.** If selection across genera were the primary driver, LEPA, Brassica, and Arabidopsis would share dominant residues at most HV columns (functionally fixed across 25+ My of Brassicaceae evolution). We observe the opposite — the dominant residue within LEPA is independent of the dominant residue in either reference genus at 74 % of HV columns. The within-LEPA 100 % concordance reflects the simpler probabilistic outcome that drift fixes the most-common ancestral allele in every bottlenecked lineage — and Synonymy group 1 (Allele_050 + Allele_057 + 11 HV-identical relatives, 103 of 158 BL-assigned AAAA individuals) was already at high enough ancestral frequency that drift fixed it independently in every BL.

The 16 % of HV columns where LEPA = Brassica = Arabidopsis dominant residue mark the **truly deeply conserved SI-recognition positions under genus-spanning selection**. The remaining 84 % are positions where each genus's S-allele pool has been shaped by its own demographic history. For LEPA, that history is **drift-driven fixation of a small ancestral allele family across all five independent bottleneck lineages** — the molecular signature of the C3 cascade Stage 1, observed at residue resolution.

This finding closes the loop with the BL stratification: the pan-BL fixation of Synonymy group 1 documented by the per-BL accumulation curves (Q2), the allele-sharing UpSet ([Figure 6](#figure-6)), and the TP2 reproductive-effort analyses ([Figure 10](#figure-10)) is precisely the same signal recovered here at the residue level. **Drift acting independently in each BL on a shared ancestral allele pool produces the appearance of cross-BL convergence at the residue level (every BL fixes the most-common ancestral allele) without invoking pan-Brassicaceae selection.**

### Conservation implication

The mechanism dictates the intervention strategy. Increasing census population size within an EO alone is insufficient to reverse the AAAA trajectory — if the prevailing AAAA frequency exceeds the reproductive dead-end threshold, within-population growth simply produces more AAAA offspring (Stage 3–4 dynamic). **Allele introduction via inter-BL crosses is the primary lever.** Introducing B, C, and D S-alleles from the 13 AABC individuals (Q4) into crosses with AAAA individuals simultaneously restores SI compatibility AND re-engages balancing selection — the mechanism that, once functional, will favour the spread of introduced S-alleles through subsequent generations (reversing Stage 2). The 13 AABC seed parents are therefore the only endogenous genetic resource capable of reversing the C3 cascade.

---

## Q6 — Cross-based hypothesis testing for informed breeding {#q6-cross-based-hypothesis-testing-for-informed-breeding}

### Why this question matters

The bioinformatics in Q1–Q5 produces *predictions* about which sequence-defined allele bins represent which functional SI specificities, and identifies the diversity reservoir (BL4) and the priority seed parents (13 AABC individuals). To move from prediction to breeding-programme practice, those predictions must be **validated by controlled crosses**. This question lays out the experimental framework, the underlying decision logic, and the resulting cross plan.

### Hypervariable position identification — cross-Brassicaceae validation

The cross design rests on a single critical bioinformatic step: identifying which positions in the S-domain alignment actually discriminate alleles for SI specificity. Distances on these hypervariable (HV) positions are biologically meaningful for recognition specificity, unlike full-domain distances which are dominated by conserved structural residues.

The HV positions are detected by an identical sliding-window Shannon-entropy scan applied separately to LEPA, Brassica (22 alleles), and Arabidopsis (10 alleles) SRK alignments, then validated by cross-Brassicaceae permutation testing ([Figure 16](#figure-16)). The LEPA HV set comprises **73 columns in four regions** of the S-domain (alignment cols 187–199, 240–253, 276–290, 358–388). All three pairwise HV-region overlaps are highly significant (10 000-permutation test): LEPA ↔ Brassica obs = 17 cols (null mean 8.2, p = 0.0005); LEPA ↔ Arabidopsis 22 cols (p = 0.0024); Brassica ↔ Arabidopsis 43 cols (p < 0.0001) — far more co-located than expected by chance.

As **independent structural validation**, all 11 of the 12 mappable SCR9-contact residues from the *B. rapa* eSRK9–SCR9 crystal structure (Ma et al. 2016, PDB 5GYY) fall within these shared HV regions when mapped to LEPA alignment coordinates — direct evidence that the entropy scan is detecting the SI-specificity surface itself.

<a name="figure-16"></a>

![Figure 16: S-domain variability landscape across Brassicaceae. Top panel: smoothed per-column Shannon entropy for LEPA (blue, n = 63 alleles), Brassica (red, n = 22), and Arabidopsis (green, n = 10), each with its own mean + 1 SD threshold (dotted lines). Tracks below mark each species' HV regions detected by identical sliding-window criterion. Bottom track: 11 SCR9-contact residues from the Ma et al. 2016 *B. rapa* eSRK9–SCR9 crystal structure (PDB 5GYY) mapped to LEPA alignment columns. Permutation test: LEPA ↔ Brassica HV overlap = 17 columns (null mean 8.2, p = 0.0005, 10 000 permutations).](figures/SRK_variability_landscape.png)

### Synonymy network — collapsing 63 sequence bins toward functional specificities

Pairwise distances recomputed on the 73 canonical HV positions reveal a strongly bimodal structure: 126 allele pairs are HV-identical (distance = 0); all other within-class pairs are below 0.087; and one allele (Allele_061) is separated from all others by a distance of ~0.96. UPGMA clustering automatically detects this gap, splitting the 63 allele bins into Class I (62 alleles) and Class II (1 allele) — corresponding to the documented Brassicaceae phylogenetic split.

Within Class I, the 126 HV-identical allele pairs form **nine tight synonymy groups** ranging from 2 to 13 alleles per group ([Figure 17](#figure-17)), accounting for 40 of the 62 Class I alleles. The remaining 22 Class I alleles and the single Class II allele are isolated. **If HV identity reliably predicts shared recognition specificity, the 63 sequence bins collapse to 32 functional specificities** (9 synonymy groups + 23 isolated alleles). The largest synonymy group (Synonymy group 1) comprises **13 alleles observed in 103 AAAA individuals** — including Allele_050 and Allele_057, the two pan-BL fixed alleles documented in Q4 — making it the most prevalent specificity in the dataset and the highest priority for synonymy confirmation by crossing.

<a name="figure-17"></a>

![Figure 17: Allele synonymy network — HV-identical synonymy groups. The nine synonymy groups (coloured nodes connected by red HV-identical edges) and 23 isolated alleles (grey nodes) collapse 63 sequence-defined bins into 32 candidate functional specificities, pending experimental confirmation by Synonymy_test crosses. Synonymy group 1 (13 alleles, 103 AAAA individuals) is the highest-priority synonymy target.](figures/SRK_synonymy_network_groups.png)

### The cross plan — five nested hypotheses with explicit genotype constraints

The cross plan is generated by `srk_cross_plan.py` (Step 22e of the bioinformatic pipeline). It combines two independent axes of evidence to assign each candidate cross to a hypothesis level:

- **Axis 1 — Sequence-based cross category** (from the Step 22b cross-design table): Incompatible (HV = 0, same Class) / Synonymy_test (0 < HV < 0.04, same Class) / Compatible_within (HV ≥ 0.04, same Class) / Compatible_cross (different Class).
- **Axis 2 — Genotype-based feasibility** (from the Step 12 zygosity TSV): which parents are AAAA (clean specificity) and which are heterozygous (require polyploid pollen-segregation accounting + paired controls).

Every cross row in the output TSVs is fully traceable from these two axes back to the underlying upstream files; the complete provenance and assumption checklist is in [`Bioinformatics_pipeline.md`](Bioinformatics_pipeline.md) (Step 22e) and [`index.Rmd`](https://svenbuerki.github.io/SRK_bioinformatics/) (`#crossplan`).

The five hypotheses, with their genotype requirements:

| Level | Question | Maternal | Paternal | Allele constraint |
|---|---|---|---|---|
| **H0** | Does SI rejection actually work? | AAAA | AAAA | both alleles in **same** synonymy group |
| **H1a** | Within-Class compatible baseline | AAAA | AAAA | different synonymy groups, no Synonymy_test edge (HV ≥ 0.04) |
| **H1b** | Between-Class compatible baseline | AAAA Class I | Heterozygous (AAAB / AABB) carrier of Allele_061 (Class II) | mother allele ≠ father's other Class I alleles |
| **H2** | Synonymy bin boundaries | AAAA | AAAA | different synonymy groups, with Synonymy_test edge (0 < HV < 0.04) |
| **H3** | Hidden bins (no AAAA representative) | AAAA carrier of allele M known Compatible_within with father's main allele | Heterozygous carrier of the **hidden** allele | M ≠ father's main allele AND M ≠ hidden allele; **requires paired AAAA × AAAA control** |

<a name="figure-18"></a>

![Figure 18: Step 22e hypothesis-testing cross plan summary. Left: cross counts per hypothesis level (101 unique crosses, 1 323 cross attempts at the recommended replicate counts). Right: decision tree per phase outcome, showing how each hypothesis's result updates the validated functional S-allele table. The plan is constrained by the genotype distribution in the current 272-individual dataset — H1b is limited to 1 cross because Allele_061 has only one heterozygous carrier; H3 covers 20 of 29 hidden bins because the remaining 9 lack AAAA mothers with a Compatible_within partner allele.](figures/SRK_cross_plan_summary.png)

### Cross plan results

The phased structure and cross counts are summarised in [Figure 18](#figure-18).

| Hypothesis | Question | Crosses | Replicates each | Total attempts |
|---|---|---:|---:|---:|
| H0 | Does SI rejection work? | 21 | 9 | 189 |
| H1a | Within-Class compatible baseline | 10 | 9 | 90 |
| H1b | Between-Class compatible baseline | **1** | 9 | **9** |
| H2 | Synonymy bin boundaries | 32 | 15 | 480 |
| H3 | Hidden bins (heterozygous donor + paired control) | 37 (= 20 main + 17 paired) | 15 | 555 |
| **Total** | | **101** | — | **1 323** |

**Sample-size constraints.** H1b is sample-limited to 1 cross because Allele_061 has a single carrier in the dataset (Library008_barcode08, AABB Allele_050+Allele_061; AABB pollen segregates 17 % AA + 67 % AB + 17 % BB → 84 % carries the Class II specificity directly). H3 covers 20 of 29 hidden bins; the remaining 9 lack AAAA mothers carrying a Compatible_within partner allele. Adding additional Allele_061 carriers and additional AAAA mothers would directly relax both constraints. The carrier-detection logic reads `Allele_composition` strings authoritatively, so the cross plan automatically expands as new individuals are added — re-run Steps 22a → 22b → 22e after new genotyping data lands.

### Decision tree → validated functional S-allele table

- **H0 pass** (≥ 80 % of Incompatible crosses produce 0 seeds) → SI is functional → proceed. **H0 fail** → SI broken at the Stage 5 polyploid breakdown level → revise the framework.
- **H1a + H1b** outcomes calibrate the seed-yield scale.
- **H2 outcomes** decide which synonymy groups should be **merged** (0 seeds → same SI specificity) or **kept separate** (yield ≥ H1a baseline → distinct specificities). Most extreme case: every Synonymy_test cross is incompatible → 63 bins collapse to **32 functional specificities** (9 synonymy groups + 23 isolated alleles).
- **H3 outcomes** classify each of the 20 testable hidden alleles as compatible / incompatible with at least one major Class I specificity, by comparing total seed yield to the paired-control AA-only baseline.

The combined H2 + H3 outcomes produce a **validated functional S-allele table** that supersedes the sequence-only allele bin definitions and serves as the input to operational seed-orchard design.

---

## Conservation Recommendations

The six questions converge on three immediate priorities for the *Lepidium papilliferum* recovery programme:

1. **Cross all 13 AABC individuals this season** (Q4). They are the only endogenous source of GFS = 0.833 gametes and the only individuals that can simultaneously contribute B, C, and D S-alleles to AAAA recipients — directly reversing the C3 Stage 5 trajectory (Q5). EO67 (5 AABC), EO27 (3), and EO70 (2) are the priority maternal sites; EO76 (0 AABC) is the most degraded and the most in need of allele importation.
2. **Prioritise inter-BL allele transfers, especially BL4 → other BLs** (Q3). BL4 is the species' diversity reservoir (31 of 54 alleles, 14 of which are BL4-private). With no between-EO pollinator connections anywhere in the dataset (Q1), inter-BL crosses are the *only* mechanism for redistributing private alleles. Within-BL or within-EO crossing alone re-segregates the same drift-purged pool and cannot restore the SI system (Q3 TP1 finding).
3. **Schedule the Step 22e cross plan with Synonymy group 1 as the first H2 priority** (Q6). Synonymy group 1 contains 13 HV-identical alleles spanning 103 AAAA individuals (40 % of all AAAA species-wide). An incompatible Synonymy_test cross would consolidate the largest single block of allele bins into one functional specificity, dramatically clarifying the breeding-pool design. Conversely, a compatible Synonymy_test cross would reveal that even within Synonymy group 1 there is genuine functional diversity worth preserving.

The full execution plan is in `SRK_cross_plan_H0_SI_validation.tsv` through `SRK_cross_plan_H3_hidden_bin_tests.tsv` (Step 22e outputs). Outcomes feed back into Step 23 of the bioinformatic pipeline (`SRK_cross_result_analysis_HV.pdf`), which formally tests whether the bioinformatic compatibility predictions match the observed seed yields and produces the **validated functional S-allele table** that anchors operational seed-orchard design.

---

## References

Aguilar, R., Ashworth, L., Galetto, L., & Aizen, M. A. (2006). Plant reproductive susceptibility to habitat fragmentation: Review and synthesis through a meta-analysis. *Ecology Letters*, 9(8), 968–980. https://doi.org/10.1111/j.1461-0248.2006.00927.x

Brennan, A. C., Harris, S. A., Tabah, D. A., & Hiscock, S. J. (2002). The population genetics of sporophytic self-incompatibility in *Senecio squalidus* L. (Asteraceae) I: S allele diversity in a natural population. *Heredity*, 89(6), 430–438. https://doi.org/10.1038/sj.hdy.6800154

Buerki, S., Geisler, M., Martinez, P., Beck, J., & Robertson, I. (2026). Genomics and phylogenetics inform a species recovery plan for a threatened allopolyploid plant. *Botanical Journal of the Linnean Society*, boaf116. https://doi.org/10.1093/botlinnean/boaf116

Busch, J. W., & Schoen, D. J. (2008). The evolution of self-incompatibility when it is costly. *Trends in Plant Science*, 13(3), 128–135. https://doi.org/10.1016/j.tplants.2008.01.001

Byers, D. L., & Meagher, T. R. (1992). Mate availability in small populations of plant species with homomorphic sporophytic self-incompatibility. *Heredity*, 68(4), 353–359. https://doi.org/10.1038/hdy.1992.49

Castric, V., & Vekemans, X. (2004). Plant self-incompatibility in natural populations: A critical assessment of recent theoretical and empirical advances. *Molecular Ecology*, 13(10), 2873–2889. https://doi.org/10.1111/j.1365-294X.2004.02296.x

Comai, L. (2005). The advantages and disadvantages of being polyploid. *Nature Reviews Genetics*, 6(11), 836–846. https://doi.org/10.1038/nrg1711

Elam, D. R., Ridley, C. E., Goodell, K., & Ellstrand, N. C. (2007). Population size and relatedness affect fitness of a self-incompatible invasive plant. *Proceedings of the National Academy of Sciences*, 104(2), 549–552. https://doi.org/10.1073/pnas.0607259104

Ma, R., Han, Z., Hu, Z., Lin, G., Gong, X., Zhang, H., Nasrallah, J. B., & Chai, J. (2016). Structural basis for specific self-incompatibility response in *Brassica*. *Cell Research*, 26(12), 1320–1329. https://doi.org/10.1038/cr.2016.129

Mable, B. K. (2004). Polyploidy and self-compatibility: Is there an association? *New Phytologist*, 162(3), 803–811. https://doi.org/10.1111/j.1469-8137.2004.01055.x

Mable, B. K., Beland, J., & Di Berardo, C. (2004). Inheritance and dominance of self-incompatibility alleles in polyploid *Arabidopsis lyrata*. *Heredity*, 93(5), 476–486. https://doi.org/10.1038/sj.hdy.6800526

Schierup, M. H., Vekemans, X., & Christiansen, F. B. (1997). Evolutionary dynamics of sporophytic self-incompatibility alleles in plants. *Genetics*, 147(2), 835–846. https://doi.org/10.1093/genetics/147.2.835

Vekemans, X., & Slatkin, M. (1994). Gene and allelic genealogies at a gametophytic self-incompatibility locus. *Genetics*, 137(4), 1157–1165. https://doi.org/10.1093/genetics/137.4.1157

Willi, Y., Van Buskirk, J., & Fischer, M. (2005). A threefold genetic Allee effect: Population size affects cross-compatibility, inbreeding depression and drift load in the self-incompatible *Ranunculus reptans*. *Genetics*, 169(4), 2255–2265. https://doi.org/10.1534/genetics.104.034553

Wright, S. (1939). The distribution of self-sterility alleles in populations. *Genetics*, 24(4), 538–552.
