# Polyploid Crossing Model — Overview (v2, 2026-05-24)

A computer model that supports conservation of **slickspot peppergrass** (*Lepidium papilliferum*, hereafter **LEPA**), a tetraploid plant native to southwestern Idaho. The model recommends cross-pollination strategies for managers, predicts how those strategies affect the species' self-incompatibility (SI) system over generations, and stratifies recommendations by the spatial structure of remaining populations.

This document is the plain-language overview of the model after the **2026-05 GFS-aware refactor**. The biology and analytical framework are inherited from the upstream report by [`svenbuerki/SRK_bioinformatics`](https://github.com/svenbuerki/SRK_bioinformatics); this model is the Python implementation that operationalises that framework for cross planning.

---

## What's changed since v1

This is a substantial revision. If you read v1, the headline shifts are:

| Aspect | v1 (pre-2026-05) | v2 (this version) |
|---|---|---|
| **Allele representation** | 94 raw SRK proteins (binary detection matrix) | 27 *biological* alleles after synonymy collapse (19 isolated + 8 synonymy groups derived from upstream Step 22b) |
| **Dataset size** | 124 ingroup individuals, 24 populations | 335 ingroup individuals, 26 populations, 16 Element Occurrences (EOs) across 5 Bottleneck Lineages (BLs) |
| **Genotype tiers** | computed at protein level only | recomputed at biological level after synonymy collapse — many AABC reclassified to AAAB / AABB / AAAA |
| **Conservation diagnostic** | allele-frequency variance only | now reports **TP1** (mating-pool functionality: Shannon evenness J × P_compat) and **TP2** (genotype quality: mean GFS × prop_AAAA) per the LEPA report's two-axis framework |
| **Strategies tracked** | Random, Optimized, Preservation, Demographic | Random, Optimized, Preservation, **GFS-Preservation** (new) — plus optional inter-EO rescue mode |
| **Inheritance mode** | tetrasomic (implicit, hard-coded) | tetrasomic *and* disomic-averaged, with a sensitivity comparison in notebook 09b |
| **SI rule** | strict SI only | strict SI **plus optional empirical leakage** parameter matching Sven's L_hat = 0.18 |
| **Elite-preservation bug** | AAAA carriers of rare alleles were *over*-promoted into elites (the rarity score weighted them up by 4×) | new GFS filter: AAAA → weight 0, AAAB → 0.5, AABB+ → 1.0 |
| **Rare-allele identification** | frequency threshold (0.05) | carrier-count threshold (≤ 2 plants) — the frequency rule broke at 27 alleles where uniform = 0.037 |
| **Allele-loss safety net** | always on (masked the actual dynamics) | opt-in only, default off for honest comparisons |
| **Stratification** | by population (`Pop`) only | by **Pop, EO, and BL** — inter-BL transfers can be modeled and recommended |
| **Open scientific questions** | implicit / unstated | documented in `Questions_for_Sven.md` (5 open questions) |
| **Critical review** | not done | `Critical_Review.md` (5-section balanced-but-steel-man critique) |

The bottom-line empirical consequence: **all 5 BLs and all 6 focus EOs are TP2-CRITICAL at the biological-allele level**, the AABC "rescue lever" shrinks from 15 protein-level individuals to **2** biological-level individuals, and the LEPA report's "inter-BL transfers are the only recovery mechanism" claim is empirically validated in the simulation.

---

## Glossary of terms and acronyms

This document uses domain terminology from plant self-incompatibility genetics and the upstream SRK_bioinformatics framework. Every acronym is defined here on first use.

| Term | Definition |
|---|---|
| **AAAA / AAAB / AABB / AABC / ABCD** | Tetraploid genotype tiers, ordered by descending homozygosity. A tetraploid carries 4 allele copies; the label encodes the copy-count pattern of distinct alleles. AAAA = 1 distinct allele × 4 copies (the "reproductive dead-end" tier). ABCD = 4 distinct alleles × 1 copy each (the rarest, highest-fitness tier). |
| **Allopolyploid** | A polyploid species whose chromosome sets came from different ancestral species. LEPA is described as allopolyploid. Allopolyploids typically (but not always) show disomic inheritance. |
| **BL** | **Bottleneck Lineage.** The LEPA species range partitions into 5 BLs (BL1–BL5), identified by hierarchical clustering of population centroids on a 500 m pollinator-dispersal threshold. Each BL is an evolutionarily independent demographic unit; S-allele erosion has proceeded independently in each. |
| **C3 cascade** | **Compatibility Collapse Cascade** — the 5-stage mechanism (described in the upstream report Q6) by which habitat fragmentation drives a polyploid SI population toward AAAA dominance via leaky self-fertilization, even without any loss-of-function mutation in the SI machinery. |
| **EO** | **Element Occurrence** — a discrete population locality. Each EO belongs to exactly one BL. The 6 *focus EOs* (>= 30 individuals each) are EO18, EO25, EO27, EO67, EO70, EO76. |
| **GFS** | **Genotypic Fitness Score.** Per-individual metric: `GFS = 1 - Σ n_k(n_k - 1) / 12` where `n_k` is the copy count of allele *k*. Equals the proportion of an individual's diploid gametes that carry two distinct alleles. Range 0 (AAAA, 0/6 heterozygous gametes) to 1 (ABCD, 6/6). |
| **inheritance mode** | How the 4 chromosomes pair during gamete formation. **Tetrasomic** = random pairing across all 4 chromosomes (C(4,2) = 6 equally likely gametes). **Disomic** = the 4 chromosomes form two homeologous subgenome pairs that segregate independently (4 gametes per individual, structured by which alleles are on which subgenome). The model supports both, with disomic averaged uniformly over subgenome partitions since per-individual subgenome assignment is unobservable from the current genotyping data. |
| **J** | **Shannon evenness** of allele frequencies: `J = H / ln(k_observed)` where H is Shannon entropy and k is the number of observed alleles. **TP1 x-axis.** J = 1 at NFDS equilibrium; J → 0 as one allele dominates. |
| **L_hat** | Empirical SI leakage rate inferred from observed AAAA proportions in the LEPA dataset. Sven's pair-level `L_hat = 0.18` means ~18% of strict-incompatible plant pairs effectively produce seed via aggregated leakage mechanisms (Lewis 1947 competitive interaction, partial loss-of-function, recognition errors, environmental SI failure). The model uses a per-gamete equivalent `L_HAT_GAMETE = 0.033` (since 1 − (1 − 0.033)^6 ≈ 0.18 for the 6 gametes of tetrasomic inheritance). |
| **NFDS** | **Negative Frequency-Dependent Selection** — the balancing mechanism that protects rare alleles under SI by giving rare-allele carriers more compatible mates. At equilibrium NFDS drives allele frequencies toward equal. |
| **P_compat** | **Compatible-pair fraction.** Fraction of random plant pairs that are SI-compatible (at least one paternal gamete passes the SI check on the maternal stigma). **TP1 y-axis.** A value of 0.40 means ~40% of pairs can produce seed. |
| **prop_AAAA** | Fraction of individuals in the AAAA tier (reproductive dead-ends — AA pollen, rejected wherever A is in the mother). **TP2 y-axis.** |
| **SCR** | **S-locus Cysteine-Rich** protein — the pollen-side recognition ligand in Brassicaceae sporophytic SI. Coats the pollen surface; recognised by SRK on the stigma. |
| **SI** | **Self-Incompatibility** — the molecular mechanism that prevents self-fertilization by rejecting pollen carrying S-alleles shared with the maternal plant. |
| **SRK** | **S-locus Receptor Kinase** — the stigma-side receptor for SCR. The locus that defines S-allele identity in this study. The S-locus is *the* gene this entire project orbits around. |
| **Synonymy group** | A set of SRK proteins with identical hypervariable-region residues (upstream Step 22b). Synonymous proteins are predicted to share SI recognition specificity and so are *collapsed* into a single biological allele in this model. The 8 synonymy groups in LEPA (plus 19 isolated alleles) yield 27 biological alleles from 58 protein bins. |
| **Tetraploid** | An organism with 4 copies of each chromosome (vs 2 in diploids). LEPA is tetraploid. |
| **TP1** | **Tipping Point 1** — mating-pool functionality, the (J, P_compat) diagnostic. Quadrants map to conservation actions: MONITOR, AUGMENT-evenness, AUGMENT-richness URGENTLY, BIOBANK + RESTORE. |
| **TP2** | **Tipping Point 2** — genotype quality. CRITICAL if `mean_GFS < 0.667` AND `prop_AAAA > 0.30`; AT_RISK if one threshold is breached; OK otherwise. |

---

## What the model is for

LEPA is endangered: the species exists in small, fragmented populations across the Snake River Plain of Idaho, USA. Genetic diversity at the SRK (S-locus receptor kinase) gene controls whether plants can find compatible mates — and the surviving populations show severe erosion of that diversity.

The model serves four purposes:

1. **Diagnosis** — quantify how badly each population's SI system is degraded, using the two-axis TP1/TP2 framework from the upstream LEPA report.
2. **Strategy comparison** — simulate four different management strategies (no intervention, optimized crosses, allele preservation, GFS-aware preservation) on each population and compare outcomes.
3. **Cross-plan generation** — produce a ranked list of recommended crosses for conservation managers per EO and per BL, including the operationally distinct *inter-EO rescue crosses* that use the species' 2 biological-level AABC carriers to rescue AAAA-dominated populations elsewhere.
4. **Sensitivity analysis** — quantify how much the recommendations depend on unverified biological assumptions (inheritance mode, SI leakage rate, etc.) so that uncertainty is documented for any conservation decision.

---

## The biology in brief

### Tetraploidy and the S-locus

LEPA is **tetraploid** — each plant carries 4 copies of every chromosome (vs 2 for humans). At the S-locus, each plant therefore holds 4 S-alleles. In a healthy population there are many S-alleles, individuals carry diverse combinations, and self-fertilization is impossible because any pollen carrying a matching S-allele is rejected by the maternal stigma.

### Self-incompatibility (SI)

When pollen lands on a flower, the stigma checks whether the pollen carries any S-allele that matches the maternal plant. If yes → rejected. If no → accepted, fertilization proceeds. This forces plants to mate with genetically different partners and is the species' defence against inbreeding.

In Brassicaceae (LEPA's family), SI is *sporophytic with competitive interaction* — multiple SCR proteins on the pollen surface can interfere with each other, occasionally allowing self-fertilization to bypass the strict rule. The model captures this empirically through a `leakage` parameter calibrated to L_hat ≈ 0.18.

### Negative frequency-dependent selection (NFDS)

Plants with rare alleles can mate with more partners (their pollen is less often rejected) and therefore leave more offspring. This balancing force pushes allele frequencies toward equal — the *NFDS equilibrium*. In small populations, drift overpowers NFDS and rare alleles are lost.

### The C3 cascade

The upstream LEPA report identifies a 5-stage mechanism by which fragmented populations drift toward AAAA dominance: bottleneck → balancing-selection breakdown → mate limitation → drift dominates → leaky-SI selfing produces more AAAA. **Once a population crosses ~30% AAAA, the cascade becomes self-reinforcing and irreversible without intervention.** All 6 focus EOs in LEPA are now past that threshold.

### Genotypic Fitness Score (GFS)

A per-individual metric: the proportion of an individual's diploid gametes that are heterozygous (carry two distinct alleles). Ranges 0 (AAAA, reproductive dead-end) to 1 (ABCD, fully heterozygous gametes). GFS is the right fitness signal for this model because under SI with competitive interaction, heterozygous gametes have more places to land than homozygous ones.

---

## How the model represents the data

### Biological alleles, not raw proteins

Upstream sequencing identifies 49 distinct SRK protein sequences in LEPA. **But many of these proteins share identical hypervariable-region residues** — they are predicted to have the same SI recognition specificity. The upstream `srk_allele_hypotheses.py` analysis groups these into 8 *synonymy groups* (plus 19 isolated singletons), yielding **27 biological alleles** in total, of which 21 are observed in the current dataset.

The model operates on the biological-allele view. Genotype tiers (AAAA, AAAB, AABB, AABC, ABCD) are recomputed at the biological level — and many individuals reclassify (for example, an individual carrying two proteins from Synonymy_group_1 was protein-level AABB but is biological-level AAAA). The protein-level view is preserved in the export for traceability.

### Spatial stratification: BL and EO

Each individual is annotated with its EO (the locality where it was collected) and BL (the regional bottleneck lineage). These annotations let the model:

- Report TP1/TP2 metrics per EO and per BL
- Distinguish *within-EO* crosses (achievable by natural pollinators) from *inter-EO* / *inter-BL* crosses (requiring human-mediated pollen or seed transfer)
- Generate stratified cross plans tailored to each EO

The BL ordering and colour palette match the upstream framework exactly (RColorBrewer Set1: BL1 purple, BL2 blue, BL3 red, BL4 orange, BL5 green) so figures are visually consistent with the upstream report.

### Inheritance mode

The model supports both **tetrasomic** (uniform pairing across all 4 chromosomes) and **disomic-averaged** (two homeologous subgenomes segregating independently). The default is tetrasomic; the disomic-averaged option exists so that the inheritance-mode question can be sensitivity-tested (notebook 09b). Whether LEPA actually shows tetrasomic or disomic inheritance at the S-locus is open question Q1 for Sven Buerki.

### SI leakage

Strict SI (`leakage = 0`) reproduces the textbook rule: any pollen sharing an allele with the mother is rejected. With `leakage > 0`, each rejected paternal gamete has probability `leakage` of being accepted anyway — the empirical aggregate of all SI-escape mechanisms. The model defines `L_HAT_GAMETE = 0.033` (matching Sven's pair-level L_hat = 0.18 for tetrasomic inheritance) as the recommended default for realistic simulations.

---

## How the model recommends crosses

### Step 1: Build the compatibility matrix

For every pair of plants `(i, j)`, the model computes:
- Whether `j` can fertilize `i` under the current SI rule and inheritance mode
- The expected offspring genotype distribution if the cross succeeds
- The expected **GFS** of those offspring (used as the genotype-quality term in the optimizer)
- The expected per-allele frequency contribution to the next generation

### Step 2: Optimize cross weights

Given the compatibility matrix, the optimizer chooses a *weighting* over compatible crosses that minimises a combined objective:

- **Allele-frequency term** (always on): how far the next generation's allele frequencies will be from uniform (the NFDS target)
- **Rare-allele preservation term** (when rare alleles are present): quadratic floor penalty that prevents the optimizer from sacrificing rare alleles
- **Genotype-quality term** (GFS-Preservation strategy only): penalty proportional to `mean(1 - expected_offspring_GFS)` weighted by the chosen crosses

### Step 3: Mandate rare-allele crosses

Before the general optimizer, the model identifies alleles held by only 1-2 plants and mandates one cross per rare allele — pairing the carrier with a partner of tier AABB or higher (per Crossing Plan §3 step 7). This converts an AAAA carrier of a rare allele from "dead end" to "rescue source."

### Step 4: Apply GFS-aware elite filtering

For the Preservation and GFS-Preservation strategies, the top ~10% of individuals are retained unchanged across generations as "elites." The elite score is `(sum of 1/freq for each allele the individual carries) × quality_weight`, where quality_weight is 0 for AAAA, 0.5 for AAAB, and 1.0 for AABB+. This **fixes the v1 bug** where AAAA carriers of rare alleles received the highest possible score (4 × 1/freq) and dominated the elite set — the model was preserving its worst possible carriers.

### Step 5: Inter-EO rescue (optional)

For populations whose within-EO crossing cannot reach a recovered state (the LEPA report's BIOBANK + RESTORE designation), the model can search for inter-EO rescue opportunities — pairing AABC carriers in one EO with AAAA individuals in *other* EOs whose single A allele is not in the carrier's allele set. Each such cross converts an AAAA dead-end lineage into an AABC offspring lineage in one generation.

---

## The four strategies compared

| Strategy | What it does | When to use it |
|---|---|---|
| **Random mating** | No intervention — plants mate at random under whatever SI rule applies. | Baseline / null model. Shows what happens without management. |
| **Optimized** | Mathematically optimal crosses to minimise allele-frequency variance, but no allele preservation. | Best convergence speed; risks losing rare alleles to optimization aggressiveness. |
| **Preservation** | Optimized + elite carryover + mandatory rare crosses + rare-allele floor penalty. | The default for managed conservation programs that want to prevent allele loss. |
| **GFS-Preservation** | Preservation + GFS term in the optimizer + GFS filter on the elite set. | The recommended strategy from this work. Trades a small chi-squared cost for a substantial improvement in mean GFS and prop_AAAA. |

A fifth optional mode — **inter-EO rescue** — augments any of the above by adding the species' 2 biological-level AABC carriers to recipient EOs that lack rescue capacity internally.

---

## Key empirical results

From the integrated analysis in `notebooks/09_genotype_quality.ipynb` (LEPA 2026-05-11 dataset, 335 individuals, 27 biological alleles):

| # | Finding |
|---|---|
| **F1** | Synonymy collapse reduces the AABC "rescue lever" from 15 protein-level individuals to **2 biological-level individuals** — one in EO18 (BL5), one in EO67 (BL4). **Neither is in EO70**, contradicting the protein-level Crossing Plan claim. |
| **F2** | All 5 BLs and all 6 focus EOs are **TP2-CRITICAL** at the biological-allele level. EO70 is worst (mean_GFS = 0.012, prop_AAAA = 0.982); EO27 is best (mean_GFS = 0.201, prop_AAAA = 0.660). |
| **F3** | **EO27 stress test:** GFS-Preservation reaches `mean_GFS = 0.999 ± 0.001` after 5 generations from a starting 0.201. Optimized = 0.919, Preservation = 0.901, Random = 0.765. |
| **F4** | **EO70 lever test:** within-EO crossing fails entirely under random mating (the population's reproductive-collapse signature). With inter-EO rescue (adding the 2 AABC carriers), GFS-Preservation reaches `mean_GFS = 0.997` from starting 0.040. This is the empirical validation of the LEPA report's "inter-BL transfers are the only recovery mechanism" claim. |
| **F5** | **237 feasible inter-EO rescue crosses** with up to 83% AABC offspring per cross at the top of the ranking. EO76 (BL3, 66 mothers reachable) and EO70 (BL2, 55 mothers reachable) are the priority intervention targets. |
| **F6** | The inheritance-mode sensitivity check (notebook 09b) shows the effect of disomic-averaged vs tetrasomic inheritance on mean_GFS is **+0.013 ± 0.020** — within trial noise and much smaller than between-strategy differences. The strategy ranking is robust. |

---

## What's not yet resolved

Five open questions remain that require Sven Buerki's input (documented in `Questions_for_Sven.md`):

1. **Q1 — Inheritance mode.** Is LEPA's S-locus tetrasomic or disomic? Sensitivity analysis shows this is unlikely to change strategy recommendations but does shift per-individual GFS values by ±15%.
2. **Q2 — Per-subgenome haplotype data.** Even if Q1 confirms disomic, the model cannot simulate disomic correctly without per-individual subgenome assignments, which are not in the current data.
3. **Q3 — Leakage mechanism mix.** Is L_hat = 0.18 entirely Lewis 1947 competitive interaction, or aggregates other SI-failure modes? Affects whether the implementation should be mechanistic or empirical.
4. **Q4 — SI-escape candidate status.** Have the 7 phenotyping samples been tested? If confirmed self-compatible, EO76 (which holds 5/7) needs separate management framing.
5. **Q5 — AAAB pollen-donor rule.** The Crossing Plan recommends "AAAB never as pollen donor"; under realistic leakage AAAB carriers of rare alleles might be useful in selective inter-EO transfers. Worth Sven's intuition before locking the rule.

For the full reviewer-style critique of the model (assumptions, mathematical correctness, refactor scope, implementation traps, upstream gaps) see `Critical_Review.md`.

---

## Practical takeaway

For a conservation manager working with LEPA, the model now provides:

1. **A TP1 + TP2 diagnostic** per EO and per BL, directly comparable to the figures in the upstream LEPA report.
2. **A ranked list of recommended within-EO crosses** per EO, generated by the GFS-Preservation strategy.
3. **A short list of high-yield inter-EO rescue crosses** — the species' 2 biological AABC carriers paired with 237 candidate AAAA mothers across other EOs, ranked by expected AABC offspring fraction.
4. **A clear classification of EOs by intervention type**: EO27 and EO25 respond to within-EO crossing; EO18 / EO67 / EO76 benefit from a mix of within-EO and inter-EO; **EO70 requires inter-EO rescue + seed banking and cannot recover by within-EO means alone**.

The strategy comparison shows clearly that GFS-Preservation outperforms Random and the older Preservation strategy on every population tested, especially when realistic SI leakage is included. This is the recommended strategy for any deployed conservation program.

---

## Where things live in this repo

| File / directory | What it contains |
|---|---|
| `CLAUDE.md` | Project conventions, upstream framework, naming, copy-count inference, BL palette |
| `Genotype_Quality_Crossing_Plan.md` | The 2026-04-30 GFS refactor plan (authored against pre-synonymy data; some specifics now stale, the direction stands) |
| `Critical_Review.md` | 5-section critical review of the model (biology, math, refactor plan, implementation traps, upstream gaps) |
| `Questions_for_Sven.md` | Five open questions only Sven can resolve, with our default assumptions documented |
| `MODEL_OVERVIEW.md` | *(this file)* plain-language overview |
| `src/polyploid_utils.py` | Shared utility module — every notebook imports from here |
| `src/bl_constants.py` | BL palette, ordering, EO ↔ BL lookup (mirrored from upstream) |
| `notebooks/00_load_data_Salleles.ipynb` | Data loader: protein matrix + synonymy + zygosity → biological-allele population, with EO/BL annotations |
| `notebooks/01_genotypes.ipynb` … `05_visualize.ipynb` | Synthetic concept demos (toy 8-allele populations) — useful for teaching the mechanics, not used in the deliverable analyses |
| `notebooks/08.2_real_analysis_Salleles.ipynb` | Original real-data analysis (predates the GFS refactor; kept for reference) |
| `notebooks/09_genotype_quality.ipynb` | **The integrated GFS-aware analysis.** Run this for the headline figures, the cross plan, the EO27 stress test, the EO70 lever test, and the synthetic regression test. |
| `notebooks/09b_inheritance_sensitivity.ipynb` | Sensitivity analysis: tetrasomic vs disomic-averaged, strict SI vs empirical leakage. Robustness check for Sven Q1. |
| `notebooks/archive/` | Pre-refactor notebooks (00 loader, 06 demography, 07 collapse prediction, 08 real analysis, 08.1 Pop 27) — kept for traceability, no longer maintained |
| `data/salleles/` | The current S-allele data files synced from upstream `tables/`, including `SRK_synonymy_groups.csv` and the BL summary files |
| `data/population.pkl` | The processed biological-allele population (regenerated by notebook 00) |
| `data/revised/` | Legacy 94-protein-level data (read by archived notebooks) |
| `outputs/figures/` | Generated figures from the real-data notebooks |

If you read just one document beyond this one, read `notebooks/09_genotype_quality.ipynb` — it contains the headline analysis, the figures, and the practical cross recommendations.
