# Polyploid Crossing Model — Slickspot Peppergrass

## Project Overview

This project models allopolyploid genetics for **slickspot peppergrass** (*Lepidium papilliferum*), a tetraploid plant (4 allele copies per locus). The goal is to enumerate all possible crossing outcomes between genotypes, then use optimization to accelerate convergence toward **NFDS equilibrium** (equal S-allele frequencies) across a population.

## Biology Context

- **Ploidy**: Tetraploid (4 allele copies per individual at the S-locus)
- **Allopolyploid**: The genome contains subgenomes from different ancestral species, meaning alleles from distinct origins coexist
- **S-locus alleles**: The population carries a pool of SRK proteins (94 detected, 92 in ingroup). Each individual holds exactly 4 copies (with possible duplicates per the zygosity pattern).
- **Genotype notation**: A sorted tuple/frozenset of 4 integer or string allele IDs (e.g., (1, 5, 23, 47)). Order does not matter — (1,5,23,47) == (47,23,5,1). Always store in sorted canonical form.
- **Equivalence**: Genotypes are unordered multisets of 4 alleles. An individual can carry duplicate alleles (e.g., (1,1,5,23)) but may also carry 4 unique alleles.
- **Crossing**: Each parent contributes 2 alleles (diploid gametes) to offspring. For a tetraploid, gamete formation involves choosing 2 of 4 alleles (C(4,2) = 6 gamete types per parent).

## Self-Incompatibility (SI) & Negative Frequency-Dependent Selection (NFDS)

This project focuses on the **S-locus** controlling self-incompatibility:

- **Self-incompatibility**: A plant rejects pollen that shares S-alleles with the maternal plant. The more S-alleles shared between pollen donor and recipient, the lower the crossing compatibility.
- **Negative frequency-dependent selection (NFDS)**: Rare S-alleles confer a fitness advantage because their carriers are compatible with a larger fraction of the population. Common alleles are at a disadvantage — more potential mates will reject pollen carrying them.
- **Equilibrium target**: Under NFDS, the expected equilibrium is **equal frequency of all S-alleles** in the population. This is driven by frequency-dependent compatibility — rare alleles are naturally favored — not by the assumptions of classical Hardy-Weinberg (random mating, no selection). The optimization accelerates what NFDS already does naturally, rather than overriding it.
- **Practical goal**: Given a real population with known genotypes, determine which crosses to prioritize to most efficiently drive allele frequencies toward equilibrium (equalize rare vs. common alleles).
- **Scale**: With 92 SRK proteins in the ingroup and 4 copies per individual, the number of unique genotypes is C(n+3, 4) (multiset combinations). The real population of 124 individuals spans 24 populations and contains a small subset of possible genotypes.
- **Population grouping**: Individuals are grouped by the **Pop** field (biological populations / element occurrences), not by sequencing library. Analysis supports both within-population and cross-population crossing strategies.

## Core Goals

1. **Enumerate all unique tetraploid genotypes** for a given set of S-alleles (configurable pool size, 50–100+)
2. **Determine crossing compatibility** — for every pair of genotypes, compute SI compatibility and which gamete combinations are allowed
3. **Compute all pairwise crosses** — for compatible pairs, determine the complete set of offspring genotypes and their probabilities
4. **Build a crossing outcome matrix** showing offspring genotype distributions for all parent combinations
5. **Measure distance from equilibrium** — compute current allele frequencies and their deviation from uniform (equal-frequency) target
6. **Optimize crossing strategy** — use gradient descent (or similar optimization) on a fitness landscape to find which individuals to cross to reach allele-frequency equilibrium fastest
7. **Visualize results** — produce plots and tables showing:
   - Current vs. target allele frequency distributions
   - The fitness landscape of allele frequency space
   - Optimal crossing paths toward equilibrium
   - Recommended crosses ranked by effectiveness
   - Compatibility network between genotypes

## Technical Stack

- **Language**: Python 3
- **Key libraries** (expected):
  - `itertools` — combinatorial enumeration of genotypes and gametes
  - `numpy` — matrix operations for crossing tables and frequency calculations
  - `scipy.optimize` — gradient descent / optimization
  - `matplotlib` / `seaborn` — visualization of fitness landscapes and crossing outcomes
  - `pandas` — tabular output of crossing matrices and results
  - `networkx` (optional) — graph representation of crossing paths

## Key Concepts

### Gamete Formation (Tetraploid)
A tetraploid parent with genotype (S1, S5, S23, S47) produces gametes by choosing 2 of 4 alleles:
- C(4,2) = 6 possible gamete types: (S1,S5), (S1,S23), (S1,S47), (S5,S23), (S5,S47), (S23,S47)
- Each gamete equally likely under random segregation (probability 1/6 each)

### Crossing Compatibility (SI)
A cross between parent A (maternal) and parent B (pollen donor) is filtered by SI:
- Pollen gametes sharing any S-allele with the maternal parent may be rejected (depending on the dominance/interaction model)
- Compatible gametes proceed; incompatible ones are blocked
- The offspring distribution is computed only from compatible gamete combinations

### Crossing Two Parents
Offspring genotype = union of one maternal gamete + one compatible paternal gamete. Combine all compatible gamete pairs, sort alleles to canonical form, and tally probabilities.

### Equilibrium Under NFDS
- Under NFDS at the S-locus, equilibrium = **equal allele frequencies** across all S-alleles in the population
- **Distance metric**: Deviation from uniform frequency distribution (e.g., chi-squared, KL divergence, or variance of allele frequencies)
- Real populations may be far from equilibrium due to drift, bottlenecks, or founder effects

### Fitness Landscape & Optimization
- **State**: Current population allele frequency distribution
- **Fitness function**: Negative distance from equal-frequency equilibrium
- **Decision variables**: Which crosses to prioritize (selection weights for specific parent pairs)
- **Compatibility constraint**: Only SI-compatible crosses are valid
- **Optimization**: Gradient descent (or discrete optimization) over crossing strategy space to minimize generations to equilibrium

### Random Mating vs. Optimized Crossing: Model Behavior

The simulation compares two strategies: **random mating** (uniform random parent pairing) and **optimized crossing** (gradient-descent weighted crosses). Under random mating, the model correctly projects that common alleles increase in frequency over time — the opposite of the NFDS equilibrium goal. This is expected behavior, not an artifact, for the following reasons:

1. **SI filtering does not compensate for frequency skew under random mating.** The compatibility check rejects pollen gametes sharing any allele with the maternal plant, but random parent selection still disproportionately picks parents carrying common alleles. Offspring are therefore enriched for those alleles.

2. **Genetic drift in small populations** amplifies this effect. With 15–31 individuals per population, random sampling variance is high. Common alleles are statistically more likely to persist while rare alleles can be lost entirely through drift.

3. **No NFDS mating advantage in the random model.** In nature, NFDS gives rare-allele carriers a reproductive advantage (they are compatible with more mates). The random mating model does not weight by compatibility advantage — parents are chosen uniformly, so rare-allele carriers get no preferential selection. This represents a worst-case scenario: reproduction without natural or managed selection pressure favoring rare alleles.

4. **The optimized strategy counteracts this** by weighting crosses that boost rare alleles and suppress overrepresented ones, explicitly accelerating what NFDS does naturally.

The divergence between random (blue) and optimized (red) trajectories in the simulation plots is the core result of the model — it quantifies the benefit of strategic crossing management over unmanaged reproduction.

## Real Data: SRK Protein Genotyping (Revised)

Five revised data files in `data/revised/` provide SRK (S-locus Receptor Kinase) genotyping results for 128 *L. papilliferum* individuals across 6 sequencing libraries. The revised data includes a **zygosity file** that directly specifies each individual's genotype pattern, eliminating the need for imputation.

### Revised Data Files

| File | Description |
|---|---|
| `SRK_individual_zygosity.tsv` | N_proteins, Zygosity (Homozygous/Heterozygous), Genotype pattern (AAAA/AABB/AABC/ABCD) per individual |
| `SRK_individual_genotypes.tsv` | Binary presence/absence matrix (128 individuals x 94 SRK proteins) |
| `SRK_individual_protein_table.tsv` | Long-format read-level data (1,079 rows); used to resolve AABC doubled allele |
| `SRK_individual_status_report.tsv` | QC classification and filtering status (154 individuals, 128 passed) |
| `sampling_metadata.csv` | Population (Pop), Ingroup/Outgroup flag, Library, sample identifiers (156 samples) |

### Legacy Data Files (in `data/deprecated/`)

The original data files and their generated outputs have been moved to `data/deprecated/` and are **no longer used** in the pipeline.

### Data Quality

- **94 distinct SRK proteins** detected across all individuals
- **92 proteins** detected in ingroup individuals; 39 are singletons (1 individual only), 53 are shared (2+ individuals)
- **128 individuals** passed QC into genotype data; 124 are ingroup *L. papilliferum*
- **4 excluded**: 3 outgroup species (2 *L. montanum*, 1 *L. philonitron*) + 1 with no metadata (SRK_BEA)
- **Population structure**: 24 distinct populations (Pop field); 4 major populations (25, 76, 67, 27) contain 95/124 ingroup individuals

## Genotype Assignment (No Imputation)

The zygosity file directly specifies genotype patterns, enabling **deterministic assignment** for all individuals. No imputation or proportional read-count allocation is needed.

| Pattern | N proteins | Rule | Ambiguity | Count |
|---|---|---|---|---|
| AAAA | 1 | 4 copies of the single protein | None | 64 |
| AABB | 2 | 2 copies of each protein | None | 49 |
| ABCD | 4 | 1 copy of each protein | None | 6 |
| AABC | 3 | Highest read-count protein gets 2 copies; other 2 get 1 each | Resolved by read depth | 5 |

119/124 genotypes are fully deterministic. Only 5 AABC individuals require read-depth resolution from the protein table.

## Project Structure

```
polyploid-model/
  CLAUDE.md                    # This file
  environment.yml              # Mamba/conda environment specification
  data/
    revised/                          # Current data (used by pipeline)
      SRK_individual_zygosity.tsv     #   Zygosity pattern per individual
      SRK_individual_genotypes.tsv    #   Binary presence/absence matrix (128 x 94)
      SRK_individual_protein_table.tsv#   Long-format read-level data
      SRK_individual_status_report.tsv#   QC classification per individual
      sampling_metadata.csv           #   Population, Ingroup flag, sample IDs
    deprecated/                        # Legacy data files (no longer used)
      SRK_individual_genotypes.tsv    #   Old presence/absence matrix (119 x 172)
      SRK_individual_allele_table.tsv #   Old read-level data
      imputed_population.pkl          #   Old imputed genotypes pickle
      imputed_genotypes.tsv           #   Old imputed genotypes TSV
    population.pkl                    # Generated by notebook 00: assigned genotypes + pop groups
    assigned_genotypes.tsv            # Generated by notebook 00: human-readable genotype assignments
  notebooks/
    00_load_data.ipynb         # Data loading, EDA, genotype assignment, export (run first)
    01_genotypes.ipynb         # Genotype enumeration, canonical forms, allele pool setup
    02_crossing.ipynb          # Gamete formation, SI compatibility, pairwise crosses
    03_equilibrium.ipynb       # Allele frequency analysis, distance from equilibrium
    04_optimize.ipynb          # Fitness landscape, gradient descent, crossing strategy
    05_visualize.ipynb         # Comprehensive visualizations and final output
    06_real_analysis.ipynb     # Full pipeline on real population data
```

## Conventions

- **Development environment**: Jupyter Lab notebooks
- Use sorted tuples as canonical genotype representation (e.g., `(1, 5, 23, 47)`)
- Alleles are represented as integers (S1, S2, ..., Sn) for computational efficiency
- All genotype keys should be normalized (sorted) before storage or comparison
- Each notebook should be self-contained but may import shared utility functions defined in earlier notebook cells or a shared module
- Markdown cells should explain the genetics logic at each step
- Verify key results with small worked examples (e.g., 3–4 alleles) before scaling up

## Running

```bash
mamba env create -f environment.yml   # Create the environment
mamba activate polyploid-model        # Activate it
jupyter lab                           # Launch Jupyter Lab

# Conceptual walkthrough (demo data):
# Open notebooks in order: 01 → 02 → 03 → 04 → 05

# Real data analysis:
# 1. Run 00_load_data.ipynb first (generates data/population.pkl)
# 2. Run 06_real_analysis.ipynb (loads pickle, runs full pipeline)
```
