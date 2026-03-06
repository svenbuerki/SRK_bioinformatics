# Polyploid Crossing Model — Slickspot Peppergrass

## Project Overview

This project models allopolyploid genetics for **slickspot peppergrass** (*Lepidium papilliferum*), a tetraploid plant (4 allele copies per locus). The goal is to enumerate all possible crossing outcomes between genotypes, then use optimization to accelerate convergence toward **NFDS equilibrium** (equal S-allele frequencies) across a population.

## Biology Context

- **Ploidy**: Tetraploid (4 allele copies per individual at the S-locus)
- **Allopolyploid**: The genome contains subgenomes from different ancestral species, meaning alleles from distinct origins coexist
- **S-locus alleles**: The population carries a large pool of possible S-alleles (50–100+ distinct alleles). Each individual holds exactly 4 of these.
- **Genotype notation**: A sorted tuple/frozenset of 4 integer or string allele IDs (e.g., (1, 5, 23, 47)). Order does not matter — (1,5,23,47) == (47,23,5,1). Always store in sorted canonical form.
- **Equivalence**: Genotypes are unordered multisets of 4 alleles. An individual can carry duplicate alleles (e.g., (1,1,5,23)) but may also carry 4 unique alleles.
- **Crossing**: Each parent contributes 2 alleles (diploid gametes) to offspring. For a tetraploid, gamete formation involves choosing 2 of 4 alleles (C(4,2) = 6 gamete types per parent).

## Self-Incompatibility (SI) & Negative Frequency-Dependent Selection (NFDS)

This project focuses on the **S-locus** controlling self-incompatibility:

- **Self-incompatibility**: A plant rejects pollen that shares S-alleles with the maternal plant. The more S-alleles shared between pollen donor and recipient, the lower the crossing compatibility.
- **Negative frequency-dependent selection (NFDS)**: Rare S-alleles confer a fitness advantage because their carriers are compatible with a larger fraction of the population. Common alleles are at a disadvantage — more potential mates will reject pollen carrying them.
- **Equilibrium target**: Under NFDS, the expected equilibrium is **equal frequency of all S-alleles** in the population. This is driven by frequency-dependent compatibility — rare alleles are naturally favored — not by the assumptions of classical Hardy-Weinberg (random mating, no selection). The optimization accelerates what NFDS already does naturally, rather than overriding it.
- **Practical goal**: Given a real population with known genotypes, determine which crosses to prioritize to most efficiently drive allele frequencies toward equilibrium (equalize rare vs. common alleles).
- **Scale**: With 50–100 S-alleles and 4 copies per individual, the number of unique genotypes is C(n+3, 4) (multiset combinations). For n=50, this is ~316,251 possible genotypes — but real populations will only contain a subset.

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

## Real Data: SRK Allele Genotyping

Two real data files in `data/` contain SRK (S-locus Receptor Kinase) genotyping results for 119 *L. papilliferum* individuals across 5 library/population groups:

- **`SRK_individual_genotypes.tsv`** — Binary presence/absence matrix (119 individuals x 172 alleles). Each cell is 0 or 1 indicating whether that allele was detected in that individual.
- **`SRK_individual_allele_table.tsv`** — Long-format read-level data (648 rows). Each row is one sequencing read mapping an allele to an individual. Multiple rows per (individual, allele) pair reflect read depth.

### Data Quality Limitations

- **Incomplete genotyping**: Only 6/119 individuals have exactly 4 detected alleles. Most have 1-3 due to low sequencing depth.
- **Singletons**: ~126/172 alleles appear in only 1 individual, which may reflect sequencing artifacts or genuinely rare alleles.
- **Low read depth**: ~35% of allele calls are supported by a single read, making them unreliable.
- **Population structure**: Individuals come from 5 library groups (Library001-005) plus 1 hybrid, representing distinct populations.

## Imputation Strategy

Since *L. papilliferum* is tetraploid (4 allele copies per individual), all individuals must have exactly 4 alleles assigned. The imputation rules below convert variable-depth sequencing data into tetraploid genotypes.

**WARNING: All real-data results depend on these imputed genotypes. Interpret with caution.**

| Detected alleles | Rule | Example |
|---|---|---|
| 1 | Homozygous: (A,A,A,A) | 1 allele, any reads -> 4 copies |
| 2 | Proportional to read counts, fill 4 slots | A(6) B(2) -> (A,A,A,B) |
| 3 | Proportional to read counts, fill 4 slots | A(5) B(2) C(1) -> (A,A,B,C) |
| 4 | Direct: one copy each | A,B,C,D -> (A,B,C,D) regardless of reads |
| >4 | **Strategy A:** top 4 by reads, proportional fill | Drops lowest-read alleles |
| >4 | **Strategy B:** top 3 fixed, rotate 4th through remaining | Creates multiple genotype entries |

### Largest Remainder Method (Slot Allocation)

1. Compute each allele's proportion of total reads
2. Multiply by 4 to get ideal (fractional) slot counts
3. Floor each value to get guaranteed slots
4. Distribute remaining slots to alleles with largest fractional remainders
5. Ties broken alphabetically (by allele name)

For exactly 4 detected alleles, skip proportional allocation and assign 1 copy each.

## Project Structure

```
polyploid-model/
  CLAUDE.md                    # This file
  environment.yml              # Mamba/conda environment specification
  data/
    SRK_individual_genotypes.tsv      # Binary presence/absence matrix (raw)
    SRK_individual_allele_table.tsv   # Long-format read-level data (raw)
    imputed_population.pkl            # Generated by notebook 00: imputed genotypes + metadata
    imputed_genotypes.tsv             # Generated by notebook 00: imputed genotypes (human-readable)
  notebooks/
    00_load_data.ipynb         # Data loading, EDA, imputation, export (run first for real data)
    01_genotypes.ipynb         # Genotype enumeration, canonical forms, allele pool setup
    02_crossing.ipynb          # Gamete formation, SI compatibility, pairwise crosses
    03_equilibrium.ipynb       # Allele frequency analysis, distance from equilibrium
    04_optimize.ipynb          # Fitness landscape, gradient descent, crossing strategy
    05_visualize.ipynb         # Comprehensive visualizations and final output
    06_real_analysis.ipynb     # Full pipeline on real imputed population data
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
# 1. Run 00_load_data.ipynb first (generates data/imputed_population.pkl)
# 2. Run 06_real_analysis.ipynb (loads pickle, runs full pipeline)
```
