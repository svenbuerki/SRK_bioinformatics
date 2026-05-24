# Polyploid Crossing Model — Slickspot Peppergrass

## Project Overview

This project models allopolyploid genetics for **slickspot peppergrass** (*Lepidium papilliferum*), a tetraploid plant (4 allele copies per locus). The goal is to enumerate all possible crossing outcomes between genotypes, then use optimization to accelerate convergence toward **NFDS equilibrium** (equal S-allele frequencies) across a population.

## Upstream repo & source of truth

This project is downstream of **[`svenbuerki/SRK_bioinformatics`](https://github.com/svenbuerki/SRK_bioinformatics)** (local clone at `~/Repos/SRK_bioinformatics/`). That repo defines the analytical framework (GFS, TP1, TP2, BL stratification, allele synonymy) and produces the data tables this project consumes. **When a metric, threshold, label, or naming convention is in doubt, the upstream repo is canonical.**

Sibling repo: **[`svenbuerki/LEPA_EO_spatial_clustering`](https://github.com/svenbuerki/LEPA_EO_spatial_clustering)** — defines the BL framework consumed by both projects.

`~/Repos/SRK_bioinformatics/modeling/` is an embedded mirror of this repo (sync via the `polyploid` remote + `sync-from-srk` branch). Don't edit there — polyploid-model is canonical.

## Biology Context

- **Ploidy**: Tetraploid (4 copies at the S-locus). Allopolyploid (subgenomes from different ancestors).
- **S-locus alleles** (2026-05-11 dataset): 49 alleles observed, predicted species total ~60 (MM + Chao1 consensus). After synonymy collapse: **27 effective biological alleles** (19 isolated + 8 synonymy groups).
- **Genotype notation**: Sorted tuple of 4 allele IDs, e.g. `("Allele_017", "Allele_017", "Allele_023", "Allele_044")`. Unordered multiset — duplicates allowed.
- **Crossing**: Each parent contributes 2 alleles (diploid gametes). C(4,2) = 6 gamete types per parent.

## Naming conventions (mirror upstream)

- **Allele names**: `Allele_001`..`Allele_058` (3-digit zero-padded). Use these exact strings, not arbitrary integer IDs. A name-to-int mapping is fine internally as long as the string label is preserved for I/O.
- **Sample IDs / Individual**: `LibraryNNN_barcodeMM` (e.g., `Library001_barcode02`). Joinable to `SampleID` in `sampling_metadata.csv`.
- **EO labels**: `EO18`, `EO25`, `EO27`, `EO67`, `EO70`, `EO76` (six focus EOs). Composite EO entries in upstream files use `"EO118; EO76"` form (semicolon-separated).
- **BL labels**: `BL1`–`BL5` (five bottleneck lineages). Colors are **LOCKED** to RColorBrewer Set1: BL1 purple `#984EA3`, BL2 blue `#377EB8`, BL3 red `#E41A1C`, BL4 orange `#FF7F00`, BL5 green `#4DAF4A`. Source: `~/Repos/SRK_bioinformatics/Scripts/srk_bl_constants.py`.
- **BL orderings**: `BL_ORDER` = area-then-connectivity = `BL4, BL5, BL3, BL1, BL2` (for bars/panels/tables). `BL_ORDER_NUMERIC` = `BL1`..`BL5` (for scatter legends).
- **Synonymy groups**: `Isolated` for 19 stand-alone alleles; `"Synonymy group 1"` .. `"Synonymy group 8"` for the 8 grouped sets. Group 1 (15 alleles, includes `Allele_050`/`Allele_051`) is the dominant fixation driver.
- **Genotype tiers**: `AAAA`, `AAAB`, `AABB`, `AABC`, `ABCD` (computed at the **biological** level after synonymy collapse — the `Genotype` column in the zygosity TSV is computed at the **protein** level and must be recomputed).

## Self-Incompatibility (SI) & NFDS

- **SI**: Pollen rejected if it shares any S-allele with the maternal plant.
- **NFDS**: Rare alleles confer fitness advantage (compatible with more mates). Equilibrium = **equal frequency of all S-alleles**.
- **Practical goal**: Optimize crossing strategy to drive allele frequencies toward equilibrium fastest.
- **Scale (2026-05-11)**: 335 ingroup samples, 17 populations, 16 EOs across all 5 BLs. Six focus EOs: EO18, EO25, EO27, EO67, EO70, EO76.

## GFS, TP1, TP2 (canonical, from upstream Scripts/SRK_individual_GFS.R)

**GFS** (per-individual Genotypic Fitness Score): `GFS_i = 1 − Σ_k n_k(n_k−1) / 12` where `n_k` is the copy count of allele *k*.

| Tier | Copy pattern | GFS |
|---|---|---|
| ABCD | (1,1,1,1) | 1.000 |
| AABC | (2,1,1) | 0.833 |
| AABB | (2,2) | 0.667 |
| AAAB | (3,1) | 0.500 |
| AAAA | (4) | 0.000 |

**TP2 (Tipping Point 2)** — EO/BL flagged **CRITICAL** when both `mean_GFS < 0.667` AND `prop_AAAA > 0.30`. Use these exact constants (`TP2_MEAN_GFS = 0.667`, `TP2_PROP_AAAA = 0.30`).

**TP1** — Evenness J (Shannon H / ln k) × P_compat (fraction of random pairs that are SI-compatible). Strict SI (L=0) vs leaky SI (L=0.25; empirical L̂≈0.18). Quadrants → MONITOR / AUGMENT-evenness / AUGMENT-richness / BIOBANK + RESTORE.

## Crossing Strategies

Three strategies compared in simulations:

1. **Random mating** (blue): Uniform random parent pairing. Common alleles increase due to sampling bias — SI filtering alone doesn't compensate. Worst-case scenario without managed selection.
2. **Optimized crossing** (red): L-BFGS-B gradient descent minimizes chi-squared distance from uniform allele frequencies. Recomputed each generation.
3. **Optimized + Preservation** (green): Adds allele loss prevention via three mechanisms:
   - **Elitism** (`select_elites`): Top ~10% individuals by inverse-frequency rarity score retained unchanged
   - **Mandatory rare-allele crosses** (`get_mandatory_rare_crosses`): One best cross per rare allele (freq < 5%) performed first
   - **Optimizer penalty** (`compute_optimal_weights` with `rare_allele_indices`): Quadratic penalty `preservation_weight * max(0, target - expected)^2` prevents sacrificing rare alleles

**Preservation parameters**: `preserve_rare=True`, `rare_threshold=0.05`, `elite_frac=0.1`, `preservation_weight=10.0`

**Tradeoff**: Preservation slightly slows convergence but prevents irreversible allele loss.

4. **Optimized + Preservation + Demography** (green, diamond): Extends strategy 3 with realistic population dynamics:
   - **Logistic growth** (`logistic_n_offspring`): `N_{t+1} = N + N*r*(1 - N/K)` with Poisson noise
   - **Effective population size** (`effective_population_size`): Ne estimated from SI compatibility structure
   - **Harmonic mean Ne** (`ne_harmonic_mean`): Long-term Ne dominated by bottleneck generations

**Demographic parameters**: `carrying_capacity=K`, `growth_rate=0.5`, `demographic_stochastic=True`

## Technical Stack

Python 3: `itertools`, `numpy`, `scipy.optimize`, `matplotlib`/`seaborn`, `pandas`, `networkx`

## Key Concepts

- **Gamete formation**: C(4,2)=6 diploid gametes per tetraploid parent, equally likely
- **SI filtering**: Pollen gametes sharing any allele with maternal plant are rejected
- **Offspring**: maternal gamete (2) + compatible paternal gamete (2) = 4-allele offspring in canonical form
- **Equilibrium**: Equal allele frequencies. Distance measured by variance, chi-squared, KL divergence, extinct/endangered allele counts
- **Fitness landscape**: Optimization over crossing weight space to minimize chi-squared distance from uniform

## Real Data

**Current path (`data/salleles/`)** — S-allele representation synced from `svenbuerki/SRK_bioinformatics` upstream `tables/`:

- `SRK_individual_allele_genotypes.tsv` — 335 individuals × 49 alleles. **Values are observed PROTEIN counts, not tetraploid COPY counts** (upstream README mislabels them). Row sums are NOT 4 in general — e.g., AAAA individuals often have sum=1 or 2. Copy counts must be inferred from the zygosity TSV's `Genotype` + `Allele_composition` columns.
- `SRK_individual_zygosity.tsv` — 335 rows. Columns: `Individual`, `N_distinct_alleles`, `N_total_proteins`, `Zygosity`, `Genotype` (AAAA/AAAB/AABB/AABC/ABCD), `Allele_composition` (e.g., `Allele_018(1)+Allele_050(1)+Allele_049(2)`). The (N) in `Allele_composition` is protein count, not copy count. The `Genotype` value is computed at the **protein level**.
- `sampling_metadata.csv` — 409 rows. Columns include `SampleID`, `Library`, `barcode`, `Pop`, `EO_w_sub`, `Ingroup` (1=analyzed), `OccurrenceID`. Join key: `SampleID = Individual`.
- `SRK_synonymy_groups.csv` — 58 alleles: 19 `Isolated` + 8 `Synonymy group N` (Group 1: 15 alleles incl. dominant `Allele_050`/`Allele_051`; Group 2: 9; Group 3: 4; Group 4: 3; Groups 5–8: 2 each). Yields **27 effective biological alleles**.

**Copy-count inference** (from protein counts + zygosity Genotype + Allele_composition):
- AAAA: 4 copies of the single named allele
- AAAB: smaller protein count → B (1 copy); larger → A (3 copies)
- AABB: 2 + 2 copies
- AABC: protein count 2 → doubled allele; the other two are singletons
- ABCD: 1 + 1 + 1 + 1

**Synonymy collapse** — at the biological level, alleles in the same synonymy group are a single SI specificity. Example: `Allele_050(1)+Allele_049(2)` (both Group 1) collapse to a single biological allele with 3 copies, shifting a protein-level AABC to a biological-level AAAB. **The crossing model operates on biological alleles**, so genotype tiers and GFS must be recomputed after synonymy collapse — do not use the upstream `Genotype` column verbatim.

**Join chain** (Individual → BL):
```
SRK_individual_allele_genotypes.tsv  (Individual)
    └─ join SampleID = Individual ─→  sampling_metadata.csv  (Pop, EO_w_sub)
                                          └─ join on EO ─→  EO_group_BL_summary.csv  (BL, Drift_index)
                                                              [mirrored from sibling spatial-clustering repo]
```

**Legacy path (`data/revised/`)** — 94-protein binary matrix used by archived notebooks (`notebooks/archive/`). Files: `SRK_individual_genotypes.tsv` (128×94 binary), `SRK_individual_zygosity.tsv` (older schema), `SRK_individual_protein_table.tsv`, `SRK_individual_status_report.tsv`, `sampling_metadata.csv`.

## Project Structure

```
polyploid-model/
  CLAUDE.md                    # This file
  Genotype_Quality_Crossing_Plan.md  # GFS / TP2 refactor plan (2026-04-30)
  environment.yml              # Mamba/conda environment
  src/
    polyploid_utils.py         # Shared utility functions (all notebooks import from here)
  data/
    salleles/                  # Current S-allele data (dosage matrix + synonymy groups)
    revised/                   # Legacy 94-protein binary matrix (read by archived notebooks)
    filespopulation27/         # Pop 27 subset (read by archived 08.1)
    deprecated/                # Pre-revised legacy data
    population.pkl             # Generated by 00_load_data_Salleles.ipynb
    assigned_genotypes.tsv     # Generated by 00_load_data_Salleles.ipynb
  notebooks/
    00_load_data_Salleles.ipynb     # Load S-allele data, assign genotypes, export pickle (run first)
    01_genotypes.ipynb              # Genotype enumeration, canonical forms (synthetic concepts)
    02_crossing.ipynb               # Gamete formation, SI compatibility (synthetic concepts)
    03_equilibrium.ipynb            # Allele frequency analysis (synthetic concepts)
    04_optimize.ipynb               # Fitness landscape, gradient descent (synthetic concepts)
    05_visualize.ipynb              # Comprehensive visualizations (synthetic concepts)
    08.2_real_analysis_Salleles.ipynb  # Full pipeline on real S-allele data
    09_genotype_quality.ipynb       # (planned) GFS-aware crossing strategy on biological alleles
    archive/                        # Pre-refactor notebooks on the old data path
      README.md                     # Explains what was archived and why
      00_load_data.ipynb            # Old loader for data/revised/ (94-protein binary matrix)
      06_demography.ipynb           # Demography on old population.pkl
      07_collapse_prediction.ipynb  # Collapse cascade on old population.pkl
      08_real_analysis.ipynb        # Old real-data pipeline
      08.1_real_analysis_pop27.ipynb  # Pop 27 subset on old data
```

## Conventions

- Sorted tuples as canonical genotype representation: `(1, 5, 23, 47)`
- Integer allele IDs for computational efficiency
- All notebooks import from `src/polyploid_utils.py` via `sys.path.insert(0, "../src")`
- Jupyter Lab for development; markdown cells explain genetics logic
- **Always run Python commands within the mamba environment**: `mamba run -n polyploid-model python3 ...`
- Use `layout="constrained"` in `plt.subplots()` instead of `plt.tight_layout()`

## Running

```bash
mamba env create -f environment.yml   # Create the environment
mamba activate polyploid-model        # Activate it
jupyter lab                           # Launch Jupyter Lab

# Conceptual walkthrough (synthetic data): notebooks 01 → 02 → 03 → 04 → 05
# Real data (S-allele dosage): run 00_load_data_Salleles, then 08.2_real_analysis_Salleles
# Archived (old protein-level path): see notebooks/archive/README.md
```
