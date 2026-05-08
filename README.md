# SRK Bioinformatics Pipeline

## Overview

This repository contains a comprehensive bioinformatics pipeline for analyzing S-receptor kinase (SRK) haplotype diversity in self-incompatible plant species using Oxford Nanopore long-read sequencing. The pipeline is specifically designed for threatened plant species research, enabling assessment of self-incompatibility system integrity in populations where genetic diversity loss may compromise reproductive success.

## Authors

**Sven Buerki, Ph.D.** (he/him)\
Associate Professor\
Dr. Christopher Davidson Endowed Chair in Botany\
Department of Biological Sciences\
Boise State University\
Boise, Idaho, USA

📧 [svenbuerki\@boisestate.edu](mailto:svenbuerki@boisestate.edu)

**Jim Beck, Ph.D.**\
Manager/HPC Engineer IV\
Research Computing\
Boise State University\
Boise, Idaho, USA

📧 [jimbeck\@boisestate.edu](mailto:jimbeck@boisestate.edu)

## Key Features

-   **Multi-replicate CANU assembly** optimized for rare haplotype recovery
-   **Polyploid-aware variant calling** and read-backed phasing with WhatsHap
-   **Reference-guided sequence orientation** and gap-filling procedures
-   **SRK-specific protein validation** with optional domain analysis
-   **Abundance-based filtering** to distinguish authentic alleles from technical artifacts
-   **Tetraploid zygosity analysis** with allele copy-count dosage inference (AAAA, AAAB, AABB, AABC, ABCD) for complex polyploid genetics
-   **Individual Genotypic Fitness Score (GFS)** quantifying heterozygous gamete proportion per individual — differentiates AABB (GFS = 0.667) from AAAB (GFS = 0.500) to guide seed parent selection and EO prioritisation
-   **Population genetic dataset generation** for demographic inference
-   **Species allele richness estimation** using Michaelis-Menten, Chao1, and iNEXT estimators to quantify total SRK diversity and population deficits relative to the species optimum

## Pipeline Workflow

The pipeline consists of **23 main steps organized into four phases**, progressing from within-library sequence assembly to cross-library integration, population genetic analyses, and experimental validation of S-allele hypotheses through controlled crosses.

## Phase 1: SRK Amplicon Sequence Assembly

1.  **Prepare Canonical Sequences** – Preparation of reference SRK sequences used for downstream assembly and validation.
2.  **Nanopore Amplicon Assembly and Phasing Pipeline** – Multi-CANU assembly of Nanopore amplicons followed by haplotype phasing to reconstruct allelic sequences.
3.  **Haplotype Orientation and Consolidation Pipeline** – Standardization of sequence orientation and consolidation of phased haplotypes.
4.  **FASTA Sequence Filtering and Reference Integration Pipeline** – Quality filtering of sequences and integration of validated reference alleles.
5.  **Multiple Sequence Alignment** – Alignment of nucleotide sequences using MAFFT.
6.  **Exon Extraction from Multiple Sequence Alignments Pipeline** – Extraction of coding regions based on AUGUSTUS gene annotations.
7.  **MSA Gap Backfilling and Terminus Processing Pipeline** – Reference-guided backfilling of terminal alignment gaps and sequence terminus correction.
8.  **DNA to Amino Acid Translation Pipeline** – Frame-specific translation of nucleotide sequences into proteins. **Note:** This step is optional and may be skipped if protein sequences are not required.

## Phase 2: Functional Proteins, S Alleles and Genotyping

9.  **SRK Protein Translation, Alignment, and Abundance Filtering Pipeline** – Translation of SRK alleles into proteins, alignment, and filtering based on abundance thresholds to retain functional candidates.
10. **Distance-Based SRK S-Allele Definition** – Grouping of functional proteins into allele bins using pairwise amino acid p-distance on the S-domain ectodomain, with sensitivity analysis and amino acid variation heatmaps.
11. **SRK S-Allele Genotyping Pipeline** – Assignment of distance-defined S-allele bins to individuals, producing long-format allele tables (with allele copy counts) and count-based genotype matrices. Each cell in the wide matrix records the number of distinct proteins from that individual assigned to that allele bin, capturing dosage information for downstream zygosity inference and allele frequency analysis.
12. **SRK Zygosity Analysis Pipeline for Tetraploid Species** – Classification of individuals as homozygous or heterozygous using the allele count matrix. Infers tetraploid genotype classes (AAAA, AAAB, AABB, AABC, ABCD) from both the number of distinct allele bins and their copy counts, and outputs an `Allele_composition` string (e.g. `Allele_044(2)+Allele_047(2)`) for direct use in cross planning.

**⚠️ Step 12b (In Progress) — SRK Class Assignment and Dominance Prediction** – Phylogenetic classification of SRK proteins into Class I (dominant) or Class II (recessive) using maximum-likelihood inference (IQ-TREE3). Automated class assignment by longest internal branch exploits the ancient divergence between SRK classes maintained by balancing selection. Per-individual dominance predictions guide cross-compatibility decisions in a conservation breeding program. **This analysis is under active development and may not function correctly on all datasets — results require validation against published SRK reference sequences with known class assignments.**

## Phase 3: Data Analyses

13. **Element Occurrence and Bottleneck Lineage Integration** – Bridges the SRK individual-level dataset to the spatial framework defined in the sibling project [LEPA_EO_spatial_clustering](https://github.com/svenbuerki/LEPA_EO_spatial_clustering). For each post-QC individual, joins EO label, geographic group, bottleneck lineage (BL1–BL5), and drift index from `EO_group_BL_summary.csv`. Produces `SRK_individual_BL_assignments.tsv`, the join key used by all subsequent Phase 3 steps to add BL stratification (EO sorted within BL for management; BL aggregated for inferential power on small localities).
14. **Population Genetics Statistics** – Estimation of population-level diversity metrics, including heterozygosity, mean alleles per individual, total allele counts, and effective allele numbers. Allele frequencies are based on copy counts summed across individuals (from the count matrix), giving a proper tetraploid frequency estimate.
15. **Allele Accumulation Curves** – Rarefaction-based analysis of SRK allele discovery across individuals to evaluate patterns consistent with negative frequency-dependent selection versus genetic drift. Includes estimation of total species allele richness using Michaelis-Menten asymptote fitting, Chao1, and iNEXT estimators. Outputs an empirical species optimum used as a baseline in steps 16 and 17.
16. **Allele Frequency Analysis** – Species- and population-level χ² tests of allele frequency distributions to assess deviations from equal-frequency expectations under NFDS. The estimated species allele richness from step 15 is used as the optimum, quantifying how many alleles each population is missing relative to the species pool.
17. **TP1 Tipping Point Analysis** – Diagnostic scatter plot synthesising Steps 14–16, positioning each EO by the proportion of the species optimum retained (x axis) and by allele frequency evenness (Ne/N, y axis). EOs breaching both thresholds (< 50% of optimum; Ne/N < 0.80) are flagged CRITICAL.
18. **Allele Composition Comparison Across Element Occurrences** – UpSet plot and pairwise sharing heatmap quantifying how S-allele sets partition across Element Occurrences, identifying private alleles and alleles shared across all populations. Requires only standard Python dependencies (pandas, numpy, matplotlib).
19. **Individual Genotypic Fitness Score (GFS)** – Per-individual metric quantifying the proportion of heterozygous diploid gametes a tetraploid can produce. Differentiates dosage-imbalanced genotypes (AABB vs AAAB) invisible to zygosity analysis alone. Outputs per-individual GFS values and ranked seed parent lists per EO.
20. **TP2 Tipping Point Analysis** – EO-level assessment placing mean GFS and proportion of AAAA individuals in interaction. EOs breaching both thresholds simultaneously (mean GFS < 0.667; proportion AAAA > 30%) are flagged CRITICAL, AT RISK, or OK.
21. **Reproductive Effort Support per Element Occurrence** – Horizontal proportional bar chart showing the fraction of individuals at each GFS tier per EO.

## Phase 4: Testing S-allele Hypotheses

22. **HV-Based Allele Hypothesis Testing and Crossing Design** – Moving-window variability scan of the S-domain alignment identifies hypervariable (HV) positions where alleles actually diverge. Pairwise distances restricted to those HV positions are used for UPGMA clustering, which automatically detects the Class I / Class II phylogenetic split. Within the majority class, three cross categories are assigned by HV distance: W (HV-identical alleles, d = 0 — expected incompatible), N (small HV difference, d < threshold — synonymy test), and P_within (substantial within-class HV divergence — expected compatible). Between-class crosses are labelled P_cross (guaranteed compatible positive controls). An allele similarity heatmap and cross design summary figure are also produced.
23. **Cross Result Analysis** – Reads completed crossing records and tests whether the W / N / P category predicts seed yield, validating the sequence-based allele definitions against experimental cross-compatibility data. Activated by setting `CROSS_TSV` in the script once crossing data are available.

## Requirements

### Software Dependencies

-   **CANU** (v2.0+) - Long-read assembly
-   **RACON** (v1.4+) - Assembly polishing
-   **minimap2** (v2.17+) - Read alignment
-   **samtools** (v1.10+) - SAM/BAM processing
-   **FreeBayes** (v1.3+) - Variant calling
-   **WhatsHap** (v1.0+) - Haplotype phasing
-   **bcftools** (v1.10+) - VCF processing
-   **MAFFT** (v7.0+) - Multiple sequence alignment
-   **seqkit** (v0.12+) - FASTA manipulation

### Python Dependencies

``` bash
pip install biopython pandas numpy
```

### R Dependencies

``` r
install.packages(c("dplyr", "ggplot2", "readr", "tidyr", "forcats", "scales", "bookdown"))

# Optional — required for iNEXT richness estimation in step 15:
install.packages("iNEXT")

# Optional — required for labelled scatter plot in step 18:
install.packages("ggrepel")
```

## Installation

1.  Clone the repository:

``` bash
git clone https://github.com/svenbuerki/SRK_bioinformatics.git
cd SRK_bioinformatics
```

2.  Install required software (see Requirements section)

3.  Make scripts executable:

``` bash
chmod +x scripts/*.sh
chmod +x scripts/*.py
```

## Usage

### Quick Start

1.  **Prepare your data structure:**

```         
project/
├── Library001/
│   ├── barcode01/
│   │   └── *.fastq.gz files
│   ├── barcode02/
│   └── ...
└── reference_sequences.fasta
```

2.  **Phase 1 — SRK Amplicon Sequence Assembly (Steps 1–8):**

``` bash
# Step 2: Assembly and phasing
./scripts/01_nanopore_assembly_phasing.sh

# Step 3: Orientation
./scripts/02_haplotype_orientation.sh

# Continue with remaining steps...
```

3.  **Phase 2 — Functional Proteins, S Alleles and Genotyping (Steps 9–12):**

``` bash
# Steps 9-12: Protein filtering, allele definition, genotyping, zygosity
python scripts/09_srk_protein_filtering.py
Rscript scripts/10_srk_allele_definition.R
Rscript scripts/11_srk_allele_genotyping.R
Rscript scripts/12_zygosity_analysis.R
```

4.  **Phase 3 — Data Analyses (Steps 13–21):**

``` bash
# Step 13: EO + bottleneck lineage integration (run first — produces the join key for all downstream steps)
python3 SRK_BL_integration.py

# Step 14: Population genetics statistics
Rscript SRK_population_genetic_summary.R

# Step 15: Allele accumulation curves and richness estimation
# (must run before steps 16 and 17)
Rscript SRK_allele_accumulation_analysis.R

# Step 16: Allele frequency analysis (reads SRK_species_richness_estimates.tsv)
Rscript SRK_chisq_species_population.R

# Step 17: TP1 tipping point analysis (reads Steps 14 and 15 outputs)
Rscript SRK_TP1_tipping_point.R

# Step 18: Allele composition comparison across Element Occurrences
python SRK_allele_sharing_EOs.py

# Steps 19–20: Individual Genotypic Fitness Score (Step 19) and TP2 tipping point analysis (Step 20)
Rscript SRK_individual_GFS.R
```

5.  **Phase 4 — Testing S-allele Hypotheses (Steps 22–23):**

``` bash
# Step 22: HV-based allele hypothesis testing and crossing design
python srk_allele_hypotheses.py

# Step 23: Cross result analysis
# Set CROSS_TSV = "<your_cross_results_file>" in the script, then re-run:
python srk_allele_hypotheses.py
```

### Detailed Usage

See the [Pipeline Documentation](https://svenbuerki.github.io/SRK_bioinformatics/) for comprehensive step-by-step instructions.

For a concise step-by-step protocol (scripts, inputs, outputs, key parameters), see [Bioinformatics_pipeline.md](Bioinformatics_pipeline.md).

## Input Data Requirements

-   **Nanopore amplicon sequencing data** (FASTQ format)
-   **Canonical reference sequences** (FASTA format)
-   **Gene structure annotations** (AUGUSTUS CSV format)
-   **Organized directory structure** by library and barcode

## Output Files

### Key Outputs

-   `SRK_functional_proteins.fasta` - Validated SRK protein sequences
-   `SRK_individual_allele_genotypes.tsv` - Allele copy-count matrix (individuals × allele bins; values = number of distinct proteins per individual per allele bin)
-   `SRK_individual_allele_table.tsv` - Long-format allele assignment table (Individual, Protein, Allele, Count)
-   `SRK_individual_zygosity.tsv` - Zygosity classifications with columns: `N_distinct_alleles`, `N_total_proteins`, `Zygosity`, `Genotype` (AAAA/AAAB/AABB/AABC/ABCD), `Allele_composition`
-   `SRK_self_compatible_candidates.txt` - Potentially self-compatible individuals

-   `SRK_TP1_summary.tsv` - Per-EO allele richness, frequency evenness, and TP1 status
-   `SRK_TP1_tipping_point.pdf` - TP1 diagnostic scatter plot (richness retained × frequency evenness)

### Population Genetic Outputs (Phase 3)

#### Step 13 — EO + Bottleneck Lineage integration

-   `SRK_individual_BL_assignments.tsv` - Per-individual EO, Group, BL (BL1–BL5), Drift_index, BL_status (Assigned / Inferred / Unassigned). Join key consumed by every Phase 3/4 R script that adds a BL stratification.

#### Step 14 — Population genetics

-   `SRK_population_genetic_summary.tsv` - EO-level summary with `BL` and `Drift_index` columns; rows sorted by `BL → EO`
-   `SRK_population_genetic_summary_BL.tsv` - **NEW**. BL-level aggregate (5 rows) with `N_EOs`, `Mean_drift`, plus per-group metrics
-   `SRK_population_genetic_summary.pdf` - 2-page PDF: page 1 EO bars sorted by BL with Dark2 BL color palette; page 2 BL aggregate bars

#### Step 15 — Allele accumulation

-   `SRK_allele_accumulation_curves.pdf` - Species + per-EO + per-BL pages with MM/Chao1/iNEXT asymptote reference lines
-   `figures/SRK_allele_accumulation_species.png` - Standalone species curve (54 alleles observed in 272 ingroup individuals; MM = 69, Chao1 = 73)
-   `figures/SRK_allele_accumulation_combined.png` - All qualifying EO curves on shared axes, **colored by parent BL**
-   `figures/SRK_allele_accumulation_BL_combined.png` - **NEW**. The 5 BL aggregate curves on shared axes
-   `figures/SRK_allele_accumulation_drift_erosion.png` - EO stacked bars sorted by BL with BL color strip on the x-axis baseline
-   `figures/SRK_allele_accumulation_BL_drift_erosion.png` - **NEW**. BL aggregate stacked bars (BL4 = 32% lost; BL1/BL2 ≥ 86% lost)
-   `SRK_allele_accumulation_stats.tsv` - Per-level curve statistics (Species + 6 EOs + 5 BLs), MM/Chao1/iNEXT estimates, sampling adequacy targets
-   `SRK_species_richness_estimates.tsv` - Consensus species allele richness estimate (input to Step 16). Computed from all 272 ingroup individuals.

#### Step 16 — Allele frequency χ²

-   `SRK_chisq_species_population.tsv` - χ² test statistics at three levels (Species, EO, BL); new `BL` column; rows ordered Species → EOs (sorted by BL) → BLs
-   `SRK_chisq_species_population_frequency_plots.pdf` - 12 pages (1 species + 6 EOs sorted by BL + 5 BLs); EO/BL pages colored by parent BL palette

#### Steps 17–21 (other Phase 3 outputs)

-   `SRK_TP1_summary.tsv` - Per-EO allele richness, frequency evenness, and TP1 status (Step 17)
-   `SRK_TP1_tipping_point.pdf` - TP1 diagnostic scatter plot (richness retained × frequency evenness) (Step 17)
-   `SRK_allele_upset_EOs.pdf` - UpSet plot of all pairwise and higher-order allele set intersections across Element Occurrences (Step 18)
-   `SRK_allele_sharing_heatmap_EOs.pdf` - Pairwise allele sharing heatmap between Element Occurrences (Step 18)
-   `SRK_individual_GFS.tsv` - Per-individual Genotypic Fitness Score, genotype class (AAAA/AAAB/AABB/AABC/ABCD), and EO assignment (Step 19)
-   `SRK_EO_GFS_summary.tsv` - EO-level mean GFS, genotype class proportions, and Tipping Point 2 status (CRITICAL / AT RISK / OK) (Step 20)
-   `SRK_GFS_plots.pdf` - Four diagnostic plots: stacked composition bars, individual GFS jitter with mean, TP2 tipping point map, and absolute count bars (Steps 19–20)

### Phase 4 Outputs (Testing S-allele Hypotheses)

-   `SRK_variability_landscape.pdf` / `figures/SRK_variability_landscape.png` - Per-column S-domain variability profile with HV regions shaded (Step 22)
-   `SRK_HV_allele_distances.tsv` - 63×63 pairwise distance matrix computed on 75 HV positions (Step 22)
-   `SRK_functional_allele_groups.tsv` - Allele bin → phylogenetic class assignment, AAAA count, cross power (Step 22)
-   `SRK_synonymy_candidates.tsv` - All within-class allele pairs with HV distance and testability flag (Step 22)
-   `SRK_allele_similarity_heatmap.pdf` / `figures/SRK_allele_similarity_heatmap.png` - 63×63 HV similarity heatmap ordered by UPGMA with class strips (Step 22)
-   `SRK_AAAA_cross_design_HV.tsv` - All AAAA × AAAA pairs ranked by category (W / N / P_within / P_cross) with HV distance and expected outcome (Step 22)
-   `SRK_cross_design_summary.pdf` / `figures/SRK_cross_design_summary.png` - Three-panel figure: HV distance distribution, cross category schematic, N-cross interpretation (Step 22)
-   `SRK_HV_cluster_figure.pdf` / `figures/SRK_HV_cluster_figure.png` - UPGMA dendrogram (HV distances) coloured by class + AAAA availability bar chart (Step 22)
-   `SRK_cross_result_analysis_HV.pdf` - Seed yield distributions and success rates by cross category, with Kruskal-Wallis and Mann-Whitney U tests (Step 23; requires cross data)

### Quality Control Reports

-   `SRK_individual_status_report.tsv` - Comprehensive processing summary
-   `SRK_zygosity_distribution.pdf` - Genotype frequency visualization

## Applications

This pipeline is designed for:

-   **Population genetic analysis** of self-incompatibility systems
-   **Conservation genetics** assessment in threatened species
-   **Demographic inference** from S-locus diversity patterns
-   **Self-compatibility evolution** studies
-   **Polyploid genetics** research

## Case Study: *LEPA* (Threatened Brassicaceae)

For a full worked example of results produced by this pipeline, see [LEPA_SRK_report.md](LEPA_SRK_report.md).

Key findings from **272 individuals across five Element Occurrences** (EO25, EO27, EO67, EO70, EO76; Libraries 001–009):
- **54 observed S-allele bins** species-wide, with a predicted total of **71** (consensus of Michaelis-Menten = 69 and Chao1 = 73)
- All Element Occurrences retain only 9–32% of the species-level SI repertoire (EO70: 9%; EO27: 32%)
- Allele frequencies are significantly skewed from NFDS expectations at every level (χ² *p* < 10⁻⁷)
- 60% of individuals carry an AAAA genotype (single allele, four copies), flagging a high risk of reduced SI function
- A further 16% carry an AAAB genotype (GFS = 0.500) — dosage-imbalanced individuals that appear heterozygous but produce fewer diverse gametes than AABB individuals (GFS = 0.667)
- All five EOs are flagged **CRITICAL** for Tipping Point 2 (mean GFS 0.18–0.32; proportion AAAA 53–69%); EO67 is least degraded; EO76 is most degraded with no AABC individuals
- S-allele sets are largely private to each EO: only **2 alleles** (Allele_050 and Allele_057) are shared across all five EOs; EO27 holds the most private alleles (12), EO70 the fewest (2)
- Managed crossing simulations show 77–95% variance reduction within one generation vs. random mating, with zero allele loss under optimised preservation strategies

## Citation

TBD

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

For questions or issues:

-   Open an [issue](https://github.com/svenbuerki/SRK_bioinformatics/issues)
-   Contact: [[svenbuerki\@boisestate.edu](mailto:svenbuerki@boisestate.edu){.email}]

## Acknowledgments

TBD

------------------------------------------------------------------------

**Note:** This pipeline is specifically optimized for SRK analysis in Brassicaceae species but may be adaptable to other self-incompatibility systems with appropriate modifications.
