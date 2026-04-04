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
-   **Population genetic dataset generation** for demographic inference
-   **Species allele richness estimation** using Michaelis-Menten, Chao1, and iNEXT estimators to quantify total SRK diversity and population deficits relative to the species optimum

## Pipeline Workflow

The pipeline consists of **15 main steps organized into three phases**, progressing from within-library sequence assembly to cross-library integration and final population genetic analyses.

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

13. **Population Genetics Statistics** – Estimation of population-level diversity metrics, including heterozygosity, mean alleles per individual, total allele counts, and effective allele numbers. Allele frequencies are based on copy counts summed across individuals (from the count matrix), giving a proper tetraploid frequency estimate.
14. **Allele Accumulation Curves** – Rarefaction-based analysis of SRK allele discovery across individuals to evaluate patterns consistent with negative frequency-dependent selection versus genetic drift. Includes estimation of total species allele richness using Michaelis-Menten asymptote fitting, Chao1, and iNEXT estimators. Outputs an empirical species optimum used as a baseline in step 15.
15. **Allele Frequency Analysis** – Species- and population-level χ² tests of allele frequency distributions to assess deviations from equal-frequency expectations under NFDS. The estimated species allele richness from step 14 is used as the optimum, quantifying how many alleles each population is missing relative to the species pool.
16. **Allele Composition Comparison Across Element Occurrences** – UpSet plot and pairwise sharing heatmap quantifying how S-allele sets partition across Element Occurrences, identifying private alleles and alleles shared across all populations. Requires only standard Python dependencies (pandas, numpy, matplotlib).

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
install.packages(c("dplyr", "ggplot2", "bookdown"))

# Optional — required for iNEXT richness estimation in step 14:
install.packages("iNEXT")
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

2.  **Run the within-library processing:**

``` bash
# Step 2: Assembly and phasing
./scripts/01_nanopore_assembly_phasing.sh

# Step 3: Orientation
./scripts/02_haplotype_orientation.sh

# Continue with remaining steps...
```

3.  **Run cross-library analysis (Phase 2):**

``` bash
# Steps 9-12: Protein filtering, allele definition, genotyping, zygosity
python scripts/09_srk_protein_filtering.py
Rscript scripts/10_srk_allele_definition.R
Rscript scripts/11_srk_allele_genotyping.R
Rscript scripts/12_zygosity_analysis.R
```

4.  **Run population genetic analyses (Phase 3):**

``` bash
# Step 13: Population genetics statistics
Rscript SRK_popgen_statistics.R

# Step 14: Allele accumulation curves and richness estimation
# (must run before step 15)
Rscript SRK_allele_accumulation_analysis.R

# Step 15: Allele frequency analysis (reads SRK_species_richness_estimates.tsv)
Rscript SRK_chisq_species_population.R

# Step 16: Allele composition comparison across Element Occurrences
python SRK_allele_sharing_EOs.py
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

### Population Genetic Outputs (Phase 3)

-   `SRK_allele_accumulation_curves.pdf` - Species- and population-level accumulation curves with estimated asymptotes
-   `SRK_allele_accumulation_stats.tsv` - Curve statistics including MM, Chao1, and iNEXT richness estimates per level
-   `SRK_species_richness_estimates.tsv` - Consensus species allele richness estimate (input to step 15)
-   `SRK_chisq_species_population.tsv` - χ² test statistics at species and population levels
-   `SRK_chisq_species_population_frequency_plots.pdf` - Allele frequency plots showing observed distribution, NFDS expectation, and missing alleles relative to estimated optimum
-   `SRK_allele_upset_EOs.pdf` - UpSet plot of all pairwise and higher-order allele set intersections across Element Occurrences
-   `SRK_allele_sharing_heatmap_EOs.pdf` - Pairwise allele sharing heatmap between Element Occurrences

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

Key findings from **189 individuals across five Element Occurrences** (EO25, EO27, EO67, EO70, EO76):
- **47 observed S-allele bins** species-wide, with a predicted total of **75** (consensus of Michaelis-Menten = 65 and Chao1 = 84)
- All Element Occurrences retain only 11–40% of the species-level SI repertoire (EO70: 11%; EO27: 40%)
- Allele frequencies are significantly skewed from NFDS expectations at every level (χ² *p* < 10⁻⁷)
- 56% of individuals carry an AAAA genotype (single allele, four copies), flagging a high risk of reduced SI function
- S-allele sets are largely private to each EO: only **2 alleles** (Allele_044 and Allele_048) are shared across all five EOs; EO27 holds the most private alleles (10), EO70 the fewest (3)
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
