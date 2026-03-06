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
-   **Tetraploid zygosity analysis** for complex polyploid genetics
-   **Population genetic dataset generation** for demographic inference

## Pipeline Workflow

The pipeline consists of **14 main steps organized into three phases**, progressing from within-library sequence processing to cross-library integration and final population genetic analyses.

## Phase 1: Within-Library Processing

1.  **Prepare Canonical Sequences** – Preparation of reference SRK sequences used for downstream assembly and validation.
2.  **Nanopore Amplicon Assembly and Phasing Pipeline** – Multi-CANU assembly of Nanopore amplicons followed by haplotype phasing to reconstruct allelic sequences.
3.  **Haplotype Orientation and Consolidation Pipeline** – Standardization of sequence orientation and consolidation of phased haplotypes.
4.  **FASTA Sequence Filtering and Reference Integration Pipeline** – Quality filtering of sequences and integration of validated reference alleles.
5.  **Multiple Sequence Alignment** – Alignment of nucleotide sequences using MAFFT.
6.  **Exon Extraction from Multiple Sequence Alignments Pipeline** – Extraction of coding regions based on AUGUSTUS gene annotations.
7.  **MSA Gap Backfilling and Terminus Processing Pipeline** – Reference-guided backfilling of terminal alignment gaps and sequence terminus correction.
8.  **DNA to Amino Acid Translation Pipeline** – Frame-specific translation of nucleotide sequences into proteins. **Note:** This step is optional and may be skipped if protein sequences are not required.

## Phase 2: Cross-Library Analysis

9.  **SRK Protein Translation, Alignment, and Abundance Filtering Pipeline** – Translation of SRK alleles into proteins, alignment, and filtering based on abundance thresholds to retain functional candidates.
10. **SRK Functional Protein Genotyping Pipeline** – Construction of individual-level genotype matrices and detection of functional SRK proteins across samples.
11. **SRK Zygosity Analysis Pipeline for Tetraploid Species** – Classification of individuals as homozygous or heterozygous and estimation of population-level zygosity statistics in tetraploid species.

## Phase 3: Data Analyses

12. **Population Genetics Statistics** – Estimation of population-level diversity metrics, including heterozygosity, mean protein counts per individual, total allele counts, and effective allele numbers.
13. **Allele Accumulation Curves** – Rarefaction-based analysis of SRK allele discovery across individuals to evaluate patterns consistent with negative frequency-dependent selection versus genetic drift.
14. **Allele Frequency Analysis** – Species- and population-level χ² tests of allele frequency distributions to assess deviations from equal-frequency expectations under negative frequency-dependent selection.

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
install.packages(c("dplyr", "ggplot2"))
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

3.  **Run cross-library analysis:**

``` bash
# Steps 9-11: Protein analysis and genotyping
python scripts/09_srk_protein_filtering.py
python scripts/10_functional_genotyping.py
Rscript scripts/11_zygosity_analysis.R
```

### Detailed Usage

See the [Pipeline Documentation](https://svenbuerki.github.io/SRK_bioinformatics/) for comprehensive step-by-step instructions.

## Input Data Requirements

-   **Nanopore amplicon sequencing data** (FASTQ format)
-   **Canonical reference sequences** (FASTA format)
-   **Gene structure annotations** (AUGUSTUS CSV format)
-   **Organized directory structure** by library and barcode

## Output Files

### Key Outputs

-   `SRK_functional_proteins.fasta` - Validated SRK protein sequences
-   `SRK_individual_genotypes.tsv` - Population genetic matrix
-   `SRK_individual_zygosity.tsv` - Zygosity classifications
-   `SRK_self_compatible_candidates.txt` - Potentially self-compatible individuals

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
