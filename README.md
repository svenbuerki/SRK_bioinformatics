# SRK Bioinformatics Pipeline

## Overview

This repository contains a comprehensive bioinformatics pipeline for analyzing S-receptor kinase (SRK) haplotype diversity in self-incompatible plant species using Oxford Nanopore long-read sequencing. The pipeline is specifically designed for threatened plant species research, enabling assessment of self-incompatibility system integrity in populations where genetic diversity loss may compromise reproductive success.

## Author

**Sven Buerki, Ph.D.** (he/him)\
Associate Professor\
Dr. Christopher Davidson Endowed Chair in Botany\
Department of Biological Sciences\
Boise State University\
Boise, Idaho, USA

📧 [svenbuerki\@boisestate.edu](mailto:svenbuerki@boisestate.edu)

## Key Features

-   **Multi-replicate CANU assembly** optimized for rare haplotype recovery
-   **Polyploid-aware variant calling** and read-backed phasing with WhatsHap
-   **Reference-guided sequence orientation** and gap-filling procedures
-   **SRK-specific protein validation** with optional domain analysis
-   **Abundance-based filtering** to distinguish authentic alleles from technical artifacts
-   **Tetraploid zygosity analysis** for complex polyploid genetics
-   **Population genetic dataset generation** for demographic inference

## Pipeline Workflow

The pipeline consists of 11 main steps divided into two phases:

### Phase 1: Within-Library Processing (Steps 1-8)

1.  **Prepare Canonical Sequences** - Reference sequence preparation
2.  **Nanopore Amplicon Assembly and Phasing** - Multi-CANU assembly with polyploid phasing
3.  **Haplotype Orientation and Consolidation** - Standardize sequence orientation
4.  **FASTA Sequence Filtering and Reference Integration** - Quality control and reference addition
5.  **Multiple Sequence Alignment** - MAFFT-based alignment
6.  **Exon Extraction** - Extract coding sequences using AUGUSTUS annotations
7.  **MSA Gap Backfilling** - Reference-guided gap filling for terminal regions
8.  **DNA to Amino Acid Translation** - Frame-specific protein translation

### Phase 2: Cross-Library Analysis (Steps 9-11)

9.  **SRK Protein Translation and Filtering** - Comprehensive protein validation and abundance filtering
10. **Functional Protein Genotyping** - Generate population genetic matrices
11. **Zygosity Analysis** - Tetraploid genotype classification and demographic inference

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

- Open an [issue](https://github.com/svenbuerki/SRK_bioinformatics/issues) 
- Contact: [[svenbuerki\@boisestate.edu](mailto:svenbuerki@boisestate.edu){.email}]

## Acknowledgments

TBD

------------------------------------------------------------------------

**Note:** This pipeline is specifically optimized for SRK analysis in Brassicaceae species but may be adaptable to other self-incompatibility systems with appropriate modifications.
