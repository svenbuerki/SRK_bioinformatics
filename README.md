# SRK Bioinformatics Pipeline

## Overview

This repository contains a comprehensive bioinformatics pipeline for analyzing S-receptor kinase (SRK) haplotype diversity in self-incompatible plant species using Oxford Nanopore long-read sequencing. The pipeline is specifically designed for threatened plant species research, enabling assessment of self-incompatibility system integrity in populations where genetic diversity loss may compromise reproductive success.

### Roadmap — what the pipeline does, at a glance

The five-phase pipeline takes Nanopore amplicon reads through assembly, allele definition, population-genetic diagnostics, hypothesis-driven cross design, and per-individual SI-system characterisation. Each box below summarises what that phase does, the steps it covers, and the canonical output it produces — read left-to-right to follow the data flow.

![SRK pipeline overview — 5-phase roadmap from Nanopore amplicons through functional S-alleles to conservation diagnostics, cross plan, and forward-time simulation. Produced by `SRK_pipeline_overview_figure.R`.](figures/SRK_pipeline_overview.png)

Detailed step-by-step protocol in [`Bioinformatics_pipeline.md`](Bioinformatics_pipeline.md); full output organisation under [Output Organisation](#output-organisation-2026-06-14-refactor).

## Related Repositories

This pipeline integrates with a sibling project that provides the spatial framework for S-allele interpretation:

| Repository | Role |
|------------|------|
| **[`svenbuerki/LEPA_EO_spatial_clustering`](https://github.com/svenbuerki/LEPA_EO_spatial_clustering)** | Quantifies habitat fragmentation across the species range using germplasm collection events and a 500 m pollinator-dispersal threshold; identifies five **independent bottleneck lineages (BL1–BL5)** by Ward's D2 hierarchical clustering of geographic-group centroids. Produces `EO_group_BL_summary.csv`, the cross-reference key consumed by Step 13 of this pipeline. |
| **`svenbuerki/SRK_bioinformatics`** (this repo) | Reconstructs S-allele genotypes from Nanopore amplicons and stratifies all Phase 3/4 analyses by the BL framework defined above. |

The two repositories share a locked **RColorBrewer Set1 BL colour palette** (BL1 purple, BL2 blue, BL3 red, BL4 orange, BL5 green) so figures are visually consistent across both projects: the dendrogram, BL geographic-context map, and BL drift panel in the spatial-clustering repo, and the BL accumulation curves, drift erosion bars, TP1/TP2 scatters, UpSet plots, and per-BL entropy heatmap here.

### Shared ordering and palette config

All BL / EO ordering and colour decisions are centralised in two mirrored modules — **[`srk_bl_constants.R`](srk_bl_constants.R)** and **[`srk_bl_constants.py`](srk_bl_constants.py)** — sourced/imported by every pipeline script that orders or colours BLs. Both modules derive the orderings at load time from CSVs mirrored from the sibling project, so the entire pipeline reorders automatically when sampling grows — only the CSVs in `tables/` need updating.

Two BL orderings are exposed:

- **`BL_ORDER`** — high-to-low by total habitat area, with within-BL connectivity as the secondary tie-break (currently **BL4, BL5, BL3, BL1, BL2**). Sort key: `total_area_ha` (desc), then connectivity = `n_locations − n_groups` (desc), then BL name (asc). Area is the primary *Ne* proxy (carrying capacity → drift floor); connectivity is the secondary stratification for BLs of similar size. Used for **bars, faceted panels, horizontal-bar y-axes, and tables** — where BL appears as an axis category and the order itself tells a story.
- **`BL_ORDER_NUMERIC`** — alphanumeric BL1→BL5. Used for **scatter-plot legends** (TP1, TP2, accumulation-curve legends) where the axes are not BL and a readable legend matters more than the inferential order.

EOs are ordered within their parent BL by ascending mean Drift_index (lower DI = more connected = first), with EO name as deterministic tie-break, via `get_eo_order_within_bl()`. The same formula is implemented identically in R and Python so the two languages produce identical orderings.

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
-   **Targeted QC filters** addressing two failure modes observed in Library 010 processing: (i) a **reference-similarity (BLAST coverage) filter** that removes SRK paralogs carrying a localised ~500 bp insertion (98–100 % identical to canonical SRK where they align, so identity-based filtering fails); (ii) an **N-content filter** that removes contigs with N-rich termini introduced by Canu *de novo* assembly when Nanopore reads are shorter than the ~3.5 kb amplicon — the upstream cause of spurious gap blocks in downstream AA alignments.
-   **SRK-specific protein validation** with optional domain analysis
-   **Abundance-based filtering** to distinguish authentic alleles from technical artifacts
-   **Tetraploid zygosity analysis** with allele copy-count dosage inference (AAAA, AAAB, AABB, AABC, ABCD) for complex polyploid genetics
-   **Individual Genotypic Fitness Score (GFS)** quantifying heterozygous gamete proportion per individual — differentiates AABB (GFS = 0.667) from AAAB (GFS = 0.500) to guide seed parent selection and EO prioritisation
-   **Population genetic dataset generation** for demographic inference
-   **Species allele richness estimation** using Michaelis-Menten, Chao1, and iNEXT estimators to quantify total SRK diversity and population deficits relative to the species optimum

## Pipeline Workflow

The pipeline consists of **28 main steps organized into five phases**, progressing from within-library sequence assembly to cross-library integration, population genetic analyses, and experimental validation of S-allele hypotheses through controlled crosses.

## Phase 1: SRK Amplicon Sequence Assembly

1.  **Prepare Canonical Sequences** – Preparation of reference SRK sequences used for downstream assembly and validation.
2.  **Nanopore Amplicon Assembly and Phasing Pipeline** – Multi-CANU assembly of Nanopore amplicons followed by haplotype phasing to reconstruct allelic sequences.
3.  **Haplotype Orientation and Consolidation Pipeline** – Standardization of sequence orientation and consolidation of phased haplotypes.
4.  **FASTA Sequence Filtering and Reference Integration Pipeline** – Quality filtering of sequences by length (min 3250 bp, max 4000 bp) and integration of validated reference alleles.
4b. **Reference-Similarity Filter (BLAST + coverage)** – Removes SRK paralogs that pass the length filter but contain a localised ~500 bp insertion. Paralogs are 98–100 % identical to canonical SRK *where they align*, so percent identity does not discriminate — the filter uses BLAST query coverage (≥0.90) against the canonical LEPA SRK reference set. Without this filter, paralogs inflate the Step 5 MAFFT alignment and produce spurious gap blocks downstream.
5.  **Multiple Sequence Alignment** – Alignment of nucleotide sequences using MAFFT.
6.  **Exon Extraction from Multiple Sequence Alignments Pipeline** – Extraction of coding regions based on AUGUSTUS gene annotations.
7.  **MSA Gap Backfilling and Terminus Processing Pipeline** – Reference-guided backfilling of terminal alignment gaps and sequence terminus correction.
7b. **Ambiguity (N-content) Filter** – Removes contigs with N-rich termini. **Origin of the `N`s:** when Nanopore reads are shorter than the ~3.5 kb SRK amplicon (DNA fragmentation during library prep, partial pore translocation, or end-of-read truncation), Canu *de novo* assembly cannot span the full amplicon and produces contigs shorter than the canonical reference. Step 7 then pads the missing terminus with `N`. Downstream, those `N` codons translate to `X` residues in Step 8 and cause MAFFT to insert spurious gap blocks on the BEA canonical reference rows. Filtering at the DNA stage prevents the cascade.
8.  **DNA to Amino Acid Translation Pipeline** – Frame-specific translation of nucleotide sequences into proteins. **Note:** This step is optional and may be skipped if protein sequences are not required.

## Phase 2: Functional Proteins, S Alleles and Genotyping

9.  **SRK Protein Translation, Alignment, and Abundance Filtering Pipeline** – Translation of SRK alleles into proteins, alignment, and filtering based on abundance thresholds to retain functional candidates.
10. **Distance-Based SRK S-Allele Definition** – Grouping of functional proteins into allele bins using pairwise amino acid p-distance on the S-domain ectodomain. Calibration is a separate first step: `find_allele_plateau.py` applies three independent estimators (Kneedle elbow detection, slope-magnitude minimum, longest-persistent-run plateau) to the sensitivity curve and reports a recommended `N_ALLELES` on stdout, which is then set manually in `define_SRK_alleles_from_distance.py` before clustering. Step 10b (`SRK_AA_mutation_heatmap.py`) produces amino acid variation heatmaps for the resulting allele bins.
11. **SRK S-Allele Genotyping Pipeline** – Assignment of distance-defined S-allele bins to individuals, producing long-format allele tables (with allele copy counts) and count-based genotype matrices. Each cell in the wide matrix records the number of distinct proteins from that individual assigned to that allele bin, capturing dosage information for downstream zygosity inference and allele frequency analysis.
12. **SRK Zygosity Analysis Pipeline for Tetraploid Species** – Classification of individuals as homozygous or heterozygous using the allele count matrix. Infers tetraploid genotype classes (AAAA, AAAB, AABB, AABC, ABCD) from both the number of distinct allele bins and their copy counts, and outputs an `Allele_composition` string (e.g. `Allele_044(2)+Allele_047(2)`) for direct use in cross planning.

**⚠️ Step 12b (In Progress) — SRK Class Assignment and Dominance Prediction** – Phylogenetic classification of SRK proteins into Class I (dominant) or Class II (recessive) using maximum-likelihood inference (IQ-TREE3). Automated class assignment by longest internal branch exploits the ancient divergence between SRK classes maintained by balancing selection. Per-individual dominance predictions guide cross-compatibility decisions in a conservation breeding program. In the current LEPA dataset (2026-05-11) the procedure assigns Allele_055 to Class II; all other observed alleles fall in Class I. **This analysis is under active development and may not function correctly on all datasets — results require validation against published SRK reference sequences with known class assignments.**

**Step 12c — Data Quality Evaluation (Phase 3 QC gate)** – Cross-tabulates every metadata sample against its terminal pipeline outcome, producing five mutually exclusive categories (Functional, Partial translation failure, SI-escape candidate, Re-PCR, Re-DNA-extraction). Writes three deliverables: the per-sample categorisation, the lab follow-up sheet (Re-PCR + Re-DNA-extraction), and the SI-escape candidate list for phenotyping. Per-BL and per-EO stacked-bar figures use the connectivity order. Implemented in `evaluate_data_quality.py`. **Execution order (revised 2026-06-13):** 12c.ii + 12c.iii run **after Step 13** (BL integration) so per-BL stratification is available; 12c.i (library-effect tests) runs **after Step 19** (GFS) since it consumes `tables/Phase3/step19_individual_GFS.tsv`.

## Phase 3: Data Analyses

13. **Element Occurrence and Bottleneck Lineage Integration** – Bridges the SRK individual-level dataset to the spatial framework defined in the sibling project [LEPA_EO_spatial_clustering](https://github.com/svenbuerki/LEPA_EO_spatial_clustering). For each post-QC individual, joins EO label, geographic group, bottleneck lineage (BL1–BL5), and drift index from `EO_group_BL_summary.csv`. Produces `tables/Phase3/step13_individual_BL_assignments.tsv`, the join key used by all subsequent Phase 3 steps to add BL stratification (EO sorted within BL for management; BL aggregated for inferential power on small localities).
13b. **Sampling overview map** – Plots all catalogued populations of the species against their BL membership, distinguishing populations successfully sampled by the SRK pipeline (filled triangles) from those not yet sampled (filled circles). Focus EOs used in most downstream analyses (EO18, EO25, EO27, EO67, EO70, EO76) are labelled in white-fill boxes. JAR-region germplasm sub-codes (no resolved coordinates) are not plotted. Implemented in `SRK_sampling_map.R`; reads `tables/Phase2/step11_individual_allele_genotypes.tsv`, `tables/sampling_metadata.csv`, `tables/Phase3/step13_individual_BL_assignments.tsv`, and `EO_location_groups.csv` mirrored from the spatial-clustering project.
14. **Population Genetics Statistics** – Estimation of population-level diversity metrics, including heterozygosity, mean alleles per individual, total allele counts, and effective allele numbers. Allele frequencies are based on copy counts summed across individuals (from the count matrix), giving a proper tetraploid frequency estimate.
15. **Allele Accumulation Curves** – Rarefaction-based analysis of SRK allele discovery across individuals to evaluate patterns consistent with negative frequency-dependent selection versus genetic drift. Includes estimation of total species allele richness using Michaelis-Menten asymptote fitting, Chao1, and iNEXT estimators. Outputs an empirical species optimum used as a baseline in steps 16 and 17.
16. **Allele Frequency Analysis** – Species- and population-level χ² tests of allele frequency distributions to assess deviations from equal-frequency expectations under NFDS. The estimated species allele richness from step 15 is used as the optimum, quantifying how many alleles each population is missing relative to the species pool.
17. **TP1 — Mating-pool functionality** – Reframed diagnostic answering two complementary conservation questions: (i) *is random mating still viable in this population?* via the **compatible-pair fraction** P_compat (the fraction of random plant pairs that are cross-compatible under tetraploid sporophytic SI), rendered as a red/amber/green traffic-light figure per focal EO with bootstrap 95 % CIs; and (ii) *which intervention does the data support?* via P_compat × Depletion Index (DI = 1 − k_group / k_species, against the MM consensus species optimum k_species = 59) with quadrants labelled by the three breeding strategies in `tables/SRK_breeding_strategies.csv` (HEALTHY / INFORMED BREEDING (frequency skew) / INFORMED BREEDING + ALLELE INJECTION (preventive / urgent)). Scripts: `SRK_TP1_compatibility_metrics.py` (computes per-EO + per-BL metrics with bootstrap CIs), `SRK_P_compat_traffic_light.R` (stakeholder traffic-light figure), `SRK_depletion_ranking.R` (P_compat × DI conservation-ranking figure with `_observed` / `_predicted` panels and `_blank` / `_EOs` / `_all` variants for layered presentation builds).
18. **Allele Composition Comparison Across Element Occurrences** – UpSet plot and pairwise sharing heatmap quantifying how S-allele sets partition across Element Occurrences, identifying private alleles and alleles shared across all populations. Requires only standard Python dependencies (pandas, numpy, matplotlib).
19. **Individual Genotypic Fitness Score (GFS)** – Per-individual metric quantifying the proportion of heterozygous diploid gametes a tetraploid can produce. Differentiates dosage-imbalanced genotypes (AABB vs AAAB) invisible to zygosity analysis alone. Outputs per-individual GFS values and ranked seed parent lists per EO.
20. **TP2 Tipping Point Analysis** – EO-level assessment placing mean GFS and proportion of AAAA individuals in interaction. EOs breaching both thresholds simultaneously (mean GFS < 0.667; proportion AAAA > 30%) are flagged CRITICAL, AT RISK, or OK.
21. **Reproductive Effort Support per Element Occurrence** – Horizontal proportional bar chart showing the fraction of individuals at each GFS tier per EO.

## Phase 4: Per-individual SI status, null-aware genotypes, and forward simulation (Steps 25–28)

*Phase 4 must precede Phase 5: cross design is only meaningful once we know each individual's SI status (full SI, partial SI, or self-compatible), because a broken-SI individual cannot serve as a reliable parent in a cross-compatibility prediction.*

25. **Per-individual SI / pSI / SC reconstruction** – `SRK_individual_SI_status.py` (Step 25b) re-aggregates per-haplotype OK / REMOVED calls from the Step-7 `_frame1_stopcodon_log.tsv` files to classify each individual as **SI** (all 4 copies functional), **pSI** (1–3 copies broken — partial SI loss), **SC** (all 4 copies broken — full self-compatibility), or **Insufficient_data** (too few haplotypes after chimera filtering to call status). `SRK_null_allele_assignment.py` (Step 25a) optionally augments this by mapping REMOVED haplotypes to their nearest functional allele via AA-distance for chimera-aware identity assignment. `SRK_SI_status_figures.R` produces species-, BL-, and EO-level stacked-bar figures with `_full` (7-tier QC view) and `_robust` (5-tier defensible-call view) variants. **Current dataset: 247 SI / 15 pSI / 1 SC / 234 Insufficient_data.**
26. **Null-aware genotype rebuild** – `SRK_genotype_null_alleles.py` rewrites the Step 11 allele × individual matrix so that pSI individuals carry explicit `Allele_NULL` copies (proportional to `copies_nonfunctional`) and the single SC individual is represented as `Allele_NULL × 4`; rows sum to 4 across the 49 functional alleles + NULL. Outputs `tables/Phase4/step26_individual_allele_genotypes_with_nulls.tsv` + a redo flag list (`step26_samples_for_redo.tsv`) for the 234 Insufficient_data individuals. Triggers the **null-aware cascade re-runs**: Step 14b (pop gen), Step 17b (TP1 P_compat + DI ranking), Step 19b + 20b (GFS + TP2) — same scripts as the canonical Steps 14/17/19/20 but reading the null-aware genotype matrix.
27. **Forward-time tetraploid inheritance simulator** – `SRK_inheritance_simulator.py` runs a forward-time multi-replicate simulation per BL under multiple scenarios (baseline / high_drift / rescue_high / rescue_low) starting from the empirical Step 26 genotypes. Tracks per-copy NULL frequency, SC frequency, and time-to-50%-SC trajectories. `SRK_inheritance_figures.R` produces the trajectory + SC-progression + time-to-SC figures.
28. **AA-distance null-allele assignment + donor ranking** – `SRK_injection_donor_ranking.py` ranks candidate seed parents for allele-injection crosses to rescue specific BLs / EOs. `SRK_donor_ranking_figure.R` produces the donor-recovery-ladder figure (`figures/Phase4/step28_donor_recovery_ladder.{pdf,png}`).

## Phase 5: Testing S-allele Hypotheses and Cross Design (Steps 22–23)

*Phase 5 uses the per-individual SI status established in Phase 4: only individuals classified as SI (or pSI · low confidence) are valid parents for the cross plan, ensuring that experimental compatibility predictions rest on functioning recognition machinery.*

22. **HV-Based Allele Hypothesis Testing and Crossing Design** – Multi-script workflow. **Step 22a** (`srk_variability_landscape.py`): combined alignment of LEPA + *Brassica rapa*/*B. oleracea* + *Arabidopsis lyrata*/*A. halleri* SRK alleles; identical sliding-window Shannon-entropy scan applied to each species; per-species HV regions called at mean + 1×SD; anchored structurally by 11 of 12 mappable SCR9-contact residues from Ma et al. 2016 (PDB 5GYY); canonical 66 LEPA HV columns written to `tables/Phase5/step22a_LEPA_HV_positions.tsv` (in the current post-QC dataset, the cross-genus HV-overlap permutation is no longer significant for LEPA pairs — interpreted as drift-eroded standing variation, with Brassica↔Arabidopsis remaining highly significant). **Step 22b** (`srk_allele_hypotheses.py`): reads the canonical HV columns, computes HV-only pairwise distances, UPGMA-clusters into Class I / Class II, builds the synonymy network (HV-identical synonymy groups + Synonymy_test bridges), and generates the Incompatible / Synonymy_test / Compatible_within / Compatible_cross cross design — *restricted to parents whose Phase 4 SI status is `SI`*. **Steps 22c and 22d** are optional mechanism diagnostics. **Step 22e** (`srk_cross_plan.py`) is the operational deliverable — a phased cross plan testing five nested hypotheses (H0 SI validation; H1a within-Class baseline; H1b between-Class baseline; H2 synonymy bin boundaries; H3 hidden bins via heterozygous donors) with explicit AAAA / AAAB / AABB genotype requirements + paired controls for AAAB-mediated tests.
23. **Cross Result Analysis** – Reads completed crossing records and tests whether the cross category (Incompatible / Synonymy_test / Compatible_within / Compatible_cross) predicts seed yield, validating the sequence-based allele definitions against experimental cross-compatibility data. Activated by setting `CROSS_TSV` in the script once crossing data are available.

## Output Organisation (2026-06-14 refactor)

Every canonical pipeline script (Step 9 → 28) writes to a **Phase- and step-prefixed path**:

```
FASTA/                  ← external reference FASTAs (Brassica + Arabidopsis SRK)
tables/sampling_metadata.csv   ← canonical sample registry
tables/Phase{2,3,4,5}/stepXX_*.{tsv,csv,fasta}
figures/Phase{2,3,4,5}/stepXX_*.{pdf,png}
```

Every plot produces both `.pdf` (vector) and `.png` (raster) except inherently multi-page PDFs (Step 15 curves, Step 16 frequency, Step 19 GFS) where the page-level PNGs are companions. See `Bioinformatics_pipeline.md` for the full directory tree.

**Phase 1 wrapper:** `Scripts/run_within_library_suite.sh` chains Steps 2 → 8 for every library listed in a TSV config (auto-generated from `Library*/barcode??` if absent), pipes each script's prompts automatically, skips already-completed samples, and continues past per-library failures. Env vars: `CHIMERA_FILTER_MODE`, `MIN_UNIFORMITY` (default 0.2 — Library 005 calibration), `FORCE_SAMPLE`, `FORCE_CANU`, `START_AT_STEP`, `STOP_AT_STEP`.

**Phase 4 LEPA-only filter (2026-06-14):** `pad_representatives.py` now restricts Step 22 inputs to the 49 LEPA-ingroup-observed alleles (drops 6 outgroup-only alleles from *L. montanum* / *L. freemontii* / *L. philonitron*) so the cross-Brassicaceae variability analysis isn't contaminated by non-LEPA sequences.

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

> **Recommended:** run the wrapper `run_within_library_suite.sh`, which
> chains Steps 2 → 8 for every library listed in `libraries.tsv` (auto-generated
> from `Library*/barcode??` if absent), feeds each script's prompts automatically,
> skips already-completed samples / steps, and continues past per-library failures.
> See `Bioinformatics_pipeline.md` for env vars and resume semantics.
>
> ```bash
> ./run_within_library_suite.sh                         # report-mode chimera survey
> CHIMERA_FILTER_MODE=filter ./run_within_library_suite.sh   # production
> ```

The individual scripts (below) are still callable directly if you prefer to run a single step.

``` bash
# Step 2: Assembly and phasing (now failure-tolerant, with embedded
#   coverage-based chimera filter, Canu cache resume, and per-sample skip)
./nanopore_assembly_pipeline_barcode_range.sh

# Step 3: Orientation
./scripts/02_haplotype_orientation.sh

# Step 4: Length filter + canonical reference integration
./filter_and_add_references.sh

# Step 4b: Reference-similarity filter (BLAST + coverage)
#   Removes SRK paralogs (~500 bp insertion) that pass length filter
#   but inflate the MAFFT alignment downstream. New for Library 010+.
./filter_by_reference_similarity.sh

# Step 5: MAFFT multiple sequence alignment
mafft --auto --adjustdirection <_blastfilt.fasta>  > <_blastfilt_aligned.fasta>

# Step 6: Exon extraction from MSA
python extract_exons_with_annotation.py

# Step 7: Backfill terminal gaps from canonical reference
python backfill_alignment_ends.py

# Step 7b: Ambiguity (N-content) filter
#   Removes contigs with N-rich termini introduced by Step 7 when
#   Canu assembly cannot span the full amplicon (short Nanopore reads).
#   These would translate to X residues in Step 8 and cause spurious
#   gap blocks on BEA canonical references. New for Library 010+.
./filter_by_ambiguity.sh

# Step 8 (optional): DNA → AA translation, stop-codon filter, MAFFT
python translate_filter_align_AA.py
```

3.  **Phase 2 — Functional Proteins, S Alleles and Genotyping (Steps 9–12):**

``` bash
# Step 9: Protein filtering (define functional proteins from per-library backfilled DNA)
python define_functional_proteins.py

# Step 10a-pre: Calibrate N_ALLELES (Kneedle elbow on the sensitivity curve)
#   Re-run whenever the dataset changes. Prints the recommended N_ALLELES on stdout.
python find_allele_plateau.py

# Step 10a: Allele clustering at calibrated N_ALLELES (manually edit N_ALLELES at top of script)
python define_SRK_alleles_from_distance.py

# Step 10b: Amino acid variation heatmaps
python SRK_AA_mutation_heatmap.py

# Steps 11–12: Genotyping and tetraploid zygosity
python genotype_SRK_from_alleles.py
python compute_zygosity.py

# Step 12c: Data quality evaluation (end-of-Phase-2 gate)
#   Reads the exclusion audit, categorises every metadata sample, writes the
#   lab follow-up sheet (Re-PCR + Re-DNA-extraction) and the SI-escape list.
python evaluate_data_quality.py
```

4.  **Phase 3 — Data Analyses (Steps 13–21):**

``` bash
# Step 13: EO + bottleneck lineage integration (run first — produces the join key for all downstream steps)
python3 SRK_BL_integration.py

# Step 13b: Sampling overview map (sampled triangles vs unsampled circles, coloured by BL)
#   Requires Step 13 output + EO_location_groups.csv from the sibling spatial-clustering project
Rscript SRK_sampling_map.R

# Step 14: Population genetics statistics
Rscript SRK_population_genetic_summary.R

# Step 15: Allele accumulation curves and richness estimation
# (must run before steps 16 and 17)
Rscript SRK_allele_accumulation_analysis.R

# Step 16: Allele frequency analysis (reads SRK_species_richness_estimates.tsv)
Rscript SRK_chisq_species_population.R

# Step 17: TP1 — mating-pool functionality (two scripts)
#   metrics.py computes per-EO + per-BL evenness J, P_compat (L=0..0.5),
#   bootstrap 95% CIs on P_compat, L_hat, k_rarefied30, and writes
#   tables/Phase3/step17_EO_allele_richness.tsv;
#   traffic_light.R renders the stakeholder-facing red/amber/green
#   priority figure (full + blank variants);
#   depletion_ranking.R renders the P_compat × DI conservation-ranking
#   figure (observed + predicted, blank / EOs / all variants).
python SRK_TP1_compatibility_metrics.py
Rscript SRK_P_compat_traffic_light.R
Rscript SRK_depletion_ranking.R

# Step 18: Allele composition comparison across Element Occurrences
#   Python UpSet plots + pairwise heatmaps (EO + BL); R companion produces
#   an area-proportional 5-set Euler diagram of BL allele sets for
#   stakeholder presentations. Both scripts read the same Step 11 + Step 13
#   inputs, so regenerate together whenever the dataset grows.
python SRK_allele_sharing_EOs.py
Rscript SRK_allele_eulerr_BLs.R

# Steps 19–20: Individual Genotypic Fitness Score (Step 19) and TP2 tipping point analysis (Step 20)
Rscript SRK_individual_GFS.R

# Step 21: Reproductive effort support per EO and BL (with AAAA allele identity panels)
Rscript SRK_GFS_reproductive_effort.R
```

5.  **Phase 4 — Testing S-allele Hypotheses (Steps 22–23):**

``` bash
# Step 22a: cross-Brassicaceae variability landscape (Shannon entropy + permutation test)
#
# Rebuild chain (run every time Step 10a re-runs, e.g. new library, new N_ALLELES):
python3 srk_fetch_reference_alleles.py            # ONCE per project (re-run only if reference set changes)
python3 pad_representatives.py                    # pads LEPA reps to uniform length for mafft --add
mafft --add all_reference_SRKs_dedup.fasta SRK_protein_allele_representatives_padded.fasta \
      > SRK_combined_alignment.fasta              # rebuilds combined alignment
python3 srk_brassica_hv_mapping.py                # remaps Ma 2016 SCR9 contacts to new LEPA coordinates
# Then run the variability landscape:
python3 srk_variability_landscape.py

# Step 22b: HV-based allele hypothesis testing and crossing design
python3 srk_allele_hypotheses.py

# Step 22c (optional diagnostic): does Synonymy group redundancy explain LEPA's low Shannon entropy?
python3 srk_wgroup_collapse_test.py

# Step 22d (optional diagnostic): drift vs selection at LEPA HV cols (per-BL entropy + cross-genera dominant residues)
python3 srk_perBL_entropy_test.py

# Step 22e: hypothesis-testing cross plan generator (H0/H1a/H1b/H2/H3 phased protocol with genotype constraints)
python3 srk_cross_plan.py

# Step 23: Cross result analysis
# Set CROSS_TSV = "<your_cross_results_file>" in the script, then re-run:
python srk_allele_hypotheses.py

# Step 25: per-individual SI system status (SI / pSI / SC / Insufficient_data)
# Answers: "What is the status of the SI system?" at individual, BL, and EO levels.
# Re-aggregates per-haplotype OK/REMOVED calls from the Step-7 stop-codon logs
# (which Phase 1-2 drops before genotyping) without invalidating the canonical
# 49-allele catalogue. Applies a Canu-noise filter (n_REMOVED < 2 → SI) and a
# pSI confidence flag (high if frac_NF ≥ 0.25).
python3 SRK_individual_SI_status.py
Rscript SRK_SI_status_figures.R

# Step 26: null-aware tetraploid genotype rebuild + downstream reruns
# Propagates Step 25 SI status into the genotype tables: pSI_high gets explicit
# Allele_NULL copies, SC included with all-null genotype, Insufficient_data
# flagged for re-sequencing. Original Phase 1-2 tables stay frozen.
python3 SRK_genotype_null_alleles.py
Rscript SRK_population_genetic_summary_with_nulls.R       # Step 14b
/Users/sven/anaconda3/bin/python SRK_TP1_compatibility_metrics_with_nulls.py  # Step 17b
Rscript SRK_individual_GFS_with_nulls.R                   # Steps 19b + 20b
Rscript SRK_P_compat_traffic_light_with_nulls.R           # stakeholder traffic-light figure
Rscript SRK_depletion_ranking_with_nulls.R                # DI x P_compat ranking

# Step 27: forward-time inheritance simulator + figures.
# Answers: "Where is each BL heading on the SI→SC erosion axis under
# current conditions, and what conservation lever stops the trajectory?"
# Tetraploid Wright-Fisher with NULL-aware sporophytic SI; 4 scenarios
# (baseline, rescue_low, rescue_high, high_drift) × 5 BLs × 30 replicates.
/Users/sven/anaconda3/bin/python SRK_inheritance_simulator.py --n_generations 100 --n_replicates 30
Rscript SRK_inheritance_figures.R
```

### Detailed Usage

See the [Pipeline Documentation](https://svenbuerki.github.io/SRK_bioinformatics/) for comprehensive step-by-step instructions.

For a concise step-by-step protocol (scripts, inputs, outputs, key parameters), see [Bioinformatics_pipeline.md](Bioinformatics_pipeline.md).

## Input Data Requirements

-   **Nanopore amplicon sequencing data** (FASTQ format)
-   **Canonical reference sequences** (FASTA format)
-   **Gene structure annotations** (AUGUSTUS CSV format)
-   **Organized directory structure** by library and barcode

## Curated tables (`tables/`)

The `tables/` folder is the single source of truth for the most-shared CSV/TSV deliverables — downloadable directly by modeling, lab, and conservation collaborators without re-running the pipeline.

### Sample registry & population key

-   **[`tables/sampling_metadata.csv`](tables/sampling_metadata.csv)** — canonical sample registry. Links each `SampleID` (e.g., `Library001_barcode02`, joinable to the `Individual` column of every genotyping table below) to its `Pop`, `EO_w_sub`, `Library`, `barcode`, `FlowCell_ID`, `Ingroup` flag (1 = analyzed, 0 = excluded), and `OccurrenceID`.
-   **[`tables/EO_group_BL_summary.csv`](tables/EO_group_BL_summary.csv)** — EO → geographic group → bottleneck lineage (BL1–BL5) → drift index cross-reference. Generated by the sibling repository [LEPA_EO_spatial_clustering](https://github.com/svenbuerki/LEPA_EO_spatial_clustering); consumed by Step 13 of this pipeline.

### Genotyping deliverables (for modeling collaborators)

-   **[`tables/Phase2/step11_individual_allele_genotypes.tsv`](tables/Phase2/step11_individual_allele_genotypes.tsv)** — wide allele×individual count matrix (367 individuals × 49 alleles; integer copy counts). Ready-to-use design matrix for allele frequency estimation, drift/erosion modeling, simulations.
-   **[`tables/Phase2/step12_individual_zygosity.tsv`](tables/Phase2/step12_individual_zygosity.tsv)** — per-individual genotype summary with columns `N_distinct_alleles`, `N_total_proteins`, `Zygosity` (Homozygous / Heterozygous), `Genotype` (AAAA / AAAB / AABB / AABC / ABCD), `Allele_composition` (packed allele-with-counts string). Human-readable companion to the wide matrix above.
-   **[`tables/Phase5/step22b_synonymy_groups.csv`](tables/Phase5/step22b_synonymy_groups.csv)** — allele → synonymy group mapping with group size and observation counts. Lets modelers collapse putatively identical alleles into a single functional unit when needed.

**Recommended join chain for `individual → allele counts → population → bottleneck lineage`:**

```
tables/Phase2/step11_individual_allele_genotypes.tsv  (Individual)
       │ join on SampleID = Individual
       ▼
tables/sampling_metadata.csv  (Pop, EO_w_sub)
       │ join on EO
       ▼
tables/EO_group_BL_summary.csv  (EO → BL, Drift_index)
```

### Quality categorisation

-   **[`tables/Phase2/step12c_data_quality_categories.tsv`](tables/Phase2/step12c_data_quality_categories.tsv)** — master per-sample categorisation (all 433 metadata samples) with full upstream diagnostic columns preserved (Step 12c).

### Lab and phenotyping deliverables

-   **[`tables/Phase2/step12c_samples_redo.csv`](tables/Phase2/step12c_samples_redo.csv)** — single merged lab follow-up sheet. `Lab_action` discriminates **Re-PCR** (existing DNA OK, re-amplify) from **Re-DNA-extraction** (DNA contaminated, re-isolate tissue); rows sorted by `Lab_action → EO → Sample_ID` for field-collection grouping (Step 12c).
-   **[`tables/Phase2/step12c_SI_escape_candidates.csv`](tables/Phase2/step12c_SI_escape_candidates.csv)** — phenotyping deliverable: SI-escape candidates with the recommended controlled-selfing test instruction (Step 12c).

## Output Files

### Key Outputs

-   `tables/Phase2/step9_functional_proteins.fasta` - Validated SRK protein sequences
-   `tables/Phase2/step11_individual_allele_genotypes.tsv` - Allele copy-count matrix (individuals × allele bins; values = number of distinct proteins per individual per allele bin)
-   `tables/Phase2/step11_individual_allele_table.tsv` - Long-format allele assignment table (Individual, Protein, Allele, Count)
-   `tables/Phase2/step12_individual_zygosity.tsv` - Zygosity classifications with columns: `N_distinct_alleles`, `N_total_proteins`, `Zygosity`, `Genotype` (AAAA/AAAB/AABB/AABC/ABCD), `Allele_composition`
-   `tables/Phase2/step9_self_compatible_candidates.txt` - Potentially self-compatible individuals

-   `tables/Phase3/step17_EO_allele_richness.tsv` - Per-EO + per-BL evenness J, rarefied k (N=30, 1000 perms), P_compat at L ∈ {0, 0.10, 0.25, 0.50}, L̂_from_AAAA, geno_mix (AAAA:AAAB:AABB:AABC:ABCD), raw amplicon copy recovery

### Population Genetic Outputs (Phase 3)

#### Step 13 — EO + Bottleneck Lineage integration

-   `tables/Phase3/step13_individual_BL_assignments.tsv` - Per-individual EO, Group, BL (BL1–BL5), Drift_index, BL_status (Assigned / Inferred / Unassigned). Join key consumed by every Phase 3/4 R script that adds a BL stratification.

#### Step 13b — Sampling overview map

-   `figures/Phase3/step13b_sampling_map.png` / `SRK_sampling_map.pdf` - Geographic map of all catalogued *L. papilliferum* populations against their parent BL. Sampled populations drawn as filled triangles; not-yet-sampled populations as filled circles; both coloured by parent BL and sized by habitat area. The 6 focus EOs (EO18, EO25, EO27, EO67, EO70, EO76) carry bold white-fill labels with sample counts; other sampled EOs carry smaller plain labels. Headline: **359 BL-resolved samples across 17 populations in 16 EOs** + **8 JAR samples / 8 germplasm sub-codes** (not plotted) = 367 total successful samples.
-   `figures/Phase3/step13b_sampling_map_presentation.png` / `SRK_sampling_map_presentation.pdf` - Same plot with title/subtitle/caption stripped for slide decks.

#### Step 14 — Population genetics

-   `tables/Phase3/step14_population_genetic_summary.tsv` - EO-level summary with `BL` and `Drift_index` columns; rows ordered by parent BL (total habitat area, then within-BL connectivity → BL4, BL5, BL3, BL1, BL2) then by EO mean Drift_index
-   `tables/Phase3/step14_population_genetic_summary_BL.tsv` - BL-level aggregate (5 rows) with `N_EOs`, `Mean_drift`, plus per-group metrics
-   `figures/Phase3/step14_population_genetic_summary.pdf` - 2-page PDF: page 1 EO bars in connectivity order with the locked Set1 BL palette; page 2 BL aggregate bars

#### Step 15 — Allele accumulation

-   `figures/Phase3/step15_allele_accumulation_curves.pdf` - Species + per-EO + per-BL pages with MM/Chao1/iNEXT asymptote reference lines
-   `figures/Phase3/step15_allele_accumulation_species.png` - Two-panel species figure: stacked bar (observed + predicted-undetected = MM) on the left, rarefaction curve with MM/Chao1 asymptote reference lines on the right; both panels share the bar's y-axis. Same colour key (dark blue = observed, light blue = predicted undetected) as the BL/EO drift-erosion bars. Current dataset 2026-05-11: 49 observed + 10 predicted = MM 59, Chao1 = 61, consensus = 60
-   `figures/Phase3/presentation/step15_allele_accumulation_species_observed.png` / `_predicted.png` - **NEW**. Slide build-up frames: frame 1 shows the bar with only the observed segment (no curve, no MM line); frame 2 shows the full bar + curve + MM dashed line (no Chao1/iNEXT). Same canvas size across both frames for clean animation
-   `figures/Phase3/step15_allele_accumulation_combined.png` - All qualifying EO curves on shared axes, **colored by parent BL**
-   `figures/Phase3/step15_allele_accumulation_BL_combined.png` - **NEW**. The 5 BL aggregate curves on shared axes
-   `figures/Phase3/step15_allele_accumulation_drift_erosion.png` - EO stacked bars sorted by BL with BL color strip on the x-axis baseline
-   `figures/Phase3/step15_allele_accumulation_BL_drift_erosion.png` - **NEW**. BL aggregate stacked bars (BL4 = 32% lost; BL1/BL2 ≥ 86% lost)
-   `tables/Phase3/step15_allele_accumulation_stats.tsv` - Per-level curve statistics (Species + 6 EOs + 5 BLs), MM/Chao1/iNEXT estimates, sampling adequacy targets
-   `tables/Phase3/step15_species_richness_estimates.tsv` - Consensus species allele richness estimate (input to Step 16). Computed from all 367 ingroup individuals.

#### Step 16 — Allele frequency χ²

-   `tables/Phase3/step16_chisq_species_population.tsv` - χ² test statistics at three levels (Species, EO, BL); new `BL` column; rows ordered Species → EOs (sorted by BL) → BLs
-   `figures/Phase3/step16_chisq_species_population_frequency_plots.pdf` - 12 pages (1 species + 6 EOs sorted by BL + 5 BLs); EO/BL pages colored by parent BL palette

#### Step 17 — TP1 mating-pool functionality

-   `tables/Phase3/step17_EO_allele_richness.tsv` - Per-EO (6 rows, N ≥ 15) + per-BL (5 rows) metric table: `level`, `group`, `BL`, `N`, `geno_mix`, `prop_AAAA`, `L_hat_from_AAAA`, `raw_copy_recovery`, `k_observed`, `k_rarefied30_mean/sd`, `shannon_H`, `evenness_J`, `chi2_uniform_p`, `P_compat_uniform_4x`, `P_compat_L0/L0.10/L0.25/L0.50`. Inferred tetraploid genotypes from `SRK_zygosity_from_genotype.R` (every individual filled to 4 S-allele copies). Exact P_compat formula assumes tetrasomic inheritance.
-   `figures/Phase3/step17_P_compat_traffic_light_EO.{png,pdf}` - Stakeholder-facing traffic-light figure: six focal EOs ranked worst-to-best by strict-SI P_compat, coloured red (< 0.20, failed) / amber (0.20–0.40, struggling) / green (≥ 0.40, sustainable), with bootstrap 95 % CIs as horizontal error bars and BL colour squares on the left edge. Headline result drops out at a glance: EO70 (BL2) is the only RED population. Produced by `SRK_P_compat_traffic_light.R`.
-   `figures/Phase3/step17_P_compat_traffic_light_EO_blank.{png,pdf}` - Same framework with no data drawn — for predictions slide; reveal the full figure after explaining the threshold logic.
-   `figures/Phase3/step17_depletion_ranking_{observed,predicted}_{all,EOs,blank}.{png,pdf}` - P_compat × Depletion Index conservation-ranking figures with quadrants labelled HEALTHY / INFORMED BREEDING (frequency skew) / INFORMED BREEDING + ALLELE INJECTION (preventive / urgent). `_observed` uses `k_rarefied30` (conservative); `_predicted` uses the MM asymptote. `_all` plots both BLs and EOs; `_EOs` is EO-only; `_blank` has zones only. Produced by `SRK_depletion_ranking.R`.
-   Legacy `SRK_TP1_tipping_point.R` archived under `archive/`.

#### Step 18 — Allele composition / sharing

-   `figures/Phase3/step18_allele_upset_EOs.{pdf,png}` - EO UpSet (sorted by parent BL; bars/dots colored by parent BL)
-   `figures/Phase3/step18_allele_sharing_heatmap_EOs.{pdf,png}` - Pairwise allele sharing heatmap between EOs
-   `figures/Phase3/step18_allele_upset_BLs.{pdf,png}` - **NEW**. BL-level UpSet — the direct test of independent bottlenecks
-   `figures/Phase3/step18_allele_sharing_heatmap_BLs.{pdf,png}` - **NEW**. 5×5 BL pairwise sharing heatmap (BL4↔BL5 share 12 alleles, the highest off-diagonal value)
-   `figures/Phase3/step18_allele_eulerr_BLs.{pdf,png}` - **NEW (2026-05-22)**. Area-proportional 5-set Euler diagram of BL allele sets — same data as the BL UpSet, rendered as familiar Venn-style ellipses for stakeholder presentations. Produced by the R companion script `SRK_allele_eulerr_BLs.R` (regenerate together with the Python UpSet whenever new data are added).

PNGs land in `figures/`; PDFs in the project root. **Headline finding (2026-05-11 refresh):** 28 of 49 alleles (53 %) are private to a single BL; only Allele_043 and Allele_046 are shared across all 5 BLs.

#### Steps 19–20 — Individual GFS + TP2 tipping point

-   `tables/Phase3/step19_individual_GFS.tsv` - Per-individual Genotypic Fitness Score, genotype class (AAAA/AAAB/AABB/AABC/ABCD), EO, **BL, Drift_index** (Step 19)
-   `tables/Phase3/step20_EO_GFS_summary.tsv` - EO-level mean GFS, genotype class proportions, TP2 status, **BL column** (Step 20)
-   `tables/Phase3/step20_BL_GFS_summary.tsv` - **NEW**. BL-level mean GFS, class proportions, TP2 status (5 rows; all CRITICAL)
-   `figures/Phase3/step19_GFS_plots.pdf` - Multi-page PDF: composition (proportional + counts), per-individual jitter, TP2 scatter — EOs sorted by parent BL with BL-coloured x-axis labels (Steps 19–20)
-   `figures/Phase3/step20_TP2_tipping_point.png` - TP2 scatter; **EOs as circles, BLs as triangles, both Set1-coloured by parent BL**. All 5 BLs and 5/6 plotted EOs CRITICAL (EO67 only AT RISK)

#### Step 21 — Reproductive effort support

-   `figures/Phase3/step21_GFS_reproductive_effort.pdf` - 2-page PDF: EO panel + BL panel — GFS tier composition with reproductive effort annotations
-   `figures/Phase3/step21_GFS_AAAA_allele_composition.pdf` - 2-page PDF: EO panel + BL panel — allele identity of AAAA individuals; Synonymy group 1 alleles highlighted
-   `tables/Phase3/step21_BL_reproductive_effort_summary.tsv` - **NEW**. BL-level reproductive support summary
-   `figures/Phase3/step21_GFS_reproductive_effort_{EO,BL}.png` - Per-EO and per-BL proportional bars
-   `figures/Phase3/step21_GFS_AAAA_allele_composition_{EO,BL}.png` - Per-EO and per-BL AAAA allele identity. **Headline:** Allele_043 + Allele_046 are pan-BL across all 5 BLs and pan-EO in 5/6 focus EOs — confirms shared Synonymy group 1 fixation despite independent bottlenecks

### Phase 4 Outputs (Testing S-allele Hypotheses)

-   `figures/Phase5/step22a_variability_landscape.pdf` / `figures/Phase5/step22a_variability_landscape.png` - **Multi-panel cross-Brassicaceae figure**: per-species smoothed Shannon entropy (LEPA + Brassica + Arabidopsis), per-species HV-region tracks, and Ma 2016 SCR9-contact residue markers (Step 22a)
-   `tables/Phase5/step22a_HV_regions_per_species.tsv` - LEPA / Brassica / Arabidopsis HV runs (start–end, length, threshold) (Step 22a)
-   `tables/Phase5/step22a_LEPA_HV_positions.tsv` - **66 canonical LEPA HV columns** (consumed by Step 22b) (Step 22a)
-   `SRK_HV_overlap_permutation.tsv` - Three pairwise HV-overlap permutation tests (10 000 perms). In the current post-QC LEPA dataset the Brassica↔Arabidopsis HV overlap remains highly significant, but the two LEPA pairs are no longer significant — interpreted as drift-eroded standing variation at LEPA HV columns (Step 22a)
-   `SRK_brassica_hv_mapping.tsv` - 12 Ma 2016 SCR9-contact residues mapped to LEPA columns (helper output)
-   `SRK_HV_allele_distances.tsv` - 58×58 pairwise distance matrix computed on the 66 canonical HV columns (Step 22b)
-   `SRK_functional_allele_groups.tsv` - Allele bin → phylogenetic class assignment, AAAA count, cross power (Step 22b)
-   `SRK_synonymy_candidates.tsv` - All within-class allele pairs with HV distance and testability flag (Step 22b)
-   `SRK_synonymy_groups.csv` - Per-allele synonymy group membership (8 synonymy groups + 19 isolated → 27 effective bins) (Step 22b)
-   `figures/Phase5/step22b_allele_similarity_heatmap.{pdf,png}` - 49×49 HV similarity heatmap ordered by UPGMA with class strips (Step 22b)
-   `SRK_AAAA_cross_design_HV.tsv` - All AAAA × AAAA pairs ranked by category (Incompatible / Synonymy_test / Compatible_within / Compatible_cross) with HV distance and expected outcome (Step 22b)
-   `figures/Phase5/step22b_cross_design_summary.{pdf,png}` - Three-panel figure: HV distance distribution, cross category schematic, Synonymy_test cross interpretation (Step 22b)
-   `figures/Phase5/step22b_HV_cluster_figure.{pdf,png}` - UPGMA dendrogram (HV distances) coloured by class + AAAA availability bar chart (Step 22b)
-   `SRK_synonymy_network_groups.{pdf,png}` / `SRK_synonymy_network_tests.{pdf,png}` - Synonymy network: 9 HV-identical synonymy groups + N-connectivity condensed graph (Step 22b)
-   `SRK_LEPA_synonymy_group_representatives.tsv` - Synonymy group → representative allele mapping (kept vs redundant) (Step 22c)
-   `SRK_wgroup_collapse_entropy_summary.tsv` - Per-species entropy before/after Synonymy group collapse (Step 22c)
-   `SRK_variability_landscape_wgroup_collapsed.{pdf,png}` - Side-by-side LEPA full / LEPA collapsed / Brassica / Arabidopsis entropy comparison (Step 22c)
-   `SRK_perBL_HV_residue_table.tsv` - Per-(HV col × BL) dominant residue + frequency + Shannon entropy among AAAA individuals (Step 22d)
-   `SRK_perBL_HV_concordance_summary.tsv` / `..._per_hv_col.tsv` - Within-LEPA concordance + cross-genera (Brassica, Arabidopsis) dominant-residue match per HV col (Step 22d)
-   `SRK_perBL_entropy_figure.{pdf,png}` - Per-BL entropy heatmap + dominant-residue comparison across BL1–BL5 + Brassica + Arabidopsis at LEPA HV cols (Step 22d)
-   `SRK_cross_plan_H0_SI_validation.tsv` / `..._H1a_within_class_baseline.tsv` / `..._H1b_between_class_baseline.tsv` / `..._H2_synonymy_tests.tsv` / `..._H3_hidden_bin_tests.tsv` - Phased hypothesis-testing cross plan with mother / father IDs, predicted compatibility, replicate counts, and decision rules (Step 22e)
-   `SRK_cross_plan_summary.tsv` / `SRK_cross_plan_summary.{pdf,png}` - Per-phase cross counts + decision-tree figure (Step 22e)
-   `SRK_cross_result_analysis_HV.pdf` - Seed yield distributions and success rates by cross category, with Kruskal-Wallis and Mann-Whitney U tests (Step 23; requires cross data)

### Phase 5 Outputs (Per-Individual SI System Status — Step 25)

> Answers the question: *"What is the status of the SI system?"* at the individual, BL, and EO levels.

-   `tables/Phase4/step25b_individual_SI_status.tsv` — One row per ingroup individual: `n_haps_OK`, `n_haps_REMOVED`, `frac_nonfunctional`, `copies_functional`, `copies_nonfunctional`, `SI_status` (SI / pSI / SC / Insufficient_data), `pSI_confidence` (high / low). Re-aggregates per-haplotype OK/REMOVED calls from the Step-7 `_frame1_stopcodon_log.tsv` files **without** invalidating the canonical 49-allele catalogue (Step 25)
-   `figures/Phase4/step25b_SI_status_species_{full,robust}.png` — Species-level stacked bar. `_full` = all 497 individuals (7-tier inc. pSI_low + Insufficient_data); `_robust` = 263 with defensible calls (5-tier SI/pSI_1nf/pSI_2nf/pSI_3nf/SC) (Step 25)
-   `figures/Phase4/step25b_SI_status_by_BL_{full,robust}.png` — Per-BL stacked %, BLs in `BL_ORDER`. Same full/robust split (Step 25)
-   `figures/Phase4/step25b_SI_status_by_EO_{full,robust}.png` — Per-EO stacked % faceted by BL, within-BL order by ascending drift index. Same full/robust split (Step 25)

**Rule:** Step 25a aligns every Step-7 REMOVED haplotype to the 49 functional reps and flags chimeric Canu artefacts (`AA_distance > 0.05` or `n_stops > 10`). Step 25b counts only real broken haplotypes: `copies_nonfunctional = round(4 × n_REMOVED_real / (n_OK + n_REMOVED_real))`. **Snapshot (n = 401 ingroup):** SI = 247 (62 %), pSI = 15 (10 high / 5 low), SC = 1 (Library010_barcode53, EO76), Insufficient_data = 138 (heavily concentrated in EO76, flagged for re-sequencing).

### Phase 5b Outputs (Null-Aware Genotype Rebuild — Step 26)

> Propagates Step 25 SI status into the canonical genotype tables so downstream metrics (Step 14 pop-gen, Step 17 TP1, Steps 19+20 GFS/TP2) reflect the broken-copy reality. Phase 1–2 outputs stay frozen as the functional-only reference.

-   `tables/Phase4/step26_individual_allele_genotypes_with_nulls.tsv` — 234 individuals × 49 functional alleles + `Allele_NULL` + `Genotype_class_flag`. pSI_high gets explicit nulls (`copies_nonfunctional`); SC = 4 × NULL; rows sum to 4. pSI_low is promoted to SI (Canu-noise floor) (Step 26)
-   `tables/Phase4/step26_individual_zygosity_with_nulls.tsv` — Genotype labels like `AB00` (1 functional + 1 functional + 2 null) or `0000` (SC); 12-tier null-aware scheme (Step 26)
-   `tables/Phase4/step26_samples_for_redo.tsv` — 234 Insufficient_data individuals flagged for re-sequencing, with EO / BL / hap counts / priority (Step 26)
-   `SRK_population_genetic_summary_with_nulls.tsv` / `_BL_with_nulls.tsv` / `.pdf` — Null-aware EO and BL pop-gen with new columns: `Frac_nonfunctional_alleles`, `Mean_null_copies`, `Prop_pSI`, `N_SC` (Step 14b)
-   `tables/Phase4/step17b_EO_allele_richness_with_nulls.tsv` — Null-aware TP1 metrics (P_compat with Allele_NULL in the frequency vector); BL5 P_compat at L=0 drops 0.472 → 0.392 (Step 17b)
-   `SRK_individual_GFS_with_nulls.tsv` / `SRK_EO_GFS_summary_with_nulls.tsv` / `SRK_BL_GFS_summary_with_nulls.tsv` — Null-aware GFS (`GFS_func × n_func/4`); all 5 BLs CRITICAL on TP2 (Steps 19b + 20b)
-   `figures/Phase4/step19b_GFS_with_nulls_composition.png` + `figures/Phase4/step20b_TP2_with_nulls_scatter.png` — Stacked-bar genotype tiers per BL + null-aware TP2 scatter (Steps 19b + 20b)
-   `figures/Phase4/step17b_P_compat_traffic_light_with_nulls.{png,pdf}` (and `_blank` for presentations) — Null-aware random-mating viability per EO. Headline shifts vs canonical: EO76 slips from green into amber (0.36), EO27 / EO25 sit right at the 0.40 floor (0.40 / 0.42), EO70 stays red (0.13). Script: `SRK_P_compat_traffic_light_with_nulls.R`
-   `figures/Phase3/step17_depletion_ranking_observed_with_nulls_{all,EOs,blank}.{png,pdf}` / `SRK_depletion_ranking_predicted_with_nulls_*` — Null-aware DI × P_compat ranking (observed and MM-predicted). Every BL and every focal EO now sits in the right half (DI ≥ 0.5, severe depletion). Script: `SRK_depletion_ranking_with_nulls.R`

**Headline finding:** functional-only P_compat overestimated mating-pool size by 0.03–0.08 across BL3–BL5; null-aware GFS shows 60–75 % of individuals in every BL fall to GFS = 0 (AAAA-equivalent or worse). The conclusion that LEPA is reproductively compromised is *strengthened*, not weakened, by the null-aware rebuild.

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

Key findings from the **335 successful ingroup samples** — 325 BL-resolved across **17 populations in 16 Element Occurrences** spanning all five bottleneck lineages BL1–BL5 (with the 6 focus EOs EO18, EO25, EO27, EO67, EO70, EO76 driving most downstream analyses); plus 10 samples in the JAR region with unresolved coordinates; plus **59 metadata samples reserved for lab follow-up** (Libraries 001–010). See `figures/Phase3/step13b_sampling_map.png` for the sampling overview:

- **49 observed S-allele bins** species-wide, with a predicted total of **60** (consensus of Michaelis-Menten = 59 and Chao1 = 61) — the species retains ~82% of its expected SI repertoire, but the retention is unevenly distributed across BLs.
- The five bottleneck lineages have eroded their S-allele pools very unevenly: **BL4 retains ~68%** of its expected richness while **BL1 and BL2 have lost ≥86%**, consistent with their lower within-BL connectivity.
- Allele frequencies are significantly skewed from NFDS expectations at every level (Species, EO, BL; χ² *p* < 10⁻⁷).
- **28 of 49 alleles (53%) are BL-private**; only **Allele_043 + Allele_046** (Synonymy group 1) are pan-BL, and these are the alleles driving most AAAA fixation across the species.
- A majority of individuals carry an **AAAA genotype** (single allele, four copies, GFS = 0), flagging high risk of reduced SI function; a further substantial fraction carry an **AAAB** genotype (GFS = 0.500) — dosage-imbalanced individuals that appear heterozygous but produce fewer diverse gametes than AABB individuals (GFS = 0.667).
- All five BLs and 5/6 focus EOs are flagged **CRITICAL** for Tipping Point 2 (TP2); EO67 (BL4) is the only AT RISK EO.
- A formal **Step 12c Data Quality Evaluation** identifies **5 of 7 Complete_loss SI-escape candidates concentrated in EO76 (BL3)**, framing EO76 as a candidate SI → SC transitional population under the Igić–Lande–Kohn framework.
- A hypothesis-testing **cross plan (Step 22e)** issues 104 phased crosses (H0 SI validation; H1a/b within/between-Class baselines; H2 synonymy-bin tests; H3 hidden-bin tests via heterozygous donors), with every cross traceable to its sequence-based category and its genotype-based feasibility from Step 12.

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
