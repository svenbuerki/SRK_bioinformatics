# SRK Bioinformatics Pipeline — Protocol

> **Purpose:** Step-by-step protocol for running the SRK haplotype reconstruction and population analysis pipeline. For full scientific background, see [index.Rmd](index.Rmd) (rendered at https://svenbuerki.github.io/SRK_bioinformatics/).

---

## Environment

Activate the conda environment before running any bioinformatics scripts:

```bash
conda activate nanopore_pipeline
```

For Python scripts that require BioPython/pandas:

```bash
source srk_env/bin/activate
```

---

## Phase 1: SRK Amplicon Sequence Assembly

### Step 1 — Prepare Canonical Sequences

> **Run ONCE only.** Produces the reference sequences used throughout the pipeline.

**Script:** `reverse_complement_startcodon_fasta.py`

**Command:**
```bash
python reverse_complement_startcodon_fasta.py
# Prompted input: path to SRK_canonical_haplotype_sequences.fasta
```

**Input:**
- `SRK_canonical_haplotype_sequences.fasta` — raw SRK paralogs from genome mining

**Key parameters:**
- Sequences are reverse-complemented (gene is on the negative strand)
- Start codon (ATG) is corrected if missing from the first 1–2 positions

**Outputs:**
- `SRK_canonical_haplotype_sequences_revcomp_ATGfixed.fasta` — corrected reference sequences
- `SRK_canonical_haplotype_sequences_revcomp_ATGfixed_startcodon_fix.log` — log of corrections

**Gene annotation** (run separately, lab Linux computer):
```bash
augustus --species=arabidopsis --strand=both --genemodel=partial --protein=on \
  --codingseq=on --introns=on --start=on --stop=on --cds=on --gff3=on \
  SRK_canonical_haplotype_sequences_revcomp_ATGfixed.fasta \
  > augustus_SRK_canonical_haplotype_sequences_revcomp_ATGfixed_arabidopsis.gff
```

**Output:**
- `augustus_SRK_canonical_haplotype_sequences_revcomp_ATGfixed_arabidopsis.csv` — exon annotations used in Step 6

---

### Step 2 — Nanopore Amplicon Assembly and Phasing

**Script:** `nanopore_assembly_pipeline_barcode_range.sh`

**Command:**
```bash
./nanopore_assembly_pipeline_barcode_range.sh
# Prompted inputs: LibraryID (e.g. Library0001), start barcode, end barcode
```

> **Important:** LibraryID format must be `Library0001` (no underscores).

**Expected directory structure:**
```
Library001/
├── barcode01/
│   └── *.fastq.gz
├── barcode02/
└── ...
```

**Key parameters:**
- 4 CANU replicates (`REPS=4`) to detect low-frequency haplotypes
- `genomeSize=5k`, `correctedErrorRate=0.15`
- 3 rounds of RACON polishing
- FreeBayes variant calling with `-p 4` (tetraploid ploidy)
- WhatsHap polyphase with `--ploidy 4`
- Up to 4 phased haplotypes exported per individual via `bcftools consensus`

**Outputs per barcode:**
- `combined_unique_contigs.fasta` — deduplicated contigs from all 4 CANU runs
- `polished_r3.fasta` — final polished assembly
- `aligned.bam` / `.bai` — read alignments
- `variants.vcf.gz` / `phased.vcf.gz` — variants and phased variants
- `${Library}_${barcode}_hap[1-4].fasta` — individual haplotype FASTA files
- `${Library}_${barcode}_Phased_haplotypes.fasta` — combined haplotype file

---

### Step 3 — Haplotype Orientation and Consolidation

**Script:** `orient_haplotypes_library.sh`

**Command:**
```bash
./orient_haplotypes_library.sh
# Prompted inputs:
#   Library name (e.g. Library001)
#   Canonical reference FASTA (full path)
```

**Inputs:**
- Library folder containing all `barcode*/` subdirectories with `*_Phased_haplotypes.fasta`
- Canonical reference: `SRK_canonical_haplotype_sequences_revcomp_ATGfixed.fasta`

**Key behaviour:**
- Each barcode's contigs are aligned to the reference with minimap2
- Contigs with **no primary alignment** to the reference are discarded (off-target assemblies). A WARNING is printed to the terminal for any barcode where this occurs.
- Contigs mapping to the reverse strand are reverse-complemented
- This filter was added after Library009 produced off-target contigs in barcodes 72, 75, 94, 40, 38, and 28 that stretched the MAFFT alignment from ~4,500 bp to 10,779 bp

**Outputs:**
- `${Library}_${barcode}_Phased_haplotypes.oriented.fasta` — on-target oriented haplotypes per barcode
- `all_${Library}_Phased_haplotypes.fasta` — consolidated oriented haplotypes for the entire library

---

### Step 4 — FASTA Sequence Filtering and Reference Integration

**Script:** `filter_and_add_references.sh`

**Command:**
```bash
./filter_and_add_references.sh
# Prompted inputs:
#   Merged haplotype FASTA (e.g. all_Library001_Phased_haplotypes.fasta)
#   Canonical reference FASTA (full path)
#   Minimum sequence length (bp)
#   Maximum sequence length (bp)
```

**Key parameters:**
- Minimum length: **3250 bp** (for ~4500 bp expected amplicon)
- Maximum length: **4000 bp** — removes mis-assembled or chimeric contigs that stretch the MSA. Library009 (96 samples) produced outlier assemblies of 4134–4841 bp in barcodes 32, 40, and 72 that caused MAFFT to produce an 18,212 bp-wide alignment instead of the expected ~4,500 bp. The max filter eliminates these. Reference sequences are prepended before filtering and are always retained regardless of this threshold.

**Input:**
- `all_Library00X_Phased_haplotypes.fasta` — output of Step 3

**Output:**
- `all_Library00X_Phased_haplotypes_filtered_min3250_max4000.fasta` — length-filtered sequences including canonical references

---

### Step 5 — Multiple Sequence Alignment

**Tool:** MAFFT

**Command:**
```bash
mafft --auto --adjustdirection all_Library00X_Phased_haplotypes_filtered_min3250_max4000.fasta \
  > all_Library00X_Phased_haplotypes_filtered_min3250_max4000_aligned.fasta
```

> Run this command for each library. `--adjustdirection` handles any residual orientation issues.

**Input:** `all_Library00X_Phased_haplotypes_filtered_min3250_max4000.fasta`

**Output:** `all_Library00X_Phased_haplotypes_filtered_min3250_max4000_aligned.fasta`

---

### Step 6 — Exon Extraction from Multiple Sequence Alignments

**Script:** `extract_exons_with_annotation.py`

**Command:**
```bash
python extract_exons_with_annotation.py
# Prompted inputs:
#   MSA FASTA path (e.g. all_Library006_Phased_haplotypes_filtered_min3250_max4000_aligned.fasta)
#   Canonical reference ID (exact header, e.g. SRK_BEA_hybrid_bp_hap1_p_ctg_fa|amp_1|h1tg000019l|1700582|3407)
#   AUGUSTUS CSV annotation file (full path)
#   Strand of annotation for canonical reference (+)
#   Minimum exon length (0)
```

**Inputs:**
- Aligned MSA FASTA from Step 5
- Canonical reference ID (must match exact FASTA header)
- `augustus_SRK_canonical_haplotype_sequences_revcomp_ATGfixed_arabidopsis.csv` — from Step 1
- Strand: `+`
- Min exon length: `0`

**Output:**
- `all_Library00X_Phased_haplotypes_filtered_min3250_max4000_aligned_exons.fasta`

---

### Step 7 — MSA Gap Backfilling and Terminus Processing

**Script:** `backfill_alignment_ends.py`

**Command:**
```bash
python backfill_alignment_ends.py
# Prompted inputs:
#   MSA FASTA file (exon FASTA from Step 6)
#   Canonical reference ID (exact FASTA header)
```

**Key parameter:**
- `WINDOW = 25` bp backfilled from reference at each terminus (hardcoded)
- Leading/trailing terminal gaps beyond the backfill window are replaced with `N`

**Input:** `all_Library00X_Phased_haplotypes_filtered_min3250_max4000_aligned_exons.fasta`

**Output:** `all_Library00X_Phased_haplotypes_filtered_min3250_max4000_aligned_exons_backfilled.fasta`

---

### Step 8 — DNA to Amino Acid Translation (Optional)

> **Note:** This step is optional and applies to individual libraries. It may be skipped if protein-level per-library analysis is not required; protein translation for the full cross-library dataset is performed in Step 9.

**Script:** `translate_filter_align_AA.py`

**Command:**
```bash
python translate_filter_align_AA.py
# Prompted inputs:
#   Aligned DNA FASTA file (backfilled exon FASTA from Step 7)
#   Translation frame (1)
```

**Key parameters:**
- Frame: `1`
- Genetic code: standard nuclear (table 1)
- Only sequences without internal stop codons are retained in the final filtered output

**Input:** `all_Library00X_Phased_haplotypes_filtered_min3250_max4000_aligned_exons_backfilled.fasta`

**Key output:**
- `all_Library00X_Phased_haplotypes_filtered_min3250_max4000_aligned_exons_backfilled_frame1_AA_filtered_aligned.fasta`

---

## Phase 2: Functional Proteins, S Alleles and Genotyping

> Steps 9–12 operate on the **merged cross-library dataset** (all `*_exons_backfilled.fasta` files in the working directory). Run these steps from the directory containing all library outputs.

### Step 9 — SRK Protein Translation, Alignment, and Abundance Filtering

**Script:** `define_functional_proteins.py`

**Command:**
```bash
python define_functional_proteins.py
# Prompted inputs (after overwrite check):
#   Translation frame: 1
#   Minimum sequences per protein (recommended 2–5): 5
#   Maximum functional proteins per individual (e.g. 4 for allotetraploid): 4
#   Use enhanced SRK domain validation? (y/n): n
```

**Inputs (auto-detected by glob):**
- All `*_exons_backfilled.fasta` files in the working directory

**Key parameters:**
- Input glob pattern: `*_exons_backfilled.fasta`
- `frame = 1`
- `min_length = 100` aa (hardcoded)
- `min_count = 5` — abundance filter (minimum times a protein must appear)
- `max_alleles = 4` — individuals with >4 distinct proteins are excluded (tetraploid ploidy)
- Enhanced SRK domain validation: `n`

**Outputs:**
- `SRK_proteins_raw.fasta` — all translated proteins before filtering
- `SRK_proteins_aligned.fasta` — MAFFT-aligned proteins
- `SRK_functional_proteins.fasta` — final abundance-filtered proteins
- `SRK_functional_protein_key.tsv` — mapping: Original_sequence_ID → Protein (overloaded individuals excluded)
- `SRK_individual_status_report.tsv` — per-individual processing summary
- `SRK_self_compatible_candidates.txt` — individuals with no/low functional SRK
- `SRK_too_many_alleles.txt` — individuals excluded for exceeding `max_alleles`

---

### Step 10 — Distance-Based SRK S-Allele Definition

**Scripts:** `define_SRK_alleles_from_distance.py` and `SRK_AA_mutation_heatmap.py`

> Run Step 10a first so allele assignments are available for Step 10b row ordering.

#### Step 10a — Allele clustering

**Command:**
```bash
python define_SRK_alleles_from_distance.py
```

**Key parameters (edit at top of script):**
| Parameter | Default | Notes |
|-----------|---------|-------|
| `INPUT_FASTA` | `SRK_functional_proteins_aligned.fasta` | Aligned proteins from Step 9 |
| `N_ALLELES` | `63` | Option A: fix number of alleles directly |
| `DIST_THRESHOLD` | `0.01` (1%) | Option B: fix p-distance cutoff; use sensitivity plot elbow to calibrate |
| `DOMAIN_REGION` | `(31, 430)` | S-domain columns (1-based); adjust if species-specific annotation available |

> **Calibration:** inspect `SRK_protein_distance_analysis.pdf` sensitivity curve; choose the elbow (typically 1–2% for Nanopore data). Use `N_ALLELES` (Option A) if the elbow is easier to read on the y-axis than on the distance axis. For the current dataset (Libraries 001–009, 308 individuals, 302 functional proteins), the sensitivity curve showed a plateau between 50–100 alleles centred at 63 (implied threshold ≈ 0.005), which was used as `N_ALLELES`.

**Outputs:**
- `SRK_protein_allele_assignments.tsv` — Protein → Allele mapping
- `SRK_protein_allele_representatives.fasta` — one representative sequence per allele
- `SRK_protein_distance_analysis.pdf` — sensitivity curve + pairwise distance heatmap

#### Step 10b — Amino acid variation heatmaps

**Command:**
```bash
python SRK_AA_mutation_heatmap.py
```

**Key parameters (edit at top of script):**
| Parameter | Default |
|-----------|---------|
| `INPUT_FASTA` | `SRK_functional_proteins_aligned.fasta` |
| `ALLELE_TSV` | `SRK_protein_allele_assignments.tsv` |
| `MAX_GAP_FREQ` | `0.20` |
| `MAX_POSITIONS` | `300` |

**Outputs:**
- `SRK_AA_mutation_heatmap.pdf` — AA physicochemical class heatmap (proteins × variable positions)
- `SRK_AA_frequency_heatmap.pdf` — AA frequency heatmap per variable position
- `SRK_AA_variable_positions.tsv` — variable positions with entropy and AA composition

---

### Step 11 — SRK S-Allele Genotyping

**Script:** `genotype_SRK_from_alleles.py`

**Command:**
```bash
python genotype_SRK_from_alleles.py
# Prompted: overwrite check for output files
```

**Key settings (edit at top of script):**
```python
PROTEIN_KEY_FILE  = "SRK_functional_protein_key.tsv"   # from Step 9
ALLELE_TSV_FILE   = "SRK_protein_allele_assignments.tsv"  # from Step 10a
METADATA_FILE     = "sampling_metadata.csv"   # set None to skip ingroup filter
INGROUP_ONLY      = True  # set False to include all Library* individuals
```

**Inputs:**
- `SRK_functional_protein_key.tsv` — from Step 9
- `SRK_protein_allele_assignments.tsv` — from Step 10a
- `sampling_metadata.csv` — must contain `Library`, `barcode`, `Ingroup` columns

**Outputs:**
- `SRK_individual_allele_table.tsv` — long format: Individual | Protein | Allele | Count
- `SRK_individual_allele_genotypes.tsv` — wide count matrix: individuals × allele bins (values = allele copy counts)

---

### Step 12 — SRK Zygosity Analysis (Tetraploid)

**Script:** `SRK_zygosity_from_genotype.R`

**Command:**
```bash
Rscript SRK_zygosity_from_genotype.R
```

**Input:** `SRK_individual_allele_genotypes.tsv` — from Step 11

**Genotype classification logic:**
| Genotype | Criterion |
|----------|-----------|
| AAAA | 1 distinct allele bin |
| AAAB | 2 bins, dominant count ≥ 2× minor count |
| AABB | 2 bins, roughly equal counts |
| AABC | 3 bins |
| ABCD | ≥ 4 bins |

**Outputs:**
- `SRK_individual_zygosity.tsv` — columns: Individual, N_distinct_alleles, N_total_proteins, Zygosity, Genotype, Allele_composition (e.g. `Allele_044(2)+Allele_047(2)`)
- `SRK_zygosity_distribution.pdf` — bar plot of genotype class frequencies

---

### Step 12b — SRK Class Assignment and Dominance Prediction ⚠️ IN PROGRESS

> ⚠️ **Work in Progress — this step may not function correctly.** The script has not been fully validated on this dataset and the automated class assignment requires confirmation against published SRK references with known class assignments. Do not use outputs for conservation management decisions without independent validation. See also `tmp_SRK_classes_prediction.Rmd` (working draft, possibly outdated).

**Biological rationale:** SRK alleles in Brassicaceae are phylogenetically structured into two classes with a strict dominance hierarchy — Class I alleles are dominant over Class II in heterozygous individuals. This directly determines cross-compatibility: individuals sharing a Class I allele are incompatible regardless of Class II alleles. Knowing which class each individual carries is therefore essential for planning crosses in a conservation breeding program.

**Script:** `srk_phylogenetic_analysis.py`

**Command:**
```bash
python srk_phylogenetic_analysis.py
# IQ-TREE3 executable auto-detected; prompted for path if not found
```

**Inputs:**
- `SRK_functional_proteins_aligned.fasta` — from Step 9
- `SRK_individual_genotypes.tsv` — ⚠️ **Note:** the script expects a binary presence/absence matrix; verify compatibility with the current count-matrix format (`SRK_individual_allele_genotypes.tsv`) before running
- `sampling_metadata.csv` — must contain `Library`, `barcode`, `Ingroup` columns

**Key parameters:**
- IQ-TREE model set: `LG,WAG,JTT`; UFBoot + SH-aLRT: 1000 replicates
- Class split: longest internal branch (automated — validate in FigTree/iTOL)

**Outputs:**
- `SRK_ingroup_proteins_aligned.fasta` — filtered ingroup alignment
- `SRK_phylogeny.treefile` — ML tree (unrooted)
- `SRK_phylogeny_midpoint_rooted.treefile` — for visual validation in FigTree/iTOL
- `SRK_phylogeny.iqtree` — full IQ-TREE log
- `SRK_protein_classes.tsv` — Protein → Class I / Class II
- `SRK_individual_class_genotype.tsv` — per-individual dominance prediction (Class_I_dominant / Class_II_only_codominant / no_functional_SRK)

---

## Phase 3: Data Analyses

> All Phase 3 scripts read from the shared outputs of Phase 2. Run from the same working directory.

Phase 3 evaluates the evolutionary status and conservation implications of the SI system by integrating population genetic theory with the genotype data generated in Phases 1 and 2. Three interconnected scientific questions structure the analyses.

**Q1 — What are the population genetic parameters of the SRK system?**
Before interpreting SI diversity in an evolutionary framework, it is necessary to characterise each population's fundamental genetic state: allele richness, effective allele number, heterozygosity, and mean alleles per individual (Step 13). These descriptors provide the empirical foundation for all downstream comparisons.

**Q2 — What is the total allele richness of the SI system at the species level?**
Under negative frequency-dependent selection (NFDS), rare S-alleles confer a reproductive advantage because pollen bearing a rare allele can fertilise a larger proportion of compatible stigmas. This self-reinforcing dynamic drives the accumulation of allelic diversity and, at evolutionary equilibrium, pushes allele frequencies towards equality. The *species allele richness* — estimated from accumulation curves and asymptotic richness estimators (Michaelis-Menten, Chao1, iNEXT) applied to the full dataset (Step 14) — approximates this NFDS equilibrium. It represents the total pool of functionally distinct S-alleles maintained across the species and serves as the reference optimum against which population-level deficits are measured. This question is foundational: without a species-level baseline, it is impossible to determine how much allelic diversity individual populations have lost relative to the evolutionary expectation.

**Q3 — How has genetic drift eroded the SI system at the population level?**
In small or isolated populations, genetic drift counteracts NFDS by stochastically eliminating rare alleles, skewing allele frequencies away from the equal-frequency expectation, and increasing the prevalence of homozygous genotypes. Drift erodes the SI system at two hierarchical levels, each corresponding to a conservation tipping point. **Tipping Point 1 (TP1)** is breached when a population has lost so many S-alleles relative to the species optimum that inter-population allele transfers are required to restore richness. **Tipping Point 2 (TP2)** is breached when the distribution of the remaining alleles across individuals has degraded to the point that managed crossing within the population cannot recover reproductive fitness without external introductions. Phase 3 quantifies both levels through four complementary analyses:

*TP1 — Allele richness and frequency analyses (Steps 14–17):*

- **Allele accumulation curves** (Step 14) reveal how many S-alleles each population retains relative to the species optimum, measuring the depth of drift-driven allele loss.
- **Allele frequency analysis** (Step 15) tests whether population-level frequencies deviate from the equal-frequency NFDS expectation via χ² goodness-of-fit, using the species richness estimate from Step 14 as the optimum.
- **TP1 tipping point analysis** (Step 16) synthesises Steps 13–15 into a single diagnostic plot, positioning each EO by the proportion of the species optimum it retains (x axis) and by allele frequency evenness — the ratio of effective allele number to observed allele count (Ne/N) — as a measure of departure from NFDS equal-frequency expectations (y axis). EOs breaching both thresholds (< 50% of optimum and Ne/N < 0.80) are flagged CRITICAL.
- **Allele composition comparison** (Step 17) identifies which alleles are private to single populations versus shared across Element Occurrences, distinguishing populations that retain unique allelic diversity from those impoverished by isolation or founder effects.

*TP2 — Genotypic fitness analysis (Steps 18–20):*

- **Genotypic Fitness Score** (Step 18) translates allele-level diversity into individual reproductive fitness. In a tetraploid, each diploid gamete is formed by randomly sampling 2 of the 4 allele copies at the SRK locus, yielding C(4,2) = 6 equally probable gamete combinations. GFS is the proportion of those combinations that carry two distinct alleles. This reveals dosage asymmetries invisible to zygosity classification alone: an AABB individual (2+2 copies) produces 4 of 6 heterozygous gametes (GFS = 0.667), while an AAAB individual (3+1 copies) produces only 3 of 6 (GFS = 0.500), even though both carry exactly two distinct alleles. In self-incompatible plants, an individual's value as a breeding partner depends not only on which key/lock types it carries, but also on the number of distinct alleles present across its genome copies. Individuals with more key/lock types can participate in more compatible crosses, making them especially valuable in a managed breeding program.
- **TP2 Tipping Point Analysis** (Step 19) places individual-level GFS in interaction with the population-level proportion of AAAA genotypes. EOs breaching both criteria (mean GFS < 0.667 and > 30% AAAA) are flagged **CRITICAL**; those breaching either one are flagged **AT RISK**; the remainder are **OK**.
- **Reproductive effort support** (Step 20) visualises, per EO, the proportion of individuals at each GFS tier, quantifying the fraction of the population capable of contributing allelic diversity to compatible crosses.

Together, these analyses trace the consequences of demographic decline from species-level diversity baselines, through population-level allele erosion and frequency skew, to individual-level reproductive fitness — providing a quantitative framework for prioritising conservation interventions.

---

### Step 13 — Population Genetics Statistics

**Script:** `SRK_population_genetic_summary.R`

**Command:**
```bash
Rscript SRK_population_genetic_summary.R
```

**Inputs:**
- `SRK_individual_allele_genotypes.tsv` — from Step 11
- `SRK_individual_zygosity.tsv` — from Step 12
- `sampling_metadata.csv` — must contain `SampleID`, `Pop`, `Ingroup` columns

**Key metrics computed per population:**
- N_alleles — allele richness
- Effective_alleles_Ne = 1/Σpᵢ² (based on `colSums` copy counts)
- Prop_heterozygous / Prop_homozygous
- Mean_alleles per individual

**Outputs:**
- `SRK_population_genetic_summary.tsv`
- `SRK_population_genetic_summary.pdf` — 4-panel plot (Ne, heterozygosity, allele richness, N)

---

### Step 14 — Allele Accumulation Curves and Richness Estimation

> **Run before Step 15.** Produces the species richness estimate used as the NFDS optimum in Step 15.

**Script:** `SRK_allele_accumulation_analysis.R`

**Command:**
```bash
Rscript SRK_allele_accumulation_analysis.R
```

**Inputs:**
- `SRK_individual_allele_genotypes.tsv` — from Step 11
- `sampling_metadata.csv` — must contain `SampleID`, `Pop`, `Ingroup` columns

**Key parameters:**
- Permutations: 1000 (default)
- Populations with <5 individuals excluded from population-level curves
- Three richness estimators: Michaelis-Menten (MM), Chao1, iNEXT (optional — `install.packages("iNEXT")`)

**Outputs:**
- `SRK_allele_accumulation_curves.pdf` — species- and population-level curves with asymptote reference lines
- `SRK_allele_accumulation_combined.png` — single-panel combined plot showing all EO accumulation curves on the same axes; end-of-curve labels show observed allele count and MM-predicted total (e.g. `EO27 (15/26)`); legend entries include `obs=` and `MM=` values
- `SRK_allele_accumulation_drift_erosion.png` — stacked bar chart decomposing the allele deficit per EO into observed alleles, alleles predicted-undetected (EO MM − observed), and alleles lost to genetic drift (species MM − EO MM); quantifies the irreversible erosion of SI diversity driven by drift (60–89% of the 65-allele species optimum lost per EO)
- `SRK_allele_accumulation_stats.tsv` — curve statistics per level including MM, Chao1, iNEXT estimates and sampling adequacy targets
- `SRK_species_richness_estimates.tsv` — **consensus species allele richness; required as input for Step 15**

---

### Step 15 — Allele Frequency Analysis

> Requires `SRK_species_richness_estimates.tsv` from Step 14.

**Script:** `SRK_chisq_species_population.R`

**Command:**
```bash
Rscript SRK_chisq_species_population.R
```

**Inputs:**
- `SRK_species_richness_estimates.tsv` — from Step 14 (MM estimate used as NFDS optimum)
- `SRK_individual_allele_genotypes.tsv` — from Step 11
- `sampling_metadata.csv`

**Key analysis:**
- χ² goodness-of-fit test vs. equal-frequency NFDS expectation at species and population levels
- Populations with <5 individuals excluded
- All frequency plots share x-axis scaled to species optimum (MM estimate); missing-allele zone shown

**Outputs:**
- `SRK_chisq_species_population.tsv` — χ² statistics (X2, df, p-value) per level
- `SRK_chisq_species_population_frequency_plots.pdf` — ranked allele frequency plots with NFDS expectation, LOESS trend, and missing-allele zone

---

### Step 16 — TP1 Tipping Point Analysis

> Requires outputs from Steps 13 and 14. Run after both have completed.

**Script:** `SRK_TP1_tipping_point.R`

**Command:**
```bash
Rscript SRK_TP1_tipping_point.R
```

**Inputs:**
- `SRK_population_genetic_summary.tsv` — from Step 13
- `SRK_species_richness_estimates.tsv` — from Step 14

**Key metrics per EO:**
- `prop_optimum` = N_alleles / MM_estimate — proportion of the species-level S-allele pool retained, using the Michaelis-Menten asymptote as the species optimum.
- `evenness` = Ne / N_alleles — frequency evenness index. Ne (effective allele number) = 1/Σpᵢ², where pᵢ is the frequency of allele i; it answers *how many equally frequent alleles would produce the same level of diversity as observed*. When all alleles are at equal frequency (NFDS ideal), Ne = N and the ratio = 1.0. Genetic drift pushes the ratio downward by creating dominant alleles and marginalising rare ones. An EO with evenness = 0.50, for example, has allele frequencies so uneven that only half the observed alleles are effectively contributing to SI function — the remainder are present in such low copy numbers that they are rarely expressed in crosses and are at elevated risk of loss through drift.

**TP1 thresholds (adjustable at top of script):**

| Parameter | Default | Criterion |
|-----------|---------|-----------|
| `TP1_PROP_OPTIMUM` | 0.50 | EO retains < 50% of species optimum |
| `TP1_EVENNESS` | 0.80 | Ne/N_alleles below 0.80 |

EOs breaching both criteria are flagged **CRITICAL**; those breaching either one are flagged **AT RISK**; the remainder are **OK**.

**Outputs:**
- `SRK_TP1_summary.tsv` — per-EO N_alleles, prop_optimum, evenness, and TP1 status
- `SRK_TP1_tipping_point.png` — scatter plot positioning each EO by richness retained (x) and frequency evenness (y), with colour-coded zone polygons for all three categories (OK, AT RISK, CRITICAL) and threshold lines
- `SRK_TP1_tipping_point_blank.png` — same plot without data points, for presentation use

---

### Step 17 — Allele Composition Comparison Across Element Occurrences

**Script:** `SRK_allele_sharing_EOs.py`

**Command:**
```bash
# Default: reads from current directory, compares EOs 25 27 67 70 76
python SRK_allele_sharing_EOs.py

# Custom paths / EO selection
python SRK_allele_sharing_EOs.py \
    --genotypes SRK_individual_allele_genotypes.tsv \
    --metadata sampling_metadata.csv \
    --pops 25 27 67 70 76 \
    --min-samples 5 \
    --outdir figures/
```

**Inputs:**
- `SRK_individual_allele_genotypes.tsv` — from Step 11
- `sampling_metadata.csv` — must contain `Pop` and `Ingroup` columns

**Key parameters:**

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--pops` | `25 27 67 70 76` | EO IDs to compare (space-separated; matched as strings) |
| `--min-samples` | `5` | Minimum individuals per EO |
| `--outdir` | `.` | Output directory |

**Key analysis:**
- For each non-empty intersection of EO sets, counts alleles exclusive to exactly that combination
- UpSet plot sorted by intersection size; private alleles colour-coded per EO
- Pairwise heatmap shows raw shared-allele counts for every EO pair (diagonal = EO allele richness)
- No additional dependencies beyond the `polyploid-model` environment (pandas, numpy, matplotlib)

**Outputs:**
- `SRK_allele_upset_EOs.pdf` — UpSet plot of all allele set intersections
- `SRK_allele_sharing_heatmap_EOs.pdf` — pairwise sharing heatmap

---

### Step 18 — Individual Genotypic Fitness Score (GFS)

**Script:** `SRK_individual_GFS.R`

**Command:**
```bash
Rscript SRK_individual_GFS.R
```

**Inputs:**

| File | Description |
|------|-------------|
| `SRK_individual_zygosity.tsv` | Per-individual tetraploid genotype pattern from Step 12 (`AAAA`/`AAAB`/`AABB`/`AABC`/`ABCD`) |
| `sampling_metadata.csv` | Must contain `SampleID`, `Pop`, and `Ingroup` columns |

**Key analysis:**

Computes the **Genotypic Fitness Score** — the proportion of heterozygous diploid gametes that a tetraploid individual can produce. GFS is derived directly from the genotype pattern called in Step 12, avoiding dependence on partial haplotype recovery (many individuals have fewer than 4 phased proteins):

| Genotype | GFS |
|----------|-----|
| ABCD | 1.000 |
| AABC | 0.833 |
| AABB | 0.667 |
| AAAB | 0.500 |
| AAAA | 0.000 |

This metric differentiates genotype classes that zygosity analysis treats as equivalent (e.g., AABB = 0.667 vs AAAB = 0.500).

**Key parameters (top of script):**

| Parameter | Default | Notes |
|-----------|---------|-------|
| `EO_MAP` | 5-EO lookup | Maps sub-population codes to EO identifiers; extend for new sub-populations |
| `TP2_MEAN_GFS` | `0.667` | Adjustable threshold for mean GFS (used in Step 19) |
| `TP2_PROP_AAAA` | `0.30` | Adjustable threshold for proportion AAAA (used in Step 19) |

**Outputs:**

| File | Content |
|------|---------|
| `SRK_individual_GFS.tsv` | Per-individual GFS, genotype class, EO, and method |

---

### Step 19 — TP2 Tipping Point Analysis

> Requires `SRK_individual_GFS.tsv` from Step 18. Performed by the same script (`SRK_individual_GFS.R`); TP2 flags are computed after per-individual GFS values are aggregated to the EO level.

**Script:** `SRK_individual_GFS.R` (same run as Step 18)

**Key analysis:**

Places individual-level GFS in interaction with population-level proportion of AAAA genotypes to evaluate whether an EO's allele distribution has degraded beyond the reach of managed crossing. EO-level summaries are evaluated against two thresholds:

| Threshold | Default | Criterion |
|-----------|---------|-----------|
| `TP2_MEAN_GFS` | 0.667 | Mean GFS below AABB level |
| `TP2_PROP_AAAA` | 0.30 | >30% AAAA individuals |

EOs breaching both thresholds simultaneously are flagged **CRITICAL**; those breaching either one are flagged **AT RISK**; the remainder are **OK**.

**Outputs:**

| File | Content |
|------|---------|
| `SRK_EO_GFS_summary.tsv` | EO-level GFS statistics, class proportions, TP2 status |
| `SRK_GFS_plots.pdf` | Four diagnostic plots: stacked bars (proportional + count), individual jitter, TP2 tipping point map |
| `SRK_GFS_plots_p3_TP2_tipping_point.png` | Standalone TP2 tipping point map with colour-coded zone polygons for all three categories (OK, AT RISK, CRITICAL) |
| `SRK_GFS_plots_p3_TP2_tipping_point_blank.png` | Same plot without data points, for presentation use |

---

### Step 20 — Reproductive Effort Support per Element Occurrence

> Requires `SRK_individual_GFS.tsv` from Step 18.

**Script:** `SRK_GFS_reproductive_effort.R`

**Command:**
```bash
Rscript SRK_GFS_reproductive_effort.R
```

**Input:** `SRK_individual_GFS.tsv` — per-individual GFS scores from Step 18

**Key analysis:**

Visualises, for each of the five focus EOs, the proportion of individuals at each GFS tier as a horizontal proportional bar chart. In a self-incompatible plant, an individual's value as a breeding partner depends not only on which key/lock types it carries, but also on the number of distinct alleles present across its genome copies. Individuals with more key/lock types can participate in more compatible crosses, making them especially valuable in a managed breeding program. The figure makes this explicit by distinguishing the fraction of each population that can contribute allelic diversity to gametes (GFS > 0) from those that cannot (AAAA, GFS = 0), and annotates each EO with the proportion of supporting individuals and mean GFS.

EOs are ordered by mean GFS (ascending); the TP2 AAAA threshold (30%) is marked as a dashed vertical line. Annotations to the right of each bar show: proportion supporting (%), count (n / N), and mean GFS.

**Key parameters (top of script):**

| Parameter | Default | Notes |
|-----------|---------|-------|
| `EO_FOCUS` | 5-EO vector | Element Occurrences to include |
| `TP2_PROP_AAAA` | `0.30` | Threshold line position (matches TP2 threshold from Step 19) |

**Outputs:**

| File | Content |
|------|---------|
| `SRK_GFS_reproductive_effort.pdf` | Horizontal proportional bar chart — GFS tier composition per EO with reproductive effort annotations |
| `SRK_GFS_reproductive_effort.png` | Same figure as PNG at 200 dpi |

---

---

## Phase 4: Testing S-allele Hypotheses

> Phase 4 uses the allele bin definitions from Phase 2 and the individual GFS data from Phase 3 to design and analyse controlled crossing experiments. All scripts run from the same working directory as Phase 2 and Phase 3.

### Step 21 — HV-Based Allele Hypothesis Testing and Crossing Design

**Script:** `srk_allele_hypotheses.py`

**Command:**
```bash
python srk_allele_hypotheses.py
```

**Inputs:**

| File | From step | Description |
|------|-----------|-------------|
| `SRK_protein_allele_representatives.fasta` | Step 10a | One representative sequence per allele bin (pre-aligned) |
| `SRK_individual_zygosity.tsv` | Step 12 | Per-individual tetraploid genotype pattern |
| `SRK_individual_allele_table.tsv` | Step 11 | Individual → allele → copy count table |

**Key parameters (edit at top of script):**

| Parameter | Default | Notes |
|-----------|---------|-------|
| `DOMAIN_REGION` | `(31, 430)` | S-domain columns; must match Step 10 |
| `WINDOW_SIZE` | `20` | Sliding window width (aa) for smoothing the variability profile |
| `PEAK_SD_FACTOR` | `1.0` | Positions above mean + k×SD of the smoothed profile are flagged HV; lower to include more HV positions |
| `WITHIN_CLASS_THRESHOLD` | `0.04` | HV p-distance boundary between N (synonymy test) and P_within (compatible) within a class |
| `N_GROUPS` | `None` | Override auto class detection with a fixed number of groups |
| `DISTANCE_THRESHOLD` | `None` | Override auto class detection with a fixed HV distance cutoff |
| `CROSS_TSV` | `None` | Set to cross results filename to activate Step 22 |

> **Auto class detection:** by default the script cuts the UPGMA tree at the largest gap in merge heights, which automatically identifies the Class I / Class II phylogenetic split. Inspect `SRK_HV_cluster_figure.pdf` and the printed gap table to verify the cut is biologically sensible before adjusting `WITHIN_CLASS_THRESHOLD`.

**Script workflow:**

| Part | Description |
|------|-------------|
| 1 | Moving-window per-column variability scan → identify HV positions |
| 2 | HV-only pairwise distances + UPGMA clustering → phylogenetic class split |
| 2b | Allele similarity heatmap (colour scale spans within-class range) |
| 3 | Functional group table, synonymy candidates, AAAA cross design matrix |
| 3b | Cross design summary figure (distance distribution, category schematic, N-cross interpretation) |
| 3c | Synonymy network: W-only groups figure + N-connectivity condensed figure; per-allele synonymy group CSV |
| 4 | UPGMA dendrogram + AAAA availability bar chart |
| 5 | Cross result analysis (activated when `CROSS_TSV` is set) |

**Outputs:**

| File | Content |
|------|---------|
| `SRK_variability_landscape.pdf` | Per-column and smoothed S-domain variability profile; HV regions shaded |
| `SRK_HV_allele_distances.tsv` | Pairwise HV-only distance matrix (alleles × alleles) |
| `SRK_allele_similarity_heatmap.pdf` | 63×63 HV similarity heatmap ordered by UPGMA; colour scale spans within-class range |
| `SRK_functional_allele_groups.tsv` | Allele bin → class assignment, AAAA count, cross power (full / singleton / none) |
| `SRK_synonymy_candidates.tsv` | All within-class allele pairs with HV distance, cross tier, and testability flag |
| `SRK_AAAA_cross_design_HV.tsv` | All AAAA × AAAA pairs ranked W → N → P_within → P_cross, with HV distance and expected outcome |
| `SRK_cross_design_summary.pdf` | Three-panel figure: HV distance distribution, cross category schematic, N-cross interpretation logic |
| `SRK_HV_cluster_figure.pdf` | UPGMA dendrogram (HV distances) coloured by class + AAAA availability bar chart |
| `SRK_synonymy_network_W.pdf` | Two-panel figure: HV-identical W-groups (left, individual alleles); isolated alleles with no HV-identical partner (right) |
| `SRK_synonymy_network_N.pdf` | Condensed super-node graph: each node = one W-group or isolated allele; edges = N-bridges between groups; node size ∝ individuals observed |
| `SRK_synonymy_groups.csv` | Per-allele synonymy group membership: group ID, group size, observed individual count (all genotypes), AAAA count, group total observed |
| `figures/*.png` | PNG copies of all figures at 200 dpi |

**Cross categories:**

| Category | HV distance | Definition | Expected outcome |
|----------|-------------|------------|-----------------|
| W | d = 0 | HV-identical alleles | No seeds — incompatibility predicted by sequence identity |
| N | 0 < d < threshold | Small HV differences, same class | Unknown — synonymy test: does this substitution change specificity? |
| P_within | d ≥ threshold, same class | Substantial within-class HV divergence | Seeds expected — within-class positive control |
| P_cross | different class | Between phylogenetic classes | Seeds expected — guaranteed compatible (Class I × Class II) |

> **Interpreting N-cross outcomes:** incompatible N cross → allele bins share SI specificity → merge bins (synonymous alleles). Compatible N cross → small HV difference is functionally real → bin boundary confirmed.

---

### Step 22 — Cross Result Analysis

> Requires completed crossing records. Activate by setting `CROSS_TSV` in the script, then re-run `srk_allele_hypotheses.py`.

**Script:** `srk_allele_hypotheses.py` (Part 5)

**Command:**
```bash
# Edit the script: set CROSS_TSV = "<your_cross_results_file>"
python srk_allele_hypotheses.py
```

**Inputs:**

| File | Description |
|------|-------------|
| `SRK_functional_allele_groups.tsv` | Class assignments from Step 21 |
| `SRK_HV_allele_distances.tsv` | HV distances for cross category assignment |
| `SRK_individual_allele_table.tsv` | Used to assign alleles to non-AAAA cross parents |
| Cross results file | One row per cross; columns `Mother`, `Father`, and a seed count column (auto-detected) |

**Statistical tests:**

| Test | Purpose |
|------|---------|
| Kruskal-Wallis H | Overall test of seed yield differences across W / N / P_within / P_cross |
| Mann-Whitney U (pairwise) | All category pairs |

**Outputs:**

| File | Content |
|------|---------|
| `SRK_cross_result_analysis_HV.pdf` | Seed yield box/strip plot + success rate bar chart by cross category |

---

## Key Files at Each Phase Boundary

| After step | Critical file(s) for downstream |
|------------|----------------------------------|
| Step 1 | `SRK_canonical_haplotype_sequences_revcomp_ATGfixed.fasta`, `augustus_*.csv` |
| Step 3 | `all_${Library}_Phased_haplotypes.fasta` |
| Step 7 | `all_${Library}_*_exons_backfilled.fasta` |
| Step 9 | `SRK_functional_proteins.fasta`, `SRK_functional_protein_key.tsv` |
| Step 10a | `SRK_protein_allele_assignments.tsv` |
| Step 11 | `SRK_individual_allele_genotypes.tsv` |
| Step 12 | `SRK_individual_zygosity.tsv` |
| Step 14 | `SRK_species_richness_estimates.tsv` |
| Step 16 | `SRK_TP1_summary.tsv` |
| Step 18 | `SRK_individual_GFS.tsv` |
| Step 19 | `SRK_EO_GFS_summary.tsv` |
| Step 20 | `SRK_GFS_reproductive_effort.pdf`, `SRK_GFS_reproductive_effort.png` |
| Step 21 | `SRK_functional_allele_groups.tsv`, `SRK_AAAA_cross_design_HV.tsv`, `SRK_synonymy_candidates.tsv`, `SRK_synonymy_groups.csv` |
| Step 22 | `SRK_cross_result_analysis_HV.pdf` |

---

## Metadata File Requirements (`sampling_metadata.csv`)

| Column | Used in steps | Description |
|--------|---------------|-------------|
| `Library` | 12 (optional), Class step | Library number (int or `Library001` format) |
| `barcode` | 12 (optional), Class step | Barcode number (int or `barcode01` format) |
| `SampleID` | 13, 14, 15, 17 | Individual ID matching `${Library}_${barcode}` format |
| `Pop` | 13, 14, 15, 17 | Population identifier |
| `Ingroup` | 11–15, 17 | `1` = include, `0` = exclude (outgroup) |
