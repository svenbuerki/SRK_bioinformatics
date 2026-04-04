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

**Outputs:**
- `${Library}_${barcode}_Phased_haplotypes.oriented.fasta` — oriented haplotypes per barcode
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
```

**Key parameters:**
- Minimum length: **3250 bp** (for ~4500 bp expected amplicon)
- Reference sequences are prepended before filtering so they are always retained

**Input:**
- `all_Library00X_Phased_haplotypes.fasta` — output of Step 3

**Output:**
- `all_Library00X_Phased_haplotypes_filtered_min3250.fasta` — length-filtered sequences including canonical references

---

### Step 5 — Multiple Sequence Alignment

**Tool:** MAFFT

**Command:**
```bash
mafft --auto --adjustdirection all_Library00X_Phased_haplotypes_filtered_min3250.fasta \
  > all_Library00X_Phased_haplotypes_filtered_min3250_aligned.fasta
```

> Run this command for each library. `--adjustdirection` handles any residual orientation issues.

**Input:** `all_Library00X_Phased_haplotypes_filtered_min3250.fasta`

**Output:** `all_Library00X_Phased_haplotypes_filtered_min3250_aligned.fasta`

---

### Step 6 — Exon Extraction from Multiple Sequence Alignments

**Script:** `extract_exons_with_annotation.py`

**Command:**
```bash
python extract_exons_with_annotation.py
# Prompted inputs:
#   MSA FASTA path (e.g. all_Library006_Phased_haplotypes_filtered_min3250_aligned.fasta)
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
- `all_Library00X_Phased_haplotypes_filtered_min3250_aligned_exons.fasta`

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

**Input:** `all_Library00X_Phased_haplotypes_filtered_min3250_aligned_exons.fasta`

**Output:** `all_Library00X_Phased_haplotypes_filtered_min3250_aligned_exons_backfilled.fasta`

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

**Input:** `all_Library00X_Phased_haplotypes_filtered_min3250_aligned_exons_backfilled.fasta`

**Key output:**
- `all_Library00X_Phased_haplotypes_filtered_min3250_aligned_exons_backfilled_frame1_AA_filtered_aligned.fasta`

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
| `N_ALLELES` | `None` | Option A: fix number of alleles directly |
| `DIST_THRESHOLD` | `0.01` (1%) | Option B: fix p-distance cutoff; use sensitivity plot elbow to calibrate |
| `DOMAIN_REGION` | `(31, 430)` | S-domain columns (1-based); adjust if species-specific annotation available |

> **Calibration:** inspect `SRK_protein_distance_analysis.pdf` sensitivity curve; choose the elbow (typically 1–2% for Nanopore data). Use `N_ALLELES` (Option A) if the elbow is easier to read on the y-axis than on the distance axis.

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

### Step 16 — Allele Composition Comparison Across Element Occurrences

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

---

## Metadata File Requirements (`sampling_metadata.csv`)

| Column | Used in steps | Description |
|--------|---------------|-------------|
| `Library` | 12 (optional), Class step | Library number (int or `Library001` format) |
| `barcode` | 12 (optional), Class step | Barcode number (int or `barcode01` format) |
| `SampleID` | 13, 14, 15 | Individual ID matching `${Library}_${barcode}` format |
| `Pop` | 13, 14, 15 | Population identifier |
| `Ingroup` | 11–15 | `1` = include, `0` = exclude (outgroup) |
