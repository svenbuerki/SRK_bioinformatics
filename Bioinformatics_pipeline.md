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

Phase 3 evaluates the evolutionary status and conservation implications of the SI system by integrating population genetic theory with the genotype data generated in Phases 1 and 2. Before any inferential question can be addressed at the right evolutionary scale, individuals must be assigned to the demographic units that share — and have shared — a gene-flow history. Element Occurrences (EOs) are the unit of management, but they are nested within larger evolutionary units: groups of EOs that descend from a common ancestral colonisation event and have evolved independently since fragmentation severed connectivity between them. These independent demographic units are termed **bottleneck lineages (BL1–BL5)** and are derived in a sibling spatial-clustering project (`LEPA_EO_spatial_clustering`) from a 500 m pollinator-dispersal threshold, hull-to-hull connectivity analysis, and Ward's D2 hierarchical clustering of group centroids (silhouette-optimal k = 5). Phase 3 therefore opens with **Step 13 — Element Occurrence and Bottleneck Lineage Integration**, which joins each individual to its EO, geographic group, BL, and drift index. All subsequent steps stratify in parallel by **EO sorted within BL** (the management view) and by **BL** (the evolutionary view), the latter recovering the small-locality samples that EO-level statistics cannot use. Three interconnected scientific questions then structure the analyses.

**Q1 — What are the population genetic parameters of the SRK system?**
Before interpreting SI diversity in an evolutionary framework, it is necessary to characterise each population's fundamental genetic state: allele richness, effective allele number, heterozygosity, and mean alleles per individual (Step 14). These descriptors provide the empirical foundation for all downstream comparisons.

**Q2 — What is the total allele richness of the SI system at the species level?**
Under negative frequency-dependent selection (NFDS), rare S-alleles confer a reproductive advantage because pollen bearing a rare allele can fertilise a larger proportion of compatible stigmas. This self-reinforcing dynamic drives the accumulation of allelic diversity and, at evolutionary equilibrium, pushes allele frequencies towards equality. The *species allele richness* — estimated from accumulation curves and asymptotic richness estimators (Michaelis-Menten, Chao1, iNEXT) applied to the full dataset (Step 15) — approximates this NFDS equilibrium. It represents the total pool of functionally distinct S-alleles maintained across the species and serves as the reference optimum against which population-level deficits are measured. This question is foundational: without a species-level baseline, it is impossible to determine how much allelic diversity individual populations have lost relative to the evolutionary expectation.

**Q3 — How has genetic drift eroded the SI system at the population level?**
In small or isolated populations, genetic drift counteracts NFDS by stochastically eliminating rare alleles, skewing allele frequencies away from the equal-frequency expectation, and increasing the prevalence of homozygous genotypes. Drift erodes the SI system at two hierarchical levels, each corresponding to a conservation tipping point. **Tipping Point 1 (TP1)** is breached when a population has lost so many S-alleles relative to the species optimum that inter-population allele transfers are required to restore richness. **Tipping Point 2 (TP2)** is breached when the distribution of the remaining alleles across individuals has degraded to the point that managed crossing within the population cannot recover reproductive fitness without external introductions. Phase 3 quantifies both levels through four complementary analyses:

*TP1 — Allele richness and frequency analyses (Steps 15–18):*

- **Allele accumulation curves** (Step 15) reveal how many S-alleles each population retains relative to the species optimum, measuring the depth of drift-driven allele loss.
- **Allele frequency analysis** (Step 16) tests whether population-level frequencies deviate from the equal-frequency NFDS expectation via χ² goodness-of-fit, using the species richness estimate from Step 15 as the optimum.
- **TP1 tipping point analysis** (Step 17) synthesises Steps 14–16 into a single diagnostic plot, positioning each EO by the proportion of the species optimum it retains (x axis) and by allele frequency evenness — the ratio of effective allele number to observed allele count (Ne/N) — as a measure of departure from NFDS equal-frequency expectations (y axis). EOs breaching both thresholds (< 50% of optimum and Ne/N < 0.80) are flagged CRITICAL.
- **Allele composition comparison** (Step 18) identifies which alleles are private to single populations versus shared across Element Occurrences, distinguishing populations that retain unique allelic diversity from those impoverished by isolation or founder effects.

*TP2 — Genotypic fitness analysis (Steps 19–21):*

- **Genotypic Fitness Score** (Step 19) translates allele-level diversity into individual reproductive fitness. In a tetraploid, each diploid gamete is formed by randomly sampling 2 of the 4 allele copies at the SRK locus, yielding C(4,2) = 6 equally probable gamete combinations. GFS is the proportion of those combinations that carry two distinct alleles. This reveals dosage asymmetries invisible to zygosity classification alone: an AABB individual (2+2 copies) produces 4 of 6 heterozygous gametes (GFS = 0.667), while an AAAB individual (3+1 copies) produces only 3 of 6 (GFS = 0.500), even though both carry exactly two distinct alleles. In self-incompatible plants, an individual's value as a breeding partner depends not only on which key/lock types it carries, but also on the number of distinct alleles present across its genome copies. Individuals with more key/lock types can participate in more compatible crosses, making them especially valuable in a managed breeding program.
- **TP2 Tipping Point Analysis** (Step 20) places individual-level GFS in interaction with the population-level proportion of AAAA genotypes. EOs breaching both criteria (mean GFS < 0.667 and > 30% AAAA) are flagged **CRITICAL**; those breaching either one are flagged **AT RISK**; the remainder are **OK**.
- **Reproductive effort support** (Step 21) visualises, per EO, the proportion of individuals at each GFS tier, quantifying the fraction of the population capable of contributing allelic diversity to compatible crosses.

Together, these analyses trace the consequences of demographic decline from species-level diversity baselines, through population-level allele erosion and frequency skew, to individual-level reproductive fitness — providing a quantitative framework for prioritising conservation interventions.

---

### Step 13 — Element Occurrence and Bottleneck Lineage Integration

**Script:** `SRK_BL_integration.py`

**Command:**
```bash
python3 SRK_BL_integration.py
```

**Purpose.** Bridges the SRK individual-level dataset to the spatial framework defined in the sibling project `LEPA_EO_spatial_clustering`. For each post-QC individual, the script resolves a complete sampling lineage record: which Element Occurrence (EO) it was collected from, which connected geographic group(s) that EO occupies, which independent bottleneck lineage (BL1–BL5) the group belongs to, and the relative drift index of its habitat patch. The output table is the join key used by every subsequent Phase 3 step to add the BL stratification.

**Inputs:**
- `sampling_metadata.csv` — raw sample registry (resolved to ingroup individuals only).
- `SRK_individual_zygosity.tsv` — canonical post-QC individual list (currently 272 rows). Used to filter the sampling registry to the individuals that actually completed Phase 2.
- `EO_group_BL_summary.csv` — authoritative EO → Group → BL → Drift_index mapping, produced by `LEPA_EO_spatial_clustering/EO_spatial_clustering.R`. Path is auto-resolved at `../../LEPA_EO_spatial_clustering/data/EO_group_BL_summary.csv` relative to the pipeline directory.

**Outputs:**
- `SRK_individual_BL_assignments.tsv` — one row per post-QC individual with columns `Individual`, `Pop`, `EO`, `Group`, `BL`, `Drift_index`, `BL_status`. The `BL_status` column flags each row as `Assigned` (clean numeric `Pop` matched to a known EO), `Inferred` (resolved by hyphen-prefix or manual override) or `Unassigned` (germplasm sub-code that requires user resolution against the spatial-clustering source CSV).
- A console summary showing per-BL individual counts and the EO breakdown, plus a list of any unresolved `Pop` codes.

**Key implementation notes:**
- `Pop` values that are pure numeric (e.g. `27`, `118`) are matched to the canonical EO label (e.g. `EO27`, `EO118`) only if that EO is present in the BL summary; otherwise they are flagged Unassigned.
- Compound codes such as `26-3` or `24-14` are resolved by splitting on `-` and matching the prefix to a known EO.
- The `POP_TO_EO_OVERRIDE` dictionary at the top of the script records manual mappings for irregular codes; extend it as germplasm sub-codes are cross-referenced against `Peggy_EOs_Germplasm_w_lat_long_from_Events_30Apr2026.csv` (the restricted input to the spatial-clustering project).
- The script reads all TSV/CSV inputs with `encoding="utf-8-sig"` to strip Excel-prepended BOM characters that otherwise corrupt the first column name.

**Why this step exists.** Without it, downstream analyses can only stratify by EO, which (a) misses the larger evolutionary structure — EOs are nested within BLs that represent demographically independent replicates of the bottleneck-and-drift process — and (b) is forced to discard the small localities (n = 1–3 each) for which EO-level statistics are unreliable. By aggregating those small localities into their parent BL, this step rescues those individuals for analysis at the lineage scale while keeping the EO label available for management-level reporting.

**Conservation implication.** The five BLs define the minimum number of geographic regions that must be sampled to capture the full landscape-wide allele pool. Comparing allele richness, evenness, and composition across BLs directly tests whether *independent* bottlenecks have shaped S-allele diversity in the species, and quantifies how much of the deficit measured at EO level is shared across lineages versus private to a single one.

---

### Step 14 — Population Genetics Statistics

**Script:** `SRK_population_genetic_summary.R`

**Command:**
```bash
Rscript SRK_population_genetic_summary.R
```

**Inputs:**
- `SRK_individual_allele_genotypes.tsv` — from Step 11
- `SRK_individual_zygosity.tsv` — from Step 12
- `SRK_individual_BL_assignments.tsv` — **from Step 13 (required for BL stratification)**

The script restricts analyses to BL-assigned individuals (262 of 272; the 10 unresolved germplasm sub-codes are dropped because they cannot be placed into an EO/BL).

**Key metrics computed per group (EO and BL):**
- N_alleles — allele richness
- Effective_alleles_Ne = 1/Σpᵢ² (based on `colSums` copy counts)
- Prop_heterozygous / Prop_homozygous
- Mean_alleles per individual
- Top_alleles — three most abundant alleles with copy counts

**BL color palette** (RColorBrewer Set1, matching the sibling LEPA_EO_spatial_clustering project, locked at this step and reused by Steps 15–21):
`BL1 = #E41A1C` (red), `BL2 = #377EB8` (blue), `BL3 = #4DAF4A` (green), `BL4 = #984EA3` (purple), `BL5 = #FF7F00` (orange). Unassigned individuals (if retained) use `#999999`.

**Outputs:**
- `SRK_population_genetic_summary.tsv` — EO-level table sorted by `BL → EO`. Columns: `BL`, `EO`, `Drift_index`, `N_individuals`, `N_individuals_with_alleles`, `N_alleles`, `Effective_alleles_Ne`, `Prop_heterozygous`, `Prop_homozygous`, `Mean_alleles`, `Top_alleles`, `Population` (= EO label, kept for backwards compatibility).
- `SRK_population_genetic_summary_BL.tsv` — **NEW**. BL-level aggregate table (5 rows). Columns: `BL`, `N_EOs`, `Mean_drift`, plus the same per-group metrics.
- `SRK_population_genetic_summary.pdf` — 2-page, 4-panel plot. Page 1: EO bars sorted by BL with bars colored by parent BL (Set1 palette + legend). Page 2: BL aggregate bars colored by BL.

---

### Step 15 — Allele Accumulation Curves and Richness Estimation

> **Run before Step 16.** Produces the species richness estimate used as the NFDS optimum in Step 16.

**Script:** `SRK_allele_accumulation_analysis.R`

**Command:**
```bash
Rscript SRK_allele_accumulation_analysis.R
```

**Inputs:**
- `SRK_individual_allele_genotypes.tsv` — from Step 11
- `sampling_metadata.csv` — must contain `SampleID`, `Pop`, `Ingroup` columns
- `SRK_individual_BL_assignments.tsv` — from Step 13

**Split-filter design.** This is the only step where the sample set differs by analysis level: the **species-level curve and richness estimators use ALL 272 ingroup individuals** (so the species pool baseline includes alleles unique to the 10 individuals whose germplasm sub-codes are not yet resolved to a BL — those alleles still count toward the species pool). The **EO-level and BL-level curves use only the 262 BL-assigned individuals**, where geographic provenance is required. Without the split, the species pool would shrink from 54 observed / MM=69 down to 50 observed / MM=63 — losing four alleles unique to the unresolved samples.

**Key parameters:**
- Permutations: 1000 (default)
- EOs with <5 individuals excluded from EO-level curves; all 5 BLs always included
- EO curves are sorted by parent BL and colored using the locked Set1 BL palette (Step 14)
- Three richness estimators: Michaelis-Menten (MM), Chao1, iNEXT (optional — `install.packages("iNEXT")`)

**Outputs:**
- `SRK_allele_accumulation_curves.pdf` — species + per-EO + per-BL curves on separate pages, each with MM/Chao1/iNEXT asymptote reference lines.
- `figures/SRK_allele_accumulation_species.png` — **NEW**. Standalone species curve with MM/Chao1 asymptote reference lines, y-axis extended to accommodate the asymptotes.
- `figures/SRK_allele_accumulation_combined.png` — EO accumulation curves on shared axes, **colored by parent BL** (Set1 palette); end-of-curve labels show `EOXX (observed/MM)`; legend grouped by BL.
- `figures/SRK_allele_accumulation_BL_combined.png` — **NEW**. The 5 BL aggregate accumulation curves on shared axes, colored by BL. Provides the headline visual test of "are independent bottleneck lineages approaching saturation independently?"
- `figures/SRK_allele_accumulation_drift_erosion.png` — stacked bar per EO (observed / predicted-undetected / lost-to-drift) with EOs sorted by parent BL and a **thin BL color strip along the x-axis baseline** that visually groups EOs within each BL.
- `figures/SRK_allele_accumulation_BL_drift_erosion.png` — **NEW**. Same stacked-bar decomposition aggregated to the BL level. Headline finding: BL4 has lost only ~32% of the species pool to drift; BL1 and BL2 have each lost ≥86%.
- `SRK_allele_accumulation_stats.tsv` — per-level curve statistics: rows for `Species`, then 6 EO rows (sorted by parent BL), then 5 BL rows. Columns include MM, Chao1, iNEXT estimates and sampling adequacy targets.
- `SRK_species_richness_estimates.tsv` — **consensus species allele richness; required as input for Step 16**. Computed from the all-ingroup species curve; contract unchanged from before the BL re-frame.

---

### Step 16 — Allele Frequency Analysis

> Requires `SRK_species_richness_estimates.tsv` from Step 15.

**Script:** `SRK_chisq_species_population.R`

**Command:**
```bash
Rscript SRK_chisq_species_population.R
```

**Inputs:**
- `SRK_species_richness_estimates.tsv` — from Step 15 (MM estimate used as NFDS optimum)
- `SRK_individual_allele_genotypes.tsv` — from Step 11
- `sampling_metadata.csv`
- `SRK_individual_BL_assignments.tsv` — from Step 13

**Same split-filter design as Step 15:** the species-level χ² test uses ALL 272 ingroup individuals; EO and BL χ² tests use only the 262 BL-assigned individuals.

**Key analysis:**
- χ² goodness-of-fit test vs. equal-frequency NFDS expectation
- Three levels: Species (1 row), EOs with N ≥ 5 (sorted by parent BL), BLs (all 5 always tested)
- All frequency plots share x-axis scaled to species optimum (MM estimate); missing-allele zone shown
- Each EO/BL frequency page is colored by its parent BL (Set1 palette) so the BL identity is visually obvious without checking the legend

**Outputs:**
- `SRK_chisq_species_population.tsv` — χ² statistics (X2, df, p-value) per level. Columns: `N_individuals`, `N_alleles`, `X2`, `df`, `p_value`, `Level` (Species/EO/BL), `Population`, **`BL`** (new). Rows ordered Species → EOs (sorted by BL) → BLs.
- `SRK_chisq_species_population_frequency_plots.pdf` — 12 pages: 1 species + 6 EOs (sorted by BL) + 5 BLs. Bars colored by parent BL; species page uses neutral grey.

---

### Step 17 — TP1 Tipping Point Analysis

> Requires outputs from Steps 14 and 15. Run after both have completed.

**Script:** `SRK_TP1_tipping_point.R`

**Command:**
```bash
Rscript SRK_TP1_tipping_point.R
```

**Inputs:**
- `SRK_population_genetic_summary.tsv` — from Step 14 (EO-level)
- `SRK_population_genetic_summary_BL.tsv` — from Step 14 (BL aggregates)
- `SRK_species_richness_estimates.tsv` — from Step 15

**Two analysis levels overlaid on the same scatter:**
- **EO** — circles, one per EO with N ≥ 5 individuals, colored by parent BL (Set1 palette, matching Steps 14–16 and the LEPA_EO_spatial_clustering project)
- **BL** — triangles, one per bottleneck lineage (BL1–BL5), colored by BL

**Key metrics per group:**
- `prop_optimum` = N_alleles / MM_estimate — proportion of the species-level S-allele pool retained, using the Michaelis-Menten asymptote as the species optimum.
- `evenness` = Ne / N_alleles — frequency evenness index. Ne (effective allele number) = 1/Σpᵢ², where pᵢ is the frequency of allele i; it answers *how many equally frequent alleles would produce the same level of diversity as observed*. When all alleles are at equal frequency (NFDS ideal), Ne = N and the ratio = 1.0. Genetic drift pushes the ratio downward by creating dominant alleles and marginalising rare ones. An EO/BL with evenness = 0.50, for example, has allele frequencies so uneven that only half the observed alleles are effectively contributing to SI function — the remainder are present in such low copy numbers that they are rarely expressed in crosses and are at elevated risk of loss through drift.

**TP1 thresholds (adjustable at top of script):**

| Parameter | Default | Criterion |
|-----------|---------|-----------|
| `TP1_PROP_OPTIMUM` | 0.50 | group retains < 50% of species optimum |
| `TP1_EVENNESS` | 0.80 | Ne/N_alleles below 0.80 |
| `EO_MIN_N` | 5 | minimum sample size to plot an EO point (BLs are always plotted) |

Groups breaching both criteria are flagged **CRITICAL**; those breaching either one are flagged **AT RISK**; the remainder are **OK**.

**Outputs:**
- `SRK_TP1_summary.tsv` — per-EO TP1 status with new `BL` column; rows arranged by parent BL → prop_optimum
- `SRK_TP1_summary_BL.tsv` — **NEW**, per-BL TP1 status (5 rows)
- `figures/SRK_TP1_tipping_point.png` — combined EO + BL scatter
- `SRK_TP1_tipping_point.pdf` — same as PNG, root-level (backwards compat)
- `figures/SRK_TP1_tipping_point_blank.png` — empty zones for presentations

**Headline result (current dataset, 2026-05-07):** every BL and 5 of 6 EOs are CRITICAL. Only EO18 (N = 5) escapes CRITICAL, due to small-N evenness inflation. BL4 has the highest prop_optimum (0.45) but the lowest evenness (0.33) — the diversity reservoir of the species is CRITICAL through frequency skew rather than richness loss.

**Backwards-compatibility note (2026-05-07):** the old `EO_MAP` lookup that mapped raw `Pop` values (`"27"` → `"EO27"`) was removed. Step 14 now writes `Population` already populated with the canonical EO label, so the script reads it directly.

---

### Step 18 — Allele Composition Comparison

**Script:** `SRK_allele_sharing_EOs.py`

**Command:**
```bash
python3 SRK_allele_sharing_EOs.py
```

**Inputs:**
- `SRK_individual_allele_genotypes.tsv` — from Step 11
- `SRK_individual_BL_assignments.tsv` — from Step 13 (provides EO and BL labels)

**Two analysis levels in parallel:**
- **EO** — UpSet + sharing heatmap for EOs with N ≥ 5 individuals, sorted by parent BL
- **BL** — UpSet + sharing heatmap for the 5 bottleneck lineages (no minimum N)

EO bars and dots are colored by parent BL using the locked Set1 palette (Steps 14–17). BL bars are colored by their own BL color.

**Key parameters:**

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--genotypes` | `SRK_individual_allele_genotypes.tsv` | Allele count matrix |
| `--bl-assignments` | `SRK_individual_BL_assignments.tsv` | Bridge table from Step 13 |
| `--min-samples` | `5` | Minimum individuals per EO; BLs always included |
| `--pdf-dir` | `.` | Output directory for PDFs |
| `--png-dir` | `figures` | Output directory for PNGs |

**Key analysis:**
- For each non-empty intersection of group sets, counts alleles exclusive to exactly that combination
- UpSet plot sorted by intersection size
- Pairwise heatmap shows raw shared-allele counts for every group pair (diagonal = group allele richness)
- No additional dependencies beyond the `polyploid-model` environment (pandas, numpy, matplotlib)

**Outputs (PDFs at root, PNGs in `figures/`):**
- `SRK_allele_upset_EOs.{pdf,png}` — EO UpSet (sorted by parent BL, dots colored by parent BL)
- `SRK_allele_sharing_heatmap_EOs.{pdf,png}` — pairwise EO heatmap
- `SRK_allele_upset_BLs.{pdf,png}` — **NEW**, BL-level UpSet
- `SRK_allele_sharing_heatmap_BLs.{pdf,png}` — **NEW**, 5×5 BL heatmap

**Headline result (current dataset, 2026-05-07):** 33 of 54 alleles (61%) are private to a single BL; only Allele_050 and Allele_057 are shared across all 5 BLs. The near-disjoint BL allele sets are the **direct test of independent bottlenecks** — a single shared species-level bottleneck would predict overlapping losses, not the observed lineage-private alleles. BL4 holds 14 private alleles (45% of its 31-allele complement), confirming its role as the species' diversity reservoir.

---

### Step 19 — Individual Genotypic Fitness Score (GFS)

**Script:** `SRK_individual_GFS.R`

**Command:**
```bash
Rscript SRK_individual_GFS.R
```

**Inputs:**

| File | From step | Description |
|------|-----------|-------------|
| `SRK_individual_zygosity.tsv` | Step 12 | Per-individual tetraploid genotype pattern (`AAAA`/`AAAB`/`AABB`/`AABC`/`ABCD`) |
| `SRK_individual_BL_assignments.tsv` | Step 13 | EO and BL assignments (262 BL-assigned ingroup individuals) |

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

**BL-stratified design:** EO and BL summaries are computed in parallel using the Step 13 bridge. EOs are sorted by parent BL throughout (Set1 palette, matches `LEPA_EO_spatial_clustering`); BL aggregates pool all BL-assigned individuals.

**Key parameters (top of script):**

| Parameter | Default | Notes |
|-----------|---------|-------|
| `EO_MIN_N` | `5` | Minimum sample size for an EO to appear in plots; BLs always plotted |
| `BL_PALETTE` | Set1 | `BL1=#E41A1C`, `BL2=#377EB8`, `BL3=#4DAF4A`, `BL4=#984EA3`, `BL5=#FF7F00` |
| `TP2_MEAN_GFS` | `0.667` | Adjustable threshold for mean GFS (used in Step 20) |
| `TP2_PROP_AAAA` | `0.30` | Adjustable threshold for proportion AAAA (used in Step 20) |

**Outputs:**

| File | Content |
|------|---------|
| `SRK_individual_GFS.tsv` | Per-individual GFS, genotype class, EO, **BL, Drift_index**, and method |

---

### Step 20 — TP2 Tipping Point Analysis

> Requires `SRK_individual_GFS.tsv` from Step 19. Performed by the same script (`SRK_individual_GFS.R`); TP2 flags are computed at both the EO level and the BL level.

**Script:** `SRK_individual_GFS.R` (same run as Step 19)

**Key analysis:**

Places individual-level GFS in interaction with the proportion of AAAA genotypes to evaluate whether allele distributions have degraded beyond the reach of managed crossing. Summaries are evaluated against two thresholds at both EO and BL levels:

| Threshold | Default | Criterion |
|-----------|---------|-----------|
| `TP2_MEAN_GFS` | 0.667 | Mean GFS below AABB level |
| `TP2_PROP_AAAA` | 0.30 | >30% AAAA individuals |

Groups breaching both thresholds simultaneously are flagged **CRITICAL**; those breaching either one are flagged **AT RISK**; the remainder are **OK**. The TP2 figure plots EOs as circles and BLs as triangles on the same axes, both coloured by parent BL using the locked Set1 palette.

**Headline:** all 5 BLs are CRITICAL on TP2 (BL3 worst at mean GFS = 0.196, prop AAAA = 0.662); 5 of 6 plotted EOs are CRITICAL (EO67 AT RISK only).

**Outputs:**

| File | Content |
|------|---------|
| `SRK_EO_GFS_summary.tsv` | EO-level GFS statistics, class proportions, TP2 status, **BL column** |
| `SRK_BL_GFS_summary.tsv` | **NEW** — BL-level GFS statistics, class proportions, TP2 status |
| `SRK_GFS_plots.pdf` | Multi-page PDF: composition (proportional + counts), jitter, TP2 scatter — EOs sorted by parent BL with BL-coloured x-axis labels |
| `figures/SRK_GFS_plots_p1_composition_proportional.png` | Stacked-proportional bar (EOs sorted by BL) |
| `figures/SRK_GFS_plots_p2_individual_jitter.png` | Per-individual GFS jitter with EO means |
| `figures/SRK_GFS_plots_p3_TP2_tipping_point.png` | TP2 scatter — EOs (circles) + BLs (triangles), both Set1-coloured |
| `figures/SRK_GFS_plots_p3_TP2_tipping_point_blank.png` | Same plot without data points, for presentation use |
| `figures/SRK_GFS_plots_p4_composition_counts.png` | Stacked-count bar (EOs sorted by BL) |

---

### Step 21 — Reproductive Effort Support per Element Occurrence and Bottleneck Lineage

> Requires `SRK_individual_GFS.tsv` from Step 19 (now carries the BL column from the Step 13 bridge).

**Script:** `SRK_GFS_reproductive_effort.R`

**Command:**
```bash
Rscript SRK_GFS_reproductive_effort.R
```

**Inputs:**

| File | From step | Description |
|------|-----------|-------------|
| `SRK_individual_GFS.tsv` | Step 19 | Per-individual GFS scores, genotype classes, **EO + BL + Drift_index** |
| `SRK_individual_zygosity.tsv` | Step 12 | Per-individual allele composition used to identify the single allele carried by AAAA individuals |

**Key analysis:**

Visualises, for each focus EO and each BL, the proportion of individuals at each GFS tier as a horizontal proportional bar chart. In a self-incompatible plant, an individual's value as a breeding partner depends on the number of distinct alleles present across its genome copies — individuals with more key/lock types can participate in more compatible crosses. The figures distinguish individuals that can contribute allelic diversity to gametes (GFS > 0) from those that cannot (AAAA, GFS = 0), and annotate each row with proportion supporting, count, and mean GFS.

**BL-stratified design:** EOs are sorted by parent BL then by mean GFS within BL; y-axis labels are coloured by parent BL (Set1 palette). The BL panel pools all BL-assigned individuals into 5 lineage rows, sorted by mean GFS. The TP2 AAAA threshold (30%) is marked on both panels.

**AAAA allele identity panels:** unpack the AAAA bar by showing which alleles each AAAA individual carries. *Pan-BL* alleles are present in AAAA individuals across every BL; *pan-EO* alleles are present in AAAA individuals in ≥80% of focus EOs (relaxed threshold accommodates EO18, the smallest focus EO with only 4 AAAA individuals carrying a different drift signature). Allele_050 and Allele_057 (W-group 1, HV-identical, likely the same SI specificity) appear as pan-BL across all 5 BLs and pan-EO in 5 of 6 focus EOs — confirming that independent bottlenecks have all converged on the same fixed SI specificity.

**Headline:** BL3 has the lowest reproductive effort support (34%, mean GFS = 0.196); BL2 has the highest W-group 1 saturation among AAAA individuals (69%); EO18 is the only focus EO without W-group 1 alleles in its AAAA pool (4 AAAA individuals, 4 distinct non-W-group-1 alleles).

**Key parameters (top of script):**

| Parameter | Default | Notes |
|-----------|---------|-------|
| `EO_MIN_N` | `5` | Minimum sample size for an EO to appear in plots (matches Steps 17/19) |
| `BL_PALETTE` | Set1 | `BL1=#E41A1C`, `BL2=#377EB8`, `BL3=#4DAF4A`, `BL4=#984EA3`, `BL5=#FF7F00` |
| `TP2_PROP_AAAA` | `0.30` | Threshold line position (matches TP2 threshold from Step 20) |
| `PAN_EO_FRAC` | `0.80` | Allele must appear in AAAA individuals in ≥80% of focus EOs to be flagged pan-EO |

**Outputs:**

| File | Content |
|------|---------|
| `SRK_GFS_reproductive_effort.pdf` | 2-page PDF: EO panel + BL panel — GFS tier composition with reproductive effort annotations |
| `SRK_GFS_AAAA_allele_composition.pdf` | 2-page PDF: EO panel + BL panel — allele identity of AAAA individuals; W-group 1 alleles (Allele_050, Allele_057) highlighted |
| `SRK_BL_reproductive_effort_summary.tsv` | **NEW** — BL-level reproductive effort summary (n, n_supporting, prop_supporting, prop_AAAA, mean_GFS) |
| `figures/SRK_GFS_reproductive_effort_EO.png` | EO-level proportional bar (PNG) |
| `figures/SRK_GFS_reproductive_effort_BL.png` | **NEW** — BL-level proportional bar (PNG) |
| `figures/SRK_GFS_AAAA_allele_composition_EO.png` | EO-level AAAA allele identity (PNG) |
| `figures/SRK_GFS_AAAA_allele_composition_BL.png` | **NEW** — BL-level AAAA allele identity (PNG) |

---

---

## Phase 4: Testing S-allele Hypotheses

> Phase 4 uses the allele bin definitions from Phase 2 and the individual GFS data from Phase 3 to design and analyse controlled crossing experiments. All scripts run from the same working directory as Phase 2 and Phase 3.

### Step 22 — HV-Based Allele Hypothesis Testing and Crossing Design

> Step 22 is now a **two-script workflow**. First, the variability landscape is built from a multi-species alignment (LEPA + Brassica + Arabidopsis SRK alleles) and HV columns are detected by Shannon entropy with a permutation-tested cross-Brassicaceae overlap. The validated LEPA HV columns are then consumed by the allele-hypothesis testing script (HV-distance matrix → UPGMA classes → cross design → synonymy network).

#### Step 22a — Cross-Brassicaceae S-domain Variability Landscape

**Scripts:**
- `srk_fetch_reference_alleles.py` — one-time NCBI fetch of *Brassica rapa*, *B. oleracea*, *Arabidopsis lyrata*, and *A. halleri* SRK reference proteins (≈22 + 10 sequences after S-haplotype dedup).
- `srk_brassica_hv_mapping.py` — maps the 12 SCR9-contact residues from Ma et al. 2016 (PDB 5GYY, *B. rapa* eSRK9) to LEPA alignment coordinates via `mafft --add`.
- `srk_variability_landscape.py` — primary script: combined alignment + per-species Shannon entropy + per-species HV-region calls + permutation tests + figure.

**Commands:**
```bash
# One-time reference acquisition + alignment (skip if already done)
python3 srk_fetch_reference_alleles.py
mafft --add all_reference_SRKs_dedup.fasta SRK_protein_allele_representatives_padded.fasta \
      > SRK_combined_alignment.fasta
python3 srk_brassica_hv_mapping.py

# Variability landscape (re-run whenever LEPA representatives change)
python3 srk_variability_landscape.py
```

**Inputs:**

| File | From | Description |
|------|------|-------------|
| `SRK_combined_alignment.fasta` | mafft --add | LEPA representatives + Brassica + Arabidopsis SRK references |
| `SRK_brassica_hv_mapping.tsv` | `srk_brassica_hv_mapping.py` | Ma 2016 SCR9-contact residues mapped to LEPA columns |

**Key parameters (top of `srk_variability_landscape.py`):**

| Parameter | Default | Notes |
|-----------|---------|-------|
| `DOMAIN_REGION` | `(31, 430)` | S-domain columns (must match Step 10) |
| `WINDOW_SIZE` | `20` | Sliding-window width for smoothing |
| `PEAK_SD_FACTOR` | `1.0` | HV threshold = mean + k × SD on smoothed Shannon entropy |
| `MIN_HV_RUN` | `3` | Minimum consecutive HV columns to call a region (drops singletons) |
| `PERM_N` | `10000` | Permutations for HV-overlap significance test |

**Outputs:**

| File | Content |
|------|---------|
| `SRK_variability_landscape.pdf` / `figures/SRK_variability_landscape.png` | Multi-panel: per-species smoothed Shannon entropy; per-species HV-region tracks; Ma 2016 SCR9-contact residue markers |
| `SRK_HV_regions_per_species.tsv` | LEPA / Brassica / Arabidopsis HV runs (start–end, length, threshold) |
| `SRK_LEPA_HV_positions.tsv` | **Canonical LEPA HV columns** (one row per column, 1-based alignment coordinate). Consumed by Step 22b. |
| `SRK_HV_overlap_permutation.tsv` | Pairwise HV-overlap permutation tests (LEPA↔Brassica, LEPA↔Arabidopsis, Brassica↔Arabidopsis) |

**Headline:** All three pairwise HV overlaps are highly significant (LEPA↔Brassica obs=17 cols, p=0.0005; LEPA↔Arabidopsis 22 cols, p=0.0024; Brassica↔Arabidopsis 43 cols, p<0.0001). All 11 of the 12 mappable Ma 2016 SCR9-contact residues fall within the cross-species shared HV regions — independent structural validation that the variability scan is detecting the SI-specificity surface. LEPA additionally has a unique HV peak at LEPA cols 358–388 not present in Brassica/Arabidopsis (candidate Lepidium-specific specificity site). Wu-Kabat sanity check: Jaccard with Shannon-entropy HV cols = 0.45 (LEPA), 0.88 (Brassica), 0.75 (Arabidopsis) — strong concordance under both metrics.

---

#### Step 22b — Allele Hypothesis Testing and Crossing Design

**Script:** `srk_allele_hypotheses.py`

**Command:**
```bash
python3 srk_allele_hypotheses.py
```

**Inputs:**

| File | From step | Description |
|------|-----------|-------------|
| `SRK_protein_allele_representatives.fasta` | Step 10a | One representative sequence per allele bin (pre-aligned) |
| `SRK_individual_zygosity.tsv` | Step 12 | Per-individual tetraploid genotype pattern |
| `SRK_individual_allele_table.tsv` | Step 11 | Individual → allele → copy count table |
| `SRK_LEPA_HV_positions.tsv` | Step 22a | **Canonical HV columns (overrides internal scan when present)** |
| `SRK_brassica_hv_mapping.tsv` | `srk_brassica_hv_mapping.py` | Ma 2016 markers for the variability-landscape page (optional) |

**Key parameters (edit at top of script):**

| Parameter | Default | Notes |
|-----------|---------|-------|
| `DOMAIN_REGION` | `(31, 430)` | S-domain columns; must match Step 10 |
| `WINDOW_SIZE` | `20` | Sliding-window width for the legacy internal scan (used only as fallback) |
| `PEAK_SD_FACTOR` | `1.0` | Threshold for the legacy internal scan |
| `WITHIN_CLASS_THRESHOLD` | `0.04` | HV p-distance boundary between N (synonymy test) and P_within (compatible) |
| `N_GROUPS` | `None` | Override auto class detection with a fixed number of groups |
| `DISTANCE_THRESHOLD` | `None` | Override auto class detection with a fixed HV distance cutoff |
| `CROSS_TSV` | `None` | Set to cross results filename to activate Step 23 |
| `LEPA_HV_POSITIONS_FILE` | `SRK_LEPA_HV_positions.tsv` | When this file exists, its HV columns OVERRIDE the internal sliding-window scan |

> **HV-column source of truth:** Parts 2–5 (HV-distance matrix, UPGMA, cross design, synonymy network) operate on the columns listed in `SRK_LEPA_HV_positions.tsv`. The internal sliding-window scan in Part 1 still runs (and produces the legacy figure) but its HV calls are *overridden* whenever the canonical file exists. This guarantees that every downstream test uses the cross-Brassicaceae-validated HV columns from Step 22a.

> **Auto class detection:** by default the script cuts the UPGMA tree at the largest gap in merge heights, which automatically identifies the Class I / Class II phylogenetic split. Inspect `SRK_HV_cluster_figure.pdf` and the printed gap table to verify the cut is biologically sensible before adjusting `WITHIN_CLASS_THRESHOLD`.

**Script workflow:**

| Part | Description |
|------|-------------|
| 1 | Sanity-check internal variability scan (overridden by `SRK_LEPA_HV_positions.tsv` if present) |
| 2 | HV-only pairwise distances + UPGMA clustering → phylogenetic class split |
| 2b | Allele similarity heatmap (colour scale spans within-class range) |
| 3 | Functional group table, synonymy candidates, AAAA cross design matrix |
| 3b | Cross design summary figure (distance distribution, category schematic, N-cross interpretation) |
| 3c | Synonymy network: W-only groups figure + N-connectivity condensed figure; per-allele synonymy group CSV |
| 4 | UPGMA dendrogram + AAAA availability bar chart |
| 5 | Cross result analysis (activated when `CROSS_TSV` is set; → Step 23) |

**Outputs:**

| File | Content |
|------|---------|
| `SRK_HV_allele_distances.tsv` | Pairwise HV-only distance matrix (alleles × alleles) computed on the canonical HV columns |
| `SRK_allele_similarity_heatmap.pdf` | HV similarity heatmap ordered by UPGMA; colour scale spans within-class range |
| `SRK_functional_allele_groups.tsv` | Allele bin → class assignment, AAAA count, cross power (full / singleton / none) |
| `SRK_synonymy_candidates.tsv` | All within-class allele pairs with HV distance, cross tier, and testability flag |
| `SRK_AAAA_cross_design_HV.tsv` | AAAA × AAAA pairs ranked W → N → P_within → P_cross, with HV distance and expected outcome |
| `SRK_cross_design_summary.pdf` | Three-panel figure: HV distance distribution, cross category schematic, N-cross interpretation |
| `SRK_HV_cluster_figure.pdf` | UPGMA dendrogram coloured by class + AAAA availability bar chart |
| `SRK_synonymy_network_W.pdf` | HV-identical W-groups + isolated singletons |
| `SRK_synonymy_network_N.pdf` | Condensed super-node graph (W-groups + N-bridge edges; node size ∝ individuals) |
| `SRK_synonymy_groups.csv` | Per-allele W-group membership and counts |
| `figures/*.png` | PNG copies of all figures at 200 dpi |

**Headline (current dataset, 73 canonical HV columns):** UPGMA splits 63 alleles into 2 functional groups (FG01 = 62 alleles, FG02 = Allele_061 outlier). Synonymy network: **9 W-groups + 23 isolated singletons → 32 effective bins** (down from 63), driven by W-group 1 (13 HV-identical alleles incl. Allele_050 + Allele_057 = the pan-BL fixed S-specificity). 5663 W pairs, 6365 N pairs (synonymy tests), 1338 P_within pairs in the cross-design matrix.

**Cross categories:**

| Category | HV distance | Definition | Expected outcome |
|----------|-------------|------------|-----------------|
| W | d = 0 | HV-identical alleles | No seeds — incompatibility predicted by sequence identity |
| N | 0 < d < threshold | Small HV differences, same class | Unknown — synonymy test: does this substitution change specificity? |
| P_within | d ≥ threshold, same class | Substantial within-class HV divergence | Seeds expected — within-class positive control |
| P_cross | different class | Between phylogenetic classes | Seeds expected — guaranteed compatible (Class I × Class II) |

> **Interpreting N-cross outcomes:** incompatible N cross → allele bins share SI specificity → merge bins (synonymous alleles). Compatible N cross → small HV difference is functionally real → bin boundary confirmed.

---

### Step 23 — Cross Result Analysis

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
| `SRK_functional_allele_groups.tsv` | Class assignments from Step 22 |
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
| Step 15 | `SRK_species_richness_estimates.tsv` |
| Step 17 | `SRK_TP1_summary.tsv` |
| Step 19 | `SRK_individual_GFS.tsv` |
| Step 20 | `SRK_EO_GFS_summary.tsv` |
| Step 21 | `SRK_GFS_reproductive_effort.pdf`, `SRK_GFS_AAAA_allele_composition.pdf` |
| Step 22 | `SRK_functional_allele_groups.tsv`, `SRK_AAAA_cross_design_HV.tsv`, `SRK_synonymy_candidates.tsv`, `SRK_synonymy_groups.csv` |
| Step 23 | `SRK_cross_result_analysis_HV.pdf` |

---

## Metadata File Requirements (`sampling_metadata.csv`)

| Column | Used in steps | Description |
|--------|---------------|-------------|
| `Library` | 12 (optional), Class step | Library number (int or `Library001` format) |
| `barcode` | 12 (optional), Class step | Barcode number (int or `barcode01` format) |
| `SampleID` | 13, 14, 15, 17 | Individual ID matching `${Library}_${barcode}` format |
| `Pop` | 13, 14, 15, 17 | Population identifier |
| `Ingroup` | 11–15, 17 | `1` = include, `0` = exclude (outgroup) |
