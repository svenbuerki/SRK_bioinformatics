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

## Output organisation (refactored 2026-06-14)

Every Phase 2 and Phase 3 script reads and writes from a fixed Phase- and step-prefixed path. This makes every input/output unambiguously traceable to the step that produced it, eliminates root-level clutter, and avoids stale-copy bugs.

```
Tables/
├── sampling_metadata.csv          ← canonical metadata (symlinked from root)
├── Phase2/
│   ├── step9_*.{fasta,tsv,txt}    ← functional protein outputs
│   ├── step10a_*.tsv,fasta        ← allele clustering
│   ├── step10b_*.tsv              ← variable positions
│   ├── step11_*.tsv               ← genotyping
│   ├── step12_*.tsv               ← zygosity
│   └── step12c_*.tsv,csv          ← data-quality audit + lab deliverables
└── Phase3/
    ├── step13_*.tsv               ← BL integration
    ├── step14_*.tsv               ← population genetics
    ├── step15_*.tsv               ← accumulation curves + species richness
    ├── step16_*.tsv               ← χ² allele frequency
    ├── step17_*.tsv               ← TP1 / P_compat / DI metrics
    ├── step19_*.tsv               ← individual GFS
    ├── step20_*.tsv               ← TP2 EO + BL summaries
    └── step21_*.tsv               ← reproductive effort BL summary

figures/
├── Phase2/
│   ├── step10a_protein_distance_analysis.{pdf,png}
│   ├── step10b_AA_mutation_heatmap.{pdf,png}
│   ├── step10b_AA_frequency_heatmap.{pdf,png}
│   ├── step12_zygosity_distribution.{pdf,png}
│   ├── step12c_library_effect_summary.{pdf,png}
│   ├── step12c_data_quality_per_BL.{pdf,png}
│   └── step12c_data_quality_per_EO.{pdf,png}
└── Phase3/
    ├── step14_population_genetic_summary.{pdf,png}
    ├── step15_allele_accumulation_curves.pdf  (multi-page)
    ├── step15_allele_accumulation_species.{pdf,png}
    ├── step15_allele_accumulation_{combined,BL_combined,drift_erosion,BL_drift_erosion}.png
    ├── step16_chisq_species_population_frequency_plots.pdf  (multi-page)
    ├── step17_P_compat_traffic_light_EO{,_blank}.{pdf,png}
    ├── step17_depletion_ranking_{observed,predicted}_{blank,EOs,all}.{pdf,png}
    ├── step18_allele_{upset,sharing_heatmap}_{EOs,BLs}.{pdf,png}
    ├── step18_allele_eulerr_BLs.{pdf,png}
    ├── step19_GFS_plots.pdf + p1/p2/p4 PNGs
    ├── step20_TP2_tipping_point{,_blank}.png
    ├── step21_GFS_reproductive_effort{,_EO,_BL}.{pdf,png}
    ├── step21_GFS_AAAA_allele_composition{,_EO,_BL}.{pdf,png}
    └── presentation/
        └── step15_allele_accumulation_species_{observed,predicted}.png
```

**Conventions:**
- Every plot script produces BOTH `.pdf` and `.png` (vector + raster) except multi-page PDFs (Step 15 curves, Step 16 frequency, Step 19 GFS) where the page-level PNGs are companions.
- Step prefix `stepN_` (lowercase, integer). Sub-step letters (`step10a`, `step12c`, etc.) retain their letter suffix.
- Every script's READ path includes the upstream step's prefix; this makes the dependency chain visible directly in the source.
- `Tables/sampling_metadata.csv` is the **single canonical metadata file**; a symlink at `./sampling_metadata.csv` preserves backward compatibility for legacy scripts.

**Phase 4 (Steps 25–28 + Step 14b/17b/19b/20b null-aware reruns) refactored 2026-06-14; promoted from Phase 5 to Phase 4 on 2026-06-15 because per-individual SI status is a prerequisite for cross design.** Tables in `Tables/Phase4/`, figures in `figures/Phase4/`. Includes:

- `step25a_null_allele_assignments.tsv` + intermediate FASTAs/alignments
- `step25b_individual_SI_status.tsv` + `step25b_SI_status_{species,by_BL,by_EO}_{full,robust}.png`
- `step26_individual_{allele_genotypes,zygosity}_with_nulls.tsv` + `step26_samples_for_redo.tsv`
- `step14b_population_genetic_summary{,_BL}_with_nulls.tsv` + `step14b_population_genetic_summary_with_nulls.pdf`
- `step17b_EO_allele_richness_with_nulls.tsv` + `step17b_P_compat_traffic_light_with_nulls{,_blank}.{pdf,png}` + `step17b_depletion_ranking_{observed,predicted}_with_nulls_{blank,EOs,all}.{pdf,png}`
- `step19b_individual_GFS_with_nulls.tsv` + `step19b_GFS_with_nulls_composition.png` + `step20b_TP2_with_nulls_scatter.png`
- `step27_inheritance_{trajectories,time_to_sc}.tsv` + `step27_inheritance_{pNULL_trajectories,SC_progression,time_to_sc}.png`
- `step28_injection_donor_ranking.tsv` + `step28_donor_recovery_ladder.{pdf,png}`

**Phase 4 cascade dependency reminder.** The Phase 4 chain depends on Phase 3 outputs (`Tables/Phase3/step13_individual_BL_assignments.tsv`, `Tables/Phase2/step12_individual_zygosity.tsv`). Run order: 25b → 25a → 25b refresh → 25b figures → 26 → {14b, 17b, 19b/20b} → 27 → 27 figures → 28 → 28 figures.

**Phase 5 (Steps 22a–22e + 23, cross design) refactored 2026-06-14; demoted from Phase 4 to Phase 5 on 2026-06-15 because cross design must consume the SI-status-validated parent list from Phase 4.** External reference FASTAs (Brassica + Arabidopsis SRK alleles) live in **`FASTA/`** at the repo root. Intermediate Phase 5 FASTAs and step-prefixed TSVs land in `Tables/Phase5/`; figures in `figures/Phase5/`. `pad_representatives.py` filters Step 10a representatives to LEPA-ingroup-observed alleles only (49 of 55) — eliminates outgroup contamination from the variability-landscape analysis.

```
FASTA/                  ← external reference FASTAs (Brassica + Arabidopsis SRK)
Tables/Phase5/          ← step22a_*.{fasta,tsv}, step22b_*.{tsv,csv}, step22c_*.tsv, step22d_*.tsv, step22e_*.tsv
figures/Phase5/         ← step22a_variability_landscape.{pdf,png}, step22b_*.{pdf,png},
                          step22c_*.{pdf,png}, step22d_*.{pdf,png}, step22e_cross_plan_summary.{pdf,png}, step23_*.pdf
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

### Pipeline wrapper (recommended) — `run_within_library_suite.sh`

**Script:** `run_within_library_suite.sh`

Steps 2 → 8 form a strictly sequential within-library suite. The wrapper runs them in order for every library in a config TSV, feeds each script's prompted inputs automatically, captures per-library logs, and continues past failures.

**Command:**
```bash
./run_within_library_suite.sh                       # report-mode chimera survey
CHIMERA_FILTER_MODE=filter ./run_within_library_suite.sh   # production
```

**Config file (`libraries.tsv`)** — tab-separated `library<TAB>start_barcode<TAB>end_barcode`. If absent, the wrapper auto-generates it by scanning `Library*/barcode??` and listing each library with its first/last barcode. Lines starting with `#` are skipped (use to disable a library on a given run).

**Resume / idempotency:**

| Granularity | Default behaviour | Override |
|---|---|---|
| Per-step (Steps 3–8 outputs) | Skip if expected output file exists and is non-empty | `FORCE=1` |
| Per-sample (Step 2 `*_Phased_haplotypes.fasta`) | Skip whole barcode if final FASTA exists | `FORCE_SAMPLE=1` |
| Per-Canu-rep (inside Step 2) | Skip rep if `CANU_rep<N>/*.contigs.fasta` exists | `FORCE_CANU=1` |

Together these let an interrupted run resume cleanly: just re-invoke the wrapper, and only the missing outputs are computed.

**Env vars:**

| Variable | Default | Purpose |
|---|---|---|
| `CHIMERA_FILTER_MODE` | `report` | `report` = pass-through + per-contig stats TSV; `filter` = drop chimeric contigs |
| `MIN_MEAN_COV` | `20` | mean depth floor for chimera filter |
| `MIN_UNIFORMITY` | `0.2` | `min_window_cov / median_cov` floor (calibrated 2026-06-12 on Library 005) |
| `CHIMERA_WINDOW` | `200` | sliding-window size (bp) for `min_window_cov` |
| `START_AT_STEP` / `STOP_AT_STEP` | `2` / `8` | restrict the step range |
| `FORCE`, `FORCE_SAMPLE`, `FORCE_CANU` | `0` | bypass the respective skip layers |

**Outputs (per run):**
- `logs/wrapper_<timestamp>/<Library>.log` — full per-library log capturing every step's stdout/stderr
- `logs/wrapper_<timestamp>/master_summary.tsv` — `library / status / last_completed_step / failed_step / log_file`

Step-by-step descriptions below remain authoritative; Step 2's "Key behaviours" block has been augmented with the chimera filter, Canu cache resume, and sample-skip logic.

---

### Step 2 — Nanopore Amplicon Assembly and Phasing

**Script:** `nanopore_assembly_pipeline_barcode_range.sh`

**Command (when running individually outside the wrapper):**
```bash
./nanopore_assembly_pipeline_barcode_range.sh
# Prompted inputs: LibraryID (e.g. Library001), start barcode, end barcode
```

> **Important:** LibraryID format must be `Library001` (no underscores).

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

**Embedded coverage-based chimera filter** (helper: `chimera_coverage_filter.py`)

After Canu rep merging + `seqkit rmdup` and before RACON polishing, all unique contigs are reference-aligned with minimap2, per-base depth is read with `samtools depth -a`, and per-contig statistics are computed:

| Metric | Definition |
|---|---|
| `mean_cov` | mean read depth across the contig |
| `median_cov` | median read depth |
| `min_window_cov` | minimum of the mean depth across overlapping 200 bp sliding windows |
| `uniformity` | `min_window_cov / median_cov` — chimeric contigs collapse at a junction, depressing this value |

A contig is **DROP** when `mean_cov < MIN_MEAN_COV` OR `uniformity < MIN_UNIFORMITY`.

In `report` mode (default during calibration) every contig is retained and the per-contig TSV is written to `logs/${library}_${run_id}_${barcode}.chimera_stats.tsv`. In `filter` mode DROP contigs are removed before RACON.

Library 005 calibration (2026-06-12): 5 barcodes, 54 contigs, 14 DROP (26 %). Uniformity values form a clean gap at 0.20 (lowest KEEP = 0.202, highest DROP = 0.182), validating the default threshold. Two extreme DROPs (uniformity < 0.10) showed textbook depth-crash signatures at the chimeric junction.

**Failure tolerance:**
The script wraps every per-sample command in a `try_step` helper that captures stdout / stderr to a per-sample log and stops only the affected sample on failure — the loop continues with the next barcode. A summary TSV (`logs/${library}_${start}-${end}_${run_id}.summary.tsv`) records `sample / status / failed_step / log_file`.

**Resume layers:**
- **Per-sample skip** — if `${library}_${barcode}_Phased_haplotypes.fasta` already exists, the whole barcode is skipped (override with `FORCE_SAMPLE=1`).
- **Canu cache** — each `CANU_rep<N>` directory is reused if its `*.contigs.fasta` exists (override with `FORCE_CANU=1`).

**Outputs per barcode:**
- `combined_unique_contigs.fasta` — deduplicated contigs from all 4 CANU runs
- `combined_unique_contigs_chimera_filtered.fasta` — chimera filter output (identical to input in `report` mode)
- `chimera_depth.tsv` — per-base depth for the chimera filter
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

### Step 4b — Reference-Similarity Filter (BLAST + coverage)

**Script:** `filter_by_reference_similarity.sh`

**Rationale:** Length-based filtering alone does not separate true SRK alleles from SRK paralogs. Library 010 (2026-05-11) revealed a bimodal length distribution: a true SRK cluster (~3300 bp) and a paralog cluster (~3900 bp) carrying a ~500 bp insertion. Both clusters are 98–100% identical where they align, so identity does not discriminate. The discriminator is BLAST **query coverage**: true SRK aligns over ≥90% of query length, paralogs over ≤85% (the inserted region does not align to canonical SRK).

Without this filter, paralogs inflate the Step 5 MSA (Library 010: 5690 bp wide; 40–45% gaps per row) and produce spurious gap blocks in the Step 8 AA alignment.

**Command:**
```bash
./filter_by_reference_similarity.sh
# Prompted inputs:
#   Length-filtered FASTA (e.g. all_Library010_Phased_haplotypes_filtered_min3250_max4000.fasta)
#   Canonical SRK reference FASTA (e.g. Canonical_sequences/SRK_canonical_haplotype_sequences_revcomp_ATGfixed.fasta)
#   Minimum query coverage [default 0.90]
#   Minimum percent identity [default 95]
#   BLAST threads [default 4]
```

**Key parameters:**
- `min_coverage = 0.90` — paralogs with the ~500 bp insertion bottom out at ~85% coverage; 0.90 cleanly separates them from true SRK alleles (which cluster at 95–100%).
- `min_identity = 95` — rejects spurious low-identity hits but does not attempt to use identity for paralog discrimination (paralogs are 98–100% identical where they align).
- Canonical reference IDs are always retained regardless of BLAST result.

**Input:**
- `all_Library00X_Phased_haplotypes_filtered_min3250_max4000.fasta` — output of Step 4
- `Canonical_sequences/SRK_canonical_haplotype_sequences_revcomp_ATGfixed.fasta` — 4 LEPA canonical SRK haplotypes (BEA hybrid)

**Output:**
- `all_Library00X_Phased_haplotypes_filtered_min3250_max4000_blastfilt.fasta` — paralog-purged FASTA, input to Step 5
- `all_Library00X_Phased_haplotypes_filtered_min3250_max4000_blastfilt_log.tsv` — per-query decision log (`query_id`, `query_length`, `pct_identity`, `coverage`, `decision`)

**Validation (Library 010):**
| Metric | Before Step 4b | After Step 4b |
|---|---:|---:|
| Sequences | 2260 | 1516 (–744 paralog/chimeric) |
| Step 5 MSA width | 5690 bp | 3735 bp |
| Per-row gap fraction | 40–45% | <20% |

Runtime: ~3 s on 2260 sequences (4 threads).

---

### Step 5 — Multiple Sequence Alignment

**Tool:** MAFFT

**Command:**
```bash
mafft --auto --adjustdirection all_Library00X_Phased_haplotypes_filtered_min3250_max4000_blastfilt.fasta \
  > all_Library00X_Phased_haplotypes_filtered_min3250_max4000_blastfilt_aligned.fasta
```

> Run this command for each library. `--adjustdirection` handles any residual orientation issues.

**Input:** `all_Library00X_Phased_haplotypes_filtered_min3250_max4000_blastfilt.fasta` (output of Step 4b)

**Output:** `all_Library00X_Phased_haplotypes_filtered_min3250_max4000_blastfilt_aligned.fasta`

---

### Step 6 — Exon Extraction from Multiple Sequence Alignments

**Script:** `extract_exons_with_annotation.py`

**Command:**
```bash
python extract_exons_with_annotation.py
# Prompted inputs:
#   MSA FASTA path (e.g. all_Library006_Phased_haplotypes_filtered_min3250_max4000_blastfilt_aligned.fasta)
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
- `all_Library00X_Phased_haplotypes_filtered_min3250_max4000_blastfilt_aligned_exons.fasta`

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

**Input:** `all_Library00X_Phased_haplotypes_filtered_min3250_max4000_blastfilt_aligned_exons.fasta`

**Output:** `all_Library00X_Phased_haplotypes_filtered_min3250_max4000_blastfilt_aligned_exons_backfilled.fasta`

---

### Step 7b — Ambiguity (N-content) Filter

**Script:** `filter_by_ambiguity.sh`

**Rationale (causal chain):**

1. **Upstream origin — short Nanopore fragments.** A subset of barcodes in any given library produces Nanopore reads that are shorter than the expected ~3.5 kb SRK amplicon (fragmentation during library prep, partial pore translocation, or end-of-read truncation). When Canu performs *de novo* assembly on these reads, it cannot reconstruct the full amplicon span and instead produces a contig that is shorter than the canonical reference at one (usually the 3′) terminus.

2. **Step 7 backfill writes `N`s.** Step 7 (`backfill_alignment_ends.py`) handles short contigs by **padding the missing terminal positions with `N`** so that all sequences share the same alignment width as the canonical reference. This is correct behaviour at the DNA stage — the `N`s preserve column registry without inventing nucleotides.

3. **Step 8 translates `N` codons to `X`.** Any codon containing `N` translates to the ambiguity character `X`. The 3′ N-block therefore appears as a contiguous block of `X` residues at the C-terminus of the predicted protein.

4. **MAFFT inflates the AA alignment.** When the AA MSA is built, MAFFT cannot match `X` to a defined residue, so it inserts gap columns to accommodate the `X` block. Those gaps appear on the **BEA canonical reference rows** (which have proper residues there), producing the visible trailing-gap artefact.

Filtering at the DNA stage (this step), before translation, removes the affected contigs cleanly and prevents the downstream cascade.

Library 010 (2026-05-11) illustrates the pattern: 128 of 1516 sequences (~8%) carried 10–100+ trailing `N`s, distributed across 5 barcodes with reproducibly short 3′ contigs. The N-count distribution is sharply bimodal (0–9 `N`s vs. 10–100+), enabling a clean threshold cut.

**Command:**
```bash
./filter_by_ambiguity.sh
# Prompted inputs:
#   FASTA to filter (Step 7 backfilled output)
#   Max N count per sequence [default 9]
#   Max N fraction [0-1, default 0.005]
```

**Key parameters:**
- `max_N = 9` — Library 010 distribution is sharply bimodal (clean cluster 0–9 N's vs. problematic 10–100+ N's). Threshold of 9 cleanly separates them.
- `max_N_fraction = 0.005` — secondary safety check (~13 N's for a 2530 bp sequence).

**Input:**
- `all_Library00X_..._blastfilt_aligned_exons_backfilled.fasta` — output of Step 7

**Output:**
- `all_Library00X_..._blastfilt_aligned_exons_backfilled_Nfilt.fasta` — N-purged DNA, input to Step 8
- `all_Library00X_..._blastfilt_aligned_exons_backfilled_Nfilt_log.tsv` — per-query decision log

**Validation (Library 010):** 1516 → 1388 (–128 sequences dropped); resolves the trailing-gap artefact on BEA canonical AA references in Step 8 output. Runtime <1 s.

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

**Input:** `all_Library00X_Phased_haplotypes_filtered_min3250_max4000_blastfilt_aligned_exons_backfilled_Nfilt.fasta`

**Key output:**
- `all_Library00X_Phased_haplotypes_filtered_min3250_max4000_blastfilt_aligned_exons_backfilled_Nfilt_frame1_AA_filtered_aligned.fasta`

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

**Scripts:** `find_allele_plateau.py` (calibration), `define_SRK_alleles_from_distance.py` (clustering), `SRK_AA_mutation_heatmap.py` (variation heatmaps)

> Run Step 10a-pre once per dataset to choose `N_ALLELES`, then Step 10a, then Step 10b.

#### Step 10a-pre — Calibrate `N_ALLELES` (Kneedle elbow detection)

Whenever the input dataset changes (new library, re-run after a filter change), run the calibration utility to mathematically identify the elbow of the sensitivity curve.

**Command:**
```bash
python find_allele_plateau.py
```

**Methods reported (three independent estimators):**
1. **Kneedle** (max perpendicular distance from the chord connecting curve endpoints) — the canonical "elbow" detection method (Satopaa et al. 2011). This is the recommended estimator and matches the visually-identified plateau in well-behaved curves.
2. **Slope-magnitude minimum** (smoothed `|dN/dt|`) — finds the flattest region within a constrained allele-count range. May lock onto secondary plateaus lower on the curve.
3. **Longest persistent run** (contiguous threshold steps where `N` stays within ±3 alleles of a candidate centre) — robust measure of plateau width.

The script also prints a per-threshold table around the Kneedle elbow so the curve shape can be sanity-checked numerically.

**Recommended action:** use the Kneedle value unless the per-threshold table shows it landing on a knife-edge transition rather than a stable level. Inspect `SRK_protein_distance_analysis.pdf` (Step 10a output) for visual confirmation.

**Outputs:** stdout report only — no files written. The chosen value of `N_ALLELES` is then manually entered at the top of `define_SRK_alleles_from_distance.py`.

#### Step 10a — Allele clustering

**Command:**
```bash
python define_SRK_alleles_from_distance.py
```

**Key parameters (edit at top of script):**
| Parameter | Default | Notes |
|-----------|---------|-------|
| `INPUT_FASTA` | `SRK_functional_proteins_aligned.fasta` | Aligned proteins from Step 9 |
| `N_ALLELES` | `58` | Option A: fix number of alleles directly; set from Step 10a-pre Kneedle output |
| `DIST_THRESHOLD` | `0.01` (1%) | Option B: fix p-distance cutoff; used only when `N_ALLELES = None` |
| `DOMAIN_REGION` | `(31, 430)` | S-domain columns (1-based); adjust if species-specific annotation available |

> **Calibration history:** for the current dataset (Libraries 001–010, 380 individuals, 376 functional proteins after Step 4b + Step 7b filtering), `find_allele_plateau.py` returned a Kneedle elbow at N = 58 (implied threshold ≈ 0.0055), confirmed visually on the sensitivity curve. Earlier (Libraries 001–009, 302 functional proteins) the value was N ≈ 63 (implied threshold ≈ 0.005). The shift reflects (i) the addition of Library 010 and (ii) the cleaner protein set produced by Steps 4b/7b.

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

### Step 12c — Data Quality Evaluation (end of Phase 2)

A formal data-quality gate at the end of Phase 2 that answers three sequential questions before downstream population-genetic analyses (Phase 3) begin:

1. **Is there a library effect impacting data interpretation?** (technical bias check)
2. **Which samples failed sequencing/assembly and need to be re-sequenced?** (lab deliverable)
3. **Which samples have non-functional SRK proteins (full to partial) and may have escaped self-incompatibility?** (phenotyping deliverable)

Outputs include per-Bottleneck-Lineage and per-EO stacked-bar plots showing the proportion of each outcome category, and two CSVs that can be handed directly to the wet-lab team.

#### Step 12c.i — Library effect tests

**Script:** `test_library_effect.py`

**Tests:**
- **Test 1 — Library × `SI_functional_status` (global χ²)**: tests whether Functional / Partial_translation_failure / Complete_loss proportions are independent of library.
- **Test 2 — Library × `Dominant_failure_mode` (global χ²)**: tests whether `premature_stop` / `ambiguous_aa` / `mixed` failure-mode classification is independent of library.
- **Tests 3 + 4 — Within-EO Library 009 / 010 vs other libraries (Fisher's exact)**: for each focus EO with multi-library sampling, tests whether the focus library's AAAA proportion (Test 3) or Complete_loss proportion (Test 4) differs from the pooled "other libraries" in the same EO.

**Outputs:**
- `SRK_library_effect_tests.tsv` — chi-square statistics, df, p-values, standardised residuals
- `SRK_library_effect_summary.pdf` — 4-panel diagnostic (residual heatmaps + within-EO bar plots)

**Interpretation rule:** a *significant* `Library × Dominant_failure_mode` association is expected if libraries differ in pipeline processing — e.g., Library 010 (Step 7b applied) should show enrichment in `premature_stop` (real LoF) and depletion in `mixed`/`ambiguous_aa` (data-quality noise). This is a *favourable* library effect — the cleanest library produces the most biologically interpretable failure pattern. A significant `Library × SI_functional_status` association, by contrast, would indicate that *the rate of biological failure itself depends on library*, which would be a confound. The within-EO tests are the cleanest test for this.

#### Step 12c.ii — Sample exclusion audit

**Script:** `audit_sample_exclusions.py` (documented in detail in *Sample Exclusion Audit* section below).

Produces `SRK_sample_exclusion_audit.tsv` consumed by Step 12c.iii.

#### Step 12c.iii — Integrate into outcome categories, CSV deliverables, and proportion plots

**Script:** `evaluate_data_quality.py`

Reads the audit and library-effect outputs and categorises every metadata sample into one of five outcome categories (mutually exclusive; samples with `Ingroup = 0` are flagged separately as Outgroup and not shown in the plots):

| Outcome category | Criterion | Lab action |
|------------------|-----------|------------|
| **Functional** | `SI_functional_status = Functional` (included) | None — sample contributes to all downstream analyses |
| **Partial_translation_failure** | `SI_functional_status = Partial_translation_failure` (included) | None — sample contributes. **NEUTRAL framing: this is most likely a technical artefact (chimeric Canu haplotypes inflating the failure denominator), not biological SI loss. The current dataset provides no tangible evidence supporting biological SI loss in this category — make no inference about SI status from this label alone.** |
| **SI_escape_candidate** | `SI_functional_status = Complete_loss` AND `Dominant_failure_mode = premature_stop` (excluded) | Phenotype via controlled selfing test (no-pollen-deposition vs self-pollen vs cross-pollen seed-set comparison) |
| **Re_PCR** | Stage 3 / Step 4 length / Step 4b BLAST paralog / Step 7b / dropped_abundance_filter / no_functional_proteins with ambiguous_aa or mixed failures | **Re-amplify SRK from the existing DNA stock.** The DNA is fine; the PCR product is the problem (no assembly produced, fragmented amplicon, paralog amplification, short product yielding N-rich Canu contigs, low yield, or dirty product translating to X residues) |
| **Re_DNA_extraction** | Step 9 too_many_alleles (>4 distinct functional proteins per individual) | **Re-isolate tissue from a single plant + re-extract DNA before any further PCR.** The DNA stock itself is contaminated (mixed sample or barcode bleed-through); re-PCR with the same DNA will preserve the contamination |

**Why the Re_PCR / Re_DNA_extraction split matters operationally.** These are *different* lab tasks with different costs and turnaround times. Re_PCR samples have viable DNA on hand — a same-day operation that just re-amplifies. Re_DNA_extraction samples require fresh tissue (potentially field collection), then DNA extraction, then PCR — a much longer pipeline. Surfacing them as separate CSVs lets the wet-lab team triage and schedule appropriately.

**Step 4b paralog reasoning.** All 24 of the Library 010 paralog samples are routed to Re_PCR because PCR product competition is stochastic — a fresh PCR may amplify the canonical SRK locus rather than the paralog. If the *same* paralog signature persists across re-PCR replicates, that signals a need for primer re-design, which is a different operation than re-PCR; in that case re-classify those samples as `Excluded_paralog` in the audit and exclude from the lab CSV.

**Outputs:**

- `Tables/SRK_data_quality_categories.tsv` — master per-sample categorisation with all upstream diagnostic columns preserved.
- **`Tables/SRK_samples_redo.csv`** — single merged lab deliverable for samples needing follow-up. First column `Lab_action` discriminates **Re-PCR** (existing DNA OK, re-amplify) from **Re-DNA-extraction** (DNA contaminated, re-isolate tissue). Rows sorted by `Lab_action → EO → Sample_ID` for field-collection grouping. Other columns: `Sample_ID`, `Library`, `Barcode`, `EO`, `EO_normalised`, `BL_inferred`, `n_raw_haps`, `Proteins_in_final_data` (informative for Re-DNA-extraction rows where it shows the contaminated extraction's distinct-protein count), `Exclusion_stage`, `Recommended_action` (stage-specific lab instruction).
- **`Tables/SRK_SI_escape_candidates.csv`** — phenotyping deliverable: only SI-escape candidates, with the recommended controlled-selfing test instruction.

**Modeling deliverables also mirrored into `Tables/`** (produced upstream of Step 12c; bundled here as a single download location for downstream collaborators):

- **`Tables/SRK_individual_allele_genotypes.tsv`** — wide allele×individual count matrix (335 ingroup individuals × 49 alleles, integer copy counts), produced by Step 11. Ready-to-use design matrix for modeling.
- **`Tables/SRK_individual_zygosity.tsv`** — per-individual genotype summary produced by Step 12. Same information as the wide matrix above, in packed-string form.
- **`Tables/SRK_synonymy_groups.csv`** — allele → synonymy group mapping produced by Step 22b. Lets modelers collapse putatively identical alleles into one functional unit.
- **`Tables/sampling_metadata.csv`** — canonical sample registry (`SampleID`, `Pop`, `EO_w_sub`, `Ingroup`, `OccurrenceID`); the population key joining the genotyping tables to `Tables/EO_group_BL_summary.csv` (EO → BL → drift index).
- `figures/SRK_data_quality_per_BL.{pdf,png}` — stacked-bar chart showing the proportion of each outcome category for each Bottleneck Lineage (BL1–BL5). X-tick labels include sample count `N`; bars are stacked Functional → Partial_translation_failure → SI_escape_candidate → Re_PCR → Re_DNA_extraction with absolute counts annotated inside each non-trivial slice.
- `figures/SRK_data_quality_per_EO.{pdf,png}` — same chart for the 6 focus EOs (N ≥ 5). X-tick labels are coloured by parent BL using the locked Set1 palette so each EO is visually placed within its lineage.

**Why this matters operationally:** the per-BL and per-EO proportion plots make it trivial to spot populations where Re-sequence + SI_escape_candidate categories are over-represented — the populations where conservation action is most urgent and where additional lab effort is best invested.

**When to run (revised 2026-06-13):** the three sub-scripts now have hard dependencies on Phase 3 outputs that did not exist when this step was first written. The correct execution order is:

1. **Step 12c.ii** (`audit_sample_exclusions.py`) — run **after Step 13** (BL integration) so the per-sample BL column is populated. Strictly the script tolerates a missing BL file, but the resulting audit is more useful with BL columns present.
2. **Step 12c.iii** (`evaluate_data_quality.py`) — run **after Step 13** and after 12c.ii. The per-BL / per-EO stacked-bar plots require BL assignments; the script will run without them but with degraded output.
3. **Step 12c.i** (`test_library_effect.py`) — run **after Step 19** (Individual GFS), since one of the four tests reads `SRK_individual_GFS.tsv`. This sub-step is therefore the last 12c sub-script to execute, not the first.

The full chain `Step 13 → 12c.ii → 12c.iii → ... → Step 19 → 12c.i` is library-agnostic, auto-detects per-library log files, and the QC sub-steps themselves take < 5 s each.

---

## Phase 3: Data Analyses

> All Phase 3 scripts read from the shared outputs of Phase 2. Run from the same working directory.

Phase 3 evaluates the evolutionary status and conservation implications of the SI system by integrating population genetic theory with the genotype data generated in Phases 1 and 2. Before any inferential question can be addressed at the right evolutionary scale, individuals must be assigned to the demographic units that share — and have shared — a gene-flow history. Element Occurrences (EOs) are the unit of management, but they are nested within larger evolutionary units: groups of EOs that descend from a common ancestral colonisation event and have evolved independently since fragmentation severed connectivity between them. These independent demographic units are termed **bottleneck lineages (BL1–BL5)** and are derived in a sibling spatial-clustering project (`LEPA_EO_spatial_clustering`) from a 500 m pollinator-dispersal threshold, hull-to-hull connectivity analysis, and Ward's D2 hierarchical clustering of group centroids (silhouette-optimal k = 5). Phase 3 therefore opens with **Step 13 — Element Occurrence and Bottleneck Lineage Integration**, which joins each individual to its EO, geographic group, BL, and drift index. All subsequent steps stratify in parallel by **EO sorted within BL** (the management view) and by **BL** (the evolutionary view), the latter recovering the small-locality samples that EO-level statistics cannot use. Three interconnected scientific questions then structure the analyses.

**Q1 — What are the population genetic parameters of the SRK system?**
Before interpreting SI diversity in an evolutionary framework, it is necessary to characterise each population's fundamental genetic state: allele richness, effective allele number, heterozygosity, and mean alleles per individual (Step 14). These descriptors provide the empirical foundation for all downstream comparisons.

**Q2 — What is the total allele richness of the SI system at the species level?**
Under negative frequency-dependent selection (NFDS), rare S-alleles confer a reproductive advantage because pollen bearing a rare allele can fertilise a larger proportion of compatible stigmas. This self-reinforcing dynamic drives the accumulation of allelic diversity and, at evolutionary equilibrium, pushes allele frequencies towards equality. The *species allele richness* — estimated from accumulation curves and asymptotic richness estimators (Michaelis-Menten, Chao1, iNEXT) applied to the full dataset (Step 15) — approximates this NFDS equilibrium. It represents the total pool of functionally distinct S-alleles maintained across the species and serves as the reference optimum against which population-level deficits are measured. This question is foundational: without a species-level baseline, it is impossible to determine how much allelic diversity individual populations have lost relative to the evolutionary expectation.

**Q3 — How has genetic drift eroded the SI system at the population level?**
In small or isolated populations, genetic drift counteracts NFDS by stochastically eliminating rare alleles, skewing allele frequencies away from the equal-frequency expectation, and increasing the prevalence of homozygous genotypes. Drift erodes the SI system at two hierarchical levels, each corresponding to a conservation tipping point. **Tipping Point 1 (TP1)** evaluates whether the mating pool is still functioning at the population level — i.e. whether enough cross-compatible plant pairs remain to sustain reproduction without external introductions. **Tipping Point 2 (TP2)** is breached when the distribution of the remaining alleles across individuals has degraded to the point that managed crossing within the population cannot recover reproductive fitness without external introductions. Phase 3 quantifies both levels through four complementary analyses:

*TP1 — Allele erosion, evenness, and mating-pool functionality (Steps 15–18):*

- **Allele accumulation curves** (Step 15) reveal how many S-alleles each population retains relative to the species optimum, measuring the depth of drift-driven allele loss.
- **Allele frequency analysis** (Step 16) tests whether population-level frequencies deviate from the equal-frequency NFDS expectation via χ² goodness-of-fit, using the species richness estimate from Step 15 as the optimum.
- **TP1 — mating-pool functionality** (Step 17) takes the depletion established by Steps 15–16 as given and asks the complementary question: *is the depleted population still mating?* Per group (EO and BL), the script computes evenness J (Shannon H / ln k) and tetraploid sporophytic-SI compatibility P_compat under increasing SI-leakage rates (L ∈ {0, 0.10, 0.25, 0.50}). Four diagnostic panels (EO and BL × strict L=0 and leaky L=0.25) place each group into one of four conservation-action quadrants: MONITOR, AUGMENT-evenness, AUGMENT URGENTLY, or BIOBANK + RESTORE. The empirical leakage rate L̂ ≈ 0.18 inferred from observed AAAA proportions reflects C3 Stage 5 directly.
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
- `SRK_individual_zygosity.tsv` — canonical post-QC individual list (currently 335 rows). Used to filter the sampling registry to the individuals that actually completed Phase 2.
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

### Step 13b — Sampling Overview Map

**Script:** `SRK_sampling_map.R`

**Command:**
```bash
Rscript SRK_sampling_map.R
```

**Purpose.** Single-panel geographic map showing every catalogued population of the species against its parent BL, with a clear visual distinction between populations successfully sampled by the SRK pipeline (filled triangles) and not-yet-sampled populations (filled circles). The figure is the recommended *first* visual in any presentation or report: it sets the geographic baseline before any per-locus analysis and makes the gap between current sampling and the full species range immediately obvious.

**Inputs:**
- `Tables/SRK_individual_allele_genotypes.tsv` — canonical successful-sample roster (= 335 individuals)
- `Tables/sampling_metadata.csv` — `Pop` and `EO_w_sub` lookup
- `SRK_individual_BL_assignments.tsv` — from Step 13 (provides per-sample EO → BL resolution)
- `EO_location_groups.csv` — coordinates for every catalogued population, mirrored from the sibling `LEPA_EO_spatial_clustering` project (path resolution: NSF data dir first, then `Tables/` fallback)
- `Tables/EO_group_BL_summary.csv` — provides parent BL for *every* catalogued EO, including ones that are catalogued but not yet sampled (so unsampled populations still get their BL colour)
- `srk_bl_constants.R` — locked Set1 palette + `BL_ORDER`

**Visual conventions:**
- **Shape**: filled triangle = sampled; filled circle = not yet sampled.
- **Colour**: parent BL (Set1 palette, matches every BL-stratified figure in the project).
- **Size**: habitat area (ha) — `sqrt`-scaled so small populations remain visible.
- **Labels**: for each sampled EO, the EO code (with `n =` sample count) is placed on the **largest** location within that EO (the convention for compound germplasm sub-codes that cannot be pinned to a specific sub-location). Focus EOs (EO18, EO25, EO27, EO67, EO70, EO76) carry bold white-fill labels; other sampled EOs carry smaller plain labels. Unsampled populations are left unlabelled.
- **Landscape backdrop**: same DEM topo + Snake River + reference cities as the sibling project's `EO_BL_geographic_context_map.png`, so the two figures cross-reference visually.

**Outputs:**

| File | Content |
|------|---------|
| `figures/SRK_sampling_map.png` / `SRK_sampling_map.pdf` | Captioned version for the report / pipeline doc |
| `figures/SRK_sampling_map_presentation.png` / `SRK_sampling_map_presentation.pdf` | Title / subtitle / caption stripped, for slide decks |

**Headline numbers (current dataset, 2026-05-21).** **325 BL-resolved samples across 17 populations in 16 EOs** (BL1: 5 samples / 4 EOs / 4 Pops · BL2: 59 / 2 / 2 · BL3: 80 / 4 / 4 · BL4: 89 / 2 / 2 · BL5: 92 / 4 / 5) + **10 samples / 9 germplasm sub-codes in the JAR region** (no resolved coordinates, not plotted) = **335 total successful ingroup samples**.

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

The script restricts analyses to BL-assigned individuals (325 of 335; the 10 unresolved germplasm sub-codes are dropped because they cannot be placed into an EO/BL).

**Key metrics computed per group (EO and BL):**
- N_alleles — allele richness
- Effective_alleles_Ne = 1/Σpᵢ² (based on `colSums` copy counts)
- Prop_heterozygous / Prop_homozygous
- Mean_alleles per individual
- Top_alleles — three most abundant alleles with copy counts

**BL palette and ordering — single source of truth.** Every Phase 3/4 script that orders or colours BLs sources/imports the mirrored modules `srk_bl_constants.R` and `srk_bl_constants.py`. Both derive their orderings at load time from CSVs mirrored in `Tables/` (`EO_BL_summary.csv`, `EO_group_BL_summary.csv`), so adding new sampling only requires (i) re-running `LEPA_EO_spatial_clustering`, (ii) refreshing those CSVs — the entire pipeline reorders automatically. The module exposes:

- `BL_COLORS` — RColorBrewer Set1 palette mapped to BL by cluster-index, matching the rendered figures in `LEPA_EO_spatial_clustering`: `BL1 = #984EA3` (purple), `BL2 = #377EB8` (blue), `BL3 = #E41A1C` (red), `BL4 = #FF7F00` (orange), `BL5 = #4DAF4A` (green). Unassigned individuals (if retained) use `#999999`.
- `BL_ORDER` — BLs sorted high-to-low by **total habitat area (ha)**, with **within-BL connectivity** as the secondary tie-break, then BL name as a deterministic final tie-break. Connectivity = `n_locations − n_groups`. For the current dataset: **BL4, BL5, BL3, BL1, BL2**. Area is the primary *Ne* proxy (carrying capacity → drift floor); connectivity is the secondary stratification for BLs of similar size. Used for **bars, faceted panels, horizontal-bar axes, and tables** wherever BL appears as a category axis and the ordering itself tells a story.
- `BL_ORDER_NUMERIC` — alphanumeric BL1→BL5. Used for **scatter-plot legends** (TP1, TP2, accumulation-curve legends) where the axes are not BL and a readable legend matters more than the inferential order.
- `get_eo_order_within_bl(eo_codes)` — EOs ordered first by parent `BL_ORDER`, then by ascending mean Drift_index within BL, with EO name as deterministic tie-break.

Composite EO entries in `EO_group_BL_summary.csv` of the form `"EO118; EO76"` (one geographic group spanning two EOs) are split inside the helper; each EO inherits the group's mean Drift_index. The R and Python implementations apply the same composite-split + tie-break rules and therefore produce identical orderings.

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

**Split-filter design.** This is the only step where the sample set differs by analysis level: the **species-level curve and richness estimators use ALL ingroup individuals** (current dataset 2026-05-11: 335 individuals; species pool: 49 observed / MM = 59 / Chao1 = 61 / consensus = 60) so the baseline includes alleles unique to the 10 individuals whose germplasm sub-codes are not yet resolved to a BL — those alleles still count toward the species pool. The **EO-level and BL-level curves use only the BL-assigned individuals** (325 of 335 in current dataset; was 262 of 272 prior to Library 010), where geographic provenance is required.

**Key parameters:**
- Permutations: 1000 (default)
- EOs with <5 individuals excluded from EO-level curves; all 5 BLs always included
- EO curves are sorted by parent BL and colored using the locked Set1 BL palette (Step 14)
- Three richness estimators: Michaelis-Menten (MM), Chao1, iNEXT (optional — `install.packages("iNEXT")`)

**Outputs:**
- `SRK_allele_accumulation_curves.pdf` — species + per-EO + per-BL curves on separate pages, each with MM/Chao1/iNEXT asymptote reference lines.
- `figures/SRK_allele_accumulation_species.png` — Two-panel species figure (stacked bar on the left + rarefaction curve on the right, sharing the bar's y-axis). The bar decomposes the species MM ceiling into the observed alleles (dark blue) and the predicted-undetected alleles (light blue, = MM − observed), using the **same colour key as the BL/EO drift-erosion bars** so a reader can read the species, BL, and EO panels together. MM/Chao1/iNEXT asymptote lines are drawn on the curve panel and labelled on its right margin; a combined legend at bottom-right covers both the bar colour key and the asymptote line styles.
- `figures/presentation/SRK_allele_accumulation_species_observed.png` / `_predicted.png` — **NEW**. Slide build-up frames of the figure above. Frame 1 (`_observed`) renders only the bar's observed segment with its y-axis (no curve, no MM line). Frame 2 (`_predicted`) renders the full bar (observed + predicted undetected) + curve + MM dashed line, with Chao1/iNEXT suppressed. Same canvas size (12 × 6.5 in) across both frames so they animate cleanly in a deck.
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

**Same split-filter design as Step 15:** the species-level χ² test uses ALL 335 ingroup individuals; EO and BL χ² tests use only the 325 BL-assigned individuals.

**Key analysis:**
- χ² goodness-of-fit test vs. equal-frequency NFDS expectation
- Three levels: Species (1 row), EOs with N ≥ 5 (sorted by parent BL), BLs (all 5 always tested)
- All frequency plots share x-axis scaled to species optimum (MM estimate); missing-allele zone shown
- Each EO/BL frequency page is colored by its parent BL (Set1 palette) so the BL identity is visually obvious without checking the legend

**Outputs:**
- `SRK_chisq_species_population.tsv` — χ² statistics (X2, df, p-value) per level. Columns: `N_individuals`, `N_alleles`, `X2`, `df`, `p_value`, `Level` (Species/EO/BL), `Population`, **`BL`** (new). Rows ordered Species → EOs (sorted by BL) → BLs.
- `SRK_chisq_species_population_frequency_plots.pdf` — 12 pages: 1 species + 6 EOs (sorted by BL) + 5 BLs. Bars colored by parent BL; species page uses neutral grey.

---

### Step 17 — TP1 Mating-pool Functionality

> Requires outputs from Steps 13, 14, and the zygosity TSV from Step 12. Run after all have completed.

**Scripts (run in order):**
1. `SRK_TP1_compatibility_metrics.py` — computes per-EO + per-BL metrics (including bootstrap 95 % CIs on P_compat **and the Depletion Index** `DI = 1 − k / k_species`, both observed and MM-predicted flavours) and writes `Tables/SRK_EO_allele_richness.tsv`
2. `SRK_P_compat_traffic_light.R` — stakeholder-facing red/amber/green priority figure (full + blank variants); answers the question *is random mating still viable in this population?* at a glance
3. `SRK_depletion_ranking.R` — **conservation-ranking figure** placing each group on P_compat × DI axes with quadrants labelled by the three breeding strategies in `Tables/SRK_breeding_strategies.csv` (HEALTHY / INFORMED BREEDING (frequency skew) / INFORMED BREEDING + ALLELE INJECTION (preventive / urgent)). Two panels (`_observed`, `_predicted`) × three variants (`_blank`, `_EOs`, `_all`) = 6 figures for layered presentation builds

**Command:**
```bash
python3 SRK_TP1_compatibility_metrics.py
Rscript SRK_P_compat_traffic_light.R
Rscript SRK_depletion_ranking.R
```

**Inputs:**
- `Tables/SRK_individual_allele_genotypes.tsv` — Step 11 allele copy-count matrix
- `Tables/sampling_metadata.csv` — EO assignment per individual
- `SRK_individual_BL_assignments.tsv` — Step 13 BL labels
- `Tables/SRK_individual_zygosity.tsv` — Step 12 inferred tetraploid genotype classes (the metrics script mirrors the same inference rule directly from the raw matrix; the zygosity TSV is the validated reference)

**Why TP1 was reframed (2026-05-22):** Depletion against the species optimum is already shown by the Step 16 erosion barplots; the previous TP1 (`prop_optimum × Ne/N_alleles`) re-told that story. The reframed TP1 asks the complementary question — *given that a population is depleted, is its mating pool still functioning, and what conservation intervention does the data support?* Axes are therefore evenness J (frequency-shape diagnostic) × P_compat (direct demographic compatibility under tetraploid sporophytic SI).

**Two analysis levels on separate panels (four figures total):**
- **EO** — circles, six EOs with N ≥ 15 (EO18, EO25, EO27, EO67, EO70, EO76), coloured by parent BL
- **BL** — triangles, all five lineages (BL1–BL5), pooling every BL-assigned ingroup individual

**Key metrics per group:**
- `k_observed` — distinct S-alleles present after tetraploid inference
- `k_rarefied30_mean ± sd` — subsample 30 individuals, 1000 perms (sample-size-corrected richness; encoded as point size in the figures)
- `k_species` — species-level optimum from Step 15 (MM consensus = 59 in the current dataset)
- `k_predicted_MM` — per-group MM asymptote from Step 15 (upper bound on richness achievable from existing sampling without further intervention)
- **`DI_observed`** = `1 − k_rarefied30 / k_species` — **Depletion Index** using sample-size-corrected richness. 0 = at species equilibrium; 1 = no alleles retained. The conservative read because it does not extrapolate beyond observed sampling.
- **`DI_predicted`** = `1 − k_predicted_MM / k_species` — Depletion Index using the MM asymptote per group. Gaps between `DI_observed` and `DI_predicted` flag populations where more sampling would likely recover additional alleles (e.g. EO27 obs 0.71 / pred 0.47; BL4 obs 0.71 / pred 0.37) versus those already near-asymptote (e.g. EO70 obs 0.92 / pred 0.88).
- `evenness_J` = Shannon H / ln(k) — frequency-shape diagnostic; J = 1 at NFDS equilibrium
- `prop_AAAA` — fraction of individuals inferred AAAA (tetraploid homozygous at SRK)
- `L_hat_from_AAAA` = prop_AAAA / 3.5 — empirical SI-leakage estimate (population-averaged tetrasomic correction; upper bound because amplicon under-recovery also produces apparent AAAA)
- `P_compat_L0/L0.10/L0.25/L0.50` — the **compatible-pair fraction**: the fraction of randomly drawn plant pairs that are cross-compatible under tetraploid sporophytic SI with co-dominance, at increasing SI-leakage levels. A value of 0.40 means roughly 40 % of random pairs in the group can produce seed. Exact multinomial formula: `Σ_{a,b} p_a p_b (1 - p_a - p_b·I[a≠b])^4 + L·(1 - that)`.
- `P_compat_L{0,0.10,0.25,0.50}_lo/_hi` — bootstrap 95 % confidence interval (2.5th / 97.5th percentile from 1000 resamples of individuals with replacement). Captures sampling uncertainty in the group's allele-frequency estimate; tight CIs at large N (e.g. EO70 N = 56: CI ±0.02), wide CIs at small N (e.g. BL1 N = 5: CI ±0.04 on a 0.07 point estimate).

**Inheritance-mode caveat:** the P_compat formula assumes **tetrasomic** inheritance with random pairing across all four chromosomes. If LEPA's SRK locus shows **disomic** inheritance (two homoeologous subgenomes segregating independently), the system behaves more like two superimposed diploid SI systems and is generally more permissive — current P_compat is then a conservative lower bound. Documented in the script docstring.

**TP1 thresholds:**

| Parameter | Default | Conservation reading |
|-----------|---------|----------------------|
| `TP1_J` | 0.80 | evenness floor (frequencies skewed from NFDS expectation) |
| `TP1_P_COMPAT` | 0.40 | mating-pool floor (fewer than 40% of random pairs compatible) |

**Quadrants map to conservation actions:**
- top-right — **MONITOR** + augment for long-term sustainability
- top-left — **AUGMENT** to restore evenness
- bottom-right — **AUGMENT URGENTLY** (too few compatible mates)
- bottom-left — **BIOBANK + RESTORE** (mating pool collapsed)

**Depletion-ranking quadrants** (`SRK_depletion_ranking.R`; axes P_compat × DI; thresholds DI = 0.50, P_compat = 0.40; quadrant labels drawn from `Tables/SRK_breeding_strategies.csv`):

| DI | P_compat | Quadrant label / intervention |
|---|---|---|
| < 0.50 | ≥ 0.40 | **HEALTHY** (random pollination) |
| ≥ 0.50 | ≥ 0.40 | **INFORMED BREEDING + ALLELE INJECTION** (preventive) |
| < 0.50 | < 0.40 | **INFORMED BREEDING** (frequency skew) |
| ≥ 0.50 | < 0.40 | **INFORMED BREEDING + ALLELE INJECTION** (urgent) |

In the current dataset every group except EO27 falls in one of the two ALLELE INJECTION quadrants. EO27 sits in INFORMED BREEDING (frequency skew) on the `_predicted` panel and on the threshold of `_observed`.

**Outputs:**
- `Tables/SRK_EO_allele_richness.tsv` — 5 BL rows + 6 EO rows × 31 columns (identifiers + depletion + genotype mix + richness + evenness + leakage ladder with bootstrap CIs)
- `Tables/SRK_breeding_strategies.csv` — 3-row reference table (Random pollination / Informed breeding / Informed breeding + allele injection) describing how parents are paired and what each strategy produces; supplies the vocabulary used in the depletion-ranking quadrants
(P_compat traffic-light + DI ranking figures listed below — no separate J × P_compat scatter is produced.)
- `figures/SRK_P_compat_traffic_light_EO.{png,pdf}` — **NEW (2026-05-22)**, stakeholder-facing traffic-light figure: six focal EOs ranked by strict-SI P_compat, coloured red / amber / green against the 0.20 and 0.40 thresholds, with bootstrap 95 % CIs as horizontal error bars and BL colour squares on the left edge. The dedicated random-mating-viability diagnostic; foregrounded in the conservation report and any stakeholder presentation
- `figures/SRK_P_compat_traffic_light_EO_blank.{png,pdf}` — companion blank variant (zones + EO labels + BL strip + legend, no data drawn); use as the predictions slide before revealing the full figure
- `figures/SRK_depletion_ranking_observed_{blank,EOs,all}.{png,pdf}` — P_compat × DI conservation-ranking figure using `DI_observed` (sample-size-corrected k). Three variants (zones only / focal EOs only / BLs + EOs together) for layered presentation builds.
- `figures/SRK_depletion_ranking_predicted_{blank,EOs,all}.{png,pdf}` — same figure using `DI_predicted` (MM asymptote per group); represents the upper bound on achievable richness without further intervention. Differences from the `_observed` panel highlight populations where additional sampling would likely recover alleles.

**Headline result (current dataset, 2026-05-22):**
- Strict SI: **EO70 (BL2)** sits alone in BIOBANK + RESTORE (P_compat = 0.08, J = 0.67, k_rare30 = 4.9). EO67 is also BIOBANK (P_compat = 0.23). EO76, EO18 sit in AUGMENT URGENTLY. EO25, EO27 in MONITOR. BL1 + BL2 in BIOBANK; BL3 borderline; BL4 + BL5 in MONITOR.
- Leaky SI (L = 0.25, ≈ empirical L̂): EO76, EO18, BL3 rescued into MONITOR; **EO70 and BL2 stay in BIOBANK under both views** — the cleanest signal that SI leakage cannot save these populations and seed banking is the priority.

**Connection to the C3 cascade:** the leaky panel visualises Stage 5 directly — populations that look catastrophic under strict SI are demographically alive *only because* SI breakdown allows leaky selfing. The cost is homozygosity accumulation (the geno_mix column shows 55–73% AAAA across the focal EOs), which feeds back into Stage 2 of C3.

**Legacy:** the prior implementation `SRK_TP1_tipping_point.R` is archived under `archive/`.

---

### Step 18 — Allele Composition Comparison

**Scripts:**
1. `SRK_allele_sharing_EOs.py` — UpSet plots + pairwise heatmaps for EOs and BLs
2. `SRK_allele_eulerr_BLs.R` — area-proportional 5-set Euler diagram of BL allele sets (stakeholder-facing companion to the BL UpSet)

**Command:**
```bash
python3 SRK_allele_sharing_EOs.py
Rscript SRK_allele_eulerr_BLs.R
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
- `figures/SRK_allele_eulerr_BLs.{png,pdf}` — **NEW (2026-05-22)**, area-proportional 5-set Euler diagram of BL allele sets. Same data as the BL UpSet, rendered as familiar Venn-style ellipses for stakeholder presentations. R companion script (`SRK_allele_eulerr_BLs.R`); regenerate together with the Python UpSet whenever new data are added. Requires the `eulerr` R package.

**Headline result (current dataset, 2026-05-11):** 26 of 49 alleles (53 %) are private to a single BL; only Allele_050 and Allele_051 are shared across all 5 BLs. The near-disjoint BL allele sets are the **direct test of independent bottlenecks** — a single shared species-level bottleneck would predict overlapping losses, not the observed lineage-private alleles. BL4 holds 10 private alleles (37 % of its 27-allele complement), confirming its role as the species' diversity reservoir.

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

**AAAA allele identity panels:** unpack the AAAA bar by showing which alleles each AAAA individual carries. *Pan-BL* alleles are present in AAAA individuals across every BL; *pan-EO* alleles are present in AAAA individuals in ≥80% of focus EOs (relaxed threshold accommodates EO18, the smallest focus EO with only 4 AAAA individuals carrying a different drift signature). Allele_050 and Allele_051 (Synonymy group 1, HV-identical, likely the same SI specificity) appear as pan-BL across all 5 BLs and pan-EO in 5 of 6 focus EOs — confirming that independent bottlenecks have all converged on the same fixed SI specificity.

**Headline:** BL3 has the lowest reproductive effort support (34%, mean GFS = 0.196); BL2 has the highest Synonymy group 1 saturation among AAAA individuals (69%); EO18 is the only focus EO without Synonymy group 1 alleles in its AAAA pool (4 AAAA individuals, 4 distinct non-Synonymy-group-1 alleles).

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
| `SRK_GFS_AAAA_allele_composition.pdf` | 2-page PDF: EO panel + BL panel — allele identity of AAAA individuals; Synonymy group 1 alleles (Allele_050, Allele_051) highlighted |
| `SRK_BL_reproductive_effort_summary.tsv` | **NEW** — BL-level reproductive effort summary (n, n_supporting, prop_supporting, prop_AAAA, mean_GFS) |
| `figures/SRK_GFS_reproductive_effort_EO.png` | EO-level proportional bar (PNG) |
| `figures/SRK_GFS_reproductive_effort_BL.png` | **NEW** — BL-level proportional bar (PNG) |
| `figures/SRK_GFS_AAAA_allele_composition_EO.png` | EO-level AAAA allele identity (PNG) |
| `figures/SRK_GFS_AAAA_allele_composition_BL.png` | **NEW** — BL-level AAAA allele identity (PNG) |

---

---

## Phase 4: Per-individual SI status, null-aware genotypes, and forward simulation

> Phase 4 establishes the per-individual self-incompatibility status (full SI, partial SI, self-compatible, insufficient_data) using the per-haplotype OK / REMOVED calls from Step 7 and the AA-distance broken-allele assignments from Step 25a. This information is *required* before Phase 5 (cross design) — a cross plan is only meaningful when the parents have a verified functioning SI machinery.

Phase 4 contains Steps 25 (per-individual SI / pSI / SC), 26 (null-aware genotype rebuild + cascade re-runs of Steps 14b / 17b / 19b+20b), 27 (forward-time tetraploid inheritance simulator), and 28 (donor ranking). All scripts read from `Tables/Phase4/` and write to `Tables/Phase4/` + `figures/Phase4/`.

(Step 22 + Step 23 — the experimental cross design — appear in **Phase 5** below, where they consume the SI-status-validated parent list produced here.)

---

## Phase 5: Testing S-allele Hypotheses and Cross Design

> Phase 5 uses the allele bin definitions from Phase 2, the individual GFS data from Phase 3, *and the per-individual SI status from Phase 4* to design and analyse controlled crossing experiments. All scripts run from the same working directory as the earlier phases. Cross design (Step 22e) is restricted to parents whose Phase 4 SI status is `SI` — broken-SI individuals cannot serve as reliable compatibility-prediction parents.

### Step 22 — HV-Based Allele Hypothesis Testing and Crossing Design

> Step 22 is now a **two-script workflow**. First, the variability landscape is built from a multi-species alignment (LEPA + Brassica + Arabidopsis SRK alleles) and HV columns are detected by Shannon entropy with a permutation-tested cross-Brassicaceae overlap. The validated LEPA HV columns are then consumed by the allele-hypothesis testing script (HV-distance matrix → UPGMA classes → cross design → synonymy network).

#### Step 22a — Cross-Brassicaceae S-domain Variability Landscape

**Scripts (in execution order):**
1. `srk_fetch_reference_alleles.py` — one-time NCBI fetch of *Brassica rapa*, *B. oleracea*, *Arabidopsis lyrata*, and *A. halleri* SRK reference proteins (≈22 + 10 sequences after S-haplotype dedup). Output: `all_reference_SRKs_dedup.fasta`. **Re-run only if the reference set changes** (rare).
2. `pad_representatives.py` — pads `SRK_protein_allele_representatives.fasta` (Step 10a output) to uniform length. Output: `SRK_protein_allele_representatives_padded.fasta`. **Re-run every time Step 10a re-runs** (new library, re-calibrated `N_ALLELES`, or any filter change upstream).
3. `mafft --add` — builds the combined LEPA + Brassica + Arabidopsis alignment. **Re-run every time pad_representatives runs.**
4. `srk_brassica_hv_mapping.py` — maps the 12 SCR9-contact residues from Ma et al. 2016 (PDB 5GYY, *B. rapa* eSRK9) to LEPA alignment coordinates. **Re-run every time the combined alignment is rebuilt** (LEPA column coordinates shift when LEPA representatives change).
5. `srk_variability_landscape.py` — primary Step 22a script: per-species Shannon entropy + per-species HV-region calls + permutation tests + figure.

**Commands:**
```bash
# (Once per project, or when reference set changes)
python3 srk_fetch_reference_alleles.py

# Rebuild the combined alignment whenever Step 10a re-runs:
python3 pad_representatives.py
mafft --add all_reference_SRKs_dedup.fasta SRK_protein_allele_representatives_padded.fasta \
      > SRK_combined_alignment.fasta
python3 srk_brassica_hv_mapping.py

# Variability landscape
python3 srk_variability_landscape.py
```

**Why the alignment must be rebuilt after every Step 10a re-run:**
The LEPA allele set defines the *columns* of the combined alignment that downstream Step 22 scripts treat as canonical (HV column indices in `SRK_LEPA_HV_positions.tsv`, SCR9-contact residue mappings in `SRK_brassica_hv_mapping.tsv`). Re-running Step 10a with a different `N_ALLELES` produces a different representative set (different sequences, possibly different lengths after trailing-gap trimming), which means the `mafft --add` profile shifts and column coordinates renumber. Skipping the rebuild would silently mix Step 22 outputs that reference different LEPA coordinate systems — a hard-to-debug failure mode.

**Validation step:** after rebuild, compare the printed HV region spans against the previous run. Substantial shifts (>10 columns) indicate the LEPA representative set has changed enough to invalidate downstream Step 22 outputs computed against the old combined alignment.

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

**Headline (current post-QC dataset, 2026-05-11).** The cross-Brassicaceae HV-overlap permutation now resolves the three pairs differently. **Brassica ↔ Arabidopsis** remains highly significant (obs = 43 cols, *p* < 0.0001), confirming the conserved selection signature between the two well-sampled genera. The two **LEPA pairs are no longer significant** — LEPA ↔ Brassica obs = 10 cols (*p* = 0.18); LEPA ↔ Arabidopsis obs = 15 cols (*p* = 0.16). This is interpreted as **drift-eroded standing variation** at LEPA HV columns: 11 of the 12 mappable Ma 2016 SCR9-contact residues still fall within the LEPA HV regions (structural validation is intact), and the 66 canonical LEPA HV columns are unchanged — but allele richness inside each column has collapsed to the point where the permutation null can no longer reject random overlap. LEPA additionally retains a unique HV peak at LEPA cols 358–388 not present in Brassica/Arabidopsis (candidate Lepidium-specific specificity site). Wu-Kabat sanity check: Jaccard with Shannon-entropy HV cols = 0.45 (LEPA), 0.88 (Brassica), 0.75 (Arabidopsis) — strong concordance under both metrics. See Step 22d for the per-BL entropy decomposition that quantifies the same drift signature.

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
| `WITHIN_CLASS_THRESHOLD` | `0.04` | HV p-distance boundary between Synonymy_test (uncertain) and Compatible_within (compatible) |
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
| 3b | Cross design summary figure (distance distribution, category schematic, Synonymy_test cross interpretation) |
| 3c | Synonymy network: synonymy-group groups figure + N-connectivity condensed figure; per-allele synonymy group CSV |
| 4 | UPGMA dendrogram + AAAA availability bar chart |
| 5 | Cross result analysis (activated when `CROSS_TSV` is set; → Step 23) |

**Outputs:**

| File | Content |
|------|---------|
| `SRK_HV_allele_distances.tsv` | Pairwise HV-only distance matrix (alleles × alleles) computed on the canonical HV columns |
| `SRK_allele_similarity_heatmap.pdf` | HV similarity heatmap ordered by UPGMA; colour scale spans within-class range |
| `SRK_functional_allele_groups.tsv` | Allele bin → class assignment, AAAA count, cross power (full / singleton / none) |
| `SRK_synonymy_candidates.tsv` | All within-class allele pairs with HV distance, cross tier, and testability flag |
| `SRK_AAAA_cross_design_HV.tsv` | AAAA × AAAA pairs ranked Incompatible → Synonymy_test → Compatible_within → Compatible_cross, with HV distance and expected outcome |
| `SRK_cross_design_summary.pdf` | Three-panel figure: HV distance distribution, cross category schematic, Synonymy_test cross interpretation |
| `SRK_HV_cluster_figure.pdf` | UPGMA dendrogram coloured by class + AAAA availability bar chart |
| `SRK_synonymy_network_groups.pdf` | HV-identical synonymy groups + isolated singletons |
| `SRK_synonymy_network_tests.pdf` | Condensed super-node graph (synonymy groups + synonymy-test bridge edges; node size ∝ individuals) |
| `SRK_synonymy_groups.csv` | Per-allele Synonymy group membership and counts |
| `figures/*.png` | PNG copies of all figures at 200 dpi |

**Headline (current dataset, 2026-05-11, 66 canonical HV columns):** UPGMA splits 58 alleles into 2 functional groups (FG01 = 57 alleles, FG02 = Allele_055 outlier). Synonymy network: **8 synonymy groups + 19 isolated singletons → 27 effective bins** (down from 58), driven by Synonymy group 1 (15 HV-identical alleles incl. Allele_050 + Allele_051 = the pan-BL fixed S-specificity, 134 AAAA individuals). 9700 Incompatible pairs, 10343 Synonymy_test pairs, 2112 Compatible_within pairs in the AAAA × AAAA cross-design matrix (22 155 total).

**Cross categories:**

| Category | HV distance | Definition | Expected outcome |
|----------|-------------|------------|-----------------|
| Incompatible | d = 0 | HV-identical alleles | No seeds — incompatibility predicted by sequence identity |
| Synonymy_test | 0 < d < threshold | Small HV differences, same class | Unknown — synonymy test: does this substitution change specificity? |
| Compatible_within | d ≥ threshold, same class | Substantial within-class HV divergence | Seeds expected — within-class positive control |
| Compatible_cross | different class | Between phylogenetic classes | Seeds expected — guaranteed compatible (Class I × Class II) |

> **Interpreting Synonymy_test cross outcomes:** incompatible Synonymy_test cross → allele bins share SI specificity → merge bins (synonymous alleles). Compatible Synonymy_test cross → small HV difference is functionally real → bin boundary confirmed.

---

#### Step 22c — Synonymy group Collapse Diagnostic (optional sanity check)

> Optional diagnostic. Run AFTER Steps 22a and 22b. Tests whether LEPA's low per-column Shannon entropy (vs Brassica / Arabidopsis) is driven by Synonymy group redundancy (multiple HV-identical alleles dragging down per-column diversity) or by other mechanisms (drift-purged rare residues, shallow phylogeny, polyploidy-relaxed selection). Re-runnable safely after every 22a → 22b cycle when new data arrive.

**Script:** `srk_wgroup_collapse_test.py`

**Command:**
```bash
python3 srk_wgroup_collapse_test.py
```

**Inputs (all from earlier Step 22 outputs — no new dependencies):**

| File | From | Description |
|------|------|-------------|
| `SRK_synonymy_groups.csv` | Step 22b | Synonymy group membership per allele |
| `SRK_combined_alignment.fasta` | Step 22a | LEPA + Brassica + Arabidopsis aligned |
| `SRK_LEPA_HV_positions.tsv` | Step 22a | Canonical LEPA HV columns (validated; not modified by this script) |

**Method:** for each synonymy group, picks one representative allele (the one with the most AAAA individuals; ties broken by lowest allele ID for reproducibility); keeps every isolated allele as-is. This collapses LEPA from 63 sequences to **32** (= 9 Synonymy group representatives + 23 isolated). The same Shannon-entropy scan as Step 22a is then run on the collapsed LEPA set, side-by-side with the original LEPA / Brassica / Arabidopsis profiles.

**Outputs (do NOT overwrite Step 22a):**

| File | Content |
|------|---------|
| `SRK_LEPA_synonymy_group_representatives.tsv` | One row per LEPA allele: kept (representative / isolated) or dropped (redundant), with Synonymy group and chosen representative |
| `SRK_wgroup_collapse_entropy_summary.tsv` | Numeric before/after summary: n, median/mean/max H, threshold, HV runs, HV cols per species |
| `SRK_variability_landscape_wgroup_collapsed.pdf` / `figures/SRK_variability_landscape_wgroup_collapsed.png` | Side-by-side variability landscape: LEPA full + LEPA collapsed + Brassica + Arabidopsis |

**Headline (current dataset):** LEPA full median entropy = 0.118 bits (n = 63); LEPA collapsed = 0.201 bits (n = 32, **+0.083 bits = +70%**); Brassica = 0.700 bits (n = 22). Synonymy group redundancy contributes to LEPA's low entropy (collapsing nearly doubles per-column entropy) but the **gap with Brassica is only ~12% closed** by collapse alone — the dominant mechanism is biological (drift-purged rare residues + shallow within-LEPA phylogeny + possibly polyploidy-relaxed selection on individual SRK alleles), not just sequence redundancy. HV-column count is robust to collapse: 73 → 59 (only −19%), so the HV-region *locations* are not artefacts of Synonymy group redundancy.

**Pipeline-order requirement:** Step 22a → 22b → 22c. Re-running 22a or 22b alone after new data lands does NOT regenerate the diagnostic; re-run 22c afterwards if you want the updated comparison.

---

#### Step 22d — Per-BL Entropy Decomposition: Drift vs Selection (optional sanity check)

> Optional diagnostic. Run AFTER Steps 13 and 22a (does NOT depend on 22b or 22c). Tests whether LEPA's residual low Shannon entropy at HV cols reflects **independent drift per Bottleneck Lineage** (each BL fixing a different residue) or **selection convergence** (all BLs sharing the same residue under common functional constraint), and whether that within-LEPA convergence corresponds to pan-Brassicaceae conservation or LEPA-specific drift on a shared ancestral pool.

**Script:** `srk_perBL_entropy_test.py`

**Command:**
```bash
python3 srk_perBL_entropy_test.py
```

**Inputs (all from existing pipeline outputs — no new dependencies):**

| File | From | Description |
|------|------|-------------|
| `SRK_combined_alignment.fasta` | Step 22a | LEPA + Brassica + Arabidopsis aligned |
| `SRK_LEPA_HV_positions.tsv` | Step 22a | Canonical LEPA HV columns |
| `SRK_individual_BL_assignments.tsv` | Step 13 | BL membership per individual |
| `SRK_individual_allele_table.tsv` | Step 11 | AAAA-allele assignments |
| `SRK_individual_zygosity.tsv` | Step 12 | Genotype patterns (filter to AAAA) |

**Method:** restrict to AAAA individuals (164 species-wide; 158 BL-assigned). For each LEPA HV column:
1. Compute the dominant residue and its frequency in each BL's AAAA gene pool.
2. Compute per-BL Shannon entropy (residue diversity within each BL).
3. Compute within-LEPA concordance (how many of the 5 BLs share the same dominant residue).
4. Compute the dominant residue in **Brassica** and **Arabidopsis** at the same alignment column (cross-genera comparison).

**Two-axis interpretation matrix:**

| Within-LEPA (5/5 BLs concordant?) | Cross-genera (LEPA = Brassica = Arabidopsis dominant?) | Interpretation |
|---|---|---|
| Yes (all 5 BLs same) | Yes (LEPA = both genera) | Pan-Brassicaceae selection — residues conserved across 25+ My of evolution |
| Yes (all 5 BLs same) | No (LEPA-specific residue) | **Drift on a shared LEPA ancestral pool** — the most-common ancestral allele was already at high frequency pre-bottleneck; drift fixed it independently in every BL |
| No (BLs discordant) | n/a | Independent drift per lineage — each BL stochastically fixed a different residue |

**Outputs:**

| File | Content |
|------|---------|
| `SRK_perBL_HV_residue_table.tsv` | Per-(HV col × BL) row: n_AAAA, dominant residue, dominant frequency, Shannon H, full residue counts |
| `SRK_perBL_HV_concordance_summary.tsv` | Headline statistics: within-LEPA concordance distribution, cross-genera match counts, verdict |
| `SRK_perBL_HV_concordance_summary_per_hv_col.tsv` | Per-HV-col detail: LEPA consensus residue, n BLs concordant, Brassica dominant, Arabidopsis dominant, match flags |
| `SRK_perBL_entropy_figure.pdf` / `figures/SRK_perBL_entropy_figure.png` | Two-panel: per-BL entropy heatmap (top) + dominant-residue heatmap with Brassica + Arabidopsis comparison rows (bottom) |

**Headline (current dataset):**
- Within-LEPA: **all 5 BLs share the same dominant residue at 66/66 HV cols (100 % concordance)** in the current dataset (2026-05-11); per-BL mean entropy = 0.000–0.042 bits.
- Cross-genera match: LEPA dominant = Brassica dominant at only **16/73 cols (22%)**; LEPA = Arabidopsis at 15/73 (21%); LEPA = both at **12/73 (16%)**.
- LEPA-specific (LEPA all-BL concordant + ≠ both Brassica AND Arabidopsis): **54/73 cols (74%)**.
- **Verdict: drift on a shared LEPA ancestral pool, NOT pan-Brassicaceae selection convergence.** Every BL fixed the same dominant allele family (Synonymy group 1 / Allele_050+057+relatives) because it was the most common ancestral haplotype pre-bottleneck. The dominant residues at 74% of HV cols are LEPA-specific (different from Brassica AND Arabidopsis), ruling out pan-genera selection as the primary driver. Only 16% of HV cols are conserved across all three genera — these represent the deeply conserved SI-recognition surface under selection across Brassicaceae.

**Combined Step 22a + 22c + 22d biological narrative**: LEPA's low Shannon entropy at HV cols (Step 22a) is mostly explained by independent drift in every BL converging on the same ancestral Synonymy group 1 alleles (Step 22d), with synonymy-group sequence redundancy contributing a smaller share (Step 22c). The 66 LEPA HV cols are real — anchored structurally by 11 of 12 Ma 2016 SCR9-contact residues — but in the current post-QC dataset (2026-05-11) the cross-genera HV-overlap permutation is no longer significant for LEPA pairs (LEPA ↔ Brassica p = 0.18; LEPA ↔ Arabidopsis p = 0.16); Brassica ↔ Arabidopsis remains highly significant (p < 0.0001). The honest reading: drift has eroded LEPA's HV signal beyond statistical detectability against an unbiased baseline — the same conclusion drawn from the per-BL entropy decomposition.

**Pipeline-order requirement:** Step 22a (HV cols) + Step 13 (BL bridge) → 22d. Optional whether 22b/22c have run.

---

#### Step 22e — Hypothesis-Testing Cross Plan Generator

> Operational deliverable. Translates the four-way cross categorisation from Step 22b (Incompatible / Synonymy_test / Compatible_within / Compatible_cross) into a phased experimental protocol that tests S-allele specificity hypotheses with explicit genotype constraints. Run AFTER Steps 22a and 22b; uses Step 13 BL assignments and Step 11/12 genotype + composition outputs. Designed to be re-run when new individuals are added.

**Script:** `srk_cross_plan.py`

**Command:**
```bash
python3 srk_cross_plan.py
```

**Hypothesis-testing framework — rationale**

Each sequence-defined allele bin is treated as a *hypothesis*: "the proteins clustered into this bin share one SI specificity." Crosses test these hypotheses by checking whether the predicted compatibility category (Incompatible / Synonymy_test / Compatible_within / Compatible_cross) matches the observed seed yield. The plan structures this as a chain of nested hypotheses, each level depending on the previous one's outcome:

| Level | Question | Maternal genotype | Paternal genotype | Allele constraint | Why |
|-------|----------|-------------------|-------------------|-------------------|-----|
| **H0** | Does SI rejection actually work? | AAAA | AAAA | both alleles in the **same** synonymy group | AAAA × AAAA gives unambiguous SI specificity on both sides — pistil expresses one identity, all pollen carries one identity. Any seed set is direct evidence of SI breakdown. |
| **H1a** | What is the within-Class compatible-cross seed-yield baseline? | AAAA | AAAA | alleles in **different** synonymy groups, NO Synonymy_test edge (HV ≥ 0.04, same Class) | Both genotypes unambiguous; predicted compatible. Anchors the upper baseline of seed yield within a Class. |
| **H1b** | What is the between-Class (max) seed-yield baseline? | AAAA Class I | Heterozygous (AAAB / AABB / AABC) carrier of a Class II allele | Mother allele ≠ any of father's other Class I alleles | Class II allele is between-class compatible (P_cross-style). Heterozygous father is required because Class II has zero AAAA carriers in the dataset. Pollen segregation (e.g., AABB → 17% AA / 67% AB / 17% BB) means most pollen carries the Class II specificity. |
| **H2** | Do bins separated by small HV differences correspond to distinct specificities? | AAAA | AAAA | alleles in **different** synonymy groups, WITH Synonymy_test edge (0 < HV < 0.04, same Class) | Same genotype rigour as H0. Outcome interpretation: 0 seeds → MERGE the two synonymy groups (synonymous specificities); H1a-baseline yield → bin boundary functionally real. |
| **H3** | Do the 29 bins lacking AAAA representatives follow the same compatibility hierarchy? | AAAA carrier of an allele M known to be Compatible_within with the father's main allele | Heterozygous (AAAB / AABB) carrier of the **hidden** allele | M ≠ father's main allele AND M ≠ hidden allele | Mother AAAA so pistil specificity is unambiguous; father heterozygous because the hidden bin doesn't exist as AAAA. **Requires a paired control** — same mother × an AAAA father carrying the father's main allele only — to measure the AA-pollen baseline yield. The hidden allele's effect is the *additional* yield above that baseline. |

**Provenance — every hypothesis assignment is traceable**

Each cross in the plan is assigned to its hypothesis level by combining two independent axes of evidence:

*Axis 1 — Sequence-based cross category (Step 22a + 22b):*

| Evidence | Source script | Source file | Threshold rationale |
|---|---|---|---|
| HV columns | `srk_variability_landscape.py` | `SRK_LEPA_HV_positions.tsv` | Shannon entropy > mean + 1×SD on smoothed profile (window = 20 aa); min run = 3 cols. Anchored structurally by overlap with Ma 2016 SCR9-contact residues (11 of 12 contacts mapped). In the current cleaner post-QC dataset (2026-05-11), the cross-Brassicaceae HV-overlap permutation is no longer significant for LEPA pairs (LEPA↔Brassica p = 0.18, LEPA↔Arabidopsis p = 0.16); Brassica↔Arabidopsis remains highly significant (p < 0.0001) — interpretation: drift has eroded LEPA's HV signal beyond statistical detectability. |
| Pairwise HV distance | `srk_allele_hypotheses.py` Part 2 | `SRK_HV_allele_distances.tsv` | p-distance computed only on the 66 canonical HV columns. |
| Class assignment | `srk_allele_hypotheses.py` Part 2 (UPGMA) | `SRK_functional_allele_groups.tsv` | Auto-detected at the largest gap in UPGMA merge heights (0.087 → 0.96 in current data) — corresponds to the well-documented Brassicaceae Class I / Class II split. |
| Synonymy group membership | `srk_allele_hypotheses.py` Part 3c | `SRK_synonymy_groups.csv` | Connected components of the graph where allele pairs are linked by HV distance = 0. |
| Cross category | `srk_allele_hypotheses.py` Part 3 | `SRK_AAAA_cross_design_HV.tsv` | `different Class → Compatible_cross`; `d = 0 → Incompatible`; `0 < d < WITHIN_CLASS_THRESHOLD → Synonymy_test`; `d ≥ WITHIN_CLASS_THRESHOLD → Compatible_within`. Default threshold 0.04 ≈ median within-class HV distance. |

*Axis 2 — Genotype-based feasibility (this script):*

| Evidence | Source | File | Why it matters |
|---|---|---|---|
| Per-individual genotype | Step 12 | `SRK_individual_zygosity.tsv` (`Genotype` column) | Determines whether an individual is a clean specificity donor (AAAA) or a heterozygous carrier (AAAB/AABB/AABC/ABCD). |
| Per-individual allele composition | Step 12 | `SRK_individual_zygosity.tsv` (`Allele_composition` column) | Authoritative source for which allele is the major (3-copy) vs minor (1-copy) in heterozygous individuals. The cross-plan parser reads this string directly (not the per-protein `SRK_individual_allele_table.tsv`). |
| BL membership | Step 13 | `SRK_individual_BL_assignments.tsv` | Tiebreaker — between-BL pairings preferred for genetic-diversity benefit. |

*Axis combination* — every cross row is the result of one cell in this matrix:

| Cross category (Axis 1) | AAAA × AAAA (clean) | AAAA × heterozygous |
|---|---|---|
| Incompatible (HV = 0, same Class) | **H0** | excluded — heterozygous father has ambiguous specificity |
| Compatible_within (HV ≥ 0.04, same Class) | **H1a** | **H3** (when AAAA mother + heterozygous donor of hidden bin; requires paired AAAA × AAAA control) |
| Compatible_cross (different Class) | not possible — no AAAA in Class II | **H1b** (AAAA Class I × heterozygous carrier of Class II allele) |
| Synonymy_test (0 < HV < 0.04, same Class) | **H2** | deferred — AAAA × AAAA is cleaner |

**Assumptions made explicit (for reviewer-facing justification)**

1. **The 66 HV columns are the SI-recognition surface.** Justified by structural overlap with the Ma 2016 *B. rapa* eSRK9–SCR9 crystal structure (PDB 5GYY, 11 of 12 contact residues map within or adjacent to LEPA HV regions). In the current cleaner post-QC dataset (2026-05-11), the cross-Brassicaceae HV-overlap permutation test is no longer significant for LEPA pairs (LEPA ↔ Brassica p = 0.18; LEPA ↔ Arabidopsis p = 0.16) but Brassica ↔ Arabidopsis remains highly significant (p < 0.0001) — the LEPA HV signal is interpreted as drift-eroded standing variation on a shared LEPA ancestral pool (see Step 22d). The structural anchor from Ma 2016 carries the validation.
2. **HV-only p-distance approximates SI-specificity divergence.** This is the working hypothesis the H2 synonymy tests will *empirically validate or refute* — it is not assumed correct.
3. **The UPGMA largest-gap cut correctly identifies Class I vs Class II.** Validated by the strongly bimodal HV-distance distribution and consistency with the documented Brassicaceae class split (current data 2026-05-11: Allele_055 is the single Class II allele, separated from all Class I alleles by a much larger HV distance).
4. **`WITHIN_CLASS_THRESHOLD = 0.04` separates small from substantial HV differences.** Empirical, derived from the within-class HV-distance distribution. Sensitivity analysis is straightforward (re-run Step 22b with a different value).
5. **AAAA × AAAA crosses give unambiguous specificity assignment.** Direct consequence of polyploid SI biology: AAAA individuals carry one specificity in pistil and four identical copies in pollen.
6. **Polyploid pollen segregation follows the C(4,2) = 6 combination rule.** Standard tetraploid meiotic assumption. H1b and H3 cross interpretations explicitly account for the resulting AA / AB / BB pollen distributions per genotype (AAAB → 50/50 AA/AB; AABB → 17/67/17 AA/AB/BB; etc.).
7. **Heterozygous-donor seed yields are interpretable only when paired with AAAA-only baseline controls.** This is *why* H3 includes paired controls — without them, the AAAB-pollen yield is confounded by the AA-pollen baseline.

**Inputs** (all from upstream Step 22 / Step 13 / Step 11–12 outputs):

| File | From | Purpose |
|------|------|---------|
| `SRK_AAAA_cross_design_HV.tsv` | Step 22b | Pre-computed cross category for all AAAA × AAAA pairs |
| `SRK_synonymy_groups.csv` | Step 22b | Synonymy-group membership per allele |
| `SRK_HV_allele_distances.tsv` | Step 22b | Full 63 × 63 HV-distance matrix (used to evaluate AAAA × heterozygous pairings) |
| `SRK_functional_allele_groups.tsv` | Step 22b | Allele → Class I / Class II |
| `SRK_individual_BL_assignments.tsv` | Step 13 | BL per individual (tiebreaker for between-BL preference) |
| `SRK_individual_zygosity.tsv` | Step 12 | `Genotype` and `Allele_composition` per individual |
| `SRK_individual_allele_table.tsv` | Step 11 | Per-protein allele assignments (cross-checked against `Allele_composition`) |

**Key parameters (top of script):**

| Parameter | Default | Notes |
|-----------|---------|-------|
| `WITHIN_CLASS_THRESHOLD` | 0.04 | Must match Step 22b. |
| `REPLICATES["H0"]` | 9 | 3 mothers × 3 flowers — enough to detect any SI leakage. |
| `REPLICATES["H1a"]` | 9 | Establishes the within-Class baseline tightly. |
| `REPLICATES["H1b"]` | 9 | Establishes the between-Class baseline. |
| `REPLICATES["H2"]` | 15 | 5 mothers × 3 flowers — higher because H2 outcomes are bimodal (uncertain). |
| `REPLICATES["H3"]` | 15 | Same as H2; plus paired controls per H3 cross. |
| `N_CROSSES["H0_within_allele"]` | 5 | Same-allele AAAA × AAAA in the largest synonymy group. |
| `N_CROSSES["H0_within_group"]` | 8 | Different-alleles-same-group within Synonymy group 1. |
| `N_CROSSES["H0_other_groups"]` | 12 | Spread across the remaining 8 synonymy groups (≈2 per group). |
| `N_CROSSES["H1a_far"]` | 10 | Maximally HV-divergent within-Class pairs. |
| `N_CROSSES["H1b"]` | 5 | Capped (sample-limited by Class II carriers). |
| `N_CROSSES["H2_per_pair"]` | 1 | One representative cross per synonymy-group-pair Synonymy_test bridge. |
| `N_CROSSES["H3_per_bin"]` | 1 | One main cross + one paired control per hidden bin. |

**Outputs:**

| File | Content |
|------|---------|
| `SRK_cross_plan_H0_SI_validation.tsv` | H0 — per-cross row: Mother / Father IDs, genotypes, alleles, BL pair, predicted Incompatible, replicate count, decision rule |
| `SRK_cross_plan_H1a_within_class_baseline.tsv` | H1a — within-Class compatible-baseline crosses (max-HV-divergent AAAA × AAAA) |
| `SRK_cross_plan_H1b_between_class_baseline.tsv` | H1b — between-Class baseline (Compatible_cross via heterozygous Class II carrier) |
| `SRK_cross_plan_H2_synonymy_tests.tsv` | H2 — one cross per inter-synonymy-group bridge, with merge / keep-separate decision rule |
| `SRK_cross_plan_H3_hidden_bin_tests.tsv` | H3 — main cross (AAAA mother × heterozygous hidden-bin carrier) + paired control (AAAA × AAAA, same mother + father main allele) |
| `SRK_cross_plan_summary.tsv` | Phase counts: N crosses, replicates per cross, total attempts, cumulative |
| `SRK_cross_plan_summary.pdf` / `figures/SRK_cross_plan_summary.png` | Two-panel figure: per-phase cross counts (left) + decision tree per phase outcome (right) |

**Headline (current dataset):**
- **101 unique crosses, 1 323 cross attempts** total.
- H0: 21 crosses (189 attempts). H1a: 10 (90). H1b: **1** (sample-limited — only 1 Class II carrier in dataset). H2: 32 (480 — one per synonymy-group-pair bridge). H3: 37 (555 — 20 hidden-bin tests + 17 paired AAAA × AAAA controls).
- 20 of 29 hidden bins are testable; the remaining 9 lack AAAA mothers carrying a Compatible_within partner allele and require additional genotyped plants.

**Pipeline-order requirement:** Step 13 + Step 22a + Step 22b → 22e. Re-run 22e after re-running upstream steps when new data arrive (the carrier-detection logic reads the authoritative `Allele_composition` strings, so the plan automatically expands as new individuals are added).

**Operational seed-orchard design** is downstream of 22e and 23 — it consumes the *validated* functional S-allele table that emerges from H2 + H3 outcomes, and is documented separately.

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
| Kruskal-Wallis H | Overall test of seed yield differences across Incompatible / Synonymy_test / Compatible_within / Compatible_cross |
| Mann-Whitney U (pairwise) | All category pairs |

**Outputs:**

| File | Content |
|------|---------|
| `SRK_cross_result_analysis_HV.pdf` | Seed yield box/strip plot + success rate bar chart by cross category |

---

### Step 25 — Per-Individual SI System Status (SI / pSI / SC)

> **Question answered:** *What is the status of the SI system?* — at the individual, bottleneck-lineage (BL), and population (EO) levels.

Step 25 reconstructs, for every ingroup individual, how many of the four expected tetraploid SRK copies produce a functional protein. It uses two inputs: the per-haplotype `OK` / `REMOVED` calls written by Step 7 (`translate_filter_align_AA.py`) to the per-library `*_frame1_stopcodon_log.tsv` files, and the per-haplotype broken-allele identity assignments produced by Step 25a (`SRK_null_allele_assignment.py`). Step 25a aligns every REMOVED haplotype's AA sequence to the 49 canonical functional allele representatives and flags chimeric Canu assembly artefacts (haplotypes that fail to identify with any functional allele by AA distance). Step 25b then counts only real broken haplotypes against the four-copy tetraploid expectation. The canonical 49-allele catalogue and all Phase 1–2 outputs are unchanged.

**Step 25a — broken-allele identity assignment.** For every haplotype tagged `REMOVED` by Step 7, the AA sequence is added to the canonical reference alignment via `mafft --add --keeplength` and p-distance to each of the 49 functional reps is computed. Confidence flag:

- **high**: `AA_distance ≤ 0.005` AND `n_stops ≤ 3` (one premature stop on an otherwise clean copy of the assigned allele)
- **medium**: `AA_distance ≤ 0.05` AND `n_stops ≤ 10`
- **low**: otherwise — chimeric Canu artefact, no biological signal

Output: `Tables/SRK_null_allele_assignments.tsv` (Sequence_ID, assigned_allele, AA_distance, first_stop_position, stop_codon_DNA, confidence).

**Step 25b — per-individual SI categorisation.** Real broken haplotypes are those flagged `high` or `medium` in 25a; chimeric haplotypes (`low`) are recorded in `n_haps_chimeric` for transparency but excluded from the tetraploid-copy signal pool:

    n_total              = n_haps_OK + n_haps_REMOVED_real
    copies_nonfunctional = round(4 × n_haps_REMOVED_real / n_total)
    SI_status:
      SI                copies_NF == 0   OR   n_haps_REMOVED_real < 2
      pSI               1 ≤ copies_NF ≤ 3   AND   n_haps_REMOVED_real ≥ 2
                        (pSI_confidence: high if frac_NF ≥ 0.25, low otherwise)
      SC                copies_NF == 4
      Insufficient_data n_total < 4

The `SRK_BEA` Brassica reference haplotypes are excluded before aggregation. The single-`REMOVED` floor on the SI side (`n_haps_REMOVED_real < 2 → SI`) is a residual noise margin after chimera filtering.

**Scripts:**

| Script | Role |
|---|---|
| `SRK_null_allele_assignment.py` | Step 25a — broken-allele identity assignment by AA p-distance to the 49 functional reps. |
| `SRK_individual_SI_status.py` | Step 25b — per-individual SI categorisation using chimera-filtered counts. |
| `SRK_SI_status_figures.R` | Renders the three stacked-bar figures (species, BL, EO) in `_full` and `_robust` variants. |

**Commands:**
```bash
/Users/sven/anaconda3/bin/python SRK_null_allele_assignment.py    # Step 25a
python3 SRK_individual_SI_status.py                                # Step 25b
Rscript SRK_SI_status_figures.R
```

**Inputs:**

| File | Description |
|---|---|
| `all_Library*_*_frame1_stopcodon_log.tsv` (×9) | Per-haplotype OK/REMOVED calls from Step 7. |
| `all_Library*_*_frame1_AA_raw.fasta` (×9) | AA sequences for every haplotype (Step 25a only). |
| `all_Library*_*_aligned_exons_backfilled.fasta` (×9) | Aligned DNA — used for the DNA codon at the first stop position (Step 25a). |
| `SRK_protein_allele_representatives.fasta` | 49 canonical functional allele reps (Step 25a). |
| `Tables/SRK_individual_allele_genotypes_with_nulls.tsv` | Used by Step 25a only to identify the canonical 49-allele set; the matrix itself is Step 26's output. |
| `Tables/SRK_data_quality_categories.tsv` | Ingroup flag + `EO_normalised` + `BL_inferred` (Step 12c). |

**Outputs:**

| File | Content |
|---|---|
| `Tables/SRK_null_allele_assignments.tsv` | One row per REMOVED haplotype: assigned functional allele, AA distance, first stop position + DNA codon, confidence. |
| `Tables/SRK_individual_SI_status.tsv` | One row per ingroup individual: `n_haps_OK`, `n_haps_REMOVED`, `n_haps_chimeric`, `n_haps_total`, `frac_nonfunctional`, `copies_functional`, `copies_nonfunctional`, `SI_status`, `pSI_confidence`. |
| `figures/SRK_SI_status_species_full.png` / `_robust.png` | Species-level stacked bar. `_full` (all 401 ingroup, 7-tier incl. pSI_low + Insufficient_data); `_robust` (drops Insufficient_data + pSI_low). |
| `figures/SRK_SI_status_by_BL_full.png` / `_robust.png` | Per-BL stacked %. Same full / robust split. |
| `figures/SRK_SI_status_by_EO_full.png` / `_robust.png` | Per-EO stacked %, faceted by BL with within-BL drift-index order. Same full / robust split. |

**Current species-level snapshot (n = 401 ingroup):** SI = 247 (62 %), pSI = 15 (3.7 %; 10 high-confidence + 5 low), SC = 1 (in EO76), Insufficient_data = 138 (heavily concentrated in EO76, flagged for re-sequencing in `Tables/SRK_samples_for_redo.tsv`). The molecular SI machinery is broadly intact; pSI is sparse and geographically diffuse; the only confirmed SC individual sits in EO76 (the SI → SC transitional hotspot identified in Q2).

---

### Step 26 — Null-Aware Tetraploid Genotype Rebuild (Phase 4)

> **Question answered:** *Once we know which individuals are pSI / SC, what does the population genetic picture look like when we stop treating their broken copies as functional?* — the integration step that propagates Step 25's per-individual SI status through to every downstream population-genetic metric (Steps 14, 17, 19, 20).

**Why this step exists.** The canonical Phase-1–2 genotype tables count only *functional* SRK copies and pad under-recovered individuals to four slots by homozygosity assumption (Step 12). For SI individuals this is correct. For pSI individuals the homozygosity padding fabricates functional copies that biologically do not exist; for SC individuals the entire genotype is missing because they fail at Step 9. Step 26 propagates Step 25's per-individual SI status into a parallel null-aware genotype matrix so downstream Phase-3 metrics can be recomputed with explicit broken alleles. The canonical Phase-1–2 tables stay frozen as the functional-only reference.

**Categorisation rule:**

| Source SI_status | Action in null-aware tables | Rationale |
|---|---|---|
| `SI` (n=247) | Pad to 4 functional copies; `Allele_NULL = 0`. | Truly SI — no nulls. |
| `pSI` + `pSI_confidence == low` (n=5) | Promote to SI; pad to 4 functional copies; `Allele_NULL = 0`. | `frac_NF < 0.25` is too dilute to assign nulls confidently. |
| `pSI` + `pSI_confidence == high` (n=10) | Scale existing allele counts down to `copies_functional`; set `Allele_NULL = copies_nonfunctional`. | Genuine partial SI — broken copies are real and enter the allele-frequency denominators. |
| `SC` (n=1) | All four slots → `Allele_NULL`. | Self-compatible escape (Library010_barcode53, EO76 / BL3). |
| `Insufficient_data` (n=138) | Kept in canonical set if they reached genotyping; flagged for re-sequencing in `SRK_samples_for_redo.tsv` either way. | Genotype call is data-thin; redo to resolve. |

The scaling uses **largest-remainder proportional rounding** so every row in the new genotype matrix still sums to four. Original Phase-1–2 tables stay frozen as the functional-only reference.

**Scripts:**

| Script | Role |
|---|---|
| `SRK_genotype_null_alleles.py` | Builds `SRK_individual_allele_genotypes_with_nulls.tsv`, `SRK_individual_zygosity_with_nulls.tsv`, and `SRK_samples_for_redo.tsv`. |
| `SRK_population_genetic_summary_with_nulls.R` | Step 14b — null-aware allele frequencies, He, Ho, Ne, plus new columns `Frac_nonfunctional_alleles`, `Mean_null_copies`, `Prop_pSI`, `N_SC`. |
| `SRK_TP1_compatibility_metrics_with_nulls.py` | Step 17b — null-aware P_compat. Treats Allele_NULL as a regular allele in the frequency vector (heuristic; conservative for high-p_NULL populations). |
| `SRK_individual_GFS_with_nulls.R` | Steps 19b + 20b — null-aware GFS using `GFS = GFS_func × (n_func / 4)`. |
| `SRK_P_compat_traffic_light_with_nulls.R` | Stakeholder-facing null-aware random-mating traffic-light per EO (red / amber / green at the 0.20 and 0.40 P_compat thresholds). |
| `SRK_depletion_ranking_with_nulls.R` | Null-aware DI × P_compat conservation ranking (observed + MM-predicted variants, EO + BL points). |

**Commands:**
```bash
python3 SRK_genotype_null_alleles.py
Rscript SRK_population_genetic_summary_with_nulls.R
/Users/sven/anaconda3/bin/python SRK_TP1_compatibility_metrics_with_nulls.py
Rscript SRK_individual_GFS_with_nulls.R
Rscript SRK_P_compat_traffic_light_with_nulls.R
Rscript SRK_depletion_ranking_with_nulls.R
```

**Outputs (key tables and figures):**

| File | Content |
|---|---|
| `Tables/SRK_individual_allele_genotypes_with_nulls.tsv` | 336 individuals × 49 functional alleles + Allele_NULL + Genotype_class_flag + EO + BL + SI_status. Rows sum to 4. |
| `Tables/SRK_individual_zygosity_with_nulls.tsv` | Per-individual genotype label (e.g. `AABC`, `AB00`, `0000`) + N_distinct_functional + N_null_copies + SI_status. |
| `Tables/SRK_samples_for_redo.tsv` | 138 Insufficient_data samples + EO + BL + hap counts + redo priority. |
| `SRK_population_genetic_summary_with_nulls.tsv` / `_BL_with_nulls.tsv` / `.pdf` | EO- and BL-level null-aware pop-gen tables and 2-page PDF figure. |
| `Tables/SRK_EO_allele_richness_with_nulls.tsv` | Null-aware TP1 metrics (P_compat at L=0/0.10/0.25/0.50 with bootstrap CIs). |
| `SRK_individual_GFS_with_nulls.tsv` / `SRK_EO_GFS_summary_with_nulls.tsv` / `SRK_BL_GFS_summary_with_nulls.tsv` | Null-aware GFS + TP2 status. |
| `figures/SRK_GFS_with_nulls_composition.png` | Stacked bars: genotype-tier composition per BL with null-aware tiers. |
| `figures/SRK_GFS_with_nulls_TP2_scatter.png` | Null-aware TP2 scatter: mean GFS × prop_zero with BL/EO markers. |

**Null-aware GFS formula:**

    GFS_func        = 1 - sum_k[ n_k(n_k - 1) ] / (n_func × (n_func - 1))
    GFS_null_aware  = GFS_func × (n_func / 4)

where `n_func = 4 - N_null_copies`. The 12-tier ordering is:

| GFS | Genotypes |
|---:|---|
| 0.000 | 0000, A000, AA00, AAA0, AAAA |
| 0.500 | AAAB, AB00, AAB0 |
| 0.667 | AABB |
| 0.750 | ABC0 |
| 0.833 | AABC |
| 1.000 | ABCD |

SC individuals and severely null-loaded genotypes collapse to the zero tier. Compare to the canonical 5-tier scheme (AAAA = 0 / AAAB = 0.5 / AABB = 0.667 / AABC = 0.833 / ABCD = 1.0).

**Canonical vs null-aware deltas:**

| Metric | Canonical (functional-only) | Null-aware | Δ |
|---|---:|---:|---:|
| BL1 P_compat (L=0) | 0.068 | 0.091 | +0.023 |
| BL2 P_compat (L=0) | 0.106 | 0.112 | +0.006 |
| BL3 P_compat (L=0) | 0.394 | 0.406 | +0.012 |
| BL4 P_compat (L=0) | 0.417 | 0.422 | +0.005 |
| BL5 P_compat (L=0) | 0.472 | 0.481 | +0.009 |
| Mean GFS per BL | 0.17 – 0.37 | 0.17 – 0.37 | ≤ 0.01 |
| Sample size (canonical set) | 335 | 336 (+1 SC) | SC included with all-null genotype |

Per-BL broken-allele frequencies after chimera filtering are small (`frac_null` = 0.005 – 0.05), so the null-aware metrics **corroborate** the canonical Phase-3 conclusions rather than revise them. The null-aware companions are kept as a defensibility check.

---

### Step 27 — Forward-Time Inheritance Simulator (Phase 4)

> **Question answered:** *Where is each BL heading on the SI → SC erosion axis under current conditions, and what conservation lever stops the trajectory?* — the dynamic companion to Step 25's snapshot.

**Why this step exists.** Step 25 tells us the per-individual SI / pSI / SC distribution today. Step 26 tells us what that distribution implies for current population-genetic metrics. Neither tells us where the species is *heading*. A forward-time tetraploid Wright–Fisher simulator initialised from the empirical Step 26 state projects each BL forward in generation time under (drift × mutation × SI rejection × selfing × inbreeding depression × optional migration), so the conservation lever can be quantified rather than guessed.

**Model.** Each individual is a length-4 vector of allele identities, with 0 = `Allele_NULL` and positive integers indexing functional alleles. Per generation:

1. Each copy → NULL with probability μ (default 1 × 10⁻⁴).
2. Pick a mother uniformly at random. If she has ≥ 3 NULL stigma copies AND `U(0,1) < s`, she selfs (default `s = 0.5`). Otherwise she attempts outcrossing for up to 20 attempts.
3. Each outcross attempt: pick a random father, sample one of his gametes by tetrasomic random pairing (with double-reduction probability α = 0.10). Reject if any functional pollen S-allele matches any functional stigma S-allele (sporophytic SI). NULL pollen alleles never trigger rejection; NULL stigma copies never reject.
4. If no compatible father is found in 20 attempts, fall back to selfing.
5. Form the new individual from the mother's gamete + father's gamete.
6. Selfed offspring have fitness 1 − δ (default δ = 0.5). All offspring resampled with weights proportional to fitness (Wright–Fisher with selection).
7. Optional inter-BL migration: with probability *m* per individual, replace it with one drawn uniformly from a pool of all other BLs' current populations.

**Scenarios (default ladder):**

| Scenario | μ | s | δ | m | Ne_scaling |
|---|---|---|---|---|---|
| baseline | 1e-4 | 0.5 | 0.5 | 0 | 1.0 |
| rescue_low | 1e-4 | 0.5 | 0.5 | 0.001 | 1.0 |
| rescue_high | 1e-4 | 0.5 | 0.5 | 0.01 | 1.0 |
| high_drift | 1e-4 | 0.5 | 0.5 | 0 | 0.5 |

**Scripts:**

| Script | Role |
|---|---|
| `SRK_inheritance_simulator.py` | Runs all scenarios × all BLs × `--n_replicates` replicates × `--n_generations` generations. Outputs trajectories TSV + first-passage-time TSV. CLI flags: `--n_generations`, `--n_replicates`, `--mu`, `--scenarios`, `--seed`, `--quick`. |
| `SRK_inheritance_figures.R` | Renders three figures from the trajectories TSV. |

**Commands:**
```bash
/Users/sven/anaconda3/bin/python SRK_inheritance_simulator.py     # default: 200 gens × 50 reps
# Faster scaled-down run:
/Users/sven/anaconda3/bin/python SRK_inheritance_simulator.py --n_generations 100 --n_replicates 30
Rscript SRK_inheritance_figures.R
```

**Inputs:**

| File | Description |
|---|---|
| `Tables/SRK_individual_allele_genotypes_with_nulls.tsv` | Step 26 augmented genotype matrix — initial state per BL. |

**Outputs:**

| File | Content |
|---|---|
| `Tables/SRK_inheritance_trajectories.tsv` | One row per (scenario × BL × replicate × generation): n_SI, n_pSI1, n_pSI2, n_pSI3, n_SC, mean_p_NULL. |
| `Tables/SRK_inheritance_time_to_sc.tsv` | First-passage time to 50 % SC frequency per (scenario × BL × replicate). |
| `figures/SRK_inheritance_pNULL_trajectories.png` | Per-BL p_NULL trajectories — thin lines per replicate, bold line per scenario median. |
| `figures/SRK_inheritance_SC_progression.png` | SC-frequency trajectories with IQR ribbons + 50 % threshold line. |
| `figures/SRK_inheritance_time_to_sc.png` | Median + IQR time-to-50%-SC bars per BL × scenario. |

**Current results snapshot (100 gens × 30 reps).** Median generations to 50 % SC frequency:

| Scenario | BL5 | BL4 | BL3 | BL2 | BL1 |
|---|---:|---:|---:|---:|---:|
| baseline | 62 | 66 | 41 | 26 | > 100 |
| rescue_high (m=0.01) | > 100 | 72 | 54 | 41 | 13 |
| high_drift (N × 0.5) | > 100 | 94 | 41 | 65 | > 100 |

**Three biological readings:** (1) Erosion is real but slow — the fastest BL (BL2) crosses 50 % SC in ~26 generations under baseline, the larger BLs take 40 – 65, and BL1 (n = 7) is dominated by stochasticity. (2) Migration does not behave as a uniform rescue: it accelerates or decelerates the trajectory depending on whether the donor pool's broken-allele frequency is above or below the recipient's. (3) **The simulator-endorsed rescue lever is targeted SI-mother / SI-father crosses** (Step 22e cross plan, filtered through the Step 25 SI status table), not random inter-BL migration.

---

## Sample Exclusion Audit (cross-phase QC report)

**Script:** `audit_sample_exclusions.py`

**Purpose:** Trace every barcode in `sampling_metadata.csv` through the SRK pipeline and report, for each sample that did **not** make it to the final genotyped dataset, the **earliest** pipeline stage at which it was lost plus the specific reason. Distinguishes samples that need to be re-sequenced (poor DNA quality, fragmented assembly) from samples that have an intrinsic data-integrity problem (SRK paralogs / off-targets / contamination) where re-sequencing will not help, and from intentional exclusions (outgroup samples).

**Command:**
```bash
python3 audit_sample_exclusions.py
```

**Inputs (auto-detected in CWD):**
- `sampling_metadata.csv` — ground truth of expected samples
- `all_Library*_Phased_haplotypes.fasta` — Step 3 outputs (per-library raw haplotypes)
- `all_Library*_*_min3250*.fasta` — Step 4 length-filtered outputs
- `all_Library*_*_blastfilt_log.tsv` — Step 4b BLAST decisions per haplotype (Library 010+)
- `all_Library*_*_backfilled_Nfilt_log.tsv` — Step 7b N-filter decisions per haplotype (Library 010+)
- `SRK_individual_status_report.tsv` — Step 9 classifications
- `SRK_individual_BL_assignments.tsv` — Step 13 BL status
- `SRK_individual_zygosity.tsv` — final included individuals

**Output:** `SRK_sample_exclusion_audit.tsv` — one row per sample with columns:

| Column | Notes |
|--------|-------|
| `Sample_ID`, `Library`, `Barcode`, `EO`, `Ingroup` | from metadata |
| `n_raw_haps` | from Step 3 Phased haplotypes |
| `n_after_step4` | from Step 4 length-filtered FASTA |
| `n_after_step4b` | from Step 4b log (blank when library has no log) |
| `n_after_step7b` | from Step 7b log (blank when library has no log) |
| `Step9_classification` | from `SRK_individual_status_report.tsv` |
| `Proteins_in_final_data` | from same |
| `SI_functional_status` | Functional / Partial_translation_failure / Complete_loss / NA — see biological columns below |
| `Dominant_failure_mode` | premature_stop / ambiguous_aa / mixed / other — parsed from Step 9 Top_failure_reasons |
| `BL_status` | Assigned / Inferred / Unassigned |
| `Final_status` | Included / Excluded |
| `Exclusion_stage` | earliest stage at which the sample was lost (see cascade below) |
| `Exclusion_reason` | human-readable explanation + re-sequencing recommendation |

**Biological columns — SI functional status (added 2026-05-12):**

Two columns flag the *molecular* SI functional status of each sample independent of the QC cascade. These surface samples where sequencing was successful but no functional SRK proteins could be recovered — a candidate signal of ongoing self-compatibility evolution.

`SI_functional_status` categorisation (samples that reached Step 9):

| Value | Criterion | Interpretation |
|---|---|---|
| `Functional` | `Functional_rate ≥ 0.5` | SI system intact at the molecular level |
| `Partial_translation_failure` | `0 < Functional_rate < 0.5` (low_functional_rate) | included in dataset, but the label is **NEUTRAL** — does NOT imply biological SI loss. Most likely a technical artefact: Canu *de novo* assembly produces chimeric haplotypes whose spliced junctions introduce premature stops, which inflate the failure denominator without reflecting real biology. Step 9's abundance filter excludes chimeric proteins from the numerator (`Functional_proteins`) but not from `Total_sequences`. The current dataset provides no tangible evidence supporting biological SI loss in this category. **Do not make biological inferences from this label alone.** |
| `Complete_loss` | `Total_sequences > 0` AND `Functional_proteins == 0` (no_functional_proteins classification) | sequencing succeeded but no functional protein recovered — candidate SI-escape signal |
| `NA` | did not reach Step 9 (excluded earlier in pipeline) | not applicable |

`Dominant_failure_mode` parses Step 9's `Top_failure_reasons` string (e.g. `"premature_stop(28);ambiguous_aa(12)"`) and returns:

| Value | Meaning | Biological signal |
|---|---|---|
| `premature_stop` | internal stop codons dominate (≥ 2× the second most common failure) | **real loss-of-function mutations in SRK — strong SI breakdown candidate** |
| `ambiguous_aa` | X residues from N-rich data dominate | **data quality artefact, NOT a biological SI escape** (typical of libraries that bypassed Step 7b) |
| `mixed` | second most common failure ≥ 50 % of the top one | interpret cautiously — both LoF and quality issues contribute |
| `other` | different failure pattern (e.g., `too_short_protein`) | rare; investigate manually |
| (blank) | no failure data | sample did not reach Step 9 |

**Stdout summary** prints the `SI_functional_status` breakdown plus an explicit list of candidate SI-escape samples (those with `Complete_loss` + `premature_stop` dominated failures). These are the strongest molecular candidates for ongoing self-compatibility evolution and warrant follow-up phenotyping (controlled selfing tests).

**Cascade (earliest exclusion wins):**
1. **Stage 3 — no assembly:** sample in metadata but no haplotypes assembled by Canu → re-sequencing recommended.
2. **Step 4 — length filter:** raw haps exist, 0 passed 3250–4000 bp → re-sequencing recommended (fragmented assembly).
3. **Step 4b — BLAST coverage filter** (Library 010+): all haps dropped at coverage 0.90 / identity 95 → **NOT a re-sequencing issue** — sample contains SRK paralogs / off-targets; sequencing more reads will not help.
4. **Step 7b — N-content filter** (Library 010+): all haps dropped at max_N = 9 → re-sequencing *may* help if more reads yield longer 3′ contigs.
5. **Step 9 — translation / abundance / ploidy filters:**
    - `no_functional_proteins` — all sequences failed translation/validation → may benefit from re-sequencing.
    - `dropped_abundance_filter` — functional proteins present but below `min_count = 5` → more sequencing may push singletons over threshold.
    - `too_many_alleles` — >4 distinct proteins per individual → **NOT a re-sequencing issue** — likely contamination / mixed sample; clean re-extraction needed.
6. **Step 11 — outgroup filter:** `Ingroup = 0` in metadata → **intentional exclusion** (NOT a quality issue).

**Auxiliary flag:** `BL_status = Unassigned` (10 individuals with germplasm sub-codes) — these are STILL INCLUDED in the zygosity matrix; they just cannot be placed into a BL for Phase 3 BL-stratified analyses.

**When to run:** every time a new library is added, after Step 11 (genotyping) completes. The audit takes < 1 s. Use the output to (i) decide which barcodes to flag for re-sampling at the next field season, (ii) identify systematic failures suggesting a library-prep or assembly issue.

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
| Step 17 | `Tables/SRK_EO_allele_richness.tsv` |
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
