# Re-sequencing calibration — interpretation (post lab-machine diagnostics)

Built from the lab-machine run of `reseq_calibration.py` over **500 on-disk samples** (Libraries 001–011), joined to the 135-sample `Re-sequence (SI status uncertain)` cohort from `tables/Phase2/step12c_samples_redo.csv`, with new per-stage pipeline diagnostics + raw-read quality columns.

Source data: `tables/Phase4/step22c_reseq_calibration_per_sample.tsv`

---

## Headline finding — the bottleneck moved

Earlier interpretation (dev workstation, 157 samples) suggested most failures were *Canu refusing to assemble*. The full 500-sample diagnostics on the lab machine show that was wrong. The actual breakdown of the 135 Re-sequence cohort:

| failure_mode | n | What the diagnostics show |
|---|---:|---|
| **`pipeline_downstream_filter`** | **122 (90 %)** | Canu produced ~8 contigs / barcode. seqkit rmdup retained ~7. Step 2 chimera filter retained ~6. **The contigs survive everything up to and including the chimera filter — but then disappear before reaching the SI status call.** Read quality is fine (N50 ≈ 925, mean Phred ≈ 29.8). |
| `partial_recovery` | 7 | 1–3 clean haps, just under the 4-hap SI threshold. |
| `chimera_dominated` | 6 | Lots of chimeric haps despite adequate clean contigs surviving Step 2. |

So **the chimera defence (Step 2 coverage filter + Step 22a AA-distance filter) is not the bottleneck for any of the 122 main-failure samples** — confirming your structural reasoning. The bottleneck is somewhere in the **protein-level processing between the chimera filter and the SI status call**: Steps 4 (length filter), Step 4b (BLAST coverage / paralog filter), Step 7 (gap backfill), Step 7b (N-content filter), or Step 9 (abundance filter).

---

## The smoking gun — Libraries 009 and 011

100 % of the 122 `pipeline_downstream_filter` samples come from just two libraries:

| Library | pipeline_downstream_filter samples | % of failures |
|---|---:|---:|
| **Library 009** | **83** | 68 % |
| **Library 011** | **39** | 32 % |
| Library 010 | 0 | 0 % |
| All older libraries (001–008) | 0 | 0 % |

Library 010 was the first to be processed under the new Step 4b (BLAST coverage / paralog) + Step 7b (N-content) filters, and it shows zero downstream failures. Library 009 (older — pre-filter design) and Library 011 (newer, added after Library 010) are the failure cohort. That's a strong hint that the new filters interact poorly with the read profiles of those two libraries.

---

## Per-EO distribution

| EO | Total flagged | pipeline_downstream_filter | partial_recovery | chimera_dominated |
|---|---:|---:|---:|---:|
| **EO76** | 47 | 43 | 1 | 3 |
| EO70 | 30 | 26 | 3 | 1 |
| EO27 | 20 | 18 | 2 | 0 |
| **EO25** | 19 | **19** | 0 | 0 |
| EO08 | 6 | 6 | 0 | 0 |
| EO18 | 6 | 5 | 0 | 1 |
| EO67 | 5 | 5 | 0 | 0 |
| EO29 | 1 | 0 | 0 | 1 |
| EO97 | 1 | 0 | 1 | 0 |

The per-EO pattern is library-driven, not population-driven. **EO25 in particular is 100 % downstream-filter failures**, which would be impossible if the problem were population-level template difficulty — it's almost certainly that EO25's samples sit predominantly in Libraries 009 / 011.

### Critical implication for the LEPA report

The "EO76 = SI → SC transitional hotspot" framing in the LEPA report Q2 / Q6 / Conservation Recommendations was anchored on EO76 carrying the only confirmed SC individual + the heaviest Insufficient_data load. The latter looks substantially inflated by the Library 009/011 downstream-filter issue rather than real biology: 43 of EO76's 47 Insufficient_data calls are downstream-filter failures, not biological signal.

**Until the downstream-filter root cause is identified and fixed, the EO76 transitional-hotspot finding should be caveated as preliminary.** Once those 43 individuals get an SI call, EO76's profile may look substantially less broken.

---

## Quality data confirms the chimera / Canu reading

Per-sample quality summary, by failure mode:

| failure_mode | n | median raw_reads | median N50 | median Phred / read | median long reads (≥ 3 kb) | median Canu contigs | median unique contigs | median chimera KEEP |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `none — SI call achieved` | 263 | 20 185 | 980 | 29.8 | 228 | 9 | 8 | 6 |
| `pipeline_downstream_filter` | 122 | 16 543 | 925 | 29.8 | 222 | 8 | 7 | 6 |
| `PCR_amplicon_failure` | 88 | 6 149 | 839 | 27.6 | 69 | 8 | 4 | 3 |
| `partial_recovery` | 7 | 24 360 | 1 017 | 30.2 | 184 | 4 | 1 | 1 |
| `chimera_dominated` | 6 | 17 766 | 820 | 29.1 | 180 | 9 | 7 | 5 |
| `DNA_contamination` | 6 | 24 049 | 901 | 29.9 | 334 | 10 | 9 | 8 |

The contrast between **pass** and **`pipeline_downstream_filter`** is the headline: the two rows are statistically indistinguishable on raw reads, N50, Phred, long reads, Canu contigs, unique contigs, and chimera KEEP counts. Whatever filter is killing the Library 009/011 samples is operating downstream of the chimera filter on contigs that look just like the passing-sample contigs.

`PCR_amplicon_failure` (88 samples in the Re-PCR cohort) is genuinely under-sequenced (median 6 149 raw reads, 69 long reads, Phred 27.6) — those samples really do need a fresh PCR or higher coverage.

---

## Revised summary for the team

Instead of "we need to re-sequence 135 SI-uncertain samples," the data now supports:

> **13 samples** need physical lab work for an SI call:
> - 7 partial-recovery samples — re-sequence on existing DNA + library at higher coverage; high probability of success.
> - 6 chimera-dominated samples — need different amplicon prep (primer redesign, longer-fragment library prep); more reads alone won't fix.
>
> **122 samples need a bioinformatics fix** — their reads, Canu assemblies, and chimera-filtered contigs are indistinguishable from the 263 passing samples. They're being killed by a downstream filter (Step 4 / 4b / 7 / 7b / 9) that interacts badly with Libraries 009 and 011. **No lab re-work is needed for these — the bioinformatics pipeline needs investigation.** Once the root cause is found, most of these 122 should produce SI calls on existing data.

That's roughly a **10×** reduction in the lab workload (135 → 13) and a clear bioinformatics task to identify the Step 4–9 culprit.

---

## Suggested next bioinformatics investigation

The 122 Library 009/011 samples are indistinguishable from the passing samples at the chimera-filter output. So the right per-sample diagnostic is to **count contigs / proteins surviving each subsequent step**:

| Stage | Expected output file |
|---|---|
| Step 3 — orientation | `Library###/barcode##/Phased_haplotypes.fasta` (or equivalent) |
| Step 4 — length filter (3250–4000 bp) | post-length output |
| Step 4b — BLAST coverage / paralog | per-sample paralog-filtered FASTA |
| Step 7 — gap backfilling | `*_backfilled.fasta` |
| Step 7b — N-content filter | filtered backfilled FASTA |
| Step 9 — translation + abundance | per-library protein outputs |

If we extend `reseq_calibration.py` to count records at each of those checkpoints for the 122 samples, the failure mode will refine into something specific like `failed_at_step_4b_blast_coverage` or `failed_at_step_7b_n_content` — and that's the bug to fix.

Want me to add those per-Step contig / protein counts to the script for the next round?

---

## Classification thresholds used

- `LOW_RAW_READS = 10 000` — below this, classify as coverage-limited if Canu also produced nothing.
- Chimera-dominated: `n_haps_chimeric ≥ 2 × max(clean haps, 1)` AND `n_haps_chimeric ≥ 4`.
- Partial-recovery: `1 ≤ clean haps < 4`.
- `pipeline_chimera_filter_too_aggressive`: ≥ 80 % of contigs dropped at Step 2 AND `n_unique_contigs > 0`.
- `canu_assembly_failure_low_quality`: mean per-read Phred < 10.
- `canu_assembly_failure_fragmented`: N50 < 800 bp.
- `canu_assembly_failure_no_long_reads`: < 50 reads ≥ 3 kb.

None of the `canu_assembly_failure_*` sub-categories fired in the lab-machine run — the read quality + length signals are uniformly healthy across all 500 samples (median N50 ≈ 925 bp, median Phred ≈ 30 per read).

---

## Files

- `tables/Phase4/step22c_reseq_calibration_per_sample.tsv` — 500 rows × 37 columns (sort by `failure_mode` or `recommended_action` for the lab).
- `tables/Phase4/step22c_reseq_calibration_summary.tsv` — bin-level dose-response.
- `figures/Phase4/step22c_reseq_calibration.pdf` / `.png` — 3-panel figure.
- `reseq_calibration.py` — the script (project root on this workstation; same path on the lab machine).
