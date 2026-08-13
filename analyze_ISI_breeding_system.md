# ISI Breeding-System Analysis — Documentation

Companion documentation for `analyze_ISI_breeding_system.py`. Explains the biological question, the statistical method, the design choices made, the input and output specifications, and how to run the analysis.

---

## 1. Purpose

Classify the breeding system of *Lepidium papilliferum* (slickspot peppergrass) from experimental crosses of selfed (S) versus outcrossed (O) plants. Fruit-set data (percent fruit set per treatment × slick spot × site) are converted into an Index of Self-Incompatibility (ISI) with bootstrap-based 95 % confidence intervals, and mapped to a three-tier breeding-system classification (self-compatible / partially self-incompatible / self-incompatible).

The analysis provides **phenotypic** evidence for the LEPA breeding-system question, complementary to the **molecular** evidence produced by Phase 4 Step 22 of the SRK amplicon pipeline (per-individual SI-status categorisation from SRK haplotype data).

---

## 2. The ISI index

The Index of Self-Incompatibility, following Zapata & Arroyo (1978) and subsequent botanical convention, is defined as

```
ISI = 1 − mean(S fruit set) / mean(O fruit set)
```

where **S** = fruit set on selfed flowers and **O** = fruit set on outcrossed flowers. Interpretation:

- ISI = 0 → selfed and outcrossed give equal fruit set (fully self-compatible)
- ISI = 1 → selfed gives zero fruit set (fully self-incompatible)
- Intermediate values → partial self-incompatibility

The metric is a **ratio of population means**, not a per-plant statistic.

### Classification thresholds

| ISI range | Category |
|---|---|
| ISI < 0.2 | Self-compatible |
| 0.2 ≤ ISI < 0.8 | Partially self-incompatible |
| ISI ≥ 0.8 | Self-incompatible |

These thresholds match the shaded bands (green / amber / red) in the output figure.

---

## 3. Bootstrap procedure

The point estimate `ISI = 1 − mean(S) / mean(O)` provides no measure of uncertainty. A **nonparametric bootstrap** is used to estimate the sampling distribution of ISI given the observed data.

Each of `B = 10 000` iterations:

1. Draw a resample of the S observations *with replacement*, size `n_S`.
2. Independently draw a resample of the O observations *with replacement*, size `n_O`.
3. Compute mean of each resample.
4. Compute `ISI_i = 1 − mean(S_resample) / mean(O_resample)`.

Summary quantities reported per level (species-wide and per EO):

- **Point ISI** — from raw sample means (no resampling)
- **Median ISI** — median of the 10 000 bootstrap replicates
- **95 % CI** — 2.5 % and 97.5 % percentiles of the bootstrap distribution
- **Prop_SC / Prop_Partial / Prop_SI** — fraction of bootstrap replicates falling in each classification band. These are the probabilities that the true species-level ISI sits in each band given the observed data. **They always sum to 1** across the three bands.
- **Classification** — assigned from the median (not the point) so it aligns with the middle of the bootstrap distribution.

### Why bootstrap of means, not bootstrap of single-pair ratios

An alternative that sometimes gets proposed is to pick one S value and one O value per iteration, compute `1 − S_i / O_j`, and repeat. This measures **variability across individual plant pairings**, not uncertainty in the species-level ISI parameter. The classification thresholds (0.2 / 0.8) are defined on the species-level metric, so the mean-based bootstrap is the one that lines up with how the classification is applied. The single-pair alternative would also be dominated by outliers and produce a jagged, mostly-uninformative distribution given this dataset's zero-heavy S values.

---

## 4. Design choices

Four data-handling choices ship as defaults. Each is a deliberate decision worth documenting.

### 4.1 Matched-sites restriction (default: on)

The analysis restricts to sites (`Site Acronym`) that carry **both** S and O observations. This ensures `mean(S)` and `mean(O)` are computed on the same set of sites — an apples-to-apples comparison rather than a pooled one.

**Why.** In the raw dataset, 5 sites (GC, SE, SF, SG, SR) have S observations but no O controls. If pooled across all sites, the S mean receives contributions from those sites while the O mean does not. Where the S-only sites happen to be systematically low (as they are here, mostly 0–10 % fruit set), pooling drags `mean(S)` down and inflates ISI upward — a real but modest bias (~0.02 in this dataset).

**Effect on this dataset.** Restricting to the 6 matched sites (KB, MHSE, NP, SRC, TMC, WG) drops the point ISI from 0.759 → 0.741 without changing the classification band.

**Reversibility.** Pass `--include-unmatched` to run the pooled analysis instead.

### 4.2 Independent resampling of S and O

S and O are resampled independently at their own sample sizes. This is standard when the two treatments are not paired at the plant level. In this dataset, S and O observations were collected at the slick-spot × treatment level (not on the same individual plant), so independent bootstrap is the correct choice.

If crosses were paired on individual plants, a paired bootstrap (draw the same plant's S and O together) would be more appropriate.

### 4.3 Zero-replacement epsilon (1 × 10⁻⁶)

Zero values in `PercentFruitSet` are replaced with `EPSILON = 1e-6` before any computation. Purpose: prevent `mean(O) = 0` in edge-case bootstrap resamples (which would produce a division by zero). The replacement is small enough that means are shifted by at most 1e-6, well below the precision of the raw data.

### 4.4 Negative-ISI clipping to 0

Bootstrap ISI values below 0 (which occur when `mean(S) > mean(O)` in a resample) are clipped to 0 before quantile computation. Rationale: negative ISI has no biological interpretation in the self-incompatibility framework — it would mean the plant sets *more* fruit when selfed than outcrossed, which is not what the classification bands describe. Clipping keeps the reported CI within the biologically interpretable range [0, 1].

---

## 5. Input data

### File

`Billinge and Robertson outcrossing data_Percent Fruit Set.xls`, located at `/Users/sven/Documents/Current_projects/SRK_bioinformatics/` (one directory above `Canu_amplicon/`).

### Sheet: `outcrossing`

131 rows × 5 columns:

| Column | Description |
|---|---|
| `EO Name` | Element Occurrence label with EO code in parentheses (e.g. `Kuna Butte (EO18-7)`). Populated on the first row of each site; forward-filled on load. |
| `Site Acronym` | Short site code (e.g. `KB`). Used as the site key for the matched-sites restriction. |
| `slick spot` | Slick-spot identifier within a site. Preserved but not used in the current analysis. |
| `treatment` | Treatment code. Values used: `S` (selfing), `AS` (outcrossing — renamed to `O` on load). Other codes (`Bsite`, `Bss`, `NN`, `C`) are filtered out. |
| `percent fruit set` | Response variable, 0 – 100 scale. |

### Observations used after filtering

- Total S+O rows: **44** in the raw file (27 S + 17 O)
- After matched-sites restriction: **37** (20 S + 17 O) across 6 sites

---

## 6. Output

### 6.1 Table: `Tables/ISI_breeding_system_summary.tsv`

Long-format summary, one row per (Level, Group). Species-wide row plus per-EO rows for EOs that meet the minimum-N filter (`n_S ≥ 3` AND `n_O ≥ 3`).

| Column | Description |
|---|---|
| `Level` | `Species` or `EO`. |
| `Group` | Species name for the species row; EO label for EO rows. |
| `N_S` | Number of S observations used. |
| `N_O` | Number of O observations used. |
| `ISI_point` | Point estimate from raw means. |
| `ISI_median` | Median of the bootstrap distribution. |
| `CI_lo`, `CI_hi` | 2.5 % and 97.5 % percentiles of the bootstrap distribution. |
| `Prop_SC` | Fraction of bootstrap replicates with ISI < 0.2. |
| `Prop_Partial` | Fraction with 0.2 ≤ ISI < 0.8. |
| `Prop_SI` | Fraction with ISI ≥ 0.8. |
| `Classification` | Assigned from the median: `Self-compatible`, `Partially self-incompatible`, `Self-incompatible`, or `Undefined`. |
| `N_boot_finite` | Number of bootstrap replicates with finite ISI (usually equal to `B`). |

### 6.2 Figure: `figures/ISI_breeding_system.{pdf,png}`

Single-panel horizontal-density figure:

- **X-axis** — Index of Self-incompatibility (ISI), range 0 – 1
- **Y-axis** — Bootstrap density (from Gaussian KDE of the 10 000 replicates)
- **Filled grey curve** — bootstrap density of the species-wide ISI
- **Solid vertical line** — median ISI
- **Horizontal bar near the baseline** — 95 % confidence interval
- **Three shaded bands with in-band labels**
  - Green (0 – 0.2): "Self-compatible"
  - Amber (0.2 – 0.8): "Partially self-incompatible"
  - Red (0.8 – 1.0): "Self-incompatible"
- **Dotted verticals** at 0.2 and 0.8 mark the classification boundaries
- **Title** — includes point ISI, 95 % CI, classification verdict, and sample sizes

PDF is vector for publication; PNG is 200 dpi raster for slides / previews.

### 6.3 Blank companion figure: `figures/ISI_breeding_system_blank.{pdf,png}`

Same axes, same shaded bands, same in-band labels, same dotted boundary lines — density curve, median line and CI whisker removed. Title reduced to "Species-wide ISI — classification framework" (no ISI value, no CI, no sample sizes). Produced automatically alongside the data figure on every run.

**Use case.** Show the blank version first as a "predictions" slide in talks: audience sees the framework (three bands + thresholds) before the data are revealed, so they understand what they are about to compare against. Reveal the data version on the next slide. Follows the same convention as `SRK_P_compat_traffic_light_EO_blank.{png,pdf}` in the LEPA TP1 figure suite.

---

## 7. Usage

### Default run (matched sites, B = 10 000)

```bash
cd /Users/sven/Documents/Current_projects/SRK_bioinformatics/Canu_amplicon
python3 analyze_ISI_breeding_system.py
```

Produces:

- `Tables/ISI_breeding_system_summary.tsv`
- `figures/ISI_breeding_system.pdf`
- `figures/ISI_breeding_system.png`

and prints the summary table + list of matched vs dropped sites to stdout.

### Command-line options

| Flag | Default | Purpose |
|---|---|---|
| `--input PATH` | Billinge & Robertson XLS in `SRK_bioinformatics/` | Alternative input file. |
| `--tables-dir DIR` | `Tables` | Output directory for the summary TSV. |
| `--figures-dir DIR` | `figures` | Output directory for the PDF + PNG figures. |
| `-B N`, `--bootstrap-samples N` | `10000` | Number of bootstrap iterations. 10 000 is the recommended sweet spot; smaller values produce unstable CI edges, larger values give diminishing returns given the sample size. |
| `--min-n N` | `3` | Minimum `n_S` AND `n_O` per EO required to include that EO in per-EO output. EOs below this threshold are dropped from the summary. |
| `--seed N` | `20260624` | RNG seed for deterministic reproducibility. Change to see how much bootstrap CIs jitter between reruns. |
| `--include-unmatched` | off | When set, includes S observations from sites that have no O control. Produces the pooled analysis (n_S=27 instead of 20 for this dataset). Use for sensitivity checks; matched sites is the recommended default. |

### Reproducibility

Same script + same seed + same input file → byte-identical outputs. The seed defaults to `20260624`.

---

## 8. Interpretation

### Species-level (current dataset, matched sites, 2026-06-24)

- Point ISI = **0.74**, 95 % CI = **0.54 – 0.90**
- Classification: **Partially self-incompatible** (by the median rule)
- **Prop_SI = 27 %** — the probability that the true species-level ISI sits in the SI band given these data
- **Prop_Partial = 73 %**, Prop_SC = 0 %

**Biological reading.** LEPA sits close to the Partial ↔ SI boundary. The point estimate is 0.06 below the SI threshold, but the upper 95 % CI (0.90) crosses into SI territory and 27 % of the bootstrap posterior lands there. This is consistent with the molecular Q2 finding of Phase 4 Step 22 (~94 % of the robust subset genotypically SI). The species is functionally SI in the majority of plants, with enough phenotypic variability (from real biology + small-sample noise) that the phenotypic classification hedges to "Partially self-incompatible."

### Per-EO

Three EOs meet the `min_n ≥ 3` per-treatment filter. Only Kuna Butte (EO18-7, `n_S=5, n_O=5`) gives an unambiguous per-EO call (`ISI = 0.95`, `CI = 0.89 – 1.00`, formally SI). The other two (Tenmile Creek EO32; New Plymouth EO66) have wide CIs and shouldn't drive classification claims on their own.

---

## 9. Caveats and limitations

- **Small overall sample.** With n_S = 20 and n_O = 17 across 6 sites, the bootstrap CI (~0.36 wide) is dominated by the smaller O sample. Additional outcross data would tighten the CI more efficiently than additional S data.
- **Independent bootstrap assumes exchangeability within treatment.** If some sites systematically produce lower fruit set for reasons unrelated to self-compatibility (e.g. pollinator scarcity, drought), the bootstrap CI understates the true uncertainty. A hierarchical bootstrap that resamples sites first, then observations within each drawn site, is a more conservative alternative — not currently implemented, but straightforward to add if required.
- **Point estimate uses raw means; the classification uses the bootstrap median.** These differ by ~0.005 in this dataset (0.741 vs 0.746) — the discrepancy is negligible here but can be non-trivial in strongly skewed bootstrap distributions.
- **Zero replacement (1e-6).** Numerically defensible; scientifically it means "we treat a reported fruit set of 0 as effectively zero, not exactly zero, so the ratio arithmetic stays defined."
- **Negative-ISI clipping.** This makes the reported CIs conservative for populations where mean(S) approaches mean(O). Without clipping, one CI in this dataset (New Plymouth's per-EO row) extended to −0.85 in the unclipped bootstrap. That value has no biological meaning under the ISI framework — clipping to 0 reflects the interpretive convention rather than a data manipulation.
- **Site coding.** The XLS uses `Site Acronym` (not `EO Name`) as the matching key for the matched-sites filter. If the raw data ever gets restructured so that one site acronym maps to multiple EOs (or vice versa), this filter needs to be revisited.

---

## 10. Dependencies

Python 3.10+ with:

- `numpy`
- `pandas`
- `matplotlib`
- `scipy` (for `scipy.stats.gaussian_kde`)
- `xlrd` ≥ 2.0.1 (for reading the `.xls` input; not needed if the file is re-saved as `.xlsx`)

Install with:

```bash
python3 -m pip install --user --break-system-packages numpy pandas matplotlib scipy 'xlrd>=2.0.1'
```

---

## 11. File locations

| Path | Role |
|---|---|
| `Canu_amplicon/analyze_ISI_breeding_system.py` | Analysis script |
| `Canu_amplicon/analyze_ISI_breeding_system.md` | This documentation |
| `SRK_bioinformatics/Billinge and Robertson outcrossing data_Percent Fruit Set.xls` | Input data |
| `Canu_amplicon/Tables/ISI_breeding_system_summary.tsv` | Summary output |
| `Canu_amplicon/figures/ISI_breeding_system.{pdf,png}` | Figure output (data) |
| `Canu_amplicon/figures/ISI_breeding_system_blank.{pdf,png}` | Blank companion figure for talks |

---

## 12. Reference

Zapata, T. R. & Arroyo, M. T. K. (1978). Plant reproductive ecology of a secondary deciduous tropical forest in Venezuela. *Biotropica* 10(3): 221–230.
