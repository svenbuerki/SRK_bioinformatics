# Phase 5 — Mate limitation and fragmentation in *Lepidium papilliferum*

**A compact overview for collaborators.** Framed on a single causal
chain (**fragmentation → genetic drift → mate limitation**); each
step below is presented as Question → Approach → Result. Full
technical detail (formulas, code, output filenames) lives in the long
companion doc [`Phase5_SRK_sampling_and_prediction.md`](Phase5_SRK_sampling_and_prediction.md).

---

## The central hypothesis — a causal chain

*Lepidium papilliferum* (slickspot peppergrass, LEPA) is a tetraploid
Brassicaceae with **sporophytic self-incompatibility (SI)** — a plant
rejects pollen carrying any of the SRK alleles the pollen parent
expresses on its own stigma. In small, patchy populations we
hypothesise a single causal chain that drives reproductive failure:

> **Fragmentation → genetic drift → mate limitation.**
>
> Habitat fragmentation shrinks the effective mating pool at each
> location (fewer plants within ~50 m pollinator flight of each
> other). Small effective mating pools intensify genetic drift on
> SRK, which erodes local SRK diversity and skews local Fg
> composition. The eroded and skewed local pool reduces the fraction
> of pollen a mother is compatible with — **mate limitation** — which
> reduces per-mother seed set.

The pipeline evaluates the chain in order — **fragmentation first**,
then its drift consequences, then its mate-limitation consequences —
because fragmentation is the physical driver upstream of everything
else. It is where location size enters the model (via
`N_fertile_effective`); a census of 1 000 plants split across
20 disconnected slickspots behaves like a location of ~50 plants for
every downstream calculation.

| Link | What we measure | Steps that produce it |
|---|---|---|
| **1. Fragmentation** | Spatial isolation of adult plants; effective mating pool size per location | Step 29b, Step 30b (§ A.5, § B.2 of long doc) |
| **2. Genetic drift on SRK** | Predicted local Fg pool size and Fg frequency composition, driven by `N_fertile_effective` | Step 30 Phase A (§ A.6) |
| **3. Mate limitation** | Predicted per-location random-mating pollen compatibility `P_compat` under sporophytic Class I / II SI | Step 30 Phase A (§ A.7, § A.8); Step 30c empirical validation on adult SRK genotypes (§ C.0) |
| **4. Reduced seed set** | Observed per-mother seed set regressed on predicted `P_compat` | Step 30 Phase B mate-limitation regression (§ C.1) — needs seed genotypes |

**Not part of this framework.** Detecting cases where SI has broken
down entirely (self-compatibility escape) is *not* an outcome of
Phase 5. Such individuals are identified during the Canu-amplicon
SRK genotyping (Phase 4 Step 22b) and enter Phase 5 as prior
information via the empirical zygosity distribution, not as an
experimental target.

---

## Scope and data

- **Species.** Tetraploid (2n = 4x). Every plant carries 4 SRK
  allele copies. Sporophytic Class I / Class II SI: Class I strictly
  dominant within a plant; between-class crosses always compatible;
  within-class crosses compatible only when expressed sets are disjoint.
- **Zygosity.** Canu-amplicon Step 23 shows ~66 % of adult plants
  carry a single distinct functional SRK identity (homozygote-like),
  32 % carry two, ~2 % carry three. This empirical distribution is
  the mother- and father-drawing prior in every P_compat calculation.
- **Data source.** `LEPA_SQL.db` — 2025 wild in-situ occurrences,
  filtered to records with coordinates inside LEPA's Idaho range.
  **39 locations, 3 140 events (individual slickspots), 765 mothers
  with seed records.**
- **Species-wide SRK prior (P1).** 32 functional groups (Fgs) built
  from the Canu-amplicon L1 carrier inventory (49 alleles collapsed
  into 32 Fgs). Dominant Fg = FG001 (41 %); Class I = 6 Fgs (~65 %),
  Class II = 26 Fgs (~35 %).

---

## Step 28 — How many seeds per mother?

**Question.** How many seeds must we genotype from each mother to
observe most of the SRK alleles she was actually crossed with, given
that she can only be pollinated by adults within ~50 m?

**Approach.** Each seed contributes 2 paternal SRK alleles under
tetraploid sporophytic inheritance. Under a uniform-pollen coupon
collector, the probability of missing any one of K local paternal
Fgs after 2·n seed alleles is `(1 − 1/K)^(2n)`. Setting
`(1 − 1/K)^(2n) = 0.10` gives the seeds needed for a 90 % chance to
see every Fg. Rule 2 caps this at **n = 15 seeds/mother** — the point
where the return per additional seed is negligible.

**Result.**

| Reachable donor plants around a mother | K (local Fgs) | Seeds/mother for 90 % detection | Delivered at Rule 2 (n = 15) |
|---|---|---|---|
| 1–2 | ~4 | 3 | 4.0 of 4 (100 %) |
| 3–5 | ~12 | 9 | 11.1 of 12 (93 %) |
| 6–10 | ~28 | 18 | 18.6 of 28 (66 %) |
| 11–20 | ~32 | 20 | 19.7 of 32 (62 %) |
| > 50 | 32 (species ceiling) | 24 | 19.7 of 32 (62 %) |

A single mother sitting in a large event cannot saturate the local
32-Fg pool from 15 seeds alone — but the shortfall closes cleanly
when several mothers' seed lots aggregate at the location scale
(Step 29).

![Figure 1 — Step 28 coverage curves. Panel A: per-mother detection of the local Fg pool as a function of seeds genotyped, one curve per event-size bin; the vertical red line at 15 marks the tetraploid Rule 2 cap. Panel B: aggregation of coverage across mothers at a location (15 seeds × M). Every event size reaches its local ceiling by 5 mothers.](figures/Phase5/step28_coverage_curves.png)

---

## Step 29 — How many mothers per location?

**Question.** How many mothers must we sample at each location to
observe every SRK allele physically present in the location's mating
pool?

**Approach.** Under tetraploid sampling each adult contributes
`PLOIDY × N_fertile = 4·N` allele copies to the location's pool. Set
the coupon-collector target: 90 % chance of observing every Fg in the
local pool at the total delivered allele draws
`A_delivered = M × (4 + 2 × 15) = 34·M`. Add a private-allele floor:
**at least one mother per event** (an isolated slickspot's private
Fg cannot be recovered from any other event).

**Result.** 39 of 39 locations reach 100 % local coverage under the
recommended M. Distribution of the field-team recipe:

- **Green tier (≤ 15 seeds/mother):** 19 / 39 locations.
- **Amber tier (16–100 seeds/mother):** 16 / 39 locations.
- **Red tier (> 100 seeds/mother, unrealistic):** 4 / 39 locations
  — all single-plant or two-plant slickspots in the EO24 group,
  where the census is simply too small to characterise fully.

---

## Step 29b — Which adults actually pollinate which?

**Question.** How do we count the effective mating pool at each
location, given that a pollinator's flight radius is finite?

**Approach.** Build a within-location graph in which two adults are
connected if they are within `R_primary = 50 m` of each other; extract
the largest connected component; define
`N_fertile_effective_50m = N_census × (largest-component share)`. The
50 m primary radius is validated in the sensitivity sweep (Figure 2):
connectivity plateaus at ≥ 75 m; sampling cost stabilises at 75–100 m;
predicted P_compat is radius-independent under empirical zygosity.

**Result.** Median location retains **77 % of its census adults** in
the largest 50 m connected component. Some locations drop to 20 %,
because their census is spread across several slickspots more than
50 m apart. This connectivity factor is what pulls `N_fertile`
downstream in the P_compat and fragmentation calculations.

![Figure 2 — Pollinator-radius sensitivity sweep. Four panels showing how connectivity, fragmentation-aware sampling cost, predicted P_compat, and the sustainable-band fraction of locations change across radii from 10 to 200 m. Connectivity plateaus at ≥ 75 m and P_compat is radius-independent under empirical zygosity, justifying 50 m as the primary radius.](figures/Phase5/step29_location_connectivity_radius_sensitivity.png)

---

## Step 29c — Fragmentation-aware sampling

**Question.** If a location is fragmented into several disconnected
50 m mating pools, how many mothers do we need per pool?

**Approach.** Decompose each location into its 50 m connected
components. For each component apply Step 29's coupon-collector,
plus a maternal floor of ≥ 1 mother per event. Sum the per-component
requirements to obtain the location's honest mother allocation
`M_frag`.

**Result.**

- **43 / 52 locations require MORE mothers** under fragmentation-aware
  allocation than under the pooled Step-29 recommendation.
- **9 / 52 locations unchanged** — single 50 m component already
  covered.
- **Total effort scales from 748 → 1 712 mothers** (2.3× increase).

This is the most consequential design shift: the old pooled
allocation systematically under-sampled fragmented locations.
`M_frag` is now the authoritative field-team recipe (see
[`step29c_sampling_frag_aware_per_event.tsv`](tables/Phase5/step29c_sampling_frag_aware_per_event.tsv)).

---

## Step 30 Phase A — Per-location predictions

### 30.1 Predicted SRK diversity per location — unbiased truth vs sampling

**Two questions, kept strictly separate.** SRK diversity at a
location has two very different meanings and the pipeline predicts
both:

1. **What Nature actually holds at this location** (unbiased truth) —
   how many SRK alleles are physically present. Depends only on
   `N_fertile_effective`; not conditioned on sampling. This is the
   Link 2 output of the causal chain and the raw material on which
   the P_compat prediction (§ 30.2) operates.
2. **What our seed sampling will recover** (sampling-inferred) —
   expected detection under M mothers × 15 seeds each. Under
   tetraploid Rule 2 each mother's seed lot contributes 4 maternal +
   2·15 = 30 paternal allele observations. **Every seed's 2 paternal
   alleles are direct samples of the local pollen donor pool** — the
   seed genotyping is literally a pollen-pool characterisation
   experiment.

**Approach.** For each location: (i) simulate the true local pool by
drawing `PLOIDY × N_fertile_effective` alleles from P1 — this gives
the unbiased local Fg diversity and its 95 % credible interval;
(ii) evaluate the coupon-collector expected detection under
`A_delivered = 34 × M` draws from that local pool.

**Result.**

| Location | Effective N (50 m) | M mothers | Nature holds | Sampling detects | Coverage |
|---|---|---|---|---|---|
| EO30-1 | 420 | 23 | ~32 | ~28 | **89 %** ⚠ |
| EO29 | 417 | 22 | ~32 | ~28 | **89 %** ⚠ |
| EO76 (largest BL3) | 416 | 62 | ~32 | ~31 | 98 % |
| EO32 (well-sampled BL5) | 327 | 38 | ~31 | ~30 | 95 % |
| EO27-1 (large BL4 pilot) | 116 | 33 | ~27 | ~26 | 98 % |
| EO67 (small BL4 pilot) | 6 | 4 | ~8 | ~8 | ~100 % |
| EO24-2 (1 plant, BL5 tail) | 1 | 1 | ~3 | ~3 | ~100 % |

- **What Nature holds** (Panel A of Figure 3) is the biological
  signal. Large well-buffered locations hold ~31–32 of the 32
  species-wide Fgs; the BL5 tail (EO24 group) holds only ~3–6
  because drift has already collapsed the local pool there. Small
  isolated locations are drift-limited, not sampling-limited.
- **What our sampling detects** (Panel B) tracks Nature closely at
  small and mid-sized locations, where the sampling design is
  comfortably above the coupon-collector threshold for the true
  local pool.
- **Two under-sampled outliers** — **EO30-1 and EO29**, both at
  89 % coverage — are the only locations where the permit-realistic
  M does not clear the 90 % target (Panel C). Adding ~5–10 more
  mothers at either would restore coverage to > 95 %. Everywhere
  else the field-team recipe delivers ≥ 95 % coverage.

![Figure 3 — Predicted SRK allele diversity per LEPA location under the tetraploid P1 finite-population model. **Panel A** — What Nature actually holds (unbiased truth, driven by `N_fertile_effective` alone; feeds the P_compat prediction). **Panel B** — What our sampling detects (M mothers × 15 seeds each; 30 paternal-allele samples of the local pollen donor pool per mother). **Panel C** — Coverage = Panel B ÷ Panel A; dotted line = 90 % target. Two locations (EO30-1, EO29) sit below the target — candidates for adding more mothers. Locations sorted by unbiased pool size; dot colour = Bottleneck Lineage; error bars = 95 % credible interval.](figures/Phase5/step30_A_diversity_unbiased_vs_sampling.png)

### 30.2 Predicted pollen compatibility per location

**Question.** Under the sporophytic Class I / II + empirical-zygosity
model, what fraction of pollen would a mother at each location be
compatible with under random mating?

**Approach.** For each location, on each simulation replicate:
(1) draw local Fg pool from P1; (2) simulate mothers from that pool
under the empirical zygosity distribution; (3) apply the sporophytic
Case-A / Case-B analytical formulas per mother; (4) average across
mothers. Report the posterior mean and 95 % credible interval across
replicates.

**Result.**

- **Species-mean P_compat = 0.68** (sustainable band by construction).
- **Traffic-light bands** (1/3 and 2/3 of species mean): failed
  < 0.23, struggling 0.23–0.45, sustainable ≥ 0.45.
- **BL5 tail (EO24 group)** is the primary conservation concern —
  predicted P_compat 0.42–0.51 with wide credible intervals
  entering the struggling band.
- All other BLs sit in the sustainable band at the mean; the tightest
  credible intervals belong to the largest locations (EO76, EO61).

![Figure 4 — Predicted per-location pollen compatibility under sporophytic Class I / II + empirical LEPA zygosity. One dot per location, error bars = 95 % credible interval, panelled by Bottleneck Lineage. Traffic-light bands: red = failed (< 0.23), amber = struggling (0.23–0.45), green = sustainable (≥ 0.45). Species mean = 0.68. BL5 tail slips into the struggling band; all other locations sit in the sustainable band at the mean.](figures/Phase5/step30_A_prediction_fecundation.png)

### 30.3 Fragmentation index

**Question.** Independent of drift on allele frequencies, how
spatially fragmented is each location's mating environment?

**Approach.** Two purely spatial indices, no allele frequencies:
`F_event = 1 − K_spatial_50m / 32` at the event scale (per slickspot),
and `F_location = 1 − within-location connectivity at 50 m` at the
location scale. Both range 0 (no fragmentation) → 1 (fully isolated).

**Result.** F_event and F_location are strongly correlated at BL2 /
BL5 (uniformly small, isolated events → high fragmentation on both
axes). BL3 and BL1 show the widest spread — some locations have
tight event-scale connectivity but many disconnected components at
the location scale, so their mating environment is layered rather
than uniformly fragmented.

![Figure 5 — Event-scale × location-scale fragmentation scatter, one dot per location, coloured by Bottleneck Lineage. Pure spatial indices — no allele frequencies enter. Diagonal locations have matched fragmentation at both scales; off-diagonal locations reveal layered structure (e.g. tight event-scale connectivity but multiple disconnected components at the location scale).](figures/Phase5/step30_A_fragmentation_index.png)

---

## Step 30 Part C § C.0 — Empirical validation of the P_compat model

**Question.** Before we invest in seed genotyping, does the
sporophytic + empirical-zygosity P_compat model reproduce what we
already observe in the adult population from the Canu-amplicon
pipeline?

**Approach.** For each of the 6 EOs with n ≥ 10 individuals having
functional SRK genotypes (338 individuals total: EO76, EO70, EO27,
EO25, EO18, EO67), compute two per-mother P_compat quantities using
the same set of observed mothers and identical Class I / II rules,
changing only the father-drawing distribution:

- **Observed P_compat** — fathers drawn from **observed local Fg
  frequencies at that EO**.
- **Predicted P_compat** — fathers drawn from the **species-wide P1
  prior** (ignoring per-EO drift).

Comparison isolates the effect of local Fg frequency drift on
random-mating compatibility, holding the mothers themselves fixed.

**Result.**

- **4 of 6 EOs sit on the 1:1 diagonal within their 95 % CI** —
  the model reproduces the data at those EOs, confirming local
  frequency drift is not moving them off species-wide expectation.
- **EO70 is a striking outlier**: observed 0.41 (struggling band) vs
  predicted 0.60 (sustainable). This is a genuine biological signal —
  local Fg pool skew makes mothers overlap far more with neighbouring
  fathers than P1 would predict. Not a model failure.
- **EO18 and EO67** show milder deviations in the same direction
  (observed ~ 0.58 vs predicted ~ 0.66). Both are candidates for
  future conservation intervention.
- **Observed zygosity distributions per EO** track the species-wide
  66 / 32 / 2 % well, except EO25 (49 % 2-distinct vs 32 %
  species-wide) — consistent with its perfect-match position on
  the diagonal.

![Figure 6 — EO-level empirical validation of the sporophytic P_compat model. Panel A: observed vs predicted mean pollen compatibility per EO (n ≥ 10 individuals), with 1:1 diagonal and traffic-light bands. 4/6 EOs sit on the diagonal within their 95 % CI; EO70 is a striking outlier (observed 0.41 in the struggling band vs predicted 0.60 sustainable) — genuine local Fg pool skew, not model failure. Panel B: distinct-identity distribution per EO vs the species-wide reference (66 / 32 / 2 %).](figures/Phase5/step30_C_pcompat_observed_vs_predicted.png)

**Caveat.** P1 was built from these same individuals, so this
comparison does not test absolute calibration but robustness to
per-EO drift. EO-scale is intermediate; the finer per-location test
will happen in Phase B once seed genotypes are available.

---

## Phase B — Closing the causal chain with seed data

**Test 1 — Mate-limitation regression (§ C.1 of the long doc).**
Once seed genotypes are back, per-location observed P_compat is
computed from real mother-father pairings. The regression
`per-mother seed set ~ β₁ · P_compat + β₂ · mating_neighbourhood + …`
decomposes the fragmentation effect into two pathways:

- **β₁ (P_compat effect)** measures the strength of the **complete
  drift → mate-limitation branch** — the reproductive consequence of
  the fragmentation-driven Fg loss and skew that Phase A predicted.
- **β₂ (mating-neighbourhood effect)** measures whether fragmentation
  has **direct effects on seed set that are NOT mediated through Fg
  composition** (e.g. fewer neighbouring adults means fewer pollinator
  visits, independent of which alleles they carry).

Together the two coefficients tell us how much of fragmentation's
reproductive cost flows through drift, and how much bypasses it. A
large β₁ + small β₂ = the drift chain fully explains the fragmentation
effect; a large β₂ + small β₁ = fragmentation reduces seed set through
pollinator-behaviour pathways we haven't modelled; both large = both
channels contribute. All three interpretations can be pre-registered
before seed data arrive.

---

## Recommended field pilot — BL4 (EO67 + EO27-1)

Before running Phase B across all 39 locations, we recommend a
**two-location pilot within Bottleneck Lineage 4 (BL4)** — one small
(**EO67**, ~10 fertile plants, 4 mothers in DB) and one large
(**EO27-1**, ~395 fertile plants, 33 mothers in DB).

**Why BL4.** Spans the full range of LEPA slickspot sizes and
internal connectivity; holds the second-biggest species SRK diversity
share; both pilot locations sit in the same BL so the pilot is a
within-BL contrast free of between-BL confounds.

**What the pilot tests.**

- **EO67** (drift-limited regime) — sporophytic + empirical-zygosity
  P_compat mean 0.67, 95 % CI [0.42, 0.87]. Any Phase B observation
  dropping EO67 out of sustainable = decisive evidence of drift-driven
  mate limitation.
- **EO27-1** (aggregation regime) — P_compat 0.69, tight CI [0.59,
  0.77]. Solid anchor point in the sustainable band; if EO27-1 comes
  in below prediction, something bigger than drift is going on.

**Pilot cost.** ≤ 555 seed genotypes (< 5 % of the full 2025 recipe
under Step 29; ≤ 1 200 seeds and < 5 % under the fragmentation-aware
allocation). Two-location within-BL contrast with a well-anchored
prediction on each end.

---

## Deliverables and where to look

| What | File | Notes |
|---|---|---|
| P1 prior over 32 Fgs | [`tables/Phase5/step30_A_prediction_prior_frequencies.tsv`](tables/Phase5/step30_A_prediction_prior_frequencies.tsv) | Species-wide Fg frequencies + Dirichlet α |
| Predicted per-location diversity | [`tables/Phase5/step30_A_prediction_location_diversity.tsv`](tables/Phase5/step30_A_prediction_location_diversity.tsv) | Species-wide + local coverage, 95 % CI |
| Predicted per-location P_compat | [`tables/Phase5/step30_A_prediction_location_pcompat.tsv`](tables/Phase5/step30_A_prediction_location_pcompat.tsv) | Sporophytic + empirical zygosity |
| Fragmentation indices | [`tables/Phase5/step30_A_fragmentation_per_event.tsv`](tables/Phase5/step30_A_fragmentation_per_event.tsv), [`_per_location.tsv`](tables/Phase5/step30_A_fragmentation_per_location.tsv) | Pure spatial, no allele frequencies |
| Field-team recipe | [`tables/Phase5/step29c_sampling_frag_aware_per_event.tsv`](tables/Phase5/step29c_sampling_frag_aware_per_event.tsv) | `M_frag` per event, authoritative |
| Empirical validation (§ C.0) | [`tables/Phase5/step30_C_pcompat_validation_at_eo.tsv`](tables/Phase5/step30_C_pcompat_validation_at_eo.tsv) | 6 EOs, observed vs predicted |

**Key figures.**

- [`step29_location_connectivity_radius_sensitivity.png`](figures/Phase5/step29_location_connectivity_radius_sensitivity.png)
  — why 50 m is the right primary radius.
- [`step30_A_diversity_unbiased_vs_sampling.png`](figures/Phase5/step30_A_diversity_unbiased_vs_sampling.png)
  — three-panel: what Nature holds (unbiased) vs what our sampling detects vs coverage. Flags the two under-sampled locations (EO30-1, EO29).
- [`step30_A_prediction_diversity.png`](figures/Phase5/step30_A_prediction_diversity.png)
  — legacy species-wide-conditioned view, kept for continuity.
- [`step30_A_prediction_fecundation.png`](figures/Phase5/step30_A_prediction_fecundation.png)
  — predicted P_compat per location, traffic-light bands.
- [`step30_A_fragmentation_index.png`](figures/Phase5/step30_A_fragmentation_index.png)
  — event-scale × location-scale fragmentation scatter.
- [`step30_C_pcompat_observed_vs_predicted.png`](figures/Phase5/step30_C_pcompat_observed_vs_predicted.png)
  — the § C.0 empirical validation.

---

## For technical detail

Formulas, code, output schemas, and design decisions are documented
in the long companion doc [`Phase5_SRK_sampling_and_prediction.md`](Phase5_SRK_sampling_and_prediction.md):

- § A.3 — species-wide prior + two-generation model.
- § A.5 — fragmentation of pollen flow, event vs location scales.
- § A.7 — sporophytic Class I / Class II SI biology.
- § A.8 — finite-population P_compat formulas (Case A / Case B).
- § B.4.2 — fragmentation-aware sampling derivation.
- § C.0 — EO-level empirical validation (the § C.0 above is the
  compact version of this section).
- § C.1 — the mate-limitation regression (Phase B).
