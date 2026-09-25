# Phase 5 — SRK sampling & prediction framework (Steps 28–30)

## Contents

The document is organised in **three parts** that map to the natural
flow of the study — setting up the model, generating predictions, and
testing those predictions against real SRK data:

- [Scientific goals](#scientific-goals) — the four hypotheses this framework tests
- **Part A** — [Foundation: setting up the model](#part-a--foundation-setting-up-the-model) · data scope, population-size definitions, within-location pollen connectivity, and the per-mother / per-location sampling design (Steps 28–29). *Everything a reviewer needs to see BEFORE any modelling begins.*
- **Part B** — [Predictions before seed genotyping (Phase A)](#part-b--predictions-before-seed-genotyping-phase-a) · the two-generation trick, prior structure, finite-population prediction of SRK diversity and pollen compatibility per location, and the diversity × compatibility cross-plot (Step 30 Phase A outputs).
- **Part C** — [Testing predictions with observed SRK data (Phase B)](#part-c--testing-predictions-with-observed-srk-data-phase-b) · mate-limitation regression, SI-escape rate test, script behaviour, and how the two tests together classify locations.
- **[Summary of 2025 Phase A results + Phase B data plan](#summary-of-2025-phase-a-results-and-phase-b-data-plan)** — one-page results synthesis with take-home messages and the exact seed-data specs needed to run Part C. *Read this first if you're short on time.*
- [Output map](#output-map--quick-reference-grouped-by-phase) — filenames organised by phase.

**Naming.** Parts **A / B / C** refer to *sections of this document*.
Output filenames use `step30_A_*` for Phase A prediction artefacts (built
without seed data) and `step30_B_*` for Phase B artefacts (built with
observed seed data); `step30_B_DEMO_*` marks synthetic Phase B for
pipeline validation. The doc's Part A + Part B cover the Phase A
artefacts; Part C covers the Phase B artefacts.

## Scientific goals

The purpose of this framework is not diversity estimation for its own sake.
It is a location-level SRK sampling and inference pipeline that lets us test
four connected hypotheses about the reproductive fate of small, isolated
LEPA populations:

1. **Mate-limitation test.** Do locations with fewer compatible pollen
   donors — driven by small mate pool, skewed Fg frequencies, or both —
   show reduced per-mother seed set? Formally: does per-mother
   `germplasmQuantityEstimate` decline with predicted per-mother
   $P_{\text{compat}} = 1 - f_a - f_b$?
2. **SI-escape test.** Under strict SI, no seed can carry a paternal Fg
   matching either of its mother's Fgs. A location where seeds with
   maternal-matching paternal Fgs appear above the strict-SI null is a
   candidate for partial-SI transition (a breakdown of the SI machinery).
3. **Fragmentation × drift decomposition.** Habitat fragmentation depresses
   K (the pollen-donor Fg pool a mother can access); genetic drift skews
   Fg frequencies at small isolated locations. Both reduce
   $P_{\text{compat}}$ but through different channels; the framework
   separates them by using the spatial K (a fragmentation metric) and the
   posterior Fg-frequency skew (a drift metric) as distinct predictors.
4. **Phenotype cross-validation (deferred).** SRK-based predictions can be
   cross-checked against per-site ISI / fruit set in the Genetic-Rescue-DB
   repository. This adds independent lines of evidence but does not shape
   the sampling design and is not modelled here in v1.

The framework has two design consequences: **(1)** the sampling protocol
must scale with each mother's real mate-availability context (Steps 28–29),
and **(2)** the inference layer must produce testable predictions under
random-mating and strict-SI nulls that can be compared to observed seed
genotypes (Step 30).

Everything in this document is built to serve those four goals; the
sampling design (Part 1) is the vehicle, and the prediction/comparison
(Part 2) is the payload.

---

## Part A — Foundation: setting up the model

Part A is everything a reviewer needs to understand **before any
modelling begins** — the data-scope filters, the population-size
definitions, the biological measurement of within-location pollen
connectivity, and the per-mother / per-location sampling design that
turns real field data into a reproducible seed-genotyping recipe.

### A.1 Two operational phases

The pipeline runs in two operational phases, corresponding to the state of
the seed genotype dataset. **All steps use LEPA_SQL.db and the P1
empirical prior in both phases**; what changes is whether observed seed
genotypes are available:

| Phase | State of data | Uses | Produces |
|---|---|---|---|
| **A — Preliminary** (before SRK genotyping) | Field data only (`Germplasm.germplasmQuantityEstimate`, `Events.organismQuantityFertile`, event coordinates), plus the P1 species-wide prior | Steps 28, 29, and Step 30 in `prediction` mode | Per-mother sampling recipe + per-location seed-count recommendations, **plus** predicted SRK diversity, P_compat, fecundation failure |
| **B — Post-genotyping** (after SRK data are back) | Everything above + observed seed genotypes (from Phase A's sampling) | Step 30 in `comparison` mode + Tests 1 & 2 | Observed vs predicted SRK diversity, mate-limitation regression, SI-escape rate test — all only for the locations that have observed data |

**Filename conventions make the phase — and its data provenance —
unambiguous.** Every Step 30 output uses one of three prefixes:

| Prefix | Meaning | Produced by |
|---|---|---|
| `step30_A_*` | **Phase A** — prior-based prediction; always safe to produce | Default `python step30_...py` |
| `step30_B_*` | **Phase B** — result from *real* observed seed genotypes | `--seed-genotypes real.tsv --mother-genotypes real.tsv` |
| `step30_B_DEMO_*` | **Phase B, synthetic** — simulated seeds for pipeline validation | `--demo` |

Additionally, every Phase B figure produced under `--demo` carries a
diagonal **"DEMO" watermark** and an explicit `[DEMO — synthetic data]`
suffix in its title, so screenshots and slide captures cannot be
mistaken for real analysis. The `--demo` mode remains valuable for
pipeline validation without ever risking file-name collision with real
Phase B deliverables.

Steps 28 and 29 outputs are all Phase A by construction (they do not
depend on seed genotypes) and keep their existing `step28_*` / `step29_*`
names.

If observed data exist for only a few locations, Phase B outputs report
on those few; the Phase A prediction outputs still cover every location
the field data know about.

### A.2 Within-location pollen connectivity (foundational input)

The LEPA "location" is a curated grouping of nearby slickspot events;
it is **not** by itself a mating unit. Whether the events at a location
actually exchange pollen depends on their spatial arrangement relative
to the pollen-flight radius. This is a foundational Phase A input
because **connectivity is a direct predictor of both realised SRK
diversity and random-mating compatibility** at the location scale:

- Fewer connected adults → smaller effective mating pool → more drift →
  fewer distinct SRK alleles.
- Fewer connected adults → higher per-mother probability that the pollen
  she sees carries her own alleles → lower random-mating compatibility.

**How it is computed** ([`step29b_location_connectivity.py`](step29b_location_connectivity.py)).
For each location, we build a graph on its events using haversine
distance and draw an edge between events whose coordinates lie within
R metres. Connected components are found by breadth-first search.
Per-location outputs at R = 10 m, **25 m (primary)**, and 50 m:

- **`connected_share_{R}m`** — fraction of adults in a component of
  more than one event (i.e. that exchange pollen with at least one other
  event under the R-radius assumption).
- **`largest_component_share_{R}m`** — fraction of adults in the
  location's largest connected component.
- **`n_components_{R}m`** — number of disconnected mating units the
  location is broken into.

**How it feeds Phase A prediction.** The per-location `largest_component_share_10m`
is used as an **effective-N multiplier**: the finite-population
compatibility model in Step 30 draws the local pollen pool from
`total_n_fertile × largest_component_share_10m` alleles instead of
`total_n_fertile`, so drift acts on the mating unit that actually
exchanges pollen. The mate-limitation regression in Phase B adds
connectivity as an explicit third predictor alongside pollen
compatibility and mating-neighbourhood size, letting the model
attribute reduced seed set to allele-frequency drift vs pollen-flow
fragmentation as distinct causal channels.

**Current LEPA reading** (2025 field only, 39 locations):

| Pollen-flight radius | Locations ≥ 90 % adults connected | Locations < 50 % |
|---|---|---|
| 10 m (small-bee patch — conservative)  | 0 / 39 | 35 / 39 |
| **25 m (extended foraging — primary)** | 8 / 39 | 21 / 39 |
| 50 m (long-flight — optimistic) | 18 / 39 | 11 / 39 |

At the primary radius (25 m) only **8 / 39 locations** behave as a
single mating unit. This is the biological reason every downstream
Phase A prediction has to run the finite-population model, not a
species-wide random-mating limit — and why the mate-limitation story
is a fragmentation × drift story rather than a pure drift story.

**Why 25 m as the primary?** 10 m is a small-bee foraging patch —
conservative; 50 m is a long-flight assumption — optimistic. The
literature on halictid and small solitary bee foraging in the Snake
River Plain does not pin a single distance for LEPA's pollinators,
so we present **25 m as the intermediate honest default**. All three
radii are computed and stored in the connectivity table; the primary
figure and Step 30's finite-population multiplier both use 25 m.

<a id="fig-1"></a>
![Figure 1: Within-location pollen connectivity across the 39 LEPA locations at the **25 m primary pollen-flight radius**. One bar per location, panelled by Bottleneck Lineage (BL4 orange, BL5 green, BL3 red, BL1 purple, BL2 blue — Set1 palette shared with the LEPA_EO_spatial_clustering project). Bar length = fraction of the location's adults that sit in a connected component containing more than one event. Vertical guides: 50 % (orange dotted) and 90 % (green dotted) thresholds. Row labels list the location code, number of events, and total adult census. **At 25 m only 8 / 39 locations reach 90 % within-location connectivity; 21 / 39 sit below 50 %.** Sensitivity views at 10 m (conservative small-bee patch) and 50 m (optimistic long-flight) are stored alongside as `step29_location_connectivity_10m.png/pdf` and `step29_location_connectivity_50m.png/pdf`. Source: `step29b_location_connectivity.py`. Data: [`step29_location_connectivity.tsv`](tables/Phase5/step29_location_connectivity.tsv).](figures/Phase5/step29_location_connectivity.png)

<a id="fig-2"></a>
![Figure 2: Radius sensitivity of within-location pollen connectivity across all 39 LEPA locations, and the biological justification for adopting **25 m as the primary radius**. Orange bars = fraction of locations reaching ≥ 50 % adults connected at each radius; green bars = fraction reaching ≥ 90 %. At **10 m** (conservative small-bee patch) **0 / 39** locations are fully connected — the framework would flag every LEPA location as fragmented, an over-strong claim. At **50 m** (optimistic long-flight) **46 %** are fully connected — the framework would under-flag fragmentation. **At 25 m, 21 % of locations reach 90 % connectivity and 46 % reach 50 %** — an intermediate reading that is neither too pessimistic nor too optimistic, and that leaves plenty of variation across BLs to detect fragmentation-driven mate limitation. Every downstream analysis in the doc uses this 25 m default; 10 m and 50 m are always computed and stored as sensitivity checks. Source: `step29b_location_connectivity.py`.](figures/Phase5/step29_location_connectivity_radius_sensitivity.png)

### A.3 Census N_fertile vs permit-realistic M sampled

The pipeline carries **two distinct notions of population size** that must
not be conflated in downstream analysis:

- **`total_n_fertile`** — the field census: total number of fertile plants
  counted at a location (or an event). This determines the *pollen SRK
  allele pool* a mother is exposed to: `K = 2 × (N_fertile − 1)`. Pollen
  is contributed by every fertile plant whether or not the seed team was
  allowed to collect from it, so K correctly uses this census number.
- **`M_mothers_in_db`** — the permit-realistic reality: number of
  mothers for which the seed team was actually allowed to collect
  seeds, i.e. the count of `germplasmID` records tied to that location.
  For a large slickspot with 543 fertile plants under a 10 % permit,
  `M_mothers_in_db` is around 55, not 543. This is the number that
  drives Step 30's location-level SRK diversity **prediction** and every
  Phase B output (comparison, Tests 1 & 2).
- **`M_achievable_ceiling`** — the theoretical design ceiling assuming
  the permit allowed sampling every fertile plant. Kept in the output
  tables as a reference upper bound; **not** used by the prediction.

The prior version of the pipeline used `M_achievable_ceiling` in the
prediction, which over-estimated detectable SRK diversity at large
slickspots (e.g. EO61 predicted 29 distinct alleles at `M = 543`, versus
21 at the permit-realistic `M = 55`). The current pipeline uses
`M_mothers_in_db` throughout.

### A.4 Data scope — what enters the pipeline

Every SQL query in Steps 28 and 29 applies two mandatory filters and one
optional filter. **Both scripts (`step28_seed_sampling_per_mother.py` and
`step29_event_location_sampling.py`) accept a `--year YYYY` flag with the
same semantics**; downstream Step 30 inherits whatever Steps 28-29
produced.

- **Mandatory · Wild, field-collected only.** Occurrences and germplasm
  records tied to ex-situ material (nursery accessions, seed-increase
  plots, in-vitro cultures) are excluded from every table and figure. The
  two SQL predicates are:
  - `Germplasm.biologicalStatus = 'Wild'` — keeps 765 records, drops 20
    ex-situ.
  - `Occurrences.provenance = 'in situ'  OR  IS NULL` — keeps 2 419 + 552
    records, drops 808 ex-situ and 1 in-vitro.
  This is not a per-run switch — greenhouse plants have no place in a
  wild-population mate-limitation and SI-escape analysis, so the filter
  is baked in.

- **Mandatory · Coordinates present and inside LEPA's known bounding
  box.** Rows with NULL or free-text `eventDecimalLatitude` /
  `eventDecimalLongitude` are dropped, and remaining rows are clipped to
  30–55 °N × −125 to −100 °E to prevent a mistyped decimal from pushing
  the spatial neighbourhood off-planet.

- **Optional · Year filter (`--year YYYY`).** Restricts to events whose
  `eventDate` (stored as `MM-DD-YYYY`) belongs to the given year. When
  set, it applies to *both* the germplasm records (Step 28) and the
  events feeding the spatial neighbourhood (so K_spatial never mixes
  survey years — see § 1.5). The current DB contains 765 wild
  seed-bearing mothers in 2025 and no other year; running without
  `--year` and with `--year 2025` therefore give identical results
  today, but the flag is what will keep the two seasons cleanly
  separated once 2026 field data arrive.

The framework has three steps, each with a distinct role and a distinct
output family:

| Step | Question | Kind of output |
|---|---|---|
| **Step 28** | How many seeds per mother? | Sampling design (per-mother) |
| **Step 29** | How many mothers per event, how many events per location? | Sampling design (per-event / per-location) |
| **Step 30** | Given priors from the preliminary study, what SRK diversity and fecundation failure should we predict per location — and how do observed seed genotypes compare? | Prediction + comparison |

Steps 28 and 29 are the **sampling-design outputs** — they answer *"how much
seed do we need to collect and genotype, and where?"* Step 30 is the
**prediction + comparison output** — it answers *"what should we expect to
find, and does what we found match?"* We keep these output families in
separate files and separate figure folders so that the sampling protocol can
be handed to the field team without any of the downstream inference machinery
attached to it.

---

### A.5 Sampling design — rationale (why the protocol must scale)

An LEPA slickspot with two flowering plants is a fundamentally different
sampling target than one with a hundred. In the two-plant case, every seed a
mother makes must have been sired by the single other plant, so a handful of
seeds already tells you the whole story. In the hundred-plant case, the pool
of potential pollen SRK alleles is up to two orders of magnitude larger, and
the same sampling effort characterises only a small fraction of it.

Any uniform "genotype *n* seeds per mother" protocol wastes effort in one
regime and underdelivers in the other. Worse, it also fails to be defensible
in the regime that matters most for conservation: **small, isolated slickspots
where mate limitation is the exact mechanism we are trying to detect**. In
those spots, a coarse uniform protocol would either over-sample (irrelevant)
or under-sample (missing the very mothers whose fecundation is at risk).

The design must therefore **scale with the mate-availability context of each
mother**, and by extension of each event and each location. The information
we need for this scaling is already in the LEPA field data:

- `Events.organismQuantityFertile` — number of flowering (potentially
  compatible) plants at each slickspot.
- `Germplasm.germplasmQuantityEstimate` — per-mother seed budget.
- `Locations → Events → Occurrences` — hierarchy needed to roll the
  design up from mother to event to location.

### A.6 Sampling design — methodology (Steps 28 & 29 formulas)

We treat mating as a **coupon-collector problem on SRK alleles**. Under a
random-mating null with equal pollen contribution from every compatible mate,
each seed's paternal SRK allele is one independent multinomial draw from a
pool of

$$K \;=\; 2 \times N_{\text{compatible}}$$

potential pollen-donor alleles (two alleles per diploid mate). We do not yet
apply a self-incompatibility filter because the SRK genotypes of every plant
in every event are not yet known; when they are, `N_compatible` will drop
further and the sample sizes will drop with it.

> **One K, two scopes.** The letter *K* is used only for **the size of
> the pollen-donor SRK-allele pool a mother is exposed to** — a single
> biological quantity. What changes is the *spatial scope* of that pool:
>
> - **K (event-only)** — pollen alleles from her event alone,
>   `2 × (N_fertile − 1)`. Used in the coupon-collector formulas below.
> - **K^(R) (spatial neighbourhood at radius R metres)** — same
>   quantity, extended to include all events within *R* m of hers. Used
>   later as the fragmentation predictor in the Part C mate-limitation
>   regression (a small K^(R) means the mother sees few pollen donors
>   because she is physically isolated, not because her local Fg pool
>   is skewed).
>
> Column names in the TSVs follow the same convention: `K_event`,
> `K_spatial_10m`, `K_spatial_25m`, `K_spatial_50m`. There is no
> separate "fragmentation K" — the fragmentation predictor *is* the
> spatial-scope version of the same K.
> A third symbol — `K_fg = 32` — appears only in Part B and refers to
> the **number of species-wide SRK allele classes** (functional groups)
> in the P1 empirical prior; it is a count of allele identities, not
> a count of alleles in a pool, and is unrelated to K / K^(R).

**Spatial mating-neighbourhood extension.** The event-only pool underestimates
K when a mother's event has other LEPA events physically close by — a small
bee will forage across such events indiscriminately. Using
`eventDecimalLatitude` / `eventDecimalLongitude` from the `Events` table
(haversine distances in metres), we compute for each event and radius R the
number of neighbouring events within R m and sum their `N_fertile`:

$$N_{\text{compatible}}^{(R)} \;=\; \Bigl(\sum_{e \in \text{neighbours}_{\leq R}} N_{\text{fertile}}(e)\Bigr) - 1$$

$$K^{(R)} \;=\; 2 \times N_{\text{compatible}}^{(R)}$$

R = **25 m** is the primary assumption — the intermediate, honest
default between the conservative small-bee patch (10 m) and the
optimistic long-flight assumption (50 m). Two sensitivity radii —
**10 m** and **50 m** — are also computed. In the current LEPA data
235 / 704 events (33 %) have a neighbour within 10 m, 439 / 704
(62 %) within 25 m, 579 / 704 (82 %) within 50 m — so the spatial
extension is a real effect, not a formality.

Per-mother recommendations at each spatial radius are computed alongside the
event-only version and stored as parallel columns
(`K_spatial_{R}m`, `n_expected_cov_90_spatial_{R}m`,
`n_achievable_exp_spatial_{R}m`, `achieved_coverage_exp_spatial_{R}m`).

> **Dual role of the spatial-neighbourhood K^(R).** Beyond sizing the
> sampling design at the mother's own event, K^(R) doubles as the
> **fragmentation predictor** for goal 3 of this framework: a small
> K^(R) at 10 m identifies a mother whose accessible pollen-donor
> pool is constrained by physical isolation. In Part C the same
> quantity enters the mate-limitation regression as a predictor,
> alongside the location's allele-frequency skew, so *fragmentation*
> and *drift* can be separated as distinct causal channels. It is the
> same K as in the coupon-collector formula above, just measured at a
> wider spatial scope.

Two sampling rules are computed side by side. Both come from the same
multinomial identity: an allele of frequency *p* is missed after *n*
independent seed draws with probability $(1-p)^n$.

**Rule 1 — expected coverage (event-size-dependent).**

Choose the smallest *n* such that the expected fraction of the K potential
alleles observed exceeds a target *c* = 0.90:

$$n_{\text{exp}} \;=\; \left\lceil \frac{\log(1-c)}{\log(1 - 1/K)} \right\rceil$$

This is the target when it is *reachable* — typically at small and mid-sized
events.

**Rule 2 — miss-probability guarantee (event-size-independent).**

Choose the smallest *n* such that any allele contributing at least $p_{\min}
= 0.10$ of siring events is detected with probability at least $1 - \alpha =
0.95$:

$$n_{\text{miss}} \;=\; \left\lceil \frac{\log \alpha}{\log(1 - p_{\min})} \right\rceil \;=\; 29$$

This is a **K-independent floor**. It is the honest target when Rule 1 becomes
impractical at very large events, where 90 %-coverage requires many hundreds
of seeds.

Both rules are wrapped with a **simulation check** (multinomial resample of
pollen alleles) so the analytical curve is bracketed by an empirical 95 % CI.

**Accounting for uneven seed production.** Mothers do not all make the same
number of seeds. The recommended *n* above is a *target*; the achievable *n*
for a given mother is capped by her real seed budget *S*
(`germplasmQuantityEstimate`):

$$n_{\text{use}} \;=\; \min(n_{\text{recommended}},\ S)$$

We therefore report, alongside the recommended target, **the coverage each
mother actually delivers with her real budget**:

$$\text{achieved coverage} \;=\; 1 - \left(1 - \tfrac{1}{K}\right)^{n_{\text{use}}}$$

This makes seed-production heterogeneity visible per mother instead of only
flagging shortfalls. A related quantity, useful as a biological ceiling, is
the **expected number of distinct pollen SRK alleles physically present in
her complete seed lot** — an upper bound no amount of genotyping can exceed:

$$\text{expected distinct alleles in seed lot} \;=\; K \times \left[1 - \left(1 - \tfrac{1}{K}\right)^S\right]$$

Both quantities appear as new columns in the per-mother output
(`achieved_coverage_exp`, `achieved_coverage_miss`,
`expected_distinct_alleles_in_seed_lot`).

**Aggregation across mothers rescues the target.** The 29-seed cap
does *not* mean we give up on the 90 %-coverage story. It is a
per-mother cap; when we pool across the M mothers at a location, each
mother contributes **31 allele draws** to the location's total pool
(2 alleles from her own genotype + 29 from her genotyped seeds). The
species-wide 90 %-of-32-SRK-alleles target under the P1 prior requires
**A = 704 allele draws pooled at the location** (§ A.8), so:

$$M \;\geq\; \frac{704}{31} \;\approx\; 23 \text{ mothers per location}$$

is enough — under the 29-seed-per-mother cap — to reach 90 % of the
species-wide SRK allele diversity at that location. This is why a
large slickspot like EO76 (62 mothers × 31 = 1 922 allele draws) far
exceeds the target even though each individual mother's 29-seed
sample only reveals ~11 % of her own local pollen pool. **Pooling
mother-level samples turns per-mother partial coverage into
location-level near-complete coverage of the species SRK diversity.**

Two levels of the same coverage story, made explicit:

| View | Question answered | Metric | What 29 seeds/mother buys |
|---|---|---|---|
| **Step 28 (per-mother, Figure below)** | How much of ONE mother's local pollen pool do her 29 seeds reveal? | Fraction of the mother's K pollen-donor alleles detected | 100 % at 2-plant events → 11 % at > 50-plant events |
| **Step 29 (per-location, § A.8)** | How much of the SPECIES-WIDE 32-SRK-allele pool does the location as a whole reveal after pooling across mothers? | Fraction of the 32 Fgs observed at the location | ≥ 90 % once M × 31 ≥ 704 — i.e. from ~23 sampled mothers upward |

Under the 2025 field data **15 / 39 locations already clear the
90 %-species target with the current 29-seed cap + realised M**; the
remaining shortfall is at the small slickspots where no per-mother
effort can compensate for having only 1–10 mothers (§ A.8, R.2).

![Figure — two-panel SRK allele detection on the **absolute-allele scale**, with the local pool capped at the species-wide ceiling of 32 Fgs. Shaded bands = 95 % simulation CI. **Panel A — per mother.** x = seeds genotyped, one curve per event-size bin; the vertical red line at 29 marks the Rule 2 operational cap (never ask any single mother for more than 29 seeds). The 29-seed dots show what each event size **delivers per mother**: 2.0 of 2 (1–2 plants), 6.0 of 6 (3–5), 12.4 of 14 (6–10), 18.8 of 30 (11–20), 19.3 of 32 (21–50), 19.3 of 32 (>50). One mother's 29 seeds cannot saturate a 32-allele pool — this is the coupon-collector limit for a single sampler, not undersampling. **Panel B — aggregation across mothers at a location.** x = number of mothers sampled at the location (29 seeds each = 31 allele draws per mother: 2 maternal + 29 paternal). At the 5-mother benchmark (green dotted line, 145 cumulative seeds): **2.0 of 2**, **6.0 of 6**, **14.0 of 14**, **29.8 of 30**, **31.8 of 32**, **31.8 of 32** — every event size reaches its local ceiling. Because Panel B is a coupon-collector simulation continuous in M, any real location can read off its own coverage by locating (its event-size bin, its actual M) on the correct curve. The "gap at large events" in Panel A closes cleanly at the location scale — see § A.8. Source: `step28_seed_sampling_per_mother.py`. Data: [`step28_coverage_curves_by_Nfertile.tsv`](tables/Phase5/step28_coverage_curves_by_Nfertile.tsv), [`step28_aggregation_curves_by_Nfertile.tsv`](tables/Phase5/step28_aggregation_curves_by_Nfertile.tsv).](figures/Phase5/step28_coverage_curves.png)

### A.7 Step 28 — per-mother sampling protocol

**Script.** [`step28_seed_sampling_per_mother.py`](step28_seed_sampling_per_mother.py)

**Input.** `LEPA_SQL.db` (`Germplasm` × `Occurrences` × `Events` × `Taxonomy`)
filtered to LEPA and to mothers with a seed-count estimate and a known
`organismQuantityFertile`. In addition, the *Data scope* filters above are
applied by default (wild-only, field-collected, coordinates present); the
optional `--year YYYY` flag restricts to a single survey year.

**Sampling outputs** (kept together in `tables/Phase5/` and
`figures/Phase5/`):

- `step28_seed_sampling_per_mother.tsv` — one row per mother: `n_fertile`,
  `n_compatible`, `K_potential_sire_alleles`, `seeds_est`,
  `n_expected_cov_90`, `n_miss_prob_10pct_95`, `n_achievable_exp`,
  `n_achievable_miss`, boolean shortfall flags, plus the
  seed-production-aware columns `achieved_coverage_exp`,
  `achieved_coverage_miss`, and `expected_distinct_alleles_in_seed_lot`.
  This is the field team's per-mother recipe.
- `step28_coverage_curves_by_Nfertile.tsv` — analytical curves + simulation
  CIs binned by event size.
- `step28_coverage_curves.pdf/png` — one line per event-size bin, target
  band at 90 % coverage.
- `step28_per_mother_budget.pdf/png` — per-mother scatter of *recommended n*
  vs *seed budget*, coloured by event-size bin. Points below the diagonal =
  budget-limited mothers.

### A.8 Step 29 — event / location rollup

**Script.** [`step29_event_location_sampling.py`](step29_event_location_sampling.py)

**Rationale.** Step 28 tells us *how many seeds per mother*. It does not tell
us *how many mothers to sample per event* or *how many events to sample per
location*. Without those two extra layers, per-mother numbers cannot be
aggregated into a location-level estimate with a stated coverage guarantee.

**Method — uniform-prior version.** We apply the same coupon-collector logic
one level up: instead of sampling pollen SRK alleles, we sample **maternal
SRK alleles** across mothers and events. The pool has size K = 2 × N_fertile
at the event level and K_loc = 2 × ΣN_fertile at the location level. Under
a uniform-per-plant allele draw, the expected coverage after M mothers
sampled (2M allele draws) is $1 - (1 - 1/K)^{2M}$.

For location-level design we allocate mother slots to events **proportionally
to N_fertile** (larger events get more sampled mothers), capped at each
event's actual N_fertile.

**Method — empirical-prior version (P1 species-wide).** The uniform pool is a
simplification: LEPA has a very skewed Fg frequency distribution (FG001
alone = 40.7 %). Under a skewed prior, common Fgs saturate quickly but rare
Fgs need much more sampling. We use the same species-wide P1 prior as
Step 30 and compute the expected fraction of the 32 Fgs observed after A
total allele draws at the location:

$$\text{expected Fg coverage} \;=\; \frac{1}{K_{\text{fg}}}\sum_{j=1}^{K_{\text{fg}}}\left[1 - (1-f_j)^{A}\right]$$

The 90 %-of-Fgs target is the smallest A for which this quantity reaches
0.90; it has no closed form and is solved by bracketed search. Under the
current P1 prior this target is **A = 704 total allele draws / location**.

At the location level, each sampled mother contributes both her own
genotype and the paternal alleles of any seeds we genotype from her:

$$A_{\text{delivered}} \;=\; M \times (2 + n_{\text{seeds}})$$

so the target A can be reached either by more mothers, more seeds per
mother, or both. When M is fixed by field logistics, the seeds-per-mother
needed is

$$n_{\text{seeds}} \;=\; \left\lceil \frac{A_{\text{target}}}{M} - 2 \right\rceil$$

For very small locations (≤ 3 plants total) this equation returns
unrealistic seeds-per-mother values (e.g. 350) — the honest signal that a
location with only a few plants cannot host 32 Fgs and therefore cannot be
characterised for the species-wide target regardless of seed effort.

**Sampling outputs** (kept in `tables/Phase5/` and `figures/Phase5/`):

- `step29_sampling_per_event.tsv` — per event: `locationID`, `eventID`,
  `n_fertile`, `M_recommended_event`, and the achievable value after capping.
- `step29_sampling_per_location.tsv` — per location. Uniform-prior columns:
  `locationID`, `n_events`, `total_n_fertile`, `M_recommended_location`,
  `M_achievable_location`, `expected_coverage_achieved`. Empirical-prior
  (P1) columns: `A_target_90pct_P1`, `A_delivered_maternal_only`,
  `A_delivered_with_seeds`, `exp_Fg_cov_maternal_only_P1`,
  `exp_Fg_cov_with_seeds_P1`, `seeds_per_mother_for_90pct_P1`.
- `step29_location_coverage_curves.tsv` — analytical curves + simulation CIs
  binned by total location size.
- `step29_location_coverage_curves.pdf/png` — coverage vs mothers-sampled,
  one line per location-size bin.

The outputs stay in the **sampling** family. Nothing here yet uses seed
genotype data.

### A.9 Two-year design

LEPA fieldwork spans multiple years and the LEPA DB indexes each event by
its `eventDate`. Every input to Steps 28–29 (`n_fertile`, `seeds_est`,
`M_actual`, spatial neighbourhood) is a per-year quantity: adult plants,
seed yields, and mating context all change between seasons. The correct
procedure is therefore to **run the sampling design once per year**, using
the same formulas — no methodological change is required.

Two implementation rules keep the two-year design honest:

- **Spatial neighbourhood is year-scoped.** When computing K_spatial for a
  mother sampled in year Y, restrict the neighbourhood sum to events of the
  same year Y. This prevents double-counting adults that appear in
  multiple survey years at the same slickspot (which are geographically
  identical and would otherwise inflate K_spatial).
- **Pooling for inference is question-specific** (see Part 2, § 2.2):
  standing SRK diversity per location is pooled across years (a stationary
  quantity); mate-limitation and SI-escape tests are *not* pooled but
  treated as repeated measures per location (with year as a fixed effect
  and location as a random effect).

The `--year YYYY` flag (documented under *Data scope*) is the single
interface for running either Step 28 or Step 29 year by year:

```
python step28_seed_sampling_per_mother.py --year 2025
python step29_event_location_sampling.py --year 2025
```

For the two-year analysis, run each step once per year and roll up in
Step 30 according to the pooling rules above.

---

---

## Part B — Predictions before seed genotyping (Phase A)

Part B turns the foundational inputs of Part A into **quantitative
predictions**: predicted SRK diversity per location, predicted
random-mating pollen compatibility per location, and the cross-plot
that ties them together. All Part B outputs are produced *before* seed
genotyping and use the P1 species-wide prior from the Canu-amplicon
preliminary study. Filenames: `step30_A_*` (Phase A prediction
artefacts).

### B.1 Rationale — the two-generation trick

Once the sampling protocol from Steps 28–29 is executed, we will have a batch
of seed DNA per mother. Each seed is diploid, carrying **one maternal SRK
allele and one paternal SRK allele**. That means one seed lot per mother
recovers, in a single extraction batch, two independent samples:

| Generation | What we recover | How |
|---|---|---|
| **Parental (G0)** | The mother's own SRK genotype | The invariant / 50 %-frequency allele across her sibs |
| **Filial (G1) — as read through paternity** | A sample of the pollen SRK allele pool she was exposed to | The variable allele across her sibs |

Aggregating across mothers of a location gives us **two independent
estimates of the location-level SRK allele frequency spectrum**: the
*maternal* one (who is standing there) and the *paternal* one (who is
actually contributing pollen). Under random mating with panmictic pollen
dispersal, the two spectra are indistinguishable. Departures flag biased
contribution, cryptic SI filtering, or immigrant pollen — themselves useful
signal.

But we do not want to wait until the seed data are in to know what we
expect. The Canu-amplicon preliminary study already characterised 32
functional SRK allele groups (Fgs) in 263 individuals — a strong empirical
**prior**. We can use that prior *now* to publish predicted diversity and
predicted fecundation failure per location, with credible intervals. When the
seed data land, prediction and observation are compared side by side, and
the locations where they disagree are the actionable ones.

### B.2 Prior structure and prediction methodology

We build three nested Dirichlet priors over the 32 Fgs:

| Prior | Source | Use |
|---|---|---|
| **P0 — Uninformative** | Dirichlet(α = 1) over the 32 Fgs | Reference baseline; what "no prior knowledge" looks like |
| **P1 — Species-wide** | 32-Fg carrier counts from `step26i_L1_carrier_inventory.tsv` | Workhorse prior for locations we haven't visited |
| **P2 — Spatially informed** | Fg counts pooled from spatial neighbours (BL / EO) | Location-specific prior when we know the neighbourhood |

Each prior has an *effective sample size* controlling how strongly it
constrains inference; larger ESS means the prior is only shifted by large
seed datasets, smaller ESS means observed genotypes dominate quickly.

**Prediction outputs — before any seed data.**

For each location and prior we compute, using Dirichlet posterior draws:

- **Expected SRK diversity per location** — expected number of distinct
  Fgs seen after 2M allele draws from the prior. Uses `M = M_mothers_in_db`
  (the permit-realistic count of mothers in the seed bank), not the
  design ceiling.
- **Expected P_compat per location — finite-population model.** Rather
  than integrating the mother's compatibility over the species-wide P1
  prior (which gives ~0.63 at every location), we simulate the location's
  finite mating pool: for each replicate we draw 2 × N_fertile alleles
  from P1, compute the realised local Fg frequencies, and evaluate each
  sampled mother's $P_{\text{compat}}(m) = 1 - f_{a_m}^{\text{local}} -
  f_{b_m}^{\text{local}}$ against those local frequencies. Large
  slickspots converge to the P1 expectation with tight CrI; small
  slickspots (N_fertile ≤ 5) drift stochastically — sometimes into
  "struggling" or "failed" traffic-light bands — with wide CrI capturing
  the founder-effect uncertainty. This is the honest Phase A signal of
  fragmentation-driven mate limitation under a random-mating null.
- **Cross-plot: predicted SRK diversity vs P_compat** — the direct
  visualisation of the fragmentation → drift → SRK diversity loss →
  mate limitation chain. Locations sit on a scatter whose reference
  curve is $E[P_{\text{compat}}] = 1 - 2/k_{\text{eff}}$ (uniform-$f$
  reference).
- **BL grouping.** Every per-location prediction figure is panelled by
  Bottleneck Lineage using `BL_ORDER` and `BL_COLORS` from
  [`srk_bl_constants.py`](srk_bl_constants.py) — a Python mirror of the
  R constants file mirrored from `LEPA_EO_spatial_clustering`. This keeps
  ordering and palette in lock-step with every other LEPA figure. The
  `locationCode → BL` mapping is applied through `base_eo()` which
  normalises EO codes to the CSV's convention (e.g. `EO8 → EO08`,
  `EO27RT → EO27`).

These are the **prediction outputs** and are kept in `tables/Phase5/`
under filenames prefixed `step30_prediction_*` and figures under
`figures/Phase5/step30_prediction_*`. They can be produced today, without
any seed data at all.

**Comparison outputs — once seed genotypes exist.**

For each location, using the seed-genotype TSV, we compute:

- Posterior Fg frequencies (prior × observed) with 95 % credible interval.
- Maternal vs paternal spectrum comparison, tested by permutation.
- Predicted vs observed per-mother seed count under the prior's
  $P_{\text{compat}}$; residuals highlight mate-limited mothers.

These are the **comparison outputs** and are kept in `tables/Phase5/`
under filenames prefixed `step30_comparison_*` and figures under
`figures/Phase5/step30_comparison_*`.

The prediction and comparison families **do not share filenames** so nothing
in the sampling protocol is mistaken for a result, and nothing in the
prediction is mistaken for observed data.

---

## Part C — Testing predictions with observed SRK data (Phase B)

Part C describes the two statistical tests and the comparison outputs
that fire once **real seed genotypes** are available. Every Part C
artefact is filename-prefixed `step30_B_*` (real data) or
`step30_B_DEMO_*` (synthetic Phase B for pipeline validation, always
carrying a diagonal DEMO watermark). This part closes the loop on the
four scientific goals of the framework.

### C.1 Test 1 — Mate-limitation regression (goals 1 + 3)

Under strict SI + random mating, a mother of Fg genotype $(a, b)$ has a
predicted proportion of ovules that will encounter compatible pollen equal
to her per-mother compatibility

$$P_{\text{compat}}(m) \;=\; 1 - f_a - f_b$$

where $f_j$ is the Fg frequency at her mating neighbourhood (the location's
posterior for the pooled-year analysis; the year-specific spatial
neighbourhood for the repeated-measures analysis). Her expected seed set
is proportional to $P_{\text{compat}}(m)$:

$$E[\text{seeds}_m] \;\propto\; \text{ovules}_m \times P_{\text{compat}}(m)$$

The **test** is a mixed-effects regression of observed
`germplasmQuantityEstimate` on predicted $P_{\text{compat}}$:

$$\text{seeds}_m \;=\; \beta_0 + \beta_1 \cdot P_{\text{compat}}(m) + \beta_2 \cdot K^{(10\text{m})}(m) + u_{\text{location}(m)} + u_{\text{year}(m)} + \epsilon_m$$

where $K^{(10\text{m})}$ is the mother's mating-neighbourhood
pollen-donor pool at the 10 m radius (§ A.6). Two coefficients, two
distinct causal channels:

- $\beta_1 > 0$ with 95 % CI excluding 0 → **evidence of mate
  limitation driven by allele-frequency drift** (skewed local
  frequencies reduce compatibility).
- $\beta_2 > 0$ conditional on $\beta_1$ → **fragmentation effect
  independent of drift** — spatial isolation reduces seed set above
  and beyond what allele skew alone explains.
- A significant $\beta_1$ with $\beta_2 \approx 0$ → **drift-dominant**
  mate limitation (small locations look fragmented but their
  mating-neighbourhood pool gains them a comparable pollen-donor
  count; the harm is the skewed local frequencies).

The regression uses `location` as a random intercept (to absorb time-
invariant site effects) and `year` as a fixed effect (to control for
across-year climate / phenology). This is where the two-year design
becomes a design advantage rather than a nuisance: within-location
year-to-year variation in mate context is what identifies $\beta_1$
cleanly.

**Predictable output columns** (per location and per mother):
`predicted_P_compat`, `predicted_seed_set`, `beta_1_estimate`,
`beta_2_estimate`, 95 % CI, and a per-mother residual
`obs_minus_pred_seeds`.

### C.2 Test 2 — Self-incompatibility escape rate test (goal 2)

Under **strict SI**, the paternal SRK allele in any seed of mother $(a, b)$
cannot equal $a$ or $b$. The rate of self-matching paternal alleles per
location is therefore expected to be

$$\pi_{\text{self-match}}^{H_0} \;=\; 0$$

Under **partial SI** or SI breakdown, some fraction of seeds carry a
paternal Fg matching one of the mother's Fgs:

$$\hat{\pi}_{\text{self-match}}(\ell) \;=\; \frac{\bigl|\{\text{seeds at loc.}\;\ell : \text{paternal Fg} \in (a_m, b_m)\}\bigr|}{\bigl|\text{seeds at loc.}\;\ell\bigr|}$$

The **test** is a permutation of the strict-SI null: for each location,
permute paternal alleles across seeds while preserving the marginal Fg
frequency vector; the p-value is the fraction of permutations reaching
$\hat{\pi}_{\text{self-match}}$ at least as extreme as observed.

Locations with $\hat{\pi}_{\text{self-match}} > 0$ significantly are
**candidate SI-escape sites**. Their status is then a validated hypothesis
for follow-up phenotyping in Genetic-Rescue-DB (goal 4, deferred).

**Predictable output columns** (per location): `n_seeds_scored`,
`n_self_matching`, `pi_self_match`, `p_permutation`, `q_bh_fdr`.

### C.3 Why the two tests together tell the story

- Locations that pass Test 1 (mate-limited) but pass Test 2 (strict SI
  preserved) → classical small-population reproductive failure,
  compounded by fragmentation and/or drift, but the SI system still works.
- Locations that pass Test 2 (SI-escape) → those are populations where the
  reproductive-assurance response has already kicked in: the SI barrier
  has partially broken down. These are where the ISI phenotypic data in
  Genetic-Rescue-DB should show elevated fruit set relative to their
  SRK-predicted $P_{\text{compat}}$.
- Locations that pass **both** Tests → the highest-priority conservation
  targets: reproductively failing *and* undergoing a mating-system
  transition.

### C.4 Step 30 script behaviour (CLI + phase filenames)

**Script.** [`step30_srk_diversity_prediction_vs_observed.py`](step30_srk_diversity_prediction_vs_observed.py)

The script always produces the prediction family. If a `--seed-genotypes`
TSV is provided (columns `locationID`, `eventID`, `germplasmID`, `seed_id`,
`maternal_Fg`, `paternal_Fg`), it additionally produces the comparison
family. If no such TSV exists but the flag `--demo` is passed, the script
simulates a plausible seed-genotype dataset from the prior itself — useful
to verify the comparison pipeline end-to-end before real seed data arrive.

### C.5 The clean payoff

- **Before seed data are back**: publishable predicted SRK diversity,
  predicted per-mother $P_{\text{compat}}$, and predicted
  fecundation-failure rates per location — with honest credible intervals.
  These can be pre-registered.
- **When seed data land**: the mate-limitation regression (§ 2.2.1) and
  the SI-escape rate test (§ 2.2.2) both fire from the same seed-genotype
  input. Together they classify each location into a 2 × 2 matrix (mate-
  limited yes/no × SI-escaped yes/no) that is directly interpretable as a
  conservation prioritisation.
- **Fragmentation × drift decomposition**: the two coefficients
  $\beta_1$ (pollen-compatibility effect) and $\beta_2$ (mating-neighbourhood-size effect) from the
  mate-limitation regression separate the two mechanisms even though they
  both reduce reproductive success — a distinction unreachable with any
  single-year, single-location analysis.
- **Two-generation efficiency**: every mother's seed lot pays double — it
  certifies her genotype (maternal inventory) and samples her pollen
  environment (paternal inventory) in one experiment.
- **Phenotype cross-validation is a natural next step**: the two-by-two
  matrix from § 2.2.1 + § 2.2.2 predicts what ISI / fruit set from
  Genetic-Rescue-DB should look like at each location, which can be
  overlaid without changing any of the code paths in this framework.

---

---

## Summary of 2025 Phase A results and Phase B data plan

This section collates the key results that live inline in Parts A, B,
and C, with cross-references so a reviewer short on time can start
here and drill into any method that catches their eye. Every result
is a **Phase A prediction** — the seed genotypes that will confirm or
refute them do not yet exist and Part C describes how to generate
them.

Every result below is a **Phase A prediction** — built from the LEPA
field census, the permit-realistic sampling record, and the
Canu-amplicon 32-Fg species-wide prior. No seed genotypes yet exist.
These are pre-registered predictions that **Phase C tests (with real
seed data)** will confirm or refute. Numbers refer to the 2025 field
season with the wild, in-situ data-scope filters: 765 mother plants
across 39 slickspots (5 Bottleneck Lineages) sampled at 5–12 % of the
census on average. Results are cross-linked back to the section that
describes the underlying method.

### R.1 Within-location pollen connectivity is the biological floor

**What we did.** For each location we built a spatial graph of its
events and computed the fraction of adults connected via pollen flow
at three foraging radii (10 m, 25 m, 50 m). Full method: [foundational
subsection](#within-location-pollen-connectivity-foundational-phase-a-input).

**Result.** At the small-bee 10 m assumption, **no LEPA location behaves
as a single mating unit** — 0 / 39 locations reach ≥ 90 % of adults
connected, and 35 / 39 sit below 50 %. At 25 m the picture eases to
8 / 39 fully connected; at 50 m to 18 / 39. Every BL5 tiny slickspot
(EO24 group) is a single-event location with no within-location
connectivity at any radius; even large BL3 locations like EO76
(23 events, 517 adults) achieve only 58 % internal connectivity at 10 m.

**Take-home for reviewers.** Fragmentation is not an abstract threat —
it is measurable at the within-location scale, before any genotyping.
This is the biological reason every downstream Phase A prediction has
to run a finite-population model, not a species-wide random-mating
limit. Connectivity enters the pipeline as an **effective-N multiplier**
for the compatibility model and as an **explicit predictor** in the
Phase B mate-limitation regression.

See [**Figure 1**](#fig-1) (per-location bars at 25 m) and [**Figure 2**](#fig-2)
(radius-sensitivity aggregate) in the foundational subsection above.

### R.2 A sampling protocol scaled to each mother's mate context (see A.5–A.8)

**What we did.** Every mother plant received an individually
calibrated seed-genotyping recommendation based on the number of
fertile plants at her slickspot and the number of seeds she produced.
Two decision rules run side-by-side: a 90 %-coverage target (fires
where the mate pool is small enough to characterise cheaply) and a
29-seed detection floor (guarantees identification of any pollen
allele siring at least 10 % of her offspring).

**Result.** For 44 of 52 LEPA locations the protocol reaches the
90 %-of-species-wide-SRK-alleles target with fewer than 30 seeds per
mother. The remaining 8 locations are tiny slickspots (≤ 10 fertile
plants) where no sampling intensity can compensate for the small
census — a mathematical certainty, not a design failure.

**Take-home for reviewers.** Field sampling effort is auditable at
the per-mother level. The field team receives one number per
germplasmID ([`step29_field_team_sampling_recipe.tsv`](tables/Phase5/step29_field_team_sampling_recipe.tsv)).
The design is
robust to permit-driven under-sampling and honest about what small
populations cannot reveal.

<a id="fig-3"></a>
![Figure 3: Recommended seeds per mother plant for each LEPA location under the P1 empirical prior. One horizontal bar per location, coloured by achievability tier: **green** (≤ 29 seeds/mother, fits the Step 28 Rule 2 miss-probability floor), **amber** (30–100 seeds/mother, achievable with focused effort), **red** (> 100 seeds/mother, unrealistic — the location's census is too small to characterise). Vertical guides: 29 (Rule 2 floor) and 100 (practical ceiling). Bars > 300 are capped for display, with the true value annotated at the right. **44 / 52 locations sit in the green tier; 5 in amber; 3 in red (EO24-1, EO24-2, EO24-7 — all single-plant or two-plant slickspots).** Source: `step29_event_location_sampling.py`.](figures/Phase5/step29_recommended_seeds_per_mother_P1.png)

### How Part B results are computed — short summary for readers

Every Part B result below (R.3 – R.5) comes from the same
**finite-population simulation model**, applied per location, with
the **25 m connectivity radius** as the biological scope of pollen
movement. Four steps:

1. **Simulate the local SRK pool (drift signature).**
   Each location's mating population has *N_fertile* × 25 m-connectivity
   plants ; those contribute 2 × N_fertile allele copies. We **draw
   those alleles from the Canu-amplicon species-wide prior (P1)** —
   the empirical frequency distribution of the 32 SRK allele groups
   across the 263 preliminary genotyped individuals. Small locations
   lose rare alleles to drift; large locations approach the species-wide
   32.

2. **Predicted SRK diversity per location** = expected number of
   *distinct* alleles in that local pool.

3. **Predicted pollen compatibility per location.** Within each
   replicate we simulate *M* sampled mother genotypes (drawn from the
   local pool, not from P1) and compute each mother's random-mating
   compatibility as **1 − f_a − f_b** using her two alleles' *local*
   frequencies. Under strict self-incompatibility, this is the
   expected fraction of her ovules that meet compatible pollen.
   Averaged over the M mothers, then over replicates.

4. **Uncertainty.** Both quantities are recomputed on 1 000 – 4 000
   independent replicates of the local pool. Reported means and 95 %
   credible intervals are the mean and 2.5 % / 97.5 % percentiles
   across replicates. Small locations have wide intervals (founder-
   effect uncertainty is large); large locations have tight intervals
   (local pool converges to P1).

**Sampling effort — what the numbers assume we will see.** *M* mothers
× 29 seeds each = **A = M × 31 allele draws** from the local pool
(2 alleles per mother from her own genotype + 29 paternal alleles per
seed). We report coverage in two flavours:

- **Species-wide coverage** — fraction of the 32 P1 alleles detected.
  Biased low for small locations (drift already removed most).
- **Local coverage** — fraction of the alleles *actually at the
  location* that we detect. The biologically honest metric; ~100 % at
  every location under the current design (see R.3 table).

**One-line take-home:** SRK diversity is what drift has left, pollen
compatibility is how well a random mother matches the neighbours she
can reach at 25 m, and both are simulated per location under the
same finite-population draw from the species-wide prior.

### R.3 Predicted SRK diversity captured per location (Part B result)

**What we did.** Under the species-wide 32-Fg prior, we predicted how
many distinct SRK alleles Phase B seed genotyping should recover at
each location, using the permit-realistic count of mothers actually
sampled (55 at EO61, 62 at EO76, down to 1–4 at the smallest
slickspots). Predictions are panelled by Bottleneck Lineage
(BL4 → BL5 → BL3 → BL1 → BL2).

**Result.** Large, well-buffered slickspots in BL3 (EO76) and BL1
(EO61) are predicted to recover 21–22 distinct alleles out of the
species-wide 32. BL5 is highly bimodal: the well-sampled locations
(EO32, EO18-7 group) reach 15–19 alleles, while the BL5 tail
(EO24, EO24-1, EO24-2, EO24-7) can recover only 2–5 alleles from
their tiny local populations.

**Take-home for reviewers.** The predicted-diversity gradient across
BLs is the biological signal of habitat fragmentation and drift
operating at different intensities across the species. BL5's tail
is where drift has already collapsed local SRK diversity to the
point of predicted mate-limitation — the primary conservation
concern.

**Species-wide vs location-local coverage — two honest metrics.**
The numbers above answer *"how many of the 32 species-wide SRK alleles
will our seed genotyping recover at each location?"* — a metric biased
against small populations because it counts alleles that drift already
removed. The **location-local coverage** metric answers the more
biologically relevant question *"how many of the SRK alleles that are
actually at this location will we detect?"*. Both quantities are
computed for every location (see the new
`predicted_local_pool_size_mean` and `predicted_local_coverage_mean`
columns in [`step30_A_prediction_location_diversity.tsv`](tables/Phase5/step30_A_prediction_location_diversity.tsv)):

| Location | Effective N (25 m) | Local pool (alleles present) | Local coverage at 29 seeds × M | Species coverage (of 32) |
|---|---|---|---|---|
| EO24-2 (1 plant) | 1 | ~2 | **~100 %** | ~6 % |
| EO67 (small pilot) | 6 | ~6 | **~100 %** | ~18 % |
| EO27-1 (large pilot) | 116 | ~22 | **~100 %** | ~55 % |
| EO32 (well-sampled BL5) | 224 | ~26 | ~98 % | ~58 % |
| EO76 (largest BL3) | 416 | ~30 | ~98 % | ~67 % |

**Every LEPA location — including the tiny BL5 slickspots — is
essentially fully characterised at the local level (~98–100 %).** The
species-wide coverage number remains useful as a cataloguing metric,
but for Phase B mate-limitation and SI-escape tests the local coverage
is what matters: we are testing reproductive dynamics on the alleles
that are physically present, not attempting a species-wide inventory.

<a id="fig-4"></a>
![Figure 4: Predicted number of distinct SRK alleles detected per LEPA location under the P1 species-wide prior (32 Fgs from the Canu-amplicon preliminary study). Horizontal layout, one dot per location, panelled by Bottleneck Lineage in the standard BL_ORDER (BL4, BL5, BL3, BL1, BL2). Dot size ∝ √M (permit-realistic count of mothers with seed records in the LEPA DB); error bars = 95 % credible interval from Dirichlet posterior draws. Y-tick labels give `locationCode (n = mothers sampled)`. Vertical dashed line marks the species-wide SRK allele ceiling of 32. **BL3 (EO76, EO38) and BL1 (EO61) predict ~21 distinct alleles; the BL5 tail (EO24 group) predicts 2–5.** Source: `step30_srk_diversity_prediction_vs_observed.py`.](figures/Phase5/step30_A_prediction_diversity.png)

### R.4 Finite-population prediction of pollen compatibility (Part B result)

**What we did.** For each location we simulated the local mating
pool by drawing 2 × N_fertile alleles from the species-wide prior,
computed local allele frequencies, and evaluated the expected
random-mating compatibility of a sampled mother against those local
frequencies. This captures the *founder effect* on small
populations directly: a slickspot of two plants literally has
four SRK alleles and no others.

**Result.** Large locations (N_fertile ≥ 50) converge to the
species-wide expected compatibility of ~0.63 (the "sustainable"
traffic-light band) with tight 95 % credible intervals. Small BL5
slickspots collapse into "struggling" (EO24-1, EO24-7 at ~0.31)
or "failed" (EO24-2 at exactly 0 — a single plant has no
compatible partner). Credible intervals widen accordingly, honestly
reflecting founder-effect uncertainty.

**Take-home for reviewers.** Fragmentation-driven drift is
predicted, quantitatively, to depress pollen compatibility below
the "failed" threshold at ≥ 3 Idaho slickspots before we
even open a seed lot. These are the specific sites where
mate-limitation is the working hypothesis.

<a id="fig-5"></a>
![Figure 5: Predicted per-mother pollen compatibility under random mating for each LEPA location. Horizontal layout, one dot per location, panelled by Bottleneck Lineage. Traffic-light background bands mark **failed** (compatibility < 0.20, red), **struggling** (0.20–0.40, orange) and **sustainable** (≥ 0.40, green) — matching the wording used elsewhere in the SRK random-mating framework. Dot position = mean predicted compatibility from the finite-population model (2 × N_fertile local alleles drawn from P1, N_fertile scaled by the location's within-10 m connectivity share); error bars = 95 % credible interval across simulation replicates; dot size ∝ √M (mothers with seed records in DB). **Large slickspots converge on ~0.63 with tight CI; BL5 tiny slickspots collapse into "struggling" or "failed" bands with wide CI reflecting founder-effect uncertainty.** Source: `step30_srk_diversity_prediction_vs_observed.py`.](figures/Phase5/step30_A_prediction_fecundation.png)

### R.5 Cross-plot: SRK diversity vs pollen compatibility (Part B result)

**What we did.** We plotted per-location predicted SRK allele
diversity against predicted pollen compatibility, with locations
coloured by BL and the theoretical mean-compatibility curve
$1 - 2/k_{\text{eff}}$ overlaid.

**Result.** The BL5 tail sits below the P1 expectation on both
axes — the direct signature of the fragmentation → drift → SRK
diversity loss → mate-limitation chain. Larger slickspots
converge on the P1 mean.

**Take-home for reviewers.** This single figure is the causal-chain
summary of Objective 3: fragmentation and drift are not abstract
threats to LEPA; they translate into a measurable, testable drop
in the number of compatible mates a mother can access.

<a id="fig-6"></a>
![Figure 6: Predicted SRK allele diversity (x) vs predicted random-mating pollen compatibility (y) per LEPA location, coloured by Bottleneck Lineage (Set1 palette: BL1 purple, BL2 blue, BL3 red, BL4 orange, BL5 green). Error bars on both axes come from the same Dirichlet posterior draws that produced the two single-quantity Phase A figures. Dot size ∝ √M (mothers with seed records in DB). Dashed grey reference curve = the theoretical mean-compatibility relation `1 − 2 / k_eff` under uniform allele frequencies — locations sit ON this curve when Fg diversity is even, below it when the local pool is skewed toward one or two common alleles (drift signature). **The BL5 tail (EO24, EO24-1, EO24-2) sits at the extreme low-diversity / low-compatibility corner; large BL3 (EO76) and BL1 (EO61) sit near the top-right where both quantities saturate near P1.** This is the single-figure summary of the fragmentation → drift → SRK diversity loss → mate-limitation chain that anchors Objective 3. Source: `step30_srk_diversity_prediction_vs_observed.py`.](figures/Phase5/step30_A_diversity_vs_pcompat.png)

### R.6 How Part C closes the loop: data needed to test the predictions

Every Part B prediction above becomes testable once we have
**seed-DNA genotypes** at SRK. The two-generation trick makes each
mother's seed lot pay double — the mother's own genotype falls out
of her sib set (constant / 50 %-frequency alleles), and the paternal
alleles of those seeds sample the pollen environment that fertilised
her ovules. One seed lot per mother = one maternal genotype +
n paternal alleles.

**Data-generation pipeline.**

1. **Seed extractions.** For each mother in the field-team recipe
   ([`step29_field_team_sampling_recipe.tsv`](tables/Phase5/step29_field_team_sampling_recipe.tsv)),
   genotype the recommended
   number of seeds (typically 29, capped at her real seed budget).
   Total per-year effort: 765 mothers × ~29 seeds ≈ 22 000 seed
   genotypes.
2. **SRK amplicon sequencing.** Any protocol that resolves the 32
   Fgs. The Canu-amplicon pipeline is already validated on adult
   plants and applies directly to seed material.
3. **Fg assignment.** Each observed SRK allele is mapped to one of
   the 32 Fgs from the preliminary study
   ([[project_functional_srk_definition]]). Alleles that do not map
   are flagged as candidates for post-hoc expansion of the Fg set.
4. **Two output tables**, which the Phase B pipeline consumes
   without any code changes:
   - `real_seed_genotypes.tsv` — one row per seed:
     `locationID`, `eventID`, `germplasmID`, `seed_id`, `maternal_Fg`,
     `paternal_Fg`.
   - `real_mother_genotypes.tsv` — one row per mother:
     `locationID`, `germplasmID`, `mother_Fg_a`, `mother_Fg_b`,
     `K_spatial`, `seeds_est`. Everything except the two `mother_Fg_*`
     columns is already in the Step 28 output; the two Fg columns
     come directly from the inferred maternal genotype in step 1.

**Running Phase B.**

```
python step30_srk_diversity_prediction_vs_observed.py \
    --year 2025 \
    --seed-genotypes  real_seed_genotypes.tsv \
    --mother-genotypes real_mother_genotypes.tsv \
    --match-seed-count
```

This produces the Phase B outputs — `step30_B_*` files, no `DEMO`
watermark — restricted to the locations for which real seed data
exist. The tests that fire automatically:

- **Test 1 — Mate-limitation regression.** Does per-mother
  observed seed count decline with predicted random-mating
  compatibility? A positive slope with 95 % CI excluding zero =
  mate-limitation confirmed at population level.
- **Test 2 — Self-incompatibility escape rate.** Do paternal
  alleles match either of the mother's alleles at rates above
  the strict-SI null? False-discovery-corrected across locations.
  A significantly positive rate at any location = candidate
  partial-SI transition.

Below is what the two tests **will look like** once real seed
genotypes exist — these are DEMO renders produced under
`--demo` (simulated Phase B data), watermarked so they cannot
be confused with real analysis.

<a id="fig-7"></a>
![Figure 7 (DEMO): Mate-limitation regression preview. One dot per LEPA location, colour = "sustainable" band. X = predicted random-mating pollen compatibility (mean across sampled mothers at that location); Y = mean observed seeds per mother. Traffic-light background bands (failed / struggling / sustainable). Dashed line = weighted OLS fit, slope + p-value printed in the legend. In this DEMO the simulator baked in a direct causal link (seed set ∝ compatibility), so the slope is highly significant. **With real data the same figure will test whether observed seed set actually declines with predicted compatibility — a positive slope with 95 % CI excluding 0 confirms mate limitation at population level.** Source: `step30_srk_diversity_prediction_vs_observed.py --demo`.](figures/Phase5/step30_B_DEMO_mate_limitation.png)

<a id="fig-8"></a>
![Figure 8 (DEMO): Self-incompatibility escape preview. One horizontal bar per LEPA location, sorted by observed rate. X = observed rate of pollen alleles matching the mother's own SRK alleles (= self-incompatibility escape rate). Red bars = locations that reject the strict-SI null at 5 % false-discovery rate; grey bars = consistent with strict SI. In this DEMO the simulator baked in 8 % SI escape rate, so most locations show detectable escape. **With real data any red bar names a candidate partial-SI population — a location where the SI machinery has broken down enough that self-pollen produces seeds.** Source: `step30_srk_diversity_prediction_vs_observed.py --demo`.](figures/Phase5/step30_B_DEMO_si_escape.png)

**Two-year extension.** The `--year` flag on Steps 28, 29 and 30
makes the pipeline year-scoped. When 2026 field data arrives, run
each step twice (once per year); the pooling rules under § 1.5
tell Phase B how to combine years for standing-diversity
inference vs how to keep them separate as repeated measures for
mate-limitation and SI-escape tests.

**Downstream: cross-validation with phenotype.** The mate-limitation
and SI-escape results define a 2 × 2 classification per location
(mate-limited yes/no × SI-escaped yes/no). Phase III of the wider
project overlays the Genetic-Rescue-DB ISI / fruit-set phenotype
against these predictions to provide an independent line of evidence.

### R.7 BL4 pilot study — one small + one large location

Before running Phase B across all 39 locations, we recommend a
**pilot within Bottleneck Lineage 4 (BL4)** using two contrasting
locations. BL4 is a good choice because it spans the whole range of
LEPA slickspot sizes and internal connectivity, holds the second-
biggest species SRK diversity share, and is not the SI-escape hot
spot (BL3), so any Phase B signal in BL4 is a mate-limitation
signal by construction.

**The two pilot locations:**

| Attribute | **EO67** (pilot — small) | **EO27-1** (pilot — large) |
|---|---|---|
| Number of events (slickspots at this location) | 2 | 8 |
| Total fertile adults (census) | 10 | 395 |
| Effective N under 25 m connectivity | 6 | 116 |
| Mothers with seed records in DB | **4** | **33** |
| Seeds / mother (Rule 2 cap) | 29 (budget-limited if lower) | 29 |
| Total seeds to genotype at this location | ≤ 124 allele draws | 1 023 allele draws |
| **Predicted local pool size** (SRK alleles at this location) | ~6 alleles | ~22 alleles |
| **Predicted local coverage** (of the alleles at this location) | **~100 %** ✓ | **~100 %** ✓ |
| Predicted species-wide coverage (of 32 P1 alleles — biased against drifted-out alleles) | ~18 % | ~55 % |
| Predicted random-mating compatibility | *"struggling"* band | *"sustainable"* band |

**Reading the two coverage numbers.**
The **local coverage** (~100 % at both pilot locations) is the biologically
honest number: our seed genotyping will characterise essentially every
SRK allele that is *physically at* EO67 and EO27-1. The **species-wide
coverage** is systematically lower at small locations not because we
are under-sampling but because drift has already removed most species
alleles from those locations. For a mate-limitation pilot, the local
coverage number is the relevant one — we are testing reproductive
dynamics on the alleles that are present, not attempting to inventory
the species pool.

**Why this pair is well-chosen for a pilot.** EO67 tests the
**budget-limited / founder-effect regime** (small population, few
mothers, low predicted compatibility, high risk of mate limitation).
EO27-1 tests the **aggregation regime** (many mothers, permit-realistic
33 × 31 = 1 023 allele draws easily crossing the A_target = 704 species
threshold). Both locations sit in the same BL, so their pilot outputs
are directly comparable — a *within-BL* contrast that avoids
between-BL confounds. The two-location pilot uses ≤ 1 100 seed
genotypes total (< 6 % of the full 2025 recipe of ~18 640).

**Field-team recipe for the pilot** — download this and hand it to the
lab: [`step29_field_team_sampling_recipe.tsv`](tables/Phase5/step29_field_team_sampling_recipe.tsv).
Filter to `locationCode ∈ {EO67, EO27-1}` to isolate the pilot rows.

**Companion tables for pilot review:**

- [`step28_seed_sampling_per_mother.tsv`](tables/Phase5/step28_seed_sampling_per_mother.tsv) — analyst view (per-mother targets, achievable coverage, budget flags).
- [`step29_sampling_per_location.tsv`](tables/Phase5/step29_sampling_per_location.tsv) — per-location design table (M, target, delivered coverage).
- [`step29_location_connectivity.tsv`](tables/Phase5/step29_location_connectivity.tsv) — the 25 m primary + 10 m / 50 m sensitivity connectivity metrics.
- [`step30_A_prediction_location_diversity.tsv`](tables/Phase5/step30_A_prediction_location_diversity.tsv) — predicted SRK allele count with 95 % credible interval per location.
- [`step30_A_prediction_location_pcompat.tsv`](tables/Phase5/step30_A_prediction_location_pcompat.tsv) — predicted random-mating pollen compatibility with 95 % credible interval per location (finite-population model at 25 m).

**Downstream in Phase B.** Once EO67 and EO27-1 seed genotypes exist,
run Step 30 with the pilot subset:

```
python step30_srk_diversity_prediction_vs_observed.py \
    --year 2025 \
    --seed-genotypes  real_seeds_BL4_pilot.tsv \
    --mother-genotypes real_mothers_BL4_pilot.tsv \
    --match-seed-count
```

The two-point mate-limitation regression and SI-escape test will fire
on the two pilot locations; scaling to the full 39-location dataset
is then just a matter of adding rows to the two TSV inputs.

---

## Output map — quick reference (grouped by phase)

All Phase A tables are direct download links — click the filename to open
or right-click → *Save link as…* to pull the TSV into your local pipeline.

### Phase A — Preliminary (before SRK genotyping)

**Tables** (all TSV):

- [`step28_events_spatial_neighborhood.tsv`](tables/Phase5/step28_events_spatial_neighborhood.tsv) — event spatial reference (lat, lon, per-radius neighbourhood counts).
- [`step28_seed_sampling_per_mother.tsv`](tables/Phase5/step28_seed_sampling_per_mother.tsv) — per-mother analyst view (K, achievable coverage, budget flags).
- [`step28_coverage_curves_by_Nfertile.tsv`](tables/Phase5/step28_coverage_curves_by_Nfertile.tsv) — per-mother analytical + simulation curves by event-size bin.
- [`step28_aggregation_curves_by_Nfertile.tsv`](tables/Phase5/step28_aggregation_curves_by_Nfertile.tsv) — aggregation curves showing expected distinct alleles at a location for M = 1 … 30 mothers × 29 seeds each, one row per (event-size bin, M).
- [`step29_sampling_per_event.tsv`](tables/Phase5/step29_sampling_per_event.tsv) — per-event allocation.
- [`step29_sampling_per_location.tsv`](tables/Phase5/step29_sampling_per_location.tsv) — per-location design table with connectivity-informed columns.
- [`step29_location_coverage_curves.tsv`](tables/Phase5/step29_location_coverage_curves.tsv) — analytical curves by location size.
- [`step29_field_team_sampling_recipe.tsv`](tables/Phase5/step29_field_team_sampling_recipe.tsv) — **FIELD TEAM per-germplasmID recipe** (one row per mother with the actionable seed count).
- [`step29_location_connectivity.tsv`](tables/Phase5/step29_location_connectivity.tsv) — within-location pollen connectivity at 10 / **25 (primary)** / 50 m.
- [`step30_A_prediction_prior_frequencies.tsv`](tables/Phase5/step30_A_prediction_prior_frequencies.tsv) — the P1 species-wide Fg prior.
- [`step30_A_prediction_location_diversity.tsv`](tables/Phase5/step30_A_prediction_location_diversity.tsv) — predicted SRK allele diversity per location.
- [`step30_A_prediction_location_pcompat.tsv`](tables/Phase5/step30_A_prediction_location_pcompat.tsv) — predicted random-mating pollen compatibility per location (finite-population model at 25 m).
- [`step30_A_prediction_per_mother_fecundation.tsv`](tables/Phase5/step30_A_prediction_per_mother_fecundation.tsv) — species-wide compatibility reference distribution.

**Figures** (PNG + PDF):

- [`step28_coverage_curves.pdf`](figures/Phase5/step28_coverage_curves.pdf) / [`.png`](figures/Phase5/step28_coverage_curves.png) — per-mother coverage vs seeds, one curve per event-size bin, with Rule 2 cap.
- [`step28_per_mother_budget.pdf`](figures/Phase5/step28_per_mother_budget.pdf) / [`.png`](figures/Phase5/step28_per_mother_budget.png) — per-mother seed budget vs recommended n.
- [`step29_location_coverage_curves.pdf`](figures/Phase5/step29_location_coverage_curves.pdf) / [`.png`](figures/Phase5/step29_location_coverage_curves.png) — uniform-K location curves.
- [`step29_location_coverage_curves_P1.pdf`](figures/Phase5/step29_location_coverage_curves_P1.pdf) / [`.png`](figures/Phase5/step29_location_coverage_curves_P1.png) — P1-prior location curves with per-location bars.
- [`step29_recommended_seeds_per_mother_P1.pdf`](figures/Phase5/step29_recommended_seeds_per_mother_P1.pdf) / [`.png`](figures/Phase5/step29_recommended_seeds_per_mother_P1.png) — per-location seed-genotyping recipe (green/amber/red tiers).
- [`step29_location_connectivity.pdf`](figures/Phase5/step29_location_connectivity.pdf) / [`.png`](figures/Phase5/step29_location_connectivity.png) — **primary connectivity map at 25 m**, BL-panelled.
- [`step29_location_connectivity_10m.pdf`](figures/Phase5/step29_location_connectivity_10m.pdf) / [`.png`](figures/Phase5/step29_location_connectivity_10m.png) — sensitivity: 10 m conservative.
- [`step29_location_connectivity_50m.pdf`](figures/Phase5/step29_location_connectivity_50m.pdf) / [`.png`](figures/Phase5/step29_location_connectivity_50m.png) — sensitivity: 50 m optimistic.
- [`step29_location_connectivity_radius_sensitivity.pdf`](figures/Phase5/step29_location_connectivity_radius_sensitivity.pdf) / [`.png`](figures/Phase5/step29_location_connectivity_radius_sensitivity.png) — aggregate across radii (justification for 25 m primary).
- [`step30_A_prediction_diversity.pdf`](figures/Phase5/step30_A_prediction_diversity.pdf) / [`.png`](figures/Phase5/step30_A_prediction_diversity.png) — BL-panelled SRK allele richness per location.
- [`step30_A_prediction_fecundation.pdf`](figures/Phase5/step30_A_prediction_fecundation.pdf) / [`.png`](figures/Phase5/step30_A_prediction_fecundation.png) — BL-panelled compatibility per location, traffic-light bands.
- [`step30_A_diversity_vs_pcompat.pdf`](figures/Phase5/step30_A_diversity_vs_pcompat.pdf) / [`.png`](figures/Phase5/step30_A_diversity_vs_pcompat.png) — cross-plot: SRK diversity × pollen compatibility.

### Phase B — Post-genotyping (requires observed seed genotypes)

**Tables** (produced by `step30_srk_diversity_prediction_vs_observed.py --seed-genotypes … --mother-genotypes …`):

- `tables/Phase5/step30_B_comparison_location_diversity.tsv` — posterior vs predicted diversity per location.
- `tables/Phase5/step30_B_comparison_maternal_vs_paternal.tsv` — spectra test per location.
- `tables/Phase5/step30_B_mate_limitation_per_location.tsv` — Test 1 · location-level.
- `tables/Phase5/step30_B_mate_limitation_per_mother.tsv` — Test 1 · per-mother detail.
- `tables/Phase5/step30_B_mate_limitation_coefficients.tsv` — Test 1 · β₁, β₂, β₃ estimates with 95 % CI.
- `tables/Phase5/step30_B_si_escape_permutation.tsv` — Test 2 · per-location Binomial test + FDR.

**Figures**:

- `figures/Phase5/step30_B_comparison_diversity.pdf/png` — observed vs predicted diversity.
- `figures/Phase5/step30_B_mate_limitation.pdf/png` — Test 1 · location scatter.
- `figures/Phase5/step30_B_si_escape.pdf/png` — Test 2 · per-location bars.

### Phase B — DEMO (synthetic pipeline-validation outputs, `--demo` mode)

Same filenames as Phase B above with `_B_` → `_B_DEMO_`. Figures also
carry a "DEMO — synthetic data" title suffix and a diagonal DEMO
watermark so a demo file can never be mistaken for a real result.
