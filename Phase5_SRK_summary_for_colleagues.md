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
- **Species-wide SRK prior (P1).** The **32 functional SRK allele
  groups (Fgs)** identified across LEPA, together with their
  species-wide empirical frequencies — a 32-slot probability
  distribution built from the Canu-amplicon L1 carrier inventory
  (49 sequence alleles collapsed into 32 Fgs across 263 individuals).
  "Draw an allele from P1" means pick one of the 32 identities with
  probability equal to its share of the inventory. Dominant Fg =
  FG001 (**41 %** of P1 mass). **Class I = 26 Fgs (~35 % of P1
  mass), Class II = 6 Fgs (~65 %)** — this inverts the "common =
  Class I" Brassica frequency pattern, which assumes neutrality. In
  LEPA genetic drift elevates FG001-006 regardless of class. See
  § A.7.2 of the long doc for the drift note; a phylogenetic
  reassignment is tracked as future work.

---

## Key concepts and terminology

Every per-location quantity in this doc is derived from a nested
spatial hierarchy. The two top levels are **quoted verbatim from the
LEPA DB `Terms` table** (the canonical glossary that ships with
`LEPA_SQL.db`) so the vocabulary matches every other LEPA analysis.
Reading from the biggest unit down to the individual plant:

- **Location** (`locationID` / `locationCode`) — DB `Locations`
  table: `locationID` = "Location Unique Barcode #" (Darwin Core
  `dwc:locationID`); `locationCode` = the EO code of the sampling
  site ("Report the unique EO # where the sampling is conducted
  (e.g., EO38)"). May span several slick spots.
  - **Within-EO location split (Phase 5 refinement).** When events
    inside the same EO sit **≥ 500 m apart with no bridging events
    in between**, Phase 5 tracks the disjoint pieces as separate
    `locationCode`s (e.g. EO24 → `EO24`, `EO24-1`, `EO24-2`,
    `EO24-7`; EO27 → `EO27`, `EO27-1`, `EO27-3`, `EO27RT`), each
    with its own `locationID`. In the current data **5 EOs are
    split this way** (EO18, EO24, EO25, EO26, EO27 → 16 Phase 5
    locationCodes). 500 m is 10× LEPA's primary pollinator radius
    (50 m), so these sub-locations cannot share pollen under any
    plausible flight distance and must be modelled as independent
    drift units.
  - **Event** (`eventID` / `occurrenceID`) — DB `Events` table:
    "**an 'Event' refers to an occupied slick spot within a
    Location**" (Darwin Core `dwc:eventID`). Each event has its own
    census of fertile plants and its own coordinates.
    - **50 m connected component** (`component_id_50m`) — a group
      of events whose plants sit within 50 m pollinator-flight range
      of each other. **Plants in the same component share a pollen
      pool; plants in different components — even at the same
      location — do not.** Phase 5 derived concept, not a DB term.
      - **Mother plant** (`germplasmID`) — an individual plant
        already collected and stored in the LEPA DB, sitting at one
        specific event.

Everything else in the pipeline is a count or a derived number
sitting on top of this hierarchy:

| Concept (full English name) | Code identifier | Rooted at | Definition |
|---|---|---|---|
| Fertile plant census | `total_n_fertile` | location | All fertile plants at a location, summed across every event. The biological potential, with no spatial filtering. |
| **Effective mating pool size** (also called **N_fertile_effective**) | `N_fert_eff` | **component (primary); location (diagnostic)** | **Per component:** `component_N_fertile`, the fertile plants that share a single 50 m pollen pool — this is the drift unit for both diversity and pollen-compatibility prediction. **Per location:** `total_n_fertile × largest_component_share_50m`, a one-number *fragmentation diagnostic* used in the N_fertile_effective figure. The two agree when the whole location is one component. |
| Connectivity share | `largest_component_share_50m` | location, built from components | `N_fert_eff (location) ÷ total_n_fertile`. 1.0 = fully connected (no fragmentation); 0.3 = 70 % of the raw census is drift-irrelevant. Fragmentation diagnostic only — the Phase A prediction loops over every component. |
| Fragmentation-aware mother target | `M_frag` | event, derived from components | For each event, the number of mothers to sample so that each 50 m connected component reaches 90 % allele-detection coverage, with a ≥ 1-per-event maternal-genotype floor. Sums across events to the location-level `M_frag_aware`. |
| Rule 2 tetraploid seed cap | 15 seeds/mother | mother plant | Each seed contributes 2 paternal allele draws from the local pollen pool. 15 seeds/mother is the coupon-collector floor for a mother to see every allele in her component's pollen pool with 90 % probability. |
| Species prior | `P1` | species-wide | The **32 Fgs identified across LEPA plus their empirical species-wide frequencies** — a 32-slot probability vector that sums to 1, built from the Canu-amplicon L1 carrier inventory. Common Fgs (e.g. FG001 at 41 %) have a large slot; rare ones have a small slot. Every per-location prediction draws alleles from P1, so small locations lose the rare Fgs to drift by chance. |

**Why components matter in one sentence.** Every per-location
prediction in this doc — SRK diversity, pollen compatibility, mother
allocation — is built **component-by-component**, because a 50 m
connected component is what trades pollen. Each component gets its
own `component_N_fertile`; location-level numbers are the **set
union** across components for diversity (Fgs are a set) and the
**size-weighted mean** across components for pollen compatibility
(a continuous rate). A location with 500 fertile plants spread
across 20 isolated slickspots behaves like 20 small drift-prone
pools, not one pool of 500.

### How the five quantities chain together

1. **Location → events.** Raw census `total_n_fertile` = Σ event
   `n_fertile` across the location's 50 m-resolved events
   (step28_events_spatial_neighborhood.tsv).
2. **Events → 50 m components.** Within each location, connect
   events whose fertile plants sit within ≤ 50 m; the connected
   components are the `component_id_50m` units
   (step29c_event_to_component_50m.tsv).
3. **Components → effective mating pool.**
   `component_N_fertile_c` = Σ_{events ∈ c} `n_fertile_e` — the
   plants that actually share one 50 m pollen pool and the drift
   unit for both predictions.
4. **Effective mating pool → SRK diversity.** For each component,
   draw `4 × component_N_fertile_c` alleles from P1 → the component's
   present Fgs. Location pool size = |union of component Fg sets|.
   Sampling is simulated per component too (mothers and seeds
   distributed proportional to component size).
5. **Effective mating pool → pollen compatibility.** For each
   component, sample mothers from its local frequencies under the
   empirical LEPA zygosity and compute sporophytic Class I / II
   P_compat per mother. Location P_compat = size-weighted mean of
   per-component P_compat.

The whole causal chain **fragmentation → drift → mate limitation**
enters at step 3 and comes out at steps 4 and 5: fragmented
locations have small components → small drift-prone pools → fewer
Fgs and lower pollen compatibility than a one-pool location of the
same raw census would predict.

---

## The question this doc builds toward

**How many seeds per mother should the lab genotype for Part C
testing?** That single number is the practical output of the whole
framework, and it depends on three things — none of which can be
skipped:

- **Fragmentation** — how many adults actually share a pollen
  environment at each location.
- **SRK diversity** — how many Fg alleles Nature holds locally under
  fragmentation-driven drift.
- **Pollen compatibility (P_compat)** — how much of the local pool
  a mother is compatible with under sporophytic Class I / II SI.

The next sections build those three quantities up in order. The
sampling design — mother count per location, seed count per mother
— is derived at the end, once the causal chain is on the table.

---

## Link 1 — Fragmentation of pollen flow (Step 29b)

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

![Figure 2 — Pollinator-radius sensitivity sweep. Four panels showing how connectivity, fragmentation-aware sampling cost, predicted P_compat, and the sustainable-band fraction of locations change across radii from 10 to 200 m. Connectivity plateaus at ≥ 75 m and P_compat is radius-independent under empirical zygosity, justifying 50 m as the primary radius.](figures/Phase5/step30_A_radius_sensitivity.png)

---

## The pivotal metric — `N_fertile_effective`

**Question.** How many adults at a location actually share a pollen
environment, and therefore contribute to genetic drift on the local
SRK pool?

**Approach.** Combine the raw census with the 50 m connectivity
share from Step 29b:

```
N_fert_eff = total_n_fertile × largest_component_share_50m
```

**Result.**

- **Well-connected locations (connectivity share ≈ 1).** EO29,
  EO70, EO118, EO26-3 (small), EO67, EO24 group — every fertile
  plant sits in the same 50 m mating pool. Raw census = `N_fert_eff`.
- **Partially fragmented (share 0.5 – 0.9).** EO76 (517 → 401),
  EO61 (543 → 315), EO32 (466 → 324), EO18-7 (242 → 158) — a large
  census loses 10 – 40 % of adults to fragmentation.
- **Heavily fragmented (share < 0.5).** EO27-1 (395 → 147, 37 %),
  EO27-1 (371 → 116, 31 %), EO18-8 (123 → 52, 42 %), EO26-2
  (22 → 11, 50 %) — half or more of the raw census is drift-
  irrelevant.

**Why it matters.** `N_fert_eff` is the single connectivity-adjusted
number that every Phase A per-location prediction is built on:

- **Figure 3 (§ 30.1)** — Panel A (unbiased local Fg diversity) is a
  coupon-collector draw of `4 × N_fert_eff` alleles from P1. Small
  `N_fert_eff` → small local pool → drift-collapsed diversity.
- **Figure 4 (§ 30.2)** — random-mating pollen compatibility is
  simulated on a local Fg pool sized by `4 × N_fert_eff`. Small
  `N_fert_eff` → skewed local frequencies → mothers overlap more
  with candidate fathers.
- **Figure 5 (§ 30.3)** — event-scale and location-scale
  fragmentation indices decompose into the same shrinkage factor
  used to compute `N_fert_eff`.
- **§ B.4.2 fragmentation-aware sampling** — the reason
  22 / 39 locations need MORE mothers under the honest allocation
  is the same 50 m connectivity shrinkage exposed here.

Wherever the framework refers to `N_fertile` as a biological input,
it means `N_fert_eff`. The raw census is Nature's biological
potential; the effective count is what actually matters for
reproduction under 50 m pollinator flight.

![Figure — `N_fertile_effective` per LEPA location, panelled by Bottleneck Lineage in BL_ORDER. Y-axis labels give `locationCode (raw N, N_fert_eff)`. **Panel A** — raw census, log scale. Nature's biological potential. **Panel B** — `N_fert_eff = raw N × largest_component_share_50m`, log scale. Adults that actually share a 50 m mating pool. **Panel C** — connectivity share = Panel B ÷ Panel A. Dashed line at 1.0 = no fragmentation. This is the single input that drives Figures 3, 4, and 5.](figures/Phase5/step30_A_N_fertile_effective.png)

---

## Step 30 Phase A — Per-location predictions

### 30.1 Predicted SRK diversity per location — unbiased truth vs sampling

**Two questions, kept strictly separate.** SRK diversity at a
location has two very different meanings and the pipeline predicts
both. This section uses the **existing LEPA dataset** — actual
mothers × actual seeds/mother in the DB. The prospective design
question ("how many mothers × how many seeds do we NEED for a new
location?") is answered in Steps 28 – 29.

1. **What Nature actually holds at this location** (unbiased truth) —
   how many SRK alleles are physically present. Depends only on
   `N_fertile_effective` = fertile plants × 50 m connectivity. This
   is the Link 2 output of the causal chain and the raw material on
   which the P_compat prediction (§ 30.2) operates.
2. **What our sampling will recover** (sampling-inferred) — expected
   detection given `A_delivered = 4·M_mothers + 2·total_seeds` at each
   location. **Every seed's 2 paternal alleles are direct samples of
   the local pollen donor pool** — the seed genotyping is a
   pollen-pool characterisation experiment.

**Approach — how the simulation works, step by step.**

Think of the species prior P1 as a bag with 32 kinds of coloured
balls — one colour per Fg. The bag is **not evenly filled**: FG001
occupies 41 % of it, the next five Fgs together another ~25 %, and
the 20+ rare Fgs each fill under 2 %. Those shares are the
*species-wide empirical frequencies* measured from the Canu-amplicon
263-individual inventory. Each LEPA location only gets a limited
handful of draws from that bag, so small locations are almost
guaranteed to miss the rare colours by chance — that is **genetic
drift**. The prediction turns this intuition into numbers in five
steps, run for every 50 m component inside every location:

1. **Build the component's local pool (drift step).** For a
   component with `component_N_fertile_c` fertile plants, draw
   `4 × component_N_fertile_c` balls from P1 (because every plant is
   tetraploid and carries 4 SRK alleles). The unique colours drawn
   are the Fgs **present at that component**; the proportions in the
   draw are the component's local Fg frequency vector `f_local_c`.
   Small `N_c` → few draws → missing some rarer Fgs.

2. **Sum components to the location (union step).** The location
   holds a Fg if *any* of its components holds it, so the location's
   pool size on this replicate = |union of per-component present
   sets|. That is why diversity is unioned (set operation), not
   averaged (continuous-rate operation).

3. **Simulate the sampling (what we will actually see).** Distribute
   the location's `M_mothers` and `total_seeds` across components
   proportional to component size (largest-remainder: the component
   with the most plants gets the most sampling effort). Per
   component this gives `M_c` mothers and `seeds_c` seeds, so
   `A_delivered_c = 4·M_c + 2·seeds_c` allele draws (4 maternal
   alleles per mother + 2 paternal alleles per seed). Draw
   `A_delivered_c` balls from the component's own `f_local_c` and
   record which colours were hit — these are the Fgs **detected at
   that component**.

4. **Union detected across components, divide by union present.**
   Location detected = |union of per-component detected sets|.
   Location coverage = detected / pool size, on this replicate.

5. **Repeat 1 000 times, report mean and 95 % credible interval.**
   Each replicate uses an independent random draw from P1 at step 1
   and an independent sampling draw at step 3; the 2.5 % and 97.5 %
   percentiles across replicates give the credible interval seen as
   error bars in Figure 3.

Why the per-component view matters: a drift-collapsed small
component (say 2 plants, holding 4 Fgs) left with 0 mothers under
step 3 contributes its Fgs to the location's "present" count at
step 2 but nothing to the "detected" count at step 4 — that is the
sampling blind spot that pulls EO26-3 and EO27RT down to 97 %
coverage below.

**Result.**

| Location | BL | N_fert_eff (50 m) | n components | M_mothers | total_seeds | Nature holds | Sampling detects | Coverage |
|---|---|---|---|---|---|---|---|---|
| EO24-2 (1 plant, BL5 tail) | BL5 | 1 | 1 | 1 | 1 | ~3 | ~2.6 | **87 %** |
| EO24 (BL5 tail) | BL5 | 2 | 1 | 1 | 5 | ~4.5 | ~4.1 | **91 %** |
| EO26-3 (fragmented BL1) | BL1 | 5 | 2 | 4 | ~70 | ~7.6 | ~7.4 | **97 %** |
| EO27RT (fragmented BL4) | BL4 | 177 | 6 | 22 | ~1700 | ~29 | ~28.1 | **97 %** |
| EO67 (small BL4 pilot) | BL4 | 10 | 2 | 4 | 71 | ~11 | ~10.7 | 99 % |
| EO27-1 (large BL4 pilot) | BL4 | 116 | 2 | 33 | 2465 | ~27 | ~27 | ~100 % |
| EO76 (largest BL3) | BL3 | 517 | 6 | 62 | 11723 | ~32 | ~32 | ~100 % |

**Every location clears the 90 % target; 6/39 now sit below 99 %
coverage** (vs 2 under the previous single-pool model). The newly
visible gap is the **sampling blind spot** in fragmented locations:
when drift-collapsed small components (1–5 plants) receive 0 mothers
under proportional allocation but still carry Fgs the larger
components do not — EO26-3 and EO27RT are the clearest examples.

The existing LEPA dataset therefore characterises Nature's truth at
every location. Total seed counts in the DB span 1 (EO24-2) →
11 723 (EO76), so most locations are well provisioned.

![Figure 3 — Predicted SRK allele diversity per LEPA location, **built per 50 m component and unioned at the location level**. Panelled by Bottleneck Lineage in canonical BL_ORDER. Y-axis labels give `locationCode (N_fert_eff, M_mothers, total_seeds)` where `N_fert_eff` is the sum of component sizes across the location. **Panel A** — Fgs present somewhere in the location (union across its components, drift-only; feeds the P_compat prediction). **Panel B** — Fgs detected by the actual LEPA DB sampling, with mothers and seeds allocated across components proportional to component size; each component's sampling is simulated against its own frequency vector and the location-level detected count is the union across components. **Panel C** — Coverage = Panel B / Panel A, with the 90 % target line. Six locations now sit below 99 % coverage (vs 2 under the previous single-pool model) — the small components left with 0 mothers under proportional allocation are now visible as a coverage gap.](figures/Phase5/step30_A_diversity_unbiased_vs_sampling.png)

### 30.2 Predicted pollen compatibility per location

**Question.** Under the sporophytic Class I / II + empirical-zygosity
model, what fraction of pollen would a mother at each location be
compatible with under random mating?

**Approach — how the simulation works, step by step.**

Same per-component logic as § 30.1, but the recognition rule changes
from "did we observe this Fg at least once?" to "if this mother's
expressed Fgs are {a, b}, what fraction of pollen in her component
would she accept?" Five steps per replicate:

1. **Build the component's local pool (same as § 30.1 step 1).**
   Draw `4 × component_N_fertile_c` alleles from P1 → local Fg
   frequencies `f_local_c`. This is the drift-collapsed pollen pool
   the component's mothers see.

2. **Sample mother genotypes under the empirical LEPA zygosity.**
   Each mother gets 4 SRK alleles drawn from `f_local_c`, but
   combined into tetraploid genotypes that match the observed LEPA
   mix: 66 % of mothers end up with 1 distinct Fg (homozygotes),
   32 % with 2 distinct, 2 % with 3 distinct (from
   `srk_zygosity_empirical.tsv` — see Part C § C.0 for how this was
   measured from the 367 Canu-amplicon adults).

3. **Compute per-mother pollen compatibility under sporophytic Class
   I / II SI.** A pollen parent expresses only Class I alleles if it
   carries any, else all its Class II alleles co-dominantly. A cross
   is rejected if parents share any expressed allele. For a mother
   expressing Fg set *M*:
   - Mother has ≥ 1 Class I allele: `P_compat = (1 − p(M))⁴` where
     p(M) is the local frequency sum of her expressed Fgs.
   - Mother is pure Class II: `P_compat = 1 − (1 − p_I)⁴ +
     (1 − p_I − p(M))⁴`.

   Average across the mothers sampled at this component → component
   P_compat on this replicate.

4. **Aggregate components to the location (weighted mean, not
   union).** P_compat is a continuous rate, so a location's number
   is the component sizes weighted mean:
   P_compat_location = Σ_c (P_compat_c · component_N_fertile_c) /
   Σ_c component_N_fertile_c. Big components dominate; small
   components still count — that is where drift pulls the location
   mean down.

5. **Repeat 400 times, report mean and 95 % credible interval.**

Why no sampling term enters P_compat: § 30.2 asks "what fraction of
pollen is compatible with a random mother at this component?" — a
biological property of the pollen pool, not a function of how many
seeds we sequence. Seed counts enter § 30.1 (sampling detection) and
Phase B (the regression against observed seed set), not here.

Both tables are emitted:
[`step30_A_prediction_location_pcompat.tsv`](tables/Phase5/step30_A_prediction_location_pcompat.tsv)
(headline per location) +
[`step30_A_prediction_component_pcompat.tsv`](tables/Phase5/step30_A_prediction_component_pcompat.tsv)
(one row per component, so struggling sub-components inside a
sustainable-mean location remain visible).

**Result.**

- **Species-mean P_compat = 0.78** (sustainable band by construction).
- **Traffic-light bands** (1/3 and 2/3 of species mean): failed
  < 0.23, struggling 0.26–0.52, sustainable ≥ 0.52.
- **BL5 tail (EO24 group)** is the primary conservation concern —
  predicted P_compat 0.42–0.51 with wide credible intervals
  entering the struggling band.
- All other BLs sit in the sustainable band at the mean; the tightest
  credible intervals belong to the largest locations (EO76, EO61).

**Connection to Figure 3.** Uses the same drift mechanism as the
diversity figure, except drift is now modelled **per component**
(one pool per 50 m connected component sized by
`4 × component_N_fertile`) rather than through a single
"largest-component" proxy. Figure 3 shows that the existing LEPA
seed data recovers the true local Fg pool at every location — that
validation transfers directly to Figure 4, meaning the location-
mean P_compat from the per-component simulation is a faithful
estimator of the population-mean P_compat.

![Figure 4 — Predicted per-location pollen compatibility under sporophytic Class I / II + empirical LEPA zygosity. Y-axis labels give `locationCode (N_fert_eff, M_mothers)` — the same convention as Figure 3, without total_seeds because seed counts do not enter this prediction. One dot per location, error bars = 95 % credible interval, panelled by Bottleneck Lineage. Traffic-light bands: red = failed (< 0.26), amber = struggling (0.26–0.52), green = sustainable (≥ 0.52). Species mean = 0.78. BL5 tail slips into the struggling band; all other locations sit in the sustainable band at the mean.](figures/Phase5/step30_A_prediction_fecundation.png)

### 30.3 Per-location mating neighbourhood size

**Question.** Independent of drift on allele frequencies, how big is
each location's event-scale mating neighbourhood — and how uniform
is it across the location's events?

**Approach — purely spatial, no allele frequencies.** For every
event, count `N_reachable_50m = Σ N_fertile in other events within
50 m`: the number of pollen-donor plants a flower on that event
could reach at the 50 m pollinator range. Then **one box plot per
location** summarises its events' `N_reachable_50m` distribution,
with locations stacked vertically by Bottleneck Lineage.

**How to read it.**
- **Where a box sits on the x-axis** → typical event-scale neighbourhood size.
- **How wide the box is / how long the whiskers are** → how uniform the location's events are. A narrow box = events have similar neighbourhoods; a wide box = some events are well-connected while others sit alone.
- **The red dotted line (N = 1)** = the single-plant SI floor: a flower on this event has nobody to mate with.
- **The grey dashed line (N = 8 plants = 32 tetraploid allele copies)** = the coupon-collector floor for the 32-Fg species pool. Events to the left of it have fewer reachable allele copies than the species-wide ceiling and are therefore drift-limited at the event scale, no matter the location mean.

**Result.** BL5 tail and the smallest BL1 locations (EO26-3 and some
EO8 sub-locations) sit near or at the N = 1 floor — their events
cannot reach the coupon-collector floor on their own. BL4, BL3's
EO76, BL1's EO61 / EO29, and BL2's EO70 sit at or above the 32-copy
line across all their events — well-connected everywhere. The
widest boxes (EO18-7, EO8) are the layered locations: some events
are richly connected, others are isolated.

**Prefer the N_fertile_effective figure (Panel C) for a single
per-location fragmentation score.** This figure is the complementary
*event-level* texture.

![Figure 5 — Per-location box plots of `N_reachable_50m` (pollen-donor plants reachable within 50 m per event). One row per location, grouped by Bottleneck Lineage; row labels give `(events, adults)`. Each box summarises the location's events' reachable-neighbour counts (log x-axis). **Red dotted line: N = 1** = the single-plant SI floor (nobody to mate with on that event). **Grey dashed line: N = 8 plants = 32 tetraploid allele copies** = the coupon-collector floor for the 32-Fg species pool. Narrow boxes = uniform event neighbourhoods; wide boxes = layered location (some events well-connected, others isolated). Purely spatial — no allele frequencies enter. Source: `step30b_fragmentation_index.py`.](figures/Phase5/step30_A_fragmentation_index.png)

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
- **EO70 is a striking outlier**: observed 0.53 vs predicted 0.70 —
  a 0.17-unit shortfall, right at the struggling/sustainable boundary
  (0.52 threshold). This is a genuine biological signal — local Fg
  pool skew makes mothers overlap far more with neighbouring fathers
  than P1 would predict. Not a model failure.
- **EO18 and EO67** show milder deviations in the same direction
  (observed ~0.69–0.70 vs predicted ~0.75–0.76). Both are candidates
  for future conservation intervention.
- **Observed zygosity distributions per EO** track the species-wide
  66 / 32 / 2 % well, except EO25 (49 % 2-distinct vs 32 %
  species-wide) — consistent with its perfect-match position on
  the diagonal.

![Figure 6 — EO-level empirical validation of the sporophytic P_compat model. Panel A: observed vs predicted mean pollen compatibility per EO (n ≥ 10 individuals), with 1:1 diagonal and traffic-light bands (failed < 0.26, struggling < 0.52, sustainable ≥ 0.52). 4/6 EOs sit on the diagonal within their 95 % CI; EO70 is the clearest outlier (observed 0.53 vs predicted 0.70, right at the struggling/sustainable boundary) — genuine local Fg pool skew, not model failure. Panel B: distinct-identity distribution per EO vs the species-wide reference (66 / 32 / 2 %).](figures/Phase5/step30_C_pcompat_observed_vs_predicted.png)

**Caveat.** P1 was built from these same individuals, so this
comparison does not test absolute calibration but robustness to
per-EO drift. EO-scale is intermediate; the finer per-location test
will happen in Phase B once seed genotypes are available.

### C.0.a — Part C anchor at Phase 5 location scale (clean-overlap EOs)

**Question.** Three of the six § C.0 EOs are **1:1 with a Phase 5
locationCode** — no 500 m within-EO split — so their adult SRK
genotypes can be used **directly** as a Part C validation anchor at
the Phase 5 location scale, no event-level remap required. Does the
Phase 5 **per-location** prediction reproduce what we see at the
location scale?

**The clean overlap set.** EO67 (37 adults), EO70 (74), EO76 (76) =
**187 adults**. The other three § C.0 EOs (EO18, EO25, EO27) are
split under the 500 m rule and need an `Individual → germplasmID →
eventID → locationCode` join before their 151 adults can be used
here; that is scheduled as future work.

**Approach.** For each clean-overlap location, run `step30d_partC_clean_overlap.py`:
1. Observed per-mother P_compat under sporophytic Class I / II +
   empirical zygosity, with fathers drawn from **observed local Fg
   frequencies** (not P1).
2. Observed distinct Fgs = size of the Fg set recovered from the
   adult genotypes.
3. **No-drift upper bound on distinct Fgs** at the adult sample size:
   `E[distinct Fgs | 4·n_adults draws from P1]`. If this upper bound
   is near the species ceiling (which it is for all three locations:
   19.2 at n=37, 24.0 at n=74, 24.0 at n=76), any large gap down to
   the observed count isolates **drift**, not sampling.
4. Compare observed to the Phase 5 per-location prediction
   (`step30_A_prediction_location_pcompat.tsv` for P_compat,
   `step30_A_prediction_location_diversity.tsv` for the component-
   unioned pool size).

**Result.**

| Location | Adults | Observed P_compat (95 % CI) | Phase 5 pred P_compat (95 % CI) | Observed Fgs | No-drift upper bound | Phase 5 pred Fgs (95 % CI) |
|---|---:|---|---|---:|---:|---|
| **EO67** | 37 | 0.697 (0.660–0.737) | 0.724 (0.632–0.803) | **7** of 32 | 19.2 | 10.8 (7.0–15.0) |
| **EO70** | 74 | **0.532** (0.510–0.559) | **0.776** (0.740–0.806) | **6** of 32 | 23.8 | 28.3 (25.0–31.0) |
| **EO76** | 76 | 0.722 (0.694–0.753) | 0.775 (0.749–0.797) | **9** of 32 | 24.0 | 31.8 (31.0–32.0) |

- **EO67** — observed and Phase 5 predicted CIs overlap on P_compat;
  observed Fg count sits inside the Phase 5 CI → **model passes** at
  this small BL4 location.
- **EO70** — observed P_compat 0.24 units BELOW Phase 5 prediction
  (CIs do not overlap); observed Fg diversity 6/32 vs predicted
  28/32 vs no-drift upper bound 24. **Massive drift collapse at a
  large location** — a genuine biological signal, confirming the
  § C.0 EO-scale finding at the Phase 5 location scale.
- **EO76** — observed P_compat marginally below prediction; observed
  Fg diversity 9/32 vs predicted 32/32 vs no-drift upper bound 24.
  **Severe diversity collapse at the single largest LEPA site** —
  the P_compat model still holds up because Class II (6 Fgs,
  including FG001 at 41 %) dominates and makes most mothers
  compatible even with limited diversity.

**Outputs.**
- [`step30_B_partC_clean_overlap_per_location.tsv`](tables/Phase5/step30_B_partC_clean_overlap_per_location.tsv) — one row per clean location: observed vs predicted P_compat and SRK diversity.
- [`step30_B_partC_clean_overlap_fg_frequencies.tsv`](tables/Phase5/step30_B_partC_clean_overlap_fg_frequencies.tsv) — long form, one row per (location, Fg) with observed f and P1 f.
- [`step30_B_partC_clean_overlap.png`](figures/Phase5/step30_B_partC_clean_overlap.png) — three-panel figure (P_compat scatter, Fg count comparison, Fg frequency spectrum).

![Figure 7 — Phase 5 Part C anchor at the three clean-overlap EOs (EO67, EO70, EO76; 1:1 with a Phase 5 locationCode). **Panel A** — observed vs Phase 5 predicted pollen compatibility per location, with traffic-light bands and the 1:1 diagonal. EO67 and EO76 near the diagonal; EO70 is the drift outlier (observed 0.53 vs predicted 0.78). **Panel B** — distinct Fg count per location: **filled square = observed**, **open blue circle = no-drift upper bound at the adult sample size** (`E[distinct Fgs | 4·n_adults draws from P1]`). The large gap between open circle and filled square at EO70 and EO76 isolates drift from sampling. **Panel C** — observed Fg frequency spectrum per location vs the species-wide P1 reference (grey). FG001 dominates everywhere, consistent with drift elevating the common Fg at small-population tails. Source: `step30d_partC_clean_overlap.py`.](figures/Phase5/step30_B_partC_clean_overlap.png)

---

## Sampling design — the answer to "how many mothers × how many seeds?"

Now that fragmentation, drift and mate limitation are on the table,
the sampling design falls out of the causal chain. Every choice below
is derived from a quantity that has already been established, not
proposed independently.

### How many seeds per mother?

**Question.** How many seeds must we genotype from each mother to
characterise her local pollen environment well enough to test the
mate-limitation prediction in § C.1?

**Approach — Rule 2 tetraploid allele-detection floor.** Each seed
contributes 2 paternal SRK alleles under tetraploid sporophytic
inheritance. Under a uniform-pollen coupon collector, the probability
of missing any one of K local paternal Fgs after 2·n seed alleles is
`(1 − 1/K)^(2n)`. Setting `(1 − 1/K)^(2n) = 0.10` gives the seeds
needed for a 90 % chance to see every Fg. **Rule 2 caps this at
n = 15 seeds/mother** — the point where the return per additional
seed is negligible.

![Figure — Step 28 coverage curves. Panel A: per-mother detection of the local Fg pool as a function of seeds genotyped, one curve per event-size bin; the vertical red line at 15 marks the tetraploid Rule 2 cap. Panel B: aggregation of coverage across mothers at a location (15 seeds × M). Every event size reaches its local ceiling by 5 mothers.](figures/Phase5/step28_coverage_curves.png)

**Does 15 seeds also give enough Part C testing power?**
Coupon-collector justifies 15 as the allele-detection floor. Two
additional simulations confirm it also passes the regression
thresholds:

- **Per-mother P_compat precision** — observed P_compat has binomial
  SE `√(p·(1−p) / n_seeds)`. Drops from ~0.28 at 3 seeds to
  ~0.12 at 15 (visible plateau); marginal gains beyond.
- **§ C.1 mate-limitation regression power** — full-pipeline
  simulation of 505 mothers × 39 locations. At n_seeds = 15:
  **99.6 % power** for a medium effect (β₁ = 100 seeds per unit
  P_compat), **69 % for a small effect (β₁ = 50, below the 80 %
  target)**. Errors-in-variables attenuation of β̂₁ is ~0.42 — real
  but doesn't prevent detection at realistic effect sizes.

15 seeds is comfortably enough for medium-to-large mate-limitation
signals. The marginal power at small effect sizes comes from the same
25-mother / 6-location shortage the mother-count section flags —
the 2026 top-up addresses both.

![Figure — Justifying 15 seeds for Part C testing. Left: per-mother P_compat precision vs seed count for four true P_compat values, with a plateau visible at ~15 seeds. Right: § C.1 mate-limitation regression power vs seed count for four effect sizes; 15 seeds delivers ≥ 99 % power for β₁ ≥ 100.](figures/Phase5/step28d_matelim_power.png)

### How many mothers per location?

**Question.** How many mothers must we sample at each location to
observe every SRK allele physically present in the location's mating
pool?

**Approach.** Under tetraploid sampling each adult contributes
`PLOIDY × N_fert_eff = 4·N_fert_eff` allele copies to the location's
pool. Coupon-collector target: 90 % chance of observing every Fg in
the local pool at the total delivered allele draws
`A_delivered = M × (4 + 2 × 15) = 34·M`. Plus a private-allele floor:
**at least one mother per event** (an isolated slickspot's private
Fg cannot be recovered from any other event).

**Refinement — fragmentation-aware allocation.** If a location is
fragmented into several disconnected 50 m mating pools, we decompose
the location into its components and allocate mothers per component,
plus the ≥ 1-per-event floor. The result is `M_frag` per event.

**Result.**

- **Design target: 505 mothers across 39 locations, 234 events.**
- **Already collected in the LEPA DB: 765 mothers** — 1.5× the target.
- **33 / 39 locations are fully covered** by the current DB.
- **6 / 39 locations are short of the target** — 25 mothers total.
  The six short locations (EO27-1, EO27RT, EO26-2, EO8, EO67, EO8)
  are the field-team top-up target for the 2026 season.

`M_frag` per event is the authoritative field-team recipe (see
[`step29c_sampling_frag_aware_per_event.tsv`](tables/Phase5/step29c_sampling_frag_aware_per_event.tsv)).

### The locked recipe

With mother count and seed count fixed, the design is frozen for the
next field season:

| Quantity | Value |
|---|---|
| Primary pollinator radius | **50 m** |
| Effective mating pool metric | **N_fert_eff** = census × 50 m largest-component share |
| Mother allocation per location | Fragmentation-aware `M_frag` (per 50 m component + ≥ 1 per event) |
| Seeds per mother (tetraploid Rule 2) | **15 seeds** |
| **Total mothers across 39 locations** | **505** |
| **Total seed genotypes** | **505 × 15 = 7 575 seeds** |
| Coupon-collector target | 90 % detection per 50 m component |

**Three authoritative files:**

- **Field-team recipe** (for a new field season):
  [`step29c_sampling_frag_aware_per_event.tsv`](tables/Phase5/step29c_sampling_frag_aware_per_event.tsv) —
  one row per event, `M_frag` = number of mothers to sample there.
- **Lab recipe for Part C** (draws specific mothers from the DB):
  [`step29c_partC_germplasmID_selection.tsv`](tables/Phase5/step29c_partC_germplasmID_selection.tsv) —
  one row per SELECTED `germplasmID` already in the LEPA DB, with
  `n_seeds_to_genotype = min(15, seeds_available)`. Selection is
  **per 50 m component**: mothers within a component share their
  pollen pool, so coverage travels freely within a component; only
  the ≥ 1-mother-per-event maternal-genotype floor is a strict
  per-event rule. Delivers **431 mothers × ≤ 15 seeds = 6 459 seeds**
  for Part C from the current DB, with **76 mothers short across 35
  components** flagged for a 2026 field top-up.
- **Event → component lookup**:
  [`step29c_event_to_component_50m.tsv`](tables/Phase5/step29c_event_to_component_50m.tsv) —
  one row per event mapping (locationID, eventID) → 50 m component,
  the single canonical source for "which events share a pollen pool?".

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
  P_compat mean 0.73, 95 % CI [0.56, 0.88]. Any Phase B observation
  dropping EO67 out of sustainable = decisive evidence of drift-driven
  mate limitation.
- **EO27-1** (aggregation regime) — P_compat 0.78, tight CI [0.72,
  0.84]. Solid anchor point in the sustainable band; if EO27-1 comes
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
| Predicted per-50 m-component diversity | [`tables/Phase5/step30_A_prediction_component_diversity.tsv`](tables/Phase5/step30_A_prediction_component_diversity.tsv) | Pool size, allocated mothers/seeds, coverage per component |
| Predicted per-location P_compat | [`tables/Phase5/step30_A_prediction_location_pcompat.tsv`](tables/Phase5/step30_A_prediction_location_pcompat.tsv) | Sporophytic + empirical zygosity |
| Predicted per-50 m-component P_compat | [`tables/Phase5/step30_A_prediction_component_pcompat.tsv`](tables/Phase5/step30_A_prediction_component_pcompat.tsv) | Exposes struggling sub-components inside a sustainable-mean location |
| Fragmentation indices | [`tables/Phase5/step30_A_fragmentation_per_event.tsv`](tables/Phase5/step30_A_fragmentation_per_event.tsv), [`_per_location.tsv`](tables/Phase5/step30_A_fragmentation_per_location.tsv) | Pure spatial, no allele frequencies |
| Event → 50 m component lookup | [`tables/Phase5/step29c_event_to_component_50m.tsv`](tables/Phase5/step29c_event_to_component_50m.tsv) | Canonical "which events share a pollen pool?" |
| Field-team recipe (new field season) | [`tables/Phase5/step29c_sampling_frag_aware_per_event.tsv`](tables/Phase5/step29c_sampling_frag_aware_per_event.tsv) | `M_frag` per event, authoritative |
| **Lab recipe for Part C — germplasmIDs to sample** | [`tables/Phase5/step29c_partC_germplasmID_selection.tsv`](tables/Phase5/step29c_partC_germplasmID_selection.tsv) | **One row per selected germplasmID already in the LEPA DB; per-50 m-component allocation + `n_seeds_to_genotype`. Drives Part C.** |
| Empirical validation (§ C.0) | [`tables/Phase5/step30_C_pcompat_validation_at_eo.tsv`](tables/Phase5/step30_C_pcompat_validation_at_eo.tsv) | 6 EOs, observed vs predicted |
| **Part C anchor at Phase 5 location scale (§ C.0.a)** | [`tables/Phase5/step30_B_partC_clean_overlap_per_location.tsv`](tables/Phase5/step30_B_partC_clean_overlap_per_location.tsv) + [`_fg_frequencies.tsv`](tables/Phase5/step30_B_partC_clean_overlap_fg_frequencies.tsv) | **Clean-overlap EOs (EO67, EO70, EO76) — 187 adults ready to feed Part C now; observed vs Phase 5 pred P_compat + Fg diversity + no-drift upper bound** |

**Key figures.**

- [`step30_A_radius_sensitivity.png`](figures/Phase5/step30_A_radius_sensitivity.png)
  — why 50 m is the right primary radius.
- [`step30_A_diversity_unbiased_vs_sampling.png`](figures/Phase5/step30_A_diversity_unbiased_vs_sampling.png)
  — three-panel: what Nature holds (unbiased) vs what our sampling detects vs coverage.
- [`step30_A_prediction_fecundation.png`](figures/Phase5/step30_A_prediction_fecundation.png)
  — predicted P_compat per location, traffic-light bands.
- [`step30_A_N_fertile_effective.png`](figures/Phase5/step30_A_N_fertile_effective.png)
  — raw census vs `N_fert_eff` vs connectivity share per location. The pivotal metric.
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
