# Phase 5 — Mate limitation and fragmentation in *Lepidium papilliferum*

**A compact overview for collaborators.** Framed on a single causal
chain (**fragmentation → genetic drift → mate limitation**); each
step below is presented as Question → Approach → Result. Full
technical detail (formulas, code, output filenames) lives in the long
companion doc [`Phase5_SRK_sampling_and_prediction.md`](Phase5_SRK_sampling_and_prediction.md).

---

## Executive summary

**Approach.** We aim to test within-location connectivity in
*Lepidium papilliferum* by sweeping candidate pollinator radii
(10 – 200 m) across the 39 LEPA locations with **2025 field data**
(wild, in-situ, with coordinates). The chosen **50 m primary radius**
(Figure 1) defines a **mating pool** — the set of plants reachable
in one pollinator flight — and quantifies within-location
**fragmentation**. These two spatial descriptors drive the
predictions of **SRK allele diversity** and **random-mating pollen
compatibility** — the two metrics that jointly determine whether
mate limitation operates at a location. Locations hold 1–6 mating
pools each; **31 of 101 pools sit below the 8-plant coupon-collector
floor** for the 32-allele species pool (Figure 2).

**Methodology.** For each mating pool we simulate drift from the
species-wide 32-allele empirical prior P1, aggregate SRK diversity
by set union across pools (Figure 3), and compute pollen
compatibility as a size-weighted mean under the **sporophytic
Class I / Class II + empirical zygosity** model (Figures 3b, 4).
The Part C anchor (Figures 5–6) compares predictions against
Phase 4 adult SRK genotypes at the three EOs — **EO67, EO70,
EO76** — that map 1:1 to a Phase 5 location.

**Result + interpretation.** The diversity prediction
**over-estimates** observed counts by ~20 alleles at every clean-
overlap location (predicted 11 / 28 / 32 vs observed 7 / 6 / 9),
because P1 itself has already absorbed decades of drift and local
pools have drifted further. **Yet the pollen-compatibility
prediction tracks observation closely** — EO67 (observed 0.70 /
predicted 0.72) and EO76 (0.72 / 0.78) overlap on their 95 % CIs;
only EO70 shows a real gap (0.53 vs 0.78), a drift signal consistent
with its FG024 monoculture. The hypothesis decomposition (Figure 7)
explains why: **within-class allele frequency spread, not allele
count, drives pollen compatibility**, and tetraploid zygosity
composition actively buffers locations with many homozygous mothers.
The decisive variable is not *how many* alleles survive but *how
their frequencies and genotypes are arranged* — which is why a
huge diversity gap can coexist with an accurate pollen-compatibility
prediction.

**Next test.** Within-BL5 pilot (option A2): **EO48_7** (connected,
1 mating pool) vs **EO18-7_19** (fragmented, 3 mating pools) —
20 mothers × 25 seeds → ~300 genotyped seedlings at 60 % germination.

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
else. It is where location size enters the model (via the per-
mating-pool `component_N_fertile`); a census of 1 000 plants split
across 20 disconnected slickspots behaves like 20 small drift-prone
pools of ~50 plants each, not one pool of 1 000.

| Link | What we measure | Steps that produce it |
|---|---|---|
| **1. Fragmentation** | 50 m mating-pool structure per location (count + size per pool) | Step 29b, Step 29c, Step 29d (§ A.5, § B.2 of long doc) |
| **2. Genetic drift on SRK** | Predicted local allele pool size + frequency composition **per mating pool**, aggregated to the location | Step 30 Phase A (§ A.6) |
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
| **Mating pool size** | `component_N_fertile` | component | The fertile plants that share a single 50 m pollen pool — **this is the drift unit for both diversity and pollen-compatibility prediction**. Each 50 m connected component = one mating pool. A location can hold one mating pool (fully connected) or several (fragmented); see Figure 2. |
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

## Link 1 — Pollinator-radius choice (Step 29b)

**Question.** What pollen-flight radius best represents LEPA's
mating process, given the biology of its small-bee pollinators?

**Approach.** Sweep candidate radii from 10 m to 200 m. At each
radius, build a within-location graph in which two adults are
connected if they sit within that radius of each other, and record
(i) the fraction of adults in a multi-event pollen-flow component
(connectivity), (ii) the fragmentation-aware sampling cost
(§ Sampling design), and (iii) the predicted random-mating pollen
compatibility under the sporophytic Class I/II model.

**Result.** Connectivity plateaus at ≥ 75 m; sampling cost stabilises
at 75–100 m; predicted pollen compatibility is radius-independent
under the empirical zygosity (§ 30.2). **50 m is the sweet spot**:
within the halictid / small-bee foraging literature range,
captures 43 % of the sampling-cost reduction, and keeps meaningful
fragmentation variation across BLs. Every downstream metric in this
doc — mating pool structure (Figure 2), event-scale reachability
(Figure 2b), SRK diversity prediction (Figure 3), pollen compatibility
prediction (Figure 4) — is computed at this 50 m choice.

![Figure 1 — Pollinator-radius sensitivity sweep. Four panels showing how connectivity, fragmentation-aware sampling cost, predicted P_compat, and the sustainable-band fraction of locations change across radii from 10 to 200 m. Connectivity plateaus at ≥ 75 m and P_compat is radius-independent under empirical zygosity, justifying 50 m as the primary radius.](figures/Phase5/step30_A_radius_sensitivity.png)

---

## Mating pool structure per location — the drift unit (Step 29c)

**The 50 m choice above is the knob that controls every downstream
fragmentation metric.** Fix the radius first; the mating pool
structure (Figure 2) and the event-scale reachability texture
(Figure 2b) follow directly from it, and in turn feed the Phase A
predictions (Figures 3, 4).

**Question.** Within a location, how many pollen pools does a plant
belong to, and how big is each one? A **50 m connected component**
is the Phase 5 **mating pool**: the set of adult plants whose events
are reachable from each other at the primary pollinator radius.
Plants in the same mating pool share pollen; plants in different
mating pools at the same location do not. The mating pool is the
**drift unit** on which every Phase A prediction is built.

**Approach.** For each location, build an event-level graph in which
two events are connected if any pair of their plants sits within
≤ 50 m. Connected components of that graph are the location's
mating pools. For each mating pool *c*, `component_N_fertile_c` is
the sum of `n_fertile_e` across its events — the adult count that
drives the per-component drift simulations in § 30.1 and § 30.2.
Reference: `step29c_fragmentation_aware_sampling.py` writes the
event → mating-pool lookup (`step29c_event_to_component_50m.tsv`);
`step29d_mating_pool_structure.py` builds the display here.

**Result (dataset-wide).**

- **39 locations hold 101 distinct 50 m mating pools**
  (mean 2.6 pools per location, median 2, max 6).
- Pool sizes span the full biological range: **1 adult (SI floor) → 420 adults**;
  median pool size = 23 adults.
- **31 of 101 mating pools (31 %) sit below the N = 8 coupon-collector
  floor** (fewer than 32 tetraploid allele copies in the pool) —
  these are drift-limited for the 32-allele species pool regardless
  of location mean.
- **6 of 101 mating pools (6 %) sit at the N = 1 single-plant SI floor** —
  a single plant has nobody to mate with at the 50 m radius.
- Pattern across BLs: BL5 tail (EO24 group) carries the small-pool
  burden; BL4 has large, mostly fragmentation-robust locations; BL1
  shows the widest *within-location* spread (EO8, EO26-3 split
  across many small pools).

**Why this matters.** Every Phase A prediction below is built on
this structure: the per-component SRK diversity prediction (Figure 3)
simulates each mating pool independently and unions the Fg sets at
the location level; the per-component pollen compatibility prediction
(Figure 4) runs the sporophytic simulation per pool and reports a
size-weighted mean per location; the fragmentation-aware sampling
allocation (§ Sampling design) assigns mothers per pool, with the
≥ 1-mother-per-event floor layered on top. Reading Figure 2 first
is the fastest way to anticipate which locations will stand out in
Figures 3 and 4.

![Figure 2 — Mating-pool structure per LEPA location. Two aligned panels; one row per locationID (unified `{EOID}_{locationID}` label), rows stacked vertically and grouped by Bottleneck Lineage in canonical BL_ORDER (BL4 → BL5 → BL3 → BL1 → BL2); within each BL rows are sorted by largest-pool size. **Panel A — Within-location connectivity share** = `largest_pool_N / total_adults`, the fraction of a location's adults that sit in its biggest 50 m mating pool. Horizontal bars run 0.0 → 1.0: **1.0 = fully connected** (whole census in one mating pool); **below 0.5 = majority of adults sit outside the biggest pool** (highly fragmented). Reference dotted lines at 0.5 (red) and 1.0 (grey). This is the single-number per-location fragmentation diagnostic. **Panel B — Mating-pool sizes.** Each dot = one 50 m connected component (= one mating pool), placed at its `component_N_fertile` adult count on the log₂ x-axis; dot size scales with pool size. **Red dotted line: N = 1** = the single-plant SI floor. **Grey dashed line: N = 8 plants = 32 tetraploid allele copies** = the coupon-collector floor for the 32-allele species pool. Row labels: `locationCode_locationID (K pools, X total adults)`. 6/101 pools at or below the SI floor; 31/101 below the coupon-collector floor. Source: `step29d_mating_pool_structure.py`.](figures/Phase5/step29d_mating_pool_structure.png)

### Event-scale companion — pollen-donor reachability per event

Figure 2 above shows every mating pool at every location. The same
spatial data also carries event-scale texture: within a given
location, how uniform is each event's reachable neighbourhood?

**Approach — purely spatial, no allele frequencies.** For every
event, count `N_reachable_50m = Σ N_fertile in other events within
50 m`: the number of pollen-donor plants a flower on that event
could reach at the 50 m pollinator range. **One box plot per
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

**Figure 2 is the primary mating-pool structure; Figure 2b below is
the complementary event-scale texture** — useful when a location's
pool-level summary (Figure 2) hides event-level outliers.

![Figure 2b — Per-location box plots of `N_reachable_50m` (pollen-donor plants reachable within 50 m per event). One row per location, grouped by Bottleneck Lineage; row labels give `(events, adults)`. Each box summarises the location's events' reachable-neighbour counts (log x-axis). **Red dotted line: N = 1** = the single-plant SI floor (nobody to mate with on that event). **Grey dashed line: N = 8 plants = 32 tetraploid allele copies** = the coupon-collector floor for the 32-Fg species pool. Narrow boxes = uniform event neighbourhoods; wide boxes = layered location (some events well-connected, others isolated). Purely spatial — no allele frequencies enter. Source: `step30b_fragmentation_index.py`.](figures/Phase5/step30_A_fragmentation_index.png)

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
   how many SRK alleles are physically present. Driven **per mating
   pool** by `component_N_fertile` (Figure 2), then unioned to the
   location level. This is the Link 2 output of the causal chain and
   the raw material on which the pollen-compatibility prediction
   (§ 30.2) operates.
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

**The underlying biology — sporophytic Class I / II dominance.**
LEPA SRK alleles fall into two dominance classes. **Class I is
dominant**: a plant carrying any Class I allele expresses only its
Class I alleles on both pollen and stigma, silencing its Class II
alleles. A plant carrying only Class II alleles co-expresses all
of them. **A cross is rejected if the two parents share any
expressed allele**; between-class crosses (Class I plant × Class II
plant) are always compatible by construction. This single rule
drives the whole § 30.2 prediction — and crucially, it is what
makes the § C.0.b hypothesis decomposition work: without Class I
dominance, within-class drift concentration (EO70's FG024) and
zygosity composition (EO76's homozygous buffer) would be
indistinguishable from a plain allele-frequency effect. Figure 3b
walks through the rule in three panels before we get into the
simulation steps.

![Figure 3b — Sporophytic SI with Class I / Class II dominance in tetraploid LEPA. **Panel A — dominance within one plant.** Case A: a plant with ≥ 1 Class I allele expresses only its Class I alleles; its Class II alleles are silent (shown faded). Case B: a plant carrying only Class II alleles expresses all four Class II alleles co-dominantly. **Panel B — between-plant recognition, worked example.** Mother M carries {FG001, FG002, FG024, FG031}; her Case-A expressed set is {FG001, FG002}. Three candidate fathers: F1 shares FG001 with M → rejected; F2 is all-Class-II so between-class → always compatible; F3 shares FG002 with M → rejected. **Panel C — compatibility rule by cross type.** Class I × Class I: compatible if their expressed Class I alleles differ (shared Class II is irrelevant because Class II is silent on both sides). Class I × Class II: always compatible by construction (disjoint expressed classes). Class II × Class II: all four alleles expressed on both sides, compatible only if none are shared.](figures/Phase5/step30_A_si_model_schematic.png)

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

**Field → lab correction — 60 % germination rate.** Part C
genotypes **seedlings**, not seeds (user-confirmed design,
2026-10-03). LEPA seeds germinate at **≈ 60 %** under greenhouse
conditions, so to end up with the Rule 2 target of 15 genotyped
seedlings per mother the field protocol must germinate
`ceil(15 ÷ 0.60) = 25` seeds per mother. The totals become:

- **Collect / germinate:** 505 mothers × 25 seeds = **12 625 seeds**.
- **Expected after germination:** ~15 seedlings/mother × 505 = **~7 575 seedlings**.
- **Genotype:** up to 15 seedlings per mother → **7 575 seedling genotypes**.

The coupon-collector allele-detection guarantee still holds —
Rule 2 is a floor on *seedling* genotypes delivered to the lab, not
on seeds extracted in the field. The field-team recipe in
[`step29c_partC_germplasmID_selection.tsv`](tables/Phase5/step29c_partC_germplasmID_selection.tsv)
now carries both `n_seeds_to_germinate` and `n_seedlings_to_genotype`.

**Does 15 seedlings give enough Part C testing power?**
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
`PLOIDY × component_N_fertile = 4·component_N_fertile` allele copies
to its mating pool. Coupon-collector target: 90 % chance of
observing every allele in a given mating pool at the total delivered
allele draws `A_delivered = M × (4 + 2 × 15) = 34·M`. Plus a private-
allele floor: **at least one mother per event** (an isolated
slickspot's private allele cannot be recovered from any other event).

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
| Drift unit | **50 m mating pool** = 50 m connected component of events; `component_N_fertile` adults per pool (Figure 2) |
| Mother allocation per location | Fragmentation-aware `M_frag` (per 50 m mating pool + ≥ 1 per event) |
| Rule 2 seedling-genotype floor | **15 seedlings/mother** |
| Field-side germination assumption | **60 %** |
| Seeds to germinate per mother | **25** (= ceil(15 ÷ 0.60)) |
| **Total mothers across 39 locations** | **505** |
| **Total seeds to germinate** | **505 × 25 = 12 625** |
| **Total seedlings to genotype** | **~505 × 15 = ~7 575** |
| Coupon-collector target | 90 % allele detection per 50 m mating pool |

**Three authoritative files:**

- **Field-team recipe** (for a new field season):
  [`step29c_sampling_frag_aware_per_event.tsv`](tables/Phase5/step29c_sampling_frag_aware_per_event.tsv) —
  one row per event, `M_frag` = number of mothers to sample there.
- **Lab recipe for Part C** (draws specific mothers from the DB):
  [`step29c_partC_germplasmID_selection.tsv`](tables/Phase5/step29c_partC_germplasmID_selection.tsv) —
  one row per selected `germplasmID` already in the LEPA DB, sorted
  `EOID → locationID → component → germplasmID`. Columns include
  `n_seeds_to_germinate` = min(25, seeds_available),
  `n_seedlings_expected` = round(n_seeds_to_germinate × 0.60),
  `n_seedlings_to_genotype` = min(15, n_seedlings_expected).
  Selection is **per 50 m mating pool**: mothers within a pool share
  pollen so coverage travels freely within a pool; only the ≥ 1-
  mother-per-event maternal-genotype floor is a strict per-event
  rule. Delivers **431 mothers from the current DB → ~10 713 seeds
  to germinate → ~6 428 seedlings to genotype**, with **76 mothers
  short across 35 mating pools** flagged for a 2026 field top-up.
- **Event → component lookup**:
  [`step29c_event_to_component_50m.tsv`](tables/Phase5/step29c_event_to_component_50m.tsv) —
  one row per event mapping (locationID, eventID) → 50 m component,
  the single canonical source for "which events share a pollen pool?".

---

## Step 30 Part C § C.0 — Empirical validation of the P_compat model

The primary Part C validation runs at **Phase 5 location scale**
on the three clean-overlap EOs (§ C.0.a and § C.0.b below). The
other three EOs with Phase 4 SRK genotypes — **EO18, EO25, EO27** —
are split under the 500 m rule and cannot yet be remapped to Phase 5
locationCodes. Their observed-vs-predicted pollen compatibility is
reported here at EO scale as a background diagnostic, no figure.

**Current numbers for the three split EOs** (from `step30c_srk_validation_at_eo_level.py`,
observed fathers from local f vs predicted from species-wide P1,
same observed mothers for both):

| EO | adults | Observed pollen compatibility (95 % CI) | Predicted (95 % CI) |
|---|---:|---|---|
| EO27 | 61 | 0.77 (0.73–0.80) | 0.80 (0.77–0.84) |
| EO25 | 51 | 0.77 (0.74–0.80) | 0.80 (0.77–0.84) |
| EO18 | 39 | 0.69 (0.64–0.74) | 0.75 (0.71–0.79) |

EO27 and EO25 pass at EO scale; EO18 shows a mild drift signal
with overlapping CIs. The three clean-overlap EOs (EO67, EO70,
EO76) are reported at Phase 5 location scale in § C.0.a — see the
embedded figure below.

**Caveat.** P1 was built from these same individuals, so species-mean
alignment is guaranteed; what this tests is **robustness to per-EO
drift**, not absolute calibration.

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

![Figure 5 — Phase 5 Part C anchor at the three clean-overlap EOs (EO67, EO70, EO76; 1:1 with a Phase 5 locationCode). Panel order follows the causal chain: diversity → pollen compatibility → per-Fg drift fingerprint. **Panel A — SRK diversity.** Phase 5 predicted (x) vs observed in adults (y), one square per location with 95 % CI horizontal error bars. 1:1 diagonal + ceiling at 32 Fgs. EO67 sits on the diagonal (model passes); EO70 (6/32) and EO76 (9/32) sit far below (large drift gap). **Panel B — pollen compatibility.** Phase 5 predicted vs observed, with the traffic-light background (red = failed, amber = struggling, green = sustainable) and 1:1 diagonal. EO67 and EO76 close to the diagonal; EO70 is the clear outlier (observed 0.53 vs predicted 0.78). **Panel C — per-Fg drift residual** `f_observed − f_P1` per location (one row each, Fgs sorted left-to-right by species-wide P1 frequency, most-common → rarest). **Green bars = Fg enriched vs P1** (drift favoured it); **red bars = Fg depleted vs P1**; **× markers = Fg absent at the location** (lost entirely). The residual pops out the drift fingerprint that the raw-frequency plot blurred — EO70 shows classic FG001 drift (+21 %) + FG024 (+18 %) with 26/32 Fgs absent; EO76 shows milder enrichment of FG012 / FG010; EO67's small-population signature elevates the normally-rare FG018 and FG023 instead of the common Fgs — a founder-effect signature rather than classical drift. Source: `step30d_partC_clean_overlap.py`.](figures/Phase5/step30_B_partC_clean_overlap.png)

### C.0.b Diversity collapse → pollen compatibility: hypothesis decomposition

**Why this analysis exists.** Panel A of § C.0.a (SRK diversity
predicted vs observed) exposed a disconnect that pollen
compatibility alone would have hidden: EO70 and EO76 both lose
~23 of their predicted alleles, but observed pollen compatibility
drops by 0.24 at EO70 and barely 0.03 at EO76. EO67, with just 7
alleles against a predicted 10.8, lands at 0.70 — in the
sustainable band. **Without the diversity trigger from § C.0.a we
would not have known there was a mechanism question to ask.**

**Three competing hypotheses can produce different pollen-
compatibility responses to the same diversity loss.**

1. **Within-class allele spread.** Class I and Class II total
   masses can be identical, but all of Class I's mass may be
   locked in a single allele (drift monoculture), killing
   within-class compatibility (a Class I × Class I cross is
   rejected whenever parents share an allele).
2. **Class I / Class II mass balance.** Raising / lowering the
   Class I share changes how much "between-class rescue" is
   available to Class II mothers (between-class crosses are
   always compatible under the sporophytic model).
3. **Zygosity composition.** The fraction of mothers that carry
   1 / 2 / 3 distinct SRK identities. Multi-identity mothers
   express a bigger set and face more compatible fathers for
   Class II; the effect is more complex for Class I.

#### How the decomposition is built

The same simulation that produced Figure 4 (pollen compatibility
per location) is run **five
times per location**, keeping the observed mother genotypes fixed
and swapping ONE part of the father-drawing distribution at a time
against the species-wide reference. For each scenario, 800
candidate fathers are drawn, sporophytic Class I / II compatibility
is evaluated per mother, averaged over her trials, and then over
all her location's mothers. Bootstrap over 400 resamples of
mothers gives a 95 % CI.

The **only** thing that differs between the five bars for a given
location is the father-drawing distribution:

| Scenario (bar colour)                              | Father allele frequencies            | Father zygosity composition |
|---                                                 |---                                   |---                          |
| **Observed** (black)                               | Observed local                       | Observed local              |
| **Swap within-class spread** (blue)                | **P1 shape within each class, scaled to the observed Class I and Class II totals** | Observed local |
| **Swap Class I / II balance** (red)                | **Within-class shape kept observed, class totals rescaled to the species-wide values** | Observed local |
| **Swap zygosity** (yellow)                         | Observed local                       | **Species-wide 66 / 32 / 2 % of 1-/2-/3-distinct** |
| **Phase 5 prediction** (green)                     | Species-wide P1                      | Species-wide                |

For each single-swap bar, the **driver fraction** = (counterfactual
mean − observed mean) / (species-wide prediction mean − observed
mean). A value near **+100 %** says that one factor alone closes
the entire observed → species-wide gap. A value near **0 %** says
the factor was not involved. A **negative value** means swapping
away from the observed state moves pollen compatibility *away*
from the species-wide prediction — the observed state on that axis
is **buffering** the location.

#### How to read the figure

- **Observed black bar below green**: the location has lost
  pollen compatibility relative to the species-wide prediction.
- **A single-swap bar jumps up toward green**: that factor
  explains most of the gap.
- **A single-swap bar stays next to black**: that factor is not
  involved.
- **A single-swap bar drops below black**: that factor is
  buffering the location — observed state on that axis is more
  favourable than species-wide.

#### Walk-through — the three clean-overlap locations

| Location | Alleles obs / pred | Pollen compatibility obs → species-wide pred | Within-class spread | Class I / II balance | Zygosity composition |
|---|---|---|---:|---:|---:|
| **EO67** | 7 / 10.8 (-3.8) | 0.70 → 0.76 | **+100 %** | −7 % | −13 % |
| **EO70** | 6 / 28.3 (-22.3) | 0.54 → 0.70 | **+88 %** | −3 % | −6 % |
| **EO76** | 9 / 31.8 (-22.8) | 0.75 → 0.77 | **+206 %** | +1 % | **−115 %** |

- **EO70 — classical drift collapse.** FG024 absorbs 35 % of the
  pool, FG001 another 62 %; four other alleles sit under 1 %.
  A FG024-homozygous Class I mother sees `(1 − 0.35)⁴ ≈ 0.17`
  compatibility, vs ≈ 0.43 at EO67 where Class I is split four
  ways. The blue bar jumps up to 0.68 (closing most of the gap);
  red and yellow stay next to black. **Within-class concentration
  explains the deficit; class balance and zygosity are not
  involved.**
- **EO67 — looks bad on paper, pollen compatibility survives.**
  Only 7 of 32 alleles present, but those 7 split as **4 Class I
  alleles (FG024 / FG018 / FG023 / FG016, well spread) and 3
  Class II alleles**. Each Class I mother expresses a smaller
  share of Class I mass than at EO70, so Class I × Class I
  compatibility stays reasonable. Observed zygosity (32 %
  multi-identity) is slightly higher than species-wide — the
  yellow bar sits just below black, meaning observed zygosity is
  mildly buffering. **Mate limitation at EO67 is not as bad as
  the headline 7/32 count suggests.**
- **EO76 — predicted near-perfect, observed severely collapsed,
  pollen compatibility still fine.** 9/32 alleles — a collapse
  comparable to EO70 — but pollen compatibility drops only 0.03.
  Decomposition exposes the mechanism: **zygosity is actively
  buffering EO76** (yellow sits at 0.72, below the black 0.75,
  a −115 % move relative to the gap). The location is 76 %
  homozygous (vs species-wide 66 %), and because its dominant
  Class II allele (FG001 at 46 %) is less dominant than EO70's
  FG001 (62 %), those homozygous Class II mothers express a
  small, non-dominant set and face many compatible fathers. Blue
  overshoots the species-wide prediction (within-class diversity
  at EO76 is actually better-spread than drift-only P1 would
  deliver).

#### Methodological take-home

A location's SRK diversity gap tells us drift has happened; the
pollen-compatibility response is **not** a monotonic function of
how many alleles were lost. It depends on **which** alleles were
lost, **how dominant** the remaining alleles are, and **how multi-
identity** the mothers are. Capturing this requires an accurate
sporophytic Class I / II model plus empirical zygosity (§ A.7–A.8)
— a simpler diploid gametophytic approximation would collapse all
three channels into a single "effective diversity" number and
mis-predict EO76 outright. The hypothesis decomposition therefore
doubles as validation of the model's mechanistic structure: it
successfully resolves cases where the three channels pull in
different directions.

**Workflow placement.** This analysis is a follow-through on § C.0.a's
diversity trigger, not a stand-alone pollen-compatibility test.
When a location's observed allele diversity sits noticeably below
the Phase 5 prediction, run `step30e_pcompat_hypothesis_decomposition.py`
to see which of the three channels is driving the downstream
pollen-compatibility response.

**Outputs.**

- [`step30_B_partC_hypothesis_decomposition.tsv`](tables/Phase5/step30_B_partC_hypothesis_decomposition.tsv)
  — one row per clean-overlap location with the diversity trigger
  (observed vs Phase 5 predicted distinct allele counts), observed +
  three single-swap counterfactual pollen compatibility estimates
  with bootstrap 95 % CIs, observed Class I mass, observed zygosity
  distribution, and the three "driver fraction" columns.
- [`step30_B_partC_hypothesis_decomposition.png`](figures/Phase5/step30_B_partC_hypothesis_decomposition.png)
  — grouped bar chart per location (see caption below).

![Figure 6 — Competing-hypothesis decomposition of per-location pollen compatibility. Triggered by the SRK diversity discrepancies in § C.0.a. For each clean-overlap location (EO67, EO70, EO76) the figure shows five pollen-compatibility values from the same simulation (sporophytic Class I / II + empirical zygosity, 800 candidate fathers per mother, observed mother genotypes held fixed), each differing in which part of the father-drawing distribution is swapped to the species-wide reference. **Black — Observed:** fathers drawn from observed local allele frequencies and observed local zygosity. **Blue — Swap within-class spread:** keep observed Class I and Class II total masses, but reshape the within-class spread to match the P1 pattern (isolates the "drift monoculture" channel). **Red — Swap Class I / II balance:** keep within-class shape observed, rescale the two class totals to the species-wide values (isolates between-class rescue). **Yellow — Swap zygosity:** keep observed allele frequencies, swap father zygosity to species-wide 66/32/2 % (isolates the per-mother expressed-set channel). **Green — Phase 5 prediction:** everything swapped to species-wide. Error bars = bootstrap 95 % CI (400 mother resamples). Traffic-light bands shaded in the background; dashed grey line = species-wide pollen compatibility 0.78. **Reading rules:** a blue/red/yellow bar jumping toward green means that factor caused the gap; staying next to black means it was not involved; dropping below black means it is **buffering** the location (observed state on that axis is better than species-wide). **EO70** — blue dominates (within-class drift monoculture, FG024 at 35 % of pool); red and yellow neutral. **EO67** — small deficit, blue explains it; zygosity slightly buffering. **EO76** — small deficit despite collapsing from 32 to 9 alleles; yellow drops below black, meaning zygosity (76 % homozygous mothers) is actively buffering the location, and blue overshoots green because the specific within-class spread at EO76 is better than drift-only P1 would deliver. Source: `step30e_pcompat_hypothesis_decomposition.py`.](figures/Phase5/step30_B_partC_hypothesis_decomposition.png)

**Note — why P1 is empirical, not uniform.** A uniform-frequency P1
(Dirichlet(α = 1) over 32 Fgs, equivalent to the P0 "uninformative"
prior kept as a reference baseline) would represent a **neutral null
with no drift history**. That is not LEPA: decades of habitat loss
have already pushed the species through drift, so empirical P1 — the
observed Fg frequencies in the Canu-amplicon 263-individual inventory
— is the realistic *starting point* against which per-location
further drift is measured. Panel C makes this concrete: EO70's local
pool has moved *past* P1 toward FG001 dominance, so even the already-
skewed empirical prior underestimates how collapsed the local pool
is. A uniform P1 would start from a flat distribution and declare
every observed location "drifted", which is both less informative
(the species-wide signal is real) and less actionable (the baseline
would not reflect how LEPA actually enters the modelling frame).

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

## Recommended field pilot — BL5 (EO48_7 + EO18-7_19)

Before running Phase B across all 39 locations, we recommend a
**two-location pilot within Bottleneck Lineage 5 (BL5)** — one
**fully-connected** location (**EO48_7**, 98 adults in 1 mating
pool, 9 mothers in DB) paired with one **fragmented** location
(**EO18-7_19**, 34 adults across 3 mating pools, 11 mothers in DB).
User-selected 2026-10-03 (option A2 in
[`step29c_partC_BL5_pilot_candidates.tsv`](tables/Phase5/step29c_partC_BL5_pilot_candidates.tsv),
below).

**Why BL5.** It is the LEPA Bottleneck Lineage with the **widest
within-BL variation** in both census size and mating-pool structure
(see Figure 2): BL5 holds the drift-collapsed tail (EO24 group,
1–3 adult singletons) *and* the largest, most-connected locations
(EO32_6 at 466 adults; EO48_7 at 98 adults in a single mating pool).
A within-BL5 contrast therefore rules out between-BL noise while
testing the full span of the fragmentation × drift axis the model
predicts matters.

**What this specific pair tests.**

- **EO48_7 — fully-connected regime** (`component_N_fertile = 98`,
  single mating pool). Phase 5 predicts a sustainable location
  because `N_fert_eff = total_n_fertile`; this location is the
  cleanest test of the model's "no fragmentation → species mean"
  prediction.
- **EO18-7_19 — fragmented regime** at a similar order of magnitude
  for total adults (34) but split across 3 mating pools.
  Phase 5's per-component simulation treats this as three
  independent drift experiments; the comparison with EO48_7
  isolates the pure fragmentation effect from raw-size effects.

**Caveat on the Phase 4 adult SRK data.** The two BL5 EOs that
*do* have adult SRK genotypes from Phase 4 (EO25, EO18 at n ≥ 10)
are both split under the 500 m rule and sit in § C.0 only — they
cannot yet be used at the Phase 5 location scale. So this pilot is a
**seed-genotyping (Phase B) pilot**, not a retrospective-on-adults
pilot like § C.0.a.

**Pilot cost.** 20 mothers × 25 seeds = **500 seeds to germinate**
→ ~300 seedlings to genotype (at 60 % germination), compared to
~12 600 for the full 2026 design. **< 4 % of the full genotyping
budget** for a within-BL contrast with a clear a priori hypothesis.

**All three options kept on file.**

| Option | Robust candidate | Drift-sensitive candidate | Mothers in DB | Note |
|---|---|---|---:|---|
| A1 | **EO32_6** (5 pools, 466 adults) | **EO25-B_21** (2 pools, 10 adults) | 38 + 7 | Maximum size contrast; both have mothers; drift-sensitive is small-but-not-tiny. |
| **A2 ★** | **EO48_7** (1 pool, 98 adults) | **EO18-7_19** (3 pools, 34 adults) | 9 + 11 | **User-recommended.** Isolates the fragmentation × drift axis — similar order-of-magnitude totals, opposite pool structure. |
| A3 | **EO18-7_17** (5 pools, 242 adults) | **EO24-7_25** (1 pool, 3 adults) | 33 + 3 | Maximum biological contrast but drift-sensitive has only 3 adults — poor statistical power. |

Source table:
[`step29c_partC_BL5_pilot_candidates.tsv`](tables/Phase5/step29c_partC_BL5_pilot_candidates.tsv)
(generated by `build_bl5_pilot_candidates.py`).

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
| Background EO-scale diagnostic for split EOs (§ C.0) | [`tables/Phase5/step30_C_pcompat_validation_at_eo.tsv`](tables/Phase5/step30_C_pcompat_validation_at_eo.tsv) | EO18, EO25, EO27 (the three EOs whose Phase 5 location-scale remap is still pending); no figure |
| **Part C anchor at Phase 5 location scale (§ C.0.a)** | [`tables/Phase5/step30_B_partC_clean_overlap_per_location.tsv`](tables/Phase5/step30_B_partC_clean_overlap_per_location.tsv) + [`_fg_frequencies.tsv`](tables/Phase5/step30_B_partC_clean_overlap_fg_frequencies.tsv) | **Clean-overlap EOs (EO67, EO70, EO76) — 187 adults ready to feed Part C now; observed vs Phase 5 pred P_compat + SRK diversity + no-drift upper bound** |
| **Hypothesis decomposition (§ C.0.b)** | [`tables/Phase5/step30_B_partC_hypothesis_decomposition.tsv`](tables/Phase5/step30_B_partC_hypothesis_decomposition.tsv) | **Decomposes the pollen-compatibility deviation at each clean-overlap location into within-class spread / Class I-II balance / zygosity contributions; triggered by the SRK diversity gap** |

**Key figures** (ordered by appearance in this doc).

- **Figure 1** — [`step30_A_radius_sensitivity.png`](figures/Phase5/step30_A_radius_sensitivity.png) — why 50 m is the right primary pollinator radius.
- **Figure 2** — [`step29d_mating_pool_structure.png`](figures/Phase5/step29d_mating_pool_structure.png) — per-location mating-pool structure: one dot per 50 m connected component (= one mating pool), log x-axis for pool size.
- **Figure 2b** — [`step30_A_fragmentation_index.png`](figures/Phase5/step30_A_fragmentation_index.png) — per-location box plots of per-event pollen-donor reachability (event-scale companion to Figure 2).
- **Figure 3** — [`step30_A_diversity_unbiased_vs_sampling.png`](figures/Phase5/step30_A_diversity_unbiased_vs_sampling.png) — three-panel: what Nature holds (unbiased) vs what our sampling detects vs coverage.
- **Figure 4** — [`step30_A_prediction_fecundation.png`](figures/Phase5/step30_A_prediction_fecundation.png) — predicted pollen compatibility per location, traffic-light bands.
- **Figure 5** — [`step30_B_partC_clean_overlap.png`](figures/Phase5/step30_B_partC_clean_overlap.png) — § C.0.a Part C anchor at Phase 5 location scale (SRK diversity + pollen compatibility + per-allele drift residual).
- **Figure 6** — [`step30_B_partC_hypothesis_decomposition.png`](figures/Phase5/step30_B_partC_hypothesis_decomposition.png) — § C.0.b competing-hypothesis decomposition.

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
