# Phase 5 — SRK sampling & prediction framework (Steps 28–30)

> **Looking for a shorter overview?** A compact companion doc
> [`Phase5_SRK_summary_for_colleagues.md`](Phase5_SRK_summary_for_colleagues.md)
> covers the same framework in ~400 lines with Question → Approach →
> Result blocks per step and one figure per link of the causal chain.
> Share that with collaborators; use the doc you are currently reading
> for methods, formulas, code, and full outputs.

## Executive summary

<span style="color:#777"><strong>Context.</strong></span> *Lepidium papilliferum*
(slickspot peppergrass) is a federally-listed tetraploid
Brassicaceae whose Idaho range spans the **Snake River Plain** and
the Jarbidge Foothills. This work focuses on the **Snake River
Plain populations — 39 spatially isolated locations in documented
decline** — as prioritised by current stakeholder conservation
needs; the Jarbidge Foothills populations are out of scope here
and will be addressed separately. This work investigates the **genomic mechanism
of that decline**: habitat fragmentation shrinks the local
breeding pool, genetic drift erodes the pool of self-incompatibility
(SRK) alleles each plant needs to find a compatible mate, reduced
pollen compatibility surfaces as reduced seed set. The three
questions that follow — how to identify the breeding units, how
to predict allele diversity and pollen compatibility inside each,
and how to test those predictions with seed data — are all
operationalised at the scale of the local breeding unit (**deme**).

<span style="color:#777"><strong>Approach.</strong></span> Because the 39 locations are
already spatially isolated from each other, between-location pollen
flow is approximately zero and **the location is the natural
starting point for identifying local breeding units**. The question
we ask here is whether **within-location further subdivision is
needed** — i.e. whether a single location holds one breeding unit
or several — because SRK allele diversity (how many SRK alleles
physically live inside a breeding unit) and pollen compatibility
(what fraction of pollen × stigma combinations succeed inside it)
are both per-breeding-unit quantities. We follow a Wright-style
genetic-neighbourhood argument (Wright 1943, 1946; Levin & Kerster
1974; Vekemans & Hardy 2004 — see § References) implemented as a
hard-threshold connectivity rule: sweep pollinator radii
(10 – 200 m) across the 39 locations with **2025 field data**
(wild, in-situ, with coordinates — the data scope defined in
§ A.2), pick the **50 m primary radius** ([Figure 1](#fig-1))
where connectivity, sampling cost, and pollen compatibility all
plateau, and read off a **50 m connected component = one
operational deme**. The 50 m threshold is a step-function stand-in
for the (unmeasured) dispersal variance σ²; the resulting
partition is a geographic upper bound on realised gene flow
([Figure 1b](#fig-1b); § A.5).

<span style="color:#777"><strong>Methodology.</strong></span> For each operational deme we
use **empirical allele-frequency data from plant genotyping** to
simulate genetic drift, aggregate SRK allele diversity by set
union across demes to the location ([Figure 3](#fig-3)), and
compute pollen compatibility as a size-weighted mean under the
**sporophytic Class I / Class II model with empirical zygosity**
([Figure 4](#fig-4), [Figure 5](#fig-5)). The deme is both the
hypothesis the pipeline *uses* to generate these predictions and
the hypothesis the next seed-genotyping campaign (see **Next test**
below) will *test*.

<span style="color:#777"><strong>Result + interpretation.</strong></span> **Deme-structure
census ([Figure 1b](#fig-1b)).** The 39 locations collectively
hold **101 operational demes** at 50 m (mean 2.6 demes per
location, range 1 – 6). **13 of 39 locations are a single
connected deme**; the remaining **26 are fragmented into 2 – 6
demes** (10 locations at 2 demes, 6 at 3, 3 at 4, 4 at 5, 3 at 6).
Deme sizes span **1 – 420 adults** (median 23). **31 of 101 demes
(31 %) hold fewer than the 8 plants required to physically carry
the species-wide pool of 32 SRK alleles** (4 alleles per tetraploid
plant × 8 plants = 32 copies), and **6 of 101 (6 %) hold a single
plant with nobody to mate with at 50 m**. So within-location
subdivision is not just real but widespread — two-thirds of the
range operates as several small demes rather than one location-
wide one.

**Prediction vs observation at three locations ([Figure 10](#fig-10)).**
Comparing predictions against observed adult SRK genotypes at
**EO67, EO70, EO76**, the diversity prediction **over-estimates**
observed counts by ~20 SRK alleles everywhere (predicted 11 / 28 /
32 vs observed 7 / 6 / 9) — because the empirical allele-frequency
prior has already absorbed decades of species-wide drift, and
local demes have drifted further on top of it. **Yet the pollen-
compatibility prediction tracks observation closely** — EO67
(observed 0.70 / predicted 0.72) and EO76 (0.72 / 0.78) overlap
on their 95 % CIs; only EO70 shows a real gap (0.53 vs 0.78), a
drift signal consistent with its high FG024 frequency (**0.35
locally, vs 0.18 species-wide**). The hypothesis decomposition
([Figure 11](#fig-11)) explains why: **within-class allele
frequency spread, not allele count, drives pollen compatibility**,
and tetraploid zygosity composition actively buffers locations
with many homozygous mothers. The decisive variable is not *how
many* alleles survive but *how their frequencies and genotypes are
arranged* — which is why a huge diversity gap can coexist with an
accurate pollen-compatibility prediction. **Three mechanisms —
Class I dominance (26 of 32 Fgs), the species's current class
imbalance that keeps common Class I alleles numerous even after
drift, and tetraploid homozygosity — stack into structural
redundancy that lets a location keep breeding after severe allele
loss, reversing the standard allele-count-equals-mating-success
intuition.** **A radius sweep at
EO70 and EO76 ([Figure 13](#fig-13), § C.0.c) tests whether the
diversity gap is just an artefact of the deme being too wide at
50 m.** Rebuilding the deme partition at radii 10 – 150 m leaves
the predicted diversity flat and well above observed at every
radius, including the 10 m extreme that fragments EO76 into 17
small demes and EO70 into 3. The gap therefore **cannot be a pure
spatial-partitioning error** under the current model: the only
way it closes is if each deme's drift history has diverged from
the species-wide prior — i.e. each deme carries its own frequency
vector that the species-wide prior does not capture, which the
next seed-genotyping campaign is designed to measure. The 50 m
operational deme is defensible as a first-pass partition, and the
direction of error is favourable: over-predicting at a radius
that is already a geographic upper bound says real demes are
**at most** 50 m and could be tighter, so a future refinement can
only add demes (and therefore mothers), giving us finer-scale
information for free. (EO67 is a null — its deme partition is
invariant across the sweep and its observed count is already
inside the 95 % CI of the prediction.)

<span style="color:#777"><strong>Next test.</strong></span> We propose a between-location pilot
(§ C.6) that pairs a **big + connected anchor** with a **small +
fragmented partner**, chosen from across the full LEPA dataset so
that the pair spans the usable predicted pollen-compatibility
spectrum. Three options are on the table — all share the same
anchor (**EO29_8**, BL1: 417 adults in a single deme, 22
mothers in DB, predicted pollen compatibility 0.779):

| Option | Pair kind | Drift-sensitive partner | Predicted pollen-compatibility gap | Size ratio | Mothers (anchor + partner) | Trade-off |
|:---|:---|:---|:---:|:---:|:---:|:---|
| **B1 ★** | within-BL1 | **EO26-2_35** (22 adults, 3 demes, share 0.50, 8 mothers) | **0.039** | 19× | 22 + 8 = 30 | Strongest confound control (both BL1 → same evolutionary context) and largest within-BL predicted gap. |
| B2 | within-BL1, bigger size gap | EO26-3_34 (34 adults, 3 demes, share 0.71, 14 mothers) | 0.025 | 12× | 22 + 14 = 36 | More mothers → less sampling noise, but less fragmented partner and smaller contrast. |
| B3 | between-BL | EO25-B_21 (BL5: 10 adults, 2 demes, 7 mothers) | **0.051** | **42×** | 22 + 7 = 29 | Largest predicted gap and biggest size ratio, but confounds fragmentation × BL identity. |

**B1 is recommended** (★): it maximises within-BL contrast
sharpness while holding BL identity constant. Full pair evaluations
in [step29c_partC_between_location_candidates.tsv](Tables/Phase5/step29c_partC_between_location_candidates.tsv).
Per-mother sampling recipe (germplasmIDs to pull from the LEPA DB,
with per-deme allocation, seeds-to-germinate and seedlings-to-
genotype columns, 60 % germination correction baked in) in
[step29c_partC_germplasmID_selection.tsv](Tables/Phase5/step29c_partC_germplasmID_selection.tsv).

## Contents

The document is organised in **three parts** that follow the causal
order of the study — we first build the model and generate location-
level predictions, then derive the sampling protocol from those
predictions, then test the predictions with real seed data:

- [Scientific goals](#scientific-goals) — the central hypothesis (a causal chain: fragmentation → drift → mate limitation) and how the framework tests it
- **Part A** — [Model and predictions](#part-a--model-and-predictions) · data scope, species-wide P1 prior, finite-population model, and the Phase A predictions of SRK diversity and pollen compatibility per location (Step 30 Phase A outputs, Figures 4–6). *This is what we expect each location to look like — before we ever open a seed lot.*
- **Part B** — [Sampling protocol derived from the predictions](#part-b--sampling-protocol-derived-from-the-predictions) · within-location pollen connectivity, per-mother seed count (Step 28), per-location mother count with private-allele floor (Step 29), and two-year design. *This is what the field team must do to test the Part A predictions.*
- **Part C** — [Testing predictions with observed SRK data (Phase B)](#part-c--testing-predictions-with-observed-srk-data-phase-b) · preliminary EO-level validation with Phase 4 adult genotypes (§ C.0), mate-limitation regression, script behaviour, data-generation pipeline, and the BL5 pilot (§ C.6).
- [Appendix — Out of scope / future work](#appendix--out-of-scope--future-work) · SI-escape permutation-test scaffold retained for pipeline-validation purposes only; not part of the Phase 5 scientific analysis.
- [Output map](#output-map--quick-reference-grouped-by-phase) — filenames organised by phase.

**Naming.** Parts **A / B / C** refer to *sections of this document*.
Output filenames use `step30_A_*` for Phase A prediction artefacts (built
without seed data) and `step30_B_*` for Phase B artefacts (built with
observed seed data); `step30_B_DEMO_*` marks synthetic Phase B for
pipeline validation. Steps 28–29 outputs are all Phase A by construction
and keep their existing `step28_*` / `step29_*` names.

## Scientific goals

### The central hypothesis — a causal chain

The purpose of this framework is not diversity estimation for its own
sake. It tests a **single causal hypothesis about how habitat
fragmentation degrades reproduction in *Lepidium papilliferum***:

> **Fragmentation → genetic drift → mate limitation.**
>
> Spatial isolation of adult plants (fragmentation of pollen flow)
> shrinks the effective deme at each location. Small effective
> demes intensify genetic drift on SRK allele frequencies,
> which erodes local SRK diversity and skews local Fg composition.
> The eroded and skewed local pool reduces the probability that any
> given pollen grain is compatible with any given mother
> (**mate limitation**), which in turn reduces per-mother seed set.

The chain has four measurable links, and the framework operationalises
them in the order they must be evaluated:

| # | Link | Where measured | What it produces |
|---|---|---|---|
| 1 | **Fragmentation** — spatial isolation of adult plants at each location | § A.5 fragmentation indices; § B.2 within-location connectivity (Step 29b) | `N_fertile_effective_50m` = census × largest-connected-component share at 50 m; per-event / per-location fragmentation indices `F_event`, `F_location` |
| 2 | **Genetic drift on SRK** — small effective demes lose alleles and skew Fg frequencies | § A.6 predicted SRK diversity per location (Step 30 Phase A) | Predicted local Fg pool size and predicted per-EO local Fg frequencies (with 95 % CI), driven by `N_fertile_effective` |
| 3 | **Mate limitation** — the drifted local Fg pool determines random-mating pollen compatibility under sporophytic Class I / II SI | § A.7 SI biology, § A.8 finite-population P_compat (Step 30 Phase A) | Predicted per-location `P_compat` (with 95 % CI); § C.0 empirical validation on adult SRK genotypes |
| 4 | **Reduced seed set** — the reproductive consequence at the mother | § C.1 mate-limitation regression (Phase B, needs seed genotypes) | Observed per-mother seed set regressed on predicted `P_compat` and on the direct fragmentation term |

**Reading order — why fragmentation must be evaluated first.**
Fragmentation is the physical driver upstream of everything else. It
determines `N_fertile_effective`, which is the *only* place where
location size enters the drift-diversity-compatibility chain. This is
why Step 29b (spatial connectivity) runs before Step 30 (drift +
compatibility predictions), and why § A.5 (fragmentation) precedes
§ A.6 (diversity) and § A.7–A.8 (compatibility) in the doc. A location
with 1 000 census adults but only 60 in its largest 50 m deme
is drift-limited as if it were a 60-plant location, and the whole
downstream chain (loss of Fgs → lower P_compat → lower seed set) is
what the pipeline predicts on that reduced N.

### What the Phase B regression will test — and how it distinguishes the two mechanisms

Once seed genotypes are available (Phase B), the § C.1 regression
`seeds_per_mother ~ β₁ · predicted_P_compat + β₂ · mating_neighbourhood + …`
puts numbers on the chain:

- **β₁ (P_compat effect)** measures the strength of the
  drift-diversity-compatibility branch — the *complete chain* from
  fragmentation through seed set.
- **β₂ (mating-neighbourhood effect)** measures whether fragmentation
  has *additional* direct effects on seed set that are NOT mediated
  through Fg composition (e.g. reduced pollinator visitation because
  neighbouring adults are too sparse to draw pollinators, independent
  of which alleles they carry).

The two coefficients together decompose the fragmentation effect into
"acting through drift" (β₁) and "acting through other pathways" (β₂),
and can be pre-registered before seed data arrive.

### Complementary goal — phenotype cross-validation (deferred)

SRK-based predictions can be cross-checked against per-site ISI /
fruit set in the Genetic-Rescue-DB repository. This adds independent
lines of evidence but does not shape the sampling design and is not
modelled here in v1.

### Not addressed by this framework

Self-incompatibility escape (plants in which SI has broken down
entirely, producing viable self-seed) is *not* an outcome of Phase 5.
Such individuals are identified during the Canu-amplicon SRK
genotyping (Phase 4 Step 22b) and enter Phase 5 as prior information
via the empirical zygosity distribution (§ A.7.3a), not as an
experimental target. A DEMO SI-escape permutation-test scaffold
exists in the Step 30 code for pipeline-validation purposes only —
it is documented in the "Out of scope / future work" appendix at
the end of this doc.

### Design consequences

The causal chain drives two design decisions: **(1)** the sampling
protocol must scale with each mother's *real* mate-availability
context — her fragmentation-adjusted local deme, not the raw
census (Steps 28–29 build this into `N_fertile_effective`); and
**(2)** the inference layer must publish testable per-location
predictions of drift-loss and compatibility that can be regressed
against observed seed set (Step 30 in `prediction` mode). Part A
publishes the predictions; Part B derives the sampling protocol
needed to test them; Part C runs the test.

---

## Key concepts and terminology

Every per-location quantity in this doc is derived from a nested
spatial hierarchy. The two top levels are **quoted verbatim from the
LEPA DB `Terms` table** (the canonical glossary that ships with
`LEPA_SQL.db`) so the vocabulary is consistent with every other LEPA
analysis. Reading from the biggest unit down to the individual plant:

- **Location** (`locationID` / `locationCode`) — the DB's `Locations`
  table defines `locationID` as "Location Unique Barcode #" (Darwin
  Core `dwc:locationID` — "An identifier for the set of dcterms:Location
  information"), with `locationCode` holding the EO code of the
  sampling site ("Report the unique EO # where the sampling is
  conducted (e.g., EO38)"). Locations may span several slick spots
  across a landscape.
  - **Within-EO location split (Phase 5 refinement).** Several EOs
    contain events that sit **≥ 500 m apart with no bridging events
    in between**. Those spatially disjoint pieces of a single EO are
    tracked as **separate `locationCode`s** suffixed with a dash
    (e.g. EO24 → `EO24`, `EO24-1`, `EO24-2`, `EO24-7`; EO27 →
    `EO27`, `EO27-1`, `EO27-3`, `EO27RT`), each with its own
    `locationID`. In the current data **5 EOs are split this way**
    (EO18, EO24, EO25, EO26, EO27 → 16 Phase 5 locationCodes among
    them). The rationale: 500 m exceeds LEPA's primary pollinator
    radius (50 m) by an order of magnitude, so these sub-locations
    cannot share pollen under any plausible flight distance and must
    be modelled as independent drift units. This refinement is
    layered **on top of** the DB's canonical Location — it does not
    rewrite it, and the DB's `locationID` barcode is preserved for
    every sub-location.
  - **Event** (`eventID` / `occurrenceID`) — the DB's `Events` table
    defines an event as "**an 'Event' refers to an occupied slick
    spot within a Location**" (Darwin Core `dwc:eventID` — "A unique
    identifier for the event (e.g., a field survey collecting
    *Lepidium papilliferum* in a slickspot)"). Each event has its
    own census of fertile plants and its own coordinates in the DB.
    - **50 m connected component** (`component_id_50m`) — a group
      of events whose plants sit within 50 m pollinator-flight range
      of each other. **Plants in the same component share a pollen
      pool; plants in different components — even at the same
      location — do not.** This is a Phase 5 derived concept, not a
      DB term.
      - **Mother plant** (`germplasmID`) — an individual plant
        already collected and stored in the LEPA DB, sitting at one
        specific event.

Everything else in the pipeline is a count or a derived number
sitting on top of this hierarchy:

| Concept (full English name) | Code identifier | Rooted at | Definition |
|---|---|---|---|
| Fertile plant census | `total_n_fertile` | location | All fertile plants at a location, summed across every event. The biological potential, with no spatial filtering. |
| **Effective deme size** (also called **N_fertile_effective**) | `N_fert_eff` | **component (primary); location (as a diagnostic)** | **Per component:** `component_N_fertile`, the fertile plants that share a single 50 m pollen pool — **this is the drift unit that drives per-component SRK diversity and per-component pollen compatibility.** **Per location:** `total_n_fertile × largest_component_share_50m`, a one-number *fragmentation diagnostic* used in Figure 2 to flag locations whose raw census is split across multiple components. The two numbers agree when the whole location is one component. |
| Connectivity share | `largest_component_share_50m` | location, built from components | `N_fert_eff (location) ÷ total_n_fertile`. 1.0 = fully connected (no fragmentation); 0.3 = 70 % of the raw census is drift-irrelevant. **Fragmentation diagnostic only — the actual Phase A prediction loops over every component.** |
| Fragmentation-aware mother target | `M_frag` | event, derived from components | For each event, the number of mothers to sample so that each 50 m connected component reaches 90 % allele-detection coverage, with a ≥ 1-per-event maternal-genotype floor. Sums across events to the location-level `M_frag_aware`. |
| Tetraploid per-mother seedling floor | 15 seedlings/mother | mother plant | Each genotyped seedling contributes 2 paternal allele draws from the local pollen pool. 15 seedlings/mother is the per-mother floor that gives a 90 % chance of seeing every allele in her deme. At the project's 60 % greenhouse germination rate (§ B.3.2) the field collection target is `ceil(15 ÷ 0.60) = 25` seeds per mother. See § B.3. |
| Species prior | `P1` | species-wide | The **32 Fgs identified across LEPA plus their empirical species-wide frequencies** — a 32-slot probability vector summing to 1, built from the Canu-amplicon L1 carrier inventory. FG001 fills 41 % of the vector, the next five Fgs ~25 % together, and the 20+ rare Fgs each < 2 %. "Draw an allele from P1" means pick a Fg with probability equal to its slot in this vector. Every per-location prediction draws from P1, so small locations lose the rare Fgs to drift by chance. |

**Why components matter in one sentence.** Every per-location
prediction in this doc — SRK diversity, pollen compatibility, mother
allocation — is built **component-by-component**, because a 50 m
connected component is what trades pollen. Each component is
simulated with its own `component_N_fertile`; location-level numbers
are the **set union** across components for diversity (Fgs are a
set — present *somewhere* in the location), and the **size-weighted
mean** across components for pollen compatibility (a continuous
rate — Σ_c P_compat_c · N_c / Σ_c N_c). A location with 500 fertile
plants spread across 20 isolated slickspots behaves like 20 small
drift-prone pools, not one pool of 500.

### How the five quantities chain together

Every number in Part A is derived in the same five-step order,
starting from the spatial hierarchy above:

1. **Location → events.** A location's raw census
   `total_n_fertile` is the sum of per-event `n_fertile` counts
   across all its 50 m-resolved events ([`step28_events_spatial_neighborhood.tsv`](tables/Phase5/step28_events_spatial_neighborhood.tsv)).

2. **Events → 50 m components.** For each location, build a graph on
   its events: two events are connected if their fertile plants sit
   within ≤ 50 m. Connected components of that graph are the 50 m
   components `component_id_50m` ([`step29c_event_to_component_50m.tsv`](tables/Phase5/step29c_event_to_component_50m.tsv),
   written by Step 29c; the same connectivity logic is reported at
   location scale in Step 29b).

3. **Components → effective deme size.** For each component *c*,
   `component_N_fertile_c` = Σ_{e ∈ c} n_fertile_e. These are the
   plants that actually share one 50 m pollen pool. The location-level
   diagnostic `N_fertile_effective_50m = total_n_fertile ×
   largest_component_share_50m` is Figure 2 only; the predictions
   never collapse to it.

4. **Effective deme → SRK diversity prediction (per component).**
   Draw 4 × component_N_fertile_c SRK alleles i.i.d. from the species-
   wide prior P1 → the component's drift-collapsed Fg set.
   **Location SRK diversity = |union of per-component Fg sets|**
   (a set operation — a Fg present in *any* component is present at
   the location). Sampling coverage is built the same way: mothers
   and seeds are allocated across components proportional to
   `component_N_fertile_c`; per-component sampling is simulated
   against each component's own frequency vector; the location
   detected count is the union across components. Outputs:
   [`step30_A_prediction_location_diversity.tsv`](tables/Phase5/step30_A_prediction_location_diversity.tsv)
   + [`step30_A_prediction_component_diversity.tsv`](tables/Phase5/step30_A_prediction_component_diversity.tsv).

5. **Effective deme → pollen compatibility prediction (per
   component).** From each component's local Fg frequency vector,
   simulate mother genotypes under the empirical LEPA zygosity (66 %
   single-identity, 32 % 2-distinct, 2 % 3-distinct) and compute
   each mother's sporophytic Class I / II P_compat against candidate
   fathers drawn from the same component pool.
   **Location P_compat = size-weighted mean of per-component
   P_compat**, Σ_c P_compat_c · component_N_fertile_c / Σ_c
   component_N_fertile_c (a continuous rate, so size-weighted mean,
   not union). Outputs:
   [`step30_A_prediction_location_pcompat.tsv`](tables/Phase5/step30_A_prediction_location_pcompat.tsv)
   + [`step30_A_prediction_component_pcompat.tsv`](tables/Phase5/step30_A_prediction_component_pcompat.tsv).

The whole causal chain of this framework — fragmentation → genetic
drift → mate limitation — enters the predictions through step 3
(effective deme size per component) and comes out at steps 4
and 5. A fragmented location at step 2 produces small pools at
step 3, which drives drift-collapsed Fg sets at step 4 and lower
pollen compatibility at step 5.

---

## Part A — Model and predictions

Part A builds the model and publishes the location-level predictions
that everything else in this framework is designed to test. The chain
is: **species-wide prior → finite-population draw at each location →
predicted SRK diversity → predicted pollen compatibility → cross-plot
that summarises the fragmentation → drift → mate-limitation chain**.
No seed data are used here. These predictions can be pre-registered
today; Part B derives the sampling protocol required to test them,
and Part C runs the tests once real seed genotypes exist.

### A.1 Two operational phases

The pipeline runs in two operational phases, corresponding to the state of
the seed genotype dataset. **All steps use LEPA_SQL.db and the P1
empirical prior in both phases**; what changes is whether observed seed
genotypes are available:

| Phase | State of data | Uses | Produces |
|---|---|---|---|
| **A — Preliminary** (before SRK genotyping) | Field data only (`Germplasm.germplasmQuantityEstimate`, `Events.organismQuantityFertile`, event coordinates), plus the P1 species-wide prior | Steps 28, 29, and Step 30 in `prediction` mode | Per-mother sampling recipe + per-location seed-count recommendations, **plus** predicted SRK diversity, pollen compatibility, fecundation failure |
| **B — Post-genotyping** (after SRK data are back) | Everything above + observed seed genotypes (from Phase A's sampling) | Step 30 in `comparison` mode + Test 1 | Observed vs predicted SRK diversity, mate-limitation regression — all only for the locations that have observed data |

**Filename conventions make the phase — and its data provenance —
unambiguous.** Every Step 30 output uses one of three prefixes:

| Prefix | Meaning | Produced by |
|---|---|---|
| `step30_A_*` | **Phase A** — prior-based prediction; always safe to produce | Default `python step30_...py` |
| `step30_B_*` | **Phase B** — result from *real* observed seed genotypes | `--seed-genotypes real.tsv --mother-genotypes real.tsv` |
| `step30_B_DEMO_*` | **Phase B, synthetic** — simulated seeds for pipeline validation | `--demo` |

Every Phase B figure produced under `--demo` carries a diagonal
**"DEMO" watermark** and an explicit `[DEMO — synthetic data]` suffix
in its title, so screenshots and slide captures cannot be mistaken for
real analysis.

If observed data exist for only a few locations, Phase B outputs report
on those few; the Phase A prediction outputs still cover every location
the field data know about.

### A.2 Data scope — what enters the pipeline

Every SQL query in Steps 28–30 applies two mandatory filters and one
optional filter. **All scripts accept a `--year YYYY` flag with the
same semantics.**

- **Mandatory · Wild, field-collected only.** Occurrences and germplasm
  records tied to ex-situ material (nursery accessions, seed-increase
  plots, in-vitro cultures) are excluded from every table and figure. The
  two SQL predicates are:
  - `Germplasm.biologicalStatus = 'Wild'` — keeps 765 records, drops 20
    ex-situ.
  - `Occurrences.provenance = 'in situ'  OR  IS NULL` — keeps 2 419 + 552
    records, drops 808 ex-situ and 1 in-vitro.
  This is not a per-run switch — greenhouse plants have no place in a
  wild-population mate-limitation and fragmentation analysis, so the
  filter is baked in.

- **Mandatory · Coordinates present and inside LEPA's known bounding
  box.** Rows with NULL or free-text `eventDecimalLatitude` /
  `eventDecimalLongitude` are dropped, and remaining rows are clipped to
  30–55 °N × −125 to −100 °E to prevent a mistyped decimal from pushing
  the spatial neighbourhood off-planet.

- **Optional · Year filter (`--year YYYY`).** Restricts to events whose
  `eventDate` belongs to the given year. When set, it applies to *both*
  the germplasm records (Step 28) and the events feeding the spatial
  neighbourhood (so K_spatial never mixes survey years — see § B.5).
  The current DB contains 765 wild seed-bearing mothers in 2025 and no
  other year; running without `--year` and with `--year 2025` therefore
  give identical results today, but the flag is what will keep the two
  seasons cleanly separated once 2026 field data arrive.

The framework has three steps, each with a distinct role and a distinct
output family. The table below is **ordered by the causal flow**, not by
historical script numbering — prediction first, then the sampling
design derived from it:

| Step | Question | Kind of output |
|---|---|---|
| **Step 30 (Phase A)** | Given the species-wide prior and the per-component effective deme size, what SRK diversity and pollen compatibility do we predict per 50 m component and per location? | Prediction (per-component + per-location) |
| **Step 29** | How many events per location, and how many mothers per event, are needed so that each 50 m connected component reaches the 90 % allele-detection target under those predictions? | Sampling design (per-event / per-location) |
| **Step 28** | Given the per-mother allele exposure already set by Steps 30 and 29, how many seeds per mother must the lab genotype to recover the within-mother pollen pool at 90 % probability? | Sampling design (per-mother) |
| **Step 30 (Phase B)** | Once observed seed genotypes exist, how do they compare to the Phase A prediction, and does per-mother seed set decline with predicted pollen compatibility? | Comparison + mate-limitation regression |

### A.3 Species-wide prior P1 and ploidy

> **Ploidy note — LEPA is TETRAPLOID (2n = 4x).** Every somatic plant
> carries **4 SRK allele copies**; meiotic reduction produces **2x
> pollen grains** with 2 SRK alleles per grain, and 2x eggs with 2
> alleles per egg. Every seed therefore carries **2 maternal + 2
> paternal SRK alleles** (not 1 + 1). This is baked into every
> count-based formula in Parts A and B via a single project-wide
> constant `PLOIDY = 4` (see `step28_seed_sampling_per_mother.py`).
> **The sporophytic self-incompatibility model** with Class I / Class II
> dominance now drives § A.7's pollen compatibility formula ([`srk_si_model.py`](srk_si_model.py)
> + [`tables/Phase5/srk_fg_class.tsv`](tables/Phase5/srk_fg_class.tsv))
> and Phase C's β₁ interpretation follows the new sporophytic scale.

Every Phase A prediction below is built on a **species-wide empirical
prior** over the 32 functional SRK allele groups (Fgs) rather than a
flat / uninformative prior. The Canu-amplicon preliminary study
already characterised **32 Fgs in 263 individuals**, which gives Phase A
enough information to publish predicted diversity and predicted
fecundation failure per location, with credible intervals, *before* any
seed data land. When the seed data arrive (Part B → Part C), prediction
and observation are compared side by side and the locations where they
disagree are the actionable ones.

We build three nested Dirichlet priors over the 32 Fgs:

| Prior | Source | Use |
|---|---|---|
| **P0 — Uninformative** | Dirichlet(α = 1) over the 32 Fgs | Reference baseline; what "no prior knowledge" looks like |
| **P1 — Species-wide** | 32-Fg carrier counts from `step26i_L1_carrier_inventory.tsv` | Workhorse prior for locations we haven't visited |
| **P2 — Spatially informed** | Fg counts pooled from spatial neighbours (BL / EO) | Location-specific prior when we know the neighbourhood |

Each prior has an *effective sample size* controlling how strongly it
constrains inference; larger ESS means the prior is only shifted by large
seed datasets, smaller ESS means observed genotypes dominate quickly. **P1
is the workhorse** used by every prediction below.

### A.4 Finite-population model and prediction methodology

Every Part A prediction (A.5–A.8) comes from the same **per-50 m-
component finite-population simulation**, aggregated to the location
level in the way the quantity biologically supports. The 50 m
connected component is the drift unit — plants in different components
of the same location do not share pollen, so they do not share a drift
history.

1. **Simulate each component's local SRK pool (drift signature).**
   For each 50 m connected component *c* with `component_N_fertile`
   adults, draw `4 × component_N_fertile` alleles i.i.d. from the
   Canu-amplicon species-wide prior P1. Small components lose rare
   alleles to drift; large components approach the species-wide 32.

2. **Predicted SRK diversity.**
   - **Per component:** number of *distinct* Fgs in that component's
     simulated pool.
   - **Per location:** |union of per-component Fg sets| — a Fg
     present in *any* component is present at the location.
   - **Sampling coverage** (per component): mothers and realised seeds
     are distributed across components proportional to
     `component_N_fertile` (largest-remainder integer allocation);
     each component's `A_delivered_c = 4·M_c + 2·seeds_c` draws are
     taken against its own frequency vector and the detected Fgs
     recorded. Location coverage = |union detected| / |union present|.

3. **Predicted pollen compatibility.**
   - **Per component:** simulate mother genotypes from the
     component's local Fg frequencies under the empirical LEPA
     zygosity distribution, compute each mother's sporophytic Class
     I / II P_compat against candidate fathers drawn from the same
     component pool (see § A.7–A.8).
   - **Per location:** size-weighted mean of per-component P_compat,
     Σ_c (P_compat_c · component_N_fertile_c) / Σ_c
     component_N_fertile_c. (Weighted mean, not union — P_compat is
     a continuous quantity rather than a set.)

4. **Uncertainty.** Both diversity and P_compat are recomputed on
   400 – 1 000 independent per-component replicates. Reported means
   and 95 % credible intervals are the mean and 2.5 % / 97.5 %
   percentiles across those replicates. Small components have wide
   intervals (founder-effect uncertainty is large); large components
   have tight intervals (local pool converges to P1).

**Sampling effort — what the numbers assume we will see.** *M* mothers
× 15 seeds each (tetraploid per-mother cap) = **A = M × 34 allele draws** from
the local pool (4 alleles per mother from her own tetraploid
genotype + 2 paternal alleles × 15 seeds = 30 paternal draws). We
report coverage in two flavours:

- **Species-wide coverage** — fraction of the 32 P1 alleles detected.
  Biased low for small locations (drift already removed most).
- **Local coverage** — fraction of the Fgs *actually at the location*
  (union across its components) that we detect under the per-component
  sampling allocation. The biologically honest metric; drops below
  99 % when drift-collapsed small components receive 0 mothers under
  proportional allocation and still contribute Fgs the larger
  components do not carry (see § A.6).

**One-line take-home:** diversity is what drift has left and the
sampling has detected, pollen compatibility is how well a random
mother matches the neighbours she can reach at 50 m, and both are
simulated **per component** against the same P1 prior; location-level
numbers are the union (diversity, a set) or the size-weighted mean
(pollen compatibility, a continuous rate) of the component-level
results.

**Prediction outputs.** Kept in `tables/Phase5/` prefixed
`step30_A_prediction_*` and figures under `figures/Phase5/step30_A_*`.
They can be produced today, without any seed data at all. The
corresponding *comparison* outputs (Phase B, once seed genotypes exist)
carry a different prefix (`step30_B_*`) so prediction and comparison
cannot be mistaken for one another.

### A.5 Fragmentation of pollen flow — population, deme, event scales

**A three-level spatial hierarchy for an annual plant.**
*L. papilliferum* is an annual whose above-ground presence at a
slickspot can switch on and off between years depending on seed-
bank germination and growing-season conditions. A single-year
"location" footprint is therefore unstable: a slickspot can show
plants in 2025 and not in 2026, or vice versa, while the
**underlying seed-bank patch persists**. To get a spatial frame
that does not flicker with the above-ground emergence, Phase 5
distinguishes three nested units — a stable **population** that
spans years, a within-year **deme** that is the drift unit, and
the raw per-year **event**:

| Level | Operational definition | Biological meaning | Reference |
|---|---|---|---|
| **Species** | the 32 Fgs that make up the species-wide prior P1 | *L. papilliferum* on the Snake River Plain | this doc + § A.3 |
| **Population** (metapopulation unit) | 500 m-separated connected components, built from the **pooled 2025 + 2026 event coordinates + historical locationID centroids** | A demographically independent spatial unit. 500 m is >> contemporary pollen flight (≈ 10 × the primary pollinator radius) and >> LEPA seed dispersal (gravity + short-distance wind in low-stature slickspot vegetation), so between-population gene flow on an ecological timescale is negligible. The population is **the stable patch** over the project timescale — the above-ground emergence may be zero in a given year while the seed-bank patch persists. | Levins 1969; Hanski 1998; Freckleton & Watkinson 2002 |
| **Deme** (local breeding unit) | 50 m connected component **within one year's event set**, nested inside a population | Operational deme / Wright's genetic neighbourhood — the within-year pollen pool and drift unit on which all SRK allele diversity and pollen compatibility predictions are built | Wright 1943, 1946; Levin & Kerster 1974; Vekemans & Hardy 2004 |
| **Event** | one occupied slick spot in one year | Darwin Core `dwc:eventID` — the raw observation record | DB definition |

Every quantitative output in Phase 5 — predicted SRK allele
diversity, predicted pollen compatibility, mother-sampling
allocation — is a **per-deme** quantity before it aggregates
upward. The deme is the biological unit of inference; the
population is the stable spatial frame it nests inside; the radius
choices below are the operational knobs that delineate each level.

#### A.5.0 Populations — the stable spatial frame (Step 29a)

**Why we need a population level above the deme.** The 50 m deme
is a **within-year** unit — it uses the pollinator-flight radius to
define who could have exchanged pollen *in a given season*. For a
multi-year study of an annual plant, that leaves a gap: how do
you compare drift at the "same location" across years when the
above-ground footprint shifts? Classical metapopulation theory
answers this by separating the **patch** (geographic unit that
persists) from the **local population** (the biological community
occupying it in a given year, which can collapse to the seed
bank). Phase 5 adopts the same split: the **population** is the
patch, the **deme** is the within-year local population inside it.

**Why 500 m — three independent justifications converge.**

1. **≈ 10 × the primary pollinator radius (50 m).** Contemporary
   pollen gene flow across a 500 m gap under the halictid /
   small-bee foraging envelope is effectively zero. The 50 m
   deme-level connectivity itself plateaus by 75 m
   ([Figure 1](#fig-1)), so a gap an order of magnitude wider
   is firmly outside the ecologically-realised neighbourhood.
2. **>> LEPA seed dispersal.** Brassicaceae in low-stature arid
   vegetation disperse seeds predominantly by gravity + short-
   distance wind on the order of **metres**. 500 m is well beyond
   the realised seed shadow — long-distance dispersal events are
   rare and demographically irrelevant on the project timescale.
3. **Consistent with the existing within-EO 500 m rule.** Phase 5
   already uses 500 m as the within-EO location-splitting threshold
   (§ Spatial hierarchy — DB terms, 5 EOs split → 16 Phase 5
   locationCodes). Promoting that rule to the formal population
   threshold removes a loose end rather than adding a new parameter.

**Operational algorithm.** Note that **eventIDs are NOT stable
across years** — the field crews assign fresh barcodes each
season, so there is no explicit "same slickspot visited again"
link in the DB. We therefore need a cross-year **slickspot-matching
step** based on coordinates + a buffer that absorbs GPS drift
before the population graph is built.

1. Load all **2025 + 2026 LEPA events** (coords + `n_fertile`;
   same data-scope filters as § A.2), tagged with their source
   year.
2. Add the **historical centroid** of every DB-level locationID
   (= mean lat/lon across every event ever recorded at that
   locationID). This anchors the population on the DB-level
   location identity even when both 2025 and 2026 above-ground
   emergence happened to be zero at some slickspot.
3. **Cross-year slickspot matching** — build a **10 m** haversine
   graph on the pooled 2025 ∪ 2026 event set; connected components
   of that graph are **slickspots** (`slickspotID`). The 10 m
   buffer absorbs typical consumer-GPS error (3–5 m open sky, up
   to ~10 m in poor conditions). Each slickspot is tagged
   `both_years` / `2025_only` / `2026_only`.
4. Pool slickspot centroids + locationID centroids into the final
   point set, each row tagged with its source.
5. Build the **500 m haversine connectivity graph** on that set
   — an edge wherever two points sit within 500 m.
6. **Connected components = populations.** Assign a reproducible
   `populationID` to each component (sorted by total `n_fertile`
   across both years DESC, then westmost lon ASC).
7. Write a two-level crosswalk TSV
   (`step29a_population_crosswalk.tsv`) that maps every
   `(year, eventID) → slickspotID → populationID`, preserving the
   DB `locationID` and Phase 5 `locationCode` in extra columns so
   legacy references in older tables and figures do not break.
8. **Nested step within each population × year**: rebuild the 50 m
   deme partition on that year's events only (Step 29c logic, now
   scoped to `(populationID, year)` rather than Phase 5
   locationCode).
9. **Seed-cleaning priority ranking.** The Phase 5 collaborators
   are still cleaning 2026 seed and extracting germplasm; this
   step feeds that queue. For each population, score by (a)
   presence in **both** 2025 and 2026 above-ground samples, (b)
   the number of same-slickspot re-visits, (c) total `n_fertile`.
   The output flags **one large and one small candidate** —
   populations with ≥ 2 re-visited slickspots in both years and
   consistent above-ground trend — as the preliminary across-year
   analysis targets.

**What a population IS and ISN'T.**

- **IS** — a geographic patch that persists across the 2025 + 2026
  project window and the historical DB record; the stable frame
  for cross-year comparison; the between-population scale beyond
  which contemporary pollen and seeds do not flow.
- **ISN'T** — a genetic deme (that's the within-year 50 m level);
  nor a strictly panmictic unit (the demes inside it can still be
  internally fragmented at 50 m); nor a management unit (BLs
  remain the management stratum).

**Downstream consequences.** The old Phase 5 locationCode becomes
a legacy label kept in the crosswalk for traceability but no
longer the inferential unit. Current `(39 locationCodes, 101 demes)`
re-expresses as `(N populations, M_2025 + M_2026 demes)`. SRK
diversity and *P*<sub>compat</sub> become **year-resolved per
population** (one prediction per `populationID × year`), which
turns the drift story into a measurable **temporal** signal:
within one population, the 2025 and 2026 above-ground samples are
two independent draws from the same seed bank, and systematic
drift between them is a direct estimate of per-population drift
beyond P1 (the H2 scenario of § C.0.c, now observable in the
project's own data rather than needing a decades-long historical
baseline).

**Script.** `step29a_populations.py` builds the crosswalk +
population-structure summary; `step30g_population_year_predictions.py`
runs the per (populationID, year) SRK diversity + pollen
compatibility simulation; `step29a_population_maps.py` draws the
overview + candidate zoom figures.

**Dataset-wide results (first run, 2025 + 2026 pooled).** Site
visitation effort was equal across years; 2026 adds many more
occurrence records (mother plants collected for seed banking) that
Phase 5 does not use at this stage, so raw `N_fertile` sums are
the appropriate across-year axis.

- **44 populations** total across the Snake River Plain, built
  from **514 slickspots** (514 cross-year matched components out
  of 644 raw events; 10 m buffer) + **52 historical locationID
  centroids**.
- **22 populations present in both 2025 and 2026** (above-ground
  emergence); **11 only in 2025**; **11 only in 2026**. 0
  populations are purely historical (every population is anchored
  by at least one year's above-ground events).
- **38 of 514 slickspots (7.4 %) were re-visited in both years**;
  the other 476 are single-year — the raw signal of the annual
  plant's year-to-year above-ground flickering.
- **258 per-year demes** across the 44 populations (vs 101 under
  the previous single-year-single-locationCode framework).
- **Legacy locationCode mapping.** 21 of 30 Phase 5 locationCodes
  present in the 2-year data map 1:1 to a population. 2
  populations aggregate ≥ 2 legacy locationCodes (the clearest
  example is populationID = 1, which merges EO27-1 and EO27RT —
  the 500 m-separation rule now connects slickspots that used to
  be split, consistent with the metapopulation framing).

**Across-year prediction behaviour.** `step30g` ran the per
(populationID, year) simulation on all 44 × up-to-2 = 66
population-year combinations; the 22 both-year populations yielded
across-year comparisons:

- **13 stable** (both 2025 and 2026 95 % CIs overlap on both SRK
  diversity and pollen compatibility).
- **9 CI-mismatched** on SRK diversity — all 9 are driven by
  deme-size differences year-to-year; **0 of 22 populations show
  a pollen-compatibility CI mismatch**, reflecting the structural
  redundancy documented in § C.0.b.
- **4 crash candidates** — populations where `N_fertile_2026` ≤
  25 % of `N_fertile_2025` with negative Δdiversity. Striking
  highlight: **populationID = 10 (locationCode EO76)** — one of
  the three § C.0.a clean-overlap EOs — collapsed from 445 to 16
  plants between years. If the 2026 genotyping confirms the same
  severely-collapsed allele pool observed at EO76 in 2025 (9/32
  Fgs), that is **direct second-year evidence for the H2 scenario
  of § C.0.c** — the diversity gap is not a one-year sampling
  fluke but a persistent per-deme drift signal.
- **Preliminary LARGE + SMALL candidate pair for the across-year
  analysis** (`step30g_stable_candidate_pair.tsv`):
  - **LARGE = populationID 1** (locationCodes EO27-1 + EO27RT),
    `N_fertile` 548 → 1239, predicted diversity 31.9 → 32.0 at
    the species ceiling, predicted pollen compatibility 0.779 →
    0.776 (stable).
  - **SMALL = populationID 35** (locationCode EO67),
    `N_fertile` 10 → 12, predicted diversity 10.9 → 11.7,
    predicted pollen compatibility 0.759 → 0.752 (stable at a
    tiny scale).

**Figures.** Overview + candidate panels + crash panels + across-
year scatter are produced by `step29a_population_maps.py` and
`step30g_plot.py` — see output map at the end of this doc.

**Pollinator-radius choice.** Every fragmentation and downstream prediction in Phase 5 depends on the assumed pollen-flight radius. A dedicated sensitivity sweep across 10, 25, 50, 75, 100, 150, 200 m ([Figure 1](#fig-1)) demonstrates that **50 m is the biologically sound primary radius**: it sits within the halictid / small-bee foraging literature range, captures 43 % of the sampling-cost reduction available on the radius curve, and retains meaningful fragmentation variation across BLs (median connectivity 0.81, not yet saturated at 1.00 like at 75 m+). Every downstream analysis in this doc uses 50 m as the primary radius; 10 m and 25 m are always computed and stored as sensitivity checks.

<a id="fig-1"></a>
![Figure 1](figures/Phase5/step30_A_radius_sensitivity.png)

**Figure 1.** Full pollinator-radius sensitivity sweep across 10, 25, 50, 75, 100, 150, 200 m — the biological justification for adopting **50 m as the primary pollinator radius**. **Panel A** (landscape connectivity) shows the fraction of adults in a multi-event pollen-flow component; the median across locations rises from 0.00 at 10 m to 0.45 at 25 m to 0.81 at 50 m and plateaus at 1.00 by 75 m. **Panel B** (§ B.4.2 sampling cost) shows the total mothers required across all 39 locations; it drops steeply from 887 at 10 m to 505 at 50 m (43 % reduction) then flattens to 328 at 200 m. **Panel C** (§ A.8 pollen-compatibility prediction, empirical zygosity) is essentially flat across all radii at ~0.78 — because ~66 % of LEPA plants are single-identity homozygotes, radius-driven changes in effective N barely move P_compat. **Panel D** shows that all 39 locations stay in the "sustainable" band regardless of radius. **The two curves that DO change with radius are connectivity (Panel A) and sampling cost (Panel B); both stabilise around 75–100 m. 50 m sits within the halictid / small-bee foraging literature range (10–100 m), captures 43 % of the sampling-cost reduction available on that curve, and retains meaningful fragmentation variation across BLs (0.81 median, not yet saturated at 1.00 like at 75 m+). Source: `step29a_pollinator_radius_sensitivity.py`. Data: [`step30_A_radius_sensitivity_summary.tsv`](tables/Phase5/step30_A_radius_sensitivity_summary.tsv), [`step30_A_radius_sensitivity_per_location.tsv`](tables/Phase5/step30_A_radius_sensitivity_per_location.tsv).

**What this predicts and why it is separate from A.5–A.7.** A.5
(predicted SRK diversity) and A.7 (predicted pollen compatibility)
both live on the **drift axis**: they answer *"what alleles has drift
left at this location, and how compatible is a mother against those
local frequencies?"* This section adds the **spatial-flow axis**: how
much of the surviving diversity a mother can actually **reach** given
where LEPA events sit on the landscape and how far the pollinators
fly. Both are goal-3 predictions, but they answer different questions
— *"what alleles remain?"* (A.5) vs *"which alleles can pollinate
her?"* (A.8) — and they are **independent**: a fragmented location
can still be drift-neutral (many alleles, none reachable) and a
well-connected location can still be drift-collapsed (few alleles,
all reachable). Splitting them here lets Phase C (§ C.1) test which
channel is dominant at each location as two independent coefficients
(β₁ = drift via pollen compatibility, β₂ = fragmentation via K^(25m)).

**Why this metric is honest.** The fragmentation index uses only three
inputs: event coordinates, N_fertile per event, and the 50 m primary
pollinator radius. It uses **no allele frequencies, no priors, no
simulations** — so it can be published today, independent of any
SRK genotype work, and it depends on no modelling choice beyond the
50 m radius (already justified in § B.2 with the 10 m / 50 m
sensitivity view).

**Event scale — the mother-level view.** For each event we compute

$$K^{(25\text{m})} \;=\; 4 \times \Bigl( \sum_{e \in \text{neighbours}_{\leq 25\text{m}}} N_{\text{fertile}}(e) \,-\, 1 \Bigr)$$

the number of pollen-donor SRK allele copies a mother at that event
can reach within pollinator range. This is exactly the K used in
§ B.3 for the allele-detection formulas, viewed here as a
**fragmentation predictor** at the mother scale.
The **distribution of K^(25m) across events within each location**
is captured per location by the sum of its event-scale neighbourhood
counts and summarised at the location scale by F_location (below)
and by the deme structure in Figure 1b.
In the 2025 field, **439 / 704 events (62 %)** have at least one
other event within 50 m; the number of reachable donor plants per
event ranges from **0** (isolated singletons) to **> 250** (dense
BL3, BL4 clusters), and under tetraploid this corresponds to
K^(25m) allele copies from **0** to **> 1 000**.

**Location scale — the site-level view.** At the location scale
fragmentation is the fraction of adults NOT in a multi-event
pollen-flow component,

$$F_{\text{location}} \;=\; 1 - \text{connected-share at 50 m}$$

using the connected-share metric from § B.2. `F_location = 0` when
every adult exchanges pollen with at least one other event; `F_location = 1`
when every adult is an isolated island. In the 2025 field, **21 / 39
locations (54 %) sit at F_location ≥ 0.5** — the majority of LEPA
sites are structurally fragmented at pollinator scale. The
location-scale summary is shown in Figure 7 (§ B.2); [Figure 1b](#fig-1b)
decomposes each location into its 50 m demes so the per-
component structure is explicit.

**Patterns at the location scale** that the Figure 1b decomposition
makes visible:

- **Well-connected locations** — a single large deme well
  above the 8-plant species-pool threshold; every adult belongs to
  the same pollen pool. Examples: EO30-1, EO27-3, EO29, EO70.
- **Mixed-connectivity locations** — several demes of
  varying sizes within one location, so mothers at the same
  locationID face radically different pollen contexts. Examples:
  EO27-1, EO18-7, EO26-3, EO8 (several sub-locations). These are
  the locations where goal 3's β₂ (fragmentation) coefficient in
  Phase C will have the most within-location leverage.
- **Uniformly isolated locations** — a single deme at
  N = 1 – 3 adults; under strict SI no seed set is possible from
  that pool alone. Examples: EO24, EO24-1, EO24-2, EO24-7.
  These are the sites where fragmentation is predicted to be *fully*
  driving mate limitation.

**How A.8 hooks into Phase C.** The per-mother K^(25m) is the β₂
predictor in the mate-limitation regression (§ C.1); the per-location
F_location is available as an explicit location-level covariate. Both
are pure spatial metrics — no drift contamination — so β₂ isolates
the fragmentation channel cleanly. Locations with wide event-scale
K^(25m) spread (mixed-connectivity boxes above) contribute the most
within-location statistical leverage for identifying β₂ separately
from between-location random effects.

**Take-home for reviewers.** Fragmentation is a distinct, testable
prediction. It is what a pollinator sees, not what a SRK allele
frequency looks like. Publishing it alongside A.5–A.7 gives goal 3
two clean prediction axes (drift, fragmentation) before any seed
genotype exists — and Phase C then quantifies which axis is dominant
at each location.

#### A.5.1 `N_fertile_effective` — the pivotal metric

The fragmentation quantities above (F_event, F_location, event-scale
K^(50m)) are useful diagnostics, but they are all downstream of a
single derived population number that carries the causal chain into
drift and mate limitation:

```
N_fert_eff = total_n_fertile × largest_component_share_50m
```

**What it means.** `N_fertile_effective` is the number of adults that
actually share a pollen environment at each location, given the 50 m
primary pollinator flight range. The **raw census** at a location
tells us how many fertile plants exist there; `N_fertile_effective`
tells us how many of them are close enough to each other to trade
pollen. It is the *drift-relevant* population size: the value that
determines how strongly random loss and skew act on local SRK allele
frequencies (§ A.6; the drift step of the causal chain) and therefore how
random-mating pollen compatibility drops at fragmented locations
(§ A.8; the mate-limitation step of the causal chain).

**How it differs between locations.** [Figure 2](#fig-2) shows
`N_fert_eff` per location alongside the raw census and their ratio
(the connectivity share). Three patterns emerge:

- **Well-connected locations** (connectivity share ≈ 1). EO24 group,
  EO118, EO70, EO29, EO26-3 (small), EO67, EO72-2, EO26-4 — every
  fertile plant sits in the same 50 m deme. Raw census =
  `N_fert_eff`. Small BL5-tail locations sit here by default because
  a location of 1–3 plants cannot be fragmented further.
- **Partially fragmented locations** (connectivity share 0.5 – 0.9).
  EO76, EO32, EO18-7, EO61, EO8 clones — a large census but the
  50 m connectivity is imperfect, so 10 – 50 % of the raw census is
  invisible to drift. EO61 (raw 543 → `N_fert_eff` 315) loses 42 %
  of its adults to fragmentation.
- **Heavily fragmented locations** (connectivity share < 0.5).
  EO27-1 (raw 371 → 116, share 31 %), EO27-1 (raw 395 → 147, share
  37 %), EO18-8 (raw 123 → 52, share 42 %), EO26-2 (raw 22 → 11,
  share 50 %). Half or more of the raw census does not contribute
  to drift-relevant N. These are the locations where a census
  number would systematically overstate the *effective* deme.

**How it is used downstream.** `N_fert_eff` is the single connectivity-
adjusted number that every Phase A prediction relies on:

- **Figure 3 (§ A.6)** — Panel A (unbiased local Fg diversity) is a
  allele-detection draw of `PLOIDY × N_fert_eff` alleles from P1.
  Small N_fert_eff → small local pool → drift-collapsed diversity.
- **Figure 4 / Figure 5 (§ A.8)** — pollen compatibility is
  simulated on a local Fg pool sized by `PLOIDY × N_fert_eff`, with
  M mothers drawn from that pool. Small N_fert_eff → skewed local
  frequencies → mothers more likely to face fathers with overlapping
  expressed alleles.
- **§ C.1 mate-limitation regression** — β₁ (P_compat effect) fires
  through the entire chain rooted in `N_fert_eff`; β₂ (mating-
  neighbourhood effect) directly uses `K^(50m)` per mother, which is
  the per-mother analogue of `N_fert_eff` at the individual scale.
- **§ B.4.2 fragmentation-aware sampling** — the reason 22 / 39
  locations need MORE mothers than the pooled Step 29 recommendation
  is the same `largest_component_share_50m` shrinkage exposed in
  Panel C of Figure 2.

**Take-home.** Wherever the doc previously said "N_fertile" as a
biological input, it means `N_fertile_effective`. The raw census is
Nature's biological potential; the effective count is what actually
matters for reproduction under 50 m pollinator flight.

#### Deme structure per location — the drift unit (Step 29c)

The 50 m choice above is the knob that controls every downstream
fragmentation metric. Before summarising the per-location structure
as the single number `N_fertile_effective` (Figure 2), it is worth
seeing the raw structure that `N_fert_eff` collapses: within each
location, how many pollen pools does a plant belong to, and how big
is each one?

##### The deme as a working hypothesis

Everything the pipeline predicts rests on how a **deme** is
delineated: **SRK allele diversity is a per-deme count** (how many
Fgs physically live in the deme), and **pollen compatibility is a
per-deme rate** (what fraction of pollen × stigma combinations
succeed inside the deme). Both aggregate to the location by set
union and size-weighted mean respectively (§ A.6, § A.8). **Mate
limitation at a location — the ultimate target of this framework —
is therefore inherited directly from the deme-level predictions.
Get the deme wrong and the whole causal chain is wrong, so
delineating demes correctly is paramount.**

The theoretical scaffolding is Wright's genetic-neighbourhood
concept (Wright 1943, 1946; Levin & Kerster 1974; Vekemans & Hardy
2004), which identifies the pollen-flight scale as the biologically
meaningful partitioning scale within a spatially structured
population.

We do **not** estimate Wright's neighbourhood parameter
*N*<sub>b</sub> itself — that would require parent-offspring
dispersal distances or a fine-scale *F*<sub>ST</sub> ~ distance
curve (Hardy & Vekemans 1999; Vekemans & Hardy 2004), neither of
which is currently available for *Lepidium papilliferum*. Instead
we delineate **operational demes** using a hard-threshold
connectivity rule: two events are joined in the location graph
when any pair of their plants sits within 50 m, and each
connected component is one deme. The 50 m radius is a
step-function stand-in for the (unmeasured) dispersal variance σ²,
calibrated by the sensitivity sweep in which connectivity, cost,
and *P*<sub>compat</sub> all plateau at ≥ 75 m. The partition is
therefore a **geographic / topological upper bound** on realised
gene flow — the actual flow across a 50 m gap could be lower than
the physical bound suggests, but not higher.

**The deme definition is a testable working hypothesis.** The
Part C seed-genotyping design (§ B.3–B.4, § C.0–C.1) is built to
compare the per-location *P*<sub>compat</sub> prediction against
observed per-mother compatibility recovered from seed genotypes.
A systematic mismatch at the location level is **not a failure of
the framework** — it is information about **how good the 50 m
operational deme is as an estimator of the realised deme**. A
future phase of work could then refine the deme definition (shrink
σ² if gene flow is tighter than 50 m suggests; extend the radius
if cross-component flow is appreciable; add pollinator-behaviour
weights to the connectivity rule) and close the gap. The 50 m
operational deme is both the hypothesis the pipeline *uses* to
generate predictions and the hypothesis the Part C data will
*test* (see § References for the full literature on
fragmentation-driven drift and operational-deme delineation).

##### How the deme is built

**Definition.** A **50 m connected component** is the Phase 5
**deme**: the set of adult plants whose events are reachable
from each other at the primary pollinator radius. Plants in the same
deme share pollen; plants in different demes at the
same location do not. The deme is the **drift unit** on which
every Phase A prediction is built.

**Build.** For each location, construct an event-level graph in
which two events are connected if any pair of their plants sits
within ≤ 50 m. Connected components of that graph are the location's
demes. For each deme *c*, `component_N_fertile_c` is
the sum of `n_fertile_e` across its events — the adult count that
drives the per-component drift simulations in § A.6 and § A.8.
`step29c_fragmentation_aware_sampling.py` writes the event →
deme lookup (`step29c_event_to_component_50m.tsv`);
`step29d_mating_pool_structure.py` builds the display below.

**Dataset-wide result — how many demes do the 39 LEPA locations
hold?** The connectivity rule yields **101 operational demes
across 39 isolated locations** at 50 m. The per-location breakdown
is informative:

| Demes at a location | Locations | Interpretation |
|:---:|:---:|:---|
| 1 | **13** | single connected deme; the location-wide census and the breeding unit coincide |
| 2 | 10 | two sub-demes separated by a > 50 m internal gap |
| 3 | 6 | three sub-demes |
| 4 | 3 | four sub-demes |
| 5 | 4 | five sub-demes |
| 6 | 3 | six sub-demes — the most fragmented locations |

So **~⅓ of locations (13 / 39) operate as a single location-wide
deme**, and **~⅔ (26 / 39) require within-location subdivision**
into 2–6 operational demes to represent gene flow correctly.
Deme sizes span **1 adult (SI floor) → 420 adults** with a median
of **23 adults**. **31 of 101 demes (31 %)** sit below the 8-plant
species-pool threshold (fewer than enough plants to physically
carry 32 SRK alleles), and **6 of 101 (6 %)** sit at the N = 1
single-plant SI floor (a plant with nobody to mate with at 50 m).
The BL5 tail (EO24 group) concentrates the small demes; BL1 holds
the widest *within-location* spread (EO8, EO26-3 fragmented across
many small sub-demes). Figure 1b plots every deme and location
in this distribution.

**Conclusion for the modelling framework.** Within-location
subdivision is widespread, not an edge case. Treating each
location as a single pool would collapse 101 drift units into 39
and would misclassify two-thirds of the range — enough to drive
qualitatively different predictions of SRK allele diversity and
pollen compatibility. From § A.6 onward the drift unit is therefore
the 50 m operational deme, not the location.

**Downstream dependence.** Every per-location prediction below is
derived from this structure: Figure 1b summarises the deme
inventory per location (count and sizes), and Figures 3 and 5 run
per-component simulations whose aggregation back to the location
depends on this deme decomposition.

<a id="fig-1b"></a>
![Figure 1b](figures/Phase5/step29d_mating_pool_structure.png)

**Figure 1b.** Deme structure per LEPA location. Two aligned panels; one row per locationID (unified `{EOID}_{locationID}` label), rows stacked vertically and grouped by Bottleneck Lineage in canonical BL_ORDER (BL4 → BL5 → BL3 → BL1 → BL2); within each BL rows sorted by largest-deme size. **Panel A — Within-location connectivity share** = `largest_pool_N / total_adults`, the fraction of a location's adults that sit in its biggest 50 m deme. Horizontal bars 0.0 → 1.0: 1.0 = fully connected; below 0.5 = majority of adults outside the biggest deme. Reference dotted lines at 0.5 (red) and 1.0 (grey). Single-number per-location fragmentation diagnostic. **Panel B — Deme sizes.** Each dot = one 50 m connected component (= one deme), placed at its `component_N_fertile` adult count on the log₂ x-axis; dot size scales with deme size. Red dotted line: N = 1 = single-plant SI floor. Grey dashed line: N = 8 plants = 32 tetraploid allele copies = 8-plant species-pool threshold for the 32-allele species pool. Source: `step29d_mating_pool_structure.py`. Data: [`step29d_mating_pool_summary.tsv`](tables/Phase5/step29d_mating_pool_summary.tsv).

### A.6 Predicted SRK diversity per location

**Two questions, kept strictly separate.** SRK diversity at a
location has two very different meanings, and the pipeline predicts
both. This section uses the **whole existing LEPA dataset** — actual
M mothers × actual seeds/mother recorded in the DB — to describe
what our real data show. The prospective design question ("how many
mothers × how many seeds do we NEED?") is answered separately in
Steps 28 and 29 (§ B.3 – § B.4).

1. **What the population actually holds at this location** — the *unbiased*
   local Fg diversity. Depends only on the effective deme size
   `N_fertile_effective = N_fertile × largest-connected-component share at 50 m`
   (§ B.2). Not conditioned on sampling. This is the number that
   feeds drift → mate limitation of the causal chain: the local drifted Fg
   pool composition is what the P_compat prediction (§ A.7–A.8) is
   built on. It is what we are trying to characterise.
2. **What our sampling will recover** — the *sampling-inferred*
   detection under the existing LEPA DB.
   `A_delivered = 4·M + 2·(total seeds recorded at the location)`,
   where each seed's 2 paternal alleles are direct samples of the
   **local pollen donor pool** — the seed genotyping is a pollen-pool
   characterisation experiment. Total per-location seed counts in
   the real DB span **1 (EO24-2, 1 mother) → 11 723 (EO76, 62 mothers)**.

**Coverage = sampling / unbiased** measures how well the existing
LEPA dataset recovers the biological truth at each location.

**What the population holds (the drift step of the causal chain — unbiased).** Under
the tetraploid P1 finite-population model, the expected number of
distinct Fgs physically present at a location is a with-replacement
draw of `PLOIDY × N_fertile_effective` alleles from P1:

- Large, well-buffered slickspots (**EO30-1, EO29, EO76, EO61, EO32**
  with `N_fertile_effective` ≥ 300) hold **~31–32 of the 32
  species-wide Fgs** — drift is not yet a concern at that census
  size.
- Mid-sized locations (`N_fertile_effective` 100–200) hold
  **~26–30 Fgs** — drift has begun to erode diversity but the
  local pool still covers most of the species pool.
- Small locations (`N_fertile_effective` < 50) hold **~15–22 Fgs**.
- The BL5 tail (**EO24, EO24-1, EO24-2, EO24-7** with 1–3 plants
  each) holds **only ~3–6 Fgs**: drift has already collapsed the
  local pool to a handful of alleles. This is the drift step of the causal
  chain firing in the population — the P_compat prediction (§ A.8) at these
  locations then reflects the mate-limitation consequence of that
  collapse.

**What our sampling will recover (sampling-inferred).** Applying
`A_delivered_actual` from the existing LEPA DB to the allele-detection
calculation against each location's true local pool:

| Location | BL | N_fert_eff (50 m) | M_mothers | total_seeds | Population holds | Sampling detects | Coverage |
|---|---|---|---|---|---|---|---|
| EO24-2 (1 plant, BL5 tail) | BL5 | 1 | 1 | 1 | ~3 | ~2.6 | **88 %** ⚠ |
| EO24 (BL5 tail) | BL5 | 2 | 1 | 5 | ~4.5 | ~4.1 | **91 %** ⚠ |
| EO24-1 (BL5 tail) | BL5 | 2 | 2 | 10 | ~4.5 | ~4.5 | 99 % |
| EO26-3 (small BL1) | BL1 | 4 | 4 | 20 | ~6.6 | ~6.5 | 99 % |
| EO67 (small BL4 pilot) | BL4 | 6 | 4 | 71 | ~8 | ~8 | ~100 % |
| EO27 (BL4) | BL4 | 173 | 15 | 1251 | ~29 | ~29 | ~100 % |
| EO30-1 (BL4) | BL4 | 420 | 23 | 2249 | ~32 | ~31.5 | ~100 % |
| EO27-1 (large BL4 pilot) | BL4 | 116 | 33 | 2465 | ~27 | ~27 | ~100 % |
| EO32 (well-sampled BL5) | BL5 | 324 | 38 | 6028 | ~31 | ~31 | ~100 % |
| EO29 (BL1) | BL1 | 417 | 22 | 5419 | ~32 | ~31.5 | ~100 % |
| EO61 | BL1 | 315 | 55 | 9820 | ~31 | ~31 | ~100 % |
| EO76 (largest BL3) | BL3 | 401 | 62 | 11723 | ~32 | ~31.5 | ~100 % |

**Every LEPA location clears the 90 % target and almost all reach
≥ 99 % coverage.** The only two locations noticeably below full
recovery are single-plant slickspots in the BL5 tail — **EO24-2**
(1 plant × 1 seed = 6 allele observations against a ~3-Fg pool)
and **EO24** (1 plant × 5 seeds). Their coverage is capped by the
seed lot's size, not by mother sampling: EO24-2 has only one plant,
so no more mothers exist to add.

The existing LEPA dataset therefore characterises Nature's truth at
every location. The prospective sampling design in § B.3 – § B.4
carries forward the tetraploid per-mother cap of 15 seeds as
a floor for new field seasons.

**Take-home for reviewers.** The unbiased local Fg pool (Panel A of
Figure 3) is the biological signal — the *raw material* on which the
sporophytic P_compat prediction (§ A.8) operates and therefore the
quantity that carries the mate-limitation story. The sampling-inferred
detection (Panel B) confirms that our field-team recipe recovers
that signal at ≥ 95 % coverage at all but two locations (Panel C).
The BL5 tail's low diversity in Panel A is not an artefact of
sampling — it is what drift has already done to those slickspots.

<a id="fig-3"></a>
![Figure 3](figures/Phase5/step30_A_diversity_unbiased_vs_sampling.png)

**Figure 3.** Predicted SRK allele diversity per LEPA location under the tetraploid P1 finite-population model, **built per 50 m component and aggregated to the location by set union** (same drift unit as Figure 5). Panelled by Bottleneck Lineage in canonical BL_ORDER (BL4 → BL5 → BL3 → BL1 → BL2, top-to-bottom); within each BL row, locations sorted by unbiased pool size (small → large). Y-axis labels give `locationCode (N_fert_eff, M_mothers, total_seeds)` where `N_fert_eff` = sum of component_N_fertile across the location's 50 m components, `M_mothers` = mothers in the LEPA DB, `total_seeds` = total seeds recorded at the location across all mothers (raw count, not a mean — per-mother seed counts can vary widely at the same location). **Panel A** — What the population actually holds. For each 50 m component, draw 4 × component_N_fertile alleles from P1 and record its present Fgs; location pool size = |union across components|. Not conditioned on sampling. Feeds drift → mate limitation of the causal chain. **Panel B** — What our sampling will detect. Mothers and seeds are distributed across components proportional to component_N_fertile (largest-remainder); per-component sampling draws `A_delivered_c = 4·M_c + 2·seeds_c` alleles against the component's own frequency vector, and the location detected count = |union of per-component detected Fgs|. **Panel C** — Coverage fraction = Panel B ÷ Panel A. Dotted line = 90 % target. Error bars = 95 % credible interval across simulation replicates. Two single-plant BL5 slickspots (EO24-2, EO24) sit at ~87–91 % coverage — now joined by small-component locations where drift and proportional allocation leave 1–5-plant components with zero sampled mothers (6/39 locations <99 % coverage). Vertical dashed line in Panels A/B = species-wide ceiling of 32 Fgs. Per-component rows: [`step30_A_prediction_component_diversity.tsv`](tables/Phase5/step30_A_prediction_component_diversity.tsv). Source: `step30_srk_diversity_prediction_vs_observed.py`.


### A.7 Sporophytic self-incompatibility with Class I / Class II dominance

This section states the biological SI model; § A.8 is the
mathematical implementation and § C.1 the Phase C test that uses
its predictions. [Figure 4](#fig-4) below is a purely pedagogical
schematic that summarises the whole section in three panels — the
dominance rule inside one plant (§ A.8.3), the worked example of
one mother vs three candidate fathers (§ A.8.5), and the
compatibility rule by cross type (§ A.8.4).

<a id="fig-4"></a>
![Figure 4](figures/Phase5/step30_A_si_model_schematic.png)

**Figure 4.** Sporophytic SI with Class I / Class II dominance in tetraploid LEPA (Phase 5 § A.7). **Panel A — dominance within one plant.** Case A (plant with ≥ 1 Class I allele): only its Class I alleles are expressed on both pollen and stigma; Class II alleles are silent (shown faded). Case B (plant with only Class II alleles): all four Class II alleles are expressed co-dominantly. **Panel B — between-plant recognition, worked example.** Mother M carries `{FG001, FG002, FG024, FG031}`; her Case-A expressed set is {FG001, FG002}. Three candidate fathers: F1 shares FG001 with M → rejected; F2 is all-Class-II so between-class → always compatible; F3 shares FG002 with M → rejected. **Panel C — compatibility rule by cross type.** Class I × Class I: compatible if their expressed Class I alleles differ (Class II silent on both sides — sharing them is irrelevant). Class I × Class II: always compatible by construction (disjoint expressed classes). Class II × Class II: all four alleles expressed on both sides, compatible only if none are shared. Source: `si_model_schematic.py`.

#### A.7.1 Sporophytic recognition

Self-incompatibility rejects pollen whose SRK identity matches the
stigma's, preventing self-fertilisation. LEPA is a Brassicaceae, so
its SI is **sporophytic**: SRK recognition is determined by the
**sporophyte generation** (2n) rather than the haploid gamete (n).
The pollen coat carries proteins deposited during pollen development
by the sporophyte tapetum, so the stigma decides against the whole
pollen *parent*'s expressed genotype rather than the individual
gamete's allele. In LEPA the sporophyte is **tetraploid (2n = 4x)**,
so the pollen parent contributes 4 SRK alleles' worth of expression
information (subject to the dominance rules in § A.8.3). Compatibility
depends on both parents' expressed genotypes.

#### A.7.2 Two allelic classes

Brassicaceae S-haplotypes fall into two dominance classes:

- **Class I** — the dominant class. Large, highly diverged allele
  cluster — in Brassica species Class I typically carries the
  majority of distinct S-haplotypes.
- **Class II** — the recessive (or co-dominant among themselves)
  class. Smaller, more conserved allele cluster with fewer distinct
  members.

**Provisional Phase 5 mapping — a known shortcut.** Without a
sequence-based phylogenetic assignment against Brassica Class I /
Class II reference S-haplotypes, we assign LEPA's 32 Fgs to a class
on **count ratio alone** — 26 Class I + 6 Class II — so the mapping
matches the "Class I has more alleles" biological constraint:

| Class | Fgs | Sub-alleles per Fg | Total P1 frequency |
|---|---|---|---|
| **Class I (provisional)** | 26 (FG007–FG032) | 1 | ~35 % |
| **Class II (provisional)** | 6 (FG001–FG006) | 2 – 10 | ~65 % |

The mapping is stored as a first-class TSV
([`srk_fg_class.tsv`](tables/Phase5/srk_fg_class.tsv)); rerunning
`step30_srk_diversity_prediction_vs_observed.py` regenerates every
downstream number from the current version.

**On the P1 frequency pattern.** In Brassica species without strong
drift, Class I alleles are usually the *more common* alleles — the
dominant class carries the majority of allele frequency as well as
allele count. In LEPA that pattern breaks down: FG001 alone sits at
41 % of P1, and the six multi-allele Fgs (now Class II in the
provisional mapping) carry ~65 % of P1 frequency. The reason is
**genetic drift**: LEPA populations are small and fragmented
(§ A.5), so the "typical" Brassica frequency pattern is not the
null expectation here. Allele counts, which are preserved under
drift, are the more reliable classification signal — and
`srk_fg_class.tsv` is built on counts for exactly that reason.

**To do — proper phylogenetic reassignment.** The current mapping
should be replaced with a sequence-based classification: for each
of the 32 Fgs, align its SRK S-domain against reference Brassica
Class I (e.g. B. rapa S52) and Class II (e.g. B. rapa S29, S40)
haplotypes, assign each Fg to the class it clusters with. This is
a one-off analysis; results will supersede the provisional mapping
and every downstream P_compat prediction. Tracked as a future
work item.

#### A.7.3 Dominance within a tetraploid plant

LEPA carries 4 SRK allele copies (§ A.3). Under the classical
Brassicaceae rule:

- Class I strictly dominates Class II. A plant with ≥ 1 Class I
  allele expresses only its Class I alleles on both stigma and
  pollen; Class II alleles are silent.
- An all-Class-II plant expresses all four alleles co-dominantly.
- Within a class, alleles are co-dominant (a within-Class-I
  hierarchy for LEPA, if any, is not yet resolved).

A plant's **expressed set** is the subset of its 4 alleles actually
shown on stigma and pollen. Two cases:

- **Case A** — plant has ≥ 1 Class I allele → expressed set = its
  Class I alleles.
- **Case B** — plant has 0 Class I alleles → expressed set = all 4
  Class II alleles.

Same rule for mother and father.

#### A.7.3a How many identities does each plant display?

Because every plant carries 4 SRK allele *copies* but not all of them
are always distinct identities and not all distinct identities are
always expressed, it is useful to state clearly how many identities a
plant actually shows on its pollen and stigma. All 4 gene copies
produce SRK protein, but if any of them encode the same identity or
are Class-II-silent under Case A, the number of *distinct identities*
recognised at the stigma is smaller than 4:

| Genotype | 4-copy composition | Expressed identity set | # distinct identities |
|---|---|---|---|
| Class I homozygote | `{A, A, A, A}` (A Class I) | `{A}` | 1 |
| Class I duplex (AABB) | `{A, A, B, B}` (both Class I) | `{A, B}` | 2 |
| Class I diverse | `{A, B, C, D}` (all Class I) | `{A, B, C, D}` | 4 |
| Mixed with dominance | `{A, A, B, B}` (A Class I, B Class II) | `{A}` (Class II silent) | 1 |
| Class II homozygote | `{B, B, B, B}` (B Class II) | `{B}` | 1 |
| Class II heterozygote | `{B, B, C, C}` (both Class II) | `{B, C}` | 2 |
| Class II diverse | `{B, C, D, E}` (all Class II) | `{B, C, D, E}` | 4 |

**Consequence for compatibility.** SI recognition operates on the
*identity* level, not the copy level. A Class I homozygote AAAA
displays a single identity {A}; the fact that A occurs on 4 copies
means more SRK protein (dosage) but a single receptor phenotype at the
stigma. Under the § A.8 Case-A formula `P_compat = (1 − p(M))⁴` a
Class I homozygote has *higher* `P_compat` than a diverse heterozygote
of the same class, because `p(M)` — the local frequency mass of the
expressed set — is smaller when the set is smaller. Homozygosity is
therefore *protective* against mate limitation at the individual
scale, not a disadvantage.

**Empirical LEPA zygosity — this actually matters a lot.** The
Canu-amplicon Step 23 analysis of 367 preliminary genotyped
individuals with ≥ 1 functional SRK copy reveals the observed
distribution of distinct functional SRK identities per plant:

| # distinct functional identities | Genotype pattern | Count | Fraction |
|---|---|---|---|
| **1** (single identity — homozygous) | AAAA + AAA0 + AA00 | 241 | **65.7 %** |
| **2** (two identities) | AABB + AAAB + AAB0 | 118 | **32.2 %** |
| **3** (three identities) | AABC | 8 | **2.2 %** |
| 4 (four identities) | — | 0 | 0 % |

Source: [`srk_zygosity_empirical.tsv`](tables/Phase5/srk_zygosity_empirical.tsv).

**Two-thirds of LEPA plants carry only one distinct SRK identity** —
this is far higher than a naive independent-tetraploid-draw model
would predict (~5 % homozygosity under P1). The step30
finite-population simulation therefore draws each mother's genotype
from this empirical distribution: 65.7 % of simulated mothers are
homozygous, 32.2 % are 2-distinct, 2.2 % are 3-distinct. Candidate
fathers are drawn from the same distribution, giving a symmetric
mother/father empirical-zygosity model. The § A.8 species-mean
`P_compat` jumps from **0.16 (naive independent-draws model) →
0.69 (empirical LEPA zygosity, provisional Class I/II)** — a huge biological signal that
homozygosity is doing real protective work in this system.

#### A.7.4 The between-plant recognition rule

A cross (father → mother) is **compatible** if and only if:

**Mother's expressed set ∩ Father's expressed set = ∅.**

- **Class I × Class I** — compatible if and only if no shared expressed Class I
  allele.
- **Class I × Class II** — always compatible; the two expressed sets
  belong to disjoint classes by construction.
- **Class II × Class II** — compatible if and only if no shared expressed
  Class II allele.

The between-class rule is why sporophytic SI buffers reproduction
under drift: a location holding both classes still crosses freely
regardless of allele identity.

#### A.7.5 A worked example

Mother **M** carries `{FG001, FG002, FG024, FG031}` — two Class I
(FG001, FG002) + two Class II. Case A: her expressed set =
{FG001, FG002}; her Class II alleles are silent.

Three candidate fathers:

- **F1** = `{FG001, FG007, FG024, FG032}` — Case A, expressed
  {FG001}. Shares FG001 with M → **rejected**.
- **F2** = `{FG015, FG018, FG024, FG032}` — Case B, expressed all
  four Class II. No overlap with M's Class I set → **compatible**.
- **F3** = `{FG002, FG003, FG010, FG018}` — Case A, expressed
  {FG002, FG003}. Shares FG002 with M → **rejected**.

Under random mating, the fraction of compatible fathers is the § A.8
Case-A formula `P_compat = (1 − p(M))⁴`, where `p(M)` is the local
frequency mass of {FG001, FG002}.

### A.8 Finite-population prediction of pollen compatibility

**What we did.** Effective deme size is a property of each
**50 m connected component inside a location**, not of the location
as a whole: plants inside a component share pollen; plants in a
different component (same location but no pollen link) do not.
Phase A therefore simulates **each component independently** and
aggregates to the location level by component-size weighting:

1. For each 50 m component *c* with component_N_fertile *N_c*, draw
   4 × N_c alleles from the species-wide prior (tetraploid,
   see § A.3), compute the component's local Fg frequency vector,
   sample mother genotypes under the empirical LEPA zygosity
   distribution (§ A.6.3a) and evaluate each one's random-mating
   compatibility under the **sporophytic Class I / Class II
   dominance model** implemented in
   [`srk_si_model.py`](srk_si_model.py).
2. Component P_compat = mean over sampled mothers; posterior CI
   taken across simulation replicates.
3. **Location P_compat on replicate *k*** = size-weighted mean of
   its components, Σ_c (P_compat_{c,k} · N_c) / Σ_c N_c. Posterior
   CI across replicates.
4. Both tables are emitted —
   [`step30_A_prediction_location_pcompat.tsv`](tables/Phase5/step30_A_prediction_location_pcompat.tsv)
   (one row per location, the headline Phase A number) and
   [`step30_A_prediction_component_pcompat.tsv`](tables/Phase5/step30_A_prediction_component_pcompat.tsv)
   (one row per (location, component), so that fragmented
   locations whose headline mean is "sustainable" but which hide a
   small-component drift signal are not invisible).

The SI model itself captures three biological facts that the
diploid gametophytic approximation used in Part 1 could not:

- **Sporophytic SI.** Rejection is determined by the pollen parent's
  diploid (here, tetraploid) genotype — not by the individual pollen
  gamete's allele. Cross A × B is compatible if and only if parents A and B
  share no expressed allele.
- **Class I dominance within a plant.** A tetraploid carrying any
  Class I alleles expresses only those Class I alleles on pollen
  and stigma; a plant carrying only Class II alleles expresses all
  of them co-dominantly.
- **Between-class always compatible.** A Class I plant × Class II
  plant share no expressed alleles by construction, so their cross
  is always compatible. This is the mechanism behind the H1b
  hypothesis in the existing cross-plan
  ([`step26e_cross_plan_H1b_between_class_baseline.tsv`](tables/Phase5/step26e_cross_plan_H1b_between_class_baseline.tsv)).

**Class assignment.** The Fg → Class mapping lives in
[`tables/Phase5/srk_fg_class.tsv`](tables/Phase5/srk_fg_class.tsv)
as a first-class TSV, editable by hand. The **provisional data-driven
default** puts the 6 Fgs with ≥ 2 sub-alleles in Class I (FG001,
FG002, FG003, FG004, FG005, FG006 — together ~65 % of P1 frequency)
and the 26 single-allele Fgs in Class II (~35 %). FG024 is flagged
`REVIEW` because it is single-allele but 18 % of P1 carriers, and
its class may deserve manual reclassification once the SI literature
is checked. Rerunning `step30_srk_diversity_prediction_vs_observed.py`
regenerates every downstream number from an edited TSV.

**Analytical formulas.** Under a local Fg frequency vector *f* and
Class I total mass p_I = Σ_{j ∈ Class I} f_j, the sporophytic
pollen compatibility for a tetraploid mother with expressed set M and expressed
mass p(M) = Σ_{j ∈ M} f_j is closed-form:

- Case A — mother has ≥ 1 Class I allele (M ⊆ her Class I alleles):
  $P_{\text{compat}} = (1 - p(M))^4$
- Case B — mother has only Class II alleles (M = her 4 alleles):
  $P_{\text{compat}} = 1 - (1 - p_I)^4 + (1 - p_I - p(M))^4$

**Result under empirical LEPA zygosity (§ A.8.3a).** Species-mean
pollen compatibility is **~0.78** — very close to the Part-1 diploid
gametophytic estimate of 0.63, but reached via a completely different
mechanism. Under a naive independent-tetraploid-draw assumption the
species mean is only ~0.16, but that model over-counts heterozygosity
(it predicts ~5 % homozygotes; LEPA shows ~66 %). Once each mother
and father is drawn under the empirical zygosity distribution, ~66 %
of mothers are single-identity Case A homozygotes with a very small
expressed set → high per-mother pollen compatibility, ~32 % have
two identities → intermediate, ~2 % have three → lower. The
weighted mean lands at 0.78. Traffic-light bands are recalibrated
against this new species mean: **failed < 0.260, struggling
0.260–0.520, sustainable ≥ 0.520** (1/3 and 2/3 of species mean).
Every LEPA location's mean sits close to 0.78; **small-slickspot
uncertainty still shows up in wide credible intervals** — the BL5
tail (EO24 group, EO24-1, EO24-2, EO24-7) has 95 % CIs spanning
from "struggling" to nearly-1, honestly reflecting founder-effect
variance in the class + zygosity composition of their tiny local
pools.

**Take-home for reviewers.** The sporophytic model is more forgiving
than the Part-1 diploid gametophytic approximation. Two independent
mechanisms defend LEPA reproduction against drift: **(1)** the
between-class compatibility of Class I × Class II crosses, and
**(2)** the sheer numerical dominance of Class I (65 % of P1
frequency), which makes most mothers Class I heterozygotes with
non-zero pollen compatibility almost regardless of local drift. The mate-limited
signal Phase C's β₁ will pick up is therefore expected to be
**subtle** — driven by within-class allele skew rather than the
gross "compatibility floor" the diploid model implied.

**Traffic-light bands are recorded** in
[`step30_A_traffic_light_bands.tsv`](tables/Phase5/step30_A_traffic_light_bands.tsv)
so that any downstream analysis can join against them and reproduce
the categorical labels.

**Connection to Figure 3.** The per-location pollen compatibility
prediction is driven by the same drift mechanism the diversity
figure uses — the simulated local Fg pool — except that drift is now
modelled **per component** (one pool per 50 m connected component,
size `4 · component_N_fertile`) rather than through a single
"largest-component" proxy for the whole location (fragmentation → drift).
The location-level number in Figure 5 is the size-weighted mean
across those components. Seed counts do **not** enter the pollen
compatibility prediction. Figure 3 demonstrates that the existing
LEPA dataset recovers Nature's local Fg pool at every location
(coverage ≥ 90 % everywhere; ≥ 99 % at all but the two BL5
singletons). That validation transfers directly here: because the
sampling recovers the local Fg pool composition, the location-mean
pollen compatibility computed from the per-component simulation is
a faithful estimator of the population-mean pollen compatibility at
each location.

<a id="fig-5"></a>
![Figure 5](figures/Phase5/step30_A_prediction_fecundation.png)

**Figure 5.** Predicted per-mother pollen compatibility under the sporophytic tetraploid Class I / Class II model with **empirical LEPA zygosity** (§ A.8.3a). One dot per LEPA location, panelled by Bottleneck Lineage. Traffic-light background bands mark **failed** (pollen compatibility < 0.260, red), **struggling** (0.260–0.520, orange) and **sustainable** (≥ 0.520, green) — recalibrated against the sporophytic + empirical-zygosity species-mean of **0.780** (green dotted line). Dot position = **size-weighted mean of per-component pollen compatibility**: each 50 m connected component inside the location is simulated independently (4 × component_N_fertile alleles drawn from P1, mothers sampled under the empirical LEPA zygosity distribution — 66 % single-identity homozygotes, 32 % 2-distinct, 2 % 3-distinct, see [`srk_zygosity_empirical.tsv`](tables/Phase5/srk_zygosity_empirical.tsv) — compatibility evaluated against 300 candidate fathers drawn the same way), and the location-level number is the mean of its components weighted by component_N_fertile. Error bars = 95 % credible interval across simulation replicates; dot size ∝ √M (mothers with seed records in DB). The finer-grained per-component rows are preserved in [`step30_A_prediction_component_pcompat.tsv`](tables/Phase5/step30_A_prediction_component_pcompat.tsv) so that fragmented locations whose headline mean is "sustainable" do not hide a struggling sub-component. **Every LEPA location's mean sits close to the species mean because 66 % of mothers express only one SRK identity, which minimises their p(M) footprint under the § A.8.4 recognition rule; BL5 tiny slickspots retain wide CI reflecting founder-effect variance in class + zygosity composition.** Source: `step30_srk_diversity_prediction_vs_observed.py`. Class assignments: [`srk_fg_class.tsv`](tables/Phase5/srk_fg_class.tsv). Empirical zygosity: [`srk_zygosity_empirical.tsv`](tables/Phase5/srk_zygosity_empirical.tsv). Bands: [`step30_A_traffic_light_bands.tsv`](tables/Phase5/step30_A_traffic_light_bands.tsv).

### A.9 Cross-plot: SRK diversity vs pollen compatibility

**What we did.** We plotted per-location predicted SRK allele
diversity (x, from § A.6) against predicted sporophytic pollen
compatibility (y, from § A.8), with locations coloured by BL and the
sporophytic + empirical-zygosity species-mean pollen compatibility
(0.693) drawn as a horizontal reference line.

**Result.** Under the sporophytic Class I / II model with the
empirical LEPA zygosity distribution (§ A.8.3a), **most locations
sit near the species mean 0.78**. The BL5 tail (EO24 group) is
displaced left on the diversity axis (only ~2–5 alleles present)
but stays close to the species mean on the compatibility axis — the
combination of between-class buffering (§ A.8.4) and the ~66 %
homozygosity rate prevents even severely drift-collapsed populations
from falling into "failed" territory at the mean. What *does*
separate small locations from large ones is the **width of the 95 %
credible interval**: large slickspots have tight CIs anchored on the
species mean; small ones have wide CIs that can extend into
"struggling" territory, honestly reflecting the founder-effect risk
that a small pool draws an unlucky class + zygosity combination.

**Take-home for reviewers.** The sporophytic + empirical-zygosity
model reframes goal 3. The story is no longer "drift-collapsed sites
are predicted to fail outright" — because Class I / Class II
buffering plus 66 % single-identity homozygosity make outright failure
extremely hard. Instead, the risk is that a specific small slickspot
draws an unlucky combination (Class-II-only pool with one high-
frequency allele, or a rare 3-distinct heterozygote that ties up
its expressed footprint) and lands in "struggling" territory *within*
its own credible interval. Phase C's β₁ test picks this up as
localised within-class allele skew, which is a more subtle and more
biologically honest signal than the "global compatibility floor"
narrative that the diploid gametophytic model implied.

<a id="fig-6"></a>
![Figure 6](figures/Phase5/step30_A_diversity_vs_pcompat.png)

**Figure 6.** Predicted SRK allele diversity (x) vs predicted sporophytic pollen compatibility (y) per LEPA location, coloured by Bottleneck Lineage (Set1 palette: BL1 purple, BL2 blue, BL3 red, BL4 orange, BL5 green), under the sporophytic + empirical-zygosity model (§ A.8.3a). Error bars on both axes come from the same Dirichlet posterior draws that produced the two single-quantity Phase A figures. Dot size ∝ √M (mothers with seed records in DB). **Green dotted horizontal line = sporophytic + empirical-zygosity species-mean pollen compatibility (0.780)** — the reference every location can be read against under the § A.7 model. Most locations sit right on the species mean because ~66 % of LEPA plants express only one SRK identity, giving them a small p(M) footprint and therefore high pollen compatibility. **BL5 tail (EO24, EO24-1, EO24-2)** sits at low predicted diversity (x ≈ 2–5) but stays close to the species mean on the y-axis with much *wider* credible intervals — the honest founder-effect signal that a small slickspot can drop into "struggling" (y < 0.520) if class + zygosity composition breaks unfavourably at that specific site. This figure is the single-figure summary of goal 3's fragmentation × drift decomposition: **sporophytic + empirical zygosity flatten the mean, but small locations still carry drift-driven risk in their credible intervals**. Source: `step30_srk_diversity_prediction_vs_observed.py`.

---

## Part B — Sampling protocol derived from the predictions

Part A tells us *what SRK diversity each location should carry*. Part B
derives *what the field team must sample to test those predictions*.
The sampling design is not an arbitrary field protocol: it is engineered
to characterise, at each location, the SRK alleles that Part A's
finite-population model says are physically present, at a stated
statistical guarantee.

**Reading order.** § B is organised so the per-mother threshold —
the tetraploid per-mother cap of 15 seeds — is established first, then flows into
the mother allocation downstream:

- § B.1 – § B.2 set up the two population-size concepts and the
  50 m connectivity structure.
- § B.3 (per-mother seed count) establishes the tetraploid per-mother
  cap of 15 seeds/mother and the § B.3.1 power justification for
  Part C testing. The 15-seed answer defines the per-mother allele
  draw count A_delivered = 4 + 2·15 = 34, which the mother
  allocation inherits.
- § B.4 (per-location mother count) uses `N_fert_eff` from § A.5.1,
  the P1 prior, and the § B.3 per-mother draw count to allocate
  mothers per 50 m connected component plus a ≥ 1-per-event floor
  (fragmentation-aware allocation § B.4.2, locked recipe § B.4.3).
- § B.5 covers the two-year pooling rules.

### B.0 What one seed lot yields — the two-generation trick

Each sampled seed is **tetraploid** like its parents — 2 maternal + 2
paternal SRK alleles at each SRK locus — so one seed lot per mother
recovers, in a single extraction batch, two independent samples:

| Generation | What we recover | How |
|---|---|---|
| **Parental (G0)** | The mother's 4 SRK allele copies (her full genotype) | Invariant / 50 %-frequency alleles across her sibs |
| **Filial (G1) — as read through paternity** | The pollen SRK allele pool she was exposed to, 2 draws per seed | Variable alleles across her sibs |

Aggregating across the mothers of a location gives us **two independent
estimates of the location-level SRK allele frequency spectrum**: the
*maternal* one (who is standing there) and the *paternal* one (who is
actually contributing pollen). Under random mating with panmictic pollen
dispersal, the two spectra are indistinguishable. Departures flag biased
contribution, cryptic SI filtering, or immigrant pollen — themselves
useful signal. Every per-mother and per-location sampling target below is
engineered against this two-generation yield, so that one seed lot
delivers both halves at once.

### B.1 Census N_fertile vs permit-realistic M sampled

The pipeline carries **two distinct notions of population size** that must
not be conflated in downstream analysis:

- **`total_n_fertile`** — the field census: total number of fertile plants
  counted at a location (or an event). This determines the *pollen SRK
  allele pool* a mother is exposed to: `K = 4 × (N_fertile − 1)` under
  tetraploid LEPA (§ A.3). Pollen
  is contributed by every fertile plant whether or not the seed team was
  allowed to collect from it, so K correctly uses this census number.
- **`M_mothers_in_db`** — the permit-realistic reality: number of
  mothers for which the seed team was actually allowed to collect
  seeds, i.e. the count of `germplasmID` records tied to that location.
  For a large slickspot with 543 fertile plants under a 10 % permit,
  `M_mothers_in_db` is around 55, not 543. This is the number that
  drives the Phase A predictions above and every Phase B output.
- **`M_achievable_ceiling`** — the theoretical design ceiling assuming
  the permit allowed sampling every fertile plant. Kept in the output
  tables as a reference upper bound; **not** used by the prediction.

The prior version of the pipeline used `M_achievable_ceiling` in the
prediction, which over-estimated detectable SRK diversity at large
slickspots (e.g. EO61 predicted 29 distinct alleles at `M = 543`, versus
21 at the permit-realistic `M = 55`). The current pipeline uses
`M_mothers_in_db` throughout.

### B.2 Within-location pollen connectivity — the biological floor

The LEPA "location" is a curated grouping of nearby slickspot events;
it is **not** by itself a mating unit. Whether the events at a location
actually exchange pollen depends on their spatial arrangement relative
to the pollen-flight radius. This is a foundational Phase A input
because **connectivity is a direct predictor of both realised SRK
diversity and random-mating compatibility** at the location scale:

- Fewer connected adults → smaller effective deme → more drift →
  fewer distinct SRK alleles.
- Fewer connected adults → higher per-mother probability that the pollen
  she sees carries her own alleles → lower random-mating compatibility.

**How it is computed** ([`step29b_location_connectivity.py`](step29b_location_connectivity.py)).
For each location, we build a graph on its events using haversine
distance and draw an edge between events whose coordinates lie within
R metres. Connected components are found by breadth-first search.
Per-location outputs at R = 10 m, 25 m, and **50 m (primary)**:

- **`connected_share_{R}m`** — fraction of adults in a component of
  more than one event (i.e. that exchange pollen with at least one other
  event under the R-radius assumption).
- **`largest_component_share_{R}m`** — fraction of adults in the
  location's largest connected component.
- **`n_components_{R}m`** — number of disconnected mating units the
  location is broken into.

**How it feeds Phase A prediction.** The per-location `largest_component_share_50m`
is used as an **effective-N multiplier**: the finite-population
compatibility model in Step 30 draws the local pollen pool from
`total_n_fertile × largest_component_share_50m` alleles instead of
`total_n_fertile`, so drift acts on the mating unit that actually
exchanges pollen. The mate-limitation regression in Phase B adds
connectivity as an explicit third predictor alongside pollen
compatibility and mating-neighbourhood size, letting the model
attribute reduced seed set to allele-frequency drift vs pollen-flow
fragmentation as distinct causal channels.

**Current LEPA reading** (2025 field only, 39 locations):

| Pollen-flight radius | Locations ≥ 90 % adults connected | Locations < 50 % adults connected |
|---|---|---|
| 10 m (small-bee patch — conservative) | 0 / 39 | 35 / 39 |
| 25 m (short-flight — sensitivity) | 8 / 39 | 21 / 39 |
| **50 m (extended foraging — primary)** | 18 / 39 | 11 / 39 |
| 75 m (long-flight — sensitivity) | 39 / 39 | 0 / 39 |
| 200 m (upper-bound sensitivity) | 39 / 39 | 0 / 39 |

At the primary radius (50 m), **18 / 39 locations** already behave as
a single mating unit and **11 / 39** sit below 50 % connectivity — a
meaningful gradient across BLs that keeps the fragmentation channel
identifiable in Phase C's β₂ regression while matching the halictid /
small solitary bee foraging literature.

**Why 50 m as the primary?** A pollinator-radius sensitivity sweep
across 10, 25, 50, 75, 100, 150, 200 m
([`step30_A_radius_sensitivity_summary.tsv`](tables/Phase5/step30_A_radius_sensitivity_summary.tsv))
shows two things:

- **Connectivity plateaus at 75 m** — median location-scale
  connectivity hits 1.00 (every LEPA location fully connects into one
  pollen-flow component). Anything above 75 m brings no new adults
  into the deme.
- **§ B.4.2 sampling cost stabilises at 75–100 m** — total mothers
  required drops steeply from 887 at 10 m to 505 at 50 m (43 % drop)
  and 416 at 75 m (another 18 % drop), then flattens.

50 m is where fragmentation still has *variation across BLs* (median
connectivity 0.81, not yet saturated at 1.00), sampling effort has
already dropped **~43 %** vs 10 m, and pollen-compatibility
predictions barely change with radius (empirical zygosity dominates).
50 m sits comfortably within the published halictid + small-bee
foraging range (10 – 100 m). See the § A.5 pollinator-radius
sensitivity figure below for the full sweep.

**Take-home for reviewers.** Fragmentation is not an abstract threat —
it is measurable at the within-location scale, before any genotyping.
This is the biological reason every Part A prediction has to run a
finite-population model, not a species-wide random-mating limit.
Connectivity enters the pipeline as an **effective-N multiplier**
for the compatibility model and as an **explicit predictor** in the
Phase B mate-limitation regression. The per-location connectivity
share is visualised directly as Panel C of [Figure 2](#fig-2)
(`N_fert_eff` in § A.5.1); the full radius sensitivity that
justifies 50 m as the primary choice is in [Figure 1](#fig-1)
(§ A.5).

### B.3 Per-mother seed count — sampling design (Step 28)

**Rationale — why the protocol must scale.** An LEPA slickspot with two
flowering plants is a fundamentally different sampling target than one
with a hundred. In the two-plant case, every seed a mother makes must
have been sired by the single other plant, so a handful of seeds
already tells you the whole story. In the hundred-plant case, the pool
of potential pollen SRK alleles is up to two orders of magnitude larger,
and the same sampling effort characterises only a small fraction of it.

Any uniform "genotype *n* seeds per mother" protocol wastes effort in
one regime and underdelivers in the other. Worse, it fails to be
defensible in the regime that matters most for conservation: **small,
isolated slickspots where mate limitation is the exact mechanism we
are trying to detect**. In those spots, a coarse uniform protocol would
either over-sample (irrelevant) or under-sample (missing the very
mothers whose fecundation is at risk). The design must therefore
**scale with the mate-availability context of each mother**.

**Methodology.** We treat mating as a **allele-detection problem on
SRK alleles**. Under a random-mating null with equal pollen contribution
from every compatible mate, each seed's paternal SRK allele is one
independent multinomial draw from a pool of

$$K \;=\; \text{PLOIDY} \times N_{\text{compatible}} \;=\; 4 \times N_{\text{compatible}}$$

potential pollen-donor SRK allele copies (**4 alleles per tetraploid
mate**, see § A.3). We do not yet apply a self-incompatibility filter
because the SRK genotypes of every plant in every event are not yet
known; when they are, `N_compatible` will drop further and the sample
sizes will drop with it.

> **One K, two scopes.** The letter *K* is used only for **the size of
> the pollen-donor SRK-allele pool a mother is exposed to** — a single
> biological quantity. What changes is the *spatial scope* of that pool:
>
> - **K (event-only)** — pollen alleles from her event alone,
>   `4 × (N_fertile − 1)` under tetraploid LEPA. Used in the
>   allele-detection formulas below.
> - **K^(R) (spatial neighbourhood at radius R metres)** — same
>   quantity, extended to include all events within *R* m of hers. Used
>   later as the fragmentation predictor in the Part C mate-limitation
>   regression (a small K^(R) means the mother sees few pollen donors
>   because she is physically isolated, not because her local Fg pool
>   is skewed).
>
> Column names in the TSVs follow the same convention: `K_event`,
> `K_spatial_10m`, `K_spatial_50m`, `K_spatial_50m`. A third symbol —
> `K_fg = 32` — appears only when referring to the **number of
> species-wide SRK allele classes** (functional groups) in the P1
> empirical prior; it is a count of allele identities, not a count of
> alleles in a pool.

**Spatial mating-neighbourhood extension.** The event-only pool
underestimates K when a mother's event has other LEPA events physically
close by — a small bee will forage across such events indiscriminately.
Using `eventDecimalLatitude` / `eventDecimalLongitude` from the `Events`
table (haversine distances in metres), we compute for each event and
radius R the number of neighbouring events within R m and sum their
`N_fertile`:

$$N_{\text{compatible}}^{(R)} \;=\; \Bigl(\sum_{e \in \text{neighbours}_{\leq R}} N_{\text{fertile}}(e)\Bigr) - 1$$

$$K^{(R)} \;=\; 4 \times N_{\text{compatible}}^{(R)}$$

R = **50 m** is the primary assumption. Sensitivity radii **10 m** and
**50 m** are also computed. In the current LEPA data 235 / 704 events
(33 %) have a neighbour within 10 m, 439 / 704 (62 %) within 50 m,
579 / 704 (82 %) within 50 m — so the spatial extension is a real
effect, not a formality.

**Two decision rules on the per-mother seed count.** Both come from the
same multinomial identity: an allele of frequency *p* is missed after
*n* independent seed draws with probability $(1-p)^n$.

**Rule 1 — expected coverage (event-size-dependent).** Choose the
smallest *n* such that the expected fraction of the K potential alleles
observed exceeds a target *c* = 0.90:

$$n_{\text{exp}} \;=\; \left\lceil \frac{\log(1-c)}{\log(1 - 1/K)} \right\rceil$$

This is the target when it is *reachable* — typically at small and mid-sized
events.

**Rule 2 — miss-probability guarantee (event-size-independent).** Choose
the smallest *n* seeds such that any allele contributing at least
$p_{\min} = 0.10$ of siring events is detected with probability at
least $1 - \alpha = 0.95$. In terms of the required number of
independent paternal *allele draws*,

$$n_{\text{draws}} \;=\; \left\lceil \frac{\log \alpha}{\log(1 - p_{\min})} \right\rceil \;=\; 29$$

Under **tetraploid LEPA** each seed carries 2 paternal alleles, so the
required seed count is

$$n_{\text{seeds}} \;=\; \left\lceil \frac{n_{\text{draws}}}{\text{paternal alleles per seed}} \right\rceil \;=\; \left\lceil \frac{29}{2} \right\rceil \;=\; 15$$

This is a **K-independent floor** — the honest per-mother cap when
Rule 1 becomes impractical at very large events. (Under a diploid
model — 1 paternal allele per seed — Rule 2 would give 29 seeds. The
tetraploid cap is exactly half.) Both rules are wrapped with a
**simulation check** (multinomial resample of pollen alleles) so the
analytical curve is bracketed by an empirical 95 % CI.

**Accounting for uneven seed production.** Mothers do not all make the
same number of seeds. The recommended *n* above is a *target*; the
achievable *n* for a given mother is capped by her real seed budget *S*
(`germplasmQuantityEstimate`):

$$n_{\text{use}} \;=\; \min(n_{\text{recommended}},\ S)$$

We therefore report, alongside the recommended target, **the coverage
each mother actually delivers with her real budget**:

$$\text{achieved coverage} \;=\; 1 - \left(1 - \tfrac{1}{K}\right)^{n_{\text{use}}}$$

**Aggregation across mothers rescues the target.** The 15-seed cap
does *not* mean we give up on the 90 %-coverage story. It is a
per-mother cap; when we pool across the M mothers at a location, each
mother contributes **34 allele draws** to the location's total pool
under tetraploid — 4 alleles from her own genotype + 2 paternal
alleles per seed × 15 seeds = 30. Under the P1 species-wide prior the
90 %-of-32-species-alleles target requires **A ≈ 704 allele draws
pooled at the location** (§ B.4), so `M ≥ 704 / 34 ≈ 21` mothers is
enough to reach it. This is why a large slickspot like EO76 (62
mothers × 34 = 2 108 allele draws) far exceeds the target even though
each individual mother's 15-seed sample only reveals ~19 of her 32
local pollen alleles.

**Script.** [`step28_seed_sampling_per_mother.py`](step28_seed_sampling_per_mother.py). Outputs the
per-mother recipe TSV plus the two-panel coverage figure below.

<a id="fig-9"></a>
![Figure 9](figures/Phase5/step28_coverage_curves.png)

**Figure 9.** two-panel SRK allele detection under **tetraploid LEPA** (4 SRK copies per plant; 2 paternal alleles per seed) on the absolute-allele scale, with the local pool capped at the species-wide ceiling of 32 Fgs. Shaded bands = 95 % simulation CI. **Panel A — per mother.** x = seeds genotyped, one curve per event-size bin (each seed contributes 2 paternal allele draws); the vertical red line at **15** marks the per-mother operational cap under tetraploid (never ask any single mother for more than 15 seeds). The 15-seed dots show what each event size **delivers per mother**: **4.0 of 4** (1–2 plants), **11.1 of 12** (3–5), **18.6 of 28** (6–10), **19.7 of 32** (11–20), **19.7 of 32** (21–50), **19.7 of 32** (>50). One mother's 15 seeds cannot saturate a 32-allele pool — this is the allele-detection limit for a single sampler, not undersampling. **Panel B — aggregation across mothers at a location.** x = number of mothers sampled at the location (15 seeds each = **34 allele draws per mother: 4 maternal + 2·15 = 30 paternal**). At the 5-mother benchmark (green dotted line, 75 cumulative seeds): **4.0 of 4**, **12.0 of 12**, **27.9 of 28**, **31.9 of 32**, **31.9 of 32**, **31.9 of 32** — every event size reaches its local ceiling. Because Panel B is a allele-detection simulation continuous in M, any real location can read off its own coverage by locating (its event-size bin, its actual M) on the correct curve. The "gap at large events" in Panel A closes cleanly at the location scale — see § B.4. Source: `step28_seed_sampling_per_mother.py`. Data: [`step28_coverage_curves_by_Nfertile.tsv`](tables/Phase5/step28_coverage_curves_by_Nfertile.tsv), [`step28_aggregation_curves_by_Nfertile.tsv`](tables/Phase5/step28_aggregation_curves_by_Nfertile.tsv).

#### B.3.1 Does 15 seeds give enough Part C testing power?

Figure 9 is a allele-detection argument: 15 seeds is the tetraploid
per-mother floor for detecting every SRK allele at a mother's local
pool with 90 % probability. But **allele detection is not the same
as regression power** — the § C.1 mate-limitation test needs
per-mother observed P_compat precise enough that β̂₁ is
identifiable. Two Monte-Carlo simulations answer that separately.

**Figure 9a — per-mother P_compat precision** ([`step28d_pcompat_precision.png`](figures/Phase5/step28d_pcompat_precision.png)).
For a mother with true P_compat p, each of her seeds is an independent
Bernoulli(p) draw of paternal compatibility. Observed P_compat has
binomial SE `√(p·(1−p) / n_seeds)`. **Result:** SE drops steeply from
~0.28 at 3 seeds to ~0.12 at **15 seeds** (visible plateau), then
improves only marginally (~0.06 at 50 seeds). Precision-per-seed
returns diminish sharply above 15.

**Figure 9b — § C.1 mate-limitation regression power** ([`step28d_matelim_power.png`](figures/Phase5/step28d_matelim_power.png)).
Full-pipeline simulation of 505 mothers across 39 locations under the
model `seeds ~ intercept + β₁ · P_compat_obs + ε`, with `P_compat_obs`
noised by the seed-count-driven binomial variance. **Result at
n_seeds = 15**:

- **β₁ = 100 seeds/unit** (medium effect: half the ovule set difference
  between P_compat 0.5 and 1.0 mothers) → **99.6 % power**. Well above
  the 80 % target.
- **β₁ = 150 or 200** → 100 % power.
- **β₁ = 50** (small effect) → 69 % power. Below 80 %; would need
  ~30 seeds to reach 80 % at that effect size.

**Errors-in-variables attenuation.** Because `P_compat_obs` is a noisy
predictor, β̂₁ is biased toward zero (attenuation ≈ 0.42 at
n_seeds = 15 — β̂₁ recovers ~42 % of true β₁). This affects effect-
size *estimation* but not the ability to *detect* β₁ ≠ 0 at
realistic effect sizes. Where reporting β̂₁ as an effect size, we
apply a standard errors-in-variables correction.

**Take-home.** 15 seeds/mother clears both the allele-detection
threshold (Figure 9) and the regression-power threshold for medium-
to-large β₁ (Figure 9b). Small-effect regression power is marginal —
locations under-sampled by the fragmentation-aware allocation (the 25
mothers short at 6 locations in § B.4.3) are the most likely place
that shortfall shows up. The 2026 top-up sampling addresses both.

#### B.3.2 Field → lab correction — 60 % germination rate

Part C genotypes **seedlings**, not seeds (design confirmed
2026-10-03). LEPA seeds germinate at **≈ 60 %** under greenhouse
conditions, so the per-mother target of 15 genotyped **seedlings** per
mother translates to a field collection target of
`ceil(15 ÷ 0.60) = 25` seeds per mother. The 90 % allele-detection
floor in Figure 9 is a floor on the number of seedling genotypes
delivered to the lab, not on the number of seeds extracted in the
field.

Totals under the fragmentation-aware allocation (§ B.4.2):

- **Collect / germinate:** 505 mothers × 25 seeds = **12 625 seeds**.
- **Expected after germination (at 60 %):** ~15 seedlings/mother × 505
  = **~7 575 seedlings**.
- **Genotype:** up to 15 seedlings per mother → **~7 575 seedling genotypes**.

The germination rate is held in a single constant
`GERMINATION_RATE = 0.60` in
`step29c_fragmentation_aware_sampling.py`. The per-mother sampling
TSV `step29c_partC_germplasmID_selection.tsv` now carries both
`n_seeds_to_germinate = min(25, seeds_available)` and
`n_seedlings_to_genotype = min(15, round(n_seeds_to_germinate · 0.60))`
so the field protocol and the lab target are both auditable from one
file.

### B.4 Per-location mother count — sampling design (Step 29)

> **Note.** § B.4 and § B.4.1 below are the *reference / audit* view
> of the per-location allocation. **The operational Phase 5 allocation
> is § B.4.2** (fragmentation-aware, adopted as default). Read this
> section for the method underlying the location-scale allele-detection
> and 704-allele-draw target; use the § B.4.2 per-event `M_frag` column
> for actual field-team numbers.

**Rationale.** Step 28 tells us *how many seeds per mother*. It does
not tell us *how many mothers to sample per event* or *how many events
to sample per location*. Without those two extra layers, per-mother
numbers cannot be aggregated into a location-level estimate with a
stated coverage guarantee.

**Method — uniform-prior version.** We apply the same allele-detection
logic one level up: instead of sampling pollen SRK alleles, we sample
**maternal SRK alleles** across mothers and events. Under tetraploid
LEPA the pool has size K = 4 × N_fertile at the event level and
K_loc = 4 × ΣN_fertile at the location level. Under a
uniform-per-plant allele draw, the expected coverage after M mothers
sampled (each contributing PLOIDY = 4 maternal allele draws = 4M
draws in total) is $1 - (1 - 1/K)^{4M}$. Mother slots are allocated
to events **proportionally to N_fertile** (larger events get more
sampled mothers), capped at each event's actual N_fertile.

**Method — empirical-prior version (P1 species-wide).** The uniform
pool is a simplification: LEPA has a very skewed Fg frequency
distribution (FG001 alone = 40.7 %). Under a skewed prior, common Fgs
saturate quickly but rare Fgs need much more sampling. Using the same
species-wide P1 prior as Part A, the expected fraction of the 32 Fgs
observed after A total allele draws at the location is:

$$\text{expected Fg coverage} \;=\; \frac{1}{K_{\text{fg}}}\sum_{j=1}^{K_{\text{fg}}}\left[1 - (1-f_j)^{A}\right]$$

The 90 %-of-Fgs target is the smallest A for which this quantity
reaches 0.90. Under the current P1 prior this target is **A = 704
total allele draws / location**. Each sampled mother contributes both
her own genotype and the paternal alleles of any seeds we genotype:

$$A_{\text{delivered}} \;=\; M \times (\text{ploidy} + \text{paternal alleles per seed} \times n_{\text{seeds}}) \;=\; M \times (4 + 2 n_{\text{seeds}})$$

under tetraploid LEPA. The target A can be reached either by more
mothers, more seeds per mother, or both. When M is fixed by field
logistics, the seeds-per-mother needed is

$$n_{\text{seeds}} \;=\; \left\lceil \frac{A_{\text{target}}/M - 4}{2} \right\rceil$$

**Two levels of the same coverage story, made explicit:**

| View | Question answered | Metric | What 15 seeds/mother buys (tetraploid per-mother cap) |
|---|---|---|---|
| **Step 28 (per-mother, [Figure 9](#fig-9))** | How much of ONE mother's local pollen pool do her 15 seeds reveal? | Distinct pollen-donor alleles detected out of K local | 100 % at ≤ 5-plant events → ~19 of 32 at > 20-plant events |
| **Step 29 (per-location, this section)** | How much of the SPECIES-WIDE 32-SRK-allele pool does the location as a whole reveal after pooling across mothers? | Fraction of the 32 Fgs observed at the location | ≥ 90 % once M × 34 ≥ 704 — i.e. from ~21 sampled mothers upward |

**Result — sampling protocol scaled to each mother's mate context.**
Every mother plant received an individually calibrated seed-genotyping
recommendation based on the number of fertile plants at her slickspot
and the number of seeds she produced. Two decision rules run
side-by-side (the aspirational expected-coverage target and the
tetraploid per-mother operational cap at 15 seeds; both derived in
§ B.3.1). Under the P1 species-wide prior and the tetraploid
`M × 34 draws / mother` accounting, **46 / 52 locations (88 %) reach
the 90 %-of-32-Fgs target with their permit-realistic mother count**
(from Step 29's summary). The 6 that fall short are tiny slickspots
(≤ 10 fertile plants) where no sampling intensity can compensate for
the small census — a mathematical certainty, not a design failure.

Field sampling effort is auditable at the per-mother level. The field
team receives one number per germplasmID
([`step29_field_team_sampling_recipe.tsv`](tables/Phase5/step29_field_team_sampling_recipe.tsv)).

#### B.4.1 How many mothers per location to observe every predicted SRK allele

**Question.** Panel B of [Figure 9](#fig-9) shows *expected*
coverage as we aggregate mothers. But biologists asking "have I seen
every allele my location carries?" need a **probability**, not an
expectation — the last allele is the hardest to catch. And preliminary
Canu-amplicon data show **private alleles** across events (drift
signature), so any event skipped at sampling time carries a risk of
missing its private alleles entirely.

**Two floors combined.** For each location we compute the smallest
number of mothers `M_recommended` that meets **both**:

- **Per-mother allele-detection floor** (`M_uniform_full_detection`): under the
  uniform-frequency assumption, the smallest M such that the
  probability of observing **every one** of the location's predicted
  SRK alleles is at least **90 %**. Uses M × 34 allele draws (**4
  maternal + 30 paternal** per mother under the tetraploid per-mother cap of 15
  seeds) against the local pool K_local = min(4 · (total fertile
  plants − 1), 32 species alleles).
- **Private-allele floor** (`M_event_coverage`): at least one mother
  per event at the location. This is a *drift-aware* constraint —
  private alleles cannot be caught in an event that is not visited.

The recommendation is `M_recommended = max(per-mother floor,
private-allele floor)`, and the field team is asked to **spread the
mothers across events** so every event contributes ≥ 1.

**Result — per event-size bin under the uniform model (tetraploid):**

| Event size | K local (4·(N−1)) | M for 90 % chance to see all K |
|---|---|---|
| 1–2 fertile plants | 4 | 1 mother |
| 3–5 | 12 | 2 mothers |
| 6–10 | 28 | 5 mothers |
| 11–20 | 32 | 6 mothers |
| 21–50 | 32 | 6 mothers |
| >50 | 32 | 6 mothers |

**Result — per location combining both floors.** Median
`M_recommended = 7` mothers per location; range 1–82. **The
private-allele floor binds at 31 / 52 locations** — for the majority
of LEPA locations the real sampling constraint is not the
allele-detection maths, it is the number of events that must each be
represented. Numeric per-location results in
[`step28_mothers_for_full_detection_by_location.tsv`](tables/Phase5/step28_mothers_for_full_detection_by_location.tsv);
the visual view has been superseded by [Figure 12](#fig-12) (§ B.4.2)
which shows the honest fragmentation-aware allocation against what
the LEPA DB already holds.

**Caveat.** The uniform-frequency assumption is optimistic; when
allele frequencies are skewed (which is what drift produces), rare
alleles need substantially more mothers than the 90 % bound suggests.
The private-allele floor is the practical safeguard — it forces at
least one mother per event so no event's private alleles are missed
even when frequency skew is severe.

#### B.4.2 Fragmentation-aware sampling — comparison with the current allocation

**Why the § B.4.1 allocation under-counts fragmentation.** § B.4.1
computes the M for 90 % allele detection using the location's **aggregate**
`total_N_fertile` as a single pool. But two events sitting 300 m apart
at the same locationID are not one mating unit — no pollen crosses
between them. Treating them as one pool implicitly assumes a mother
sampled at event A pays for detection at event B, which is only true
when A and B are within pollinator range. The per-location mating-
pool decomposition in [Figure 1b](#fig-1b) (§ A.5) makes this vivid:
at locations like EO27-1 (8 events across multiple disjoint 50 m
mating units), the aggregate pool is a mathematical fiction.

**Fragmentation-aware allocation — the method.** For each location we
build the 50 m adjacency graph on its events (same graph § B.2 uses
for connectivity), find the connected components, and apply the
90 % allele-detection rule **per component** rather
than per whole location:

1. **Per component c** — total N_fertile in c → local pool
   K_c = min(4·(N_c − 1), 32 species alleles) → mothers needed
   M_c = mothers_for_full_detection(K_c, 0.90). Isolated singleton
   components → M_c = 1.
2. **Within a component**, distribute M_c across events proportional
   to N_fertile(e), with a **maternal-genotype floor of ≥ 1 mother
   per event** so every event gets visited.
3. **Location total**: M_frag_aware = Σ_c mothers allocated to c.

**Result — how does the current LEPA DB match the target?**
Under the fragmentation-aware allocation, the design asks for
**505 mothers across the 39 LEPA locations**. The LEPA DB already
holds **765 collected mothers**. Comparing target to what's already
in hand at the location scale ([Figure 12](#fig-12)):

- **33 / 39 locations are fully covered by the LEPA DB** — the
  mothers already collected meet or exceed the § B.4.2 target.
- **6 / 39 locations are short of the target** — the DB does not
  hold enough mothers to satisfy `M_frag_aware`. Those six are
  EO27-1 (short by 10), EO27RT (5), EO26-2 (4), EO8 (3), EO67 (2),
  and EO8 (1). Total shortage: **25 mothers**.
- The 25-mother shortage is the field-team top-up target for the
  2026 season — every one of the six locations is a candidate for
  focused collection.

**Per-component granularity.** The location totals hide within-
location concentration effects: a location can have enough mothers in
bulk but they may cluster at well-covered events, leaving parts of
the location under-provisioned. The Part C selector reallocates
within 50 m connected components because events inside a component
share a pollen pool (that is the biological definition of 50 m
connectivity), so coverage travels freely within a component. At the
component scale, **35 / 103 components are short by a combined
76 mothers** (see the Part C lab recipe below). The component-scale
gap is the finer-grained diagnostic; the 25-mother location-scale
gap is what the field team plans against.

**Interpretation.** The B.4.1 allocation is **optimistic** because it
assumes pollen mixes across the whole location. Fragmentation-aware
is the **honest** cost of the same 90 % guarantee under the actual
50 m mating structure. The delta is the number of mothers B.4.1 was
implicitly assuming would come "for free" through pollen sharing that
doesn't actually happen.

**Decision — we adopt fragmentation-aware allocation as the default
for Phase 5.** The 1.9× effort increase is a real permit implication,
but the alternative — publishing a 90 % coverage claim that only
holds under an unrealistic single-pool assumption — is worse for
the fragmentation × drift decomposition (Goal 2). If our sampling
protocol builds the fragmentation assumption *in* by under-allocating
at fragmented locations, we cannot then use those locations to *test*
fragmentation as a mate-limitation channel; we would be arguing from a
floor we ourselves lowered. § B.4.2 is the allocation that Phase C's
β₂ regression is entitled to inherit.

**What supersedes what.** § B.4.1's per-location `M_recommended`
column and § B.4's proportional-to-N_fertile allocation are now
**reference / audit views**, not the operational recipe. The
operational per-event allocation is
[`step29c_sampling_frag_aware_per_event.tsv`](tables/Phase5/step29c_sampling_frag_aware_per_event.tsv),
column `M_frag` — one row per event, one number per event, which
already satisfies both the per-component 90 %
guarantee and the ≥ 1 mother-per-event maternal-genotype floor.

**How this hooks into Phase C.** Phase C's mate-limitation regression
(§ C.1) inherits the new allocation automatically: the β₂ predictor is
per-mother K^(50m), which is a property of the landscape and doesn't
depend on which allocation delivered her seeds. What *does* change is
the **statistical power** — the fragmentation-aware allocation adds
mothers exactly at the locations where within-location K^(50m)
variance is highest, which is where β₂ has the most identifiability.

<a id="fig-12"></a>
![Figure 12](figures/Phase5/step29c_sampling_comparison.png)

**Figure 12.** Do we have the mothers the fragmentation-aware design asks for? Per LEPA location, panelled by Bottleneck Lineage in BL_ORDER. **Solid bar** = mothers needed under the § B.4.2 target (`M_frag_aware`). **Light bar** = mothers already collected and stored in the LEPA DB. Green "covered" annotation = DB has ≥ target; red "short by N" = DB has fewer than the target and would need a field top-up. Row labels give locationCode, number of events, number of 50 m components, census, and effective deme. **33 / 39 locations are fully covered by the current DB. 6 / 39 locations are short (25 mothers total — the 2026 field top-up).** The six short locations are EO27-1 (short by 10), EO27RT (5), EO26-2 (4), EO8 (3), EO67 (2), and EO8 (1). Source: `step29c_fragmentation_aware_sampling.py`. Data: [`step29c_sampling_comparison_per_location.tsv`](tables/Phase5/step29c_sampling_comparison_per_location.tsv), [`step29c_sampling_frag_aware_per_event.tsv`](tables/Phase5/step29c_sampling_frag_aware_per_event.tsv).

#### B.4.3 The locked field-team recipe

**Design decisions locked-in for the next field season:**

| Quantity | Value | Source |
|---|---|---|
| Primary pollinator radius | **50 m** | § A.5, [Figure 1](#fig-1), step29a sensitivity sweep |
| Drift unit | **50 m deme** = 50 m connected component of events; `component_N_fertile` adults per deme | § A.5, [Figure 1b](#fig-1b), step29c + step29d |
| Per-mother seedling-genotype floor | **15 seedlings/mother** | § B.3, step28 |
| Field-side germination assumption | **60 %** | § B.3.2, step29c `GERMINATION_RATE` |
| Seeds to germinate per mother | **25** (= ceil(15 ÷ 0.60)) | § B.3.2 |
| Mother allocation per location | Fragmentation-aware `M_frag` (per 50 m deme + ≥ 1 per event) | § B.4.2, step29c |
| **Total mothers across 39 locations** | **505** | step29c |
| **Total seeds to germinate** | **505 × 25 = 12 625 seeds** | Locked |
| **Total seedlings to genotype** | **~505 × 15 = ~7 575 seedlings** | Locked |
| Allele-detection target | 90 % probability of detecting every allele in each 50 m deme | step28, step29c |

**Rationale — why these are the minimum.** Under the tetraploid per-mother cap
(15 seeds/mother), each mother delivers 4 maternal + 30 paternal =
34 allele draws. At the biggest 50 m component (EO76's largest cluster
with ~150 fertile plants), the per-mother floor for 90 %
detection of a 32-Fg pool is roughly 6 mothers × 34 draws each — well
below 15 seeds per mother. Reducing seeds per mother below 15 would
require adding more mothers to compensate, which is a worse trade
because the paternal contribution per seed doubles the allele draws
much more efficiently than adding a fresh mother's 4-copy genotype.

**Rationale — why fragmentation-aware.** The pivotal `N_fert_eff`
metric (§ A.5.1) tells us that many locations have a raw census far
larger than their drift-relevant deme. Sampling proportional
to raw census would over-sample well-connected locations and
under-sample fragmented ones. Sampling per 50 m component (this
recipe) puts effort where mate limitation is actually expected —
the locations whose `N_fert_eff` is a small fraction of their raw
census.

**Three authoritative files** — one for the field team, one for the
lab, and one for querying the event → component mapping.

1. [`Tables/Phase5/step29c_sampling_frag_aware_per_event.tsv`](tables/Phase5/step29c_sampling_frag_aware_per_event.tsv)
   — the **field-team recipe**. One row per event, column
   `M_frag` = number of mothers to sample at that event. Used by the
   field team when planning a new season's collection.
2. [`Tables/Phase5/step29c_partC_germplasmID_selection.tsv`](tables/Phase5/step29c_partC_germplasmID_selection.tsv)
   — the **lab recipe for Part C**. One row per selected `germplasmID`
   already in the LEPA DB, drawn from the field-team recipe using the
   **50 m deme target** (not per-event). Sorted
   `EOID → locationID → component_id_within_loc → germplasmID`
   (adopted 2026-10-03). Selection rule per deme: (i) enforce
   the ≥ 1-mother-per-event maternal-genotype floor at every event
   that has germplasm in the DB; (ii) fill the remaining
   `component_M_target − floor_count` slots from the deme's
   leftover mothers sorted by `seeds_available` descending.
   Deterministic tie-break by `germplasmID`.
   **Seedling-not-seed targets** — each row carries both
   `n_seeds_to_germinate = min(25, seeds_available)` and
   `n_seedlings_to_genotype = min(15, round(n_seeds_to_germinate · 0.60))`,
   so the field protocol and the lab target are one join away.
   Columns: `EOID`, `locationCode`, `locationID`,
   `component_id_within_loc`, `germplasmID`, `occurrenceID`,
   `eventID`, `seeds_available`, `seedlings_target`,
   `n_seeds_to_germinate`, `n_seedlings_expected`,
   `n_seedlings_to_genotype`, `n_seeds_to_genotype` (deprecated alias),
   `component_M_target`, `component_n_available_in_DB`,
   `component_gap`, `selection_reason` (`floor` or `top_seeds`),
   `selection_priority_within_component`, `event_M_frag`.
3. [`Tables/Phase5/step29c_event_to_component_50m.tsv`](tables/Phase5/step29c_event_to_component_50m.tsv)
   — the **event → 50 m component lookup**. One row per event,
   columns `locationID`, `locationCode`, `eventID`, `component_id_50m`,
   `event_n_fertile`, `component_N_fertile`, `component_N_events`,
   `component_K`, `component_M_target`. The single canonical source
   for the question "which events share a pollen pool?" — used both
   by the Part C selector and available for any downstream analysis.

**What the lab recipe delivers today.** From the 765 germplasmIDs
already in the LEPA DB, the component-aware selection algorithm picks
**431 mothers × ≤ 15 seeds each = 6 459 seeds** for Part C genotyping.
Of the 505 M_frag target, the DB is short **76 mothers across 35 of
the 103 50 m components** — the events in those components either
lack germplasmIDs entirely or hold fewer than the component's
allele-detection target. Those 76 mothers are the field-team top-up
target for the 2026 season, tracked in the `component_gap` column
of the lab recipe.

### B.5 Two-year design

LEPA fieldwork spans multiple years and the LEPA DB indexes each event
by its `eventDate`. Every input to Steps 28–29 (`n_fertile`,
`seeds_est`, `M_actual`, spatial neighbourhood) is a per-year quantity:
adult plants, seed yields, and mating context all change between
seasons. The correct procedure is therefore to **run the sampling
design once per year**, using the same formulas — no methodological
change is required.

Two implementation rules keep the two-year design honest:

- **Spatial neighbourhood is year-scoped.** When computing K_spatial for
  a mother sampled in year Y, restrict the neighbourhood sum to events
  of the same year Y. This prevents double-counting adults that appear
  in multiple survey years at the same slickspot (which are
  geographically identical and would otherwise inflate K_spatial).
- **Pooling for inference is question-specific.** Standing SRK diversity
  per location is pooled across years (a stationary quantity);
  the mate-limitation test is *not* pooled but treated as repeated
  measures per location (with year as a fixed effect and location as
  a random effect).

The `--year YYYY` flag (documented under § A.2) is the single interface
for running either Step 28 or Step 29 year by year:

```
python step28_seed_sampling_per_mother.py --year 2025
python step29_event_location_sampling.py --year 2025
```

For the two-year analysis, run each step once per year and roll up in
Step 30 according to the pooling rules above.

---

## Part C — Testing predictions with observed SRK data (Phase B)

Part C describes the two statistical tests and the comparison outputs
that fire once **real seed genotypes** are available. Every Part C
artefact is filename-prefixed `step30_B_*` (real data) or
`step30_B_DEMO_*` (synthetic Phase B for pipeline validation, always
carrying a diagonal DEMO watermark). This part closes the loop on the
four scientific goals of the framework.

Section § C.0 provides a **preliminary empirical validation** of the
Part A sporophytic model using adult SRK genotypes that are already in
hand from the Canu-amplicon Phase 4 pipeline — well before any Phase B
seed data are available. It is not one of the two goal-linked tests
(those need seed genotypes), but a sanity check that the model reproduces
the data it was built on and does not fall over at any single EO.

### C.0 EO-scale pre-check — background numbers for the split EOs

The primary Part C validation runs at **Phase 5 location scale** on
the three clean-overlap EOs (§ C.0.a / § C.0.b below). The other
three Phase 4 EOs with n ≥ 10 functional-carrier adults — **EO18,
EO25, EO27** — are **split under the 500 m rule** (see § A.2 /
"Within-EO location split") and cannot yet be remapped to Phase 5
locationCodes without an `Individual → germplasmID → eventID`
join. Until that remap lands, their observed-vs-predicted pollen
compatibility is reported here at EO scale as a background
diagnostic, not a figure.

**Setup.** Per-individual SRK genotypes from
`Tables/Phase4/step23_individual_allele_genotypes_with_nulls.tsv`.
For each EO, observed per-mother pollen compatibility uses fathers
drawn from the **observed local Fg frequency vector**; predicted uses
the species-wide P1 prior. Mothers are held fixed between the two.
Both use the sporophytic Class I / II dominance rules from § A.7.3
and the empirical LEPA zygosity distribution from § A.7.3a.
Confidence intervals come from a mother-resampling bootstrap.

**Current numbers for the three split EOs** (from `step30c_srk_validation_at_eo_level.py`):

| EO | adults | Observed P_compat (95 % CI) | Predicted P_compat (95 % CI) | Interpretation |
|---|---:|---|---|---|
| EO27 | 61 | 0.77 (0.73–0.80) | 0.80 (0.77–0.84) | CIs overlap — model passes at EO scale |
| EO25 | 51 | 0.77 (0.74–0.80) | 0.80 (0.77–0.84) | CIs overlap — model passes at EO scale |
| EO18 | 39 | 0.69 (0.64–0.74) | 0.75 (0.71–0.79) | Mild drift signal, CIs touch |

The three clean-overlap EOs (EO67, EO70, EO76) are **not** reported
here — they are promoted to § C.0.a at Phase 5 location scale and
decomposed in § C.0.b.

**Caveat.** P1 was built from these same individuals, so species-mean
alignment is guaranteed; what this comparison tests is **robustness
to per-EO drift**, not absolute calibration. Once the split EOs can
be remapped to Phase 5 locationCodes, they will join § C.0.a / § C.0.b
and this background section will fold into the output map.

**Outputs.**

- Table: `tables/Phase5/step30_C_pcompat_validation_at_eo.tsv` — one
  row per EO, with observed + predicted P_compat mean and 95 % CI,
  observed Class I mass, and observed distinct-identity distribution.
- Script: `step30c_srk_validation_at_eo_level.py`.

### C.0.a Part C anchor at Phase 5 location scale (clean-overlap EOs)

Three of the § C.0 EOs are **1:1 with a Phase 5 locationCode** — no
500 m within-EO split — so their Phase 4 adult SRK genotypes can
feed Part C **directly** at the Phase 5 location scale, with no
event-level remap required. The clean set is **EO67 (37 adults),
EO70 (74), EO76 (76) = 187 adults**. The other three § C.0 EOs
(EO18, EO25, EO27) are split under the 500 m rule and need an
`Individual → germplasmID → eventID → locationCode` join before
their 151 adults can be used at the finer scale; that is scheduled
as future work.

**Script `step30d_partC_clean_overlap.py`** reads
`Tables/Phase4/step23_individual_allele_genotypes_with_nulls.tsv`,
filters to the three clean EOs, and for each location computes:

- **Observed per-mother P_compat** (sporophytic Class I/II +
  empirical zygosity, fathers drawn from the OBSERVED local Fg
  frequency vector) and **observed distinct Fg count**.
- A **no-drift upper bound on the observed Fg count** at the adult
  sample size, `E[distinct Fgs | 4·n_adults draws from P1]`. At
  n ≥ 37 adults the bound sits ≥ 19 Fgs, so any large deficit
  from the bound to the observed count isolates **drift**, not
  sampling.
- Comparison against the Phase 5 **per-location** prediction
  (`step30_A_prediction_location_pcompat.tsv` for P_compat,
  `step30_A_prediction_location_diversity.tsv` for the component-
  unioned pool size).

**Result.** EO67 — observed and Phase 5 CIs overlap on P_compat;
observed Fg count sits inside the Phase 5 CI → model passes. **EO70
— observed P_compat 0.53 vs predicted 0.78 (CIs do not overlap);
observed 6/32 Fgs vs no-drift upper bound 24 → massive drift
collapse at a large location**, confirming the § C.0 EO-scale
finding at the Phase 5 location scale. **EO76 — observed 9/32 Fgs
vs no-drift upper bound 24 → severe diversity collapse at LEPA's
single largest site**, though observed P_compat (0.72) stays close
to Phase 5's 0.78 because Class II (6 Fgs, FG001 at 41 %) dominates
and keeps most mothers compatible even with limited local diversity.

**Outputs.**

- Table: `tables/Phase5/step30_B_partC_clean_overlap_per_location.tsv`
  — one row per clean location with observed vs Phase 5 predicted
  pollen compatibility, observed vs Phase 5 predicted distinct allele
  count, and the no-drift upper bound at the adult sample size.
- Table: `tables/Phase5/step30_B_partC_clean_overlap_fg_frequencies.tsv`
  — long form (locationCode × allele) with observed f, species-wide
  P1 f, and present/absent flag.
- Figure 10: `figures/Phase5/step30_B_partC_clean_overlap.png` — see
  below.
- Script: `step30d_partC_clean_overlap.py`.

![Figure 10](figures/Phase5/step30_B_partC_clean_overlap.png)

**Figure 10.** Phase 5 Part C anchor at the three clean-overlap EOs (EO67, EO70, EO76; 1:1 with a Phase 5 locationCode). Panel order follows the biological flow: SRK diversity → pollen compatibility → per-allele drift residual. **Panel A** — predicted vs observed number of distinct SRK alleles per location (ceiling = 32), 1:1 diagonal; EO67 sits on the diagonal, EO70 and EO76 fall far below (6/32 and 9/32 observed). **Panel B** — predicted vs observed pollen compatibility under random mating, with the traffic-light background bands (red = failed < 0.26, amber = struggling < 0.52, green = sustainable ≥ 0.52). EO67 and EO76 near the diagonal; EO70 is the clear outlier (observed 0.53 vs predicted 0.78). **Panel C** — allele frequency residual per allele per location (one row each), alleles sorted by species-wide frequency most-common → rarest; **green = positive residual (observed > species-wide)**, **red = negative (depleted)**, **× = absent at this location**. The residual exposes the drift signature that the raw-frequency plot blurred — EO70 shows classical FG001 (+21 %) and FG024 (+18 %) enrichment with 26/32 alleles absent; EO76 has milder enrichment of FG012 / FG010; EO67 elevates the normally-rare FG018 and FG023 instead of the common alleles — a founder-effect signature. Source: `step30d_partC_clean_overlap.py`.

**Why P1 is empirical, not uniform.** A uniform-frequency P1
(Dirichlet(α = 1) over the 32 Fgs, equivalent to the P0 "uninformative"
prior kept as a reference baseline in § A.3) would represent a
*neutral null with no drift history*. That is not LEPA: decades of
habitat loss have already pushed the species through drift, so
**empirical P1 is the realistic starting point** against which per-
location further drift is measured. Panel C of Figure 7 makes this
concrete — EO70's local Fg pool has moved further toward FG001
dominance (0.61) than P1 itself (0.41), so even the already-skewed
empirical prior underestimates the local pool's collapse. A uniform
P1 would start from a flat distribution and declare every observed
location "drifted", losing both the species-wide signal and the
calibration readers need to interpret per-location deviation.

### C.0.b Hypothesis decomposition — diversity collapse → pollen compatibility

The § C.0.a Part C anchor exposed a disconnect that pollen
compatibility alone would have hidden: EO70 and EO76 both lose
~23 of their predicted alleles, but observed pollen compatibility
drops by 0.25 at EO70 and barely 0.02 at EO76. EO67, with just 7
alleles against a predicted 10.8, lands at 0.70 — in the
sustainable band. The script `step30e_pcompat_hypothesis_decomposition.py`
decomposes each location's observed → species-wide gap into three
competing mechanisms by swapping ONE factor at a time to the
species-wide reference while keeping observed mother genotypes
fixed:

1. **Within-class allele spread.** Same class totals, but P1-shaped
   within each class. Isolates drift concentration (e.g. FG024 at
   EO70 absorbing nearly all of Class I's mass).
2. **Class I / II mass balance.** Same within-class shape, but
   class totals rescaled to the species-wide values. Isolates
   between-class-rescue capacity.
3. **Zygosity composition.** Same local allele frequencies, but
   father zygosity drawn from the species-wide distribution.
   Isolates the per-mother expressed-set-size channel.

The driver fraction = (counterfactual − observed) / (prediction −
observed): a value near +100 % says the swapped factor by itself
explains the whole gap; a negative value says the factor is
**buffering** the location.

**Headline numbers (from the current run).**

- **EO70 — classical within-class drift collapse.** Within-class
  spread explains **+88 %** of the pollen compatibility deficit;
  class balance and zygosity contribute within noise. Mechanistic
  story: FG024 at 35 % of the pool makes a FG024-homozygous Class
  I mother see `(1 − 0.35)⁴ ≈ 0.17` compatibility vs ~0.43 at
  EO67 where Class I is split four ways.
- **EO67 — looks bad on paper, pollen compatibility survives.** 7
  of 32 alleles present, but 4 of them are Class I (FG024 / FG018
  / FG023 / FG016, well spread) and 3 are Class II. Within-class
  spread explains +100 % of the small observed → prediction gap;
  zygosity (slightly higher multi-identity than species-wide) is
  mildly buffering.
- **EO76 — predicted near-perfect, observed severely collapsed,
  pollen compatibility still fine.** 9 of 32 alleles. **Zygosity
  explains −115 % of the gap** (actively buffering): the location
  is 76 % homozygous and its dominant Class II allele FG001 sits
  at 46 % (vs EO70's 62 %), so FG001-homozygous mothers face many
  compatible fathers. Swapping to the species-wide zygosity would
  *lower* pollen compatibility at EO76 by 0.03.

**Why the sporophytic Class I / II + empirical-zygosity model
matters.** The decomposition is only meaningful because the model
can resolve the three channels separately. A simpler diploid
gametophytic approximation would collapse all three into a single
"effective diversity" number and mis-predict the EO76 outcome.
The hypothesis test therefore doubles as validation of the model's
mechanistic structure.

**Biological take-home — the deme buffers intense drift.** The
three channels above are not independent model knobs but three
**buffering mechanisms stacked on top of each other** that let a
deme keep breeding after severe allele loss:

1. **Class I dominance (26 of 32 Fgs).** Every between-class cross
   is compatible by construction, and Class I > Class II within a
   plant, so a mother carrying any Class I allele has her
   compatibility set primarily by that allele. Losing *rare*
   Class I alleles barely moves the mean because the common Class I
   alleles still carry the pool.
2. **Frequency-dependent selection is forgiving when drift
   collapses onto *common* alleles.** The classical SI catastrophe
   (Lawrence 2000; Castric & Vekemans 2004) is a deme that collapses
   onto so few alleles that most mothers share them all → crash.
   At EO70 the 6 surviving Fgs are the species's *common* ones
   (FG001 at 41 %, FG024 at 35 % locally) — the opposite of the
   worst case. Different mothers still carry different combinations,
   so the pollen pool still finds compatible targets.
3. **Tetraploid zygosity actively buffers at homozygous-rich sites.**
   A homozygous mother expresses a single Fg at her stigma, which
   pollen fathers can more easily avoid matching than a
   heterozygous mother's larger expressed set. The yellow bar at
   EO76 in Figure 11 **drops below the observed black bar** — the
   quantitative signature of this mechanism.

**What this reverses in the standard conservation-genetics
intuition.** The usual story *lose SRK alleles → mating failure*
treats allele count as the bottom line. In a **tetraploid +
sporophytic Class I / II system with the specific class imbalance
LEPA has**, the mating-level consequence of allele loss is
**actively buffered** until the pool either (i) collapses onto a
single class or (ii) homozygosity becomes so extreme that the
expressed-set channel constrains compatibility. The 32 → 6 drift
collapse at EO70 is severe by any count, and yet the pollen pool
still works at 70 % compatibility (close to the species mean
0.78). **The system has structural redundancy that diversity-
counting alone cannot see** — a central Phase 5 result and a
reason the pollen-compatibility metric, not raw allele count, is
the one that enters the § C.1 mate-limitation regression. The
B1 pilot pair in § C.6 will tell us whether any location is near
either tipping point.

**Outputs.**

- Table: `tables/Phase5/step30_B_partC_hypothesis_decomposition.tsv`
  — per-location diversity trigger + three single-swap pollen
  compatibility counterfactuals with bootstrap 95 % CIs + driver
  fractions.
- Figure 11: `figures/Phase5/step30_B_partC_hypothesis_decomposition.png`
  — see below.
- Script: `step30e_pcompat_hypothesis_decomposition.py`.

![Figure 11](figures/Phase5/step30_B_partC_hypothesis_decomposition.png)

**Figure 11.** Competing-hypothesis decomposition of per-location pollen compatibility, triggered by the SRK diversity discrepancy in Figure 10 Panel A. For each clean-overlap location (EO67, EO70, EO76) five pollen-compatibility values are shown from the same simulation (sporophytic Class I / II + empirical zygosity, 800 candidate fathers per mother, observed mother genotypes held fixed), differing only in which part of the father-drawing distribution is swapped to the species-wide reference. **Green — Predicted (species-wide):** fathers from P1 + species-wide zygosity. **Black — Observed:** fathers from observed local allele frequencies and observed local zygosity. **Blue — Swap within-class spread:** keep observed Class I and Class II total masses, reshape the within-class spread to P1 (isolates drift concentration). **Red — Swap Class I / II balance:** keep within-class shape observed, rescale class totals to species-wide (isolates between-class rescue). **Yellow — Swap zygosity:** keep observed allele frequencies, swap father zygosity to species-wide 66/32/2 % (isolates the per-mother expressed-set channel). Error bars = bootstrap 95 % CI. Traffic-light bands shaded in the background; dashed grey line = species-wide pollen compatibility 0.78. **Gold star** marks the dominant single-factor explanation of the Predicted → Observed gap at each location. Reading rules: a blue/red/yellow bar jumping toward green = that factor caused the gap; staying next to black = not involved; dropping below black = **buffering** (observed state better than species-wide). **EO70** — classical within-class drift collapse (FG024 absorbs 35 % of pool). **EO67** — small deficit explained by within-class spread; zygosity slightly buffering. **EO76** — small deficit despite 23/32 alleles lost; **zygosity actively buffers** (yellow drops below black), and within-class spread overshoots because the specific alleles present at EO76 are better-spread than drift-only P1 would deliver. Source: `step30e_pcompat_hypothesis_decomposition.py`.

### C.0.c Hypothesis test — can the diversity gap be closed by tightening the deme radius?

**Why this test exists.** Figure 10 Panel A shows the Phase 5 50 m
diversity prediction **over-estimates** observed SRK allele counts
at the clean-overlap locations: predicted 11 / 28 / 32 vs observed
7 / 6 / 9 at EO67 / EO70 / EO76. § C.0.b above explains why pollen
compatibility nevertheless tracks observation; here we address the
complementary question — **can the diversity gap at EO70 and EO76
be explained by the operational deme being too wide at 50 m?**

**Two competing scenarios.**

- **H1 — the 50 m operational deme is too wide.** Realised gene
  flow is tighter than 50 m; a smaller radius partitions the
  location into more and smaller demes with (within the current
  model) independent drift histories drawn from the species-wide
  prior P1. Smaller demes draw fewer alleles per deme, and even
  after union across more demes the predicted location pool should
  shrink. **Diagnostic signature:** predicted diversity drops
  toward observed as radius tightens.
- **H2 — per-deme drift history exceeds the species-wide prior
  P1.** Each deme has drifted *past* what P1 captures, so no
  spatial repartition at any radius closes the gap under the
  current model. **Diagnostic signature:** predicted diversity
  stays flat and well above observed at every radius.

**Scope of the test.** The sweep rebuilds the deme partition at
radii 10, 25, 50, 75, 100, 150 m and re-runs the per-deme SRK
diversity simulation (2000 replicates per combination). Each
sub-deme draws `PLOIDY × component_N_fertile` alleles from P1 —
i.e. the test varies the *spatial partition* while holding the
*drift model* fixed. **This test can only reject H1 within the P1
drift assumption; it does not simulate H2 directly** (which would
require per-deme empirical frequency vectors that do not exist
yet; the Part C seed-genotyping campaign is designed to deliver
those vectors).

**Preamble — EO67 is omitted from the sweep.** EO67 has only two
events in the 2025 record, placed > 150 m apart. The deme
partition is therefore **invariant at every tested radius** (2
demes everywhere), making the radius sweep non-informative at
EO67. The observed count at EO67 (7) also sits inside the 95 % CI
of the Phase 5 prediction — there is no gap to explain. EO67 is
retained in the TSV for completeness, but the H1-vs-H2 figure
restricts to **EO70** and **EO76**, which do show meaningful
fragmentation across the sweep (EO70: 1 → 3 demes; EO76: 5 → 17
demes as radius tightens).

**Result — H1 is rejected within the P1 drift assumption at EO70
and EO76** ([Figure 13](#fig-13)).

| Location | Observed | Predicted at r = 10 m | at 50 m | at 150 m | Deme count 10 m / 50 m / 150 m |
|---|:---:|:---:|:---:|:---:|:---:|
| **EO70** | 6 | 28.3 [25, 31] | 28.3 [25, 31] | 28.4 [25, 31] | 3 / 1 / 1 |
| **EO76** | 9 | 31.8 [31, 32] | 31.8 [31, 32] | 31.8 [31, 32] | 17 / 6 / 5 |

95 % CIs in brackets. The prediction **does not move** across the
sweep at either location: even fragmenting EO76 into 17 small
demes at 10 m, or EO70 into 3 demes at 10 m, leaves the union-of-
per-deme pools at essentially the same 28–32 Fgs. The ~22-allele
gap at EO70 and ~23-allele gap at EO76 do not close.

**Why the prediction is flat across the sweep — the mechanics.**
The location-level prediction is the **union** of per-deme Fg sets,
not a weighted sum. When a radius change re-partitions a location
into more or fewer sub-demes, the **total number of allele copies
sampled at that location stays fixed at `PLOIDY × N_fertile_total`**
— the sweep only redistributes those draws across sub-demes. Each
sub-deme on its own may miss rare Fgs under P1 (where FG001 ≈ 41 %,
five more Fgs ≈ 25 % combined, and 20+ rare Fgs are each < 2 %),
but a rare Fg at frequency *p* now has **multiple sub-demes to
appear in**, and the union captures it with probability
`1 − (1 − p)^(PLOIDY × N_total)` — which is governed by the
**location total, not the partition**. The partition therefore only
matters when (a) `N_fertile_total` is already below the P1
saturation ceiling (~50–100 adults, where the full-pool draw itself
starts missing rare Fgs — e.g. EO67's 10 adults), **or** (b)
sub-deme sizes drop into the low single digits AND the location
total is also near the ceiling. Neither applies at EO70 (161
plants, 644 copies) or EO76 (517 plants, 2068 copies): both are
far above the ceiling at every radius in the sweep. Concretely, at
EO70's 10 m extreme the 161 plants split into sub-demes of 78 / 65 /
18 that each recover ~26 / 25 / 19 Fgs independently — but their
union still delivers 28, same as the single 161-plant deme at 50 m.
This is a general feature, not a parameter artefact: **under a
drift model that is identical across sub-demes, spatial sub-
partitioning washes out at the location level once the total pool
saturates P1**. Varying the partition can therefore never close the
diversity gap under P1 — only a drift model that differs per deme
can (i.e. H2, which requires per-deme empirical priors from Part C
seed-genotyping).

**Interpretation — what we can and cannot say.**

- **Can say — the 50 m operational deme is defensible as a first-
  pass partition for Phase 5 predictions.** The sweep gives no
  evidence that a tighter radius would reconcile predicted with
  observed diversity under the current model.
- **Can say — the diversity gap at EO70 and EO76 cannot be a pure
  spatial-partitioning artefact.** The only way the gap closes
  under the current model is if drift histories *inside* each
  deme have diverged from P1, i.e. if each deme carries its own
  frequency vector that P1 does not capture. That is the H2
  scenario. **In biological terms: the sweep rules out a
  geometric explanation for the diversity gap and leaves the
  within-deme biological one — prolonged drift, founder effects,
  or local bottlenecks that have stripped alleles beyond what the
  species-wide prior encodes. The deme is still the right unit of
  inference; what the data say is that the drift history inside
  each deme is more severe than a single species-wide frequency
  vector captures, which is the central prediction of the
  fragmentation → drift → mate-limitation causal chain this
  framework was built to test.**
- **Cannot say — the sweep alone proves H2.** We have not
  simulated H2 (per-deme frequency vectors drifted past P1); the
  Part C seed-genotyping campaign is designed to measure those
  per-deme vectors empirically and complete the test.
- **Direction of error is favourable.** Over-predicting at a
  radius that is already a geographic upper bound (50 m) says real
  demes are **at most** 50 m and could be tighter. If a future
  refinement — pairwise kinship on Part C seedlings (Hardy &
  Vekemans 1999; Vekemans & Hardy 2004), or dropping the step
  function for a continuous dispersal kernel — were to shrink the
  deme, it would only *add* demes and therefore mothers per
  location, giving us finer-scale information for free. The
  opposite error (under-predicting ⇒ demes too narrow ⇒ false-
  positive fragmentation) would waste sampling on sub-demes that
  do not exist, so the current direction is the one working in
  our favour.

**Output files.**

- Table: `tables/Phase5/step30f_srk_diversity_radius_sweep.tsv`
  — rows for all three EOs × six radii with `pred_mean` /
  `pred_lo95` / `pred_hi95`, `n_components`, `component_sizes`,
  `observed`, `obs_inside_CI`.
- Figure 13: `figures/Phase5/step30f_srk_diversity_radius_sweep.png`
  — gap-focused single panel (y = predicted − observed); EO67
  omitted (invariant deme count across sweep).
- Script: `step30f_srk_diversity_hypothesis_test.py`.

<a id="fig-13"></a>
![Figure 13](figures/Phase5/step30f_srk_diversity_radius_sweep.png)

**Figure 13.** SRK diversity gap (y = predicted − observed distinct SRK alleles) vs pollinator radius at EO70 and EO76. x = radius used to rebuild the deme partition (10 → 150 m, log scale). Lines + 95 % CI bands = per-location per-radius simulation (2000 replicates per combination). Right-side labels show the gap at the largest radius; parenthetical observed counts come from the Phase 4 adult SRK genotypes (a single measurement, not a sweep output — hence shown as labels, not as horizontal lines). Zero line = perfect match. Vertical dotted line at 50 m marks the current operational deme. H1 (deme too wide) would predict the gap to shrink toward zero as radius tightens; it stays flat at +22 and +23 across all radii, including the 10 m extreme that fragments EO70 into 3 demes and EO76 into 17 — rejecting H1 within the P1 drift assumption. EO67 is omitted from the figure because its deme partition is invariant across the sweep (only two events in 2025, > 150 m apart), making the test non-informative there; see § C.0.c text. Source: `step30f_srk_diversity_hypothesis_test.py`. Data: [`step30f_srk_diversity_radius_sweep.tsv`](tables/Phase5/step30f_srk_diversity_radius_sweep.tsv).

### C.1 Test 1 — Mate-limitation regression (goals 1 + 2)

Under strict SI + random mating in a **tetraploid sporophytic** system
(§ A.3, § A.8), a mother's per-mother compatibility is given by the
closed-form § A.8 Case-A / Case-B formulas, driven by the local Fg
frequency vector *f* and the Class I / II assignment in
[`srk_fg_class.tsv`](tables/Phase5/srk_fg_class.tsv):

- **Case A** — mother has ≥ 1 Class I allele; expressed set M ⊆
  Class I; $P_{\text{compat}}(m) = (1 - p(M))^4$.
- **Case B** — mother is all Class II; expressed set M = her 4 alleles;
  $P_{\text{compat}}(m) = 1 - (1 - p_I)^4 + (1 - p_I - p(M))^4$.

The species-mean sporophytic + empirical-zygosity pollen compatibility is ~0.78 (reflecting the Class I majority: with 26 of 32 Fgs in the dominant class, between-class crosses — always compatible — dominate the mean. The
diploid gametophytic 0.63), so the regression coefficient β₁ lives on
that rescaled axis — the traffic-light bands in § A.8 anchor the
interpretation. Her expected seed set is proportional to
$P_{\text{compat}}(m)$:

$$E[\text{seeds}_m] \;\propto\; \text{ovules}_m \times P_{\text{compat}}(m)$$

The **test** is a mixed-effects regression of observed
`germplasmQuantityEstimate` on predicted $P_{\text{compat}}$:

$$\text{seeds}_m \;=\; \beta_0 + \beta_1 \cdot P_{\text{compat}}(m) + \beta_2 \cdot K^{(50\text{m})}(m) + u_{\text{location}(m)} + u_{\text{year}(m)} + \epsilon_m$$

where $K^{(50\text{m})}$ is the mother's mating-neighbourhood
pollen-donor pool at the 50 m radius (§ B.3). Two coefficients that
together decompose the fragmentation effect into a drift-mediated
pathway and a direct pathway:

- **β₁ > 0 with 95 % CI excluding 0 → the complete fragmentation → drift →
  mate-limitation chain fires.** Fragmentation reduced
  `N_fertile_effective`; small N_e drove allele-frequency drift; the
  drifted local Fg pool reduces `P_compat`; reduced `P_compat`
  reduces seed set. Under the sporophytic Class I / II model,
  between-class compatibility (Class I × Class II always compatible)
  buffers most locations against drift (§ A.8), so β₁ is expected to
  fire at locations where drift has removed a whole class from the
  local pool — this signal is subtler than a diploid gametophytic
  version would have implied.
- **β₂ > 0 conditional on β₁ → fragmentation has additional effects
  beyond the drift chain.** Spatial isolation reduces seed set above
  and beyond what allele skew alone explains — for example, sparser
  neighbouring adults may attract fewer pollinator visits regardless
  of which alleles they carry.
- **β₁ > 0 with β₂ ≈ 0 → the drift chain fully explains the
  fragmentation effect.** The fragmentation effect on seed set is
  entirely captured by its downstream consequence on `P_compat`;
  no non-drift pathway is needed.
- **β₁ ≈ 0 with β₂ > 0 → fragmentation acts entirely outside the
  drift channel.** Small locations lose seed set for reasons other
  than SRK-mediated mate limitation (e.g. pollinator behaviour).
  Would indicate a model refinement: either the SRK Fg → Class map
  needs revising, or non-genetic mechanisms need to be added.

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

### C.2 SI-escape rate test (moved to Appendix — out of scope for Phase 5)

The SI-escape permutation test that lived here in earlier drafts is
now in the [Appendix](#appendix--out-of-scope--future-work). The DEMO
scaffold is retained in the Step 30 code for pipeline-validation
purposes; the SI-escape analysis is not part of the Phase 5
mate-limitation + fragmentation framework. See the appendix for the
model, the null hypothesis, the permutation procedure, and the DEMO
figure.

### C.3 Step 30 script behaviour (CLI + phase filenames)

**Script.** [`step30_srk_diversity_prediction_vs_observed.py`](step30_srk_diversity_prediction_vs_observed.py)

The script always produces the prediction family. If a `--seed-genotypes`
TSV is provided (columns `locationID`, `eventID`, `germplasmID`, `seed_id`,
`maternal_Fg`, `paternal_Fg`), it additionally produces the comparison
family. If no such TSV exists but the flag `--demo` is passed, the script
simulates a plausible seed-genotype dataset from the prior itself — useful
to verify the comparison pipeline end-to-end before real seed data arrive.

### C.4 The clean payoff

- **Before seed data are back**: publishable predicted SRK diversity,
  predicted per-mother $P_{\text{compat}}$, and predicted
  fecundation-failure rates per location — with honest credible intervals.
  These can be pre-registered.
- **When seed data land**: the mate-limitation regression (§ C.1)
  fires from the seed-genotype input and directly tests whether
  predicted per-location pollen compatibility explains observed seed set.
- **Fragmentation × drift decomposition**: the two coefficients
  $\beta_1$ (pollen-compatibility effect) and $\beta_2$
  (mating-neighbourhood-size effect) from the mate-limitation regression
  separate the two mechanisms even though they both reduce reproductive
  success — a distinction unreachable with any single-year,
  single-location analysis.
- **Two-generation efficiency**: every mother's seed lot pays double —
  it certifies her genotype (maternal inventory) and samples her pollen
  environment (paternal inventory) in one experiment.
- **Phenotype cross-validation is a natural next step**: the per-location
  predictions can be overlaid with ISI / fruit set from Genetic-Rescue-DB
  without changing any of the code paths in this framework.

### C.5 How Phase B data are generated and consumed

Every Part A prediction becomes testable once we have **seed-DNA
genotypes** at SRK. The two-generation trick makes each mother's seed
lot pay double — the mother's own genotype falls out of her sib set
(constant / 50 %-frequency alleles), and the paternal alleles of those
seeds sample the pollen environment that fertilised her ovules. One
seed lot per mother = one maternal genotype + n paternal alleles.

**Data-generation pipeline.**

1. **Seed extractions.** For each mother in the field-team recipe
   ([`step29_field_team_sampling_recipe.tsv`](tables/Phase5/step29_field_team_sampling_recipe.tsv)),
   genotype the recommended number of seeds (**typically 15 under
   the tetraploid per-mother cap**, capped at her real seed budget; 2 paternal
   allele draws per seed). Total per-year effort: 765 mothers × ~15
   seeds ≈ 11 500 seed genotypes. (Under the fragmentation-aware
   § B.4.2 allocation the mother count roughly doubles, so the total
   is closer to 25 000 seed genotypes.)
2. **SRK amplicon sequencing.** Any protocol that resolves the 32 Fgs.
   The Canu-amplicon pipeline is already validated on adult plants and
   applies directly to seed material.
3. **Fg assignment.** Each observed SRK allele is mapped to one of the
   32 Fgs from the preliminary study
   ([[project_functional_srk_definition]]). Alleles that do not map
   are flagged as candidates for post-hoc expansion of the Fg set.
4. **Two output tables**, which the Phase B pipeline consumes without
   any code changes:
   - `real_seed_genotypes.tsv` — one row per seed: `locationID`,
     `eventID`, `germplasmID`, `seed_id`, `maternal_Fg`, `paternal_Fg`.
   - `real_mother_genotypes.tsv` — one row per mother: `locationID`,
     `germplasmID`, `mother_Fg_a`, `mother_Fg_b`, `K_spatial`,
     `seeds_est`. Everything except the two `mother_Fg_*` columns is
     already in the Step 28 output; the two Fg columns come directly
     from the inferred maternal genotype in step 1.

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
exist. Below is what the two tests **will look like** once real seed
genotypes exist — these are DEMO renders produced under `--demo`
(simulated Phase B data), watermarked so they cannot be confused with
real analysis.

<a id="fig-14"></a>
![Figure 14](figures/Phase5/step30_B_DEMO_mate_limitation.png)

**Figure 14 (DEMO).** Mate-limitation regression preview. One dot per LEPA location, colour = "sustainable" band. X = predicted random-mating pollen compatibility (mean across sampled mothers at that location); Y = mean observed seeds per mother. Traffic-light background bands (failed / struggling / sustainable). Dashed line = weighted OLS fit, slope + p-value printed in the legend. In this DEMO the simulator baked in a direct causal link (seed set ∝ compatibility), so the slope is highly significant. **With real data the same figure will test whether observed seed set actually declines with predicted compatibility — a positive slope with 95 % CI excluding 0 confirms mate limitation at population level.** Source: `step30_srk_diversity_prediction_vs_observed.py --demo`.

*(A DEMO Figure A1 for the SI-escape rate test exists at
[`step30_B_DEMO_si_escape.png`](figures/Phase5/step30_B_DEMO_si_escape.png).
It is described in the [Appendix](#appendix--out-of-scope--future-work);
the SI-escape analysis is not part of the Phase 5 mate-limitation +
fragmentation framework.)*

**Two-year extension.** The `--year` flag on Steps 28, 29 and 30
makes the pipeline year-scoped. When 2026 field data arrives, run
each step twice (once per year); the pooling rules under § B.5 tell
Phase B how to combine years for standing-diversity inference vs how
to keep them separate as repeated measures for the mate-limitation test.

**Downstream: cross-validation with phenotype.** Per-location
mate-limitation results can be overlaid with the Genetic-Rescue-DB
ISI / fruit-set phenotype to provide an independent line of evidence
under Phase III of the wider project.

### C.6 BL5 pilot study — EO48_7 (connected) + EO18-7_19 (fragmented)

Before running Phase B across all 39 locations, we recommend a
**two-location pilot within Bottleneck Lineage 5 (BL5)**: one
**fully-connected** location (**EO48_7**, 98 adults in a single
50 m deme, 9 mothers in the LEPA DB) paired with one
**fragmented** location (**EO18-7_19**, 34 adults distributed across
3 demes, 11 mothers in the DB). User-selected 2026-10-03
(option A2 in
[`Tables/Phase5/step29c_partC_BL5_pilot_candidates.tsv`](tables/Phase5/step29c_partC_BL5_pilot_candidates.tsv)).

**Why BL5.** It is the LEPA Bottleneck Lineage with the **widest
within-BL variation** in both census size and deme structure
(see § A.5 Figure 1b): BL5 holds the drift-collapsed tail (EO24
group, 1 – 3 adult singletons) *and* the largest, most-connected
locations (EO32_6 at 466 adults; EO48_7 at 98 adults in a single
deme). A within-BL5 contrast therefore rules out between-BL
noise while testing the full span of the fragmentation × drift axis
the model predicts matters.

**The two pilot locations.**

| Attribute | **EO48_7** (connected) | **EO18-7_19** (fragmented) |
|---|---|---|
| 50 m demes (`n_components_50m`) | **1** | **3** |
| Total fertile adults | 98 | 34 |
| Largest-pool `component_N_fertile` | 98 | 20 |
| Mothers with seed records in DB | **9** | **11** |
| Seedling-genotype target (§ B.3.2) | 15 seedlings/mother | 15 seedlings/mother |
| Seeds to germinate per mother (60 % germination) | 25 seeds | 25 seeds |
| Predicted per-location pollen compatibility (sporophytic + empirical zygosity, § A.8) | **sustainable** — single-pool model's cleanest test case | **sustainable mean, wider CI** — three independent small-pool drift experiments averaged by `component_N_fertile` weight |

**What this specific pair tests.**

- **EO48_7 — fully-connected regime** (single deme,
  component_N_fertile = 98). Phase 5 predicts a sustainable location
  because the location is unfragmented; the deme-based
  simulation reduces to a single pool. This location is the cleanest
  test of the model's "no fragmentation → species mean" prediction.
- **EO18-7_19 — fragmented regime** at a similar order of magnitude
  for total adults (34) but split across 3 demes (20, 12, 2
  adults — two pools above the 8-plant species-pool threshold, one
  below). The per-component simulation treats each as an independent
  drift experiment; the comparison with EO48_7 isolates the pure
  fragmentation effect from the raw-size effect.

**Caveat on retrospective adult SRK data.** The two BL5 EOs that
*do* have adult SRK genotypes from Phase 4 at n ≥ 10 (EO25 and EO18)
are both split under the 500 m rule and sit in § C.0 only — they
cannot yet be used at the Phase 5 location scale. So this pilot is a
**seed-genotyping (Phase B) pilot**, not a retrospective-on-adults
pilot like § C.0.a. The 2026 pilot can run as soon as seeds from
EO48 and EO18-7 are collected and germinated.

**Pilot cost.** 20 mothers × 25 seeds = **500 seeds to germinate**
→ ~15 seedlings/mother × 20 = **~300 seedlings to genotype**
at 60 % germination (§ B.3.2). That is **< 4 % of the full 2026
genotyping budget** (12 625 seeds / ~7 575 seedlings) for a
within-BL contrast with a clear a priori hypothesis.

**Three options kept on file** (full table in
[`Tables/Phase5/step29c_partC_BL5_pilot_candidates.tsv`](tables/Phase5/step29c_partC_BL5_pilot_candidates.tsv),
generator `build_bl5_pilot_candidates.py`):

| Option | Connected candidate | Fragmented / drift-sensitive candidate | Mothers in DB | Design note |
|---|---|---|---:|---|
| A1 | EO32_6 (5 pools, 466 adults) | EO25-B_21 (2 pools, 10 adults) | 38 + 7 | Maximum size contrast; drift-sensitive is small but not tiny. |
| **A2 ★** | **EO48_7 (1 pool, 98 adults)** | **EO18-7_19 (3 pools, 34 adults)** | **9 + 11** | **User-recommended.** Isolates fragmentation × drift at matched order-of-magnitude total adults. |
| A3 | EO18-7_17 (5 pools, 242 adults) | EO24-7_25 (1 pool, 3 adults) | 33 + 3 | Maximum biological contrast but drift-sensitive has only 3 adults — low genotyping statistical power. |

**Field-team recipe for the pilot.** Under the § B.4.2 default, the
authoritative per-event allocation is `M_frag` in
[`step29c_sampling_frag_aware_per_event.tsv`](tables/Phase5/step29c_sampling_frag_aware_per_event.tsv);
filter to the pilot `locationID`s (EO48_7 → locationID 7, EO18-7_19
→ locationID 19) to isolate the two locations. The lab recipe is
[`step29c_partC_germplasmID_selection.tsv`](tables/Phase5/step29c_partC_germplasmID_selection.tsv),
sorted `EOID → locationID → component → germplasmID`; the per-mother
columns `n_seeds_to_germinate` (= 25) and `n_seedlings_to_genotype`
(= 15 at 60 % germination) are already baked in. The older
[`step29_field_team_sampling_recipe.tsv`](tables/Phase5/step29_field_team_sampling_recipe.tsv)
is kept for reference / audit but no longer drives field effort.

**Companion tables for pilot review:**

- [`step29c_sampling_comparison_per_location.tsv`](tables/Phase5/step29c_sampling_comparison_per_location.tsv) — head-to-head of current vs fragmentation-aware M per location; the two pilot rows show how much extra effort the honest allocation asks for.
- [`step28_seed_sampling_per_mother.tsv`](tables/Phase5/step28_seed_sampling_per_mother.tsv) — analyst view (per-mother targets, achievable coverage, budget flags).
- [`step29_sampling_per_location.tsv`](tables/Phase5/step29_sampling_per_location.tsv) — per-location design table (M, target, delivered coverage) under the older § B.4 allocation, kept for audit.
- [`step29_location_connectivity.tsv`](tables/Phase5/step29_location_connectivity.tsv) — the 50 m primary + 10 m / 50 m sensitivity connectivity metrics.
- [`step30_A_prediction_location_diversity.tsv`](tables/Phase5/step30_A_prediction_location_diversity.tsv) — predicted SRK allele count with 95 % credible interval per location.
- [`step30_A_prediction_location_pcompat.tsv`](tables/Phase5/step30_A_prediction_location_pcompat.tsv) — predicted random-mating pollen compatibility with 95 % credible interval per location (finite-population model at 50 m).

**Downstream in Phase B.** Once EO67 and EO27-1 seed genotypes exist,
run Step 30 with the pilot subset:

```
python step30_srk_diversity_prediction_vs_observed.py \
    --year 2025 \
    --seed-genotypes  real_seeds_BL4_pilot.tsv \
    --mother-genotypes real_mothers_BL4_pilot.tsv \
    --match-seed-count
```

The two-point mate-limitation regression will fire on the two pilot
locations; scaling to the full 39-location dataset is then just a
matter of adding rows to the two TSV inputs.

---

## Output map — quick reference (grouped by phase)

All Phase A tables are direct download links — click the filename to open
or right-click → *Save link as…* to pull the TSV into your local pipeline.

### Phase A — Preliminary (before SRK genotyping)

**Tables** (all TSV):

- [`step28_events_spatial_neighborhood.tsv`](tables/Phase5/step28_events_spatial_neighborhood.tsv) — event spatial reference (lat, lon, per-radius neighbourhood counts).
- [`step28_seed_sampling_per_mother.tsv`](tables/Phase5/step28_seed_sampling_per_mother.tsv) — per-mother analyst view (K, achievable coverage, budget flags).
- [`step28_coverage_curves_by_Nfertile.tsv`](tables/Phase5/step28_coverage_curves_by_Nfertile.tsv) — per-mother analytical + simulation curves by event-size bin.
- [`step28_aggregation_curves_by_Nfertile.tsv`](tables/Phase5/step28_aggregation_curves_by_Nfertile.tsv) — aggregation curves showing expected distinct alleles at a location for M = 1 … 30 mothers × 15 seeds each (tetraploid per-mother cap), one row per (event-size bin, M).
- [`step28_mothers_for_full_detection_by_bin.tsv`](tables/Phase5/step28_mothers_for_full_detection_by_bin.tsv) — uniform-model M for a 90 % chance of observing every one of K_local alleles, one row per event-size bin.
- [`step28_mothers_for_full_detection_by_location.tsv`](tables/Phase5/step28_mothers_for_full_detection_by_location.tsv) — per-location `M_recommended` combining the per-mother allele-detection bound with the private-allele floor (≥ 1 mother per event).
- [`step29c_sampling_frag_aware_per_event.tsv`](tables/Phase5/step29c_sampling_frag_aware_per_event.tsv) — per-event fragmentation-aware allocation under the 50 m connected-component decomposition; columns include `component_id`, `component_K`, `component_M_target`, `M_paternal_only`, `M_frag`.
- [`step29c_partC_BL5_pilot_candidates.tsv`](tables/Phase5/step29c_partC_BL5_pilot_candidates.tsv) — § C.6 BL5 pilot candidate pairs (A1, A2 ★, A3) with per-location deme structure, mother counts in the LEPA DB, and design rationale. One row per candidate.
- [`step29d_mating_pool_summary.tsv`](tables/Phase5/step29d_mating_pool_summary.tsv) — § A.5 per-location deme structure summary (n_mating_pools, largest/smallest/median deme size, demes below SI / 8-plant species-pool thresholds). Feeds Figure 1b and the BL5 pilot candidate file.
- [`step29c_sampling_comparison_per_location.tsv`](tables/Phase5/step29c_sampling_comparison_per_location.tsv) — per-location head-to-head of current (§ B.4.1) vs fragmentation-aware (§ B.4.2) M_recommended, with `n_components_50m`, `sum_M_paternal_only`, and `delta`.
- [`step29_sampling_per_event.tsv`](tables/Phase5/step29_sampling_per_event.tsv) — per-event allocation.
- [`step29_sampling_per_location.tsv`](tables/Phase5/step29_sampling_per_location.tsv) — per-location design table with connectivity-informed columns.
- [`step29_location_coverage_curves.tsv`](tables/Phase5/step29_location_coverage_curves.tsv) — analytical curves by location size.
- [`step29_field_team_sampling_recipe.tsv`](tables/Phase5/step29_field_team_sampling_recipe.tsv) — **FIELD TEAM per-germplasmID recipe** (one row per mother with the actionable seed count).
- [`step29_location_connectivity.tsv`](tables/Phase5/step29_location_connectivity.tsv) — within-location pollen connectivity at 10 / **25 (primary)** / 50 m.
- [`step30_A_prediction_prior_frequencies.tsv`](tables/Phase5/step30_A_prediction_prior_frequencies.tsv) — the P1 species-wide Fg prior.
- [`step30_A_prediction_location_diversity.tsv`](tables/Phase5/step30_A_prediction_location_diversity.tsv) — predicted SRK allele diversity per location (union across its 50 m components).
- [`step30_A_prediction_component_diversity.tsv`](tables/Phase5/step30_A_prediction_component_diversity.tsv) — per 50 m component diversity: pool size, allocated mothers + seeds, A_delivered, Fgs detected, and coverage.
- [`step30_A_prediction_location_pcompat.tsv`](tables/Phase5/step30_A_prediction_location_pcompat.tsv) — predicted random-mating pollen compatibility per location (size-weighted mean across its 50 m components).
- [`step30_A_prediction_component_pcompat.tsv`](tables/Phase5/step30_A_prediction_component_pcompat.tsv) — per 50 m component pollen compatibility — exposes struggling sub-components hidden inside a sustainable-mean location.
- [`step30_A_prediction_per_mother_fecundation.tsv`](tables/Phase5/step30_A_prediction_per_mother_fecundation.tsv) — species-wide compatibility reference distribution under the sporophytic tetraploid model.
- [`srk_fg_class.tsv`](tables/Phase5/srk_fg_class.tsv) — Fg → dominance class (I / II) mapping used by the sporophytic pollen compatibility model in § A.8. Provisional data-driven default; editable by hand as biology is refined.
- [`srk_zygosity_empirical.tsv`](tables/Phase5/srk_zygosity_empirical.tsv) — empirical LEPA distribution of distinct functional SRK identities per plant, from Canu-amplicon Step 23 (n = 367). Used to draw mother and father genotypes under the sporophytic finite-population model (§ A.8.3a, § A.8).
- [`step30_A_traffic_light_bands.tsv`](tables/Phase5/step30_A_traffic_light_bands.tsv) — recalibrated failed / struggling / sustainable band boundaries against the sporophytic species-mean pollen compatibility.
- [`step30_A_fragmentation_per_event.tsv`](tables/Phase5/step30_A_fragmentation_per_event.tsv) — per-event fragmentation index F_event = 1 − K_spatial_50m / 32 (purely spatial, no allele frequencies).
- [`step30_A_fragmentation_per_location.tsv`](tables/Phase5/step30_A_fragmentation_per_location.tsv) — per-location fragmentation index F_location = 1 − within-location connectivity at 50 m, with median event-scale F for the same location.

**Figures** (PNG + PDF):

- [`step28_coverage_curves.pdf`](figures/Phase5/step28_coverage_curves.pdf) / [`.png`](figures/Phase5/step28_coverage_curves.png) — per-mother coverage vs seeds and aggregation across mothers, one curve per event-size bin, with the per-mother cap. Justifies the tetraploid per-mother cap of 15 seeds/mother.
- [`step28d_pcompat_precision.pdf`](figures/Phase5/step28d_pcompat_precision.pdf) / [`.png`](figures/Phase5/step28d_pcompat_precision.png) — Part C justification (§ B.3.1) — per-mother observed P_compat SE vs seed count; precision plateau at ~15 seeds.
- [`step28d_matelim_power.pdf`](figures/Phase5/step28d_matelim_power.pdf) / [`.png`](figures/Phase5/step28d_matelim_power.png) — Part C justification (§ B.3.1) — § C.1 mate-limitation regression power vs seed count; ≥ 99 % power at n_seeds = 15 for medium/large β₁, marginal at small β₁.
- [`step29c_sampling_comparison.pdf`](figures/Phase5/step29c_sampling_comparison.pdf) / [`.png`](figures/Phase5/step29c_sampling_comparison.png) — fragmentation-aware sampling target (§ B.4.2) vs mothers already collected in the LEPA DB, BL-panelled.
- - [`step30_A_diversity_unbiased_vs_sampling.pdf`](figures/Phase5/step30_A_diversity_unbiased_vs_sampling.pdf) / [`.png`](figures/Phase5/step30_A_diversity_unbiased_vs_sampling.png) — three-panel per-location diversity figure: what the population holds (unbiased, `N_fertile_effective`), what our sampling detects (from actual LEPA DB seed counts), and the coverage fraction. Flags coverage against the 90 % target.
- [`step30_A_si_model_schematic.pdf`](figures/Phase5/step30_A_si_model_schematic.pdf) / [`.png`](figures/Phase5/step30_A_si_model_schematic.png) — three-panel pedagogical schematic of the sporophytic Class I / II model (§ A.7): dominance within a plant, worked example, compatibility rule by cross type.
- [`step30_A_prediction_fecundation.pdf`](figures/Phase5/step30_A_prediction_fecundation.pdf) / [`.png`](figures/Phase5/step30_A_prediction_fecundation.png) — BL-panelled compatibility per location, traffic-light bands.
- [`step30_A_diversity_vs_pcompat.pdf`](figures/Phase5/step30_A_diversity_vs_pcompat.pdf) / [`.png`](figures/Phase5/step30_A_diversity_vs_pcompat.png) — cross-plot: SRK diversity × pollen compatibility.
- [`step30_A_N_fertile_effective.pdf`](figures/Phase5/step30_A_N_fertile_effective.pdf) / [`.png`](figures/Phase5/step30_A_N_fertile_effective.png) — three-panel per-location view of the pivotal `N_fertile_effective` metric: raw census, connectivity-adjusted `N_fert_eff`, and the connectivity share. The single input that drives Figures 3, 4, and 5.

**Phase A validation (adult SRK genotypes from Phase 4 — Part C § C.0):**

**Primary Phase A validation — § C.0.a + § C.0.b** (clean-overlap EOs at Phase 5 location scale):

- [`step30_B_partC_clean_overlap_per_location.tsv`](tables/Phase5/step30_B_partC_clean_overlap_per_location.tsv) — § C.0.a one row per clean location (EO67, EO70, EO76 — 187 adults). Observed vs Phase 5 predicted pollen compatibility, observed vs Phase 5 predicted distinct allele count, no-drift upper bound at the adult sample size.
- [`step30_B_partC_clean_overlap_fg_frequencies.tsv`](tables/Phase5/step30_B_partC_clean_overlap_fg_frequencies.tsv) — long form, one row per (clean locationCode × allele) with observed f, species-wide P1 f, and present/absent flag.
- [`step30_B_partC_clean_overlap.pdf`](figures/Phase5/step30_B_partC_clean_overlap.pdf) / [`.png`](figures/Phase5/step30_B_partC_clean_overlap.png) — Figure 10: three-panel per-location figure (SRK diversity + pollen compatibility + per-allele drift residual).
- [`step30_B_partC_hypothesis_decomposition.tsv`](tables/Phase5/step30_B_partC_hypothesis_decomposition.tsv) — § C.0.b decomposition: diversity trigger (observed vs Phase 5 predicted allele counts), observed + three single-swap counterfactual pollen compatibility estimates, driver fractions for within-class spread, Class I/II balance, zygosity.
- [`step30_B_partC_hypothesis_decomposition.pdf`](figures/Phase5/step30_B_partC_hypothesis_decomposition.pdf) / [`.png`](figures/Phase5/step30_B_partC_hypothesis_decomposition.png) — Figure 11: grouped bar chart with gold-star markers on the dominant explanation per location.
- [`step30f_srk_diversity_radius_sweep.tsv`](tables/Phase5/step30f_srk_diversity_radius_sweep.tsv) — § C.0.c SRK-diversity radius sweep: one row per (EO, radius) with `n_components`, `component_sizes`, `pred_mean` / `pred_lo95` / `pred_hi95`, `observed`, `obs_inside_CI`.
- [`step30f_srk_diversity_radius_sweep.pdf`](figures/Phase5/step30f_srk_diversity_radius_sweep.pdf) / [`.png`](figures/Phase5/step30f_srk_diversity_radius_sweep.png) — Figure 13: predicted-vs-observed diversity across radii 10 → 150 m for EO67, EO70, EO76 (deme-size vs residual-drift test).

**Background diagnostic for split EOs — § C.0** (EO-scale pending Phase 5 remap):

- [`step30_C_pcompat_validation_at_eo.tsv`](tables/Phase5/step30_C_pcompat_validation_at_eo.tsv) — per-EO observed vs predicted pollen compatibility for EO18, EO25, EO27 (clean-overlap EOs are now in the two files above instead). Observed Class I mass + distinct-identity distribution per EO. No companion figure; see § C.0 for the three-row result table.

### Phase B — Post-genotyping (requires observed seed genotypes)

**Tables** (produced by `step30_srk_diversity_prediction_vs_observed.py --seed-genotypes … --mother-genotypes …`):

- `tables/Phase5/step30_B_comparison_location_diversity.tsv` — posterior vs predicted diversity per location.
- `tables/Phase5/step30_B_comparison_maternal_vs_paternal.tsv` — spectra test per location.
- `tables/Phase5/step30_B_mate_limitation_per_location.tsv` — Test 1 · location-level.
- `tables/Phase5/step30_B_mate_limitation_per_mother.tsv` — Test 1 · per-mother detail.
- `tables/Phase5/step30_B_mate_limitation_coefficients.tsv` — Test 1 · β₁, β₂, β₃ estimates with 95 % CI.
- `tables/Phase5/step30_B_si_escape_permutation.tsv` — out-of-scope permutation-test scaffold ([Appendix](#appendix--out-of-scope--future-work)).

**Figures**:

- `figures/Phase5/step30_B_comparison_diversity.pdf/png` — observed vs predicted diversity.
- `figures/Phase5/step30_B_mate_limitation.pdf/png` — Test 1 · location scatter.
- `figures/Phase5/step30_B_si_escape.pdf/png` — out-of-scope permutation-test scaffold ([Appendix](#appendix--out-of-scope--future-work)).

### Phase B — DEMO (synthetic pipeline-validation outputs, `--demo` mode)

Same filenames as Phase B above with `_B_` → `_B_DEMO_`. Figures also
carry a "DEMO — synthetic data" title suffix and a diagonal DEMO
watermark so a demo file can never be mistaken for a real result.

---

## Appendix — Out of scope / future work

This appendix collects analyses that the Step 30 code can produce
but that **are not part of the Phase 5 scientific analysis**. They
are retained as pipeline-validation scaffolding and as a starting
point for future work under a separate study, not as a claim about
what the current data can address.

### SI-escape rate test (was § C.2 in earlier drafts)

**Why this is out of scope for Phase 5.** The Phase 5 framework tests
**mate limitation** (does drift-driven Fg loss reduce per-mother seed
set?) and **fragmentation** (do disconnected 50 m demes
compound the drift effect?). Detecting locations where the SI
machinery has broken down — self-compatibility escape — is a
separate scientific question with different data requirements
(large per-mother seed lots + genotyped mother tissue + independent
phenotypic confirmation of self-seed set). Self-compatibility calls
at the *individual* level are already produced by the Canu-amplicon
pipeline (Phase 4 Step 22b) and enter Phase 5 as prior information
via the empirical zygosity distribution (§ A.7.3a), not as an
experimental target. The permutation-test scaffold below is retained
so the pipeline can be validated end-to-end under `--demo`; it is
also a natural starting point should a separate SI-escape study be
funded.

**Model.** Under **strict SI**, the paternal SRK allele in any seed
of mother $(a, b)$ cannot equal $a$ or $b$. The rate of self-matching
paternal alleles per location is therefore expected to be

$$\pi_{\text{self-match}}^{H_0} \;=\; 0$$

Under **partial SI** or SI breakdown, some fraction of seeds carry a
paternal Fg matching one of the mother's Fgs:

$$\hat{\pi}_{\text{self-match}}(\ell) \;=\; \frac{\bigl|\{\text{seeds at loc.}\;\ell : \text{paternal Fg} \in (a_m, b_m)\}\bigr|}{\bigl|\text{seeds at loc.}\;\ell\bigr|}$$

The **test** is a permutation of the strict-SI null: for each location,
permute paternal alleles across seeds while preserving the marginal Fg
frequency vector; the p-value is the fraction of permutations reaching
$\hat{\pi}_{\text{self-match}}$ at least as extreme as observed.

Locations with $\hat{\pi}_{\text{self-match}} > 0$ significantly would
be **candidate SI-escape sites** — a validated hypothesis for
follow-up phenotyping in Genetic-Rescue-DB, under a separate study.

**Output columns** (per location, produced by the DEMO pipeline):
`n_seeds_scored`, `n_self_matching`, `pi_self_match`,
`p_permutation`, `q_bh_fdr`.

**DEMO figure.**

![Figure A1](figures/Phase5/step30_B_DEMO_si_escape.png)

**Figure A1 (DEMO).** Self-incompatibility escape preview. One horizontal bar per LEPA location, sorted by observed rate. X = observed rate of pollen alleles matching the mother's own SRK alleles (= self-incompatibility escape rate). Red bars = locations that reject the strict-SI null at 5 % false-discovery rate; grey bars = consistent with strict SI. In this DEMO the simulator baked in 8 % SI escape rate, so most locations show detectable escape. **This figure exists only to validate the pipeline end-to-end; the SI-escape analysis is out of scope for Phase 5 (see the § A.1 preamble above).** Source: `step30_srk_diversity_prediction_vs_observed.py --demo`.

---

## References

The population-genetics scaffolding for the operational-deme
framework (§ A.5), the fragmentation → drift → mate-limitation
causal chain (§ Scientific goals), and the sporophytic
self-incompatibility biology (§ A.7). BibTeX entries live in
[`Phase5_references.bib`](Phase5_references.bib).

- **Aguilar, R., Quesada, M., Ashworth, L., Herrerías-Diego, Y. & Lobo, J. (2008).** Genetic consequences of habitat fragmentation in plant populations: susceptible signals in plant traits and methodological approaches. *Molecular Ecology* 17, 5177–5188. — Fragmentation signatures in plant populations; motivation for the fragmentation → drift → mate-limitation causal chain tested here.
- **Castric, V. & Vekemans, X. (2004).** Plant self-incompatibility in natural populations: a critical assessment of recent theoretical and empirical advances. *Molecular Ecology* 13, 2873–2889. — Review of SI population-genetic theory; context for the sporophytic Class I / Class II model in § A.7.
- **Freckleton, R.P. & Watkinson, A.R. (2002).** Large-scale spatial dynamics of plants: metapopulations, regional ensembles and patchy populations. *Journal of Ecology* 90, 419–434. — Framework for distinguishing patch dynamics from true metapopulations in plants; the formal justification for separating the stable "patch" from the within-year above-ground local population used in § A.5.0.
- **Hanski, I. (1998).** Metapopulation dynamics. *Nature* 396, 41–49. — Modern synthesis of metapopulation theory; motivates the population-as-stable-patch convention used in § A.5.0.
- **Hardy, O.J. & Vekemans, X. (1999).** Isolation by distance in a continuous population: reconciliation between spatial autocorrelation analysis and population genetics models. *Heredity* 83, 145–154. — Estimator for Wright's neighbourhood from fine-scale *F*<sub>ST</sub> ~ distance; the statistical route to the *N*<sub>b</sub> we do **not** currently estimate.
- **Levins, R. (1969).** Some demographic and genetic consequences of environmental heterogeneity for biological control. *Bulletin of the Entomological Society of America* 15, 237–240. — Original formulation of the metapopulation concept; theoretical foundation for the population level introduced in § A.5.0.
- **Honnay, O. & Jacquemyn, H. (2007).** Susceptibility of common and rare plant species to the genetic consequences of habitat fragmentation. *Conservation Biology* 21, 823–831. — Rare-species context for LEPA; drift effects on SI loci are a known risk category.
- **Lawrence, M.J. (2000).** Population genetics of the homomorphic self-incompatibility polymorphisms in flowering plants. *Annals of Botany* 85 (Suppl. A), 221–226. — Allele-count theory for *S*-locus systems; background for the 32-Fg species-wide pool.
- **Levin, D.A. & Kerster, H.W. (1974).** Gene flow in seed plants. *Evolutionary Biology* 7, 139–220. — Pollen-mediated gene flow in flowering plants; the operational-deme radius is in the empirical range reported here for small-bee pollinators.
- **Schierup, M.H., Vekemans, X. & Christiansen, F.B. (1998).** Allelic genealogies in sporophytic self-incompatibility systems in plants. *Genetics* 150, 1187–1198. — Coalescent theory for sporophytic SI; drift expectations on *S*-locus allele frequencies with dominance.
- **Vekemans, X. & Hardy, O.J. (2004).** New insights from fine-scale spatial genetic structure analyses in plant populations. *Molecular Ecology* 13, 921–935. — Methods review for pollen- and seed-dispersal-driven population structure; the framework we would use to validate the 50 m operational-deme choice with genetic data.
- **Wright, S. (1943).** Isolation by distance. *Genetics* 28, 114–138. — The genetic-neighbourhood concept (*N*<sub>b</sub> = 4πσ²δ); theoretical scaffolding for why a pollinator-flight scale is biologically meaningful.
- **Wright, S. (1946).** Isolation by distance under diverse systems of mating. *Genetics* 31, 39–59. — Extension of the 1943 neighbourhood concept to mixed mating systems; relevant background for a tetraploid sporophytic SI species.
