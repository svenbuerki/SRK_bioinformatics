# Phase 5 — SRK sampling & prediction framework (Steps 28–30)

## Contents

The document is organised in **three parts** that follow the causal
order of the study — we first build the model and generate location-
level predictions, then derive the sampling protocol from those
predictions, then test the predictions with real seed data:

- [Scientific goals](#scientific-goals) — the four hypotheses this framework tests
- **Part A** — [Model and predictions](#part-a--model-and-predictions) · data scope, species-wide P1 prior, finite-population model, and the Phase A predictions of SRK diversity and pollen compatibility per location (Step 30 Phase A outputs, Figures 4–6). *This is what we expect each location to look like — before we ever open a seed lot.*
- **Part B** — [Sampling protocol derived from the predictions](#part-b--sampling-protocol-derived-from-the-predictions) · within-location pollen connectivity, per-mother seed count (Step 28), per-location mother count with private-allele floor (Step 29), and two-year design. *This is what the field team must do to test the Part A predictions.*
- **Part C** — [Testing predictions with observed SRK data (Phase B)](#part-c--testing-predictions-with-observed-srk-data-phase-b) · preliminary EO-level validation with Phase 4 adult genotypes (§ C.0), mate-limitation regression, SI-escape rate test, script behaviour, data-generation pipeline, and the BL4 pilot.
- [Output map](#output-map--quick-reference-grouped-by-phase) — filenames organised by phase.

**Naming.** Parts **A / B / C** refer to *sections of this document*.
Output filenames use `step30_A_*` for Phase A prediction artefacts (built
without seed data) and `step30_B_*` for Phase B artefacts (built with
observed seed data); `step30_B_DEMO_*` marks synthetic Phase B for
pipeline validation. Steps 28–29 outputs are all Phase A by construction
and keep their existing `step28_*` / `step29_*` names.

## Scientific goals

The purpose of this framework is not diversity estimation for its own sake.
It is a location-level SRK sampling and inference pipeline that lets us test
four connected hypotheses about the reproductive fate of small, isolated
LEPA populations:

1. **Mate-limitation test.** Do locations with fewer compatible pollen
   donors — driven by small mate pool, skewed Fg frequencies, or both —
   show reduced per-mother seed set? Formally: does per-mother
   `germplasmQuantityEstimate` decline with predicted per-mother
   pollen compatibility (sporophytic Class I / II tetraploid model,
   see § A.7 for the biology and § A.8 for the maths)?
2. **SI-escape test.** Under strict SI, no seed can carry a paternal Fg
   matching either of its mother's Fgs. A location where seeds with
   maternal-matching paternal Fgs appear above the strict-SI null is a
   candidate for partial-SI transition (a breakdown of the SI machinery).
3. **Fragmentation × drift decomposition.** Habitat fragmentation depresses
   K (the pollen-donor Fg pool a mother can access); genetic drift skews
   Fg frequencies at small isolated locations. Both reduce pollen
   compatibility but through different channels; the framework
   predicts them on two independent axes in Part A — **drift** via
   predicted SRK diversity and pollen compatibility (§ A.7–A.8) and **fragmentation**
   via the purely spatial pollen-flow index at event and location scales
   (§ A.5) — and Phase C's mate-limitation regression (§ C.1) then tests
   both simultaneously as two independent coefficients (β₁ = drift,
   β₂ = fragmentation).
4. **Phenotype cross-validation (deferred).** SRK-based predictions can be
   cross-checked against per-site ISI / fruit set in the Genetic-Rescue-DB
   repository. This adds independent lines of evidence but does not shape
   the sampling design and is not modelled here in v1.

The framework has two design consequences: **(1)** the sampling protocol
must scale with each mother's real mate-availability context (Steps 28–29),
and **(2)** the inference layer must produce testable predictions under
random-mating and strict-SI nulls that can be compared to observed seed
genotypes (Step 30). Part A publishes the predictions; Part B derives the
sampling protocol needed to test them; Part C runs the tests.

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
| **B — Post-genotyping** (after SRK data are back) | Everything above + observed seed genotypes (from Phase A's sampling) | Step 30 in `comparison` mode + Tests 1 & 2 | Observed vs predicted SRK diversity, mate-limitation regression, SI-escape rate test — all only for the locations that have observed data |

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
  wild-population mate-limitation and SI-escape analysis, so the filter
  is baked in.

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
output family:

| Step | Question | Kind of output |
|---|---|---|
| **Step 28** | How many seeds per mother? | Sampling design (per-mother) |
| **Step 29** | How many mothers per event, how many events per location? | Sampling design (per-event / per-location) |
| **Step 30** | Given the species-wide prior, what SRK diversity and pollen compatibility should we predict per location — and how do observed seed genotypes compare? | Prediction + comparison |

### A.3 The two-generation trick and species-wide prior

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

Once Part B's sampling protocol is executed, we will have a batch
of seed DNA per mother. Each seed is **tetraploid** like its parents —
2 maternal + 2 paternal SRK alleles at each SRK locus — so one seed
lot per mother recovers, in a single extraction batch, two independent
samples:

| Generation | What we recover | How |
|---|---|---|
| **Parental (G0)** | The mother's 4 SRK allele copies (her full genotype) | Invariant / 50 %-frequency alleles across her sibs |
| **Filial (G1) — as read through paternity** | The pollen SRK allele pool she was exposed to, 2 draws per seed | Variable alleles across her sibs |

Aggregating across mothers of a location gives us **two independent
estimates of the location-level SRK allele frequency spectrum**: the
*maternal* one (who is standing there) and the *paternal* one (who is
actually contributing pollen). Under random mating with panmictic pollen
dispersal, the two spectra are indistinguishable. Departures flag biased
contribution, cryptic SI filtering, or immigrant pollen — themselves useful
signal.

But we do not want to wait until the seed data are in to know what we
expect. The Canu-amplicon preliminary study already characterised **32
functional SRK allele groups (Fgs) in 263 individuals** — a strong empirical
**prior**. We can use that prior *now* to publish predicted diversity and
predicted fecundation failure per location, with credible intervals. When the
seed data land, prediction and observation are compared side by side, and
the locations where they disagree are the actionable ones.

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

Every Part A prediction (A.5–A.8) comes from the same **finite-population
simulation model**, applied per location, with the **50 m connectivity
radius** as the biological scope of pollen movement. Four steps:

1. **Simulate the local SRK pool (drift signature).**
   Each location's mating population has `N_fertile × 50 m-connectivity`
   plants; those contribute `2 × N_fertile` allele copies. We **draw
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
× 15 seeds each (tetraploid Rule 2) = **A = M × 34 allele draws** from
the local pool (4 alleles per mother from her own tetraploid
genotype + 2 paternal alleles × 15 seeds = 30 paternal draws). We
report coverage in two flavours:

- **Species-wide coverage** — fraction of the 32 P1 alleles detected.
  Biased low for small locations (drift already removed most).
- **Local coverage** — fraction of the alleles *actually at the
  location* that we detect. The biologically honest metric; ~100 % at
  every location under the current design (see A.5 table).

**One-line take-home:** SRK diversity is what drift has left, pollen
compatibility is how well a random mother matches the neighbours she
can reach at 50 m, and both are simulated per location under the
same finite-population draw from the species-wide prior.

**Prediction outputs.** Kept in `tables/Phase5/` prefixed
`step30_A_prediction_*` and figures under `figures/Phase5/step30_A_*`.
They can be produced today, without any seed data at all. The
corresponding *comparison* outputs (Phase B, once seed genotypes exist)
carry a different prefix (`step30_B_*`) so prediction and comparison
cannot be mistaken for one another.

### A.5 Fragmentation of pollen flow — event and location scales


**Pollinator-radius choice.** Every fragmentation and downstream prediction in Phase 5 depends on the assumed pollen-flight radius. A dedicated sensitivity sweep across 10, 25, 50, 75, 100, 150, 200 m ([Figure 1](#fig-1)) demonstrates that **50 m is the biologically sound primary radius**: it sits within the halictid / small-bee foraging literature range, captures 43 % of the sampling-cost reduction available on the radius curve, and retains meaningful fragmentation variation across BLs (median connectivity 0.81, not yet saturated at 1.00 like at 75 m+). Every downstream analysis in this doc uses 50 m as the primary radius; 10 m and 25 m are always computed and stored as sensitivity checks.

<a id="fig-1"></a>
![Figure 1: Full pollinator-radius sensitivity sweep across 10, 25, 50, 75, 100, 150, 200 m — the biological justification for adopting **50 m as the primary pollinator radius**. **Panel A** (landscape connectivity) shows the fraction of adults in a multi-event pollen-flow component; the median across locations rises from 0.00 at 10 m to 0.45 at 25 m to 0.81 at 50 m and plateaus at 1.00 by 75 m. **Panel B** (§ B.4.2 sampling cost) shows the total mothers required across all 39 locations; it drops steeply from 887 at 10 m to 505 at 50 m (43 % reduction) then flattens to 328 at 200 m. **Panel C** (§ A.8 pollen-compatibility prediction, empirical zygosity) is essentially flat across all radii at ~0.68 — because ~66 % of LEPA plants are single-identity homozygotes, radius-driven changes in effective N barely move P_compat. **Panel D** shows that all 39 locations stay in the "sustainable" band regardless of radius. **The two curves that DO change with radius are connectivity (Panel A) and sampling cost (Panel B); both stabilise around 75–100 m. 50 m sits within the halictid / small-bee foraging literature range (10–100 m), captures 43 % of the sampling-cost reduction available on that curve, and retains meaningful fragmentation variation across BLs (0.81 median, not yet saturated at 1.00 like at 75 m+). Source: `sensitivity_pollinator_radius.py`. Data: [`step30_A_radius_sensitivity_summary.tsv`](tables/Phase5/step30_A_radius_sensitivity_summary.tsv), [`step30_A_radius_sensitivity_per_location.tsv`](tables/Phase5/step30_A_radius_sensitivity_per_location.tsv).](figures/Phase5/step30_A_radius_sensitivity.png)

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
§ B.3 for the coupon-collector formulas, viewed here as a
**fragmentation predictor** at the mother scale.
[Figure 2](#fig-2) shows the **distribution of K^(25m) across events
within each location**, panelled by BL. Plotting the raw distribution
rather than a single summary keeps the fine-scale process visible —
which is exactly the level at which mate limitation acts, since a
mother at an N_reachable = 1 event and a mother at an N_reachable ≈
100 event of the same locationID face very different pollen contexts.
In the 2025 field, **439 / 704 events (62 %)** have at least one other
event within 50 m; the number of reachable donor plants per event
ranges from **0** (isolated singletons) to **> 250** (dense BL3, BL4
clusters), and under tetraploid this corresponds to K^(25m) allele
copies from **0** to **> 1 000**.

**Location scale — the site-level view.** At the location scale
fragmentation is the fraction of adults NOT in a multi-event
pollen-flow component,

$$F_{\text{location}} \;=\; 1 - \text{connected-share at 50 m}$$

using the connected-share metric from § B.2. `F_location = 0` when
every adult exchanges pollen with at least one other event; `F_location = 1`
when every adult is an isolated island. In the 2025 field, **21 / 39
locations (54 %) sit at F_location ≥ 0.5** — the majority of LEPA
sites are structurally fragmented at pollinator scale. The
location-scale summary is shown in Figure 7 (§ B.2); [Figure 2](#fig-2)
here shows the event-scale distribution underneath it so the two views
are complementary rather than duplicative.

**Reading the boxplot.** Three patterns dominate [Figure 2](#fig-2):

- **Well-connected locations** — a tight box entirely to the right of
  the **N = 8 species-pool floor**. Every event reaches ≥ 32 tetraploid
  allele copies; mothers experience a uniform pollen environment.
  Examples: EO30-1, EO27-3, EO29, EO70.
- **Mixed-connectivity locations** — a wide box spanning **N = 1 → ≥
  100**. Some events are isolated singletons while others sit in
  connected clusters, and mothers of the same locationID face radically
  different mate contexts. Examples: EO27-1, EO18-7, EO26-3, EO8
  (several sub-locations). These are the locations where goal 3's β₂
  (fragmentation) coefficient in Phase C will have the most
  within-location leverage.
- **Uniformly isolated locations** — a tight box at **N ≤ 1** (or a
  single dot at N = 0). Every event is a lone plant; under strict SI
  no seed set is possible. Examples: EO24, EO24-1, EO24-2, EO24-7.
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

<a id="fig-2"></a>
![Figure 2: Event-scale distribution of the pollen-donor **plant count** reachable within 50 m per LEPA location, panelled by Bottleneck Lineage in the standard BL_ORDER (BL4, BL5, BL3, BL1, BL2). X-axis = `N_reachable_50m` = Σ N_fertile in events within 50 m − 1 (donor plants, not allele copies — plotted as plants to avoid conflation with the 32-Fg species allele-class count). Vertical **red dotted line at N = 1**: single-plant SI floor; every event at or below this line has zero reachable pollen donors and no seed set is possible under strict self-incompatibility. Vertical **grey dashed line at N = 8 plants**: coupon-collector floor for the 32-Fg species pool under **tetraploid LEPA** (4 alleles per plant × 8 plants = 32 allele copies). A mother reaching ≥ 8 donor plants has enough allele copies for the species pool to be *physically* reachable in principle; whether drift preserved the diversity is A.5's question. Log₂ x-axis so the 1 → 8 range (below the species floor) and the 8 → 500 range are both legible. Row labels give the location code, number of events, and total adult census. The figure exposes the fine-scale spatial process that a location-scale metric collapses: **a tight box far right of 8** (EO30-1, EO29, EO70) = every event well-connected, uniform pollen environment; **a wide box spanning 1 → hundreds** (EO27-1, EO18-7, EO26-3, EO8 groups) = the location holds a mix of isolated singletons and connected clusters, so mothers at different events face very different mate-availability contexts under the same locationID; **a tight box at N ≤ 1** (EO24 group, EO24-1, EO24-2, EO24-7) = every event is a lone plant. The metric uses only event coordinates, N_fertile, and the 50 m primary pollinator radius — no allele frequencies, no priors. Source: `step30b_fragmentation_index.py`. Data: [`step30_A_fragmentation_per_event.tsv`](tables/Phase5/step30_A_fragmentation_per_event.tsv), [`step30_A_fragmentation_per_location.tsv`](tables/Phase5/step30_A_fragmentation_per_location.tsv).](figures/Phase5/step30_A_fragmentation_index.png)

### A.6 Predicted SRK diversity per location

**What we did.** Under the species-wide 32-Fg prior, we predicted how
many distinct SRK alleles Phase B seed genotyping should recover at
each location, using the permit-realistic count of mothers actually
sampled (55 at EO61, 62 at EO76, down to 1–4 at the smallest
slickspots). Predictions are panelled by Bottleneck Lineage
(BL4 → BL5 → BL3 → BL1 → BL2).

**Result.** Under tetraploid sampling (PLOIDY × M = 4M allele draws per
location — see § A.3 and § B.3), large, well-buffered slickspots in BL3
(EO76) and BL1 (EO61) are predicted to recover 24–25 distinct alleles
out of the species-wide 32. BL5 is highly bimodal: the well-sampled
locations (EO32, EO18-7 group) reach 20–22 alleles, while the BL5 tail
(EO24, EO24-1, EO24-2, EO24-7) can recover only 3–8 alleles from their
tiny local populations.

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
computed for every location (see the `predicted_local_pool_size_mean`
and `predicted_local_coverage_mean` columns in
[`step30_A_prediction_location_diversity.tsv`](tables/Phase5/step30_A_prediction_location_diversity.tsv)):

| Location | Effective N (50 m) | Local pool (alleles present) | Local coverage at 15 seeds × M (tetraploid Rule 2) | Species coverage (of 32) |
|---|---|---|---|---|
| EO24-2 (1 plant) | 1 | ~3 | **~100 %** | ~11 % |
| EO67 (small pilot) | 6 | ~8 | **~100 %** | ~28 % |
| EO27-1 (large pilot) | 116 | ~27 | **~98 %** | ~68 % |
| EO32 (well-sampled BL5) | 327 | ~31 | ~95 % | ~70 % |
| EO76 (largest BL3) | 416 | ~32 | ~98 % | ~77 % |

**Every LEPA location — including the tiny BL5 slickspots — is
essentially fully characterised at the local level (~95–100 %).** The
species-wide coverage number remains useful as a cataloguing metric,
but for Phase B mate-limitation and SI-escape tests the local coverage
is what matters: we are testing reproductive dynamics on the alleles
that are physically present, not attempting a species-wide inventory.

<a id="fig-3"></a>
![Figure 3: Predicted number of distinct SRK alleles detected per LEPA location under the P1 species-wide prior (32 Fgs from the Canu-amplicon preliminary study). Horizontal layout, one dot per location, panelled by Bottleneck Lineage in the standard BL_ORDER (BL4, BL5, BL3, BL1, BL2). Dot size ∝ √M (permit-realistic count of mothers with seed records in the LEPA DB); error bars = 95 % credible interval from Dirichlet posterior draws. Y-tick labels give `locationCode (n = mothers sampled)`. Vertical dashed line marks the species-wide SRK allele ceiling of 32. **BL3 (EO76, EO38) and BL1 (EO61) predict ~24–25 distinct alleles under tetraploid sampling (PLOIDY · M = 4M allele draws per location); the BL5 tail (EO24 group) predicts 3–8.** Source: `step30_srk_diversity_prediction_vs_observed.py`.](figures/Phase5/step30_A_prediction_diversity.png)

### A.7 Sporophytic self-incompatibility with Class I / Class II dominance

This section states the biological SI model; § A.8 is the
mathematical implementation and § C.1 the Phase C test that uses
its predictions. [Figure 4](#fig-4) below is a purely pedagogical
schematic that summarises the whole section in three panels — the
dominance rule inside one plant (§ A.8.3), the worked example of
one mother vs three candidate fathers (§ A.8.5), and the
compatibility rule by cross type (§ A.8.4).

<a id="fig-4"></a>
![Figure 4 — Sporophytic SI with Class I / Class II dominance in tetraploid LEPA (Phase 5 § A.7). **Panel A — dominance within one plant.** Case A (plant with ≥ 1 Class I allele): only its Class I alleles are expressed on both pollen and stigma; Class II alleles are silent (shown faded). Case B (plant with only Class II alleles): all four Class II alleles are expressed co-dominantly. **Panel B — between-plant recognition, worked example.** Mother M carries `{FG001, FG002, FG024, FG031}`; her Case-A expressed set is {FG001, FG002}. Three candidate fathers: F1 shares FG001 with M → rejected; F2 is all-Class-II so between-class → always compatible; F3 shares FG002 with M → rejected. **Panel C — compatibility rule by cross type.** Class I × Class I: compatible if their expressed Class I alleles differ (Class II silent on both sides — sharing them is irrelevant). Class I × Class II: always compatible by construction (disjoint expressed classes). Class II × Class II: all four alleles expressed on both sides, compatible only if none are shared. Source: `si_model_schematic.py`.](figures/Phase5/step30_A_si_model_schematic.png)

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

- **Class I** — older lineages, deep sub-allele polymorphism,
  dominant within a plant.
- **Class II** — younger lineages, tight single sub-allele, recessive
  or co-dominant among themselves.

The same pattern is visible in LEPA's P1 empirical prior:

| Class | Fgs | Sub-alleles per Fg | Total P1 frequency |
|---|---|---|---|
| **Class I candidates** | 6 (FG001–FG006) | 2 – 10 | ~65 % |
| **Class II candidates** | 26 (FG007–FG032 except FG024) | 1 | ~35 % |
| *Anomaly* | FG024 | 1 | 18 % — flagged `REVIEW` |

The mapping is stored as a first-class TSV
([`srk_fg_class.tsv`](tables/Phase5/srk_fg_class.tsv)); rerunning
`step30_srk_diversity_prediction_vs_observed.py` regenerates every
downstream number from the current version.

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
0.69 (empirical LEPA zygosity)** — a huge biological signal that
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

**What we did.** For each location we simulated the local mating
pool by drawing 4 × N_fertile alleles from the species-wide prior
(tetraploid, see § A.3), paired them into N_fertile tetraploid
plants, and evaluated each sampled mother's random-mating
compatibility under the **sporophytic Class I / Class II dominance
model** implemented in [`srk_si_model.py`](srk_si_model.py). The
model captures three biological facts that the diploid gametophytic
approximation used in Part 1 could not:

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
pollen compatibility is **~0.69** — very close to the Part-1 diploid
gametophytic estimate of 0.63, but reached via a completely different
mechanism. Under a naive independent-tetraploid-draw assumption the
species mean is only ~0.16, but that model over-counts heterozygosity
(it predicts ~5 % homozygotes; LEPA shows ~66 %). Once each mother
and father is drawn under the empirical zygosity distribution, ~66 %
of mothers are single-identity Case A homozygotes with a very small
expressed set → high per-mother pollen compatibility, ~32 % have
two identities → intermediate, ~2 % have three → lower. The
weighted mean lands at 0.69. Traffic-light bands are recalibrated
against this new species mean: **failed < 0.231, struggling
0.231–0.462, sustainable ≥ 0.462** (1/3 and 2/3 of species mean).
Every LEPA location's mean sits close to 0.69; **small-slickspot
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

<a id="fig-5"></a>
![Figure 5: Predicted per-mother pollen compatibility under the sporophytic tetraploid Class I / Class II model with **empirical LEPA zygosity** (§ A.8.3a). One dot per LEPA location, panelled by Bottleneck Lineage. Traffic-light background bands mark **failed** (pollen compatibility < 0.231, red), **struggling** (0.231–0.462, orange) and **sustainable** (≥ 0.462, green) — recalibrated against the sporophytic + empirical-zygosity species-mean of **0.693** (green dotted line). Dot position = mean predicted pollen compatibility from a Monte-Carlo finite-population simulation: 4 × N_fertile local alleles drawn from P1 (N_fertile scaled by within-50 m connectivity), M mothers drawn from that pool with the empirical LEPA zygosity distribution (66 % single-identity homozygotes, 32 % 2-distinct, 2 % 3-distinct — see [`srk_zygosity_empirical.tsv`](tables/Phase5/srk_zygosity_empirical.tsv)), and each mother's pollen compatibility computed against 300 candidate fathers drawn the same way. Error bars = 95 % credible interval across simulation replicates; dot size ∝ √M (mothers with seed records in DB). **Every LEPA location's mean sits close to the species mean because 66 % of mothers express only one SRK identity, which minimises their p(M) footprint under the § A.8.4 recognition rule; BL5 tiny slickspots retain wide CI reflecting founder-effect variance in class + zygosity composition.** Source: `step30_srk_diversity_prediction_vs_observed.py`. Class assignments: [`srk_fg_class.tsv`](tables/Phase5/srk_fg_class.tsv). Empirical zygosity: [`srk_zygosity_empirical.tsv`](tables/Phase5/srk_zygosity_empirical.tsv). Bands: [`step30_A_traffic_light_bands.tsv`](tables/Phase5/step30_A_traffic_light_bands.tsv).](figures/Phase5/step30_A_prediction_fecundation.png)

### A.9 Cross-plot: SRK diversity vs pollen compatibility

**What we did.** We plotted per-location predicted SRK allele
diversity (x, from § A.6) against predicted sporophytic pollen
compatibility (y, from § A.8), with locations coloured by BL and the
sporophytic + empirical-zygosity species-mean pollen compatibility
(0.693) drawn as a horizontal reference line.

**Result.** Under the sporophytic Class I / II model with the
empirical LEPA zygosity distribution (§ A.8.3a), **most locations
sit near the species mean 0.69**. The BL5 tail (EO24 group) is
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
![Figure 6: Predicted SRK allele diversity (x) vs predicted sporophytic pollen compatibility (y) per LEPA location, coloured by Bottleneck Lineage (Set1 palette: BL1 purple, BL2 blue, BL3 red, BL4 orange, BL5 green), under the sporophytic + empirical-zygosity model (§ A.8.3a). Error bars on both axes come from the same Dirichlet posterior draws that produced the two single-quantity Phase A figures. Dot size ∝ √M (mothers with seed records in DB). **Green dotted horizontal line = sporophytic + empirical-zygosity species-mean pollen compatibility (0.693)** — the reference every location can be read against under the § A.7 model. Most locations sit right on the species mean because ~66 % of LEPA plants express only one SRK identity, giving them a small p(M) footprint and therefore high pollen compatibility. **BL5 tail (EO24, EO24-1, EO24-2)** sits at low predicted diversity (x ≈ 2–5) but stays close to the species mean on the y-axis with much *wider* credible intervals — the honest founder-effect signal that a small slickspot can drop into "struggling" (y < 0.462) if class + zygosity composition breaks unfavourably at that specific site. This figure is the single-figure summary of goal 3's fragmentation × drift decomposition: **sporophytic + empirical zygosity flatten the mean, but small locations still carry drift-driven risk in their credible intervals**. Source: `step30_srk_diversity_prediction_vs_observed.py`.](figures/Phase5/step30_A_diversity_vs_pcompat.png)

---

## Part B — Sampling protocol derived from the predictions

Part A tells us *what SRK diversity each location should carry*. Part B
derives *what the field team must sample to test those predictions*.
The sampling design is not an arbitrary field protocol: it is engineered
to characterise, at each location, the SRK alleles that Part A's
finite-population model says are physically present, at a stated
statistical guarantee.

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

- Fewer connected adults → smaller effective mating pool → more drift →
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
  into the mating pool.
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
Phase B mate-limitation regression. See [**Figure 7**](#fig-7)
(per-location bars at 50 m) and [**Figure 8**](#fig-8)
(radius-sensitivity aggregate) below.

<a id="fig-7"></a>
![Figure 7: Within-location pollen connectivity across the 39 LEPA locations at the **50 m primary pollen-flight radius**. One bar per location, panelled by Bottleneck Lineage (BL4 orange, BL5 green, BL3 red, BL1 purple, BL2 blue — Set1 palette shared with the LEPA_EO_spatial_clustering project). Bar length = fraction of the location's adults that sit in a connected component containing more than one event. Vertical guides: 50 % (orange dotted) and 90 % (green dotted) thresholds. Row labels list the location code, number of events, and total adult census. **At 50 m, 18 / 39 locations reach 90 % within-location connectivity; 11 / 39 sit below 50 %.** Sensitivity views at 10 m (conservative small-bee patch) and 25 m (short-flight) are stored alongside as `step29_location_connectivity_10m.png/pdf` and `step29_location_connectivity_25m.png/pdf`. Source: `step29b_location_connectivity.py`. Data: [`step29_location_connectivity.tsv`](tables/Phase5/step29_location_connectivity.tsv).](figures/Phase5/step29_location_connectivity.png)

<a id="fig-8"></a>
![Figure 8: Radius sensitivity of within-location pollen connectivity across all 39 LEPA locations at the three sampling-recipe-relevant radii (10 / 25 / 50 m). Orange bars = fraction of locations reaching ≥ 50 % adults connected at each radius; green bars = fraction reaching ≥ 90 %. At **10 m** (conservative small-bee patch) **0 / 39** locations are fully connected — the framework would flag every LEPA location as fragmented, an over-strong claim. At **25 m** (short-flight sensitivity) **8 / 39 (21 %)** reach 90 % connectivity. **At 50 m (primary), 18 / 39 (46 %) reach 90 % connectivity and 28 / 39 (72 %) reach 50 %** — the "sweet spot" between over-fragmentation and full saturation, and the choice justified by the § A.5 pollinator-radius sweep ([Figure 1](#fig-1)). Source: `step29b_location_connectivity.py`.](figures/Phase5/step29_location_connectivity_radius_sensitivity.png)

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

**Methodology.** We treat mating as a **coupon-collector problem on
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
>   coupon-collector formulas below.
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
![Figure 9: two-panel SRK allele detection under **tetraploid LEPA** (4 SRK copies per plant; 2 paternal alleles per seed) on the absolute-allele scale, with the local pool capped at the species-wide ceiling of 32 Fgs. Shaded bands = 95 % simulation CI. **Panel A — per mother.** x = seeds genotyped, one curve per event-size bin (each seed contributes 2 paternal allele draws); the vertical red line at **15** marks the Rule 2 operational cap under tetraploid (never ask any single mother for more than 15 seeds). The 15-seed dots show what each event size **delivers per mother**: **4.0 of 4** (1–2 plants), **11.1 of 12** (3–5), **18.6 of 28** (6–10), **19.7 of 32** (11–20), **19.7 of 32** (21–50), **19.7 of 32** (>50). One mother's 15 seeds cannot saturate a 32-allele pool — this is the coupon-collector limit for a single sampler, not undersampling. **Panel B — aggregation across mothers at a location.** x = number of mothers sampled at the location (15 seeds each = **34 allele draws per mother: 4 maternal + 2·15 = 30 paternal**). At the 5-mother benchmark (green dotted line, 75 cumulative seeds): **4.0 of 4**, **12.0 of 12**, **27.9 of 28**, **31.9 of 32**, **31.9 of 32**, **31.9 of 32** — every event size reaches its local ceiling. Because Panel B is a coupon-collector simulation continuous in M, any real location can read off its own coverage by locating (its event-size bin, its actual M) on the correct curve. The "gap at large events" in Panel A closes cleanly at the location scale — see § B.4. Source: `step28_seed_sampling_per_mother.py`. Data: [`step28_coverage_curves_by_Nfertile.tsv`](tables/Phase5/step28_coverage_curves_by_Nfertile.tsv), [`step28_aggregation_curves_by_Nfertile.tsv`](tables/Phase5/step28_aggregation_curves_by_Nfertile.tsv).](figures/Phase5/step28_coverage_curves.png)

### B.4 Per-location mother count — sampling design (Step 29)

> **Note.** § B.4 and § B.4.1 below are the *reference / audit* view
> of the per-location allocation. **The operational Phase 5 allocation
> is § B.4.2** (fragmentation-aware, adopted as default). Read this
> section for the method underlying the location-scale coupon-collector
> and 704-allele-draw target; use the § B.4.2 per-event `M_frag` column
> for actual field-team numbers.

**Rationale.** Step 28 tells us *how many seeds per mother*. It does
not tell us *how many mothers to sample per event* or *how many events
to sample per location*. Without those two extra layers, per-mother
numbers cannot be aggregated into a location-level estimate with a
stated coverage guarantee.

**Method — uniform-prior version.** We apply the same coupon-collector
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

| View | Question answered | Metric | What 15 seeds/mother buys (tetraploid Rule 2) |
|---|---|---|---|
| **Step 28 (per-mother, [Figure 9](#fig-9))** | How much of ONE mother's local pollen pool do her 15 seeds reveal? | Distinct pollen-donor alleles detected out of K local | 100 % at ≤ 5-plant events → ~19 of 32 at > 20-plant events |
| **Step 29 (per-location, this section)** | How much of the SPECIES-WIDE 32-SRK-allele pool does the location as a whole reveal after pooling across mothers? | Fraction of the 32 Fgs observed at the location | ≥ 90 % once M × 34 ≥ 704 — i.e. from ~21 sampled mothers upward |

**Result — sampling protocol scaled to each mother's mate context.**
Every mother plant received an individually calibrated seed-genotyping
recommendation based on the number of fertile plants at her slickspot
and the number of seeds she produced. Two decision rules run
side-by-side (Rule 1 aspirational, Rule 2 tetraploid operational cap
at 15 seeds). Under the P1 species-wide prior and the tetraploid
`M × 34 draws / mother` accounting, **46 / 52 locations (88 %) reach
the 90 %-of-32-Fgs target with their permit-realistic mother count**
(from Step 29's summary). The 6 that fall short are tiny slickspots
(≤ 10 fertile plants) where no sampling intensity can compensate for
the small census — a mathematical certainty, not a design failure.

Field sampling effort is auditable at the per-mother level. The field
team receives one number per germplasmID
([`step29_field_team_sampling_recipe.tsv`](tables/Phase5/step29_field_team_sampling_recipe.tsv)).

<a id="fig-10"></a>
![Figure 10: Recommended seeds per mother plant for each LEPA location under the P1 empirical prior + tetraploid LEPA (2 paternal alleles per seed). One horizontal bar per location, coloured by achievability tier: **green** (≤ 15 seeds/mother, fits the Step 28 Rule 2 tetraploid floor), **amber** (16–100 seeds/mother, achievable with focused effort), **red** (> 100 seeds/mother, unrealistic — the location's census is too small to characterise). Vertical guides: 15 (Rule 2 tetraploid floor) and 100 (practical ceiling). Bars > 300 are capped for display, with the true value annotated at the right. **19 / 39 locations sit in the green tier; 16 in amber; 4 in red (all single-plant or two-plant slickspots in the EO24 group).** Source: `step29_event_location_sampling.py`.](figures/Phase5/step29_recommended_seeds_per_mother_P1.png)

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

- **Coupon-collector floor** (`M_uniform_full_detection`): under the
  uniform-frequency assumption, the smallest M such that the
  probability of observing **every one** of the location's predicted
  SRK alleles is at least **90 %**. Uses M × 34 allele draws (**4
  maternal + 30 paternal** per mother under tetraploid Rule 2 = 15
  seeds) against the local pool K_local = min(4 · (total fertile
  plants − 1), 32 species alleles).
- **Private-allele floor** (`M_event_coverage`): at least one mother
  per event at the location. This is a *drift-aware* constraint —
  private alleles cannot be caught in an event that is not visited.

The recommendation is `M_recommended = max(coupon-collector floor,
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
coupon-collector maths, it is the number of events that must each be
represented. See [Figure 11](#fig-11) and
[`step28_mothers_for_full_detection_by_location.tsv`](tables/Phase5/step28_mothers_for_full_detection_by_location.tsv).

**Caveat.** The uniform-frequency assumption is optimistic; when
allele frequencies are skewed (which is what drift produces), rare
alleles need substantially more mothers than the 90 % bound suggests.
The private-allele floor is the practical safeguard — it forces at
least one mother per event so no event's private alleles are missed
even when frequency skew is severe.

<a id="fig-11"></a>
![Figure 11: Per-location sampling target to observe every predicted SRK allele at 90 % probability under the uniform coupon-collector model, additionally requiring at least one mother per event (private-allele floor). One row per LEPA location, panelled by Bottleneck Lineage; solid coloured bar = `M_recommended` (the binding floor); open bar with the same colour outline = the uniform coupon-collector bound alone. Row labels give the location code, number of events, adult census and predicted local pool size K. Where the solid bar extends beyond the open bar, the private-allele floor is binding — the location has more events than the coupon-collector maths would ask for, and the extra mothers are needed so no event is skipped. Median `M_recommended = 7`; **31 / 52 locations** have the private-allele floor bind (i.e. `n_events > M_uniform`). Source: `step28_seed_sampling_per_mother.py`. Data: [`step28_mothers_for_full_detection_by_location.tsv`](tables/Phase5/step28_mothers_for_full_detection_by_location.tsv), [`step28_mothers_for_full_detection_by_bin.tsv`](tables/Phase5/step28_mothers_for_full_detection_by_bin.tsv).](figures/Phase5/step28_mothers_for_full_detection.png)

#### B.4.2 Fragmentation-aware sampling — comparison with the current allocation

**Why the B.4.1 allocation under-counts fragmentation.** § B.4.1
computes the coupon-collector M using the location's **aggregate**
`total_N_fertile` as a single pool. But two events sitting 300 m apart
at the same locationID are not one mating unit — no pollen crosses
between them. Treating them as one pool implicitly assumes a mother
sampled at event A pays for detection at event B, which is only true
when A and B are within pollinator range. The event-scale K^(25m)
distribution in [Figure 2](#fig-2) (§ A.5) makes this vivid: at
locations like EO8 (82 events split across 46 disjoint mating units),
the aggregate pool is a mathematical fiction.

**Fragmentation-aware allocation — the method.** For each location we
build the 50 m adjacency graph on its events (same graph § B.2 uses
for connectivity), find the connected components, and apply the
coupon-collector 90 %-full-detection rule **per component** rather
than per whole location:

1. **Per component c** — total N_fertile in c → local pool
   K_c = min(2·(N_c − 1), 32 species alleles) → mothers needed
   M_c = mothers_for_full_detection(K_c, 0.90). Isolated singleton
   components → M_c = 1.
2. **Within a component**, distribute M_c across events proportional
   to N_fertile(e), with a **maternal-genotype floor of ≥ 1 mother
   per event** so every event gets visited (records who is standing
   there, even when its paternal coverage is "free" via a
   neighbouring event's mother in the same component).
3. **Location total**: M_frag_aware = Σ_c mothers allocated to c.

**Result** (tetraploid). All 52 LEPA locations were re-allocated; the
head-to-head comparison against § B.4.1 is in [Figure 12](#fig-12):

- **43 / 52 locations need MORE mothers under fragmentation-aware
  allocation** (Δ from +1 to +114 per location).
- **9 / 52 unchanged (Δ = 0)** — locations that are one single 50 m
  component (or a single event), so per-component and per-location
  coupon-collector maths agree.
- **0 / 52 need fewer** — under the honest allocation, no location
  can reduce effort.
- **Total effort: 748 → 1 712 mothers** (a 2.3× increase in the 90 %
  guarantee's cost across the 2025 field). The tetraploid math
  amplifies this vs the diploid version because K_c per component
  doubles (4·(N−1) instead of 2·(N−1)) while the tetraploid
  per-mother draws (34) only slightly exceed the diploid (31).

**Where the extra mothers come from.** For the top-Δ locations, most
of the increase is the paternal (coupon-collector) side, not the
maternal-genotype floor. EO8 (82 events, 46 components) needs 182
mothers for the paternal 90 % guarantee summed across components; the
maternal floor adds 14 more for a final M_frag_aware = 196. The
current B.4.1 allocation of 82 is understating the honest paternal
cost of that location's fragmentation by ~114 mothers.

**Interpretation.** The B.4.1 allocation is **optimistic** because it
assumes pollen mixes across the whole location. Fragmentation-aware
is the **honest** cost of the same 90 % guarantee under the actual
50 m mating structure. The delta is the number of mothers B.4.1 was
implicitly assuming would come "for free" through pollen sharing that
doesn't actually happen.

**Decision — we adopt fragmentation-aware allocation as the default
for Phase 5.** The 2.0× effort increase is a real permit implication,
but the alternative — publishing a 90 % coverage claim that only
holds under an unrealistic single-pool assumption — is worse for
goal 3 specifically, which is the fragmentation × drift decomposition.
If our sampling protocol builds the fragmentation assumption *in* by
under-allocating at fragmented locations, we cannot then use those
locations to *test* fragmentation as a mate-limitation channel; we
would be arguing from a floor we ourselves lowered. B.4.2 is the
allocation that Phase C's β₂ regression is entitled to inherit.

**What supersedes what.** § B.4.1's per-location `M_recommended`
column and § B.4's proportional-to-N_fertile allocation are now
**reference / audit views**, not the operational recipe. The
operational per-event allocation is
[`step29c_sampling_frag_aware_per_event.tsv`](tables/Phase5/step29c_sampling_frag_aware_per_event.tsv),
column `M_frag` — one row per event, one number per event, which
already satisfies both the per-component coupon-collector 90 %
guarantee and the ≥ 1 mother-per-event maternal-genotype floor.
[`step29c_sampling_comparison_per_location.tsv`](tables/Phase5/step29c_sampling_comparison_per_location.tsv)
keeps both allocations side by side for permit and audit purposes.

**How this hooks into Phase C.** Phase C's mate-limitation regression
(§ C.1) inherits the new allocation automatically: the β₂ predictor is
per-mother K^(25m), which is a property of the landscape and doesn't
depend on which allocation delivered her seeds. What *does* change is
the **statistical power** — the fragmentation-aware allocation adds
mothers exactly at the locations where within-location K^(25m)
variance is highest (§ A.5 "mixed-connectivity" boxes), which is
where β₂ has the most identifiability. So the extra effort is not
distributed randomly across the study — it goes to the locations where
the test needs it most.

**Permit / field-team next step.** The per-event `M_frag` column is
the number the field team needs to hit at each event. The current
[`step29_field_team_sampling_recipe.tsv`](tables/Phase5/step29_field_team_sampling_recipe.tsv)
is based on the older per-mother allocation and should be regenerated
against `M_frag` before the next field season — a one-line change in
Step 29 pulling `M_frag` from the step29c per-event TSV instead of
Step 29's own proportional allocation.

<a id="fig-12"></a>
![Figure 12: Head-to-head sampling recommendation — current § B.4.1 (light bar) vs fragmentation-aware § B.4.2 (solid bar) per LEPA location under tetraploid LEPA, panelled by Bottleneck Lineage. Rows sorted within each BL by Δ = fragmentation-aware − current, ascending. Row labels give locationCode, number of events, number of 50 m components, and total adult census; Δ printed to the right of each row-pair. **43 / 52 locations need more mothers under fragmentation-aware allocation (Δ from +1 to +114); 9 / 52 unchanged (Δ = 0, single-component locations); 0 / 52 need fewer.** Total effort: current 748 mothers → fragmentation-aware 1 712 mothers (2.3×). The biggest Δ locations (EO8 +114, EO32 +81, EO27-1 +79, EO27RT +77) are all locations where a small number of large connected clusters coexist with many spatially isolated singleton events — each singleton contributing 1 mother and each cluster contributing its own coupon-collector M ≈ 6. Source: `step29c_fragmentation_aware_sampling.py`. Data: [`step29c_sampling_comparison_per_location.tsv`](tables/Phase5/step29c_sampling_comparison_per_location.tsv), [`step29c_sampling_frag_aware_per_event.tsv`](tables/Phase5/step29c_sampling_frag_aware_per_event.tsv).](figures/Phase5/step29c_sampling_comparison.png)

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
  mate-limitation and SI-escape tests are *not* pooled but treated as
  repeated measures per location (with year as a fixed effect and
  location as a random effect).

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

### C.0 Preliminary empirical validation at EO scale (Phase 4 adults)

**Question.** Before we invest in seed genotyping, does the sporophytic
Class I / Class II + empirical-zygosity P_compat model (§ A.7)
reproduce what we already observe in the adult population from the
Canu-amplicon pipeline?

**Data.** Per-individual SRK genotypes from `Tables/Phase4/step23_individual_allele_genotypes_with_nulls.tsv`
— 368 tetraploid individuals across 23 Elemental Occurrences (EOs).
Each individual carries 4 SRK copies (functional alleles plus null
copies where SI has decayed). The allele → Fg mapping is the same one
that P1 was built from (`Tables/Phase5/step26i_L1_carrier_inventory.tsv`).
Six EOs meet the inclusion threshold of n ≥ 10 individuals with at
least one functional copy: **EO76 (n = 76), EO70 (n = 74), EO27 (n = 61),
EO25 (n = 51), EO18 (n = 39), EO67 (n = 37)** — 338 individuals in total.

**Design.** For each of the six EOs, two per-mother pollen
compatibilities are computed by Step 30 Part C over the same set of
observed mothers, changing only the father-drawing distribution:

- **Observed P_compat** — real mothers × simulated fathers drawn from
  the **observed local Fg frequency vector at that EO**.
- **Predicted P_compat** — real mothers × simulated fathers drawn from
  the **species-wide P1 prior** (ignoring per-EO drift).

Both use the identical Class I / II dominance rules from § A.7.3 and
the empirical LEPA zygosity distribution from § A.7.3a. The comparison
therefore isolates the effect of **local Fg frequency drift** on
random-mating compatibility, holding the mothers themselves fixed.
Per-EO 95 % confidence intervals come from a mother-resampling
bootstrap with the father pool re-drawn each replicate.

**Result — the model reproduces four of six EOs and cleanly flags a
biological outlier.** Species-mean P_compat under the empirical
zygosity model is 0.68, matching § A.7. Four EOs — EO76, EO27, EO25,
EO67 — sit on the 1:1 diagonal within their 95 % intervals (Fig. C.0,
panel A), confirming that local frequency drift at those EOs is mild
enough not to shift random-mating compatibility away from the
species-wide expectation.

**EO70 is a striking outlier**: observed P_compat = 0.41 (struggling
band) versus predicted 0.60 (sustainable band). This is not a model
failure but a biological signal — at EO70 the local Fg pool is skewed
enough that real mothers overlap far more with their neighbours'
expressed alleles than P1 would predict. **EO18 and EO67** show a
milder drift in the same direction (observed 0.58 versus predicted
0.66). All three are candidates for future intervention: they are
places where the sporophytic model correctly diagnoses reduced
compatibility from drift alone, without needing any seed data.

**Panel B — zygosity distributions per EO.** The distinct-identity
distribution at each EO tracks the species-wide 66 / 32 / 2 %
observation from Step 23 fairly well, with one exception: EO25 shows
substantially more heterozygous individuals (49 % carry two distinct
identities, versus 32 % species-wide), which is consistent with EO25
sitting exactly on the diagonal in panel A — its Fg pool is closer to
being drift-free than the average EO.

**Interpretation and caveats.**

- P1 was built from these same individuals, so the species mean is
  guaranteed to line up on average — this comparison does not test the
  model's absolute calibration, but its **robustness to per-EO drift**.
- The scatter and the pattern of deviations (three EOs lower than
  predicted, none higher; magnitudes ordered by the intuitive severity
  of drift) is consistent with the model behaving correctly under
  realistic local skew.
- EO ≠ location. Locations sit at a finer spatial scale. Once seed
  genotypes are in hand (Phase B), § C.1 rerun at the *location* scale
  will supersede this preliminary check.
- No formal statistical test with n = 6 EOs — this is a **visual +
  quantitative validation** that the model is fit to run against seed
  data.

**Outputs.**

- Table: `tables/Phase5/step30_C_pcompat_validation_at_eo.tsv` — one
  row per EO with observed + predicted P_compat mean and 95 % CI,
  observed Class I mass, and observed distinct-identity distribution.
- Figure: `figures/Phase5/step30_C_pcompat_observed_vs_predicted.png`
  — the two-panel diagnostic described above.
- Script: `step30c_srk_validation_at_eo_level.py`.

### C.1 Test 1 — Mate-limitation regression (goals 1 + 3)

Under strict SI + random mating in a **tetraploid sporophytic** system
(§ A.3, § A.8), a mother's per-mother compatibility is given by the
closed-form § A.8 Case-A / Case-B formulas, driven by the local Fg
frequency vector *f* and the Class I / II assignment in
[`srk_fg_class.tsv`](tables/Phase5/srk_fg_class.tsv):

- **Case A** — mother has ≥ 1 Class I allele; expressed set M ⊆
  Class I; $P_{\text{compat}}(m) = (1 - p(M))^4$.
- **Case B** — mother is all Class II; expressed set M = her 4 alleles;
  $P_{\text{compat}}(m) = 1 - (1 - p_I)^4 + (1 - p_I - p(M))^4$.

The species-mean sporophytic + empirical-zygosity pollen compatibility is ~0.69 (matching, via a different mechanism, the
diploid gametophytic 0.63), so the regression coefficient β₁ lives on
that rescaled axis — the traffic-light bands in § A.8 anchor the
interpretation. Her expected seed set is proportional to
$P_{\text{compat}}(m)$:

$$E[\text{seeds}_m] \;\propto\; \text{ovules}_m \times P_{\text{compat}}(m)$$

The **test** is a mixed-effects regression of observed
`germplasmQuantityEstimate` on predicted $P_{\text{compat}}$:

$$\text{seeds}_m \;=\; \beta_0 + \beta_1 \cdot P_{\text{compat}}(m) + \beta_2 \cdot K^{(25\text{m})}(m) + u_{\text{location}(m)} + u_{\text{year}(m)} + \epsilon_m$$

where $K^{(25\text{m})}$ is the mother's mating-neighbourhood
pollen-donor pool at the 50 m radius (§ B.3). Two coefficients, two
distinct causal channels:

- $\beta_1 > 0$ with 95 % CI excluding 0 → **evidence of mate
  limitation driven by allele-frequency drift**. Under the sporophytic
  model, this fires when the drift-shifted local Fg composition shifts
  the mother out of the "sustainable" pollen compatibility band. Because between-
  class compatibility (Class I × Class II) buffers most locations
  against drift (§ A.8), the sporophytic β₁ signal is expected to be
  subtler than the diploid gametophytic version would have implied —
  the mate-limited locations are the ones where drift has removed a
  whole class from the local pool.
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
- **When seed data land**: the mate-limitation regression (§ C.1) and
  the SI-escape rate test (§ C.2) both fire from the same seed-genotype
  input. Together they classify each location into a 2 × 2 matrix (mate-
  limited yes/no × SI-escaped yes/no) that is directly interpretable as a
  conservation prioritisation.
- **Fragmentation × drift decomposition**: the two coefficients
  $\beta_1$ (pollen-compatibility effect) and $\beta_2$
  (mating-neighbourhood-size effect) from the mate-limitation regression
  separate the two mechanisms even though they both reduce reproductive
  success — a distinction unreachable with any single-year,
  single-location analysis.
- **Two-generation efficiency**: every mother's seed lot pays double —
  it certifies her genotype (maternal inventory) and samples her pollen
  environment (paternal inventory) in one experiment.
- **Phenotype cross-validation is a natural next step**: the two-by-two
  matrix from § C.1 + § C.2 predicts what ISI / fruit set from
  Genetic-Rescue-DB should look like at each location, which can be
  overlaid without changing any of the code paths in this framework.

### C.6 How Phase B data are generated and consumed

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
   tetraploid Rule 2**, capped at her real seed budget; 2 paternal
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

<a id="fig-13"></a>
![Figure 13 (DEMO): Mate-limitation regression preview. One dot per LEPA location, colour = "sustainable" band. X = predicted random-mating pollen compatibility (mean across sampled mothers at that location); Y = mean observed seeds per mother. Traffic-light background bands (failed / struggling / sustainable). Dashed line = weighted OLS fit, slope + p-value printed in the legend. In this DEMO the simulator baked in a direct causal link (seed set ∝ compatibility), so the slope is highly significant. **With real data the same figure will test whether observed seed set actually declines with predicted compatibility — a positive slope with 95 % CI excluding 0 confirms mate limitation at population level.** Source: `step30_srk_diversity_prediction_vs_observed.py --demo`.](figures/Phase5/step30_B_DEMO_mate_limitation.png)

<a id="fig-14"></a>
![Figure 14 (DEMO): Self-incompatibility escape preview. One horizontal bar per LEPA location, sorted by observed rate. X = observed rate of pollen alleles matching the mother's own SRK alleles (= self-incompatibility escape rate). Red bars = locations that reject the strict-SI null at 5 % false-discovery rate; grey bars = consistent with strict SI. In this DEMO the simulator baked in 8 % SI escape rate, so most locations show detectable escape. **With real data any red bar names a candidate partial-SI population — a location where the SI machinery has broken down enough that self-pollen produces seeds.** Source: `step30_srk_diversity_prediction_vs_observed.py --demo`.](figures/Phase5/step30_B_DEMO_si_escape.png)

**Two-year extension.** The `--year` flag on Steps 28, 29 and 30
makes the pipeline year-scoped. When 2026 field data arrives, run
each step twice (once per year); the pooling rules under § B.5 tell
Phase B how to combine years for standing-diversity inference vs how
to keep them separate as repeated measures for mate-limitation and
SI-escape tests.

**Downstream: cross-validation with phenotype.** The mate-limitation
and SI-escape results define a 2 × 2 classification per location
(mate-limited yes/no × SI-escaped yes/no). Phase III of the wider
project overlays the Genetic-Rescue-DB ISI / fruit-set phenotype
against these predictions to provide an independent line of evidence.

### C.7 BL4 pilot study — one small + one large location

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
| Effective N under 50 m connectivity | 6 | 116 |
| Mothers with seed records in DB | **4** | **33** |
| Seeds / mother (tetraploid Rule 2 cap) | 15 (budget-limited if lower) | 15 |
| Total allele draws at this location | 4 × 4 + 4 × 30 = 136 draws | 33 × 4 + 33 × 30 = 1 122 draws |
| **Predicted local pool size** (SRK alleles at this location) | ~8 alleles | ~27 alleles |
| **Predicted local coverage** (of the alleles at this location) | **~100 %** ✓ | **~100 %** ✓ |
| Predicted species-wide coverage (of 32 P1 alleles — biased against drifted-out alleles) | ~28 % | ~68 % |
| Predicted random-mating compatibility (sporophytic Class I / II + empirical zygosity, § A.7-A.7) | mean **0.67**, 95 % CI [0.42, 0.87] — *"sustainable"* band at the mean, CI stays sustainable throughout | mean **0.69**, tight CI [0.59, 0.77] — solidly *"sustainable"* |

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
**budget-limited / founder-effect regime** — small population, few
mothers; the mean sporophytic + empirical-zygosity pollen
compatibility sits near the species mean (0.67), and the 95 % CI
[0.42, 0.87] stays inside the "sustainable" band. Any Phase C
observation dropping this location out of sustainable would be
decisive evidence of drift-driven mate limitation. EO27-1 tests the
**aggregation regime** (many mothers, permit-realistic 33 mothers ×
34 tetraploid draws per mother = 1 122 allele draws, easily crossing
the A_target = 704 species threshold; the sporophytic + empirical-
zygosity CI [0.59, 0.77] anchors the "sustainable" prediction).
Both locations
sit in the same BL, so their pilot outputs are directly comparable —
a *within-BL* contrast that avoids between-BL confounds. The
two-location pilot uses ≤ 555 seed genotypes total (< 5 % of the full
2025 recipe of ~11 500 under § B.4.1, or ≤ 1 200 seeds and < 5 % of
the ~25 000 total under the § B.4.2 fragmentation-aware allocation).

**Field-team recipe for the pilot.** Under the § B.4.2 default, the
authoritative per-event allocation is `M_frag` in
[`step29c_sampling_frag_aware_per_event.tsv`](tables/Phase5/step29c_sampling_frag_aware_per_event.tsv);
filter to the pilot `locationID`s to isolate the two locations. The
older
[`step29_field_team_sampling_recipe.tsv`](tables/Phase5/step29_field_team_sampling_recipe.tsv)
is kept for reference / audit but no longer drives field effort — it
needs to be regenerated against `M_frag` before the next field
season.

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
- [`step28_aggregation_curves_by_Nfertile.tsv`](tables/Phase5/step28_aggregation_curves_by_Nfertile.tsv) — aggregation curves showing expected distinct alleles at a location for M = 1 … 30 mothers × 15 seeds each (tetraploid Rule 2), one row per (event-size bin, M).
- [`step28_mothers_for_full_detection_by_bin.tsv`](tables/Phase5/step28_mothers_for_full_detection_by_bin.tsv) — uniform-model M for a 90 % chance of observing every one of K_local alleles, one row per event-size bin.
- [`step28_mothers_for_full_detection_by_location.tsv`](tables/Phase5/step28_mothers_for_full_detection_by_location.tsv) — per-location `M_recommended` combining the coupon-collector bound with the private-allele floor (≥ 1 mother per event).
- [`step29c_sampling_frag_aware_per_event.tsv`](tables/Phase5/step29c_sampling_frag_aware_per_event.tsv) — per-event fragmentation-aware allocation under the 50 m connected-component decomposition; columns include `component_id`, `component_K`, `component_M_target`, `M_paternal_only`, `M_frag`.
- [`step29c_sampling_comparison_per_location.tsv`](tables/Phase5/step29c_sampling_comparison_per_location.tsv) — per-location head-to-head of current (§ B.4.1) vs fragmentation-aware (§ B.4.2) M_recommended, with `n_components_50m`, `sum_M_paternal_only`, and `delta`.
- [`step29_sampling_per_event.tsv`](tables/Phase5/step29_sampling_per_event.tsv) — per-event allocation.
- [`step29_sampling_per_location.tsv`](tables/Phase5/step29_sampling_per_location.tsv) — per-location design table with connectivity-informed columns.
- [`step29_location_coverage_curves.tsv`](tables/Phase5/step29_location_coverage_curves.tsv) — analytical curves by location size.
- [`step29_field_team_sampling_recipe.tsv`](tables/Phase5/step29_field_team_sampling_recipe.tsv) — **FIELD TEAM per-germplasmID recipe** (one row per mother with the actionable seed count).
- [`step29_location_connectivity.tsv`](tables/Phase5/step29_location_connectivity.tsv) — within-location pollen connectivity at 10 / **25 (primary)** / 50 m.
- [`step30_A_prediction_prior_frequencies.tsv`](tables/Phase5/step30_A_prediction_prior_frequencies.tsv) — the P1 species-wide Fg prior.
- [`step30_A_prediction_location_diversity.tsv`](tables/Phase5/step30_A_prediction_location_diversity.tsv) — predicted SRK allele diversity per location.
- [`step30_A_prediction_location_pcompat.tsv`](tables/Phase5/step30_A_prediction_location_pcompat.tsv) — predicted random-mating pollen compatibility per location (finite-population model at 50 m).
- [`step30_A_prediction_per_mother_fecundation.tsv`](tables/Phase5/step30_A_prediction_per_mother_fecundation.tsv) — species-wide compatibility reference distribution under the sporophytic tetraploid model.
- [`srk_fg_class.tsv`](tables/Phase5/srk_fg_class.tsv) — Fg → dominance class (I / II) mapping used by the sporophytic pollen compatibility model in § A.8. Provisional data-driven default; editable by hand as biology is refined.
- [`srk_zygosity_empirical.tsv`](tables/Phase5/srk_zygosity_empirical.tsv) — empirical LEPA distribution of distinct functional SRK identities per plant, from Canu-amplicon Step 23 (n = 367). Used to draw mother and father genotypes under the sporophytic finite-population model (§ A.8.3a, § A.8).
- [`step30_A_traffic_light_bands.tsv`](tables/Phase5/step30_A_traffic_light_bands.tsv) — recalibrated failed / struggling / sustainable band boundaries against the sporophytic species-mean pollen compatibility.
- [`step30_A_fragmentation_per_event.tsv`](tables/Phase5/step30_A_fragmentation_per_event.tsv) — per-event fragmentation index F_event = 1 − K_spatial_50m / 32 (purely spatial, no allele frequencies).
- [`step30_A_fragmentation_per_location.tsv`](tables/Phase5/step30_A_fragmentation_per_location.tsv) — per-location fragmentation index F_location = 1 − within-location connectivity at 50 m, with median event-scale F for the same location.

**Figures** (PNG + PDF):

- [`step28_coverage_curves.pdf`](figures/Phase5/step28_coverage_curves.pdf) / [`.png`](figures/Phase5/step28_coverage_curves.png) — per-mother coverage vs seeds and aggregation across mothers, one curve per event-size bin, with Rule 2 cap.
- [`step28_per_mother_budget.pdf`](figures/Phase5/step28_per_mother_budget.pdf) / [`.png`](figures/Phase5/step28_per_mother_budget.png) — per-mother seed budget vs recommended n.
- [`step28_mothers_for_full_detection.pdf`](figures/Phase5/step28_mothers_for_full_detection.pdf) / [`.png`](figures/Phase5/step28_mothers_for_full_detection.png) — per-location M_recommended (90 % chance to see every predicted SRK allele + ≥ 1 mother per event), panelled by BL.
- [`step29c_sampling_comparison.pdf`](figures/Phase5/step29c_sampling_comparison.pdf) / [`.png`](figures/Phase5/step29c_sampling_comparison.png) — head-to-head M_current vs M_frag_aware per location, BL-panelled, sorted by Δ.
- [`step29_location_coverage_curves.pdf`](figures/Phase5/step29_location_coverage_curves.pdf) / [`.png`](figures/Phase5/step29_location_coverage_curves.png) — uniform-K location curves.
- [`step29_location_coverage_curves_P1.pdf`](figures/Phase5/step29_location_coverage_curves_P1.pdf) / [`.png`](figures/Phase5/step29_location_coverage_curves_P1.png) — P1-prior location curves with per-location bars.
- [`step29_recommended_seeds_per_mother_P1.pdf`](figures/Phase5/step29_recommended_seeds_per_mother_P1.pdf) / [`.png`](figures/Phase5/step29_recommended_seeds_per_mother_P1.png) — per-location seed-genotyping recipe (green/amber/red tiers).
- [`step29_location_connectivity.pdf`](figures/Phase5/step29_location_connectivity.pdf) / [`.png`](figures/Phase5/step29_location_connectivity.png) — **primary connectivity map at 50 m**, BL-panelled.
- [`step29_location_connectivity_10m.pdf`](figures/Phase5/step29_location_connectivity_10m.pdf) / [`.png`](figures/Phase5/step29_location_connectivity_10m.png) — sensitivity: 10 m conservative.
- [`step29_location_connectivity_50m.pdf`](figures/Phase5/step29_location_connectivity_50m.pdf) / [`.png`](figures/Phase5/step29_location_connectivity_50m.png) — sensitivity: 50 m optimistic.
- [`step29_location_connectivity_radius_sensitivity.pdf`](figures/Phase5/step29_location_connectivity_radius_sensitivity.pdf) / [`.png`](figures/Phase5/step29_location_connectivity_radius_sensitivity.png) — aggregate across radii (justification for 50 m primary).
- [`step30_A_prediction_diversity.pdf`](figures/Phase5/step30_A_prediction_diversity.pdf) / [`.png`](figures/Phase5/step30_A_prediction_diversity.png) — BL-panelled SRK allele richness per location.
- [`step30_A_si_model_schematic.pdf`](figures/Phase5/step30_A_si_model_schematic.pdf) / [`.png`](figures/Phase5/step30_A_si_model_schematic.png) — three-panel pedagogical schematic of the sporophytic Class I / II model (§ A.7): dominance within a plant, worked example, compatibility rule by cross type.
- [`step30_A_prediction_fecundation.pdf`](figures/Phase5/step30_A_prediction_fecundation.pdf) / [`.png`](figures/Phase5/step30_A_prediction_fecundation.png) — BL-panelled compatibility per location, traffic-light bands.
- [`step30_A_diversity_vs_pcompat.pdf`](figures/Phase5/step30_A_diversity_vs_pcompat.pdf) / [`.png`](figures/Phase5/step30_A_diversity_vs_pcompat.png) — cross-plot: SRK diversity × pollen compatibility.
- [`step30_A_fragmentation_index.pdf`](figures/Phase5/step30_A_fragmentation_index.pdf) / [`.png`](figures/Phase5/step30_A_fragmentation_index.png) — event-scale × location-scale fragmentation scatter, BL-coloured (pure spatial, no drift).

**Phase A validation (adult SRK genotypes from Phase 4 — Part C § C.0):**

- [`step30_C_pcompat_validation_at_eo.tsv`](tables/Phase5/step30_C_pcompat_validation_at_eo.tsv) — per-EO observed vs predicted pollen compatibility (n ≥ 10 individuals), observed Class I mass, and observed distinct-identity distribution.
- [`step30_C_pcompat_observed_vs_predicted.pdf`](figures/Phase5/step30_C_pcompat_observed_vs_predicted.pdf) / [`.png`](figures/Phase5/step30_C_pcompat_observed_vs_predicted.png) — two-panel EO-scale diagnostic (observed vs predicted P_compat + observed zygosity distribution per EO).

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
