# Critical Review of the Polyploid Crossing Model

**Reviewer:** Claude Opus 4.6 (1M context), in collaboration with Jim Beck
**Date:** 2026-05-23
**Scope:** Source code (`src/polyploid_utils.py`, 570 LOC), the GFS refactor plan (`Genotype_Quality_Crossing_Plan.md`), and consistency with upstream framework (`svenbuerki/SRK_bioinformatics`).
**Tone:** Balanced, tilted toward steel-man — every concern stated as a reviewer would state it, with severity tags and counterpoints where warranted.

## Severity legend

- 🚨 **Blocker** — must be addressed before the planned refactor proceeds, or the refactor will encode a known defect.
- ⚠️ **Significant** — should be addressed during the refactor; defensible to ship without if explicitly acknowledged.
- 📌 **Design** — open question or tradeoff that needs an explicit decision, not necessarily a fix.
- ✅ **Validated** — checked and correct.

## Headline findings (read first)

1. 🚨 **The model's self-incompatibility rule contradicts the biological mechanism the Crossing Plan claims to model (the C3 cascade).** Brassicaceae sporophytic SI with competitive interaction allows heterozygous (AB) pollen to escape rejection on a mother carrying both A and B — this is the Lewis (1947) competitive-interaction model that the LEPA report (Q6 Stage 5) identifies as the primary driver of AAAA accumulation. The current `is_compatible` does the opposite: it rejects AB pollen on any mother containing A OR B. The Crossing Plan's GFS refactor adds GFS as a tracked metric but doesn't fix this. Without a change to `is_compatible`, the simulation cannot reproduce the C3 mechanism the conservation strategy is designed against.

2. 🚨 **The model treats inheritance as tetrasomic (uniform C(4,2)=6 gamete sampling), but LEPA is allopolyploid — the majority of allopolyploids are disomic.** Sven flags this caveat himself in the LEPA report Q4 (line 308): "If LEPA's SRK locus shows disomic inheritance ... current P_compat values are then a conservative lower bound." Under disomic inheritance the gamete distribution is structured (two homeologous subgenomes segregate independently), not uniform. This affects every prediction the model makes about offspring genotypes. Verifying the inheritance mode for the LEPA S-locus is a prerequisite for trusting the model's quantitative output.

3. 🚨 **Sven's central conservation diagnostic (TP1 four-quadrant) is not implemented in the polyploid model.** The LEPA report's most consequential finding (Q4) is the TP1 mating-pool functionality scatter, which maps each population into MONITOR / AUGMENT / BIOBANK conservation actions. The polyploid model optimizes for χ² distance to uniform — a related but distinct objective. To produce management-actionable output consistent with the LEPA framework, the model must report J × P_compat per simulation generation and per BL/EO.

4. 🚨 **The model has no concept of BL structure.** All individuals are pooled. But the central LEPA finding is: "contemporary recovery requires inter-BL allele transfers" (Q4). The model cannot recommend, evaluate, or simulate inter-BL crosses if it doesn't carry BL annotations through to the optimizer and reporting.

5. ⚠️ **`simulate_generation` includes a safety net that auto-inserts allele carriers if any allele was lost during a generation.** This is good for ensuring zero-loss reporting but it masks the actual dynamics; "preservation strategy: 0 alleles lost" results are partly the safety net at work, not the optimizer.

The rest of this document develops these and supporting findings in detail.

---

## Section A — Biological assumptions

### A.1 🚨 SI rule contradicts the C3 mechanism

**The claim under review.** `polyploid_utils.is_compatible(maternal_genotype, pollen_gamete)` returns False if any of the 2 alleles in `pollen_gamete` is in `maternal_genotype`. The `cross` function uses this to filter paternal gametes before forming offspring.

**Biological reality (from LEPA report Q6 Stage 5, lines 518–520):**

> "an AABB individual produces diploid pollen by sampling two of its four allele copies, generating three pollen types — AA (1/6), **AB (4/6 — most common)**, and BB (1/6). AB pollen carries both A-SCR and B-SCR proteins simultaneously; on an AABB pistil expressing both SRK-A and SRK-B receptors, the two competing recognition signals interfere with each other, producing a response too weak to trigger rejection. Self-fertilisation succeeds, producing predominantly AAAB offspring."

This is the Lewis (1947) competitive-interaction model, documented empirically in *Arabidopsis lyrata* (Mable et al. 2005; Mable & Adam 2007) and canonical in polyploid Brassicaceae SI (Mable 2004).

**Where the model and biology disagree.** Concrete example: father AABB × mother AABB.

| Gamete | Probability | Biology (C3 Stage 5) | Model's `is_compatible` |
|---|---|---|---|
| AA | 1/6 | Rejected (A is in mother) | Rejected (A is in mother) |
| AB | 4/6 | **Succeeds (competitive interaction)** | **Rejected (both A and B in mother)** |
| BB | 1/6 | Rejected (B is in mother) | Rejected (B is in mother) |
| **Effective selfing rate** | | **4/6 = 67%** | **0%** |

The model's rule makes selfing on AABB completely impossible. The C3 mechanism — the entire motivation for the GFS refactor — cannot occur in the simulation. The simulated population can only become AAAA via deterministic drift (sampling), not via the biological mechanism the conservation strategy is designed to counter.

**Counterpoint considered.** Sven's upstream `P_compat` formula (LEPA report Q4 Method, line 292) uses *the same rule as the model*: "pollen rejected if any of its 2 alleles matches any of the stigma's 4". He labels this "sporophytic SI with co-dominance" but the formula is gamete-vs-stigma. The two repos use the same SI rule, so the model isn't inconsistent with the upstream `P_compat` calculation.

But Sven then introduces **leakage parameter L** to model the C3 Stage 5 dynamics empirically: `P_compat(L) = P_compat_strict + L·(1 − P_compat_strict)`, with `L̂ ≈ 0.18` estimated from observed AAAA proportions. **Leakage is the upstream framework's way of capturing C3 without changing the strict SI rule.** The polyploid model has no equivalent — it operates only at L=0 (strict SI).

**Severity rationale.** The Crossing Plan's GFS refactor is motivated entirely by Section 5 of the LEPA report (C3 cascade). If the simulation cannot reproduce that mechanism, the refactor adds bookkeeping (GFS as a metric, AAAA preservation filter on elitism) without fixing the core dynamics. The strategy "prevent AAAA accumulation by refusing to use AAAA as pollen donor" is sound conservation advice but is testable only against a simulation that can produce AAAA via the biological mechanism — which the current model cannot.

**Recommended action.** Add a leakage parameter `L ∈ [0, 1]` to `is_compatible` (or, more cleanly, to `cross`) that admits otherwise-rejected pollen with probability L. Default L=0 preserves current behavior; L=0.18 reproduces the empirically observed AAAA rates. Or, more biologically faithful: implement competitive-interaction directly (AB pollen escapes when both A and B are in mother), as the LEPA report describes. Discuss with Sven which framing he prefers — the leakage parameter is simpler; the mechanistic implementation is more biologically faithful and lets the rate emerge from first principles.

### A.2 🚨 Inheritance mode is unverified

**The claim under review.** `form_gametes` returns all C(4,2)=6 ordered allele pairs uniformly. This assumes tetrasomic inheritance (random pairing across all four chromosome copies).

**Biological reality.** LEPA is allopolyploid. The general principle in the literature: most allopolyploids show disomic inheritance because the homeologous subgenomes pair preferentially with themselves, not randomly across all four copies. From the literature search (see Sources below):

> "While the majority of polyploids exhibit disomic inheritance, a few important crops are 'autopolyploids', exhibiting polysomic segregation such as alfalfa, some grass species, coffee and potato."

> "Disomic differs from polysomic inheritance in the way alleles at a locus are selected to form gametes — in an allotetraploid, instead of random selection from four chromosomes, the set is divided into two 'homeologous loci' with two chromosomes in each, and one chromosome is chosen from each locus to form gametes."

Under disomic inheritance with subgenomes (S₁, S₁') and (S₂, S₂'), the gamete distribution is the *outer product* of the two subgenome diploid pairings — only 4 of the 6 ordered pairs are accessible (S₁S₂, S₁S₂', S₁'S₂, S₁'S₂'), not all 6. The unreachable pairings are the within-subgenome ones (S₁S₁', S₂S₂').

**The empirical question.** Has LEPA's S-locus been shown to inherit tetrasomically or disomically? I don't have access to Buerki et al. 2026 (cited in the LEPA report) to confirm. Sven himself acknowledges the uncertainty (LEPA report Q4 line 308). This is a question for Sven, not something inferable from the data files alone.

**Why it matters.** If LEPA is disomic at the S-locus:

- The C(4,2)=6 gamete model overcounts heterozygous combinations and undercounts homozygous combinations.
- An AABB plant with disomic inheritance always produces AB gametes (1 from each subgenome) — never AA or BB. The C3 mechanism still applies but the rate distribution is different.
- All allele-frequency predictions, GFS values, and offspring distributions are quantitatively biased.

**Recommended action.** Add a comment + parameter to `form_gametes` flagging the tetrasomic assumption. Long-term: parameterize gamete formation by inheritance mode and provide both implementations. Raise as Sven question #6 (in addition to the existing 5 in the Crossing Plan).

### A.3 ⚠️ Dominance hierarchies among S-alleles not modeled

**Literature.** *Arabidopsis lyrata* (a close Brassicaceae tetraploid relative) shows "complicated dominance interactions among alleles in the diploid parent [that] determine self-recognition phenotypes of both pollen and stigma" (Mable et al., Heredity 2004). In sporophytic SI, only dominant SCR alleles are expressed on the pollen surface; recessive alleles are masked. The polyploid model treats all 4 alleles equally (full co-dominance).

**Upstream signal.** SRK_bioinformatics Step 12b (in progress) is assigning Class I (dominant) and Class II (recessive) by phylogenetic position. In the 2026-05-11 dataset, Allele_055 is the only observed Class II allele.

**Why it matters.** Dominance hierarchies change which alleles participate in rejection. In an AABB father where A is Class I (dominant) and B is Class II (recessive), only A-SCR is expressed → pollen carries A-SCR only → rejection depends on A in mother, not B. This changes the rejection logic substantially and may explain some of the observed paradoxes (functional SI machinery + 60% AAAA).

**Recommended action.** Until Step 12b is finalized, document the co-dominance assumption explicitly. Treat as a known limitation; revisit when Sven's class assignments are stable.

### A.4 📌 NFDS equilibrium as uniform allele frequencies — defensible but imperfect

The target frequency `1/n_alleles` for chi-squared optimization assumes uniform equilibrium. Under NFDS in infinite populations, equilibrium is uniform; in finite populations with Ne < ~50, allele frequencies fluctuate substantially around uniform (Schierup 1998; Vekemans 1998). For the LEPA dataset (Ne ≈ 15–28 per EO, per `MODEL_OVERVIEW.md`), uniform is a reasonable target but the *variance around uniform at equilibrium* is non-negligible.

**Not a bug, but a framing point.** Optimizing aggressively toward exact uniform is more stringent than necessary. A more sophisticated objective would target the expected stationary distribution under NFDS at the given Ne — but this is research-grade, not a critique of an applied tool.

### A.5 ✅ GFS formula correct

Verified against the canonical reference table:

| Tier | (n_k) | Σ n_k(n_k-1)/12 | GFS |
|---|---|---|---|
| ABCD | (1,1,1,1) | 0/12 | 1.000 ✓ |
| AABC | (2,1,1) | 2/12 | 0.833 ✓ |
| AABB | (2,2) | 4/12 | 0.667 ✓ |
| AAAB | (3,1) | 6/12 | 0.500 ✓ |
| AAAA | (4) | 12/12 | 0.000 ✓ |

The formula is the proportion of `C(4,2)=6` ordered gamete pairs that are heterozygous (carry two distinct alleles). Matches Sven's implementation in `SRK_individual_GFS.R` exactly.

### A.6 ⚠️ AB-bypass and Q5 paradox both invisible to the current model

The LEPA report (Q5 line 484): "An apparent paradox: intact SI machinery, 60 % reproductive dead-ends." The resolution: leaky SI via competitive interaction (Q6 Stage 5) — exactly the mechanism the model does not implement. The conservation strategy literally exists to counter a mechanism the model cannot simulate.

---

## Section B — Mathematics and algorithms

### B.1 ✅ Chi-squared objective is a valid distance measure

`chi_sq = Σ_a (expected_a - 1/n)² / (1/n)` is the standard Pearson chi-squared statistic between expected and uniform target. Symmetric, smooth, and proportional to variance up to scaling. KL divergence (which the code computes but doesn't optimize) is a defensible alternative; the choice matters more for interpretation than for finding the same minimum.

### B.2 ✅ Gradient correctness (verified by chain rule)

The optimizer's analytical gradient through the normalization `w_norm = w / w_sum` is correctly implemented:

$$\frac{\partial f}{\partial w_j} = \frac{1}{w_{sum}} \left( g_j - \sum_k w_{norm,k} g_k \right) \quad \text{where } g_k = (M^T r)_k$$

`grad = (g - np.dot(w_norm, g)) / w_sum` matches this exactly. No issues here.

### B.3 ⚠️ `np.abs(weights)` is redundant given the bounds

`compute_optimal_weights` uses both `bounds=[(0, 1)] * n_crosses` (which L-BFGS-B respects) and `w = np.abs(weights)` inside the objective. The abs() is redundant and introduces a discontinuity at w=0 that could complicate gradient computation if any weight rests exactly on the lower bound. Harmless in practice but should be removed for cleanliness.

### B.4 ⚠️ Preservation penalty: actively pushes UP, not just preserves

The penalty `preservation_weight * max(0, target_freq - freq_i) ** 2` is one-sided: it activates only when `freq_i < target_freq` (i.e., 1/n). For a rare allele with starting frequency f₀ < 1/n, the penalty pushes its expected frequency *toward 1/n*, not just preserving f₀. This is fine for NFDS-driven convergence but means there's no concept of "preserve the rare allele at its starting frequency without amplifying it". If a user expected "preservation = no change", they'll be surprised.

**Recommendation.** Either rename `preservation_weight` to `convergence_floor_weight` to reflect what it actually does, or add a distinct "soft floor at current frequency" option for true preservation semantics.

### B.5 ⚠️ `effective_population_size` is heuristic, not derived

```python
Ne = n_breeders * mean_compatibility
```

This is presented as a formula but it's not derived from any standard Ne formula in the population genetics literature. The standard reference for Ne under SI is Vekemans & Slatkin (1994) and Wright (1939), which give specific corrections for SI based on allele numbers and dominance. The current heuristic could be benchmarked against these (or against direct simulation with allele identity tracking).

**Counterpoint.** For a quick comparison metric the heuristic may be adequate; Ne is used in the model only for reporting, not for the optimizer.

### B.6 ⚠️ Rare-allele threshold 0.05 will misbehave post-synonymy

With 27 biological alleles, the uniform target is `1/27 ≈ 0.037`. A threshold of `0.05` for "rare" exceeds the uniform target, so at the NFDS equilibrium *every* allele would be flagged as rare. This breaks the elitism + mandatory-cross logic — everything becomes mandatory.

**Recommendation.** Recalibrate threshold for the new allele count. Two options:
- Set rare_threshold relative to target (e.g., `0.5 * target_freq`)
- Switch from frequency threshold to carrier-count threshold (the function `get_mandatory_rare_crosses` already uses `0 < c <= 2` carriers, which doesn't suffer this scaling issue)

### B.7 ⚠️ `enumerate_compatible_crosses` allele effect uses `prob / 4.0` per allele copy

`expected_freqs[allele] += prob / 4.0` for each allele in each offspring. This computes the per-cross *contribution* to allele frequency for ONE offspring drawn from that cross. When the optimizer aggregates `w_norm @ allele_effect_matrix`, it's computing expected frequency under uniform-offspring-per-cross sampling. Correct given the modeling choice, but worth noting that this treats each cross as producing exactly one offspring — it doesn't account for differential reproductive success (e.g., crosses with high compatibility producing more offspring per attempted pair).

---

## Section C — The proposed refactor (Genotype_Quality_Crossing_Plan.md)

### C.1 🚨 GFS refactor as currently described doesn't fix C3 — it instruments the symptom

The Plan §1 correctly identifies that `select_elites` over-rewards AAAA carriers of rare alleles. Adding GFS as an instrumentation metric and gating elitism on GFS prevents *amplifying* AAAA. But the *production* of AAAA via the C3 mechanism happens in the simulation itself, which is governed by `is_compatible` + `cross`. Without a change to those, the simulation cannot produce AAAA via the AB-bypass route described in Q6.

**Two possible refactor paths:**

1. **Instrument-only** (current Plan). Adds GFS tracking, gates elitism, gates mandatory rare-allele crosses. Defensible if the goal is to ensure the *managed strategy* doesn't make things worse. Cannot evaluate whether the strategy *fixes* C3.

2. **Mechanism-faithful** (extension). Modify `is_compatible` / `cross` to admit AB pollen on mothers carrying both A and B. Then GFS becomes naturally important because heterozygous gametes are the "useful" ones. The strategy becomes evaluable against the C3 mechanism.

The Plan as written is path 1. Worth deciding explicitly which one we're doing.

### C.2 ⚠️ "5 AABC individuals" is now stale (it's 15 at the protein level, fewer after synonymy collapse)

The Plan §2.2 says "5 AABC individuals exist — all in EO67 and EO70". The 2026-05-11 dataset has **15** AABC individuals at the protein level. After synonymy collapse, many will reclassify (Plan-flagged: `Allele_050(1)+Allele_049(2)` is protein-level AABC → biological-level AAAB because both are in Synonymy group 1). The final biological-level AABC count is unknown until the synonymy-aware loader runs.

**Action.** Defer Plan Phase 3 (rescue crosses, EO67/EO70 specific) until the new loader emits the post-collapse genotype-tier distribution. Then revisit feasibility.

### C.3 📌 Plan §4.1: AAAB as maternal-only — biologically conservative but worth a sensitivity check

The Plan's preferred option is "AAAB never used as pollen donor" (Recommendation b). This prevents the AAAB → AAAA slide that Q6 Stage 5 identifies as runaway. But AAAB carriers of rare alleles in BL2/BL4 (the diversity reservoirs) could be valuable inter-EO pollen sources for AAAA-dominated recipient EOs that lack the AAAB father's rare allele. Worth a sensitivity analysis comparing strategies "AAAB excluded from pollen donor" vs "AAAB allowed if partner-EO lacks the allele entirely" before locking in.

### C.4 📌 Plan §4.2: inter-EO operational feasibility is Sven's call

Plan §5 Q2 already raises this. Holding pattern — needs Sven's input.

### C.5 ⚠️ Plan §4.3 underestimates how much of the framework is in upstream

The Plan §4.3 frames Option B ("fold into Sven's notebook 08") as conditional on Sven having scaffolding. Survey of upstream confirms `notebooks/` doesn't exist there in Python form, BUT the **whole analytical framework already does**:

- GFS, TP1, TP2 metrics: `Scripts/SRK_individual_GFS.R`, `Scripts/SRK_TP1_compatibility_metrics.py`
- BL ordering and palette: `Scripts/srk_bl_constants.py` (mirrored to `.R`)
- Synonymy groups + cross categories: `Scripts/srk_allele_hypotheses.py` outputs
- Cross plan: `Scripts/srk_cross_plan.py` already generates a phased H0–H3 protocol

Option A (extend `polyploid_utils.py` + new notebook 09) is right, but the new notebook should **reuse upstream's Python module imports where possible**. Specifically: importing `srk_bl_constants` directly into the polyploid model avoids redefining BL_ORDER / BL_COLORS and guarantees consistency.

### C.6 ⚠️ Plan §3 Phase 3 inter-EO rescue depends on operational feasibility AND on the SI rule

Even if Sven confirms inter-EO transfers are operationally feasible, the rescue cross AAAA × AABC depends on:
1. The pollen-bearing parent having AABC genotype (verifiable from data)
2. The maternal AAAA parent accepting AABC pollen — which requires the AABC pollen's gametes to NOT share any allele with the AAAA mother

With AAAA mother (alleles X X X X) and AABC father (alleles Y Y Z W with W ≠ X, Y ≠ X, Z ≠ X), the AABC father produces gametes: YY (1/6), YZ (2/6), YW (2/6), ZW (1/6). All have no X → all accepted. So 100% of AABC pollen succeeds on AAAA mother (under the current SI rule). The mechanism works.

But under sporophytic SI with co-dominance (full canonical Brassicaceae model), the AABC father's pollen carries Y+Z+W SCR coatings; if mother is XXXX, all are non-X so all accepted. Same result.

✅ This particular cross is robust to SI-rule choice. Good — Phase 3 still works.

---

## Section D — Implementation traps

### D.1 ✅ Wide matrix is protein counts, not copy counts (caught and documented)

Already saved as a feedback memory and in `CLAUDE.md`. The pipeline ahead must compute copy counts via `Genotype` + `Allele_composition` parsing, never directly from the wide matrix.

### D.2 ⚠️ Synonymy collapse must happen BEFORE every downstream computation

Order of operations in the new loader:
1. Load wide protein-count matrix and zygosity.
2. Build biological-allele mapping from `SRK_synonymy_groups.csv`.
3. Collapse protein-level allele dosages to biological-level dosages.
4. **Re-derive genotype tier** (AAAA/AAAB/AABB/AABC/ABCD) at biological level — do NOT use the `Genotype` column from the zygosity TSV directly.
5. Compute copy-count tuples for each individual.
6. Build `data/population.pkl` and downstream pickle.

Skipping step 4 is the most likely silent error — the protein-level tier in the zygosity file is plausible-looking but wrong for the model.

### D.3 ⚠️ Allele ID stability across runs

If the loader assigns integer IDs by sort order over observed alleles, adding new individuals (or filtering Ingroup differently) would shift IDs. Downstream pickles produced under different filters would be incompatible.

**Recommendation.** Use the `Allele_NNN` strings as primary keys throughout; build the integer mapping from a fixed canonical pool (all 58 alleles in `SRK_synonymy_groups.csv`, ordered alphanumerically). This makes IDs stable regardless of which subset of individuals is included.

### D.4 ⚠️ BL/EO ordering and palette must come from upstream constants

Hand-rolling a BL palette in the polyploid model will produce figures that mismatch Sven's reports. The first time he opens a polyploid-model figure expecting BL4 orange and sees BL4 blue, communication breaks down.

**Recommendation.** Either:
1. Import `srk_bl_constants` directly: `sys.path.insert(0, "../../SRK_bioinformatics/Scripts")` then `from srk_bl_constants import BL_ORDER, BL_COLORS, get_eo_order_within_bl`.
2. Mirror the module into `src/polyploid_utils.py` with a comment pointing to the canonical source.

Option 1 is fragile (cross-repo path). Option 2 needs a sync discipline.

### D.5 ⚠️ `simulate_generation` safety net masks loss dynamics

Lines 337–349 (`Step 4: Final safety net — verify no allele was lost`) auto-insert allele carriers if any allele was lost. This is opt-in (only when `preserve_rare=True`), but it means the headline "preservation strategy: 0 alleles lost" metric reported in `MODEL_OVERVIEW.md` is partly the safety net at work.

**Recommendation.** Make the safety net configurable (boolean parameter, default off). Run comparisons with safety-net-off to honestly report what the optimizer alone preserves vs. what the safety net rescues. Otherwise the strategy-comparison plots are not measuring what they claim.

### D.6 📌 The strategy-comparison palette should adopt the BL palette upstream

Existing comparison palette: Random (blue/gray), Optimized (red), Preservation (blue/green), Demographic (green diamond). This collides with the locked BL palette (BL1 purple, BL2 blue, BL3 red, BL4 orange, BL5 green). Strategy and BL on the same figure (e.g., per-BL convergence plots) would produce ambiguous coloring. Pick a non-overlapping strategy palette (e.g., grayscale or sequential viridis) to keep BL color reserved.

---

## Section E — Gaps with upstream framework

### E.1 🚨 TP1 framework not implemented

Upstream's central conservation diagnostic — J × P_compat with strict and leaky panels, mapped to MONITOR / AUGMENT / BIOBANK — is not in the polyploid model. The model reports χ², variance, and KL divergence (`distance_from_equilibrium`). These tell a related but distinct story. To produce LEPA-report-compatible recommendations, the polyploid model needs:

- `evenness_J(population)`: Shannon H / ln(k_observed)
- `p_compat(population, L=0.0)`: fraction of random pairs that are SI-compatible
- A quadrant assignment given J, P_compat, and the threshold pair (J=0.80, P_compat=0.40)

These are easy to compute (one or two new functions in `polyploid_utils.py`) and let the polyploid model produce TP1 scatters per simulation generation, directly comparable to Sven's Figures 14a–d.

### E.2 ⚠️ Leaky SI not modeled

Upstream P_compat is computed at L ∈ {0, 0.10, 0.25, 0.50}; empirical L̂ ≈ 0.18. The polyploid model assumes strict SI (L=0). This is tied to A.1 — adding a leakage parameter is the minimal change that makes the model consistent with upstream and partially captures C3.

### E.3 ⚠️ Class I/Class II SRK dominance ignored

Currently affects only Allele_055 (1/49 observed alleles), so the bias is small. But Sven is actively developing Step 12b; class assignments may expand. Worth tracking.

### E.4 ✅ HV functional groups vs Synonymy groups — model uses the right one

Synonymy groups (HV-identical) are the right unit for "same SI specificity". The model should NOT use the broader HV functional groups (which are clusters of HV-similar but not HV-identical alleles). This is correctly identified in the conventions captured in CLAUDE.md.

### E.5 ⚠️ Cross categories (Step 22b) — model has only "compatible / incompatible"

Upstream defines four cross categories:
- **Incompatible** — predicted no seed
- **Synonymy_test** — same HV functional group, different synonymy bin: experimental boundary test
- **Compatible_within** — different synonymy groups, same Class
- **Compatible_cross** — different Class

The polyploid model has binary compatible/incompatible. For *simulation* this is fine. For *management recommendation*, the categories matter — a `Synonymy_test` cross may or may not produce seed; advising a manager to attempt one is a different action from a Compatible_cross. Future work, not blocking the GFS refactor.

### E.6 🚨 No BL structure in the model

The model's allele pool, population, and crossing all operate on a flat list. The LEPA report's central conservation finding ("inter-BL allele transfers are the only recovery mechanism") cannot be evaluated unless the model carries:

- A `BL` attribute per individual (loadable from the join chain Individual → EO → BL)
- A within-BL vs across-BL crossing constraint
- BL-stratified reporting (per-BL allele pool, per-BL GFS, etc.)

This needs to be designed into the synonymy-aware loader from the start, not retrofitted later.

### E.7 ⚠️ EO76 SI-escape candidates violate the SI=functional assumption

7 SI-escape candidates exist species-wide; 5 are in EO76. The model assumes SI is functional for all individuals. For EO76 specifically, the assumption is violated for ~10% of individuals. The simulation might produce wrong conservation recommendations for EO76 if these individuals are present in the simulated population without their SI-escape status being marked.

**Recommendation.** Load and propagate `SI_functional_status` from `SRK_data_quality_categories.tsv` (in `tables/`). Either exclude SI-escape candidates from the breeding pool or simulate them with appropriate behavior. At minimum, flag them when reporting per-EO results.

---

## Prioritized recommendations (in order)

Before any code is written for the GFS refactor:

1. **Decide path 1 vs path 2 (Finding A.1 / C.1).** Are we doing instrumentation-only (track GFS, gate elitism), or are we doing mechanism-faithful (add competitive interaction or leakage to `is_compatible`)? Both are defensible; the choice determines whether the simulation can reproduce C3. **My recommendation:** add a leakage parameter (small change) AND keep the strict-SI rule as default. Run simulations at both L=0 and L=0.18 (empirical) and report both.

2. **Get Sven's answer on inheritance mode (Finding A.2).** This is a question only he and Buerki et al. 2026 can answer. The Crossing Plan should adopt his answer; if disomic, gamete formation needs reworking before any quantitative result is trustworthy.

3. **Implement BL structure in the loader (Finding E.6).** Add `bl` and `eo` attributes to the population data structure now, not later. Every downstream function that reports per-individual data benefits.

4. **Implement TP1 metrics (Finding E.1).** Two functions: `evenness_J` and `p_compat`. Cheap; gives LEPA-comparable output immediately.

5. **Recalibrate rare-allele threshold (Finding B.6).** For 27 biological alleles, 0.05 is too high. Switch to carrier-count or relative-to-target threshold.

6. **Make the simulation safety-net opt-in (Finding D.5).** Run strategy comparisons with safety-net off to measure the optimizer honestly.

Then proceed with the GFS instrumentation per the Plan (Phases 1–4), with the corrections above baked in.

---

## What the model gets right (steel-man check)

The critique above is pointed; some balance is in order.

- **The conceptual architecture is sound.** Separating crossing mechanics, allele frequency targets, optimization, and preservation into composable functions is a clean design.
- **The χ² minimization framing is mathematically correct** and the gradient is correctly derived.
- **The preservation strategy (elitism + mandatory rare crosses + optimizer penalty) is a thoughtful response to the rare-allele loss problem** — even where the implementation has bugs (elitism rewarding AAAA carriers), the underlying decomposition is right.
- **GFS as a metric is correct and matches upstream.**
- **The model correctly identifies that random mating loses rare alleles fastest** — this is consistent with NFDS theory under small Ne.
- **The demographic extension (logistic growth + Poisson noise + harmonic-mean Ne tracking)** is biologically motivated and correctly implemented.

The model is a strong starting point. The critiques above identify gaps that became visible only after surveying upstream — they are not "this is broken" critiques but "this is how to align with the broader framework" critiques. The Crossing Plan's GFS refactor is the right next step; the prioritized list above just sharpens the order.

---

## Sources for external claims

- [Inheritance and dominance of self-incompatibility alleles in polyploid Arabidopsis lyrata — Heredity](https://www.nature.com/articles/6800526)
- [Self-Incompatibility in the Brassicaceae: Receptor–Ligand Signaling — PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC151257/)
- [Pollen recognition and rejection during the sporophytic self-incompatibility response — PubMed](https://pubmed.ncbi.nlm.nih.gov/14659710/)
- [Autopolyploidy exacerbates dominance masking under negative frequency-dependent selection (Arabidopsis arenosa and lyrata) — bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.06.11.659117.full.pdf)
- [Self-incompatibility in Brassicaceae crops — PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4031107/)
- [Estimation of allele frequencies in polyploids under certain patterns of inheritance — Heredity](https://www.nature.com/articles/6800728)
- [Making a functional diploid: from polysomic to disomic inheritance — New Phytologist 2010](https://nph.onlinelibrary.wiley.com/doi/10.1111/j.1469-8137.2009.03117.x)
- LEPA_SRK_report.md, sections Q2, Q4, Q5, Q6 (in `~/Repos/SRK_bioinformatics/`)
- `srk_bl_constants.py`, `SRK_individual_GFS.R`, `srk_allele_hypotheses.py`, `srk_cross_plan.py` (in `~/Repos/SRK_bioinformatics/Scripts/`)
