# Genotype-Quality-Aware Crossing Strategy — Plan

**Author:** Jim Beck (with Claude Opus 4.7)
**Date:** 2026-04-30
**Source:** Response to *SRK-Based Assessment of Self-Incompatibility in LEPA* (Sven Buerki et al.), specifically the Compatibility Collapse Cascade (C3) framework in Section 5 and the GFS / TP2 results in Section 4b.

---

## 1. Why this plan

The current `polyploid-model` preservation strategy is **allele-frequency-centric**: it minimises χ² distance from uniform allele frequencies, with rare-allele penalties to prevent extinction. It works on the TP1 axis (allele richness and evenness).

It is **blind to the TP2 axis** (genotype dosage / GFS / %AAAA). In particular:

- `select_elites` ranks individuals by `Σ 1/freq(allele)`. An AAAA carrying a rare allele X scores `4 × 1/freq(X)` — top of the ranking — and is selected as elite. The model preserves the *worst possible carrier* of our rarest alleles.
- `compute_optimal_weights` minimises χ² on allele frequencies. A cross of AAAA × ABCD that "looks good" on allele balance still locks the offspring lineage into ≥2 copies of A (since AAAA can only contribute AA gametes), capping that lineage's reachable genotype tier at AABB.
- `get_mandatory_rare_crosses` picks the "best cross" per rare allele on a metric that is not weighted toward offspring genotype quality.

Net effect: under the current preservation strategy, allele-frequency variance can converge to zero while %AAAA drifts flat or upward. The TP2 axis stays critical even when TP1 looks healed. Without addressing the dosage problem, the C3 cascade cannot be reversed.

---

## 2. Biological recap (from the LEPA report)

### 2.1 The competitive interaction model (C3 Stage 5)

A tetraploid AABB plant produces diploid pollen by sampling 2 of 4 allele copies. The gamete distribution is:

| Gamete | Probability | Notes |
|---|---|---|
| AA | 1/6 | Rejected by any pistil expressing SRK-A |
| AB | 4/6 | **Bypasses SI when both A and B are on the maternal pistil** |
| BB | 1/6 | Rejected by any pistil expressing SRK-B |

When AB pollen lands on an AABB pistil, both SCR-A and SCR-B compete for SRK binding sites simultaneously. The competing recognition signal is too weak to trigger rejection. **Self-fertilisation succeeds.** The resulting offspring are predominantly AAAB (e.g., AA egg × AB pollen).

AAAB plants then produce AA pollen 3/6 of the time (vs. 1/6 from AABB), making subsequent SI bypass even more likely. Iterate over generations → AAAA. This is how a population accumulates 56% homozygous dead-ends *while the SRK gene itself remains functional*.

### 2.2 Genotype tier ladder

| Genotype | GFS | Diploid gametes (probabilities) | Useful as pollen donor? | Useful as maternal parent? |
|---|---|---|---|---|
| AAAA | 0.000 | AA (6/6) | **Never** — pollen rejected wherever A is present | Limited — only against partners carrying ≥2 non-A alleles in one gamete; offspring locked at ≥2 copies of A |
| AAAB | 0.500 | AA (3/6), AB (3/6) | Rarely — both gamete types contain A | Limited — offspring inherit 2 A's from her |
| AABB | 0.667 | AA (1/6), AB (4/6), BB (1/6) | **AB pollen is the SI-bypass driver of selfing collapse** | Good |
| AABC | 0.833 | 5/6 heterozygous | **Excellent** — multiple non-A pollen options | **Excellent** |
| ABCD | 1.000 | All 6 heterozygous | Ideal | Ideal (none in dataset) |

In the LEPA dataset: 56% of all 189 individuals are AAAA, 21% AAAB, 21% AABB, 3% AABC, 0% ABCD. Only **5 AABC individuals exist** — all in EO67 and EO70 — and the report explicitly identifies them as "the only endogenous genetic resource capable of reversing the C3 cascade" (Section 5, conservation implication).

### 2.3 The asymmetry to exploit

AAAA individuals are **not equally useless in both roles**:

- As **pollen donor**: useless. Only AA pollen, rejected wherever A is on the pistil — i.e., almost every partner.
- As **maternal**: salvageable in narrow circumstances. Successful crosses with AABC or ABCD partners can produce AABC offspring (e.g., AAAA × AABC where the BC pollen succeeds yields AABC offspring). This converts an AAAA dead-end lineage into an AABC offspring lineage in one generation.

This asymmetry is the lever for "rescue crosses" — Phase 3 below.

---

## 3. Recommended plan

### Phase 1 — Instrumentation (no behaviour change)

Add to `src/polyploid_utils.py`:

1. `gfs(genotype: tuple) -> float` — implements `1 - Σ n_k(n_k-1) / 12` (formula on line 194 of the report).
2. `genotype_class(genotype: tuple) -> str` — returns one of `"AAAA"`, `"AAAB"`, `"AABB"`, `"AABC"`, `"ABCD"`.
3. `expected_offspring_gfs(p1: tuple, p2: tuple) -> float` — enumerates the 6×6 gamete combinations, filters to SI-compatible, weights by probability, returns mean GFS of the resulting offspring distribution.
4. Track `prop_AAAA(population)` and `mean_gfs(population)` each generation in `simulate_generation`.

### Phase 2 — Augment the optimisation objective

5. Replace the χ²-only objective in `compute_optimal_weights` with a **two-term loss**:
   - **Allele-frequency term** (existing): χ² distance from uniform.
   - **Genotype-quality term** (new): penalty proportional to `mean(1 − expected_offspring_gfs)` across the chosen crosses, weighted by a new tunable `gfs_weight` (suggest starting at 1.0 — same order as existing `preservation_weight`).
6. Modify `select_elites` to multiply rarity score by GFS:
   - AAAA rare-allele carriers score 0 (excluded from elitism).
   - AABB+ rare-allele carriers score fully.
   - AAAB carriers score at half weight.

   Rationale: elitism should preserve *breeding-capable* rare-allele carriers, not warehouses of dead alleles.
7. Modify `get_mandatory_rare_crosses` so the "best cross" for each rare allele requires the partner to be ≥AABB. AAAA-carrying-rare-allele individuals can still be *used*, but only as the maternal parent in a cross where the donor is AABC or ABCD.

### Phase 3 — New strategic mode: "Rescue crosses"

8. Add an `inter_eo_rescue_crosses` helper that pairs each of the 5 AABC individuals (EO67, EO70) against AAAA individuals in *other* EOs that lack allele B/C/D. Each such cross converts an AAAA dead-end lineage into an AABC offspring lineage in one generation. This is the only mechanism the LEPA report identifies for reversing C3 Stage 5.
9. Add a new simulation strategy `"gfs_preservation"` alongside Random / Optimized / Preservation / Demographic. Visualisations should report both:
   - existing TP1 metrics (allele-frequency variance, χ², extinct alleles)
   - new TP2 metrics (%AAAA trajectory, mean GFS trajectory)

### Phase 4 — Validation

10. **EO27 stress test.** EO27 is the most TP2-degraded (62% AAAA, mean GFS = 0.224). The new strategy should achieve comparable variance reduction to the current preservation strategy *and* a measurable drop in %AAAA after 5 generations.
11. **EO70 lever test.** EO70 retains only 6 alleles but contains 1 of the 5 species-wide AABC individuals. The optimiser should heavily prioritise crosses involving that individual.
12. **Synthetic regression test.** Construct a toy population where a rare allele is held only by an AAAA individual. The frequency-only strategy will mandate crosses using that individual; the GFS-aware strategy should either route around it or rescue it via an AABC partner first. Confirm divergent behaviour as a regression check that the genotype-quality term is actually firing.

---

## 4. Open design decisions

### 4.1 How to treat AAAB

Two options:

- **(a) Symmetric with AABB+:** penalise but allow in either role.
- **(b) Maternal-only:** AAAB never used as pollen donor (3/6 of its pollen is AA-homozygous and rejected wherever A is present anyway). This more aggressively prevents the AAAB → AAAA slide that drives Stage 5.

**Recommendation: (b)**, since it directly targets the runaway feedback loop the report identifies.

### 4.2 Within-EO vs. inter-EO

The report (Section 5 and Section 7) explicitly identifies inter-EO transfers as essential — particularly for AABC seed parents that are the only source of B/C/D alleles for AAAA-dominated EOs. But within-EO and inter-EO crosses are biologically and operationally distinct (the latter requires hand-pollination or seed transfer, not natural pollinator visits at <500 m).

**Recommendation:** Keep within-EO and inter-EO as **separate, toggleable strategies** in the simulation (`inter_eo_rescue=True/False`), defaulting to within-EO so existing comparisons stay clean. Sven's input on which inter-EO pairings are operationally feasible would shape Phase 3.

### 4.3 Where to land the code

- **Option A:** Extend `src/polyploid_utils.py` and add a new notebook `notebooks/09_genotype_quality.ipynb` that mirrors the structure of the existing `04_optimize` and `05_visualize` notebooks but on the new objective.
- **Option B:** Fold into Sven's new `notebooks/08` if his recent additions already have scaffolding for GFS / TP2 plots.

**Recommendation:** Need to inspect Sven's notebook 08 first — Option B is cleaner if his scaffolding exists; Option A is the fallback.

---

## 5. Open questions for Sven

1. **AAAB treatment.** Do you agree with (b) — AAAB never used as pollen donor in managed crosses? Or is the AB-pollen subset (3/6) worth retaining for inter-EO transfers where the maternal partner lacks both A and B?
2. **AABC seed parents — operational constraints.** The 5 AABC individuals (Library007_barcode82, Library007_barcode86 in EO67; Library008_barcode22 in EO70; the other two are the AABB-tier top-ranked at GFS = 0.667, listed in your seed-priority table) are the rescue lever. Are all five accessible for managed crossing? Are seed/pollen transfers between EO67 ↔ EO70, or from these EOs to EO25/27/76, operationally viable?
3. **GFS weight calibration.** The `gfs_weight` tunable controls the trade-off between allele-frequency convergence (TP1) and genotype-quality recovery (TP2). At very high weight the model will refuse to use AAAA individuals at all, even when they hold the only copy of a rare allele — sacrificing TP1 for TP2. Do you have intuition on where the right balance is, or should we parameter-sweep and let the data choose?
4. **Self-compatibility mutants.** The five individuals with non-functional SRK sequences (Section 5 baseline) — are these the same five AABC individuals, or distinct? If distinct, should we exclude SC mutants from the breeding pool, or treat them as a separate category?
5. **Notebook scaffolding.** Does your new `notebooks/08` already include GFS/TP2 plotting that we should build on, or is this independent work?

---

## 6. Summary

| Lever | Current model | Proposed model |
|---|---|---|
| Optimise allele frequencies (TP1) | ✅ χ² minimisation | ✅ retained |
| Penalise rare-allele extinction | ✅ preservation_weight | ✅ retained |
| Account for genotype dosage (GFS) | ❌ blind | ✅ `gfs_weight` term in objective |
| Elitism preserves rare-allele carriers | ✅ but selects AAAA | ✅ filtered to ≥AABB |
| Mandatory rare-allele crosses | ✅ but ignores partner GFS | ✅ partner must be ≥AABB |
| Use AABC seed parents to rescue AAAA lineages | ❌ | ✅ Phase 3 inter-EO rescue mode |
| Track %AAAA / mean GFS over generations | ❌ | ✅ TP2 metrics in dashboard |

The result should be a strategy that converges TP1 *and* TP2 simultaneously — moving populations toward both equal allele frequencies and a higher fraction of AABB+ individuals — instead of the current strategy that addresses only TP1 and leaves the C3 cascade in place.
