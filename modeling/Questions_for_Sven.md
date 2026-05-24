# Open Questions for Sven — Polyploid Model Refactor

**From:** Jim Beck
**Date:** 2026-05-23
**Context:** We are refactoring the `polyploid-model` Python project to operate on the new biological-allele representation (synonymy-collapsed) and to align with the LEPA report's GFS / TP1 / TP2 framework. A critical review surfaced five questions that I cannot resolve from the data alone. Each is framed below with a default assumption we'll proceed under if you don't have a definitive answer yet — but your input would let us replace defaults with calibrated parameters.

---

## Q1 — What is the inheritance mode at the LEPA S-locus?

**Question.** Is the S-locus in *Lepidium papilliferum* tetrasomic (random pairing across all 4 chromosomes → C(4,2)=6 uniform gametes) or disomic (two homeologous subgenomes segregating independently → 4 gametes per individual, structured by subgenome assignment)?

**Why it matters.** Most allopolyploids show disomic inheritance because homeologous subgenomes pair preferentially with themselves. The polyploid crossing model currently assumes tetrasomic via `form_gametes` returning all `C(4,2)=6` allele pairs uniformly. Under disomic inheritance the same genotype produces a different gamete distribution — for an AABB individual, the gamete frequencies change from `AA 1/6, AB 4/6, BB 1/6` (tetrasomic) to a partition-dependent distribution that averages `AA 1/8, AB 6/8, BB 1/8` (more heterozygous) under uniform prior. This propagates to GFS values (±50% swing per individual), C3 cascade rate, and offspring distributions from every cross.

**Has Buerki et al. 2026 established this?** You cite "LEPA's confirmed allopolyploid origin (Buerki et al., 2026)" in Q6 Stage 5 of the LEPA report, and flag the tetrasomic vs disomic uncertainty as a caveat in Q4 Layer 2 (line 308). Is the inheritance mode itself addressed in that work or in a follow-up?

**Our default if unknown.** We will implement *both* modes (`inheritance_mode='tetrasomic'` and `inheritance_mode='disomic'`) and run a sensitivity analysis (notebook `09b_inheritance_sensitivity.ipynb`) comparing strategy rankings and per-individual GFS distributions under each. If strategy ranking is preserved across modes, we proceed with tetrasomic as the primary report and disomic as a bracketed sensitivity. If ranking changes, we ask you to resolve this before we publish.

---

## Q2 — Can the upstream pipeline output per-subgenome haplotype assignments?

**Question.** The current `SRK_individual_zygosity.tsv` `Allele_composition` column tells us, for example, that `Library001_barcode02` carries `Allele_018(1)+Allele_050(1)+Allele_049(2)` — but not which alleles live on which homeologous subgenome. Can the pipeline (Step 11–12) be extended to output a per-subgenome assignment? E.g., a column `Subgenome_assignment = "S1:Allele_018+Allele_050; S2:Allele_049+Allele_049"`.

**Why it matters.** Even if you confirm disomic inheritance, the polyploid model cannot simulate disomic gamete formation correctly without knowing which alleles are on which subgenome for each individual. Without phased per-subgenome data we are forced to average over plausible configurations, which dilutes the predicted gamete distribution.

**Our default if unavailable.** Average over distinct subgenome partitions under uniform prior (the disomic-averaged mode described in Q1). Quantitatively biased toward heterozygous gametes; aggregate strategy recommendations likely still valid.

---

## Q3 — Is the empirical leakage L̂ ≈ 0.18 a mix of mechanisms, or specifically Lewis (1947) competitive interaction?

**Question.** The LEPA report Q6 Stage 5 invokes the Lewis (1947) competitive interaction model — heterozygous (AB) pollen escapes rejection on a mother carrying both A and B because the two SCRs compete and weaken the rejection signal. You quantify the aggregate rate as L̂ ≈ 0.18 from observed AAAA proportions. Is your interpretation that the entire L̂ is Lewis-style competitive interaction, or does L̂ aggregate multiple mechanisms (competitive interaction + partial LoF + environmental SI failure + chimeric assembly recognition errors)?

**Why it matters.** When I traced the dynamics of pure competitive-interaction-only SI under selfing, I found that AABB and AAAB selfing converge to a heterozygous mix (AABB ↔ AAAB), NOT to AAAA fixation. So competitive interaction alone is insufficient to explain the observed 60% AAAA. Either there's a mechanism I'm missing, or the empirical L̂ includes other leakage modes. Knowing which one drives our implementation choice: a clean competitive-interaction rule in `is_compatible`, vs an aggregate stochastic leakage parameter applied to all rejected gametes.

**Our default.** Implement L as an aggregate stochastic leakage parameter (Path 2a in the critical review) with `L=0.18` as the empirical default. This captures the observed AAAA rate regardless of mechanism mix. We will *not* implement competitive interaction as a separate biochemical rule.

---

## Q4 — Are the 7 SI-escape candidates phenotypically confirmed?

**Question.** `SRK_SI_escape_candidates.csv` flags 7 ingroup individuals for controlled-selfing tests, 5 of which cluster in EO76 (BL3). Have these tests been performed yet? If so, are any confirmed self-compatible (SC)?

**Why it matters.** The polyploid model assumes SI is functional for all individuals. If the EO76 cluster is confirmed SC, the model's recommendations for EO76 are based on a violated assumption — those individuals shouldn't be treated as breeding-pool participants under SI logic. The implementation question is whether to exclude them, model them with `leakage=1.0`, or treat them as a separate management category.

**Our default.** Load `SI_functional_status` from `SRK_data_quality_categories.tsv` and propagate per individual. Default behavior: include them in the breeding pool but flag in per-EO output. We will not change simulation logic for them until you confirm SC.

---

## Q5 — How should AAAB individuals be treated as pollen donors?

**Question.** The Crossing Plan's preferred recommendation (§4.1, option b) is "AAAB never used as pollen donor in managed crosses" — based on the reasoning that 3/6 of AAAB pollen is AA homozygous and rejected wherever A is in the mother, so the remaining 3/6 AB pollen are functionally the only useful gametes. The proposed rule prevents the AAAB → AAAA slide that drives Stage 5.

Under the strict-SI rule the polyploid model uses, this is conservative-but-defensible. Under realistic SI leakage (L̂≈0.18), AAAB pollen has access to additional inter-EO crosses where the mother doesn't carry the A allele — those crosses produce useful heterozygous offspring (e.g., AAAB × CCDD → AABC offspring, a one-generation rescue).

**Two implementation options:**
- (a) Exclude AAAB from pollen donor role entirely (Crossing Plan recommendation)
- (b) Allow AAAB as pollen donor only in inter-EO crosses where the partner-EO lacks the A allele

**Why it matters.** Affects which individuals are usable as pollen donors in the optimizer. Under (a), AAAB carriers of rare alleles in BL2/BL4 (the diversity reservoirs) are sidelined — potentially restricting recovery options. Under (b), the model is more permissive but biologically defensible.

**Our default.** Implement (b) for inter-EO crosses, (a) for within-EO crosses. Report a sensitivity comparison.

---

## Summary table

| Q | Topic | Blocks work? | Default if no answer | Best-case answer enables |
|---|---|---|---|---|
| Q1 | Inheritance mode | No (both implemented) | Tetrasomic as primary; disomic as sensitivity | Single-mode results without bracketing |
| Q2 | Subgenome haplotypes | No (average over configurations) | Uniform-prior averaging | Direct disomic simulation |
| Q3 | Leakage mechanism mix | No | Aggregate stochastic leakage L=0.18 | Mechanistic competitive-interaction rule |
| Q4 | SI-escape candidates | No (flag and include) | Include with annotation | Exclude / separate management category |
| Q5 | AAAB donor rule | No (option b default) | Inter-EO permissive, within-EO conservative | Single rule with biological backing |

---

## What we're doing in the meantime

Refactor is proceeding on the defaults above. The work plan (`Critical_Review.md` and committed task list):

1. Synonymy-aware data loader → `data/population.pkl` with biological alleles
2. BL annotations + TP1 metrics in `polyploid_utils.py`
3. `inheritance_mode` + `leakage` parameters threaded through gamete formation and SI checks
4. GFS-aware optimizer with critical-review fixes (rare threshold, safety net, elitism)
5. Inheritance-mode sensitivity comparison notebook
6. Final GFS-quality notebook with all four strategies + inter-EO rescue mode

Estimated effort: ~5 days of focused work. Your answers to Q1, Q3, and Q5 would tighten the result substantially; Q4 affects only EO76-specific recommendations.

Happy to discuss any of these by email or over a call.
