# Analysis plan and deviations log (Dark Proteome x Tobacco Signature)

Документ собран 2026-10-31 как план и журнал отклонений. Основные правила,
границы и критерии ниже были зафиксированы до основных запусков (на это
указывают даты в разделах 7, 18 и 19), но документ в целом — рабочая запись, а
не независимая предрегистрация: внешней регистрации (OSF, AsPredicted) не было,
а разделы 18.8 и 19.8 фиксируют отклонения — включая случаи, когда результат
видели до тега (баг репетиции) и до смены источника данных. Ссылки вида
«section 18.3» из кода результатов ведут в этот документ. Хэши коммитов в
файлах `results/` относятся к локальной истории разработки, которой нет в
публичном репозитории.

# Pre-registration: Analysis Plan (Dark Proteome × Tobacco Mutational Signature)

**Status:** *(historical heading; the rules below were fixed before the runs —
see dated sections and the deviations log at 18.8/19.8; this document is not an
externally registered pre-registration.)*  
**Date:** 2026-10-25  
**Project:** Dark proteome and mutational footprint of tobacco smoke

---

## 1. Hypothesis

**H1 (working).** In tumors with stronger smoking mutational signature (SBS4), the fraction of mutations falling in intrinsically disordered regions (IDR) of proteins is higher than expected given the tumor's own 96-context mutation spectrum and the distribution of those contexts between dark and ordered protein positions. The effect is expected to be stronger in SCLC than in LUAD/LUSC.

**H0.** After correction for context composition, coverage, and selection, no systematic shift remains: observed/expected ratio of mutations in IDR is consistent with 1.

---

## 2. Primary definition of "darkness"

- **Primary:** AlphaFold pLDDT < 50 (very low confidence / likely disordered).
- **Sensitivity:** pLDDT < 30 and pLDDT < 70.
- The primary threshold is fixed before analysis; sensitivity thresholds are exploratory only.

---

## 3. Mutation-to-residue mapping (blocking step)

- Transcript: GENCODE v19, protein-coding CDS.
- Translate CDS, align to UniProt sequence.
- Keep only mutations where translated residue matches UniProt residue.
- Record and report the fraction of mutations discarded due to mismatch.
- Validate on known mutations (TP53 R175H, KRAS G12C, and others from our data where reference amino acid must match).
- Check reference allele matches hg19; amino acid matches UniProt.
- Do not silently fix mismatches; count them.
- **Do not examine is_dark vs SBS4 during this step.**

---

## 4. Samples and cohorts

- Primary: SCLC (George, JHU, CLCGP).
- Secondary: TCGA LUAD and LUSC (dose-response, non-smokers).
- Check sample overlap between George, CLCGP, JHU; remove duplicates.
- Minimum mutations per sample for signature estimation: determined by bootstrap stability.

---

## 5. Signatures

- Primary predictor: **SBS4 fraction** (SBS4 / total mutations), standardized (1 SD).
- Fitting set: SBS1, SBS2, SBS4, SBS5, SBS13, SBS17, SBS40, **plus SBS31, SBS35** (platinum; amendment 2026-10-28 — treated SCLC may carry platinum signatures).
- All fits on the **coding observation space** with opportunity normalization
  (amendment 2026-10-28: original text said "WGS (George): fit on whole genome";
  only coding calls are available — "coding variants from a WGS study").
- Report bootstrap CI for SBS4.
- SBS40 is exploratory only (COSMIC: etiology unknown).

---

## 6. Primary model (fixed)

**Level: sample (Poisson with offset).** *(Amendment 2026-10-28: the original
pre-registration specified a mutation-level binomial logit with OR. Replaced by
the sample-level Poisson with offset — one single primary model; see §16.)*

O ~ offset(log E) + β1 · SBS4_std + β2 · TMB_std + β3 · wes

- O = observed dark mutations in a sample; E = expected from the tumor's own
  96-context spectrum and the position-level darkness map.
- SBS4_std = z-scored SBS4 fraction; wes = JHU+CLCGP (WES) vs George.
- **Main coefficient:** β1 — effect of 1 SD of SBS4 fraction on log O/E.
- **Primary effect measure:** exp(β1) = O/E ratio per 1 SD SBS4.
- **Primary inference:** NB2 (Var = μ + θ·μ²) with t-based CI (df = n − p) —
  dispersion-adjusted SE. Sensitivity: quasi-Poisson (Pearson χ²/df) and HC1;
  robustness variant: Poisson GEE clustered by sample.
  *(Amendment 2026-10-28, pre-run: HC1-only and quasi-Poisson failed the
  calibration criteria in §18 on simulated data; NB2-t passed.)*
- **Boundary rule (fixed 2026-10-28, pre-run):** if the NB2 dispersion estimate
  is at the boundary (alpha ≤ 1e-4) or its SEs are non-finite, the Poisson-t
  limit of NB2-t is used and labelled as such in the output.
- cancer_type and SBS4_std×cancer_type interaction enter **only after TCGA
  inclusion** (secondary analysis).
- Synonymous mutations as negative control.
- SBS1 (age) as negative control: should not predict IDR shift.

---

## 7. Equivalence boundary (fixed before analysis)

- **Minimum interesting effect:** O/E ratio = 1.2 per 1 SD SBS4.
  *(Rationale note 2026-10-31: 1.2 was chosen a priori as the smallest effect
  judged biologically interesting for a positional bias — not as a detection
  limit; the realized design constrains effects to a few percent per SD on the
  measured-predictor scale.)*
- **Equivalence interval:** (1/1.2, 1.2) = (0.833, 1.2).

**Interpretation rules:**
- If 95% CI for OR lies entirely within (0.833, 1.2) → **informative null** (effect smaller than of interest).
- If 95% CI is wider and crosses boundary → **inconclusive** (cannot exclude meaningful effect).
- If CI excludes 1 and lies outside interval → **effect detected**.

This prevents calling every narrow CI "inconclusive" and every wide CI "null".

---

## 8. O/E calculation

- Observed: mutations in IDR.
- Expected: based on tumor's own 96-context spectrum and position-level context/darkness map.
- Metric: log(O/E).
- Absolute position counts (not gene fractions).

---

## 9. Covariates (after primary analysis)

- Gene expression.
- Replication timing.
- (Chromatin deferred.)

---

## 10. Exclusion criteria (fixed)

- **Cell lines excluded** (JHU: 12 cell lines per clinical annotation SAMPLE_CLASS; tumors only, 30 samples). *(Amendment 2026-10-28: added explicitly after cell lines were found in the analysis set.)*
- Samples with fewer than N mutations (N set by bootstrap stability).
- Duplicate samples across cohorts (remove, record count).
- Mutations not mappable to protein residue (record discard fraction).
- If SBS4 is unstable (bootstrap CI includes 0 or very wide), record as limitation and proceed with caveat.

---

## 11. Sensitivity analyses

- Thresholds: pLDDT < 30, < 50, < 70.
- Signature set variations.
- Per-cohort analysis.

---

## 12. Interpretation rules

- No causal claims.
- If SBS4 CI includes zero but is narrow → informative null (per equivalence boundary).
- If SBS4 CI is wide → inconclusive.
- Report power analysis before collecting more data.

---

## 13. Synthetic data policy

- All synthetic/test files removed from analysis directory.
  *(Note 2026-10-31: synthetic **data** files are not used in the analysis; the
  generator/simulation scripts (make_synthetic_oe.py, sim_*.py) remain in the
  repository as tools.)*
- No synthetic data used in final results.

---

## 14. Commit

- This file committed to git with date before any real analysis.
- Optionally uploaded to OSF/GitHub release for public timestamp.


---

## 15. Duplicate sample rule (fixed before application)

- Compare mutation profiles between any two samples across cohorts.
- Compute **overlap coefficient** = |intersection| / min(|A|, |B|).
- If overlap coefficient ≥ 0.5 → treat as duplicate.
- Rule: keep the sample from the cohort with **larger mutation set** (George WGS > WES cohorts). If equal, keep George.
- Record all removed duplicates and their overlap coefficients.
- Apply the same rule to all pairs, including within-cohort.

---

## 16. Model deviation note

- Original pre-registration specified GLMM with random intercept on cohort.
- With only 3 cohorts, cohort variance is not estimable.
- **Deviation:** cohort as fixed effect; clustering by sample via GEE with exchangeable correlation, or GLMM with random intercept on sample only.
- This deviation is recorded here before real data analysis.
- Main coefficient unchanged: log O/E ratio per 1 SD SBS4 fraction (amended from OR per §6 amendment).

---

## 17. Disclosures (dated record of what was run and seen)

- **2026-10-28 (exploratory, pre-tag).** Model smoke tests ran on real George
  mutation labels (107 samples) with a *synthetic random* background; printed
  coefficients were seen (e.g. SBS4 O/E ratio ≈ 1.10 on synthetic background).
  These numbers carry no scientific meaning (background was uniform random) but
  are recorded here as unblinded events. A permutation smoke test on the same
  synthetic background was also seen. Annotate-validation counts on George
  (34583 → 29936 kept) were seen.
- No model has been run on real darkness background / real signatures.
  The single confirmatory run happens after the pre-registered simulations and
  a git tag (§14).

---

## 18. Confirmatory run protocol (fixed 2026-10-28, before the real run)

### 18.1 Calibration criteria (fixed before viewing results)
Simulated data, 500 runs per scenario, n = 143, overdispersion θ ∈ {0.15, 0.30}
(Gamma frailty, not modeled by the DGP):
- type-I error in **(3.5%, 6.5%)** and 95% CI coverage in **(93%, 97%)** under H0.
- Achieved by the primary method (NB2-t): θ=0.15 → type-I 5.8%, coverage 94.2%;
  θ=0.30 → type-I 5.0%, coverage 95.0%. Power for O/E = 1.2: 99.8% / 95.0%.
  (HC1: 7.4%/7.0% type-I — out of corridor; quasi-Poisson: 8.3%/8.7% — out.)

### 18.2 Run contents (one run, one output file)
1. Primary model (§6) on all 143 samples.
2. Permutation test, 10000 permutations (positions reshuffled within tumor,
   mutation types preserved).
   *(Terminology corrected 2026-10-31: implemented as a parametric
   background-resampling test — binomial draws from the background per context,
   see permutation_test.py — not a label permutation; p at the floor of
   1/(perms+1) is a bound, not an estimate.)*
3. Negative controls, same function and code path:
   a. SBS1 (age) as predictor instead of SBS4;
   b. synonymous (Silent) mutations as outcome instead of all coding.
4. Output saved to a single file with the git commit hash and SHA-256
   checksums of all input files (darkness tracks, background, opportunity,
   mutation tables, signature table).

### 18.3 Interpretation rules (fixed in advance)
- CI entirely within (0.833, 1.2) → **informative null**.
- CI crosses the boundary / too wide → **inconclusive**.
- CI excludes 1 and extends outside the interval → **effect detected**.
- **A failed negative control** (e.g. SBS1 "predicts" a shift, or the
  synonymous outcome shows a shift) → the result is recorded as
  **"under question"** regardless of the primary outcome. It does NOT license
  trying another model; any follow-up analysis is a dated deviation (§18.4).

### 18.4 After unblinding
- Any change to the model, thresholds, or mutation set is recorded as a
  separate dated deviation; the first result stays in the report unchanged.
- The confirmatory run is performed **once**.

### 18.5 Secondary analyses (fixed 2026-10-28, before the real run)
1. **Discard-bias sensitivity.** The per-sample fraction of mutations discarded
   by residue validation enters as a covariate (and/or mutations weighted by
   inverse probability of retention). If the conclusion matches the primary
   analysis, mapping-discard bias is not driving the result. (Context
   composition is already handled: E uses the retained mutations' contexts.)
2. **Overdispersion robustness.** Calibration is re-checked under lognormal
   and NB1-type (Var = μ·(1+θ·μ), non-quadratic) overdispersion. If type-I
   stays within (3.5%, 6.5%), the method is robust to overdispersion shape;
   otherwise this is recorded as a known limitation before the real run.
3. **Power caveat.** At θ = 0.30 power for O/E = 1.2 is 95%; if the real
   dispersion is higher, the "inconclusive" outcome is expected a priori.

### 18.6 Known limitations (fixed before the real run)
- The E background covers 86% of coding sites (positions without pLDDT
  annotation or unmapped transcripts are excluded); per-channel coverage is
  stable (0.861 ± 0.007 across 96 channels), i.e. channel-neutral.

### 18.7 Technical restart rule
- If the run script crashes before producing its output file, it may be fixed
  and restarted (a "technical restart"); the result of any run that completed
  is preserved and never replaced.

### 18.8 Deviations log (dated)
- **2026-10-28 — rehearsal shuffle bug (protocol deviation).** In the rehearsal
  run the `--shuffle-sbs4` flag was not passed to the primary and control model
  calls (argument omission in the script), so the rehearsal printed the REAL
  coefficients instead of shuffled ones; the real result was thereby seen
  before the pre-run tag. The bug was fixed and the fact is recorded here. The
  tagged confirmatory run was nonetheless performed exactly once (commit
  45c7c25, tag confirmatory-20261028) and its output file is preserved
  unchanged. No analysis decisions, interpretation rules or thresholds were
  altered after seeing the result (they were fixed in §6 and §18.1-18.3 before
  the rehearsal).

---

## 19. LUAD dose-response analysis (fixed 2026-10-28, before the run)

**Question (within one cancer type):** in lung adenocarcinoma, does the share
of mutations in disordered protein regions change with the strength of the
smoking signature, across its full range (never-smokers to heavy smokers)?

### 19.1 Population and data
- TCGA LUAD, primary tumors, coding SNPs, GRCh37; same mutation pipeline and
  residue validation as SCLC (match flags enforced, UniProt position with
  offset); darkness = AlphaFold pLDDT < 50 (proteome-wide track).
- **Mutation source (amendment 2026-10-28, before the completed run):**
  Firehose Legacy (`luad_tcga`) MAF instead of PanCanAtlas — the PanCanAtlas
  mutation profile contains no synonymous mutations at all, making the
  pre-registered Silent control (19.4) impossible. Firehose covers the same
  TCGA patients and the same assembly. Rules, controls and thresholds
  unchanged; see the deviations log (19.8).
- Smoking data from the Firehose Legacy clinical table for the same patients
  (TOBACCO_SMOKING_HISTORY_INDICATOR, SMOKING_PACK_YEARS).

### 19.2 Blind scale check (done before any outcome analysis) — PASSED
- On the amended (Firehose) set of record: never-smokers median SBS4 = 0.000
  (n = 32), current 0.401 (n = 45), reformed >15 y 0.293, reformed <=15 y 0.387;
  Spearman(status, SBS4) = 0.314 (p = 2e-6); Spearman(pack-years, SBS4) = 0.223
  (p = 0.005, n = 158). Firehose MAF covers 230 patients (all with smoking
  info); hypermutator rule flags 0 samples (max 1574 vs median 194).
- (Earlier PanCanAtlas check, superseded set: never 0.011 vs current 0.429.)
- Scale ordering (never < former < current) holds; proceeding.

### 19.3 Exclusion (fixed before the run)
- Hypermutated samples: coding SNPs > 10x cohort median (repair-deficiency
  proxy). Preview on counts only: 1 sample flagged (max 2020 vs median 185).

### 19.4 Primary model
- Same as SCLC (section 6): O ~ offset(log E) + SBS4_std + TMB_std;
  NB2-t with the pre-registered alpha->0 Poisson-t boundary rule;
  equivalence bounds (0.833, 1.2); same interpretation rules (18.3) and
  control-failure rule. Cancer type is constant (LUAD) and does not enter.
- Controls (same function): SBS1 as predictor; Silent mutations as outcome.
- Permutation test, 10000 permutations.

### 19.5 Secondary analyses (fixed before the run)
1. Pack-years (standardized, subset with data) as an independent smoking
   measure instead of SBS4 (removes the shared-mutations confound).
2. Smoking-status contrast: never (1) vs current (2) with SBS4 as predictor.
3. Mapping-discard fraction as covariate (as in 18.5.1).

### 19.6 Background and limitations
- Expected counts from the coding-space background (same builder as SCLC,
  channel-neutral, ratio verified); exact WES capture territory of the TCGA
  panel is unknown and recorded as a limitation.
- Pack-years missing for ~1/3 of patients; status categories cover 509/561.

### 19.7 Run discipline
- One run, one output file with git commit hash and SHA-256 of all inputs
  (run_luad.py; same code path / functions as the SCLC run).
- After unblinding: same rules as 18.4, 18.7, 18.8 (deviations dated, first
  result preserved, technical restarts allowed only before completion).
- The SCLC confirmatory result is untouched; no cross-cancer (LUAD vs SCLC)
  claims from this run — that would need a separate pre-registered comparison.
- Interpretation: the result speaks to variation WITHIN adenocarcinoma across
  the smoking range; it does not by itself establish a smokers-vs-non-smokers
  causal effect.

### 19.8 Deviations log (dated)
- **2026-10-28 — technical restart + input amendment.** The first LUAD run
  attempt crashed at CONTROL B before producing its output file: the
  PanCanAtlas mutation profile contains no synonymous mutations, so the
  pre-registered Silent control could not run. Values printed before the crash
  are disclosed here: PRIMARY O/E per 1 SD = 0.9885 (CI 0.9764..1.0009),
  CONTROL A = 0.9888 (CI 0.9643..1.0138); the hypermutator exclusion flagged
  1 sample (TCGA-17-Z031-01). Per 18.7 this is a technical restart. The
  mutation source was amended to Firehose Legacy (19.1), all pipeline stages
  were rebuilt on it, and the rules, controls and thresholds remained
  unchanged. The completed run (single, tagged) is the one of record.
