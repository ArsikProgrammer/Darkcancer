# Analysis plan and deviations log (Dark Proteome x Tobacco Signature)

This document was assembled on 2026-10-31 as an analysis plan and deviations log. The main rules, boundaries, and criteria below were fixed prior to the main runs (as indicated by the dates in sections 7, 18, and 19). However, the document as a whole is a working record, rather than an independent pre-registration: there was no external registration (OSF, AsPredicted), and sections 18.8 and 19.8 record deviations—including cases where the result was seen prior to tagging (rehearsal bug) and before the change of the data source. References such as "section 18.3" from the results code point to this document. Commit hashes in the results/ files refer to the local development history, which is not present in the public repository.

# Original analysis plan text (working title: "Pre-registration")

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
  *(Rationale note 2026-10-31: the threshold itself was fixed in this plan
  before the runs; this rationale was formulated afterwards. Honest wording:
  1.2 is the threshold of biological interest chosen at the time; the earlier
  "detection limit" explanation was wrong — the realized design constrains
  effects to a few percent per SD on the measured-predictor scale.)*
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
  *(Note 2026-10-31: this held in the development history but NOT in the public
  repository, where the plan was committed on 2026-10-31, after the runs;
  the dated sections and the deviations log are the honest record of what was
  fixed when.)*
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
  *(Clarification 2026-10-31: the medians differ by mutation set — 185 on the
  PanCanAtlas preview, 194 on the Firehose scale check, 150 in the run output
  (computed on validated mutations). The rule itself is unchanged.)*

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

---

## 20. Disorder-based darkness re-analysis (fixed 2026-11-24, before the runs)

**Motivation (external review).** AlphaFold pLDDT is a confidence measure, not
a disorder predictor, and AlphaFold is known to produce plausible-looking
structure where none exists ("hallucination"). An external reviewer recommended
repeating the analysis with specialized sequence-based disorder predictors.
The §1-19 results stand unchanged as the runs of record; this section defines
the planned re-analysis.

- **Predictor:** metapredict (pinned in requirements.txt), per-residue disorder
  propensity in [0, 1], computed on the reviewed human UniProt sequences (the
  same sequences used for residue validation).
- **State rule:** `is_disordered = score >= threshold`. **Primary threshold
  0.5; sensitivity thresholds 0.4 and 0.6** (exploratory, as §2/§11).
  Scores are computed once; thresholds are applied to the cached scores.
- **Everything else is unchanged** from §6/§18/§19: same mutation sets, mapping
  validation, cohorts, exclusions, NB2-t with log(E) offset, predictors
  (SBS4_std + TMB_std + wes), negative controls (SBS1 predictor; Silent
  outcome), equivalence bounds (0.833, 1.2), interpretation rules 18.3,
  control-failure rule, one-run discipline and output format (git commit +
  SHA-256 of all inputs).
- The static re-analysis uses the unchanged analysis code path
  (annotate_mutations_darkness -> build_background -> run_confirmatory /
  run_luad) with the disorder track in place of the pLDDT track. New runs are
  written to **new dated output files**; the §18/§19 run files are never
  replaced.

## 21. Mutation-induced disorder transitions (fixed 2026-11-24, before the runs)

**Question (external review).** Do mutations in tumours with a stronger smoke
signature differ in how they change the local disorder propensity of the
protein? Each mutation is classified into one of four transition classes by
the disorder state of the mutated position before and after the substitution:
**O-O, O-D, D-O, D-D** (O = ordered, D = disordered).

### 21.1 Per-mutation quantities
- Population: **validated missense mutations** (mapping validation exactly as
  §3; `uniprot_pos = protein_position + offset`). The alternative amino acid
  is parsed from `HGVSp_Short` (`^p?\.?([A-Z*])(\d+)([A-Z*])$`); rows without
  a parsable amino-acid change (nonsense, splice, ambiguous) are excluded and
  counted per cohort.
- `wt_score` = metapredict score of the WT UniProt sequence at `uniprot_pos`;
  `mut_score` = score at the same position after applying the substitution to
  the full sequence. `delta_pos = mut_score - wt_score`;
  `delta_window10` = mean score difference over the clipped window
  [pos-10, pos+10] (window length recorded per mutation).
- **Transition class** = `state_wt -> state_mut` at `uniprot_pos`, states from
  the §20 threshold rule (primary 0.5; sensitivity 0.4/0.6 re-derived from the
  same cached scores).

### 21.2 Silent negative control (fixed)
Silent mutations do not change the protein sequence: `delta_pos = 0` and
`delta_window10 = 0` exactly, and the class is O-O or D-D by construction.
The run asserts these invariants; any violation invalidates the run (pipeline
bug). The Silent class model uses the static disorder background
(§20 track): E[O-O] = 1 - f_dark(c), E[D-D] = f_dark(c).

### 21.3 O/E per transition class
- Per tumour t and class k: O_tk = number of its mutations of class k;
  E_tk = sum over the same mutations of p_k(context_96), where p_k(c) is the
  per-channel probability of class k. Cells with E_tk = 0 are dropped and
  counted.
- **Transition background (fixed).** p_k(c) is computed over the representative
  GENCODE v19 transcripts (longest CDS per gene, as in build_background) that
  translate to a protein **identical to the gene's UniProt sequence** (exact
  coordinates; genes counted). At every such coding site each of the 3 SNV
  alleles is translated; synonymous and stop outcomes are excluded; each
  remaining opportunity yields (context_96, alt residue) and, after mutant
  scoring, a transition class. Per channel, class fractions are estimated from
  **K = 5000 sampled opportunities** (reservoir sampling, seed 20260115);
  per-channel sampling SEs are reported. Channel coverage of the background
  is reported as in 18.6.

### 21.4 Models and rules
- Per class k (O-O, O-D, D-O, D-D) and the combined "disorder-changing" class
  (O-D + D-O, E = E_OD + E_DO): NB2-t `O ~ offset(log E) + SBS4_std + TMB_std
  (+ wes in SCLC)` with the §6 boundary rule, t-based CI, equivalence bounds
  (0.833, 1.2) and the 18.3 interpretation and control-failure rules.
- Controls (same code path): SBS1 as predictor (same outcomes); Silent
  mutations as outcome (21.2).
- Background-resampling test: per tumour and channel, multinomial draws over
  the 4 class probabilities; 10,000 runs; p reported as 1/(runs+1) bounds.
- **Continuous deltas are descriptive**: per-tumour mean `delta_pos` and
  mean `delta_window10` regressed on SBS4_std + TMB_std (+ wes) with HC1
  t-based CIs; no equivalence claim is made for score units.
- Sensitivity: thresholds 0.4/0.6 (21.1); mapping-discard covariate as 18.5.1;
  for LUAD: pack-years and never-vs-current secondaries as 19.5.

### 21.5 Cross-predictor robustness (fixed)
A second specialized predictor (IUPred2A standalone if available) scores a
seeded stratified subset of ~3000 mutations (strata: transition class x
cohort x |wt_score - 0.5| band). Reported in a separate output file:
transition-class concordance and Spearman correlation of deltas.
Pre-fixed: this check is descriptive and **cannot alter the primary results**.

### 21.6 Run discipline
One run per cohort and per analysis, one dated output file each, with git
commit and SHA-256 of all inputs. Rehearsals use the shuffled-predictor flag
and are labelled. Any post-unblinding change is a dated deviation (18.4); the
first result of each run is preserved.

---

## 22. Recovery and rebuild log (2026-11-24)

The working environment and its `data/` tree were restored from public sources
(Zenodo record 24288109, cBioPortal datahub/API, UCSC, GENCODE, AlphaFold DB,
UniProt, COSMIC). Record of what was and was not reproduced byte-for-byte:

- **Byte-identical (SHA-256 verified against the recorded run inputs):**
  `luad_smoking_data_clinical_patient.txt`, `tumor_signature_activities_refit_dedup.csv`,
  `luad_signature_activities.csv` — the signature chain (raw mutations ->
  context -> counts -> refit with seed 42) is exactly reproducible.
- **Not byte-identical (run inputs rebuilt with the public code):** the
  `*_mapping.csv`, `annotated_*.csv` and `background_pLDDT50.csv` files differ
  from the recorded hashes. Functionally the rebuild is close: mapping row
  counts match exactly (34583/8078/5190/66770), George validated rows 29927
  vs 29936 disclosed (0.03%), SCLC/LUAD annotated mutation counts within
  1.6%/4.7% of the run of record. Likely causes: UniProt sequence/mapping
  updates between 2026-10-28 and the rebuild, and code differences from the
  lost development history (commits 45c7c25 / 7f2605c7 are not in the public
  repository). The runs of record (sections 18-19) stand unchanged; nothing in
  this rebuild replaces them.
- **Bug found and fixed during recovery (2026-11-24, commit abfd734):**
  `map_mutations_to_residues.py` parsed `HGVSp_Short` without the `p.`
  prefix (`p.R175H`) used by cBioPortal MAFs, so any MAF lacking a `Codons`
  column (e.g. the Firehose `luad_tcga` MAF) validated zero rows. The fix
  accepts both notations (`^(?:p\.)?([A-Z*])(\d+)`). This is a reproduction
  bug fix; it does not alter any recorded result.
- The section 20-21 analyses run on the rebuilt pipeline (internally
  consistent inputs, same code path as sections 18-19).
