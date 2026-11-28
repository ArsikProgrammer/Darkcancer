# Dark Proteome and the Mutational Footprint of Tobacco Smoke

[![tests](https://github.com/ArsikProgrammer/Darkcancer/actions/workflows/tests.yml/badge.svg)](https://github.com/ArsikProgrammer/Darkcancer/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24288109.svg)](https://doi.org/10.5281/zenodo.24288109)

A research project on open data: I test whether the mutational signature of
tobacco smoke in the tumour genome is linked to where mutations land — in
"soft" (intrinsically disordered, without a stable fold) protein regions — or
not.

Preprint and archived record: [10.5281/zenodo.24288109](https://doi.org/10.5281/zenodo.24288109).

## Idea and hypothesis

"Soft" protein regions are segments without a stable three-dimensional shape.
In this project, "dark proteome" is the working name for such regions,
operationally defined as residues with low AlphaFold prediction confidence
(pLDDT < 50). This is a **proxy** for disorder, not a strict definition of
intrinsically disordered regions (IDRs) in the biophysical sense — a limitation
emphasised again below.

The working hypothesis (H1): after adjusting for the tumour's own 96-context
mutation spectrum — that is, its DNA damage "footprint" — the share of
mutations in soft regions increases with the strength of the smoke signature
(SBS4). The null hypothesis: no systematic shift remains after such an
adjustment. Importantly, what is tested is exactly the **residual association
after the context adjustment**; this is a narrow but correctly specified
hypothesis. No causality is claimed: there is no direct mechanism linking
"smoke" to "protein", and the study addresses coincidence only.

## Data

| What | Source | Scope |
|---|---|---|
| Small-cell lung cancer (SCLC) | George et al. 2015 (WGS), Rudin et al. 2012 (JHU), Peifer et al. 2012 (CLCGP) via cBioPortal datahub | 143 primary tumours after removing duplicates and cell lines |
| Adenocarcinoma (LUAD) | TCGA, Firehose Legacy (mutations + smoking clinical data) | 230 primary tumours, including 32 never-smokers |
| Protein structures | AlphaFold Protein Structure Database (per-residue pLDDT) | whole proteome |
| Annotations | GENCODE v19, UniProt/Swiss-Prot, COSMIC v3.3.1, hg19/GRCh37 | — |

How to download the data: see `download_tar.py` (AlphaFold archive),
`download_luad.py` (TCGA mutations and clinical data), scripts `src/01–05`
(SCLC MAFs from cBioPortal datahub), `src/download_uniprot_fasta.py` (UniProt).

## Methods

1. **Mutations → protein residue.** I translate each mutation through its
   transcript (GENCODE v19, CDS, hg19) and compare the resulting amino acid
   with UniProt. Only mutations with a confirmed residue match enter the
   analysis; discarded fractions are reported per cohort (13–17% in SCLC, 21%
   in LUAD).
2. **Annotation of soft regions.** pLDDT < 50 (AlphaFold). The threshold and
   the rules were fixed before the analysis (see `ANALYSIS_PLAN.md`).
3. **Observed and expected (O/E).** For each tumour: how many mutations landed
   in soft regions (O) and how many were expected (E) under a random placement
   that preserves the 96-context spectrum of that tumour. The background of
   "how many sites are affected by each substitution type" was computed over
   the entire protein-coding space.
4. **Model.** Poisson regression (robust variant — negative binomial) with
   log(E) as offset; the main measure is the O/E ratio per 1 standard
   deviation of the smoke signature. The threshold of an "interesting" effect
   (0.833–1.2) was fixed before the analysis (see `ANALYSIS_PLAN.md`, §7) as
   the threshold of biological interest: an association weaker than 20% per SD
   in the share of mutations in protein regions is hard to interpret
   biologically. This is **not** a detection limit: from the width of the
   confidence intervals, the design constrains effects to roughly ±3% per SD
   (with a caveat about predictor noise, see below). Robustness to
   overdispersion was checked by simulations (500 runs per scenario).
5. **Controls.** The "age-related" signature SBS1 as a negative predictor;
   synonymous mutations as a negative outcome. Failure of any control means
   the result is "under question".
6. **Background test (parametric).** Mutations are "reshuffled" onto random
   coding positions while preserving substitution types; the implementation is
   a binomial draw from the background for each context (this is a parametric
   test of the global O/E, not a classical label permutation). With 10,000
   runs the minimum attainable p ≈ 1e-4; values at this bound should be read
   as "p ≤ 1e-4".
7. **Analysis discipline.** The plan with rules and bounds was fixed before
   the runs (`ANALYSIS_PLAN.md`); the main run for each cohort was executed
   exactly once; all deviations are recorded with dates in the same document.
8. **Specialized disorder predictors (external review, plan §20).** The
   pLDDT proxy was replaced by metapredict (a sequence-based disorder
   predictor) and the full analysis was re-run on the new track
   (`results/confirmatory_disorder_run_20261124.txt`,
   `results/luad_disorder_run_20261124.txt`); thresholds 0.4/0.6 as
   sensitivity. IUPred2A cross-check on a stratified subset
   (`results/iupred_crosscheck_20261124.txt`).
9. **Mutation-induced disorder transitions (external review, plan §21).**
   Each validated missense mutation is scored on the WT and mutant protein
   sequences and classified by the disorder state of the mutated position
   before -> after: **O-O, O-D, D-O, D-D**. Per-class O/E against a
   context-preserving transition background (proteins identical to their
   UniProt sequence), NB2-t models with the same controls and equivalence
   bounds; silent mutations as a built-in negative control (delta = 0 exactly).
   Runs: `results/disorder_transitions_run_20261124.txt`,
   `results/disorder_transitions_luad_run_20261124.txt` (+ t0.4/t0.6
   sensitivity).

## Results

**Small-cell lung cancer (143 smokers).** After the context adjustment, the
smoke signature was not associated with the share of mutations in soft
regions: O/E = 1.008 per SD (95% interval 0.985–1.033) — an informative null:
effects of about 1.2 and larger are excluded.

**Adenocarcinoma (230 tumours, a range from zero in never-smokers to high
values).** The share of mutations in soft regions did not track the smoke
signature anywhere along the range: O/E = 0.9998 per SD (0.975–1.025).
Independent smoking measures agree: pack-years 0.9991, never-versus-current
smokers 1.017.

**Raw (unadjusted) association** (script `src/raw_association.py`, output in
`results/raw_association.txt`) is also absent in both cohorts: the correlation
of the soft-region mutation share with SBS4 — Pearson r = −0.12 (p = 0.16) and
Spearman ρ = −0.03 (p = 0.75) in SCLC; Pearson r = +0.08 (p = 0.25) and
Spearman ρ = +0.09 (p = 0.20) in LUAD. Thus the conclusion holds even without
the context adjustment. This is descriptive statistics without adjustment for
cohort or mutation burden; the divergence between Pearson and Spearman in SCLC
indicates outlier influence.

**Side observation.** There are somewhat fewer mutations in soft regions than
expected (O/E ≈ 0.89–0.94). The cause is not established; technical causes are
possible (sequencing coverage, expression, replication timing, selection,
filtering specifics during annotation). This is not used as a biological
conclusion.

**Bottom line.** No residual association was found between the smoke signature
and the distribution of mutations in soft protein regions: neither among
smokers with differing signature strength, nor across the full range from
never-smokers to heavy smokers. Effects ≥ 1.2 per SD are excluded; moreover,
the confidence intervals narrow plausible effects to a few percent per SD — on
the scale of the **measured** signature (due to predictor noise the true
effect could be larger than the measured one, see limitations).

**Disorder-predictor re-analysis (external review).** With metapredict in
place of the AlphaFold proxy the conclusion is unchanged: O/E = 0.9958 per SD
(0.9751–1.0170) in SCLC and 0.9952 (0.9744–1.0164) in LUAD — informative
nulls; all controls pass.

**Disorder transitions (external review).** Of 29 582 SCLC and 38 837 LUAD
validated missense mutations, 1.1–3.0% change the local disorder state
(O-D or D-O). No class tracks the smoke signature: O-D shows a
non-significant negative tendency in both cohorts (SCLC 0.909, 0.795–1.038;
LUAD 0.945, 0.824–1.083 — inconclusive by the equivalence rules), the other
classes are informative nulls. The per-tumour mean disorder change is
descriptively negative with SBS4 in LUAD (b = −0.0028 per SD, CI excludes 0)
but tiny and not pre-registered as an inferential claim. The global
depletion observed with the pLDDT track decomposes into fewer D-D and more
O-O transitions than the context-preserving background expects (both
p ≤ 1e-4), with O-D also depleted.

## Limitations

- **Absolute O/E values cannot be trusted.** The expectation model shows a
  persistent global shift (O/E ≈ 0.89–0.94) — probably due to an imperfect
  background (representative transcript per gene, expression, replication,
  selection, mutation filtering). **Relative** within-cohort comparisons are
  reliable (hence the main conclusions are drawn from associations, not from
  absolute levels).
- The smoke signature is estimated by NNLS from several signatures on hundreds
  of coding mutations — an imprecise predictor; this attenuates the observed
  effect. Pack-years, an independent measure, give the same result.
- The smoke signature and the outcome are computed from the same mutations
  (not an issue for pack-years).
- "Softness" is an AlphaFold prediction (pLDDT < 50), not an experiment and
  not a strict definition of disorder. Per external review the analysis was
  repeated with metapredict (a specialized sequence-based predictor; plan
  §20–21) with the same result; the mutant-side cross-check with a second
  predictor (IUPred2A) requires its standalone distribution and is reported
  for WT states only (82% concordance at mutated positions, Spearman 0.77).
- Exact exome capture territories are unknown; the background is computed over
  the protein-coding space.
- Conclusions concern variation within a single cancer type; the causal effect
  of "smoking versus no smoking" is not established by this design.
- The analysis plan was fixed internally (dates in `ANALYSIS_PLAN.md`), without
  external registration (OSF etc.). Files in `results/` are kept as they are:
  they use the earlier terminology (the "PERMUTATION" label) and contain commit
  hashes of the local development history, which is not in this repository;
  the current terminology is in the code and in this README.

## Reproduction

Order of steps (after downloading the data):

```
python run_pipeline.py                 # ingest -> context -> signatures
python src/map_mutations_to_residues.py ... / src/map_jhu_coordinates.py --mutations-csv ...
python src/build_alphafold_darkness.py --all-residues --output data/processed/darkness_proteome_pLDDT50.csv
python src/build_background.py --darkness data/processed/darkness_proteome_pLDDT50.csv --output data/processed/background_pLDDT50.csv
python src/build_opportunity.py --output data/processed/opportunity_96.csv
python src/refit_signatures.py --opportunity data/processed/opportunity_96.csv ...
python src/annotate_mutations_darkness.py --mapping ... --darkness ...
python src/run_confirmatory.py -- ...   # SCLC (single run)
python src/run_luad.py -- ...           # LUAD (single run)
```

Disorder-predictor analyses (plan §20–21):

```
python src/build_metapredict_disorder.py --threshold 0.5 --output data/processed/disorder_proteome_metapredict_t0.5.csv
python src/build_background.py --darkness data/processed/disorder_proteome_metapredict_t0.5.csv --output data/processed/background_metapredict_t0.5.csv
python src/annotate_disorder_transitions.py --mutations ... --mapping ... --export-subs data/processed/_subs_<cohort>.csv
python src/score_disorder_effects.py --substitutions data/processed/_subs_all.csv --output data/processed/disorder_substitution_scores.csv
python src/build_transition_background.py --threshold 0.5 --output data/processed/transition_background_metapredict_t0.5.csv
python src/annotate_disorder_transitions.py --mutations ... --mapping ... --scores data/processed/disorder_substitution_scores.csv --threshold 0.5 --output data/processed/disorder_transitions_t0.5_<cohort>.csv
python src/run_disorder_transitions.py -- ...        # SCLC (single run)
python src/run_disorder_transitions_luad.py -- ...   # LUAD (single run)
python src/crosscheck_iupred.py --mutations ... --mode rest --n 300
```

Exact arguments are in the module docstrings; the final run commands are also
duplicated in `ANALYSIS_PLAN.md`. Tests (`tests/`): 16 tests (context,
statistics, O/E logic, position parsing, disorder-transition logic) run
everywhere, including CI; 2 more — mapping validation at the known hotspots
TP53 R175H and KRAS G12C — run locally with the `data/` directory only and are
automatically skipped on a clean clone (18 in total locally).

## Repository structure

```
src/       — the pipeline: ingest (01-10), context, mapping to residues,
             darkness annotation, O/E background, signatures, models, simulations
tests/     — context and mapping checks at known mutations
results/   — final run outputs with input checksums
ANALYSIS_PLAN.md — the analysis plan, fixed before the runs
download_*.py, run_pipeline.py, requirements.txt, LICENSE
```

## References

1. George J. et al. Comprehensive genomic profiles of small cell lung cancer.
   *Nature* 524, 47–53 (2015). doi:10.1038/nature15338.
2. Rudin C.M. et al. Integrative genome analyses identify key somatic driver
   mutations in small cell lung cancer. *Nat Genet* 44, 1111–1116 (2012).
   PMID: 22941188.
3. Peifer M. et al. Comprehensive genomic analysis identifies SOX2 as a
   frequently amplified oncogene in small cell lung cancer. *Nat Genet* 44,
   1104–1110 (2012). PMID: 22941189.
4. The Cancer Genome Atlas Research Network. Comprehensive molecular profiling
   of lung adenocarcinoma. *Nature* 511, 543–550 (2014).
   doi:10.1038/nature13385.
5. Alexandrov L.B. et al. Mutational signatures associated with tobacco smoking
   in human cancer. *Science* 354, 618–622 (2016).
6. Zou H. et al. Pan-cancer assessment of mutational landscape in intrinsically
   disordered hotspots reveals potential driver genes. *Nucleic Acids Res* 50(9),
   e49 (2022). doi:10.1093/nar/gkac160.
7. Mészáros B. et al. The impact of structural disorder on protein function and
   disease. *Biomolecules* 11(3), 381 (2021). doi:10.3390/biom11030381.
8. Alexandrov L.B. et al. The repertoire of mutational signatures in human
   cancer. *Nature* 578, 94–101 (2020). doi:10.1038/s41586-020-1943-3 (COSMIC).
9. Jumper J. et al. Highly accurate protein structure prediction with AlphaFold.
   *Nature* 596, 583–589 (2021). doi:10.1038/s41586-021-03819-2.
10. Varadi M. et al. AlphaFold Protein Structure Database. *Nucleic Acids Res*
    50, D439–D444 (2022). doi:10.1093/nar/gkab1061.
11. Harrow J. et al. GENCODE: the reference human genome annotation.
    *Genome Res* 22, 1760–1774 (2012). doi:10.1101/gr.135350.111.
12. UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023.
    *Nucleic Acids Res* 51, D523–D531 (2023). doi:10.1093/nar/gkac1000.

## License

MIT, see `LICENSE`.
