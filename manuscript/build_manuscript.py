"""Build the preprint (English) as a Word document for bioRxiv submission."""
from pathlib import Path
from docx import Document
from docx.shared import Pt, Cm

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"

TITLE = ("Tobacco mutational signature SBS4 and mutation placement in predicted "
         "disordered protein regions: a context-adjusted equivalence analysis "
         "in 373 lung tumours")
AUTHOR = "Arsen Turagabulov"
AFFIL = ('10th-grade student, School No. 125 "Kemel Urpaq", Taraz, Kazakhstan')
CORRESP = "Correspondence: [email to be inserted]"


def setup(doc):
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(11)
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Cm(2)
        s.left_margin = s.right_margin = Cm(2.5)


def p(doc, text, bold=False, size=None, align=None):
    par = doc.add_paragraph()
    run = par.add_run(text)
    run.bold = bold
    if size:
        run.font.size = Pt(size)
    if align:
        par.alignment = align
    return par


def figure(doc, name, caption, width=15):
    doc.add_picture(str(FIG / name), width=Cm(width))
    doc.paragraphs[-1].alignment = 1
    cap = doc.add_paragraph()
    run = cap.add_run(caption)
    run.font.size = Pt(9)
    run.italic = True


def table(doc, header, rows, caption):
    t = doc.add_table(rows=1 + len(rows), cols=len(header))
    t.style = "Table Grid"
    for j, h in enumerate(header):
        cell = t.rows[0].cells[j]
        cell.text = h
        for r in cell.paragraphs[0].runs:
            r.bold = True
            r.font.size = Pt(9)
    for i, row in enumerate(rows, start=1):
        for j, v in enumerate(row):
            t.rows[i].cells[j].text = str(v)
            for r in t.rows[i].cells[j].paragraphs[0].runs:
                r.font.size = Pt(9)
    cap = doc.add_paragraph()
    run = cap.add_run(caption)
    run.font.size = Pt(9)
    run.italic = True


def main():
    doc = Document()
    setup(doc)

    p(doc, TITLE, bold=True, size=14, align=1)
    p(doc, AUTHOR, size=12, align=1)
    p(doc, AFFIL, size=10, align=1)
    p(doc, CORRESP, size=10, align=1)

    p(doc, "Abstract", bold=True, size=12)
    p(doc, (
        "Intrinsically disordered protein regions lack stable structure and are "
        "often functional, but whether tobacco smoke mutagenesis places its "
        "mutations preferentially in such regions is unknown. I analysed 373 "
        "lung tumours: 143 small-cell lung carcinomas (SCLC) and 230 lung "
        "adenocarcinomas (LUAD) from TCGA, including 32 never-smokers, which "
        "extends the smoking range from zero to heavy exposure. Mutations were "
        "mapped to UniProt residues and labelled as disordered when AlphaFold "
        "pLDDT < 50. For each tumour I computed the observed number of "
        "mutations in disordered regions (O) and the number expected (E) under "
        "the tumour's own 96-context mutation spectrum against a genome-wide "
        "coding-space background; association with the smoke signature SBS4 "
        "was tested by Poisson/negative-binomial regression with log(E) as "
        "offset, with equivalence bounds 0.833-1.2 per SD fixed before the "
        "runs, negative controls (SBS1 as predictor, silent mutations as "
        "outcome) and a parametric background-resampling test. No association "
        "was found: SCLC O/E = 1.008 per SD (95% CI 0.985-1.033); across the "
        "full smoking range in LUAD O/E = 0.9998 (0.975-1.025), with pack-years "
        "(0.9991) and the never-vs-current contrast (1.017) agreeing; raw "
        "correlations were also null. Effects of 1.2 per SD and larger are "
        "excluded; confidence intervals constrain effects to about +/-3% per SD "
        "on the measured-predictor scale. The strength of the tobacco "
        "mutational signature does not track mutation placement in predicted "
        "disordered protein regions."
    ))
    p(doc, "Keywords: intrinsically disordered proteins; mutational signatures; "
           "SBS4; tobacco smoking; AlphaFold; lung cancer; observed/expected "
           "analysis", size=10)

    doc.add_heading("1. Introduction", level=1)
    p(doc, (
        "A large fraction of the human proteome contains intrinsically "
        "disordered protein regions (IDRs): segments that lack a stable "
        "three-dimensional fold and often act as flexible interaction hubs "
        "(1). Because IDRs tolerate sequence change differently from folded "
        "domains, the question of where cancer mutations land relative to "
        "disorder has attracted interest: a pan-cancer analysis mapped the "
        "mutational landscape of disordered hotspots (2), and driver mutations "
        "were shown to be depleted from, but not absent of, disordered regions "
        "(3)."))
    p(doc, (
        "Tobacco smoke leaves a quantitative mutational footprint in the tumour "
        "genome: single-base substitution signature SBS4 (4, 5). Its strength "
        "can be estimated per tumour from the mutation spectrum itself, which "
        "gives a continuous exposure measure rather than a smoker/non-smoker "
        "label. Whether smoke-associated mutations are placed preferentially in "
        "predicted disordered regions is, to my knowledge, untested - and there "
        "is no direct mechanism linking a DNA damage process to a property of "
        "the translated protein."))
    p(doc, (
        "Here I test the working hypothesis (H1) that, after adjusting for the "
        "tumour's own 96-context mutation spectrum, the share of mutations "
        "falling in predicted disordered regions increases with SBS4 strength; "
        "the null hypothesis (H0) is that no systematic shift remains after "
        "such adjustment. The study is framed as a coincidence analysis, not a "
        "causal claim. Two cohorts cover complementary ranges: 143 SCLC tumours "
        "(essentially all smokers) and 230 LUAD tumours spanning never-smokers "
        "to heavy smokers. Disordered positions are defined operationally as "
        "residues with AlphaFold pLDDT < 50 - a proxy for disorder rather than "
        "a strict biophysical IDR definition, a limitation kept explicit "
        "throughout."))

    doc.add_heading("2. Materials and Methods", level=1)
    doc.add_heading("2.1 Data and cohorts", level=2)
    p(doc, (
        "SCLC: somatic coding single-nucleotide variants (SNVs) from George et "
        "al. 2015 (WGS; 107 tumours after deduplication), Rudin et al. 2012 "
        "(30 tumours; 12 cell lines excluded using sample annotations), and "
        "Peifer et al. 2012 (6 tumours after deduplication) - 143 primary "
        "tumours in total (6, 7). LUAD: TCGA Firehose Legacy somatic mutations "
        "and clinical data for 230 primary tumours, including 32 lifelong "
        "never-smokers, with pack-years available for 158 (8). Structural "
        "labels: AlphaFold Protein Structure Database per-residue pLDDT (9, "
        "10). Annotations: GENCODE v19 CDS, UniProt/Swiss-Prot, hg19/GRCh37 "
        "(11, 12). Mutational signature definitions: COSMIC v3.3.1 (13)."))
    doc.add_heading("2.2 Mutation processing and residue validation", level=2)
    p(doc, (
        "Coding SNVs in primary tumours were retained; duplicate samples across "
        "studies (overlap coefficient >= 0.5) were removed keeping the "
        "higher-coverage sample. Each mutation was translated through its "
        "GENCODE transcript and the resulting amino acid compared with the "
        "UniProt sequence; only mutations whose residue matched were analysed, "
        "and discard fractions were recorded per cohort (SCLC 13-17%, LUAD "
        "21%). Reference alleles were checked against hg19."))
    doc.add_heading("2.3 Disorder annotation", level=2)
    p(doc, (
        "A residue was labelled disordered when AlphaFold pLDDT < 50. Mutations "
        "without pLDDT coverage were excluded from both observed and expected "
        "counts (coverage 92-96% of validated mutations). The threshold and all "
        "analysis rules were fixed before the runs (Supplementary: analysis "
        "plan)."))
    doc.add_heading("2.4 Observed and expected counts", level=2)
    p(doc, (
        "For each tumour, O is the number of validated mutations in disordered "
        "residues and E is the expected number under the tumour's own 96-context "
        "spectrum: the sum over its mutations of the background probability "
        "that a site of that substitution channel is disordered. The background "
        "was computed over the whole protein-coding space (all GENCODE CDS "
        "sites with pLDDT coverage); its per-channel coverage was verified "
        "channel-neutral (ratio 0.861 +/- 0.007 across 96 channels), and O and "
        "E always use the same mutation set."))
    doc.add_heading("2.5 Statistical model", level=2)
    p(doc, (
        "Primary model: Poisson regression of O with log(E) offset, covariates "
        "SBS4 fraction (z-scored), tumour mutation burden, and - in the "
        "combined SCLC analysis - data source; the reported measure is the "
        "observed/expected ratio per 1 SD of SBS4. Overdispersion is handled by "
        "negative-binomial regression with t-based confidence intervals "
        "(NB2-t), with a pre-specified boundary rule: when the dispersion "
        "estimate is at zero, the Poisson-t limit is used. Calibration by "
        "simulation (500 runs per scenario, gamma, lognormal and NB1-type "
        "overdispersion) kept type-I error within 3.5-6.5% and coverage within "
        "93-97%. Equivalence bounds (0.833, 1.2) per SD were fixed before the "
        "runs as the threshold of biological interest (the rationale was "
        "written down afterwards; the bound is not a detection limit - the "
        "realised design constrains effects to a few percent per SD). Because O "
        "is a bounded count (binomial-like, variance E(1-p)), the Poisson "
        "variance assumption is slightly conservative here, widening rather "
        "than narrowing intervals. Raw Pearson and Spearman correlations were "
        "computed as descriptive analyses without adjustment."))
    doc.add_heading("2.6 Controls and null test", level=2)
    p(doc, (
        "Two negative controls used the identical code path: SBS1 (age-related "
        "signature) as predictor, and silent mutations as outcome. A failure of "
        "any control would mark the result as questionable. A parametric "
        "background-resampling test (10,000 runs; binomial draws from the "
        "background per context; minimal attainable p = 1e-4, reported as a "
        "bound) assessed the global O/E level."))
    doc.add_heading("2.7 Analysis transparency", level=2)
    p(doc, (
        "Analysis rules, equivalence bounds and outcome interpretation were "
        "documented before the runs; the plan and a dated deviations log are "
        "available in the code repository. Deviations include: real "
        "coefficients seen during a rehearsal before the SCLC run (a flag was "
        "not passed to the model calls), and a change of the LUAD mutation "
        "source from PanCanAtlas to Firehose Legacy after the first attempt "
        "crashed (the PanCanAtlas profile contains no synonymous mutations, "
        "making the silent control impossible). Each analysis was executed "
        "exactly once to completion; results files are preserved unchanged. "
        "Global O/E showed a consistent depletion in disordered regions "
        "(0.89-0.94) of unestablished origin (see Limitations), so conclusions "
        "rest on within-cohort associations rather than absolute O/E levels."))
    doc.add_heading("2.8 Software and reproducibility", level=2)
    p(doc, (
        "All processing and analysis code, tests (including validation at known "
        "hotspots TP53 R175H and KRAS G12C), input checksums of the analysis "
        "runs, and the analysis plan are available at "
        "https://github.com/ArsikProgrammer/Darkcancer (Python: pandas, numpy, "
        "scipy, statsmodels)."))

    doc.add_heading("3. Results", level=1)
    doc.add_heading("3.1 Cohorts and quality control", level=2)
    p(doc, (
        "After duplicate and cell-line removal the SCLC cohort comprised 143 "
        "primary tumours; the LUAD cohort comprised 230 primary tumours, 32 of "
        "them from lifelong never-smokers. Residue validation retained 79-87% "
        "of coding SNVs depending on cohort (Table 1), and 92-96% of validated "
        "mutations had pLDDT coverage. The smoking scale behaved as expected "
        "before any outcome analysis: the median SBS4 fraction was 0.000 in "
        "never-smokers (n = 32), 0.29-0.39 in former smokers, and 0.401 in "
        "current smokers; Spearman rho was 0.314 with smoking status (p = "
        "2e-6) and 0.223 with pack-years (p = 0.005, n = 158) (Figure 2)."))
    figure(doc, "fig2_smoking_scale.png",
           "Figure 2. Smoke-signature validation (LUAD). SBS4 fraction by "
           "smoking status; never-smokers show near-zero SBS4, confirming the "
           "exposure scale before outcome analysis.")
    doc.add_heading("3.2 Primary analysis: SCLC among smokers", level=2)
    p(doc, (
        "Across 143 SCLC tumours the adjusted share of mutations in disordered "
        "regions did not track SBS4 strength: O/E = 1.008 per SD (95% CI "
        "0.985-1.033, NB2-t). Per the pre-specified rules this is an "
        "informative null: effects of 1.2 per SD and larger are excluded, and "
        "the interval constrains effects to roughly -1.5% to +3.3% per SD on "
        "the measured scale. Both negative controls passed (SBS1: 0.994, "
        "0.973-1.016; silent outcome: 1.032, 0.983-1.084)."))
    doc.add_heading("3.3 Dose across the full smoking range: LUAD", level=2)
    p(doc, (
        "In LUAD the smoking range extends to zero. The dose-response slope of "
        "the adjusted disordered-region share on SBS4 was again null: O/E = "
        "0.9998 per SD (95% CI 0.975-1.025). Independent smoking measures "
        "agreed: pack-years 0.9991 (0.978-1.021, n = 158), and the "
        "never-vs-current contrast 1.017 (0.966-1.071). Adding the "
        "mapping-discard fraction as covariate left the estimate unchanged "
        "(0.9991). Controls passed (SBS1: 0.979, 0.942-1.017; silent: 1.012, "
        "0.963-1.065). Figure 3 summarises all estimates."))
    figure(doc, "fig3_forest.png",
           "Figure 3. All model estimates (O/E ratio per 1 SD, 95% CI, log "
           "scale) with the pre-specified equivalence bounds (0.833, 1.2).")
    doc.add_heading("3.4 Raw (unadjusted) associations", level=2)
    p(doc, (
        "Without the context adjustment the raw correlation of the disordered "
        "mutation share with SBS4 was also null: Pearson r = -0.12 (p = 0.16) "
        "and Spearman rho = -0.03 (p = 0.75) in SCLC; Pearson r = +0.08 (p = "
        "0.25) and Spearman rho = +0.09 (p = 0.20) in LUAD (Figure 4). These "
        "are descriptive statistics without adjustment for cohort or mutation "
        "burden; the Pearson/Spearman divergence in SCLC suggests outlier "
        "influence."))
    figure(doc, "fig4_raw_scatter.png",
           "Figure 4. Raw disordered-region mutation share versus SBS4 "
           "fraction (each point: one tumour).")
    doc.add_heading("3.5 Global level of O/E", level=2)
    p(doc, (
        "Globally, fewer mutations landed in disordered regions than expected "
        "(SCLC O/E = 0.94, LUAD 0.89). A formal significance level is not "
        "reported for this global level: the background-resampling test treats "
        "mutations as independent, whereas real mutations cluster within genes "
        "and hotspots, which inflates significance; the deviation is therefore "
        "reported descriptively. Its cause is not established; technical "
        "factors (capture coverage, expression, replication timing, selection, "
        "annotation filtering) are possible, and this observation is not "
        "interpreted biologically."))
    figure(doc, "fig1_workflow.png",
           "Figure 1. Study design.")

    doc.add_heading("4. Discussion", level=1)
    p(doc, (
        "In both cohorts - among smokers with varying signature strength and "
        "across the full range from never-smokers to heavy smokers - the "
        "placement of mutations in predicted disordered regions did not track "
        "the tobacco mutational signature. Effects of 1.2 per SD and larger are "
        "excluded by the pre-specified equivalence bounds, and the confidence "
        "intervals narrow plausible effects to a few percent per SD on the "
        "measured scale. The result matches prior expectations: there is no "
        "mechanism linking a DNA damage signature to a property of the "
        "translated protein, and the study now quantifies how large such a "
        "link could be if present."))
    p(doc, (
        "Limitations. (i) The SBS4 estimate is noisy (NNLS on a few hundred "
        "coding mutations), which attenuates observed effects toward zero; "
        "pack-years, an independent measure, agreed with the primary result. "
        "(ii) Signature and outcome are computed from the same mutations "
        "(not applicable to pack-years). (iii) Disordered positions are "
        "defined by AlphaFold pLDDT < 50, a proxy rather than experimental "
        "disorder maps. (iv) Absolute O/E levels are shifted by an "
        "unestablished background bias (global 0.94/0.89); two candidate "
        "artefacts were checked and do not explain it: recomputing the "
        "background on genes whose representative transcript translates "
        "identically to UniProt (72% of genes) moved global O/E only to "
        "0.95/0.89, and restricting both mutations and background to those "
        "genes left the depletion unchanged (0.94/0.90); requiring a "
        "three-residue validation window instead of a single residue retained "
        "98.5-99.5% of validated mutations. Conclusions therefore rest on "
        "within-cohort "
        "associations. (v) Exact exome capture territories are unknown; the "
        "background uses the protein-coding space. (vi) All SCLC tumours are "
        "from smokers; the smoking-vs-no-smoking contrast lives in LUAD. "
        "Larger cohorts, experimental disorder maps and independent signature "
        "estimates would tighten the bound further. Comparisons with the "
        "reported enrichment of mutations in disordered regions in pan-cancer "
        "settings (2) should be made cautiously: the setups differ."))
    p(doc, (
        "Conclusions. Across 373 lung tumours and the full range of tobacco "
        "exposure, the strength of the smoke mutational signature does not "
        "change where mutations land relative to predicted disordered protein "
        "regions. Effects of biological interest (>= 1.2 per SD) are excluded; "
        "smaller effects, if any, lie below a few percent per SD on the "
        "measured scale."))

    doc.add_heading("Data and code availability", level=1)
    p(doc, (
        "All source data are public: cBioPortal datahub (SCLC studies "
        "sclc_ucologne_2015, sclc_jhu, sclc_clcgp), TCGA/GDC (LUAD), AlphaFold "
        "Protein Structure Database, GENCODE v19, UniProt/Swiss-Prot, COSMIC "
        "v3.3.1. Analysis code, tests, run outputs with input checksums and the "
        "analysis plan: https://github.com/ArsikProgrammer/Darkcancer"))

    doc.add_heading("References", level=1)
    refs = [
        "Wright PE, Dyson HJ. Intrinsically disordered proteins in cellular "
        "signalling and regulation. Nat Rev Mol Cell Biol. 2015;16:18-29.",
        "Zou H, Pan Y, Chen C, et al. Pan-cancer assessment of mutational "
        "landscape in intrinsically disordered hotspots reveals potential "
        "driver genes. Nucleic Acids Res. 2022;50(9):e49.",
        "Meszaros B, Gadi I, Dosztanyi Z, et al. The impact of structural "
        "disorder on protein function and disease. Biomolecules. "
        "2021;11(3):381. doi:10.3390/biom11030381",
        "Alexandrov LB, Ju YS, Haase K, et al. Mutational signatures "
        "associated with tobacco smoking in human cancer. Science. "
        "2016;354(6312):618-622.",
        "Alexandrov LB, Nik-Zainal S, Wedge DC, et al. Signatures of "
        "mutational processes in human cancer. Nature. 2013;500:415-421.",
        "George J, Lim JS, Jang SJ, et al. Comprehensive genomic profiles of "
        "small cell lung cancer. Nature. 2015;524:47-53. "
        "doi:10.1038/nature15338",
        "Rudin CM, Durinck S, Stawiski EW, et al. Integrative genome analyses "
        "identify key somatic driver mutations in small cell lung cancer. "
        "Nat Genet. 2012;44(10):1111-1116. PMID:22941188. (Peifer M, et al. "
        "Nat Genet. 2012;44(10):1104-1110. PMID:22941189.)",
        "The Cancer Genome Atlas Research Network. Comprehensive molecular "
        "profiling of lung adenocarcinoma. Nature. 2014;511:543-550. "
        "doi:10.1038/nature13385",
        "Jumper J, Evans R, Pritzel A, et al. Highly accurate protein "
        "structure prediction with AlphaFold. Nature. 2021;596:583-589. "
        "doi:10.1038/s41586-021-03819-2",
        "Varadi M, Anyango S, Deshpande M, et al. AlphaFold Protein Structure "
        "Database. Nucleic Acids Res. 2022;50:D439-D444. "
        "doi:10.1093/nar/gkab1061",
        "Harrow J, Frankish A, Gonzalez JM, et al. GENCODE: the reference "
        "human genome annotation. Genome Res. 2012;22:1760-1774. "
        "doi:10.1101/gr.135350.111",
        "UniProt Consortium. UniProt: the Universal Protein Knowledgebase in "
        "2023. Nucleic Acids Res. 2023;51:D523-D531. doi:10.1093/nar/gkac1000",
        "Tate JG, Bamford S, Jubb HC, et al. COSMIC: the Catalogue Of "
        "Somatic Mutations In Cancer. Nucleic Acids Res. 2019;47:D941-D947.",
        "Gao J, Aksoy BA, Dogrusoz U, et al. Integrative analysis of complex "
        "cancer genomics and clinical profiles using the cBioPortal. Sci "
        "Signal. 2013;6(269):pl1.",
        "Harris CR, Millman KJ, van der Walt SJ, et al. Array programming "
        "with NumPy. Nature. 2020;585:357-362. doi:10.1038/s41586-020-2649-2",
        "Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental "
        "algorithms for scientific computing in Python. Nat Methods. "
        "2020;17:261-272. doi:10.1038/s41592-019-0686-2",
        "Seabold S, Perktold J. statsmodels: Econometric and statistical "
        "modeling with Python. Proc 9th Python in Science Conf. 2010.",
        "McKinney W. Data structures for statistical computing in Python. "
        "Proc 9th Python in Science Conf. 2010:56-61.",
    ]
    for i, r in enumerate(refs, 1):
        par = doc.add_paragraph()
        run = par.add_run(f"{i}. {r}")
        run.font.size = Pt(9)

    doc.add_heading("Tables", level=1)
    table(doc,
          ["Cohort", "Samples", "Coding SNVs", "Validated (%)", "pLDDT coverage (%)"],
          [["SCLC: George (WGS)", "107", "34,583", "29,936 (86.6)", "92.6"],
           ["SCLC: JHU (WES)", "30", "5,190", "4,403 (84.8)", "96.0"],
           ["SCLC: CLCGP (WES)", "6", "2,494", "2,070 (83.0)", "94.3"],
           ["LUAD: TCGA (WES)", "230 (32 never)", "66,770", "52,914 (79.2)", "95.1"]],
          "Table 1. Cohorts and quality control. All counts are restricted to "
          "the analysis set: 143 SCLC tumours after duplicate-sample and "
          "cell-line removal, and 230 LUAD tumours (per-cohort figures are "
          "post-filtering and sum to the analysis totals). Validated (%) is "
          "the share of coding SNVs with confirmed residue mapping; pLDDT "
          "coverage (%) is the share of validated mutations with a disorder "
          "label. Mutations entering the primary models require pLDDT "
          "coverage: SCLC 33,887 (143 tumours; equals 27,709 + 4,226 + 1,952 "
          "across cohorts) and LUAD 50,339 (230 tumours). Silent-mutation "
          "control: 7,171 (113 tumours) and 12,081 (230 tumours).")
    table(doc,
          ["Analysis", "O/E per 1 SD", "95% CI", "Interpretation"],
          [["SCLC primary (SBS4)", "1.008", "0.985-1.033", "informative null"],
           ["LUAD primary (SBS4)", "0.9998", "0.975-1.025", "informative null"],
           ["LUAD pack-years", "0.9991", "0.978-1.021", "informative null"],
           ["LUAD never vs current", "1.017", "0.966-1.071", "informative null"],
           ["LUAD discard covariate", "0.9991", "0.974-1.024", "informative null"],
           ["Control SBS1 (SCLC / LUAD)", "0.994 / 0.979", "0.973-1.016 / 0.942-1.017", "passed"],
           ["Control silent (SCLC / LUAD)", "1.032 / 1.012", "0.983-1.084 / 0.963-1.065", "passed"]],
          "Table 2. Model estimates. Equivalence bounds (0.833, 1.2) were fixed "
          "before the runs; rules for informative null / inconclusive / effect "
          "are given in the analysis plan.")

    out = HERE / "preprint.docx"
    doc.save(out)
    print("saved:", out)


if __name__ == "__main__":
    main()
