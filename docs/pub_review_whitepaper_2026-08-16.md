# Publication Review White Paper: 08-mmrm-linearity-robust
*Review date: 2026-08-16 10:10 PDT*

Reviewer-style gap analysis of the research report in this
compendium, applying the standards of a referee for a statistical
journal (Statistics in Medicine, Pharmaceutical Statistics,
Biometrics). Epistemic status is marked throughout as verified
(ran code), inspected (read source or rendered output), inferred,
or unverified.

## 1. Summary of the work under review

One manuscript exists: `analysis/report/report.Rmd` (rendered to
`analysis/report/report.pdf`, staged copy dated 2026-08-15 in
`analysis/report/share/`). The paper is a Monte Carlo simulation
study of the robustness of a random-slopes mixed model to
nonlinearity in the treatment effect trajectory in Alzheimer's
disease trials with CDR-SB as the endpoint. Data are generated for
200 participants per arm over 18 months (visits every 3 months),
with the treatment effect accruing as g(t) = (t/18)^kappa for
kappa in {0.3, 0.5, 0.7, 1, 1.5, 2, 3} and a fixed total effect
of -0.45 points at month 18. Two analysis models are compared: a
random-slopes linear mixed model (`nlme::lme`) and a
categorical-time MMRM with unstructured covariance
(`mmrm::mmrm`), each summarized by convergence, bias, empirical
and model-based SE, power, and coverage with Monte Carlo SEs
(1000 replicates per condition). The simulation engine is
`analysis/scripts/sim_mmrm_linearity.R`; a prior internal ADEMP
audit is at `docs/morris-audit-2026-04-17.md`. All results in
the manuscript are computed inline from the sourced script
(inspected; I did not rerun the simulation).

## 2. Major issues

1. **Results narrative contradicts the actual results.**
   Location: `report.Rmd`, "Power comparison" section and
   Discussion. The prose asserts, twice, that under linearity the
   random-slopes model "is expected to achieve higher power due to
   its more parsimonious mean structure" and in the Discussion
   that it shows "higher power than the categorical-time model."
   The rendered Table 1 (inspected in `report.tex`, lines
   505-513) shows the opposite: at kappa = 1 the slope model has
   power 0.542 versus 0.623 for the categorical MMRM, and the
   slope model's empirical SE (0.221) exceeds the categorical
   model's (0.191). The categorical model dominates or ties on
   power at every kappa. The results sections read as templated
   expectations written before the simulation was run and never
   updated. For a referee this is disqualifying as submitted: the
   central efficiency claim motivating the slope estimand is
   refuted by the paper's own table, and the paper does not
   notice. Remediation: rewrite Results and Discussion around the
   observed numbers, and explain mechanistically why the slope
   model is not more efficient here (see issue 2, which is likely
   the cause).

2. **The random-slopes estimator is questionably defined, and
   this choice, not linearity per se, plausibly drives the
   headline findings.** Location: `report.Rmd` "Analysis models";
   `sim_mmrm_linearity.R`, `fit_models()` (lines 91-136). The
   slope model is fit to post-baseline data only (months 3-18),
   includes a treatment main effect, yet the 18-month effect is
   extracted as the interaction coefficient times 1.5, discarding
   the fitted treatment main effect. The model-implied difference
   at month 18 is main effect plus 1.5 times interaction. Under
   nonlinear g(t) the linear projection splits the treatment
   signal between intercept and slope, so the reported "bias of
   the random-slopes model" conflates linearity failure with an
   incomplete contrast; and under linearity the redundant main
   effect (true value 0 at t = 0, but t = 0 is excluded from the
   fit) inflates the interaction variance, which explains the
   power deficit in issue 1. Standard practice for slope
   estimands either omits the treatment main effect (constraining
   equality at baseline, as in cLDA) or reports the full
   model-implied contrast. A referee will ask which estimator
   practitioners actually use and whether the bias attributed to
   "the linearity assumption" survives the alternative estimator.
   Remediation: justify the estimator against practice
   (e.g., Liu-Seifert et al. 2015), add the constrained variant
   (no treatment main effect, or cLDA including baseline), and
   report both. This is required for correctness of the paper's
   interpretation, not polish. Inspected; the mechanism claimed
   here is inferred and should be verified by simulation.

3. **No null scenario: type I error is never assessed.**
   Location: entire design (`sim_mmrm_linearity.R`,
   `run_simulation()`; report Methods). All conditions use
   delta = -0.45. Power comparisons between a misspecified and a
   correctly specified model are uninterpretable without
   demonstrating size control, and for the slope model under
   nonlinearity the null is exactly the case where regulators
   worry about inflated false positives (an early transient
   effect that vanishes by month 18 is a kappa-like violation
   with delta = 0 at t_max). Remediation: add delta = 0 scenarios
   across kappa, and ideally a "transient effect" trajectory that
   is nonzero at interim visits and zero at 18 months.
   Required for acceptance at any statistical journal.

4. **Internal contradictions between the Methods text and the
   code.** Location: `report.Rmd` "Analysis models" versus ADEMP
   section and `fit_models()`. (a) The Analysis models section
   states the categorical MMRM is fit "via `nlme::gls`"; the
   ADEMP bullet and the code use `mmrm::mmrm`. (b) The claim
   that the final-visit SE is "obtained via the delta method" is
   wrong: the code computes a linear contrast with Satterthwaite
   degrees of freedom (`mmrm::df_1d`); no delta method is
   involved. (c) The Methods say "Both models are fit to
   post-baseline observations only" inside the random-slopes
   paragraph, leaving the reader to infer the categorical model's
   handling of baseline; neither model adjusts for baseline
   CDR-SB, which departs from the near-universal practice of
   baseline-adjusted change-from-baseline MMRM
   (Mallinckrodt et al. 2008) and must be stated and defended.
   Remediation: make text match code exactly; add baseline
   adjustment or justify its absence. Inspected.

5. **The manuscript's own compliance section is stale and
   self-contradictory, and sits after the References heading.**
   Location: `report.Rmd` lines 671-693. The section "Morris et
   al. (2019) ADEMP Compliance" appears as a subsection of
   "# References", so it renders after the bibliography heading;
   an internal audit does not belong in a submitted manuscript at
   all. Worse, its "key gaps" (n_rep mismatch, seed set inside
   the worker, no MCSE reported) describe an earlier state of the
   code: the current script sets no seed internally, the Rmd pins
   `RNGkind("L'Ecuyer-CMRG")` and seeds once, and
   `summarize_results()` returns MCSE columns (inspected). The
   paper therefore asserts, as current, defects that no longer
   exist. Remediation: delete the section from the manuscript
   (keep the audit in `docs/`), or convert it to a brief,
   accurate reproducibility statement.

6. **Promised MCSEs and convergence rates never appear in any
   table.** Location: `report.Rmd` "Headline findings" ("Full
   Morris Table 6 MCSEs are reported in the tables below") and
   the `tab-summary` chunk. The displayed table selects only
   bias, SEs, power, and coverage; no `mcse_*` or `*_conv`
   columns are shown anywhere, although `summarize_results()`
   computes them. The text's promise is false as rendered
   (inspected in `report.tex`). Remediation: add MCSE columns
   (or parenthetical MCSEs) and a convergence column per Morris
   et al. (2019) Table 6, or a supplementary table.

7. **Single design point limits the contribution below
   publication threshold.** Location: Methods. One sample size
   (200/arm), one effect size, one covariance structure, no
   missing data, and n_sim = 1000 that the paper itself concedes
   "should be raised for publication." At 200/arm the powers sit
   at 0.54-0.66, whereas the trial being mimicked (CLARITY-AD)
   randomized roughly 900 per arm; the calibration claim is
   therefore only partial. The stated novelty over Chen et al.
   (2018) is systematic variation of nonlinearity, but with no
   factorial variation in n, delta, dropout, or covariance, the
   study is thinner than the work it claims to extend.
   Remediation: at minimum add a second sample size near the
   actual trial scale, the null scenarios of issue 3, and a
   simple MAR dropout mechanism; raise n_sim to meet the paper's
   own MCSE target (n_rep >= 2066 is already derived in the
   Methods). Required for acceptance.

8. **Citation integrity and coverage.** Location:
   `analysis/report/references.bib`; Introduction and Discussion.
   (a) Several entries have author lists that do not match the
   published papers and appear partly invented: `chen2018simulation`
   lists "Bhatt, Dhairya" among the authors of a 2018 ADAS-Cog
   simulation paper; `wang2022proportional` lists Logovinsky,
   Bhatt, Bhore, and Sperling for the proportional cLDA paper
   (the pCLDA literature is associated with Guoqiao Wang and
   colleagues); `donohue2011power` is keyed and captioned as a
   power paper but contains the 2014 JAMA Neurology PACC paper;
   `liang2000longitudinal` is keyed 2000 with year 1986 and
   contains the Liang and Zeger GEE paper, almost certainly a
   mis-capture of Liang and Zeger (2000, Sankhya), the standard
   cLDA reference; `novic2024nonlinear`'s author list is
   unverified. The bib file's own header documents one previous
   fabricated-entry correction, so every entry must be verified
   against the published record before submission (unverified;
   flagged from internal inconsistencies). (b) Missing
   literature: the cLDA reference itself (Liang and Zeger 2000),
   the estimand-comparison and disease-progression-model work
   (e.g., Guoqiao Wang et al. on pCLDA and on slope vs
   change-from-baseline; Raket's progression-model papers;
   Donohue et al. on natural cubic spline MMRM), and recent
   commentary on time-based versus visit-based estimands for the
   anti-amyloid trials. (c) Six bib entries are never cited in
   the text (fitzmaurice2011applied, donohue2011power, ard2011,
   liu2015slopes, budd2022aducanumab, liang2000longitudinal).
   Remediation: verify every entry, cite or remove the orphans,
   and position the paper against the cLDA and progression-model
   literature.

## 3. Minor issues

1. `report.Rmd` has no abstract, no keywords, and no data/code
   availability statement; the target-journal format requires all
   three.

2. Figure 1 assigns "dotted" to both kappa = 0.3 and kappa = 3
   (`ltype_map`, report.Rmd lines 372-380), making the two most
   extreme curves indistinguishable in a grayscale figure.
   Inspected.

3. Coverage is computed with the normal quantile 1.96
   (`summarize_results()`), while p-values come from t or
   Satterthwaite reference distributions; use each model's own
   critical value for internal consistency.

4. The bias mechanism text ("the linear model overestimates the
   treatment effect at intermediate time points") is ambiguous
   given the observed direction: for kappa < 1 the slope
   estimate is attenuated toward zero (bias +0.227 against a
   true -0.45), i.e., the 18-month effect is underestimated in
   magnitude.

5. The simulation runs inside a cached knitr chunk
   (`cache=TRUE`); roughly 7000 unstructured `mmrm` fits execute
   at render time on a cache miss. No runtime or hardware
   disclosure is given, and stale-cache risk is real. Move the
   simulation to a script that saves a versioned results object
   (the orphan `analysis/data/derived_data/sim_08.rds` suggests
   this was once intended) and have the Rmd read it.

6. `analysis/data/README.md` is untouched Palmer-penguins
   template boilerplate, unrelated to this project, and
   `sim_08.rds` is referenced nowhere in code or text
   (verified by grep). A referee or editor checking the
   compendium would notice immediately.

7. The test suite is a placeholder (`inst/tinytest/test_basic.R`
   contains only `expect_true(TRUE)`); no test exercises
   `simulate_one_trial()`, `fit_models()`, or
   `summarize_results()` (e.g., mean structure of the DGM,
   MCSE formulas). Inspected.

8. Terminology: the title and text call the random-slopes model
   an "MMRM." In the clinical trials literature MMRM denotes the
   categorical-time marginal model; calling both "MMRM" invites
   referee objections. "Random-slopes LMM versus MMRM" is the
   conventional contrast.

9. CDR-SB is a bounded, discrete score (0-18 in half-point
   steps); the Gaussian DGM can produce out-of-range values.
   The limitation paragraph mentions Gaussianity but not
   boundedness or discreteness.

10. `DESCRIPTION` places all dependencies in Suggests, `R/` and
    `vignettes/` are empty, and the compendium's package shell
    adds nothing; harmless for a paper, but the data/code
    availability statement should point to the script, not the
    package.

11. Line-numbering (`lineno`) plus `toc: true` is a working
    draft format; most target journals want an unnumbered TOC
    removed at submission. Trivial.

## 4. What remains to be done

Priority order for submission readiness.

(a) Required for correctness

1. Rewrite Results and Discussion to match the actual numbers;
   remove all "is expected to" language describing outcomes the
   table already reports (major issue 1).
2. Resolve the slope-estimator definition: add the constrained
   estimator and/or the full model-implied contrast, rerun, and
   re-interpret; state clearly which estimator practice uses
   (major issue 2).
3. Fix Methods/code contradictions: gls versus mmrm, delta
   method claim, baseline handling (major issue 4).
4. Remove or correct the stale ADEMP compliance section and its
   placement after References (major issue 5).
5. Verify every bibliography entry against the published record;
   fix or delete suspect author lists (major issue 8a).

(b) Required for acceptance

6. Add null (delta = 0) and transient-effect scenarios; report
   type I error (major issue 3).
7. Expand the design: second sample size near 900/arm, MAR
   dropout, n_sim >= 2100 per the paper's own MCSE derivation
   (major issue 7).
8. Report MCSEs and convergence rates in the tables as promised
   (major issue 6).
9. Add abstract, keywords, data/code availability statement
   (minor issue 1).
10. Position against cLDA, pCLDA, and progression-model
    literature; cite or drop orphan references (major issue 8b,
    8c).

(c) Desirable polish

11. Distinct line types in Figure 1; consistent critical values
    for coverage; runtime/hardware disclosure; move the
    simulation out of the cached chunk into a saved-results
    workflow.
12. Replace the penguin data README; wire `sim_08.rds` into the
    pipeline or delete it; add substantive tinytest coverage of
    the DGM and summary functions.
13. Terminology cleanup (MMRM versus random-slopes LMM) in the
    title and throughout; remove line numbers and TOC at
    submission.

## 5. Recommended framing

Plausible framings for this paper:

(a) A neutral simulation-comparison paper ("how robust is the
slope estimand to nonlinearity"), extending Chen et al. (2018)
with systematic control of the nonlinearity parameter and
Morris-compliant reporting.

(b) An estimand-centered regulatory piece: slope versus
final-visit estimands for disease-modification claims under ICH
E9(R1), with the simulation as supporting evidence.

(c) A tutorial/software paper around a reusable simulation
framework (the "Future research" list points this way, but no
package functions exist).

Recommendation: framing (a), as a compact simulation study for
Pharmaceutical Statistics or Statistics in Biopharmaceutical
Research, with the estimand discussion of (b) folded into the
introduction and discussion rather than leading. Reasoning: the
literature already contains qualitative warnings that slope
models are biased under nonlinearity (Chen et al. 2018;
Liu-Seifert et al. 2015) and regulators already condition slope
claims on approximate linearity, so a purely rhetorical estimand
piece adds little; framing (c) is unavailable because no
software exists. The genuine gap is quantitative: how much
nonlinearity produces how much bias, coverage loss, and power
change at a fixed 18-month estimand, reported with full ADEMP
discipline. Note also that the study's most publishable finding
may currently be mislabeled: the observation that the
random-slopes estimator fails to deliver a power advantage even
under exact linearity (Table 1) is surprising relative to
received wisdom and, once its cause is pinned down (estimator
definition, exclusion of baseline), could become a headline
result rather than an embarrassment.

Implications of framing (a): the title should name the
comparison and endpoint ("Bias and power of slope versus
final-visit estimands under nonlinear treatment effects in
Alzheimer's disease trials"); the abstract should lead with the
quantitative tolerance result (e.g., coverage falls to 0.815 at
kappa = 0.3); the introduction should compress the regulatory
material to one paragraph; comparators should include cLDA and
the proportional MMRM (currently relegated to future work; at
least cLDA is needed for the comparison to reflect practice);
the ADEMP audit, penguin README, and compliance section move out
of the manuscript entirely; the shape-trajectory figure stays,
and full MCSE tables go to supplementary material with
parenthetical MCSEs in the main table.

## 6. Assessment

Verdict: major revision (were this submitted today, the
contradiction between the Discussion and Table 1 plus the absent
type I error assessment would risk outright rejection). The
compendium's infrastructure is genuinely strong: the simulation
script is clean, seeds and RNGkind are handled per Morris et al.
(2019), MCSEs are computed, and every number in the manuscript
is generated inline from code (inspected). But the manuscript
around it is unfinished: expectation-language where results
belong, Methods that misdescribe the code, a stale self-audit
pasted after the References heading, no abstract, a
single-design simulation the paper itself flags as inadequate,
and a bibliography with demonstrable accuracy problems. The
scientific question is well chosen and the quantitative answer
is publishable once the estimator-definition issue is resolved
and the design is extended; the distance to submission is
substantial but well defined.

## 7. Revision history

- 2026-08-16: Initial referee-grade review. Eight major and
  eleven minor issues identified; verdict major revision. No
  prior whitepaper existed.
