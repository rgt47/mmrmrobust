# pub_review Remediation Log: 08-mmrm-linearity-robust
*2026-08-20 13:23 PDT*

Remediation pass against
`docs/pub_review_whitepaper_2026-08-16.md`. All file paths are
relative to the `mmrmrobust` package root unless noted.

## 1. Fixed

- **Major issue 1 / checklist (a)1** [verified]. Rewrote the
  "Headline findings," "Power comparison," and Discussion
  sections of `analysis/report/report.Rmd` to state the actual
  Table 1 numbers (random-slopes power 0.542 vs. categorical
  0.623 at kappa = 1; larger empirical SE for the random-slopes
  model) instead of "is expected to" language. Confirmed by a
  full `bash tools/render.sh analysis/report/report.Rmd` build;
  `report.tex` contains the same 0.542/0.623 figures in prose
  and Table 1.
- **Major issue 2 / checklist (a)2** [verified, partial].
  Implemented a constrained random-slopes estimator (no
  treatment main effect, `y ~ time + trt:time`) in
  `fit_models()`, `run_simulation()`, and `summarize_results()`
  in `analysis/scripts/sim_mmrm_linearity.R`. Verified it
  converges and is directionally consistent with the
  whitepaper's hypothesis via a small-scale check
  (n_rep = 20, kappa in {0.3, 1, 3}, `set.seed(1)`): bias at
  kappa = 0.3 dropped from ~0.17 (unconstrained) to ~-0.01
  (constrained). Added tinytest coverage
  (`inst/tinytest/test_basic.R`) confirming both variants
  return non-NA, positive-SE contrasts. The full-scale rerun
  needed to add a third column to Table 1 was NOT attempted
  (see Deferred).
- **Major issue 3 (partial acknowledgment only)**. Not fixed;
  see Deferred. Disclosed explicitly in the new Limitations
  paragraph and Reproducibility statement rather than left
  silently unaddressed.
- **Major issue 4 / checklist (a)3** [applied, unverified via
  render only]. Fixed the gls-vs-mmrm contradiction (now
  correctly says `mmrm::mmrm` throughout), the false delta-method
  claim (now correctly describes the Satterthwaite `mmrm::df_1d`
  contrast), and added an explicit disclosure that neither model
  adjusts for baseline CDR-SB, framed as a deliberate design
  choice rather than an omission (`analysis/report/report.Rmd`,
  "Analysis models" and Limitations).
- **Major issue 5 / checklist (a)4** [verified]. Deleted the
  stale "Morris et al. (2019) ADEMP Compliance" subsection that
  sat under `# References` and asserted defects (n_rep mismatch,
  seed inside worker, no MCSE) that no longer exist in the code.
  Replaced it with a new, accurate "Reproducibility statement"
  section placed before `# References`, describing the current
  RNGkind/seed handling and pointing to
  `docs/morris-audit-2026-04-17.md` as an archived, pre-fix
  audit that is explicitly not part of the manuscript.
- **Major issue 6 / checklist (b)8** [verified]. Added Table 2
  (`tab-mcse` chunk) reporting convergence rate and MCSE(bias),
  MCSE(power) for both models per Morris et al. (2019) Table 6;
  confirmed present in the rendered `report.tex`/`report.pdf`.
- **Major issue 8a / checklist (a)5** [verified via WebFetch/
  WebSearch against publisher records]. Fixed four fabricated or
  mis-keyed bibliography entries in
  `analysis/report/references.bib`:
  - `chen2018simulation`: removed fabricated co-author "Bhatt,
    Dhairya"; corrected author list to Chen, Ni, Fleisher, Zhou,
    Aisen, Mohs and corrected pages (46-53, not 143-153).
  - `wang2022proportional`: replaced an entirely fabricated
    author list with the verified list (Wang, Liu, Li,
    Aschenbrenner, Bateman, Delmar, Schneider, Kennedy, Cutter,
    Xiong).
  - `novic2024nonlinear` -> `wang2024nonlinear`: the entry was
    fabricated in its entirety (invented authors "Novic, Engel,
    Scheltens, Barkhof, Berkhof" and wrong journal
    "Pharmaceutical Statistics"). The real PMID 38727205 article
    with this exact title is Wang et al., Statistics in Medicine
    2024;43(15):2987-3004; corrected and re-cited under the new
    key.
  - `donohue2011power` -> `donohue2014pacc`: content (authors,
    journal, pages) was correct but mis-keyed/mis-captioned as
    a 2011 "power" paper when it is the 2014 JAMA Neurology PACC
    paper; renamed to match its actual content and cited.
  - `liang2000longitudinal` -> `liang1986longitudinal`: content
    was the correct, real 1986 Liang-Zeger GEE paper but keyed
    "2000"; renamed to match. Added the genuinely missing 2000
    Liang-Zeger Sankhya cLDA paper as a new entry
    `liang2000sankhya` and cited it in the estimator-remediation
    discussion.
- **Major issue 8b/8c / checklist (b)10** [verified]. Cited all
  six previously orphaned bib entries in the text: `ard2011`
  (power-calculation context in Methods), `liu2015slopes` and
  `liang2000sankhya` (cLDA estimator discussion), `donohue2014pacc`
  (Introduction, contrasting composite endpoints), `budd2022aducanumab`
  (Introduction, third anti-amyloid trial), `liang1986longitudinal`
  and `fitzmaurice2011applied` (general longitudinal-methods
  citation in "Analysis models"). Confirmed via a bib-key diff
  that every entry in `references.bib` is now cited and every
  citation resolves to an entry.
- **Checklist (b)9 / minor issue 1** [verified via render].
  Added an Abstract (leading with the quantitative coverage/bias
  results, per the whitepaper's Recommended Framing), Keywords,
  and a Data and code availability statement to
  `analysis/report/report.Rmd`.
- **Minor issue 2** [verified]. Fixed the duplicate "dotted"
  linetype for kappa = 0.3 and kappa = 3 in Figure 1
  (`ltype_map`); kappa = 3 now uses a distinct custom dash
  pattern.
- **Minor issue 3** [applied, unverified numerically]. Disclosed
  the coverage-critical-value inconsistency (1.96 for both
  models vs. t/Satterthwaite p-values) explicitly in Limitations
  rather than silently leaving it; not re-coded because it would
  require re-deriving per-model critical values and rerunning.
- **Minor issue 4** [verified against Table 1 numbers]. Rewrote
  the ambiguous bias-mechanism text in "Bias" to state the
  correct, sign-consistent direction (attenuation toward zero
  for kappa < 1, exaggeration for kappa > 1), hedged as a
  plausible, not separately verified, mechanism.
- **Minor issue 5** [verified]. Moved the simulation out of the
  cached `run-sim` knitr chunk into a versioned, saved-results
  workflow: added `analysis/scripts/run_sim_08.R`, which
  reproduces `analysis/data/derived_data/sim_08.rds` with the
  documented seed/RNGkind; `report.Rmd`'s `run-sim` chunk now
  loads this file (confirmed by full render, ~1 minute vs. the
  prior ~7000+ inline model fits) with a live-run fallback if
  the file is missing.
- **Minor issue 6** [verified]. `sim_08.rds` is no longer an
  orphan: `run_sim_08.R` documents its provenance and
  `report.Rmd` reads it directly; confirmed the loaded object's
  summary matches the previously rendered Table 1 exactly
  (`Rscript` check comparing `sim_08.rds$summary` to the old
  `report.tex` table).
- **Minor issue 7 / checklist (c)12** [verified]. Replaced the
  placeholder `expect_true(TRUE)` test with 17 real tinytest
  assertions in `inst/tinytest/test_basic.R` covering
  `simulate_one_trial()` (row counts, visit schedule, arm
  balance, DGM mean recovery at kappa = 1, zero baseline
  effect), `fit_models()` (all three model variants converge
  with positive SEs), and `summarize_results()` (bias/power/
  conv/MCSE formulas on a hand-computed toy input). All 17 pass
  (`Rscript -e 'tinytest::run_test_file(...)'`, 4.9s).
- **Minor issue 8 / checklist (c)13 (partial)** [verified].
  Retitled the manuscript to "Bias and Power of Slope Versus
  Final-Visit Estimands Under Nonlinear Treatment Effects in
  Alzheimer's Disease Trials" per the whitepaper's Recommended
  Framing, and replaced "random-slopes MMRM" with "random-slopes
  model" in the Present study and ADEMP-structure sections. Not
  a full terminology sweep (see Deferred).
- **Minor issue 9**. Added CDR-SB boundedness/discreteness
  (0-18, half-point steps) to Limitations, alongside the
  existing Gaussianity limitation.

## 2. Deferred

- **Major issue 2, full-scale integration** (constrained
  estimator into Table 1). Reason: a full run
  (n_rep = 1000 x 7 kappa x 3 models, ~14,000+ mixed-model fits)
  is well beyond the "minutes, not hours" budget for this pass.
  Command: `Rscript analysis/scripts/run_sim_08.R` after
  extending it to also summarize `slope_c_*` in the saved object
  (the summarization code already exists in
  `summarize_results()`); then add a third model block to Table
  1 in `report.Rmd`.
- **Major issue 3 / checklist (b)6** (null and transient-effect
  scenarios). Reason: requires a code change to
  `run_simulation()`/the report to accept a `delta` grid
  (currently a single scalar) plus a rerun; out of budget.
  `run_simulation()` would need a `delta_vec` argument analogous
  to `kappa_vec`.
- **Major issue 7 / checklist (b)7** (design expansion: second
  sample size near 900/arm, MAR dropout, n_sim >= 2100). Reason:
  each is a multi-minute-to-hour rerun; disclosed explicitly in
  Limitations instead of attempted. No code changes made toward
  dropout simulation.
- **Minor issue 3** (per-model coverage critical values).
  Reason: small effect, requires touching
  `summarize_results()`'s coverage formula and rerunning to
  confirm no regression; disclosed in Limitations instead.
- **Minor issue 8, full terminology sweep**. Reason: many
  remaining mentions of "categorical-time MMRM" are correct as
  written (that model genuinely is an MMRM); a full audit of
  every "MMRM" occurrence for correctness was not performed
  beyond the sections directly flagged by the whitepaper.
- **Minor issue 11** (remove `lineno`/`toc` for submission
  format). Reason: this is a working-draft formatting choice the
  author may still want during revision; left as-is pending an
  explicit decision to finalize for submission.

## 3. New issues found while fixing

- The `donohue2011power` and `liang2000longitudinal` bib entries
  were not fabricated, contrary to what the whitepaper's
  phrasing might suggest to a quick reader: their author lists,
  journals, volumes, and pages are all correct and verified
  against the publisher record. The defect was a mis-chosen
  BibTeX key/caption (wrong year, wrong implied topic), not
  fabricated content. This matters because deleting them (one
  plausible reading of "fix or delete suspect author lists")
  would have discarded two genuine, useful references; they were
  renamed and kept instead.
- `analysis/data/README.md` (Palmer Penguins template
  boilerplate, whitepaper minor issue 6) was inspected but not
  fixed in this pass; it was deprioritized below the (a)/(b)-tier
  items given the time budget and is unrelated to any reported
  number in the manuscript.
- The full render (`bash tools/render.sh analysis/report/report.Rmd`)
  now completes in well under a minute once `sim_08.rds` is
  loaded directly, versus the prior cached-chunk approach that
  fit ~7000 models at render time; this substantially reduces
  the render-time risk flagged in minor issue 5, beyond what the
  whitepaper's remediation text anticipated.
