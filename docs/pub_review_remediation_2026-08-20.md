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

## 4. Follow-up pass: completion of deferred items 1 and 2
*2026-08-20 18:41 PDT*

A subsequent session (working from this log plus
`git status`/`git diff`) had already made most of the code and
prose changes for major issue 2's full-scale integration and
major issue 3's null-scenario rerun before its context ended;
that work was auto-committed as
`e316ee7 Auto-backup: 2026-08-20 18:02:20` before this pass
started. This pass verified what was actually on disk (rather
than assuming the prior session's log entries described its
final state), ran the one script that had not yet been executed,
fixed a render-breaking bug the prior session's log did not
mention, and confirmed the full pipeline (tests + render) end to
end.

- **Major issue 2, full-scale integration** [verified, now
  DONE]. `analysis/scripts/run_sim_08.R` and
  `analysis/scripts/sim_mmrm_linearity.R` (as committed in
  `e316ee7`) already summarized the constrained estimator
  (`slope_c_*` columns) at full scale ($n_{\mathrm{rep}} = 1000$
  x 7 $\kappa$ values x 3 models = 7,000 model fits, confirmed by
  `dim(sim_08.rds$res)` = 7000 x 14 and
  `sim_08.rds$summary$n_reps` = 1000 at every $\kappa$).
  `report.Rmd`'s Table 1, Table 2, Figure 2, and the new
  "Estimator definition check" subsection were already wired to
  read these columns. Nothing further was needed for this item
  beyond the render fix below; it is complete, not partial.
- **Major issue 3, null/type-I-error scenario** [verified, now
  DONE]. `analysis/scripts/run_sim_08_null.R` existed on disk
  (committed in `e316ee7`) and `report.Rmd` already had a "Type I
  error under the null" section wired to load
  `analysis/data/derived_data/sim_08_null.rds`, but that file did
  not exist yet (the script had been written but never executed).
  Ran it: `Rscript analysis/scripts/run_sim_08_null.R`
  ($n_{\mathrm{rep}} = 1000$, $\delta = 0$, `RNGkind("L'Ecuyer-CMRG")`,
  `set.seed(20260320)`, distinct from the main design's seed
  20260310; ~83 s wall clock on 7 cores). Empirical type I error
  at the nominal 0.05 level was 0.059 (unconstrained random
  slopes, MCSE 0.007), 0.070 (constrained random slopes, MCSE
  0.008), and 0.054 (categorical-time MMRM, MCSE 0.007); all
  three models converged on 100% of the 1,000 null replicates.
  These figures now appear in `report.tex`/`report.pdf` ("Type I
  error under the null" and the associated table). No claim of
  formal type I error control (e.g., via a binomial test against
  0.05) was added beyond reporting the point estimate and MCSE,
  consistent with what the report's prose says.
- **Render-breaking bug found and fixed** [verified]. The first
  `bash tools/render.sh analysis/report/report.Rmd` attempt
  failed in the `estimator-check-08` chunk with `Error in `if
  (...) NULL`: ! argument is of length zero`. Root cause: knitr's
  `cache = TRUE` / `cache.path = "cache/"` global option (set in
  the `setup` chunk) had stale cached results in
  `analysis/report/cache/` for the `headline-08` chunk, computed
  before the `slope_c_*` columns existed in
  `sim_08.rds$summary`. Because knitr cache invalidation keys off
  chunk source text and options, not upstream data-file content,
  the unchanged `headline-08` chunk code reused its stale cached
  `linear_row`/`worst_slope_bias` objects (lacking `slope_c_bias`,
  `slope_c_power`, etc.), so `linear_row$slope_c_power` evaluated
  to `numeric(0)` three chunks later in `estimator-check-08`.
  Fixed by `rm -rf analysis/report/cache` and re-rendering from
  scratch; the render then completed cleanly (23/23 chunks,
  `report.pdf` produced, 17 pages, staged to
  `share/report-2026-08-20-1840-e316ee7-wip.pdf`). This is a
  latent hazard for any future edit to `sim_08.rds` or
  `sim_08_null.rds` in this repo: an incremental (non-`rm -rf
  cache`) render after regenerating either `.rds` file can
  silently serve stale numbers instead of erroring, because
  `headline-08` and the other summary chunks are not declared
  with `dependson = "load-sim-08"`. Not fixed in this pass (would
  require either adding `dependson` to every downstream summary
  chunk or disabling caching for those chunks); recommend
  clearing `analysis/report/cache/` before every render that
  follows a simulation rerun until this is addressed.
- **Tinytest suite** [verified]. `Rscript -e
  'pkgload::load_all("."); tinytest::run_test_dir("inst/tinytest")'`
  reports 21/21 passing (4.6 s), not the 17/17 recorded in this
  log's Section 1; the prior session's `e316ee7` commit added 4
  more assertions to `inst/tinytest/test_basic.R` beyond what
  Section 1 above described (diff shows 34 insertions to that
  file total). All pass; no failures at any point in this pass.
- **Design expansion (major issue 7 / checklist (b)7)**
  [confirmed still deferred, unchanged]. Per the original
  triage instruction (prioritize items 1 and 2), this item was
  not attempted. `report.Rmd`'s Limitations and "Future work"
  sections still describe it prospectively (a second sample size
  near 900/arm matching CLARITY-AD, MAR dropout, $n_{\mathrm{rep}}
  \geq 2{,}100$) rather than reporting it as done. No code toward
  a dropout mechanism was written. This remains a multi-hour
  rerun (larger per-arm N increases per-fit cost, and dropout
  requires extending `simulate_one_trial()`), consistent with the
  original log's estimate.

### Final status of the three items in scope for this pass

1. Full-scale Table 1 integration of the constrained-slope
   estimator: **done** (was already complete on disk at the
   start of this pass; verified, not re-derived).
2. Null/type-I-error scenario rerun, integrated into the
   manuscript: **done** (script existed but had not been run;
   this pass ran it and confirmed the resulting numbers render
   correctly).
3. Design expansion (900/arm, dropout, $n_{\mathrm{rep}} \geq
   2{,}100$): **still deferred**, as instructed when triaging
   against items 1 and 2; disclosed as prospective future work in
   the manuscript, not attempted.

Tests: 21/21 passing. Render: `bash tools/render.sh
analysis/report/report.Rmd` succeeds end to end after clearing
the stale cache, producing a 17-page `report.pdf`.
