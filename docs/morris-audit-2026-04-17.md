# Morris et al. (2019) ADEMP Audit: 08-mmrm-linearity-robust
*2026-04-17 09:02 PDT*

## Scope

Files audited:

- `analysis/scripts/sim_mmrm_linearity.R`
- `analysis/report/report.Rmd`

## ADEMP scorecard

| Criterion | Status | Evidence |
|---|---|---|
| Aims explicit | Partial | aims described in prose, no ADEMP heading |
| DGMs documented | Met | `sim_mmrm_linearity.R` parameterises kappa, n, visits |
| Factors varied factorially | Partial | kappa × n grid implicit, not labelled |
| Estimand defined with true value | Met | treatment-by-time interaction identified |
| Methods justified | Met | MMRM, linear-trend MMRM, ANCOVA listed |
| Performance measures justified | Partial | bias / coverage / power reported; not traced to aims |
| n_sim stated | Partial | report says `n_rep = 200` at `report.Rmd:261`; script default is `n_rep = 1000` — mismatch |
| n_sim justified via MCSE | Not met | no MCSE derivation |
| MCSE reported per metric | Not met | `summarize_results()` returns no MCSE columns |
| Seed set once | Not met | `set.seed(20260310)` lives INSIDE `run_simulation()` at `sim_mmrm_linearity.R:178` — reseeded on each call |
| RNG states stored | Not met | no `.Random.seed` capture |
| Paired comparisons | Met | all three methods applied to same simulated dataset per rep |
| Reproducibility | Partial | seed exists; RNGkind not pinned; n_rep mismatch between Rmd and script causes silent value drift |

## Overall verdict

**Partially compliant.**

## Gaps

- **n_rep mismatch** between `report.Rmd:261` (`n_rep = 200`) and the
  script default (`n_rep = 1000`) — whichever the Rmd passes wins, so the
  report is likely running with 200 reps, not the script-documented 1000.
- Seed is set inside the worker function (`sim_mmrm_linearity.R:178`), so
  `run_simulation()` can be called in other contexts that reset RNG state
  mid-pipeline. Morris §4.1: set the seed ONCE at program start.
- No MCSE on any performance metric.
- Convergence rates per kappa are not reported (Morris §5.1).
- `RNGkind("L'Ecuyer-CMRG")` not pinned.

## Remediation plan

1. Remove `set.seed()` from inside `run_simulation()` in
   `sim_mmrm_linearity.R:178`. Set the seed once at the top of
   `analysis/scripts/sim_mmrm_linearity.R` immediately before the master
   loop.
2. Align `n_rep`: pick one canonical value (1000 or higher), parameterise
   it in a single config location, and import from both the script and
   the Rmd chunk to prevent drift.
3. Derive the required n_rep from a target MCSE at the top of the script
   (e.g., power MCSE ≤ 1 pp requires n_rep ≥ 2500 at power 0.5).
4. Add `mcse_*` columns to `summarize_results()` for bias, empirical SE,
   coverage, and power per Morris Table 6.
5. Emit a `convergence_rate` column per kappa scenario.
6. Pin `RNGkind("L'Ecuyer-CMRG")` and store `.Random.seed` per replicate.
7. Add ADEMP-headed Methods section to `report.Rmd`.

## References

Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate
statistical methods. Stat Med 2019;38:2074-2102. doi:10.1002/sim.8086

---
*Source: ~/prj/res/08-mmrm-linearity-robust/mmrmrobust/docs/morris-audit-2026-04-17.md*
