# run_sim_08.R
# Driver script: runs the full-scale linearity-robustness simulation
# and saves both the raw per-replicate results and the summarized
# performance measures to analysis/data/derived_data/sim_08.rds,
# which analysis/report/report.Rmd reads directly (rather than
# re-running ~7000-10000 mixed-model fits inside a cached knitr
# chunk at render time; see pub_review whitepaper 2026-08-16, minor
# issue 5).
#
# Usage (from the mmrmrobust package root):
#   Rscript analysis/scripts/run_sim_08.R
#
# Wall-clock time: the n_rep = 1000, 7-kappa design below fits
# roughly 21,000 mixed models (unconstrained slope, constrained
# slope, categorical MMRM) and takes several minutes on a laptop.
# This is too long to run inside an interactive remediation pass;
# it must be run by the user to refresh sim_08.rds after any change
# to the data-generating model or fit_models().

here_root <- function() {
  d <- getwd()
  for (i in seq_len(8)) {
    if (file.exists(file.path(d, "DESCRIPTION")) &&
        file.exists(file.path(d, "analysis", "scripts",
          "sim_mmrm_linearity.R"))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Could not locate the mmrmrobust package root from ", getwd())
}

root <- here_root()
source(file.path(root, "analysis", "scripts", "sim_mmrm_linearity.R"))

# Morris, White, and Crowther (2019) Section 4.1: pin RNGkind and set
# the seed once at program start; matches analysis/report/report.Rmd.
RNGkind("L'Ecuyer-CMRG")
set.seed(20260310)

n_sim_reps <- 1000
sim_raw <- run_simulation(
  n_rep     = n_sim_reps,
  kappa_vec = c(0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0),
  n_per_arm = 200,
  delta     = -0.45
)
sim_summary <- summarize_results(sim_raw, true_delta = -0.45)

out_path <- file.path(root, "analysis", "data", "derived_data",
  "sim_08.rds")
saveRDS(list(res = sim_raw, summary = sim_summary), out_path)
message("Saved ", out_path)
