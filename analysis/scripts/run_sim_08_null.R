# run_sim_08_null.R
# Driver script: runs the null-scenario (delta = 0) simulation to
# report empirical type I error for all three analysis models
# (unconstrained random-slopes, constrained random-slopes,
# categorical-time MMRM), addressing pub_review whitepaper
# 2026-08-16 major issue 3 / checklist (b)6 ("Power comparisons
# between a misspecified and a correctly specified model are
# difficult to interpret fully without also demonstrating type I
# error control").
#
# Because the true treatment effect is delta = 0 at every visit
# regardless of the shape parameter kappa (g(t, kappa) * 0 = 0 for
# all kappa and t), the trajectory shape is irrelevant under the
# null and a single kappa is sufficient; kappa = 1.0 is used as a
# nominal placeholder value only (not meaningful under delta = 0).
#
# Usage (from the mmrmrobust package root):
#   Rscript analysis/scripts/run_sim_08_null.R
#
# Saves analysis/data/derived_data/sim_08_null.rds, loaded directly
# by analysis/report/report.Rmd (Results, "Type I error under the
# null").

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
# a seed distinct from the main-design run (analysis/scripts/
# run_sim_08.R uses 20260310) so the null-scenario RNG stream never
# overlaps the main design's stream.
RNGkind("L'Ecuyer-CMRG")
set.seed(20260320)

n_sim_reps <- 1000
mc_cores <- max(1, parallel::detectCores() - 1)
message("Running null-scenario simulation on ", mc_cores, " cores...")
sim_null_raw <- run_simulation_parallel(
  n_rep     = n_sim_reps,
  kappa_vec = 1.0,
  mc.cores  = mc_cores,
  n_per_arm = 200,
  delta     = 0
)
sim_null_summary <- summarize_results(sim_null_raw, true_delta = 0)

out_path <- file.path(root, "analysis", "data", "derived_data",
  "sim_08_null.rds")
saveRDS(list(res = sim_null_raw, summary = sim_null_summary), out_path)
message("Saved ", out_path)
