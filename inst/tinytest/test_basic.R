library(tinytest)

# The simulation driver lives in analysis/scripts/ (excluded from the
# built package via .Rbuildignore, per the compendium's zzcollab
# layout), so it is not installed with mmrmrobust. Locate it by
# searching upward from the working directory for the compendium
# root; skip this file gracefully (rather than fail) when the source
# tree is not available, e.g. inside a bare `R CMD check` sandbox.
find_sim_script <- function() {
  d <- getwd()
  for (i in seq_len(8)) {
    candidate <- file.path(d, "analysis", "scripts",
      "sim_mmrm_linearity.R")
    if (file.exists(candidate)) return(candidate)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  NA_character_
}

sim_script <- find_sim_script()
if (is.na(sim_script)) {
  exit_file("sim_mmrm_linearity.R not found; skipping simulation tests")
}
source(sim_script)

# --- simulate_one_trial(): data-generating mechanism -----------------

set.seed(20260819)
dat_null <- simulate_one_trial(
  n_per_arm = 50, delta = 0, kappa = 1.0
)

expect_equal(
  nrow(dat_null), 700,
  info = "simulate_one_trial returns n_per_arm * 2 * n_visits rows"
)
expect_true(
  all(c("id", "trt", "visit", "time_yr", "y", "trt_f", "visit_f") %in%
    names(dat_null)),
  info = "simulate_one_trial returns the columns fit_models() expects"
)
expect_equal(
  sort(unique(dat_null$visit)), c(0, 3, 6, 9, 12, 15, 18),
  info = "visit schedule matches the documented 3-month intervals"
)
expect_equal(
  as.numeric(table(dat_null$trt)), c(50, 50) * 7,
  info = "arms are balanced at n_per_arm each, replicated over visits"
)

set.seed(20260819)
dat_treat <- simulate_one_trial(
  n_per_arm = 4000, delta = -0.45, kappa = 1.0
)
mean_diff_18 <- with(
  dat_treat[dat_treat$visit == 18, ],
  mean(y[trt == 1]) - mean(y[trt == 0])
)
expect_true(
  abs(mean_diff_18 - (-0.45)) < 0.05,
  info = paste(
    "at kappa = 1 the simulated month-18 group mean difference",
    "recovers the true delta = -0.45 to within Monte Carlo error"
  )
)

dat_baseline <- simulate_one_trial(
  n_per_arm = 4000, delta = -0.45, kappa = 1.0
)
baseline_diff <- with(
  dat_baseline[dat_baseline$visit == 0, ],
  mean(y[trt == 1]) - mean(y[trt == 0])
)
expect_true(
  abs(baseline_diff) < 0.05,
  info = "the DGM has no treatment effect at baseline (t = 0) by construction"
)

# --- fit_models(): both models return usable contrasts ----------------

set.seed(20260819)
dat_fit <- simulate_one_trial(n_per_arm = 200, delta = -0.45, kappa = 1.0)
fits <- fit_models(dat_fit)

expect_true(
  all(c("slope", "slope_c", "cat") %in% names(fits)),
  info = paste(
    "fit_models returns the unconstrained slope, constrained",
    "slope, and categorical-time contrasts"
  )
)
expect_false(
  is.na(fits$slope["est"]),
  info = "unconstrained random-slopes model converges on well-behaved data"
)
expect_false(
  is.na(fits$slope_c["est"]),
  info = "constrained random-slopes model (no main effect) converges"
)
expect_false(
  is.na(fits$cat["est"]),
  info = "categorical-time MMRM converges on a well-behaved dataset"
)
expect_true(
  fits$slope["se"] > 0 && fits$slope_c["se"] > 0 && fits$cat["se"] > 0,
  info = "all three fitted contrasts report a strictly positive SE"
)

# --- summarize_results(): MCSE formulas and bias definition -----------

toy_res <- data.frame(
  kappa        = rep(1, 4),
  rep          = 1:4,
  slope_est    = c(-0.4, -0.5, -0.4, -0.5),
  slope_se     = c(0.2, 0.2, 0.2, 0.2),
  slope_pval   = c(0.01, 0.20, 0.01, 0.20),
  slope_conv   = c(1, 1, 1, 1),
  slope_c_est  = c(-0.45, -0.45, -0.45, -0.45),
  slope_c_se   = c(0.2, 0.2, 0.2, 0.2),
  slope_c_pval = c(0.01, 0.01, 0.01, 0.01),
  slope_c_conv = c(1, 1, 1, 1),
  cat_est      = c(-0.45, NA, -0.45, -0.45),
  cat_se       = c(0.2, NA, 0.2, 0.2),
  cat_pval     = c(0.01, NA, 0.01, 0.20),
  cat_conv     = c(1, 0, 1, 1)
)
summ <- summarize_results(toy_res, true_delta = -0.45)

expect_equal(
  summ$slope_bias, mean(c(-0.4, -0.5, -0.4, -0.5)) - (-0.45),
  info = "slope_bias is mean(estimate) - true_delta"
)
expect_equal(
  summ$slope_c_bias, 0,
  info = "constrained slope model recovers true_delta on this toy input"
)
expect_equal(
  summ$cat_conv, 3 / 4,
  info = "cat_conv is the proportion of reps with a non-missing estimate"
)
expect_equal(
  summ$cat_n_conv, 3,
  info = "cat_n_conv counts non-NA estimates, excluding the failed rep"
)
expect_equal(
  summ$slope_power, mean(c(TRUE, FALSE, TRUE, FALSE)),
  info = "slope_power is the proportion of reps with p < 0.05"
)
expect_true(
  all(c(
    "mcse_slope_bias", "mcse_slope_c_bias", "mcse_cat_bias",
    "mcse_slope_power", "mcse_cat_power",
    "mcse_slope_coverage", "mcse_cat_coverage"
  ) %in% names(summ)),
  info = "summarize_results reports MCSEs per Morris et al. (2019) Table 6"
)

# --- run_simulation_parallel(): parallel driver used for the -------
# --- full-scale reruns (analysis/scripts/run_sim_08.R, ------------
# --- run_sim_08_null.R) --------------------------------------------

RNGkind("L'Ecuyer-CMRG")
set.seed(20260819)
par_res <- run_simulation_parallel(
  n_rep = 2, kappa_vec = c(1.0, 3.0), mc.cores = 2,
  n_per_arm = 30, delta = -0.45
)

expect_equal(
  nrow(par_res), 4,
  info = "run_simulation_parallel returns one row per (kappa, rep) combination"
)
expect_equal(
  sort(unique(par_res$kappa)), c(1, 3),
  info = "run_simulation_parallel covers every requested kappa value"
)
expect_true(
  all(c(
    "slope_est", "slope_c_est", "cat_est",
    "slope_se", "slope_c_se", "cat_se"
  ) %in% names(par_res)),
  info = paste(
    "run_simulation_parallel returns the same columns as the",
    "serial run_simulation()"
  )
)
expect_true(
  all(!is.na(par_res$slope_c_est)),
  info = "run_simulation_parallel's constrained-slope estimates converge"
)
