# sim_mmrm_linearity.R
# Simulation: robustness of random-slopes MMRM to nonlinear
# treatment effects in AD trials (CDR-SB, 18-month, 3-month
# intervals)
#
# Data-generating model:
#   Y_ij = beta0 + b0_i
#          + (beta1 + b1_i) * t_j
#          + delta * g(t_j, kappa) * trt_i
#          + epsilon_ij
#
# where g(t, kappa) parameterizes the treatment effect shape:
#   kappa = 1  => linear:   g(t) = t / t_max
#   kappa < 1  => early:    g(t) = (t / t_max)^kappa
#   kappa > 1  => delayed:  g(t) = (t / t_max)^kappa
#
# The total treatment effect at t_max is always delta,
# regardless of kappa. Only the *path* changes.

library(MASS)
library(nlme)

simulate_one_trial <- function(
    n_per_arm    = 200,
    visits       = c(0, 3, 6, 9, 12, 15, 18),
    beta0        = 3.2,
    beta1_yr     = 1.5,
    delta        = -0.45,
    kappa        = 1.0,
    sd_b0        = 1.0,
    sd_b1        = 0.6,
    rho_b        = 0.3,
    sigma        = 1.2,
    ar1_phi      = 0.4
) {
  t_max <- max(visits)
  t_months <- visits
  n_visits <- length(t_months)
  n_total <- 2 * n_per_arm

  trt <- rep(c(0, 1), each = n_per_arm)

  D <- matrix(
    c(sd_b0^2, rho_b * sd_b0 * sd_b1,
      rho_b * sd_b0 * sd_b1, sd_b1^2),
    nrow = 2
  )
  b <- mvrnorm(n_total, mu = c(0, 0), Sigma = D)
  b0_i <- b[, 1]
  b1_i <- b[, 2]

  ar1_cor <- ar1_phi^abs(
    outer(seq_len(n_visits), seq_len(n_visits), "-")
  )
  R <- sigma^2 * ar1_cor

  dat_list <- vector("list", n_total)
  for (i in seq_len(n_total)) {
    eps_i <- mvrnorm(1, mu = rep(0, n_visits), Sigma = R)
    t_yr <- t_months / 12
    g_t <- (t_months / t_max)^kappa
    y <- beta0 + b0_i[i] +
      (beta1_yr + b1_i[i]) * t_yr +
      delta * g_t * trt[i] +
      eps_i
    dat_list[[i]] <- data.frame(
      id      = i,
      trt     = trt[i],
      visit   = t_months,
      time_yr = t_yr,
      y       = y
    )
  }
  dat <- do.call(rbind, dat_list)
  dat$trt_f <- factor(dat$trt, labels = c("placebo", "active"))
  dat$visit_f <- factor(dat$visit)
  dat$id <- factor(dat$id)
  dat$visit_num <- as.integer(dat$visit_f)

  dat
}

fit_models <- function(dat) {
  t_max <- max(dat$visit)
  post <- dat[dat$visit > 0, ]
  post$visit_f <- droplevels(post$visit_f)
  post$visit_num <- as.integer(post$visit_f)

  # Model 1: random-slopes (linear in time)
  fit_slope <- tryCatch({
    lme(
      y ~ trt_f * time_yr,
      random = ~ time_yr | id,
      data = post,
      control = lmeControl(
        opt = "optim",
        maxIter = 200,
        msMaxIter = 200
      )
    )
  }, error = function(e) NULL)

  # Model 2: standard MMRM -- marginal model,
  # unstructured covariance via gls
  fit_cat <- tryCatch({
    gls(
      y ~ trt_f * visit_f,
      data = post,
      correlation = corSymm(form = ~ visit_num | id),
      weights = varIdent(form = ~ 1 | visit_f),
      control = glsControl(
        opt = "optim",
        maxIter = 200,
        msMaxIter = 200
      )
    )
  }, error = function(e) NULL)

  results <- list()

  # --- Extract slope model results ---
  if (!is.null(fit_slope)) {
    cf <- summary(fit_slope)$tTable
    row_nm <- "trt_factive:time_yr"
    if (row_nm %in% rownames(cf)) {
      t_max_yr <- t_max / 12
      est <- cf[row_nm, "Value"] * t_max_yr
      se  <- cf[row_nm, "Std.Error"] * t_max_yr
      pv  <- cf[row_nm, "p-value"]
      results$slope <- c(
        est = est, se = se, pval = pv, converged = 1
      )
    } else {
      results$slope <- c(
        est = NA, se = NA, pval = NA, converged = 0
      )
    }
  } else {
    results$slope <- c(
      est = NA, se = NA, pval = NA, converged = 0
    )
  }

  # --- Extract categorical MMRM results ---
  # The coefficient trt_factive:visit_f18 is the
  # difference-in-differences: (active_18 - active_ref) -
  # (placebo_18 - placebo_ref). Since visit_f ref = 3
  # (first post-baseline), this equals the treatment
  # effect at month 18 minus the treatment effect at
  # month 3. We need the full treatment effect at 18,
  # which is trt_factive + trt_factive:visit_f18.
  if (!is.null(fit_cat)) {
    cf <- summary(fit_cat)$tTable
    last_visit_nm <- paste0("trt_factive:visit_f", t_max)
    trt_main <- "trt_factive"
    if (all(c(trt_main, last_visit_nm) %in% rownames(cf))) {
      beta_vec <- cf[, "Value"]
      vcov_mat <- vcov(fit_cat)
      idx <- c(
        which(rownames(cf) == trt_main),
        which(rownames(cf) == last_visit_nm)
      )
      est <- unname(
        beta_vec[trt_main] + beta_vec[last_visit_nm]
      )
      se <- unname(sqrt(
        vcov_mat[idx[1], idx[1]] +
          vcov_mat[idx[2], idx[2]] +
          2 * vcov_mat[idx[1], idx[2]]
      ))
      z <- est / se
      pv <- 2 * pnorm(-abs(z))
      results$cat <- c(
        est = est, se = se, pval = pv, converged = 1
      )
    } else {
      results$cat <- c(
        est = NA, se = NA, pval = NA, converged = 0
      )
    }
  } else {
    results$cat <- c(
      est = NA, se = NA, pval = NA, converged = 0
    )
  }

  results
}

run_simulation <- function(
    n_rep     = 1000,
    kappa_vec = c(0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0),
    seed      = 20260310,
    ...
) {
  set.seed(seed)

  all_results <- vector("list",
    length(kappa_vec) * n_rep
  )
  idx <- 0

  for (kappa in kappa_vec) {
    message("kappa = ", kappa)
    for (rep in seq_len(n_rep)) {
      dat <- simulate_one_trial(kappa = kappa, ...)
      fits <- fit_models(dat)

      idx <- idx + 1
      all_results[[idx]] <- data.frame(
        kappa       = kappa,
        rep         = rep,
        slope_est   = fits$slope["est"],
        slope_se    = fits$slope["se"],
        slope_pval  = fits$slope["pval"],
        slope_conv  = fits$slope["converged"],
        cat_est     = fits$cat["est"],
        cat_se      = fits$cat["se"],
        cat_pval    = fits$cat["pval"],
        cat_conv    = fits$cat["converged"],
        row.names   = NULL
      )
    }
  }

  do.call(rbind, all_results)
}

summarize_results <- function(res, true_delta = -0.45) {
  res |>
    dplyr::group_by(kappa) |>
    dplyr::summarize(
      n_reps         = dplyr::n(),
      slope_conv     = mean(slope_conv, na.rm = TRUE),
      slope_bias     = mean(slope_est, na.rm = TRUE) -
                         true_delta,
      slope_emp_se   = sd(slope_est, na.rm = TRUE),
      slope_mean_se  = mean(slope_se, na.rm = TRUE),
      slope_power    = mean(slope_pval < 0.05,
                            na.rm = TRUE),
      slope_coverage = mean(
        abs(slope_est - true_delta) <
          1.96 * slope_se,
        na.rm = TRUE
      ),
      cat_conv       = mean(cat_conv, na.rm = TRUE),
      cat_bias       = mean(cat_est, na.rm = TRUE) -
                         true_delta,
      cat_emp_se     = sd(cat_est, na.rm = TRUE),
      cat_mean_se    = mean(cat_se, na.rm = TRUE),
      cat_power      = mean(cat_pval < 0.05,
                            na.rm = TRUE),
      cat_coverage   = mean(
        abs(cat_est - true_delta) <
          1.96 * cat_se,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
}
