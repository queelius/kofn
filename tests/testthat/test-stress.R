# ===========================================================================
# Stress tests for kofn package
# ===========================================================================
# These tests push edge cases, numerical limits, and statistical properties.
# They are more expensive than the unit tests (~2-5 min total).
# ===========================================================================


# ===========================================================================
# Numerical stability: extreme parameters
# ===========================================================================

test_that("exponential loglik handles very small rates", {
  model <- kofn(k = 2, m = 2)
  gen <- rdata(model)
  set.seed(1)
  df <- gen(theta = c(0.001, 0.002), n = 50)

  ll <- loglik(model)
  val <- ll(df, c(0.001, 0.002))
  expect_true(is.finite(val))

  # Should be better than a wildly wrong guess
  val_wrong <- ll(df, c(10, 10))
  expect_true(val > val_wrong)
})

test_that("exponential loglik handles very large rates", {
  model <- kofn(k = 2, m = 2)
  gen <- rdata(model)
  set.seed(1)
  df <- gen(theta = c(50, 100), n = 50)

  ll <- loglik(model)
  val <- ll(df, c(50, 100))
  expect_true(is.finite(val))
})

test_that("Weibull loglik handles extreme shape parameters", {
  model <- kofn(k = 2, m = 2, family = "weibull")
  gen <- rdata(model)

  # Very small shape (heavy-tailed, near-exponential)
  set.seed(1)
  df_small <- gen(theta = c(0.5, 2.0, 0.5, 3.0), n = 50)
  ll <- loglik(model)
  val <- ll(df_small, c(0.5, 2.0, 0.5, 3.0))
  expect_true(is.finite(val))

  # Large shape (concentrated, near-deterministic)
  set.seed(1)
  df_large <- gen(theta = c(5.0, 2.0, 5.0, 3.0), n = 50)
  val2 <- ll(df_large, c(5.0, 2.0, 5.0, 3.0))
  expect_true(is.finite(val2))
})

test_that("loglik returns -Inf for zero or negative parameters", {
  model_exp <- kofn(k = 2, m = 2)
  gen <- rdata(model_exp)
  set.seed(1)
  df <- gen(theta = c(0.5, 0.3), n = 20)

  ll <- loglik(model_exp)
  expect_equal(ll(df, c(0, 0.3)), -Inf)
  expect_equal(ll(df, c(-1, 0.3)), -Inf)
  expect_equal(ll(df, c(0.5, -0.1)), -Inf)

  model_wei <- kofn(k = 2, m = 2, family = "weibull")
  ll_w <- loglik(model_wei)
  df_w <- data.frame(t = c(1, 2, 3))
  expect_equal(ll_w(df_w, c(0, 2, 1, 3)), -Inf)
  expect_equal(ll_w(df_w, c(1, 0, 1, 3)), -Inf)
})


# ===========================================================================
# Scaling: component count
# ===========================================================================

test_that("IE expansion works for m=1 (degenerate parallel)", {
  ie <- ie_expand(0.5)
  expect_equal(length(ie$sign), 2)  # 2^1

  # Should equal 1 - exp(-0.5*t) at t=1
  val <- sum(ie$sign * exp(-ie$rate_sum * 1))
  expect_equal(val, 1 - exp(-0.5), tolerance = 1e-12)
})

test_that("IE expansion works for m=8 (256 terms)", {
  rates <- seq(0.1, 0.8, by = 0.1)
  ie <- ie_expand(rates)
  expect_equal(length(ie$sign), 256)  # 2^8

  # Verify against direct product at t=1
  ie_val <- sum(ie$sign * exp(-ie$rate_sum * 1))
  direct_val <- prod(1 - exp(-rates * 1))
  expect_equal(ie_val, direct_val, tolerance = 1e-8)
})

test_that("exponential parallel MLE works for m=5", {
  model <- kofn(k = 5, m = 5)
  true_rates <- c(0.5, 0.3, 0.2, 0.4, 0.1)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = true_rates, n = 300)

  fitter <- fit(model)
  result <- fitter(df, n_starts = 3)
  expect_true(result$converged)
  expect_true(all(result$par > 0))

  # Sum of rates should be approximately recovered
  expect_equal(sum(result$par), sum(true_rates), tolerance = 0.3)
})

test_that("general system density works for 3-out-of-5", {
  model <- kofn(k = 3, m = 5)
  true_rates <- c(0.5, 0.3, 0.2, 0.4, 0.1)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = true_rates, n = 100)

  ll <- loglik(model)
  val <- ll(df, true_rates)
  expect_true(is.finite(val))

  # Loglik at true params should beat random params
  val_wrong <- ll(df, c(5, 5, 5, 5, 5))
  expect_true(val > val_wrong)
})


# ===========================================================================
# Cross-validation: IE vs general system density
# ===========================================================================

test_that("IE loglik matches general system density loglik for parallel", {
  # Both compute the marginal system density — they should agree exactly
  rates <- c(0.5, 0.3, 0.2)
  sys <- parallel_system(3)

  model <- kofn(k = 3, m = 3)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = rates, n = 50)

  ll_ie <- loglik(model)
  val_ie <- ll_ie(df, rates)

  val_gen <- loglik_system(df$t, sys, rates, family = "exponential")

  expect_equal(val_ie, val_gen, tolerance = 1e-8)
})


# ===========================================================================
# Cross-validation: Weibull EM vs direct MLE
# ===========================================================================

test_that("Weibull EM and direct MLE give similar estimates", {
  true_par <- c(1.5, 2.0, 2.0, 3.0)

  model_em <- kofn(k = 2, m = 2, family = "weibull", method = "em")
  model_mle <- kofn(k = 2, m = 2, family = "weibull", method = "mle")
  gen <- rdata(model_em)

  set.seed(42)
  df <- gen(theta = true_par, n = 300)

  fit_em <- fit(model_em)
  fit_mle <- fit(model_mle)

  res_em <- fit_em(df, n_starts = 3)
  res_mle <- fit_mle(df, n_starts = 3)

  # Both should converge
	expect_true(res_em$converged)
  expect_true(res_mle$converged)

  # Log-likelihoods should be close (both found the MLE)
  ll <- loglik(model_em)
  ll_em <- ll(df, res_em$par)
  ll_mle <- ll(df, res_mle$par)
  expect_equal(ll_em, ll_mle, tolerance = 0.5)
})


# ===========================================================================
# Right-censoring stress
# ===========================================================================

test_that("exponential MLE converges with heavy right-censoring", {
  model <- kofn(k = 2, m = 2)
  gen <- rdata(model)
  set.seed(42)
  # tau=1.0 with rates 0.5, 0.3 means most systems survive past tau
  df <- gen(theta = c(0.5, 0.3), n = 200, observe = observe_right_censor(tau = 1.0))

  pct_censored <- mean(df$omega == "right")
  expect_true(pct_censored > 0.3)  # should have substantial censoring

  fitter <- fit(model)
  result <- fitter(df, n_starts = 5)
  expect_true(result$converged)
  expect_true(all(result$par > 0))
})

test_that("Weibull loglik handles heavy right-censoring", {
  model <- kofn(k = 2, m = 2, family = "weibull")
  true_par <- c(1.5, 2.0, 2.0, 3.0)

  # Manually create data with ~50% right-censoring
  set.seed(42)
  gen <- rdata(model)
  df_full <- gen(theta = true_par, n = 100)

  # Apply right-censoring at median
  tau <- median(df_full$t)
  df <- df_full
  df$omega <- ifelse(df$t <= tau, "exact", "right")
  df$t[df$omega == "right"] <- tau

  ll <- loglik(model)
  val <- ll(df, true_par)
  expect_true(is.finite(val))
})


# ===========================================================================
# Monte Carlo: parameter recovery with multiple replicates
# ===========================================================================

test_that("exponential parallel MLE is approximately unbiased (sum of rates)", {
  model <- kofn(k = 2, m = 2)
  true_rates <- c(0.5, 0.3)
  true_sum <- sum(true_rates)
  gen <- rdata(model)
  fitter <- fit(model)

  n_rep <- 20
  sum_estimates <- numeric(n_rep)

  set.seed(123)
  for (r in seq_len(n_rep)) {
    df <- gen(theta = true_rates, n = 200)
    result <- fitter(df, n_starts = 2)
    if (result$converged) {
      sum_estimates[r] <- sum(result$par)
    } else {
      sum_estimates[r] <- NA
    }
  }

  valid <- !is.na(sum_estimates)
  expect_true(sum(valid) >= 15)  # at least 75% should converge

  # Mean of sum estimates should be within 15% of truth
  mean_est <- mean(sum_estimates[valid])
  expect_equal(mean_est, true_sum, tolerance = 0.15 * true_sum)
})


# ===========================================================================
# Scheme 1: periodic inspection stress
# ===========================================================================

test_that("Scheme 1 fit recovers parameters (exponential)", {
  model <- kofn(k = 2, m = 2)
  true_rates <- c(0.5, 0.3)
  s1gen <- rdata_scheme1(model)

  set.seed(42)
  df <- s1gen(theta = true_rates, n = 500, delta = 0.5)

  fitter <- fit_scheme1(model)
  result <- fitter(df, n_starts = 3)
  expect_true(result$converged)

  # Sum of rates should be approximately right
  expect_equal(sum(result$par), sum(true_rates), tolerance = 0.3)
})

test_that("Scheme 1 loglik is higher at true params than wrong params", {
  model <- kofn(k = 2, m = 2, family = "weibull")
  true_par <- c(1.5, 2.0, 2.0, 3.0)
  s1gen <- rdata_scheme1(model)

  set.seed(42)
  df <- s1gen(theta = true_par, n = 200, delta = 0.5)

  ll <- loglik_scheme1(model)
  val_true <- ll(df, true_par)
  val_wrong <- ll(df, c(0.5, 5.0, 0.5, 5.0))

  expect_true(is.finite(val_true))
  expect_true(is.finite(val_wrong))
  expect_true(val_true > val_wrong)
})


# ===========================================================================
# Coherent system: bridge and consecutive-k
# ===========================================================================

test_that("bridge system estimation recovers sorted rates", {
  br <- bridge_system()
  true_rates <- c(0.5, 0.4, 0.3, 0.2, 0.1)

  set.seed(42)
  dat <- rdata_system(br, par = true_rates, family = "exponential", n = 500)

  res <- fit_system(dat$t, br, family = "exponential", n_starts = 5)
  expect_true(res$converged)

  # Sorted estimates should roughly match sorted truth
  expect_equal(sort(res$par), sort(true_rates), tolerance = 0.3)
})

test_that("consecutive-k system works end-to-end", {
  sys <- consecutive_k_system(2, 4)
  expect_equal(sys$m, 4L)

  # Generate data and evaluate loglik
  rates <- c(0.5, 0.3, 0.2, 0.4)
  set.seed(42)
  dat <- rdata_system(sys, par = rates, family = "exponential", n = 100)
  expect_true(all(dat$t > 0))

  ll_val <- loglik_system(dat$t, sys, rates, family = "exponential")
  expect_true(is.finite(ll_val))

  # Should be better than wildly wrong rates
  ll_wrong <- loglik_system(dat$t, sys, c(10, 10, 10, 10),
                            family = "exponential")
  expect_true(ll_val > ll_wrong)
})


# ===========================================================================
# System signatures for non-trivial systems
# ===========================================================================

test_that("bridge system signature sums to 1", {
  sig <- system_signature(bridge_system())
  expect_equal(sum(sig), 1.0)
  expect_equal(length(sig), 5)
  expect_true(all(sig >= 0))
})

test_that("consecutive-k-out-of-n signature sums to 1", {
  sig <- system_signature(consecutive_k_system(2, 4))
  expect_equal(sum(sig), 1.0)
  expect_true(all(sig >= 0))
})


# ===========================================================================
# fisher_mle interface: base R stats generics
# ===========================================================================

test_that("fisher_mle result supports base R stats generics", {
  model <- kofn(k = 2, m = 2)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(0.5, 0.3), n = 200)

  fitter <- fit(model)
  result <- fitter(df, n_starts = 3)

  # coef
  expect_equal(length(coef(result)), 2)
  expect_true(all(coef(result) > 0))

  # vcov
  V <- vcov(result)
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(2, 2))
  # vcov should be positive definite (or at least have positive diagonal)
  if (all(is.finite(V))) {
    expect_true(all(diag(V) > 0))
  }

  # logLik
  ll <- logLik(result)
  expect_true(is.finite(ll))
  expect_true(ll < 0)

  # AIC / BIC
  aic <- AIC(result)
  expect_true(is.finite(aic))

  # confint
  ci <- confint(result)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 2)
  # Lower bound should be less than upper bound
  if (all(is.finite(ci))) {
    expect_true(all(ci[, 1] < ci[, 2]))
  }
})


# ===========================================================================
# Degenerate convergence guard
# ===========================================================================

test_that("multistart_mle detects degenerate convergence", {
  # Create a neg_ll that always returns the penalty value
  bad_neg_ll <- function(par) .Machine$double.xmax / 2

  # Expect warning about non-invertible Hessian from fisher_mle
  result <- suppressWarnings(
    multistart_mle(bad_neg_ll, par0 = c(1, 1), n_par = 2L,
                   n_starts = 2L, nobs = 10L)
  )
  expect_false(result$converged)
})


# ===========================================================================
# Observation mixture stress
# ===========================================================================

test_that("mixed observation scheme data generates valid observations", {
  set.seed(42)
  obs <- observe_mixture(
    observe_right_censor(tau = 5),
    observe_periodic(delta = 1.0, tau = 5),
    weights = c(0.5, 0.5)
  )

  # Apply to 100 random times
  results <- lapply(runif(100, 0.1, 10), obs)
  omegas <- vapply(results, `[[`, character(1), "omega")

  expect_true(all(omegas %in% c("exact", "right", "interval")))
  # With tau=5 and times up to 10, should have some right-censored
  expect_true(any(omegas == "right"))
  # With scheme1, should have some interval-censored
  expect_true(any(omegas == "interval"))
})


# ===========================================================================
# Weibull system density: edge cases
# ===========================================================================

test_that("weibull_f_sys is non-negative and integrates to ~1", {
  shapes <- c(1.5, 2.0)
  scales <- c(2.0, 3.0)

  # Density should be non-negative
  t_grid <- seq(0.01, 15, length.out = 200)
  f_vals <- vapply(t_grid, function(t) weibull_f_sys(t, shapes, scales),
                   numeric(1))
  expect_true(all(f_vals >= 0))

  # Trapezoidal integration should approximate 1
  dt <- diff(t_grid)
  integral <- sum((f_vals[-1] + f_vals[-length(f_vals)]) / 2 * dt)
  expect_equal(integral, 1.0, tolerance = 0.05)
})

test_that("weibull_S_sys is monotone decreasing", {
  shapes <- c(1.5, 2.0)
  scales <- c(2.0, 3.0)

  t_grid <- seq(0.1, 10, by = 0.5)
  s_vals <- vapply(t_grid, function(t) weibull_S_sys(t, shapes, scales),
                   numeric(1))

  # Should be monotone decreasing
  expect_true(all(diff(s_vals) <= 0))
  # Should start near 1 and approach 0
  expect_true(s_vals[1] > 0.5)
  expect_true(s_vals[length(s_vals)] < 0.1)
})


# ===========================================================================
# IE expansion: mathematical identities
# ===========================================================================

test_that("w_j functions sum to f_sys for exponential parallel", {
  rates <- c(0.5, 0.3, 0.2, 0.4)
  m <- length(rates)

  for (t in c(0.5, 1.0, 2.0, 5.0)) {
    # sum_j w_j(t) = f_sys(t)
    w_sum <- sum(vapply(seq_len(m), function(j) w_j_exact(t, rates, j),
                        numeric(1)))

    # f_sys via IE expansion of the full system CDF derivative
    # d/dt F_sys(t) = d/dt prod(1 - exp(-lam_i * t))
    ie <- ie_expand(rates)
    f_sys <- sum(ie$rate_sum * ie$sign * exp(-ie$rate_sum * t))
    # Note: d/dt sum(s_k * exp(-r_k * t)) = sum(-r_k * s_k * exp(-r_k * t))
    # But F_sys = prod(1-exp(-lam*t)), so f_sys = -d/dt S_sys
    # Let's use numerical derivative instead
    eps <- 1e-7
    F_t <- F_sys_exp(t, rates)
    F_t_eps <- F_sys_exp(t + eps, rates)
    f_sys_num <- (F_t_eps - F_t) / eps

    expect_equal(w_sum, f_sys_num, tolerance = 1e-5)
  }
})

test_that("w_j_integral over [0, Inf) equals 1/m for equal rates", {
  # For equal rates, by symmetry each component is equally likely
  # to be critical, so integral of w_j over [0, inf) = 1/m
  m <- 3
  rate <- 0.5
  rates <- rep(rate, m)

  for (j in seq_len(m)) {
    # Integrate from 0 to a large value (effectively infinity)
    val <- w_j_integral(0, 100, rates, j)
    expect_equal(val, 1 / m, tolerance = 1e-4)
  }
})
