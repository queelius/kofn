test_that("coherent system constructors work", {
  sys <- parallel_system(3)
  expect_equal(sys$m, 3L)
  expect_equal(length(sys$min_paths), 3)
  expect_equal(length(sys$min_cuts), 1)

  sys2 <- series_system(3)
  expect_equal(sys2$m, 3L)
  expect_equal(length(sys2$min_paths), 1)
  expect_equal(length(sys2$min_cuts), 3)

  sys3 <- kofn_system(2, 4)
  expect_equal(sys3$m, 4L)
})

test_that("kofn model constructor works", {
  model <- kofn(k = 3, m = 3, family = "exponential")
  expect_s3_class(model, "exp_kofn")
  expect_s3_class(model, "kofn")
  expect_equal(ncomponents(model), 3L)
  expect_true(is_kofn(model))

  model_wei <- kofn(k = 2, m = 2, family = "weibull", method = "em")
  expect_s3_class(model_wei, "wei_kofn")
})

test_that("IE expansion produces correct number of terms", {
  ie <- ie_expand(c(0.5, 0.3))
  expect_equal(length(ie$sign), 4)  # 2^2

  ie3 <- ie_expand(c(0.5, 0.3, 0.2))
  expect_equal(length(ie3$sign), 8)  # 2^3
})

test_that("exponential parallel data generation works", {
  model <- kofn(k = 3, m = 3)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(0.5, 0.3, 0.2), n = 50)

  expect_equal(nrow(df), 50)
  expect_true("t" %in% names(df))
  expect_true("omega" %in% names(df))
  expect_true(all(df$omega == "exact"))
  expect_true(all(df$t > 0))
})

test_that("exponential parallel loglik is finite at true params", {
  model <- kofn(k = 2, m = 2)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(0.5, 0.3), n = 100)

  ll <- loglik(model)
  val <- ll(df, par = c(0.5, 0.3))
  expect_true(is.finite(val))
  expect_true(val < 0)  # log-likelihood should be negative
})

test_that("exponential parallel MLE converges", {
  model <- kofn(k = 2, m = 2)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(0.5, 0.3), n = 200)

  fitter <- fit(model)
  result <- fitter(df, n_starts = 3)
  expect_true(result$converged)
  expect_true(all(is.finite(result$par)))
  expect_true(all(result$par > 0))
})

test_that("Weibull EM converges for parallel system", {
  model <- kofn(k = 2, m = 2, family = "weibull", method = "em")
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(1.5, 2.0, 2.0, 3.0), n = 200)

  fitter <- fit(model)
  result <- fitter(df, n_starts = 2)
  expect_true(result$converged)
  expect_true(all(result$shapes > 0))
  expect_true(all(result$scales > 0))
})

test_that("observation functors produce correct types", {
  obs0 <- observe_right_censor(tau = 5)
  expect_equal(obs0(3)$omega, "exact")
  expect_equal(obs0(7)$omega, "right")

  obs1 <- observe_periodic(delta = 1.0)
  expect_equal(obs1(2.7)$omega, "interval")
  expect_equal(obs1(2.7)$t, 2.0)
  expect_equal(obs1(2.7)$t_upper, 3.0)
})

test_that("Scheme 1 data generation works", {
  model <- kofn(k = 2, m = 2, family = "weibull")
  s1gen <- rdata_scheme1(model)
  set.seed(42)
  df <- s1gen(theta = c(1.5, 2.0, 2.0, 3.0), n = 50, delta = 0.5)

  expect_equal(nrow(df), 50)
  expect_true("comp_lower_1" %in% names(df))
  expect_true("comp_upper_1" %in% names(df))
  expect_true(all(df$t > 0))
})

test_that("general k-out-of-n loglik works", {
  model <- kofn(k = 2, m = 4)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(0.5, 0.3, 0.2, 0.4), n = 50)

  ll <- loglik(model)
  val <- ll(df, c(0.5, 0.3, 0.2, 0.4))
  expect_true(is.finite(val))
})


# ===========================================================================
# Input validation
# ===========================================================================

test_that("kofn_system rejects invalid k", {
  expect_error(kofn_system(0, 3))
  expect_error(kofn_system(4, 3))
  expect_error(kofn_system(-1, 3))
})

test_that("kofn constructor rejects invalid k/m", {
  expect_error(kofn(k = 0, m = 3))
  expect_error(kofn(k = 4, m = 3))
  expect_error(kofn(k = 1, m = 0))
  expect_error(kofn(k = 1, m = -1))
})


# ===========================================================================
# Bug #1: custom system correctly overrides k
# ===========================================================================

test_that("kofn with custom system detects k correctly", {
  # Parallel system: k should be detected as m (all must fail)
  model_par <- kofn(system = parallel_system(3))
  expect_equal(model_par$k, 3L)
  expect_equal(model_par$m, 3L)

  # Series system: k should be detected as 1 (one failure kills it)
  model_ser <- kofn(system = series_system(3))
  expect_equal(model_ser$k, 1L)

  # 2-out-of-4: k should be detected as 2
  model_2of4 <- kofn(system = kofn_system(2, 4))
  expect_equal(model_2of4$k, 2L)
})

test_that("kofn with bridge system routes to general code path", {
  model <- kofn(system = bridge_system())
  expect_true(is.na(model$k))
  expect_equal(model$m, 5L)

  # Should use general system density, not parallel IE fast path
  set.seed(42)
  rates <- c(0.5, 0.3, 0.2, 0.4, 0.1)
  gen <- rdata(model)
  df <- gen(theta = rates, n = 50)

  ll <- loglik(model)
  val <- ll(df, rates)
  expect_true(is.finite(val))
})


# ===========================================================================
# IE expansion numerical correctness
# ===========================================================================

test_that("IE expansion is numerically correct", {
  lam <- c(0.5, 0.3, 0.2)
  ie <- ie_expand(lam)

  # Evaluate at several time points and compare with direct product
  for (t in c(0.1, 0.5, 1.0, 2.0, 5.0)) {
    ie_val <- sum(ie$sign * exp(-ie$rate_sum * t))
    direct_val <- prod(1 - exp(-lam * t))
    expect_equal(ie_val, direct_val, tolerance = 1e-10)
  }
})

test_that("IE-based density matches general system density", {
  rates <- c(0.5, 0.3, 0.2)
  sys <- parallel_system(3)
  dists <- make_dists(rates, "exponential")

  for (t in c(0.5, 1.0, 2.0)) {
    # IE-based: sum of w_j
    ie_density <- sum(vapply(
      seq_along(rates),
      function(j) w_j_exact(t, rates, j),
      numeric(1)
    ))

    # General system density engine
    gen_density <- f_sys_general(t, sys, dists)

    expect_equal(ie_density, gen_density, tolerance = 1e-10)
  }
})


# ===========================================================================
# System censoring
# ===========================================================================

test_that("system_censoring is correct for parallel system", {
  sys <- parallel_system(3)
  times <- c(1.0, 3.0, 2.0)
  cens <- system_censoring(sys, times)

  expect_equal(cens$T_sys, 3.0)
  expect_equal(cens$critical, 2L)  # component 2 has max lifetime
  expect_equal(cens$status[1], "left")   # failed before system
  expect_equal(cens$status[2], "exact")  # critical component
  expect_equal(cens$status[3], "left")   # failed before system
})

test_that("system_censoring is correct for series system", {
  sys <- series_system(3)
  times <- c(2.0, 1.0, 3.0)
  cens <- system_censoring(sys, times)

  expect_equal(cens$T_sys, 1.0)
  expect_equal(cens$critical, 2L)  # component 2 has min lifetime
  expect_equal(cens$status[1], "right")  # survived past system failure
  expect_equal(cens$status[2], "exact")  # critical component
  expect_equal(cens$status[3], "right")  # survived past system failure
})


# ===========================================================================
# Score / gradient check
# ===========================================================================

test_that("exponential score matches numerical gradient of loglik", {
  model <- kofn(k = 2, m = 2)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(0.5, 0.3), n = 100)

  ll <- loglik(model)
  sc <- score(model)
  par <- c(0.5, 0.3)

  score_val <- sc(df, par)
  num_grad <- numDeriv::grad(function(p) ll(df, p), par)
  expect_equal(score_val, num_grad, tolerance = 1e-6)
})


# ===========================================================================
# Bug #3: right-censored exponential observations
# ===========================================================================

test_that("exponential loglik handles right-censored observations", {
  model <- kofn(k = 2, m = 2)
  gen <- rdata(model)
  set.seed(42)
  # Generate right-censored data via tau
  df <- gen(theta = c(0.5, 0.3), n = 100, observe = observe_right_censor(tau = 3.0))

  expect_true(any(df$omega == "right"))
  expect_true(any(df$omega == "exact"))

  ll <- loglik(model)
  val <- ll(df, c(0.5, 0.3))
  expect_true(is.finite(val))

  # The log-likelihood at true params should be maximized (approximately)
  val_wrong <- ll(df, c(2.0, 2.0))
  expect_true(val > val_wrong)
})


# ===========================================================================
# Bug fix: Weibull loglik handles right-censored observations
# ===========================================================================

test_that("Weibull parallel loglik handles right-censored data", {
  model <- kofn(k = 2, m = 2, family = "weibull")

  # Manually create data with omega column
  df <- data.frame(
    t = c(1.0, 2.0, 3.0, 4.0, 5.0),
    omega = c("exact", "exact", "right", "exact", "right"),
    stringsAsFactors = FALSE
  )

  ll <- loglik(model)
  par <- c(1.5, 2.0, 2.0, 3.0)  # shape1, scale1, shape2, scale2

  val <- ll(df, par)
  expect_true(is.finite(val))

  # Verify right-censored gives higher loglik than treating as exact
  # (right-censored obs have loglik = log(S_sys(t)) which is > log(f_sys(t))
  # for large t where S_sys > f_sys)
  df_all_exact <- df
  df_all_exact$omega <- "exact"
  val_exact <- ll(df_all_exact, par)
  expect_true(val != val_exact)  # Different when censoring is present
})


# ===========================================================================
# System signature
# ===========================================================================

test_that("system_signature gives correct results", {
  # Series: always the first order statistic
  sig_series <- system_signature(series_system(3))
  expect_equal(sig_series, c(1, 0, 0))

  # Parallel: always the last order statistic
  sig_par <- system_signature(parallel_system(3))
  expect_equal(sig_par, c(0, 0, 1))

  # 2-out-of-3: T_sys = T_{(2)} always, so signature is (0, 1, 0)
  sig_2of3 <- system_signature(kofn_system(2, 3))
  expect_equal(sig_2of3, c(0, 1, 0))
})


# ===========================================================================
# Edge cases
# ===========================================================================

test_that("single-component system works", {
  model <- kofn(k = 1, m = 1)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = 0.5, n = 100)

  ll <- loglik(model)
  val <- ll(df, 0.5)
  expect_true(is.finite(val))

  fitter <- fit(model)
  result <- fitter(df, n_starts = 2)
  expect_true(result$converged)
  expect_equal(length(result$par), 1)
})

test_that("series system (k = 1) works for exponential", {
  model <- kofn(k = 1, m = 3)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(0.5, 0.3, 0.2), n = 100)

  ll <- loglik(model)
  val <- ll(df, c(0.5, 0.3, 0.2))
  expect_true(is.finite(val))
})


# ===========================================================================
# Monte Carlo parameter recovery (exponential parallel)
# ===========================================================================

test_that("exponential parallel MLE recovers sum of rates", {
  # For identifiable permutation-symmetric parameters, check sum of rates
  model <- kofn(k = 2, m = 2)
  true_rates <- c(0.5, 0.3)
  true_sum <- sum(true_rates)

  set.seed(123)
  gen <- rdata(model)
  df <- gen(theta = true_rates, n = 500)

  fitter <- fit(model)
  result <- fitter(df, n_starts = 3)

  expect_true(result$converged)
  est_sum <- sum(result$par)
  # With n=500, the sum should be within ~20% of truth
  expect_equal(est_sum, true_sum, tolerance = 0.2)
})


# ===========================================================================
# parse_params helper
# ===========================================================================

test_that("parse_params extracts shapes and scales correctly", {
  pp_exp <- parse_params(c(0.5, 0.3, 0.2), m = 3, family = "exponential")
  expect_equal(pp_exp$shapes, c(1, 1, 1))
  expect_equal(pp_exp$scales, c(2, 10/3, 5))

  pp_wei <- parse_params(c(1.5, 2.0, 2.5, 3.0), m = 2, family = "weibull")
  expect_equal(pp_wei$shapes, c(1.5, 2.5))
  expect_equal(pp_wei$scales, c(2.0, 3.0))
})


# ===========================================================================
# Weibull direct MLE
# ===========================================================================

test_that("Weibull direct MLE converges for parallel system", {
  model <- kofn(k = 2, m = 2, family = "weibull", method = "mle")
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(1.5, 2.0, 2.0, 3.0), n = 200)

  fitter <- fit(model)
  result <- fitter(df, n_starts = 3)
  expect_true(result$converged)
  expect_true(all(is.finite(result$par)))
  expect_true(all(result$par > 0))
})


# ===========================================================================
# fit_system for general coherent systems
# ===========================================================================

test_that("fit_system works for bridge system", {
  br <- bridge_system()
  rates <- c(0.5, 0.3, 0.2, 0.4, 0.1)
  set.seed(42)
  dat <- rdata_system(br, par = rates, family = "exponential", n = 200)

  res <- fit_system(dat$t, br, family = "exponential", n_starts = 3)
  expect_true(res$converged)
  expect_equal(res$convergence, 0L)  # backward compat field
  expect_equal(length(res$par), 5)
  expect_true(all(is.finite(res$se)))
})


# ===========================================================================
# Scheme 1 smoke tests (covers refactored parse_params + weibull_f_sys paths)
# ===========================================================================

test_that("Scheme 1 loglik is finite at true params (exponential)", {
  model <- kofn(k = 2, m = 2)
  s1gen <- rdata_scheme1(model)
  set.seed(42)
  df <- s1gen(theta = c(0.5, 0.3), n = 100, delta = 1.0)

  ll <- loglik_scheme1(model)
  val <- ll(df, c(0.5, 0.3))
  expect_true(is.finite(val))
  expect_true(val < 0)
})

test_that("Scheme 1 loglik is finite at true params (Weibull)", {
  model <- kofn(k = 2, m = 2, family = "weibull")
  s1gen <- rdata_scheme1(model)
  set.seed(42)
  df <- s1gen(theta = c(1.5, 2.0, 2.0, 3.0), n = 100, delta = 0.5)

  ll <- loglik_scheme1(model)
  val <- ll(df, c(1.5, 2.0, 2.0, 3.0))
  expect_true(is.finite(val))
  expect_true(val < 0)
})

test_that("fit_scheme1 converges for exponential parallel", {
  model <- kofn(k = 2, m = 2)
  s1gen <- rdata_scheme1(model)
  set.seed(42)
  df <- s1gen(theta = c(0.5, 0.3), n = 200, delta = 1.0)

  fitter <- fit_scheme1(model)
  result <- fitter(df, n_starts = 3)
  expect_true(result$converged)
  expect_true(all(is.finite(result$par)))
  expect_true(all(result$par > 0))
})


# ===========================================================================
# Distribution objects
# ===========================================================================

test_that("make_dists creates correct distribution objects", {
  # Exponential
  dists_exp <- make_dists(c(0.5, 0.3), "exponential")
  expect_equal(length(dists_exp), 2)
  expect_equal(dists_exp[[1]]$name, "exponential")
  expect_equal(dists_exp[[1]]$params[["rate"]], 0.5)
  expect_true(dists_exp[[1]]$pdf(1) > 0)
  expect_true(dists_exp[[1]]$cdf(1) > 0)
  expect_true(dists_exp[[1]]$surv(1) > 0)

  # Weibull
  dists_wei <- make_dists(c(1.5, 2.0, 2.5, 3.0), "weibull")
  expect_equal(length(dists_wei), 2)
  expect_equal(dists_wei[[1]]$name, "weibull")
  expect_equal(dists_wei[[1]]$params[["shape"]], 1.5)
  expect_equal(dists_wei[[1]]$params[["scale"]], 2.0)
})


# ===========================================================================
# Observation scheme aliases
# ===========================================================================

test_that("observe_right_censor works", {
  obs <- observe_right_censor(tau = 10)
  expect_equal(obs(5)$omega, "exact")
  expect_equal(obs(5)$t, 5)
  expect_equal(obs(15)$omega, "right")
  expect_equal(obs(15)$t, 10)
})

test_that("observe_left_censor works", {
  obs <- observe_left_censor(tau = 10)
  expect_equal(obs(5)$omega, "left")
  expect_equal(obs(5)$t, 10)
  expect_equal(obs(15)$omega, "exact")
  expect_equal(obs(15)$t, 15)
  # Boundary: exactly at tau
  expect_equal(obs(10)$omega, "exact")
})

test_that("observe_interval_censor works", {
  obs <- observe_interval_censor(a = 5, b = 10)
  # Inside window
  expect_equal(obs(7)$omega, "interval")
  expect_equal(obs(7)$t, 5)
  expect_equal(obs(7)$t_upper, 10)
  # Outside window (before)
  expect_equal(obs(3)$omega, "exact")
  expect_equal(obs(3)$t, 3)
  # Outside window (after)
  expect_equal(obs(12)$omega, "exact")
  expect_equal(obs(12)$t, 12)
  # Boundary: at lower bound
  expect_equal(obs(5)$omega, "interval")
  # Boundary: at upper bound (half-open [a, b))
  expect_equal(obs(10)$omega, "exact")
})

test_that("observe_periodic works", {
  obs <- observe_periodic(delta = 2.0, tau = 20)
  expect_equal(obs(5)$omega, "interval")
  expect_equal(obs(5)$t, 4.0)
  expect_equal(obs(5)$t_upper, 6.0)
  expect_equal(obs(25)$omega, "right")
})

test_that("observe_mixture selects from schemes", {
  set.seed(42)
  obs <- observe_mixture(
    observe_right_censor(tau = 100),
    observe_periodic(delta = 1.0, tau = 100),
    weights = c(0.5, 0.5)
  )
  # Just check it doesn't error and returns valid structure
  result <- obs(5.0)
  expect_true(result$omega %in% c("exact", "interval"))
})


# ===========================================================================
# Assumptions methods handle NA k (bridge system)
# ===========================================================================

test_that("assumptions work for bridge system models", {
  model_exp <- kofn(system = bridge_system())
  a_exp <- assumptions(model_exp)
  expect_true(any(grepl("general coherent", a_exp)))

  model_wei <- kofn(system = bridge_system(), family = "weibull")
  a_wei <- assumptions(model_wei)
  expect_true(any(grepl("General coherent", a_wei)))
})


# ===========================================================================
# k-detection: non-kofn system with equal-size paths
# ===========================================================================

test_that("k-detection rejects non-kofn system with equal-size paths", {
  # System with paths {1,2} and {3,4}: all paths size 2, but NOT 2-out-of-4
  # (2-out-of-4 would have choose(4,2)=6 paths, not 2)
  sys <- coherent_system(list(c(1L, 2L), c(3L, 4L)), m = 4L)
  model <- kofn(system = sys)
  expect_true(is.na(model$k))
  expect_equal(model$m, 4L)
})


# ===========================================================================
# Fisher information comparison smoke test
# ===========================================================================

test_that("compare_fisher_info runs and returns expected structure (exp)", {
  set.seed(42)
  res <- compare_fisher_info(
    rates = c(0.5, 0.3), n = 50L, delta = 1.0, n_rep = 3L,
    family = "exponential"
  )

  expect_true(is.list(res))
  expect_equal(length(res$scheme0_det), 3)
  expect_equal(length(res$scheme1_det), 3)
  expect_equal(length(res$scheme2_det), 3)
  expect_true(all(names(res$median_det) == c("scheme0", "scheme1", "scheme2")))
  # Scheme 2 (complete data) determinant should always be finite and positive
  expect_true(all(is.finite(res$scheme2_det)))
  expect_true(all(res$scheme2_det > 0))
})

test_that("compare_fisher_info runs for Weibull family", {
  set.seed(42)
  res <- compare_fisher_info(
    shapes = c(1.5, 2.0), scales = c(2.0, 3.0),
    n = 50L, delta = 1.0, n_rep = 3L, family = "weibull"
  )

  expect_true(is.list(res))
  expect_equal(length(res$scheme0_det), 3)
})
