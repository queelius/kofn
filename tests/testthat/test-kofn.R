# ===========================================================================
# kofn tests after dist.structure adoption (v0.3.0+).
#
# Topology and DGP infrastructure (parallel_system, series_system,
# kofn_system, bridge_system, coherent_system, system_signature,
# system_censoring, system_lifetime, make_dists, f_sys_general,
# S_sys_general, loglik_system, fit_system, rdata_system, exp_dist,
# weibull_dist) lives in dist.structure now and is tested there. kofn
# tests focus on the inference machinery: loglik, score, fit, rdata,
# observation functors, and Scheme 1 / masked / Fisher-info routines.
# ===========================================================================


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
  expect_true(val < 0)
})


test_that("exponential parallel MLE converges", {
  model <- kofn(k = 2, m = 2)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(0.5, 0.3), n = 200)

  fitter <- fit(model)
  result <- fitter(df, n_starts = 3)
  expect_true(result$converged)
  expect_true(all(is.finite(coef(result))))
  expect_true(all(coef(result) > 0))
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

test_that("kofn constructor rejects invalid k/m", {
  expect_error(kofn(k = 0, m = 3))
  expect_error(kofn(k = 4, m = 3))
  expect_error(kofn(k = 1, m = 0))
  expect_error(kofn(k = 1, m = -1))
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


test_that("IE-based density matches dist.structure exp_kofn density", {
  rates <- c(0.5, 0.3, 0.2)
  # kofn parallel = (k_kofn = m): all components must fail. In
  # dist.structure :G convention this is k_dist = m - k_kofn + 1 = 1
  # (system functions if at least 1 component functions).
  sys <- dist.structure::exp_kofn(k = 1L, rates = rates)
  f_dist <- density(sys)

  for (t in c(0.5, 1.0, 2.0)) {
    # IE-based: sum of w_j (kofn parallel fast path)
    ie_density <- sum(vapply(
      seq_along(rates),
      function(j) w_j_exact(t, rates, j),
      numeric(1)
    ))
    expect_equal(ie_density, f_dist(t), tolerance = 1e-10)
  }
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
# Right-censored exponential observations
# ===========================================================================

test_that("exponential loglik handles right-censored observations", {
  model <- kofn(k = 2, m = 2)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(0.5, 0.3), n = 100,
            observe = observe_right_censor(tau = 3.0))

  expect_true(any(df$omega == "right"))
  expect_true(any(df$omega == "exact"))

  ll <- loglik(model)
  val <- ll(df, c(0.5, 0.3))
  expect_true(is.finite(val))

  # The log-likelihood at true params should be approximately maximized.
  val_wrong <- ll(df, c(2.0, 2.0))
  expect_true(val > val_wrong)
})


# ===========================================================================
# Weibull right-censored handling
# ===========================================================================

test_that("Weibull parallel loglik handles right-censored data", {
  model <- kofn(k = 2, m = 2, family = "weibull")

  df <- data.frame(
    t = c(1.0, 2.0, 3.0, 4.0, 5.0),
    omega = c("exact", "exact", "right", "exact", "right"),
    stringsAsFactors = FALSE
  )

  ll <- loglik(model)
  par <- c(1.5, 2.0, 2.0, 3.0)  # shape1, scale1, shape2, scale2

  val <- ll(df, par)
  expect_true(is.finite(val))

  # Right-censored gives a different loglik than treating the same
  # observations as exact.
  df_all_exact <- df
  df_all_exact$omega <- "exact"
  val_exact <- ll(df_all_exact, par)
  expect_true(val != val_exact)
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
  expect_equal(length(coef(result)), 1)
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
  # Identifiable permutation-symmetric quantity is the sum of rates.
  model <- kofn(k = 2, m = 2)
  true_rates <- c(0.5, 0.3)
  true_sum <- sum(true_rates)

  set.seed(123)
  gen <- rdata(model)
  df <- gen(theta = true_rates, n = 500)

  fitter <- fit(model)
  result <- fitter(df, n_starts = 3)

  expect_true(result$converged)
  est_sum <- sum(coef(result))
  expect_equal(est_sum, true_sum, tolerance = 0.2)
})


# ===========================================================================
# parse_params helper
# ===========================================================================

test_that("parse_params extracts shapes and scales correctly", {
  pp_exp <- kofn:::parse_params(c(0.5, 0.3, 0.2), m = 3, family = "exponential")
  expect_equal(pp_exp$shapes, c(1, 1, 1))
  expect_equal(pp_exp$scales, c(2, 10/3, 5))

  pp_wei <- kofn:::parse_params(c(1.5, 2.0, 2.5, 3.0), m = 2, family = "weibull")
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
  expect_true(all(is.finite(coef(result))))
  expect_true(all(coef(result) > 0))
})


# ===========================================================================
# Scheme 1 smoke tests
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
  expect_true(all(is.finite(coef(result))))
  expect_true(all(coef(result) > 0))
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
  expect_equal(obs(10)$omega, "exact")
})


test_that("observe_interval_censor works", {
  obs <- observe_interval_censor(a = 5, b = 10)
  expect_equal(obs(7)$omega, "interval")
  expect_equal(obs(7)$t, 5)
  expect_equal(obs(7)$t_upper, 10)
  expect_equal(obs(3)$omega, "exact")
  expect_equal(obs(3)$t, 3)
  expect_equal(obs(12)$omega, "exact")
  expect_equal(obs(12)$t, 12)
  expect_equal(obs(5)$omega, "interval")
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
  result <- obs(5.0)
  expect_true(result$omega %in% c("exact", "interval"))
})


# ===========================================================================
# Fisher information comparison smoke tests
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
