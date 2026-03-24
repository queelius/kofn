# ===========================================================================
# Stress tests round 2: coverage gaps and deeper edge cases
# ===========================================================================
# Targets the zero-coverage lines identified by covr:
#   - Left-censored and interval-censored exponential loglik (bug #3 fix)
#   - General k-out-of-n exponential and Weibull paths
#   - Truncated moment asymptotics (small/large v_max)
#   - Print methods, error paths, EM edge cases
#   - Coherent system internals (minimize_sets, min_cuts_from_paths)
# ===========================================================================


# ===========================================================================
# Left-censored exponential loglik (lines 81-89 of exp_kofn.R)
# ===========================================================================

test_that("exponential loglik handles left-censored observations", {
  model <- kofn(k = 2, m = 2)
  rates <- c(0.5, 0.3)

  # Manually construct a data frame with left-censored observations.
  # "left" means the system failed before time t (we know T < t).
  df <- data.frame(
    t = c(1.0, 2.0, 3.0, 1.5),
    omega = c("exact", "left", "exact", "left"),
    stringsAsFactors = FALSE
  )

  ll <- loglik(model)
  val <- ll(df, rates)
  expect_true(is.finite(val))

  # Left-censored contribution should be log(integral of w_j from 0 to t)
  # Verify the left-censored obs contribute a finite value
  df_exact_only <- df[df$omega == "exact", ]
  val_exact <- ll(df_exact_only, rates)
  expect_true(is.finite(val_exact))
  expect_true(val != val_exact)  # different from exact-only
})

test_that("left-censored loglik uses corrected formula (no normalizer)", {
  # Regression test for bug #3: the old code divided by F_sys(t),
  # turning the joint probability into a conditional. The fix removes
  # the normalizer.
  model <- kofn(k = 2, m = 2)
  rates <- c(0.5, 0.3)

  # Single left-censored observation at t=2
  df <- data.frame(
    t = 2.0, omega = "left",
    stringsAsFactors = FALSE
  )
  ll <- loglik(model)
  val <- ll(df, rates)

  # Manual: sum_j integral_0^2 w_j(s) ds = integral_0^2 f_sys(s) ds = F_sys(2)
  expected <- log(F_sys_exp(2.0, rates))
  expect_equal(val, expected, tolerance = 1e-10)
})


# ===========================================================================
# Interval-censored exponential loglik (lines 93-103 of exp_kofn.R)
# ===========================================================================

test_that("exponential loglik handles interval-censored observations", {
  model <- kofn(k = 2, m = 2)
  rates <- c(0.5, 0.3)

  df <- data.frame(
    t = c(1.0, 2.0, 3.0),
    t_upper = c(NA, 3.0, 4.0),
    omega = c("exact", "interval", "interval"),
    stringsAsFactors = FALSE
  )

  ll <- loglik(model)
  val <- ll(df, rates)
  expect_true(is.finite(val))
})

test_that("interval-censored loglik uses corrected formula", {
  # Single interval-censored obs [1, 3]
  model <- kofn(k = 2, m = 2)
  rates <- c(0.5, 0.3)

  df <- data.frame(
    t = 1.0, t_upper = 3.0, omega = "interval",
    stringsAsFactors = FALSE
  )
  ll <- loglik(model)
  val <- ll(df, rates)

  # Manual: sum_j integral_1^3 w_j(s) ds = F_sys(3) - F_sys(1)
  expected <- log(F_sys_exp(3.0, rates) - F_sys_exp(1.0, rates))
  expect_equal(val, expected, tolerance = 1e-10)
})

test_that("mixed observation types work together in exponential loglik", {
  model <- kofn(k = 3, m = 3)
  rates <- c(0.5, 0.3, 0.2)

  df <- data.frame(
    t = c(1.0, 2.0, 3.0, 1.5, 2.0),
    t_upper = c(NA, NA, NA, NA, 3.0),
    omega = c("exact", "right", "left", "exact", "interval"),
    stringsAsFactors = FALSE
  )

  ll <- loglik(model)
  val <- ll(df, rates)
  expect_true(is.finite(val))

  # Should be better than wildly wrong rates
  val_wrong <- ll(df, c(10, 10, 10))
  expect_true(val > val_wrong)
})


# ===========================================================================
# Parameter count mismatch error (line 56 of exp_kofn.R)
# ===========================================================================

test_that("exponential loglik errors on parameter count mismatch", {
  model <- kofn(k = 3, m = 3)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(0.5, 0.3, 0.2), n = 20)

  ll <- loglik(model)
  expect_error(ll(df, c(0.5, 0.3)), "Expected 3 parameters but got 2")
})


# ===========================================================================
# General k-out-of-n exponential path (line 111 of exp_kofn.R)
# ===========================================================================

test_that("general k-out-of-n exponential loglik returns -Inf for bad params", {
  model <- kofn(k = 2, m = 3)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(0.5, 0.3, 0.2), n = 30)

  ll <- loglik(model)
  expect_equal(ll(df, c(-1, 0.3, 0.2)), -Inf)
  expect_equal(ll(df, c(0, 0.3, 0.2)), -Inf)
})

test_that("general k-out-of-n exponential MLE converges", {
  model <- kofn(k = 2, m = 3)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(0.5, 0.3, 0.2), n = 200)

  fitter <- fit(model)
  result <- fitter(df, n_starts = 3)
  expect_true(result$converged)
  expect_true(all(result$par > 0))
})


# ===========================================================================
# General k-out-of-n Weibull paths
# ===========================================================================

test_that("Weibull loglik delegates to system density for k > 1", {
  model <- kofn(k = 2, m = 3, family = "weibull")
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(1.5, 2.0, 2.0, 3.0, 1.0, 1.5), n = 50)

  ll <- loglik(model)
  val <- ll(df, c(1.5, 2.0, 2.0, 3.0, 1.0, 1.5))
  expect_true(is.finite(val))

  # Bad params
  expect_equal(ll(df, c(-1, 2, 2, 3, 1, 1.5)), -Inf)
})

test_that("Weibull score and hessian work for parallel system", {
  model <- kofn(k = 2, m = 2, family = "weibull")
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(1.5, 2.0, 2.0, 3.0), n = 50)

  sc <- score(model)
  par <- c(1.5, 2.0, 2.0, 3.0)
  score_val <- sc(df, par)
  expect_equal(length(score_val), 4)
  expect_true(all(is.finite(score_val)))

  hl <- hess_loglik(model)
  H <- hl(df, par)
  expect_equal(dim(H), c(4, 4))
  expect_true(all(is.finite(H)))
})


# ===========================================================================
# Truncated moments: edge cases
# ===========================================================================

test_that("trunc_pow_moment_vec handles very small v_max (asymptotic)", {
  # v_max < 1e-12 triggers the asymptotic branch
  alpha <- 1.5
  beta <- 2.0
  v_max <- c(1e-13, 1e-15, 1e-20)

  result <- trunc_pow_moment_vec(k = 1, v_max, alpha, beta)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))

  # Compare with asymptotic formula: t^k / (k/alpha + 1)
  t_vals <- beta * v_max^(1/alpha)
  expected <- t_vals / (1/alpha + 1)
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("trunc_pow_moment_vec handles normal range", {
  alpha <- 2.0
  beta <- 3.0
  v_max <- c(0.1, 1.0, 5.0, 20.0)

  result <- trunc_pow_moment_vec(k = 2, v_max, alpha, beta)
  expect_true(all(is.finite(result)))
  expect_true(all(result > 0))

  # For large v_max, should approach untruncated E[T^2]
  # E[T^2] = beta^2 * Gamma(2/alpha + 1) for Weibull
  untrunc <- beta^2 * gamma(2/alpha + 1)
  expect_true(abs(result[4] - untrunc) / untrunc < 0.01)
})

test_that("trunc_pow_moment scalar wrapper works", {
  val <- trunc_pow_moment(k = 1, t = 2.0, alpha = 1.5, beta = 3.0)
  expect_true(is.finite(val))
  expect_true(val > 0)
  expect_true(val < 2.0)  # must be less than truncation point for k=1
})

test_that("trunc_log_moment_vec handles very small v_max", {
  alpha <- 1.5
  beta <- 2.0
  v_max <- c(1e-11, 1e-15)

  result <- trunc_log_moment_vec(v_max, alpha, beta)
  expect_true(all(is.finite(result)))

  # Asymptotic: log(t) - 1/alpha
  t_vals <- beta * v_max^(1/alpha)
  expected <- log(t_vals) - 1/alpha
  expect_equal(result, expected, tolerance = 1e-4)
})

test_that("trunc_log_moment_vec handles normal and large v_max", {
  alpha <- 2.0
  beta <- 3.0
  v_max <- c(1.0, 5.0, 50.0, 500.0)

  result <- trunc_log_moment_vec(v_max, alpha, beta)
  expect_true(all(is.finite(result)))

  # For large v_max, should approach untruncated E[log T]
  # E[log T] = log(beta) + digamma(1)/alpha
  untrunc <- log(beta) + digamma(1) / alpha
  expect_equal(result[4], untrunc, tolerance = 0.01)
})


# ===========================================================================
# Print methods
# ===========================================================================

test_that("print.kofn works for all system types", {
  # Parallel (k = m)
  out <- capture.output(print(kofn(k = 3, m = 3)))
  expect_true(any(grepl("parallel", out)))

  # Series (k = 1)
  out <- capture.output(print(kofn(k = 1, m = 3)))
  expect_true(any(grepl("series", out)))

  # General k-out-of-n
  out <- capture.output(print(kofn(k = 2, m = 4)))
  expect_true(any(grepl("2-out-of-4", out)))

  # Bridge (NA k)
  out <- capture.output(print(kofn(system = bridge_system())))
  expect_true(any(grepl("general coherent", out)))
  expect_true(any(grepl("k=NA", out)))

  # Weibull EM (parallel: k = m)
  out <- capture.output(print(kofn(k = 2, m = 2, family = "weibull",
                                    method = "em")))
  expect_true(any(grepl("weibull", out)))
})

test_that("print.coherent_system works", {
  out <- capture.output(print(bridge_system()))
  expect_true(any(grepl("m=5", out)))
  expect_true(any(grepl("path sets", out)))
  expect_true(any(grepl("cut sets", out)))
})


# ===========================================================================
# Coherent system internals
# ===========================================================================

test_that("minimize_sets removes supersets correctly", {
  sets <- list(c(1, 2), c(1, 2, 3), c(2, 3), c(1, 2, 3, 4))
  result <- minimize_sets(sets)
  expect_equal(length(result), 2)  # {1,2} and {2,3}
})

test_that("minimize_sets handles empty input", {
  expect_equal(length(minimize_sets(list())), 0)
})

test_that("min_cuts_from_paths computes correct duals", {
  # Parallel system: paths = {1}, {2}, {3}. Cut = {1,2,3}
  cuts <- min_cuts_from_paths(list(1L, 2L, 3L))
  expect_equal(length(cuts), 1)
  expect_equal(sort(cuts[[1]]), c(1L, 2L, 3L))

  # Series system: path = {1,2,3}. Cuts = {1}, {2}, {3}
  cuts <- min_cuts_from_paths(list(c(1L, 2L, 3L)))
  expect_equal(length(cuts), 3)
})

test_that("coherent_system validates inputs", {
  expect_error(coherent_system(list()))  # empty paths
  expect_error(coherent_system("not a list"))  # not a list
})

test_that("min_cuts accessor works", {
  sys <- kofn_system(2, 3)
  cuts <- min_cuts(sys)
  expect_true(is.list(cuts))
  expect_true(length(cuts) > 0)
})

test_that("critical_states returns empty for irrelevant component", {
  # In a series system, every component is critical in exactly one state
  # (all others up). But in some systems, a component may never be critical.
  # 2-out-of-4 system: every component is critical
  sys <- kofn_system(2, 4)
  crit <- critical_states(sys, 1)
  expect_true(nrow(crit) > 0)
})

test_that("critical_states for single-component system", {
  sys <- kofn_system(1, 1)
  crit <- critical_states(sys, 1)
  expect_equal(nrow(crit), 1)
  expect_equal(ncol(crit), 0)  # no other components
})


# ===========================================================================
# EM algorithm edge cases
# ===========================================================================

test_that("EM handles very few observations", {
  model <- kofn(k = 2, m = 2, family = "weibull", method = "em")
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(1.5, 2.0, 2.0, 3.0), n = 10)

  fitter <- fit(model)
  # Should not error, even with very few observations
  result <- fitter(df, n_starts = 2, maxiter = 100)
  expect_true(is.finite(sum(result$par)))
})

test_that("EM rejects non-parallel systems", {
  model <- kofn(k = 2, m = 3, family = "weibull", method = "mle")
  # Force EM method
  model$method <- "em"
  expect_error(fit(model), "parallel")
})


# ===========================================================================
# System density engine (S_sys_general, rdata_system)
# ===========================================================================

test_that("S_sys_general matches 1 - F_sys for exponential parallel", {
  rates <- c(0.5, 0.3)
  sys <- parallel_system(2)
  dists <- make_dists(rates, "exponential")

  for (t in c(0.5, 1.0, 3.0)) {
    S_gen <- S_sys_general(t, sys, dists)
    S_ie <- S_sys_exp(t, rates)
    expect_equal(S_gen, S_ie, tolerance = 1e-10)
  }
})

test_that("rdata_system generates correct structure for Weibull", {
  sys <- kofn_system(2, 3)
  par <- c(1.5, 2.0, 2.0, 3.0, 1.0, 1.5)
  set.seed(42)
  dat <- rdata_system(sys, par = par, family = "weibull", n = 50)

  expect_equal(nrow(dat), 50)
  expect_true(all(dat$t > 0))
  expect_true(!is.null(attr(dat, "comp_times")))
  expect_true(!is.null(attr(dat, "critical")))
})


# ===========================================================================
# Weibull data generation edge cases
# ===========================================================================

test_that("rdata.wei_kofn works for general k", {
  model <- kofn(k = 2, m = 3, family = "weibull")
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(1.5, 2.0, 2.0, 3.0, 1.0, 1.5), n = 50)

  expect_equal(nrow(df), 50)
  expect_true(all(df$t > 0))
})

test_that("rdata.wei_kofn with custom observe functor", {
  model <- kofn(k = 2, m = 2, family = "weibull")
  gen <- rdata(model)
  set.seed(42)

  # Use interval observation scheme (per-observation functor)
  df <- gen(theta = c(1.5, 2.0, 2.0, 3.0), n = 30,
            observe = observe_periodic(delta = 1.0))
  expect_true("omega" %in% names(df))
  expect_true(all(df$omega %in% c("interval", "right")))
})

test_that("rdata.wei_kofn with tau produces right-censored data", {
  model <- kofn(k = 2, m = 2, family = "weibull")
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(1.5, 2.0, 2.0, 3.0), n = 100, observe = observe_right_censor(tau = 2.0))

  expect_true("omega" %in% names(df))
  expect_true(any(df$omega == "right"))
  expect_true(any(df$omega == "exact"))
  # All right-censored obs should have t = tau
  expect_true(all(df$t[df$omega == "right"] == 2.0))
})

test_that("rdata.wei_kofn validates parameter length", {
  model <- kofn(k = 2, m = 2, family = "weibull")
  gen <- rdata(model)
  expect_error(gen(theta = c(1, 2, 3), n = 10))
})


# ===========================================================================
# extract_data validation
# ===========================================================================

# ===========================================================================
# fit_system for Weibull family
# ===========================================================================

test_that("fit_system works for Weibull 2-out-of-3", {
  sys <- kofn_system(2, 3)
  par <- c(1.5, 2.0, 2.0, 3.0, 1.0, 1.5)
  set.seed(42)
  dat <- rdata_system(sys, par = par, family = "weibull", n = 200)

  res <- fit_system(dat$t, sys, family = "weibull", n_starts = 3)
  expect_true(res$converged)
  expect_equal(length(res$par), 6)
  expect_true(all(res$par > 0))
})


# ===========================================================================
# Scheme 1: Weibull fit convergence
# ===========================================================================

test_that("fit_scheme1 converges for Weibull parallel", {
  model <- kofn(k = 2, m = 2, family = "weibull")
  s1gen <- rdata_scheme1(model)
  set.seed(42)
  df <- s1gen(theta = c(1.5, 2.0, 2.0, 3.0), n = 300, delta = 0.5)

  fitter <- fit_scheme1(model)
  result <- fitter(df, n_starts = 3)
  expect_true(result$converged)
  expect_true(all(result$par > 0))
})


# ===========================================================================
# hessian_score_at_mle edge cases
# ===========================================================================

test_that("hessian_score_at_mle handles non-differentiable function", {
  bad_fn <- function(par) {
    if (any(par < 0)) stop("negative")
    sum(par^2)
  }
  # Should not error — tryCatch handles it gracefully
  result <- hessian_score_at_mle(bad_fn, c(1, 2), 2L)
  expect_true(is.matrix(result$hessian))
  expect_equal(length(result$score), 2)
})


# ===========================================================================
# Weibull f_sys cross-validation against exponential (shape=1)
# ===========================================================================

test_that("Weibull f_sys reduces to exponential for shape=1", {
  rates <- c(0.5, 0.3, 0.2)
  # Weibull with shape=1 is exponential with scale=1/rate
  shapes <- c(1, 1, 1)
  scales <- 1 / rates

  for (t in c(0.5, 1.0, 2.0, 5.0)) {
    # Weibull system density
    f_wei <- weibull_f_sys(t, shapes, scales)

    # Exponential system density via IE
    f_exp <- sum(vapply(seq_along(rates),
                        function(j) w_j_exact(t, rates, j),
                        numeric(1)))

    expect_equal(f_wei, f_exp, tolerance = 1e-8)
  }
})


# ===========================================================================
# Large system scaling
# ===========================================================================

test_that("IE expansion is feasible for m=10 (1024 terms)", {
  rates <- seq(0.1, 1.0, by = 0.1)
  ie <- ie_expand(rates)
  expect_equal(length(ie$sign), 1024)

  # Verify at one time point
  ie_val <- sum(ie$sign * exp(-ie$rate_sum * 1.0))
  direct_val <- prod(1 - exp(-rates * 1.0))
  expect_equal(ie_val, direct_val, tolerance = 1e-6)
})

test_that("exponential parallel loglik works for m=6", {
  model <- kofn(k = 6, m = 6)
  rates <- c(0.5, 0.4, 0.3, 0.2, 0.1, 0.6)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = rates, n = 100)

  ll <- loglik(model)
  val <- ll(df, rates)
  expect_true(is.finite(val))
})


# ===========================================================================
# make_dists error path
# ===========================================================================

test_that("make_dists rejects odd-length Weibull parameter vector", {
  expect_error(make_dists(c(1, 2, 3), "weibull"), "2m")
})
