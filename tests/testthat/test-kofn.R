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
  model <- kofn(k = 1, m = 3, family = "exponential")
  expect_s3_class(model, "exp_kofn")
  expect_s3_class(model, "kofn")
  expect_equal(ncomponents(model), 3L)
  expect_true(is_kofn(model))

  model_wei <- kofn(k = 1, m = 2, family = "weibull", method = "em")
  expect_s3_class(model_wei, "wei_kofn")
})

test_that("IE expansion produces correct number of terms", {
  ie <- ie_expand(c(0.5, 0.3))
  expect_equal(length(ie$sign), 4)  # 2^2

  ie3 <- ie_expand(c(0.5, 0.3, 0.2))
  expect_equal(length(ie3$sign), 8)  # 2^3
})

test_that("exponential parallel data generation works", {
  model <- kofn(k = 1, m = 3)
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
  model <- kofn(k = 1, m = 2)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(0.5, 0.3), n = 100)

  ll <- loglik(model)
  val <- ll(df, par = c(0.5, 0.3))
  expect_true(is.finite(val))
  expect_true(val < 0)  # log-likelihood should be negative
})

test_that("exponential parallel MLE converges", {
  model <- kofn(k = 1, m = 2)
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(0.5, 0.3), n = 200)

  fitter <- fit(model)
  result <- fitter(df, n_starts = 3)
  expect_equal(result$convergence, 0L)
  expect_true(all(is.finite(result$par)))
  expect_true(all(result$par > 0))
})

test_that("Weibull EM converges for parallel system", {
  model <- kofn(k = 1, m = 2, family = "weibull", method = "em")
  gen <- rdata(model)
  set.seed(42)
  df <- gen(theta = c(1.5, 2.0, 2.0, 3.0), n = 200)

  fitter <- fit(model)
  result <- fitter(df, n_starts = 2)
  expect_equal(result$convergence, 0L)
  expect_true(all(result$shapes > 0))
  expect_true(all(result$scales > 0))
})

test_that("observation functors produce correct types", {
  obs0 <- observe_scheme0(tau = 5)
  expect_equal(obs0(3)$omega, "exact")
  expect_equal(obs0(7)$omega, "right")

  obs1 <- observe_scheme1(delta = 1.0)
  expect_equal(obs1(2.7)$omega, "interval")
  expect_equal(obs1(2.7)$t, 2.0)
  expect_equal(obs1(2.7)$t_upper, 3.0)
})

test_that("Scheme 1 data generation works", {
  model <- kofn(k = 1, m = 2, family = "weibull")
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
