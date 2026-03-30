#!/usr/bin/env Rscript
# Experiment 3: Scaling with m + Experiment 7: Weibull extension
devtools::load_all(".")
source("inst/experiments/helpers.R")

# === Exp 3: Scaling with m ===
cat("=== Exp 3: Scaling with m ===\n")
R <- 50

all_rates <- list(
  `3` = c(0.4, 0.7, 1.0),
  `4` = c(0.4, 0.6, 0.8, 1.0),
  `5` = c(0.4, 0.6, 0.8, 1.0, 1.2)
)

exp3 <- list()
for (m_str in names(all_rates)) {
  m <- as.integer(m_str)
  rates <- all_rates[[m_str]]
  rates_sorted <- sort(rates)
  k_val <- ceiling(m / 2)
  n <- 300
  model <- kofn(k = k_val, m = m, family = "exponential")

  maes_s0 <- maes_per <- maes_mask <- numeric(R)
  for (rep in seq_len(R)) {
    set.seed(2026 + rep)
    sim <- generate_system_data(n, m, model$system, rates = rates)

    # Scheme 0
    res0 <- tryCatch(fit_system(sim$sys_times, model$system, n_starts = 5),
                     error = function(e) NULL)
    # Periodic
    df_s1 <- build_scheme1_df(sim$sys_times, sim$comp_times)
    res1 <- tryCatch(fit_scheme1(model)(df_s1, n_starts = 5),
                     error = function(e) NULL)
    # Masked (full autopsy)
    gen_m <- rdata_masked(model)
    df_m <- gen_m(theta = rates, n = n, p_mask = 0)
    ll_m <- loglik_masked(model)
    neg_ll <- function(p) {
      val <- -ll_m(df_m, p)
      if (!is.finite(val)) .Machine$double.xmax / 2 else val
    }
    res_mask <- tryCatch(
      suppressWarnings(multistart_mle(neg_ll, rep(0.5, m), m, n_starts = 5, nobs = n)),
      error = function(e) NULL)

    maes_s0[rep] <- get_mae(res0, rates_sorted)
    maes_per[rep] <- get_mae(res1, rates_sorted)
    maes_mask[rep] <- get_mae(res_mask, rates_sorted)
  }

  exp3[[length(exp3) + 1]] <- data.frame(
    m = m, k = k_val,
    s0_median = median(maes_s0, na.rm = TRUE),
    periodic_median = median(maes_per, na.rm = TRUE),
    autopsy_median = median(maes_mask, na.rm = TRUE))

  cat(sprintf("m=%d k=%d: S0=%.3f, Per=%.3f, Aut=%.3f\n", m, k_val,
    median(maes_s0, na.rm = TRUE), median(maes_per, na.rm = TRUE),
    median(maes_mask, na.rm = TRUE)))
}
exp3 <- do.call(rbind, exp3)
saveRDS(exp3, "inst/precomputed/paper/exp3_scaling.rds")

# === Exp 7: Weibull extension ===
cat("\n=== Exp 7: Weibull extension ===\n")
m <- 4; k_val <- 2; n <- 300; R_w <- 50
# Weibull shape=1.5, scales chosen so mean ~ 1/rate
wei_par <- c(1.5, 2.5, 1.5, 1.7, 1.5, 1.25, 1.5, 1.0)  # interleaved shape, scale
model_w <- kofn(k = k_val, m = m, family = "weibull")
pp <- parse_params(wei_par, m, "weibull")
wei_par_sorted <- sort(wei_par)

maes_s0 <- maes_per <- maes_mask <- numeric(R_w)
for (rep in seq_len(R_w)) {
  set.seed(2026 + rep)
  sim <- generate_system_data(n, m, model_w$system, family = "weibull",
                              shapes = pp$shapes, scales = pp$scales)

  # Scheme 0
  res0 <- tryCatch(fit_system(sim$sys_times, model_w$system, family = "weibull", n_starts = 5),
                   error = function(e) NULL)
  # Periodic
  df_s1 <- build_scheme1_df(sim$sys_times, sim$comp_times)
  res1 <- tryCatch(fit_scheme1(model_w)(df_s1, n_starts = 5),
                   error = function(e) NULL)
  # Masked (full autopsy)
  gen_m <- rdata_masked(model_w, candset = "c")
  df_m <- gen_m(theta = wei_par, n = n, p_mask = 0)
  ll_m <- loglik_masked(model_w, candset = "c")
  neg_ll <- function(p) {
    val <- -ll_m(df_m, p)
    if (!is.finite(val)) .Machine$double.xmax / 2 else val
  }
  res_mask <- tryCatch(
    suppressWarnings(multistart_mle(neg_ll, rep(1, 2 * m), 2 * m, n_starts = 5, nobs = n)),
    error = function(e) NULL)

  maes_s0[rep] <- get_mae(res0, wei_par_sorted)
  maes_per[rep] <- get_mae(res1, wei_par_sorted)
  maes_mask[rep] <- get_mae(res_mask, wei_par_sorted)

  if (rep %% 10 == 0) cat(sprintf("Exp7: %d/%d done\n", rep, R_w))
}

exp7 <- data.frame(
  mechanism = c("scheme0", "periodic", "masked"),
  median_mae = c(median(maes_s0, na.rm = TRUE), median(maes_per, na.rm = TRUE),
                 median(maes_mask, na.rm = TRUE)),
  conv = c(sum(!is.na(maes_s0)), sum(!is.na(maes_per)), sum(!is.na(maes_mask))))
saveRDS(exp7, "inst/precomputed/paper/exp7_weibull.rds")

cat("Experiments 3 and 7 complete.\n")
print(exp3)
print(exp7)
