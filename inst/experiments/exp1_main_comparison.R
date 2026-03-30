#!/usr/bin/env Rscript
# Experiment 1: Main comparison of 6 symmetry-breaking mechanisms
# m=4, k=2, exponential, n=300, R=100
devtools::load_all(".")
source("inst/experiments/helpers.R")

m <- 4; k_val <- 2; n <- 300; R <- 100
rates <- c(0.4, 0.6, 0.8, 1.0)
rates_sorted <- sort(rates)
model <- kofn(k = k_val, m = m, family = "exponential")

results <- list()

for (rep in seq_len(R)) {
  set.seed(2026 + rep)

  # Generate shared component lifetimes for fair comparison
  sim <- generate_system_data(n, m, model$system, rates = rates)
  comp_times <- sim$comp_times
  sys_times <- sim$sys_times

  # --- Scheme 0 ---
  res0 <- tryCatch(fit_system(sys_times, model$system, n_starts = 5),
                   error = function(e) NULL)

  # --- Ordering constraints ---
  ll_sys <- function(par) loglik_system(sys_times, model$system, par)
  neg_ll_ord <- function(delta) {
    if (any(delta <= 0)) return(.Machine$double.xmax / 2)
    val <- -ll_sys(cumsum(delta))
    if (!is.finite(val)) .Machine$double.xmax / 2 else val
  }
  res_ord <- tryCatch(
    stats::optim(diff(c(0, rep(mean(rates), m))), neg_ll_ord,
                 method = "L-BFGS-B", lower = rep(1e-10, m)),
    error = function(e) list(par = rep(NA, m), convergence = 99))

  # --- Heterogeneous k ---
  sys_k3 <- kofn_system(3, m)
  sim2 <- generate_system_data(n / 2, m, sys_k3, rates = rates)
  sys_times_k2 <- sys_times[seq_len(n / 2)]
  sys_times_k3 <- sim2$sys_times
  neg_ll_het <- function(par) {
    if (any(par <= 0)) return(.Machine$double.xmax / 2)
    val <- -(loglik_system(sys_times_k2, model$system, par) +
             loglik_system(sys_times_k3, sys_k3, par))
    if (!is.finite(val)) .Machine$double.xmax / 2 else val
  }
  res_het <- tryCatch(
    suppressWarnings(multistart_mle(neg_ll_het, rep(0.5, m), m, n_starts = 5, nobs = n)),
    error = function(e) NULL)

  # --- Periodic inspection ---
  df_s1 <- build_scheme1_df(sys_times, comp_times)
  fitter1 <- fit_scheme1(model)
  res1 <- tryCatch(fitter1(df_s1, n_starts = 5), error = function(e) NULL)

  # --- Masked failed sets ---
  gen_m <- rdata_masked(model)
  df_m <- gen_m(theta = rates, n = n, p_mask = 0)
  ll_m <- loglik_masked(model)
  neg_ll_m <- function(p) {
    val <- -ll_m(df_m, p)
    if (!is.finite(val)) .Machine$double.xmax / 2 else val
  }
  res_mask <- tryCatch(
    suppressWarnings(multistart_mle(neg_ll_m, rep(0.5, m), m, n_starts = 5, nobs = n)),
    error = function(e) NULL)

  # --- Partial autopsy (r=2) ---
  aut_data <- generate_autopsy_data(comp_times, sys_times, r_inspect = 2, m = m)
  neg_ll_aut <- make_autopsy_negll(sys_times, aut_data$inspected,
                                   aut_data$known_failed, k_val, m)
  res_aut <- tryCatch(
    suppressWarnings(multistart_mle(neg_ll_aut, rep(0.5, m), m, n_starts = 5, nobs = n)),
    error = function(e) NULL)

  # Collect results
  get_mae_ord <- function(res) {
    if (res$convergence != 0) return(NA)
    mean(abs(cumsum(res$par) - rates_sorted))
  }

  results[[rep]] <- data.frame(
    rep = rep,
    scheme0 = get_mae(res0, rates_sorted),
    ordering = get_mae_ord(res_ord),
    heterogeneous = get_mae(res_het, rates_sorted),
    periodic = get_mae(res1, rates_sorted),
    masked = get_mae(res_mask, rates_sorted),
    partial_autopsy = get_mae(res_aut, rates_sorted),
    stringsAsFactors = FALSE
  )

  if (rep %% 10 == 0) cat(sprintf("Exp1: %d/%d done\n", rep, R))
}

exp1 <- do.call(rbind, results)
saveRDS(exp1, "inst/precomputed/paper/exp1_main_comparison.rds")
cat("Experiment 1 complete.\n")
print(round(colMeans(exp1[, -1], na.rm = TRUE), 3))
