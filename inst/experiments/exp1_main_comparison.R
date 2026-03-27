#!/usr/bin/env Rscript
# Experiment 1: Main comparison of 6 symmetry-breaking mechanisms
# m=4, k=2, exponential, n=300, R=100
devtools::load_all(".")

m <- 4; k_val <- 2; n <- 300; R <- 100
rates <- c(0.4, 0.6, 0.8, 1.0)
rates_sorted <- sort(rates)
model <- kofn(k = k_val, m = m, family = "exponential")

results <- list()

for (rep in seq_len(R)) {
  set.seed(2026 + rep)

  # Generate shared component lifetimes for fair comparison
  comp_times <- matrix(0, nrow = n, ncol = m)
  for (j in seq_len(m)) comp_times[, j] <- stats::rexp(n, rate = rates[j])

  sys_times <- vapply(seq_len(n), function(i)
    system_lifetime(model$system, comp_times[i, ]), numeric(1))

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
  comp_times2 <- matrix(0, nrow = n/2, ncol = m)
  for (j in seq_len(m)) comp_times2[, j] <- stats::rexp(n/2, rate = rates[j])
  sys_times_k2 <- sys_times[seq_len(n/2)]
  sys_times_k3 <- vapply(seq_len(n/2), function(i)
    system_lifetime(sys_k3, comp_times2[i, ]), numeric(1))
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
  s1gen <- rdata_scheme1(model)
  # Rebuild from comp_times for fair comparison
  df_s1 <- data.frame(t = sys_times)
  for (j in seq_len(m)) {
    lower_j <- floor(comp_times[, j] / 0.5) * 0.5
    df_s1[[paste0("comp_lower_", j)]] <- lower_j
    df_s1[[paste0("comp_upper_", j)]] <- lower_j + 0.5
  }
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
  inspected <- matrix(FALSE, nrow = n, ncol = m)
  known_failed <- matrix(FALSE, nrow = n, ncol = m)
  for (i in seq_len(n)) {
    idx <- sample.int(m, 2)
    inspected[i, idx] <- TRUE
    known_failed[i, idx] <- comp_times[i, idx] <= sys_times[i]
  }
  neg_ll_aut <- function(par) {
    if (any(par <= 0)) return(.Machine$double.xmax / 2)
    ll <- 0
    for (i in seq_len(n)) {
      ti <- sys_times[i]
      f_v <- stats::dexp(ti, rate = par)
      F_v <- stats::pexp(ti, rate = par)
      S_v <- 1 - F_v
      ks <- inspected[i, ] & !known_failed[i, ]
      kf <- inspected[i, ] & known_failed[i, ]
      unk <- !inspected[i, ]
      unk_idx <- which(unk)
      n_unk <- length(unk_idx)
      total <- 0
      for (mask in seq_len(2^n_unk) - 1L) {
        bits <- as.logical(as.integer(intToBits(mask))[seq_len(n_unk)])
        Fc <- c(which(kf), unk_idx[bits])
        if (length(Fc) != k_val) next
        if (any(ks[Fc])) next
        nF <- setdiff(seq_len(m), Fc)
        for (j in Fc) {
          if (f_v[j] <= 0) next
          oF <- setdiff(Fc, j)
          lt <- log(f_v[j])
          if (length(oF) > 0) { if (any(F_v[oF] <= 0)) next; lt <- lt + sum(log(F_v[oF])) }
          if (length(nF) > 0) { if (any(S_v[nF] <= 0)) next; lt <- lt + sum(log(S_v[nF])) }
          total <- total + exp(lt)
        }
      }
      if (total <= 0) return(.Machine$double.xmax / 2)
      ll <- ll + log(total)
    }
    -ll
  }
  res_aut <- tryCatch(
    suppressWarnings(multistart_mle(neg_ll_aut, rep(0.5, m), m, n_starts = 5, nobs = n)),
    error = function(e) NULL)

  # Collect results
  get_mae <- function(res) {
    if (is.null(res)) return(NA)
    est <- tryCatch(coef(res), error = function(e) res$par)
    if (any(is.na(est))) return(NA)
    mean(abs(sort(est) - rates_sorted))
  }
  get_mae_ord <- function(res) {
    if (res$convergence != 0) return(NA)
    mean(abs(cumsum(res$par) - rates_sorted))
  }

  results[[rep]] <- data.frame(
    rep = rep,
    scheme0 = get_mae(res0),
    ordering = get_mae_ord(res_ord),
    heterogeneous = get_mae(res_het),
    periodic = get_mae(res1),
    masked = get_mae(res_mask),
    partial_autopsy = get_mae(res_aut),
    stringsAsFactors = FALSE
  )

  if (rep %% 10 == 0) cat(sprintf("Exp1: %d/%d done\n", rep, R))
}

exp1 <- do.call(rbind, results)
saveRDS(exp1, "inst/precomputed/paper/exp1_main_comparison.rds")
cat("Experiment 1 complete.\n")
print(round(colMeans(exp1[, -1], na.rm = TRUE), 3))
