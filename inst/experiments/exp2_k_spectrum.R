#!/usr/bin/env Rscript
# Experiment 2: k-spectrum (Scheme 0 vs periodic across k=2,3,4)
devtools::load_all(".")

m <- 4; n <- 300; R <- 100
rates <- c(0.4, 0.6, 0.8, 1.0)
rates_sorted <- sort(rates)

results <- list()
for (k_val in c(2, 3, 4)) {
  model <- kofn(k = k_val, m = m, family = "exponential")
  for (rep in seq_len(R)) {
    set.seed(2026 + rep)
    comp_times <- matrix(0, nrow = n, ncol = m)
    for (j in seq_len(m)) comp_times[, j] <- stats::rexp(n, rate = rates[j])
    sys_times <- vapply(seq_len(n), function(i)
      system_lifetime(model$system, comp_times[i, ]), numeric(1))

    res0 <- tryCatch(fit_system(sys_times, model$system, n_starts = 5),
                     error = function(e) NULL)

    df_s1 <- data.frame(t = sys_times)
    for (j in seq_len(m)) {
      lower_j <- floor(comp_times[, j] / 0.5) * 0.5
      df_s1[[paste0("comp_lower_", j)]] <- lower_j
      df_s1[[paste0("comp_upper_", j)]] <- lower_j + 0.5
    }
    res1 <- tryCatch(fit_scheme1(model)(df_s1, n_starts = 5),
                     error = function(e) NULL)

    get_mae <- function(res) {
      if (is.null(res)) return(NA)
      est <- tryCatch(coef(res), error = function(e) NULL)
      if (is.null(est) || any(is.na(est))) return(NA)
      mean(abs(sort(est) - rates_sorted))
    }

    results[[length(results) + 1]] <- data.frame(
      k = k_val, rep = rep,
      scheme0 = get_mae(res0), periodic = get_mae(res1),
      stringsAsFactors = FALSE)
  }
  cat(sprintf("Exp2: k=%d done\n", k_val))
}

exp2 <- do.call(rbind, results)
saveRDS(exp2, "inst/precomputed/paper/exp2_k_spectrum.rds")
cat("Experiment 2 complete.\n")
for (kv in c(2,3,4)) {
  d <- exp2[exp2$k == kv, ]
  cat(sprintf("k=%d: S0=%.3f, S1=%.3f, ratio=%.1fx\n", kv,
    median(d$scheme0, na.rm=TRUE), median(d$periodic, na.rm=TRUE),
    median(d$scheme0, na.rm=TRUE) / median(d$periodic, na.rm=TRUE)))
}
