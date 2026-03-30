#!/usr/bin/env Rscript
# Experiment 2: k-spectrum (Scheme 0 vs periodic across k=2,3,4)
devtools::load_all(".")
source("inst/experiments/helpers.R")

m <- 4; n <- 300; R <- 100
rates <- c(0.4, 0.6, 0.8, 1.0)
rates_sorted <- sort(rates)

results <- list()
for (k_val in c(2, 3, 4)) {
  model <- kofn(k = k_val, m = m, family = "exponential")
  for (rep in seq_len(R)) {
    set.seed(2026 + rep)
    sim <- generate_system_data(n, m, model$system, rates = rates)

    res0 <- tryCatch(fit_system(sim$sys_times, model$system, n_starts = 5),
                     error = function(e) NULL)

    df_s1 <- build_scheme1_df(sim$sys_times, sim$comp_times)
    res1 <- tryCatch(fit_scheme1(model)(df_s1, n_starts = 5),
                     error = function(e) NULL)

    results[[length(results) + 1]] <- data.frame(
      k = k_val, rep = rep,
      scheme0 = get_mae(res0, rates_sorted),
      periodic = get_mae(res1, rates_sorted),
      stringsAsFactors = FALSE)
  }
  cat(sprintf("Exp2: k=%d done\n", k_val))
}

exp2 <- do.call(rbind, results)
saveRDS(exp2, "inst/precomputed/paper/exp2_k_spectrum.rds")
cat("Experiment 2 complete.\n")
for (kv in c(2, 3, 4)) {
  d <- exp2[exp2$k == kv, ]
  cat(sprintf("k=%d: S0=%.3f, S1=%.3f, ratio=%.1fx\n", kv,
    median(d$scheme0, na.rm = TRUE), median(d$periodic, na.rm = TRUE),
    median(d$scheme0, na.rm = TRUE) / median(d$periodic, na.rm = TRUE)))
}
