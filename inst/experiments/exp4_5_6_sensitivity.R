#!/usr/bin/env Rscript
# Experiments 4, 5, 6: Sensitivity analyses (masking noise, delta, r)
devtools::load_all(".")

m <- 4; k_val <- 2; n <- 300; R <- 50
rates <- c(0.4, 0.6, 0.8, 1.0)
rates_sorted <- sort(rates)
model <- kofn(k = k_val, m = m, family = "exponential")

# --- Exp 4: Masking noise ---
cat("=== Exp 4: Masking noise ===\n")
exp4 <- list()
for (p_mask in c(0, 0.1, 0.3, 0.5, 1.0)) {
  maes <- numeric(R)
  for (rep in seq_len(R)) {
    set.seed(2026 + rep)
    gen_m <- rdata_masked(model)
    df_m <- gen_m(theta = rates, n = n, p_mask = p_mask)
    ll_m <- loglik_masked(model)
    neg_ll <- function(p) {
      val <- -ll_m(df_m, p)
      if (!is.finite(val)) .Machine$double.xmax / 2 else val
    }
    res <- tryCatch(
      suppressWarnings(multistart_mle(neg_ll, rep(0.5, m), m, n_starts = 5, nobs = n)),
      error = function(e) NULL)
    if (!is.null(res) && !any(is.na(coef(res)))) {
      maes[rep] <- mean(abs(sort(coef(res)) - rates_sorted))
    } else {
      maes[rep] <- NA
    }
  }
  exp4[[length(exp4) + 1]] <- data.frame(
    p_mask = p_mask, median_mae = median(maes, na.rm = TRUE),
    mean_mae = mean(maes, na.rm = TRUE), conv = sum(!is.na(maes)))
  cat(sprintf("  p_mask=%.1f: median_mae=%.3f\n", p_mask, median(maes, na.rm=TRUE)))
}
exp4 <- do.call(rbind, exp4)
saveRDS(exp4, "inst/precomputed/paper/exp4_masking_noise.rds")

# --- Exp 5: Inspection delta ---
cat("\n=== Exp 5: Inspection delta ===\n")
exp5 <- list()
for (delta in c(0.1, 0.5, 1.0, 2.0)) {
  maes <- numeric(R)
  for (rep in seq_len(R)) {
    set.seed(2026 + rep)
    s1gen <- rdata_scheme1(model)
    df1 <- s1gen(theta = rates, n = n, delta = delta)
    res <- tryCatch(fit_scheme1(model)(df1, n_starts = 5),
                    error = function(e) NULL)
    if (!is.null(res) && res$converged && !any(is.na(coef(res)))) {
      maes[rep] <- mean(abs(sort(coef(res)) - rates_sorted))
    } else {
      maes[rep] <- NA
    }
  }
  exp5[[length(exp5) + 1]] <- data.frame(
    delta = delta, median_mae = median(maes, na.rm = TRUE),
    mean_mae = mean(maes, na.rm = TRUE), conv = sum(!is.na(maes)))
  cat(sprintf("  delta=%.1f: median_mae=%.3f\n", delta, median(maes, na.rm=TRUE)))
}
exp5 <- do.call(rbind, exp5)
saveRDS(exp5, "inst/precomputed/paper/exp5_inspection_delta.rds")

# --- Exp 6: Autopsy coverage ---
cat("\n=== Exp 6: Autopsy coverage ===\n")
exp6 <- list()
for (r_inspect in 0:4) {
  maes <- numeric(R)
  for (rep in seq_len(R)) {
    set.seed(2026 + rep)
    comp_times <- matrix(0, nrow = n, ncol = m)
    for (j in seq_len(m)) comp_times[, j] <- stats::rexp(n, rate = rates[j])
    sys_times <- vapply(seq_len(n), function(i)
      system_lifetime(model$system, comp_times[i, ]), numeric(1))

    if (r_inspect == 0) {
      # Scheme 0
      res <- tryCatch(fit_system(sys_times, model$system, n_starts = 5),
                      error = function(e) NULL)
    } else if (r_inspect == m) {
      # Full autopsy = masked (p_mask=0)
      gen_m <- rdata_masked(model)
      df_m <- gen_m(theta = rates, n = n, p_mask = 0)
      ll_m <- loglik_masked(model)
      neg_ll <- function(p) {
        val <- -ll_m(df_m, p)
        if (!is.finite(val)) .Machine$double.xmax / 2 else val
      }
      res <- tryCatch(
        suppressWarnings(multistart_mle(neg_ll, rep(0.5, m), m, n_starts = 5, nobs = n)),
        error = function(e) NULL)
    } else {
      # Partial autopsy
      inspected <- matrix(FALSE, nrow = n, ncol = m)
      known_failed <- matrix(FALSE, nrow = n, ncol = m)
      for (i in seq_len(n)) {
        idx <- sample.int(m, r_inspect)
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
      res <- tryCatch(
        suppressWarnings(multistart_mle(neg_ll_aut, rep(0.5, m), m, n_starts = 5, nobs = n)),
        error = function(e) NULL)
    }

    if (!is.null(res) && !any(is.na(tryCatch(coef(res), error = function(e) NA)))) {
      est <- tryCatch(coef(res), error = function(e) res$par)
      maes[rep] <- mean(abs(sort(est) - rates_sorted))
    } else {
      maes[rep] <- NA
    }
  }
  exp6[[length(exp6) + 1]] <- data.frame(
    r = r_inspect, median_mae = median(maes, na.rm = TRUE),
    mean_mae = mean(maes, na.rm = TRUE), conv = sum(!is.na(maes)))
  cat(sprintf("  r=%d: median_mae=%.3f\n", r_inspect, median(maes, na.rm=TRUE)))
}
exp6 <- do.call(rbind, exp6)
saveRDS(exp6, "inst/precomputed/paper/exp6_autopsy_coverage.rds")

cat("\nExperiments 4, 5, 6 complete.\n")
