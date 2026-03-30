# ===========================================================================
# Shared helpers for simulation experiments
# ===========================================================================
#
# Source this file from experiment scripts to avoid duplicating common
# simulation infrastructure: component time generation, Scheme 1 data
# frame construction, MAE extraction, and the partial autopsy likelihood.
# ===========================================================================


#' Generate component lifetimes and system times
#'
#' @param n Number of observations.
#' @param m Number of components.
#' @param rates Rate vector (exponential) or parsed param list (Weibull).
#' @param system A coherent_system object.
#' @param family "exponential" or "weibull".
#' @param shapes,scales Weibull parameters (used when family = "weibull").
#' @return List with comp_times (n x m matrix) and sys_times (length-n vector).
generate_system_data <- function(n, m, system, family = "exponential",
                                 rates = NULL, shapes = NULL, scales = NULL) {
  comp_times <- matrix(0, nrow = n, ncol = m)
  if (family == "exponential") {
    for (j in seq_len(m)) {
      comp_times[, j] <- stats::rexp(n, rate = rates[j])
    }
  } else {
    for (j in seq_len(m)) {
      comp_times[, j] <- stats::rweibull(n, shape = shapes[j],
                                         scale = scales[j])
    }
  }
  sys_times <- vapply(seq_len(n), function(i)
    system_lifetime(system, comp_times[i, ]), numeric(1))
  list(comp_times = comp_times, sys_times = sys_times)
}


#' Build Scheme 1 data frame from component times
#'
#' @param sys_times System lifetime vector.
#' @param comp_times Component lifetime matrix (n x m).
#' @param delta Inspection interval width.
#' @return Data frame with t and comp_lower_j, comp_upper_j columns.
build_scheme1_df <- function(sys_times, comp_times, delta = 0.5) {
  m <- ncol(comp_times)
  df <- data.frame(t = sys_times)
  comp_lower <- floor(comp_times / delta) * delta
  for (j in seq_len(m)) {
    df[[paste0("comp_lower_", j)]] <- comp_lower[, j]
    df[[paste0("comp_upper_", j)]] <- comp_lower[, j] + delta
  }
  df
}


#' Extract MAE from a fit result (sorted estimates vs sorted truth)
#'
#' @param res Fit result object (fisher_mle or list with $par).
#' @param true_sorted Sorted true parameter vector.
#' @return Scalar MAE, or NA on failure.
get_mae <- function(res, true_sorted) {
  if (is.null(res)) return(NA)
  est <- tryCatch(coef(res), error = function(e) NULL)
  if (is.null(est) || any(is.na(est))) return(NA)
  mean(abs(sort(est) - true_sorted))
}


#' Negative log-likelihood for partial autopsy
#'
#' For a k-out-of-n exponential system where r components are inspected
#' post-failure. Marginalizes over unknown component states.
#'
#' @param sys_times System lifetime vector.
#' @param inspected Logical matrix (n x m): which components were inspected.
#' @param known_failed Logical matrix (n x m): inspected components that failed.
#' @param k System failure threshold.
#' @param m Number of components.
#' @return Function(par) returning negative log-likelihood.
make_autopsy_negll <- function(sys_times, inspected, known_failed, k, m) {
  n <- length(sys_times)

  function(par) {
    if (any(par <= 0)) return(.Machine$double.xmax / 2)
    ll <- 0
    for (i in seq_len(n)) {
      ti <- sys_times[i]
      f_v <- stats::dexp(ti, rate = par)
      F_v <- stats::pexp(ti, rate = par)
      S_v <- 1 - F_v
      ks <- inspected[i, ] & !known_failed[i, ]
      kf <- inspected[i, ] & known_failed[i, ]
      unk_idx <- which(!inspected[i, ])
      n_unk <- length(unk_idx)
      total <- 0
      for (mask in seq_len(2^n_unk) - 1L) {
        bits <- as.logical(as.integer(intToBits(mask))[seq_len(n_unk)])
        Fc <- c(which(kf), unk_idx[bits])
        if (length(Fc) != k) next
        if (any(ks[Fc])) next
        nF <- setdiff(seq_len(m), Fc)
        for (j in Fc) {
          if (f_v[j] <= 0) next
          oF <- setdiff(Fc, j)
          lt <- log(f_v[j])
          if (length(oF) > 0) {
            if (any(F_v[oF] <= 0)) next
            lt <- lt + sum(log(F_v[oF]))
          }
          if (length(nF) > 0) {
            if (any(S_v[nF] <= 0)) next
            lt <- lt + sum(log(S_v[nF]))
          }
          total <- total + exp(lt)
        }
      }
      if (total <= 0) return(.Machine$double.xmax / 2)
      ll <- ll + log(total)
    }
    -ll
  }
}


#' Generate partial autopsy data
#'
#' Randomly inspects r_inspect components per observation and records
#' whether each inspected component had failed by the system failure time.
#'
#' @param comp_times Component lifetime matrix (n x m).
#' @param sys_times System lifetime vector.
#' @param r_inspect Number of components to inspect per observation.
#' @param m Number of components.
#' @return List with inspected and known_failed logical matrices.
generate_autopsy_data <- function(comp_times, sys_times, r_inspect, m) {
  n <- length(sys_times)
  inspected <- matrix(FALSE, nrow = n, ncol = m)
  known_failed <- matrix(FALSE, nrow = n, ncol = m)
  for (i in seq_len(n)) {
    idx <- sample.int(m, r_inspect)
    inspected[i, idx] <- TRUE
    known_failed[i, idx] <- comp_times[i, idx] <= sys_times[i]
  }
  list(inspected = inspected, known_failed = known_failed)
}
