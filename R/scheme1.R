# ===========================================================================
# Scheme 1: Periodic Inspection Likelihood and Fitting
# ===========================================================================
#
# Under Scheme 1 (periodic inspection at interval delta), each component's
# failure time is localized to an inspection interval [a_j-, a_j+).
# The system failure time T = max(T_1, ..., T_m) is observed exactly.
#
# The likelihood for observation i is:
#   L_i(theta) = f_sys(t_i; theta) * prod_j [F_j(a_ij+) - F_j(a_ij-)]
#
# where [a_ij-, a_ij+) is the inspection interval containing T_j.
#
# This module provides data generation, log-likelihood computation, and
# MLE fitting for Scheme 1 observations. Supports both exponential and
# Weibull component distributions.
#
# Data convention: data frames have columns t (system failure time),
# comp_lower_j, comp_upper_j (interval bounds for component j).
# ===========================================================================


# ---------------------------------------------------------------------------
# rdata_scheme1: Data generation under periodic inspection
# ---------------------------------------------------------------------------

#' Generate Scheme 1 (periodic inspection) data
#'
#' Returns a closure that generates parallel system data under periodic
#' inspection. Each observation yields the exact system failure time and
#' interval-censored component failure times.
#'
#' @param model A \code{kofn} model object (exponential or Weibull).
#' @param ... Additional arguments (currently unused).
#' @return A function \code{function(theta, n, delta)} returning a data frame
#'   with columns:
#'   \describe{
#'     \item{\code{t}}{System failure time (exact).}
#'     \item{\code{comp_lower_j}}{Lower bound of inspection interval for
#'       component j.}
#'     \item{\code{comp_upper_j}}{Upper bound of inspection interval for
#'       component j.}
#'   }
#'   The data frame has attributes \code{comp_times} (true component times),
#'   \code{delta} (inspection interval), and \code{par} (true parameters).
#'
#' @export
rdata_scheme1 <- function(model, ...) {
  m <- model$m
  lt <- model$lifetime
  family <- if (inherits(model, "wei_kofn")) "weibull" else "exponential"

  function(theta, n, delta = 1.0) {
    # Parse parameters by family
    if (family == "exponential") {
      stopifnot(length(theta) == m)
      shapes <- rep(1, m)
      scales <- 1 / theta
    } else {
      stopifnot(length(theta) == 2L * m)
      shapes <- theta[seq(1L, 2L * m, by = 2L)]
      scales <- theta[seq(2L, 2L * m, by = 2L)]
    }
    stopifnot(all(shapes > 0), all(scales > 0), delta > 0)

    # Generate component lifetimes
    comp_times <- matrix(0, nrow = n, ncol = m)
    for (j in seq_len(m)) {
      comp_times[, j] <- stats::rweibull(n, shape = shapes[j],
                                         scale = scales[j])
    }

    # System lifetime = max (parallel)
    sys_times <- apply(comp_times, 1, max)

    # Interval-censor each component to inspection grid
    comp_lower <- matrix(0, nrow = n, ncol = m)
    comp_upper <- matrix(0, nrow = n, ncol = m)
    for (j in seq_len(m)) {
      comp_lower[, j] <- floor(comp_times[, j] / delta) * delta
      comp_upper[, j] <- comp_lower[, j] + delta
    }

    # Build data frame
    df <- data.frame(sys_times)
    names(df) <- lt
    for (j in seq_len(m)) {
      df[[paste0("comp_lower_", j)]] <- comp_lower[, j]
      df[[paste0("comp_upper_", j)]] <- comp_upper[, j]
    }

    attr(df, "comp_times") <- comp_times
    attr(df, "delta") <- delta
    attr(df, "family") <- family
    attr(df, "par") <- theta

    df
  }
}


# ---------------------------------------------------------------------------
# loglik_scheme1: Log-likelihood under periodic inspection
# ---------------------------------------------------------------------------

#' Log-likelihood for Scheme 1 (periodic inspection) parallel system
#'
#' Returns a closure that computes the log-likelihood under Scheme 1
#' observation. The likelihood combines the system density at the exact
#' system failure time with interval-censored component contributions.
#'
#' @param model A \code{kofn} model object (exponential or Weibull).
#' @param ... Additional arguments (currently unused).
#' @return A function \code{function(df, par)} returning a scalar
#'   log-likelihood.
#'
#' @details
#' The log-likelihood for observation i is:
#' \deqn{\log L_i(\theta) = \log f_{sys}(t_i) + \sum_j \log[F_j(a_{ij}^+) - F_j(a_{ij}^-)]}
#'
#' where \eqn{f_{sys}} is the parallel system density and
#' \eqn{[a_{ij}^-, a_{ij}^+)} is the inspection interval containing
#' component j's failure time.
#'
#' @export
loglik_scheme1 <- function(model, ...) {
  m <- model$m
  lt <- model$lifetime
  family <- if (inherits(model, "wei_kofn")) "weibull" else "exponential"

  lower_cols <- paste0("comp_lower_", seq_len(m))
  upper_cols <- paste0("comp_upper_", seq_len(m))

  function(df, par) {
    if (any(!is.finite(par)) || any(par <= 0)) return(-Inf)
    n <- nrow(df)

    # Parse parameters by family
    if (family == "exponential") {
      shapes <- rep(1, m)
      scales <- 1 / par
    } else {
      shapes <- par[seq(1L, 2L * m, by = 2L)]
      scales <- par[seq(2L, 2L * m, by = 2L)]
    }

    ll <- 0
    for (i in seq_len(n)) {
      ti <- df[[lt]][i]

      # System density: f_sys(t_i) = sum_j f_j(t_i) * prod_{k!=j} F_k(t_i)
      wj <- numeric(m)
      for (j in seq_len(m)) {
        f_j <- stats::dweibull(ti, shape = shapes[j], scale = scales[j])
        if (f_j <= 0) next
        log_prod <- 0
        for (k in seq_len(m)) {
          if (k == j) next
          F_k <- stats::pweibull(ti, shape = shapes[k], scale = scales[k])
          if (F_k <= 0) return(-Inf)
          log_prod <- log_prod + log(F_k)
        }
        wj[j] <- f_j * exp(log_prod)
      }
      f_sys <- sum(wj)
      if (f_sys <= 0) return(-Inf)
      ll <- ll + log(f_sys)

      # Component interval contributions: prod_j [F_j(upper) - F_j(lower)]
      for (j in seq_len(m)) {
        a_lower <- df[[lower_cols[j]]][i]
        a_upper <- df[[upper_cols[j]]][i]
        F_upper <- stats::pweibull(a_upper, shape = shapes[j],
                                   scale = scales[j])
        F_lower <- stats::pweibull(a_lower, shape = shapes[j],
                                   scale = scales[j])
        prob <- F_upper - F_lower
        if (prob <= 0) return(-Inf)
        ll <- ll + log(prob)
      }
    }

    ll
  }
}


# ---------------------------------------------------------------------------
# fit_scheme1: MLE fitting under periodic inspection
# ---------------------------------------------------------------------------

#' Fit parallel system MLE under Scheme 1 (periodic inspection)
#'
#' Returns a closure that fits the model to Scheme 1 data using multi-start
#' optimization with L-BFGS-B and Nelder-Mead fallback.
#'
#' @param model A \code{kofn} model object (exponential or Weibull).
#' @param ... Additional arguments (currently unused).
#' @return A function \code{function(df, par0 = NULL, n_starts = 5L)} that
#'   returns a \code{fisher_mle} object (from the likelihood.model package).
#'
#' @details
#' The solver uses L-BFGS-B as the primary optimization method with
#' positivity constraints, falling back to Nelder-Mead on the log-parameter
#' scale if L-BFGS-B fails to converge.
#'
#' Standard errors are computed from the numerical Hessian at the MLE.
#'
#' @export
fit_scheme1 <- function(model, ...) {
  m <- model$m
  lt <- model$lifetime
  family <- if (inherits(model, "wei_kofn")) "weibull" else "exponential"
  ll_closure <- loglik_scheme1(model, ...)
  n_par_per_comp <- switch(family, exponential = 1L, weibull = 2L)
  n_par <- m * n_par_per_comp

  function(df, par0 = NULL, n_starts = 5L) {
    # Default initialization
    if (is.null(par0)) {
      mean_t <- mean(df[[lt]])
      H_m <- sum(1 / seq_len(m))
      if (family == "exponential") {
        lam0 <- H_m / max(mean_t, 0.01)
        par0 <- lam0 * seq(0.5, 1.5, length.out = m)
      } else {
        scale0 <- max(mean_t, 0.01) / H_m
        par0 <- as.numeric(rbind(
          rep(1.0, m),
          scale0 * seq(0.5, 1.5, length.out = m)
        ))
      }
    }
    stopifnot(length(par0) == n_par)

    neg_ll <- function(p) {
      val <- -ll_closure(df, p)
      if (!is.finite(val)) return(.Machine$double.xmax / 2)
      val
    }

    neg_ll_log <- function(log_p) neg_ll(exp(log_p))

    fit_from_init <- function(init_par) {
      result <- tryCatch(
        stats::optim(par = init_par, fn = neg_ll, method = "L-BFGS-B",
                     lower = rep(1e-10, n_par), upper = rep(Inf, n_par)),
        error = function(e) {
          list(par = rep(NA_real_, n_par), value = Inf, convergence = 99L)
        }
      )
      if (result$convergence != 0) {
        result2 <- tryCatch(
          stats::optim(par = log(pmax(init_par, 1e-10)), fn = neg_ll_log,
                       method = "Nelder-Mead"),
          error = function(e) {
            list(par = rep(NA_real_, n_par), value = Inf, convergence = 99L)
          }
        )
        if (result2$convergence == 0 && result2$value < result$value) {
          result <- list(par = exp(result2$par), value = result2$value,
                         convergence = result2$convergence)
        }
      }
      result
    }

    # Multi-start
    best <- fit_from_init(par0)
    if (n_starts > 1L) {
      for (s in seq_len(n_starts - 1L)) {
        init_s <- par0 * exp(stats::rnorm(n_par, sd = 0.5))
        init_s <- pmax(init_s, 1e-10)
        res_s <- fit_from_init(init_s)
        if (res_s$convergence == 0 &&
            (best$convergence != 0 || res_s$value < best$value)) {
          best <- res_s
        }
      }
    }

    # Hessian of negative log-likelihood and score at MLE
    H <- tryCatch(
      numDeriv::hessian(neg_ll, x = best$par, method = "Richardson"),
      error = function(e) NULL
    )
    hessian_mat <- if (!is.null(H) && all(is.finite(H))) H else {
      matrix(NA_real_, nrow = n_par, ncol = n_par)
    }
    score_val <- tryCatch(
      numDeriv::grad(func = function(p) ll_closure(df, p), x = best$par,
                     method = "Richardson"),
      error = function(e) rep(NA_real_, n_par)
    )

    make_mle_result(
      par       = best$par,
      loglik    = -best$value,
      hessian   = hessian_mat,
      score     = score_val,
      nobs      = nrow(df),
      converged = (best$convergence == 0)
    )
  }
}
