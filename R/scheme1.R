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
#' @examples
#' model <- kofn(k = 1, m = 2, family = "exponential")
#' gen <- rdata_scheme1(model)
#' set.seed(1)
#' df <- gen(theta = c(1, 2), n = 20, delta = 1.0)
#' head(df)
#'
#' @export
rdata_scheme1 <- function(model, ...) {
  m <- model$m
  lt <- model$lifetime
  family <- if (inherits(model, "wei_kofn")) "weibull" else "exponential"

  function(theta, n, delta = 1.0) {
    n_par_expected <- if (family == "exponential") m else 2L * m
    stopifnot(length(theta) == n_par_expected)
    pp <- parse_params(theta, m, family)
    shapes <- pp$shapes
    scales <- pp$scales
    stopifnot(all(shapes > 0), all(scales > 0), delta > 0)

    # Generate component lifetimes
    comp_times <- matrix(0, nrow = n, ncol = m)
    for (j in seq_len(m)) {
      comp_times[, j] <- stats::rweibull(n, shape = shapes[j],
                                         scale = scales[j])
    }

    # System lifetime = max (parallel)
    sys_times <- apply(comp_times, 1, max)

    # Interval-censor each component to inspection grid (vectorized)
    comp_lower <- floor(comp_times / delta) * delta
    comp_upper <- comp_lower + delta

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
#' @examples
#' model <- kofn(k = 1, m = 2, family = "exponential")
#' ll <- loglik_scheme1(model)
#' set.seed(1)
#' df <- rdata_scheme1(model)(c(1, 2), n = 30, delta = 1.0)
#' ll(df, c(1, 2))
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

    pp <- parse_params(par, m, family)
    shapes <- pp$shapes
    scales <- pp$scales

    ll <- 0
    for (i in seq_len(n)) {
      ti <- df[[lt]][i]

      # System density via shared helper
      f_sys <- weibull_f_sys(ti, shapes, scales)
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
#' @examples
#' \donttest{
#' model <- kofn(k = 1, m = 2, family = "exponential")
#' set.seed(42)
#' df <- rdata_scheme1(model)(c(1, 2), n = 50, delta = 1.0)
#' result <- fit_scheme1(model)(df)
#' coef(result)
#' }
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

    multistart_mle(neg_ll, par0, n_par = n_par, n_starts = n_starts,
                   nobs = nrow(df))
  }
}
