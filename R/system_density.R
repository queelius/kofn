# ===========================================================================
# General System Density and Likelihood Engine
# ===========================================================================
#
# Distribution-agnostic MLE engine for coherent systems under Scheme 0
# (system-level observations only: T = system lifetime).
#
# The system density is:
#   f_sys(t) = sum_j f_j(t) * sum_{x_{-j}: j critical} prod_{k!=j} p_k(x_k, t)
#
# where p_k(1, t) = S_k(t) (component k survives to t) and
#       p_k(0, t) = F_k(t) (component k fails by t).
#
# For parallel systems with exponential components this reduces to the
# IE-based formula in ie_expand.R, enabling exact cross-validation.
# ===========================================================================


#' System survival function (general coherent system)
#'
#' Computes \eqn{P(T_{sys} > t)} by summing over all \eqn{2^m} state vectors:
#' \deqn{P(T_{sys} > t) = \sum_{x: \phi(x)=1}
#'   \prod_j S_j(t)^{x_j} F_j(t)^{1-x_j}}
#'
#' This is exact for any coherent system and any component distributions.
#' Computational cost is \eqn{O(2^m)} per time point.
#'
#' @param t Numeric vector of time points.
#' @param system A [coherent_system] object.
#' @param dists List of m distribution objects (one per component), each as
#'   returned by [exp_dist()] or [weibull_dist()].
#' @return Numeric vector of survival probabilities at each `t`.
#' @seealso [f_sys_general()] for the density.
#' @export
S_sys_general <- function(t, system, dists) {
  stopifnot(inherits(system, "coherent_system"))
  m <- system$m
  stopifnot(length(dists) == m)

  # Evaluate per-component survival and CDF at each time point
  S_mat <- matrix(0, nrow = length(t), ncol = m)
  F_mat <- matrix(0, nrow = length(t), ncol = m)
  for (j in seq_len(m)) {
    S_mat[, j] <- dists[[j]]$surv(t)
    F_mat[, j] <- dists[[j]]$cdf(t)
  }

  # Enumerate all 2^m state vectors
  n_states <- 2L^m
  result <- numeric(length(t))

  for (k in seq_len(n_states) - 1L) {
    # Decode state vector: bit j of k gives component j's state (1=up, 0=down)
    x <- as.logical(as.integer(intToBits(k))[seq_len(m)])
    if (!phi(system, x)) next

    # Product of S_j^{x_j} * F_j^{1-x_j}
    log_prob_mat <- matrix(0, nrow = length(t), ncol = m)
    for (j in seq_len(m)) {
      if (x[j]) {
        log_prob_mat[, j] <- log(pmax(S_mat[, j], .Machine$double.eps))
      } else {
        log_prob_mat[, j] <- log(pmax(F_mat[, j], .Machine$double.eps))
      }
    }
    result <- result + exp(rowSums(log_prob_mat))
  }
  result
}


#' System density function (general coherent system)
#'
#' Computes the system failure density:
#' \deqn{f_{sys}(t) = \sum_j f_j(t) \sum_{x_{-j}: j \text{ critical}}
#'   \prod_{k \neq j} p_k(x_k, t)}
#'
#' Uses [critical_states()] to identify pivotal states for each component.
#' Computational cost is \eqn{O(m \cdot 2^{m-1})} per time point.
#'
#' @param t Numeric vector of time points.
#' @param system A [coherent_system] object.
#' @param dists List of m distribution objects (one per component).
#' @return Numeric vector of density values at each `t`.
#' @seealso [S_sys_general()] for the survival function.
#' @export
f_sys_general <- function(t, system, dists) {
  stopifnot(inherits(system, "coherent_system"))
  m <- system$m
  stopifnot(length(dists) == m)
  t <- as.numeric(t)
  nt <- length(t)

  # Evaluate per-component density, survival, CDF at each time point
  f_mat <- matrix(0, nrow = nt, ncol = m)
  S_mat <- matrix(0, nrow = nt, ncol = m)
  F_mat <- matrix(0, nrow = nt, ncol = m)
  for (j in seq_len(m)) {
    f_mat[, j] <- dists[[j]]$pdf(t)
    S_mat[, j] <- dists[[j]]$surv(t)
    F_mat[, j] <- dists[[j]]$cdf(t)
  }

  result <- numeric(nt)

  # Precompute critical states for all components (depends only on system)
  crit_list <- lapply(seq_len(m), function(j) critical_states(system, j))
  others_list <- lapply(seq_len(m), function(j) setdiff(seq_len(m), j))

  for (j in seq_len(m)) {
    crit <- crit_list[[j]]
    if (nrow(crit) == 0L) next

    others <- others_list[[j]]

    # For each critical state, compute the product over k != j
    state_contrib <- matrix(0, nrow = nt, ncol = nrow(crit))

    for (r in seq_len(nrow(crit))) {
      x_others <- as.logical(crit[r, ])
      log_prod <- numeric(nt)
      for (idx in seq_along(others)) {
        k <- others[idx]
        if (x_others[idx]) {
          # Component k is up: contributes S_k(t)
          log_prod <- log_prod + log(pmax(S_mat[, k], .Machine$double.eps))
        } else {
          # Component k is down: contributes F_k(t)
          log_prod <- log_prod + log(pmax(F_mat[, k], .Machine$double.eps))
        }
      }
      state_contrib[, r] <- exp(log_prod)
    }

    # Sum over critical states, then multiply by f_j(t)
    result <- result + f_mat[, j] * rowSums(state_contrib)
  }

  result
}


#' Log-likelihood for Scheme 0 (system-level observations only)
#'
#' Computes the log-likelihood \eqn{\ell(\theta) = \sum_i \log f_{sys}(t_i;
#' \theta)} for observed system lifetimes under Scheme 0, where only the
#' system failure time is observed.
#'
#' @param t_obs Numeric vector of observed system lifetimes.
#' @param system A [coherent_system] object.
#' @param par Numeric parameter vector (rates for `"exponential"`, or
#'   interleaved shape/scale for `"weibull"`).
#' @param family Character: `"exponential"` or `"weibull"`.
#' @return Scalar log-likelihood value, or `-Inf` if parameters are invalid.
#' @seealso [fit_system()] for MLE, [rdata_system()] for data generation.
#' @export
loglik_system <- function(t_obs, system, par, family = "exponential") {
  if (any(!is.finite(par)) || any(par <= 0)) return(-Inf)
  dists <- tryCatch(
    make_dists(par, family),
    error = function(e) NULL
  )
  if (is.null(dists)) return(-Inf)

  fvals <- f_sys_general(t_obs, system, dists)
  if (any(!is.finite(fvals)) || any(fvals <= 0)) return(-Inf)

  sum(log(fvals))
}


#' Fit system MLE via multi-start optimization
#'
#' Maximum likelihood estimation of component lifetime parameters from
#' system-level failure time observations (Scheme 0). Uses L-BFGS-B as the
#' primary optimizer with a Nelder-Mead (log-scale) fallback. Standard errors
#' are computed from the observed Fisher information via `numDeriv::hessian`.
#'
#' @param t_obs Numeric vector of observed system lifetimes.
#' @param system A [coherent_system] object.
#' @param family Character: `"exponential"` or `"weibull"`.
#' @param init Optional numeric starting values (length = number of
#'   parameters). If `NULL`, method-of-moments initialization with spread is
#'   used.
#' @param n_starts Integer number of random restarts (default 5).
#' @return A \code{fisher_mle} object (from likelihood.model) with
#'   backward-compatible fields \code{convergence} (integer, 0 = success)
#'   and \code{se} (standard errors). Supports \code{coef()}, \code{vcov()},
#'   \code{logLik()}, \code{AIC()}, \code{confint()}, \code{summary()}.
#' @seealso [loglik_system()] for the log-likelihood,
#'   [rdata_system()] for data generation.
#' @export
fit_system <- function(t_obs, system, family = "exponential",
                       init = NULL, n_starts = 5L) {
  stopifnot(inherits(system, "coherent_system"))
  family <- match.arg(family, c("exponential", "weibull"))
  m <- system$m

  # Number of parameters per component
  n_par_per_comp <- switch(family, exponential = 1L, weibull = 2L)
  n_par <- m * n_par_per_comp

  # Default initialization: method-of-moments spread
  if (is.null(init)) {
    mean_t <- max(mean(t_obs, na.rm = TRUE), 0.01)
    if (family == "exponential") {
      H_m <- sum(1 / seq_len(m))
      lam0 <- H_m / mean_t
      init <- lam0 * seq(0.5, 1.5, length.out = m)
    } else {
      scale0 <- mean_t * seq(0.5, 1.5, length.out = m)
      init <- as.numeric(rbind(rep(1.0, m), scale0))  # interleaved
    }
  }
  stopifnot(length(init) == n_par)

  neg_ll <- function(par) {
    val <- -loglik_system(t_obs, system, par, family)
    if (!is.finite(val)) return(.Machine$double.xmax / 2)
    val
  }

  result <- multistart_mle(neg_ll, init, n_par = n_par, n_starts = n_starts,
                           nobs = length(t_obs))

  # Backward-compatible fields for vignette / direct users
  result$convergence <- if (result$converged) 0L else 1L
  result$se <- tryCatch(
    sqrt(diag(stats::vcov(result))),
    error = function(e) rep(NA_real_, n_par)
  )
  result
}


#' Generate system lifetime data from a coherent system
#'
#' Simulates component lifetimes from the specified distributions, computes
#' system lifetimes, and records per-component censoring status.
#'
#' @param system A [coherent_system] object.
#' @param par Numeric parameter vector (rates for `"exponential"`, or
#'   interleaved shape/scale for `"weibull"`).
#' @param family Character: `"exponential"` or `"weibull"`.
#' @param n Number of observations to generate (positive integer).
#' @return A data.frame with column `t` (system lifetimes). Additional
#'   information is stored as attributes:
#'   \describe{
#'     \item{comp_times}{Numeric matrix (n x m) of component lifetimes.}
#'     \item{critical}{Integer vector of critical component indices.}
#'     \item{status}{List of character vectors giving per-component censoring
#'       status (`"exact"`, `"left"`, or `"right"`).}
#'     \item{par}{The parameter vector used for generation.}
#'   }
#' @seealso [fit_system()] for estimation, [loglik_system()] for the
#'   log-likelihood.
#' @export
rdata_system <- function(system, par, family = "exponential", n = 100L) {
  stopifnot(inherits(system, "coherent_system"))
  family <- match.arg(family, c("exponential", "weibull"))
  m <- system$m
  dists <- make_dists(par, family)

  # Generate component lifetimes
  comp_times <- matrix(0, nrow = n, ncol = m)
  for (j in seq_len(m)) {
    p <- dists[[j]]$params
    if (family == "exponential") {
      comp_times[, j] <- stats::rexp(n, rate = p["rate"])
    } else {
      comp_times[, j] <- stats::rweibull(n, shape = p["shape"],
                                         scale = p["scale"])
    }
  }

  # System lifetimes and censoring patterns
  sys_times <- numeric(n)
  critical  <- integer(n)
  status    <- vector("list", n)

  for (i in seq_len(n)) {
    cens <- system_censoring(system, comp_times[i, ])
    sys_times[i]  <- cens$T_sys
    critical[i]   <- cens$critical
    status[[i]]   <- cens$status
  }

  df <- data.frame(t = sys_times)
  attr(df, "comp_times") <- comp_times
  attr(df, "critical")   <- critical
  attr(df, "status")     <- status
  attr(df, "par")        <- par
  df
}
