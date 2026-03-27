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
#' @examples
#' sys <- parallel_system(2)
#' dists <- list(exp_dist(1), exp_dist(2))
#' S_sys_general(c(0.5, 1.0, 2.0), sys, dists)
#'
#' @seealso [f_sys_general()] for the density.
#' @export
S_sys_general <- function(t, system, dists) {
  stopifnot(inherits(system, "coherent_system"))
  m <- system$m
  stopifnot(length(dists) == m)
  t <- as.numeric(t)
  nt <- length(t)

  # Evaluate per-component log-survival and log-CDF as matrices
  log_S_mat <- matrix(0, nrow = nt, ncol = m)
  log_F_mat <- matrix(0, nrow = nt, ncol = m)
  for (j in seq_len(m)) {
    log_S_mat[, j] <- log(pmax(dists[[j]]$surv(t), .Machine$double.eps))
    log_F_mat[, j] <- log(pmax(dists[[j]]$cdf(t), .Machine$double.eps))
  }

  # Use eagerly precomputed functioning states from system$cache
  func_states <- system$cache$func

  result <- numeric(nt)
  for (s in func_states) {
    log_prod <- numeric(nt)
    if (length(s$up) > 0L) {
      log_prod <- log_prod + rowSums(log_S_mat[, s$up, drop = FALSE])
    }
    if (length(s$down) > 0L) {
      log_prod <- log_prod + rowSums(log_F_mat[, s$down, drop = FALSE])
    }
    result <- result + exp(log_prod)
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
#' Critical states are cached on the system object after first computation.
#' Computational cost is \eqn{O(m \cdot 2^{m-1})} per time point for the
#' first call; subsequent calls skip the state enumeration.
#'
#' @param t Numeric vector of time points.
#' @param system A [coherent_system] object.
#' @param dists List of m distribution objects (one per component).
#' @return Numeric vector of density values at each `t`.
#' @examples
#' sys <- parallel_system(2)
#' dists <- list(exp_dist(1), exp_dist(2))
#' f_sys_general(c(0.5, 1.0, 2.0), sys, dists)
#'
#' @seealso [S_sys_general()] for the survival function.
#' @export
f_sys_general <- function(t, system, dists) {
  stopifnot(inherits(system, "coherent_system"))
  m <- system$m
  stopifnot(length(dists) == m)
  t <- as.numeric(t)
  nt <- length(t)

  # Evaluate per-component density, survival, CDF as matrices (vectorized)
  f_mat <- matrix(0, nrow = nt, ncol = m)
  log_S_mat <- matrix(0, nrow = nt, ncol = m)
  log_F_mat <- matrix(0, nrow = nt, ncol = m)
  for (j in seq_len(m)) {
    f_mat[, j] <- dists[[j]]$pdf(t)
    S_j <- dists[[j]]$surv(t)
    F_j <- dists[[j]]$cdf(t)
    log_S_mat[, j] <- log(pmax(S_j, .Machine$double.eps))
    log_F_mat[, j] <- log(pmax(F_j, .Machine$double.eps))
  }

  # Use eagerly precomputed critical states from system$cache (environment, reference semantics)
  crit_cache <- system$cache$crit

  result <- numeric(nt)

  for (j in seq_len(m)) {
    cc <- crit_cache[[j]]
    if (is.null(cc$n_states) || cc$n_states == 0L) next

    # For each critical state, compute product via log-sum of precomputed log matrices
    state_contrib <- numeric(nt)
    for (r in seq_len(cc$n_states)) {
      log_prod <- numeric(nt)
      up_idx <- cc$up[[r]]
      down_idx <- cc$down[[r]]
      if (length(up_idx) > 0L) {
        log_prod <- log_prod + rowSums(log_S_mat[, up_idx, drop = FALSE])
      }
      if (length(down_idx) > 0L) {
        log_prod <- log_prod + rowSums(log_F_mat[, down_idx, drop = FALSE])
      }
      state_contrib <- state_contrib + exp(log_prod)
    }

    result <- result + f_mat[, j] * state_contrib
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
#' @examples
#' sys <- parallel_system(2)
#' set.seed(1)
#' df <- rdata_system(sys, par = c(1, 2), n = 30)
#' loglik_system(df$t, sys, par = c(1, 2))
#'
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
#' @examples
#' \donttest{
#' sys <- parallel_system(2)
#' set.seed(42)
#' df <- rdata_system(sys, par = c(1, 2), n = 50)
#' result <- fit_system(df$t, sys)
#' coef(result)
#' }
#'
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
  result$se <- tryCatch({
    V <- stats::vcov(result)
    if (is.null(V)) stop("NULL vcov")
    sqrt(diag(V))
  }, error = function(e) rep(NA_real_, n_par))
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
#' @examples
#' sys <- parallel_system(3)
#' set.seed(1)
#' df <- rdata_system(sys, par = c(1, 2, 3), n = 20)
#' head(df)
#' attr(df, "critical")  # which component was critical
#'
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
