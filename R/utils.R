#' @importFrom stats coef
NULL

#' Parse a flat parameter vector into shapes and scales
#'
#' Unified extraction of Weibull shape/scale parameters from the flat
#' vectors used by optimizers. For exponential family, shapes are all 1
#' and scales are 1/rate.
#'
#' @param par Numeric parameter vector.
#' @param m Number of components.
#' @param family Character: \code{"exponential"} or \code{"weibull"}.
#' @return A list with components \code{shapes} and \code{scales}.
#' @keywords internal
parse_params <- function(par, m, family) {
  if (family == "exponential") {
    list(shapes = rep(1, m), scales = 1 / par)
  } else {
    list(
      shapes = par[seq(1L, 2L * m, by = 2L)],
      scales = par[seq(2L, 2L * m, by = 2L)]
    )
  }
}


#' Default MLE solver for positive parameters
#'
#' Multi-start optimization via \pkg{compositional.mle}: L-BFGS-B
#' with positivity constraints as the primary solver, Nelder-Mead on
#' the log-parameter scale as fallback. Runs from \code{n_starts}
#' random perturbations of \code{par0} and returns the best result.
#'
#' @param neg_ll Function of \code{par} only: guarded negative
#'   log-likelihood (returns \code{.Machine$double.xmax / 2} on failure).
#' @param par0 Numeric vector of initial parameter values.
#' @param n_par Integer number of parameters.
#' @param n_starts Integer number of random restarts (default 5).
#' @param nobs Integer number of observations.
#' @return An \code{mle_numerical} result object (from algebraic.mle)
#'   with \code{coef()}, \code{vcov()}, \code{logLik()}, etc. Returns
#'   \code{NULL} if all starts fail.
#' @export
#' @examples
#' # Fit a simple exponential rate from positive data
#' x <- rexp(100, rate = 0.5)
#' neg_ll <- function(rate) -sum(dexp(x, rate, log = TRUE))
#' fit <- solve_mle(neg_ll, par0 = 1, n_par = 1, nobs = length(x))
#' coef(fit)
solve_mle <- function(neg_ll, par0, n_par, n_starts = 5L, nobs = NULL) {
  loglike <- function(par) -neg_ll(par)
  prob <- compositional.mle::mle_problem(
    loglike = loglike,
    n_obs   = nobs,
    constraint = compositional.mle::mle_constraint(
      support = function(par) all(par > 0),
      project = function(par) pmax(par, 1e-10)
    )
  )
  primary <- compositional.mle::lbfgsb(lower = rep(1e-10, n_par))
  fallback <- compositional.mle::nelder_mead()

  # Log-scale problem for Nelder-Mead fallback (ensures positivity)
  loglike_log <- function(log_par) loglike(exp(log_par))
  prob_log <- compositional.mle::mle_problem(loglike = loglike_log, n_obs = nobs)

  fit_once <- function(init) {
    res <- tryCatch(primary(prob, init), error = function(e) NULL)
    if (is.null(res) || !res$converged) {
      res2 <- tryCatch(fallback(prob_log, log(pmax(init, 1e-10))),
                        error = function(e) NULL)
      if (!is.null(res2)) {
        # Transform back to original scale
        res2$theta.hat <- exp(res2$theta.hat)
        res2$loglike <- loglike(res2$theta.hat)
        if (is.null(res) || res2$loglike > res$loglike) {
          res <- res2
        }
      }
    }
    res
  }

  best <- fit_once(par0)
  if (n_starts > 1L) {
    for (s in seq_len(n_starts - 1L)) {
      init_s <- pmax(par0 * exp(stats::rnorm(n_par, sd = 0.5)), 1e-10)
      res_s <- fit_once(init_s)
      if (!is.null(res_s) && (is.null(best) || res_s$loglike > best$loglike)) {
        best <- res_s
      }
    }
  }

  # Ensure nobs is populated for downstream BIC()/AIC() calls; the
  # compositional.mle solver does not always propagate it from mle_problem().
  if (!is.null(best) && is.null(best$nobs) && !is.null(nobs)) {
    best$nobs <- nobs
  }

  best
}
