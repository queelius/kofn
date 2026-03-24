#' Construct an MLE result object
#'
#' Creates a \code{\link[likelihood.model]{fisher_mle}} object from
#' optimization results. This gives kofn results full base R stats
#' compatibility: \code{coef()}, \code{vcov()}, \code{logLik()},
#' \code{AIC()}, \code{BIC()}, \code{confint()}, \code{summary()}.
#'
#' @param par Numeric vector of parameter estimates.
#' @param loglik Scalar log-likelihood at the MLE.
#' @param hessian Hessian matrix of the negative log-likelihood at MLE.
#' @param score Score vector at the MLE (should be near zero).
#' @param nobs Number of observations.
#' @param converged Logical indicating convergence.
#' @param ... Additional named elements stored as attributes on the result.
#' @return A \code{fisher_mle} object (from likelihood.model) with
#'   any extra fields (shapes, scales, iterations) attached.
#' @keywords internal
make_mle_result <- function(par, loglik, hessian, score, nobs,
                            converged, ...) {
  # Callers pass the Hessian of the NEGATIVE log-likelihood (observed
  # Fisher info, positive definite). fisher_mle expects the Hessian of
  # the log-likelihood itself, so we negate.
  hess_loglik <- if (all(is.finite(hessian))) -hessian else NULL
  result <- likelihood.model::fisher_mle(
    par        = par,
    loglik_val = loglik,
    hessian    = hess_loglik,
    score_val  = score,
    nobs       = nobs,
    converged  = converged
  )

  # Attach extra fields (shapes, scales, iterations, etc.)
  extras <- list(...)
  for (nm in names(extras)) {
    result[[nm]] <- extras[[nm]]
  }

  result
}


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


#' Numerical Hessian and score at a candidate MLE
#'
#' Computes the Hessian of the negative log-likelihood and the score
#' (gradient of the log-likelihood) at a parameter estimate using
#' \code{numDeriv}. Shared by \code{multistart_mle} and the EM solver.
#'
#' @param neg_ll Function of \code{par} only: negative log-likelihood.
#' @param par Numeric vector: parameter estimate.
#' @param n_par Integer: number of parameters.
#' @return A list with \code{hessian} (matrix of the negative log-likelihood)
#'   and \code{score} (gradient of the log-likelihood).
#' @keywords internal
hessian_score_at_mle <- function(neg_ll, par, n_par) {
  H <- tryCatch(
    numDeriv::hessian(neg_ll, x = par, method = "Richardson"),
    error = function(e) NULL
  )
  hessian_mat <- if (!is.null(H) && all(is.finite(H))) H else {
    matrix(NA_real_, nrow = n_par, ncol = n_par)
  }
  score_val <- tryCatch(
    -numDeriv::grad(func = neg_ll, x = par, method = "Richardson"),
    error = function(e) rep(NA_real_, n_par)
  )
  list(hessian = hessian_mat, score = score_val)
}


#' Multi-start MLE optimizer
#'
#' Common optimization infrastructure used by all fitting methods.
#' L-BFGS-B primary with Nelder-Mead log-scale fallback, multi-start
#' with log-normal perturbation. Computes Hessian and score at the MLE
#' via numDeriv and returns a \code{fisher_mle} object.
#'
#' @param neg_ll Function of \code{par} only: guarded negative
#'   log-likelihood (returns \code{.Machine$double.xmax / 2} on failure).
#' @param par0 Numeric vector of initial parameter values.
#' @param n_par Integer number of parameters.
#' @param n_starts Integer number of random restarts (default 5).
#' @param nobs Integer number of observations (for \code{fisher_mle}).
#' @param ... Additional named elements attached to the result
#'   (e.g. \code{shapes}, \code{scales}).
#' @return A \code{fisher_mle} object (from likelihood.model).
#' @keywords internal
multistart_mle <- function(neg_ll, par0, n_par, n_starts = 5L,
                           nobs = NULL, ...) {
  fit_from_init <- function(init_par) {
    result <- tryCatch(
      stats::optim(par = init_par, fn = neg_ll, method = "L-BFGS-B",
                   lower = rep(1e-10, n_par), upper = rep(Inf, n_par)),
      error = function(e) {
        list(par = rep(NA_real_, n_par), value = Inf, convergence = 99L)
      }
    )
    if (result$convergence != 0) {
      neg_ll_log <- function(log_p) neg_ll(exp(log_p))
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

  # Guard: if best value is at the penalty ceiling, the optimizer
  # "converged" to a degenerate solution — treat as non-convergence
  if (best$value >= .Machine$double.xmax / 4) {
    best$convergence <- 99L
  }

  result <- hessian_score_at_mle(neg_ll, best$par, n_par)

  make_mle_result(
    par       = best$par,
    loglik    = -best$value,
    hessian   = result$hessian,
    score     = result$score,
    nobs      = nobs,
    converged = (best$convergence == 0),
    ...
  )
}


