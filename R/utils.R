# Null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x


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


#' Extract and validate data from a system observation data frame
#'
#' Shared validation and extraction logic used by all likelihood model
#' methods. Checks that the data frame is non-empty, required columns
#' exist, decodes the candidate set matrix from prefixed Boolean columns,
#' and validates observation types.
#'
#' @param df Data frame containing system observations.
#' @param lifetime Column name for system lifetime (character).
#' @param omega Column name for observation type (character). Must contain
#'   values from \code{c("exact", "right", "left", "interval")}.
#' @param candset Column prefix for candidate set indicators (character).
#'   Columns matching \code{candset1, candset2, ...} or
#'   \code{candset.1, candset.2, ...} are extracted as a Boolean matrix.
#' @param lifetime_upper Column name for interval upper bound (character or
#'   \code{NULL}). Required when interval-censored observations are present.
#' @return A list with components:
#'   \describe{
#'     \item{t}{numeric vector of observed lifetimes}
#'     \item{omega}{character vector of observation types}
#'     \item{C}{logical matrix of candidate sets (n x m)}
#'     \item{m}{integer number of components}
#'     \item{n}{integer number of observations}
#'     \item{t_upper}{numeric vector of interval upper bounds, or \code{NULL}}
#'   }
#' @export
#' @examples
#' df <- data.frame(
#'   t = c(1.2, 3.4, 5.6),
#'   omega = c("exact", "right", "exact"),
#'   x1 = c(TRUE, FALSE, TRUE),
#'   x2 = c(TRUE, FALSE, FALSE),
#'   x3 = c(FALSE, FALSE, TRUE)
#' )
#' d <- extract_data(df, "t", "omega", "x")
#' d$m   # 3
#' d$C   # 3x3 logical matrix
extract_data <- function(df, lifetime, omega, candset,
                         lifetime_upper = NULL) {
  n <- nrow(df)
  if (n == 0L) stop("Data frame is empty")

  if (!lifetime %in% names(df)) {
    stop(sprintf("Column '%s' not found in data frame", lifetime))
  }
  if (!omega %in% names(df)) {
    stop(sprintf("Column '%s' not found in data frame", omega))
  }

  # Extract candidate set matrix from prefixed columns
  cmat <- md_decode_matrix(df, candset)
  if (is.null(cmat)) {
    stop(sprintf("No candidate set columns found with prefix '%s'", candset))
  }
  m <- ncol(cmat)

  # Validate observation types
  omega_vals <- as.character(df[[omega]])
  valid_omega <- c("exact", "right", "left", "interval")
  bad <- setdiff(unique(omega_vals), valid_omega)
  if (length(bad) > 0L) {
    stop(sprintf("Invalid omega values: %s. Must be one of: %s",
                 paste(bad, collapse = ", "),
                 paste(valid_omega, collapse = ", ")))
  }

  # Extract interval upper bounds if needed
  t_upper <- NULL
  if (any(omega_vals == "interval")) {
    if (is.null(lifetime_upper) || !lifetime_upper %in% names(df)) {
      stop("Interval-censored observations require a '",
           lifetime_upper %||% "t_upper", "' column")
    }
    t_upper <- df[[lifetime_upper]]
  }

  list(
    t = df[[lifetime]],
    omega = omega_vals,
    C = cmat,
    m = m,
    n = n,
    t_upper = t_upper
  )
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


#' Decode a matrix from prefixed columns in a data frame
#'
#' Extracts columns matching the pattern \code{var1, var2, ...} or
#' \code{var.1, var.2, ...} from a data frame and returns them as a
#' numeric matrix, ordered by index.
#'
#' @param df Data frame containing the matrix columns.
#' @param var Character prefix for the column names.
#' @return A matrix, or \code{NULL} if no matching columns are found.
#' @keywords internal
md_decode_matrix <- function(df, var) {
  stopifnot(is.data.frame(df), is.character(var))
  pat <- paste0("^", var, "\\.?(\\d+)$")
  cols <- grep(pat, colnames(df), value = TRUE)
  if (length(cols) == 0L) return(NULL)
  rank <- as.integer(sub(pat, "\\1", cols))
  as.matrix(df[cols[order(rank)]])
}
