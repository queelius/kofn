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
  pat <- paste0("^", candset, "\\.?(\\d+)$")
  cols <- grep(pat, colnames(df), value = TRUE)
  if (length(cols) == 0L) {
    stop(sprintf("No candidate set columns found with prefix '%s'", candset))
  }
  rank <- as.integer(sub(pat, "\\1", cols))
  cmat <- as.matrix(df[cols[order(rank)]])
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


#' Extract model column name defaults
#'
#' Helper function to extract default column names from a kofn model
#' object. Used by all model methods to avoid repeating the same pattern.
#'
#' @param model A \code{kofn} model object with \code{lifetime},
#'   \code{lifetime_upper}, \code{omega}, \code{candset} fields.
#' @return A list with components \code{lifetime}, \code{lifetime_upper},
#'   \code{omega}, \code{candset}.
#' @keywords internal
extract_model_defaults <- function(model) {
  list(
    lifetime = model$lifetime %||% "t",
    lifetime_upper = model$lifetime_upper %||% "t_upper",
    omega = model$omega %||% "omega",
    candset = model$candset %||% "x"
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
