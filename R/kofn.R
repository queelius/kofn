#' Create a k-out-of-n system estimation model
#'
#' Constructs a likelihood model for component lifetime estimation from
#' k-out-of-n system data. The system fails when \code{k} components have
#' failed: \code{k=1} is a parallel system, \code{k=m} is a series system.
#'
#' This model satisfies the \code{likelihood_model} concept from the
#' \code{likelihood.model} package by providing methods for
#' \code{\link[likelihood.model]{loglik}},
#' \code{\link[likelihood.model]{score}}, and
#' \code{\link[likelihood.model]{hess_loglik}}.
#'
#' The class hierarchy is:
#' \itemize{
#'   \item \code{"exp_kofn"} or \code{"wei_kofn"} (distribution-specific dispatch)
#'   \item \code{"kofn"} (shared methods)
#'   \item \code{"likelihood_model"} (generic inference infrastructure)
#' }
#'
#' @param k System parameter: system fails when k components have failed.
#'   \code{k=1} is parallel, \code{k=m} is series.
#' @param m Number of components.
#' @param family Component lifetime distribution: \code{"exponential"} or
#'   \code{"weibull"}.
#' @param method Estimation method: \code{"mle"} (direct MLE) or \code{"em"}
#'   (EM algorithm, Weibull only).
#' @param system Optional: a \code{coherent_system} object (created by
#'   \code{\link{kofn_system}}). If provided, \code{k} and \code{m} are
#'   ignored and this system is used directly.
#' @param lifetime Column name for system lifetime (default \code{"t"}).
#' @param omega Column name for observation type (default \code{"omega"}).
#' @param candset Prefix for candidate set columns (default \code{"x"}).
#' @param lifetime_upper Column name for interval upper bound (default
#'   \code{"t_upper"}).
#' @return An S3 object of class \code{c("exp_kofn"/"wei_kofn", "kofn",
#'   "likelihood_model")}.
#' @export
#' @examples
#' # Parallel system with 3 exponential components
#' model <- kofn(k = 1, m = 3, family = "exponential")
#' print(model)
#'
#' # Series system (k = m)
#' model_series <- kofn(k = 4, m = 4, family = "exponential")
#'
#' # Weibull parallel system with EM estimation
#' model_wei <- kofn(k = 1, m = 2, family = "weibull", method = "em")
kofn <- function(k = 1L, m = 2L, family = "exponential",
                 method = "mle", system = NULL,
                 lifetime = "t", omega = "omega",
                 candset = "x", lifetime_upper = "t_upper") {
  family <- match.arg(family, c("exponential", "weibull"))
  method <- match.arg(method, c("mle", "em"))

  if (method == "em" && family == "exponential") {
    warning("EM method is designed for Weibull; using MLE for exponential")
    method <- "mle"
  }

  k <- as.integer(k)
  m <- as.integer(m)

  if (is.null(system)) {
    stopifnot(k >= 1L, k <= m, m >= 1L)
    system <- kofn_system(k, m)
  } else {
    m <- system$m
    # Detect k from system structure: k-out-of-n iff all paths have equal
    # size AND the number of paths equals choose(m, k)
    path_sizes <- vapply(system$min_paths, length, integer(1))
    if (length(unique(path_sizes)) == 1L) {
      detected_k <- as.integer(path_sizes[1L])
      if (length(system$min_paths) == choose(m, detected_k)) {
        k <- detected_k
      } else {
        k <- NA_integer_
      }
    } else {
      k <- NA_integer_
    }
  }

  model <- list(
    system = system,
    k = as.integer(k),
    m = as.integer(m),
    family = family,
    method = method,
    lifetime = lifetime,
    omega = omega,
    candset = candset,
    lifetime_upper = lifetime_upper
  )

  cls <- if (family == "exponential") "exp_kofn" else "wei_kofn"
  class(model) <- c(cls, "kofn", "likelihood_model")
  model
}


#' Print method for kofn models
#'
#' Displays a human-readable summary of the k-out-of-n system model,
#' including system type, component distribution, estimation method,
#' and column name conventions.
#'
#' @param x A \code{kofn} model object.
#' @param ... Additional arguments (ignored).
#' @return The model object, invisibly.
#' @method print kofn
#' @export
#' @examples
#' print(kofn(k = 1, m = 3))
print.kofn <- function(x, ...) {
  sys_type <- if (is.na(x$k)) {
    "general coherent"
  } else if (x$k == 1L) {
    "parallel"
  } else if (x$k == x$m) {
    "series"
  } else {
    paste0(x$k, "-out-of-", x$m)
  }

  cat("k-out-of-n System Likelihood Model\n")
  cat("-----------------------------------\n")
  k_str <- if (is.na(x$k)) "NA" else as.character(x$k)
  cat("  System type:", sys_type,
      sprintf("(k=%s, m=%d)\n", k_str, x$m))
  cat("  Component distribution:", x$family, "\n")
  cat("  Estimation method:", x$method, "\n")
  cat("  Column conventions:\n")
  cat("    lifetime:", x$lifetime, "\n")
  cat("    omega:", x$omega, "\n")
  cat("    candidate set:", paste0(x$candset, "1..", x$candset, x$m), "\n")
  cat("    interval upper:", x$lifetime_upper, "\n")
  invisible(x)
}


#' Number of components in a system model
#'
#' Generic function to return the number of components in a system model.
#'
#' @param model A model object.
#' @param ... Additional arguments (currently ignored).
#' @return Integer number of components.
#' @examples
#' ncomponents(kofn(k = 1, m = 4))  # 4
#'
#' @export
ncomponents <- function(model, ...) UseMethod("ncomponents")


#' @describeIn ncomponents Method for kofn models
#'
#' Returns the number of components \code{m} in the k-out-of-n system.
#'
#' @param model A \code{kofn} model object.
#' @param ... Additional arguments (ignored).
#' @return Integer number of components.
#' @method ncomponents kofn
#' @export
#' @examples
#' ncomponents(kofn(k = 1, m = 5))
ncomponents.kofn <- function(model, ...) {
  model$m
}


#' Test if an object is a kofn model
#'
#' @param x An object to test.
#' @return Logical indicating whether \code{x} inherits from \code{"kofn"}.
#' @export
#' @examples
#' is_kofn(kofn(k = 1, m = 3))
#' is_kofn(42)
is_kofn <- function(x) {
  inherits(x, "kofn")
}
