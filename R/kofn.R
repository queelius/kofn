#' Create a k-out-of-n system estimation model
#'
#' Constructs a likelihood model for component lifetime estimation from
#' k-out-of-n system data. The system fails when \code{k} components have
#' failed: \code{k=1} is a series system (one failure kills it),
#' \code{k=m} is a parallel system (all must fail).
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
#' For non-k-of-n topologies (bridges, arbitrary coherent systems), use
#' \code{\link[dist.structure]{coherent_dist}} or one of the topology
#' shortcuts in \code{dist.structure} directly. kofn is exclusively for
#' the k-out-of-n family.
#'
#' @param k System parameter: system fails when k components have failed.
#'   \code{k=1} is series, \code{k=m} is parallel.
#' @param m Number of components.
#' @param family Component lifetime distribution: \code{"exponential"} or
#'   \code{"weibull"}.
#' @param method Estimation method: \code{"mle"} (direct MLE) or \code{"em"}
#'   (EM algorithm, Weibull only).
#' @param lifetime Column name for system lifetime (default \code{"t"}).
#' @param omega Column name for observation type (default \code{"omega"}).
#' @param lifetime_upper Column name for interval upper bound (default
#'   \code{"t_upper"}).
#' @return An S3 object of class \code{c("exp_kofn"/"wei_kofn", "kofn",
#'   "likelihood_model")}.
#' @export
#' @examples
#' # Parallel system with 3 exponential components (k = m)
#' model <- kofn(k = 3, m = 3, family = "exponential")
#' print(model)
#'
#' # Series system (k = 1)
#' model_series <- kofn(k = 1, m = 4, family = "exponential")
#'
#' # Weibull parallel system with EM estimation
#' model_wei <- kofn(k = 2, m = 2, family = "weibull", method = "em")
kofn <- function(k = 1L, m = 2L, family = "exponential",
                 method = "mle",
                 lifetime = "t", omega = "omega",
                 lifetime_upper = "t_upper") {
  family <- match.arg(family, c("exponential", "weibull"))
  method <- match.arg(method, c("mle", "em"))

  if (method == "em" && family == "exponential") {
    warning("EM method is designed for Weibull; using MLE for exponential")
    method <- "mle"
  }

  k <- as.integer(k)
  m <- as.integer(m)
  stopifnot(k >= 1L, k <= m, m >= 1L)

  model <- list(
    k = k,
    m = m,
    family = family,
    method = method,
    lifetime = lifetime,
    omega = omega,
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
#' print(kofn(k = 3, m = 3))
print.kofn <- function(x, ...) {
  sys_type <- if (x$k == 1L) {
    "series"
  } else if (x$k == x$m) {
    "parallel"
  } else {
    paste0(x$k, "-out-of-", x$m)
  }

  cat("k-out-of-n System Likelihood Model\n")
  cat("-----------------------------------\n")
  cat("  System type:", sys_type,
      sprintf("(k=%d, m=%d)\n", x$k, x$m))
  cat("  Component distribution:", x$family, "\n")
  cat("  Estimation method:", x$method, "\n")
  cat("  Column conventions:\n")
  cat("    lifetime:", x$lifetime, "\n")
  cat("    omega:", x$omega, "\n")
  cat("    interval upper:", x$lifetime_upper, "\n")
  invisible(x)
}


#' Number of components in a kofn model
#'
#' Returns the number of components \code{m} in the k-out-of-n system.
#'
#' @param model A \code{kofn} model object.
#' @param ... Additional arguments (ignored).
#' @return Integer number of components.
#' @method ncomponents kofn
#' @importFrom dist.structure ncomponents
#' @export
#' @examples
#' ncomponents(kofn(k = 5, m = 5))
ncomponents.kofn <- function(model, ...) {
  model$m
}


#' Test if an object is a kofn model
#'
#' @param x An object to test.
#' @return Logical indicating whether \code{x} inherits from \code{"kofn"}.
#' @export
#' @examples
#' is_kofn(kofn(k = 3, m = 3))
#' is_kofn(42)
is_kofn <- function(x) {
  inherits(x, "kofn")
}
