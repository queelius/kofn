# ===========================================================================
# Distribution Interface
# ===========================================================================
#
# Lightweight distribution objects for use in the system density and
# likelihood engines. Each distribution is a list with named functions:
#   pdf(t)  -- density
#   cdf(t)  -- CDF
#   surv(t) -- survival (= 1 - cdf)
#   name    -- character identifier
#   params  -- numeric vector of parameters
#
# All functions are vectorized over t.
# ===========================================================================


#' Create an exponential distribution object
#'
#' Returns a list representing an exponential distribution with the given
#' rate parameter. The list contains functions `pdf`, `cdf`, and `surv`
#' (all vectorized over `t`), plus `name` and `params` fields.
#'
#' @param rate Positive rate parameter \eqn{\lambda} (= 1/mean).
#' @return A list with elements:
#'   \describe{
#'     \item{pdf}{Function: \eqn{f(t) = \lambda e^{-\lambda t}}.}
#'     \item{cdf}{Function: \eqn{F(t) = 1 - e^{-\lambda t}}.}
#'     \item{surv}{Function: \eqn{S(t) = e^{-\lambda t}}.}
#'     \item{name}{Character: `"exponential"`.}
#'     \item{params}{Named numeric vector: `c(rate = rate)`.}
#'   }
#' @examples
#' d <- exp_dist(2)
#' d$pdf(1)   # 2 * exp(-2)
#' d$cdf(1)   # 1 - exp(-2)
#' d$surv(1)  # exp(-2)
#'
#' @export
exp_dist <- function(rate) {
  stopifnot(is.numeric(rate), length(rate) == 1L, rate > 0)
  list(
    pdf  = function(t) stats::dexp(t, rate = rate),
    cdf  = function(t) stats::pexp(t, rate = rate),
    surv = function(t) stats::pexp(t, rate = rate, lower.tail = FALSE),
    name   = "exponential",
    params = c(rate = rate)
  )
}


#' Create a Weibull distribution object
#'
#' Returns a list representing a Weibull distribution with the given shape
#' and scale parameters. The list contains functions `pdf`, `cdf`, and `surv`
#' (all vectorized over `t`), plus `name` and `params` fields.
#'
#' @param shape Positive shape parameter \eqn{k > 0}.
#' @param scale Positive scale parameter \eqn{\lambda > 0}.
#' @return A list with elements:
#'   \describe{
#'     \item{pdf}{Function: Weibull density.}
#'     \item{cdf}{Function: Weibull CDF.}
#'     \item{surv}{Function: Weibull survival.}
#'     \item{name}{Character: `"weibull"`.}
#'     \item{params}{Named numeric vector: `c(shape = shape, scale = scale)`.}
#'   }
#' @examples
#' d <- weibull_dist(shape = 2, scale = 1)
#' d$pdf(0.5)
#' d$cdf(0.5)
#'
#' @export
weibull_dist <- function(shape, scale) {
  stopifnot(is.numeric(shape), length(shape) == 1L, shape > 0)
  stopifnot(is.numeric(scale), length(scale) == 1L, scale > 0)
  list(
    pdf  = function(t) stats::dweibull(t, shape = shape, scale = scale),
    cdf  = function(t) stats::pweibull(t, shape = shape, scale = scale),
    surv = function(t) stats::pweibull(t, shape = shape, scale = scale,
                                       lower.tail = FALSE),
    name   = "weibull",
    params = c(shape = shape, scale = scale)
  )
}


#' Build distribution objects from a flat parameter vector
#'
#' Creates a list of m distribution objects (one per component) from a flat
#' numeric parameter vector.
#'
#' For `"exponential"`: `par` has m elements, one rate per component.
#' For `"weibull"`: `par` has 2m elements, interleaved as
#' `(shape_1, scale_1, shape_2, scale_2, ...)`.
#'
#' @param par Numeric parameter vector.
#' @param family Character: `"exponential"` or `"weibull"`.
#' @return List of m distribution objects (each as returned by [exp_dist()]
#'   or [weibull_dist()]).
#' @examples
#' dists <- make_dists(c(1, 2, 3), family = "exponential")
#' length(dists)  # 3
#' dists[[1]]$params  # c(rate = 1)
#'
#' @export
make_dists <- function(par, family = "exponential") {
  family <- match.arg(family, c("exponential", "weibull"))
  if (family == "exponential") {
    m <- length(par)
    lapply(seq_len(m), function(j) exp_dist(par[j]))
  } else if (family == "weibull") {
    if (length(par) %% 2 != 0) {
      stop("For Weibull family, par must have length 2m (shape_1, scale_1, ...)")
    }
    m <- length(par) / 2L
    lapply(seq_len(m), function(j) weibull_dist(par[2 * j - 1], par[2 * j]))
  }
}
