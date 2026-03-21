# ===========================================================================
# Observation Functors
# ===========================================================================
#
# Observation functors model what an experimenter actually records about
# a system failure. Each functor is a closure: function(t_true) -> list(t,
# omega, t_upper). The omega field encodes the observation type:
#
#   "exact"    -- failure time observed exactly
#   "right"    -- system survived past observed time (T > t)
#   "interval" -- failure in [t, t_upper)
#   "left"     -- failure before observed time (T < t)
#
# These compose: pass any functor as the `observe` argument to rdata().
# ===========================================================================


#' Right-censoring observation scheme
#'
#' Creates an observation functor that applies right-censoring at time
#' \code{tau}. Systems that fail before \code{tau} are observed exactly;
#' systems surviving past \code{tau} are right-censored.
#'
#' @param tau Censoring time (positive numeric, or \code{Inf} for no
#'   censoring).
#' @return Observation functor: \code{function(t_true)} returning a list
#'   with \code{t}, \code{omega} (\code{"exact"} or \code{"right"}), and
#'   \code{t_upper} (\code{NA}).
#' @export
#' @examples
#' obs <- observe_right_censor(tau = 100)
#' obs(50)   # exact:  list(t = 50, omega = "exact", t_upper = NA)
#' obs(150)  # right:  list(t = 100, omega = "right", t_upper = NA)
#'
#' # No censoring
#' obs_all <- observe_right_censor(tau = Inf)
#' obs_all(999)  # exact: list(t = 999, omega = "exact", t_upper = NA)
observe_right_censor <- function(tau = Inf) {
  force(tau)
  function(t_true) {
    if (t_true <= tau) {
      list(t = t_true, omega = "exact", t_upper = NA_real_)
    } else {
      list(t = tau, omega = "right", t_upper = NA_real_)
    }
  }
}


#' Periodic inspection observation scheme
#'
#' Creates an observation functor for periodic inspections at intervals
#' of width \code{delta}. The failure time is known only to lie within
#' the inspection interval containing it (interval-censored). Systems
#' surviving past \code{tau} are right-censored.
#'
#' @param delta Inspection interval width (positive numeric).
#' @param tau Study end time / right-censoring time (positive numeric or
#'   \code{Inf} for no right-censoring).
#' @return Observation functor: \code{function(t_true)} returning a list
#'   with \code{t} (interval lower bound or \code{tau}),
#'   \code{omega} (\code{"interval"} or \code{"right"}), and
#'   \code{t_upper} (interval upper bound or \code{NA}).
#' @export
#' @examples
#' obs <- observe_periodic(delta = 10, tau = 100)
#' obs(25)   # interval: list(t = 20, omega = "interval", t_upper = 30)
#' obs(150)  # right:    list(t = 100, omega = "right", t_upper = NA)
observe_periodic <- function(delta, tau = Inf) {
  force(delta)
  force(tau)
  function(t_true) {
    if (t_true > tau) {
      return(list(t = tau, omega = "right", t_upper = NA_real_))
    }
    lower <- floor(t_true / delta) * delta
    upper <- lower + delta
    list(t = lower, omega = "interval", t_upper = upper)
  }
}


#' Exact observation scheme (no censoring)
#'
#' Creates an observation functor that records the exact failure time
#' with no censoring. This is the trivial case — included for
#' completeness and for explicit use in Fisher information comparisons.
#'
#' @return Observation functor: \code{function(t_true)} returning a list
#'   with \code{t} (exact time), \code{omega} (\code{"exact"}), and
#'   \code{t_upper} (\code{NA}).
#' @export
#' @examples
#' obs <- observe_exact()
#' obs(42.5)  # list(t = 42.5, omega = "exact", t_upper = NA)
observe_exact <- function() {
  function(t_true) {
    list(t = t_true, omega = "exact", t_upper = NA_real_)
  }
}


#' Mixture of observation schemes
#'
#' Creates an observation functor that randomly selects from a set of
#' observation schemes for each observation. This models heterogeneous
#' monitoring environments where different units may be observed under
#' different protocols.
#'
#' @param ... Observation functors (created by \code{observe_*} functions).
#' @param weights Mixing probabilities (numeric vector). If \code{NULL},
#'   uniform weights are used. Weights are normalized to sum to 1.
#' @return Observation functor: \code{function(t_true)} that randomly
#'   selects a constituent scheme and returns its output.
#' @importFrom stats runif
#' @export
#' @examples
#' obs <- observe_mixture(
#'   observe_right_censor(tau = 100),
#'   observe_periodic(delta = 10, tau = 100),
#'   weights = c(0.7, 0.3)
#' )
#' set.seed(42)
#' obs(30)
observe_mixture <- function(..., weights = NULL) {
  schemes <- list(...)
  if (length(schemes) == 0L) {
    stop("at least one observation scheme is required")
  }
  if (is.null(weights)) {
    weights <- rep(1, length(schemes))
  }
  if (length(weights) != length(schemes)) {
    stop("weights must have the same length as the number of schemes")
  }
  weights <- weights / sum(weights)

  function(t_true) {
    idx <- sample.int(length(schemes), 1L, prob = weights)
    schemes[[idx]](t_true)
  }
}


