#' Scheme 0: System-level only (no component information)
#'
#' Creates an observation functor for Scheme 0, where only the system
#' lifetime is observed. Optionally applies right-censoring at time
#' \code{tau}. This is the "black box" scenario: we see when the system
#' fails (or is censored) but learn nothing about individual components.
#'
#' Equivalent to \code{\link{observe_right_censor}} but named for the
#' scheme hierarchy used in the k-out-of-n censoring framework.
#'
#' @param tau Right-censoring time (default \code{Inf} = no censoring).
#' @return Observation functor: a function with signature
#'   \code{function(t_true)} returning a list with components:
#'   \describe{
#'     \item{t}{observed time (system failure time or \code{tau})}
#'     \item{omega}{\code{"exact"} or \code{"right"}}
#'     \item{t_upper}{\code{NA} (not used for this scheme)}
#'   }
#' @export
#' @examples
#' obs <- observe_scheme0(tau = 100)
#' obs(50)   # exact:  list(t = 50, omega = "exact", t_upper = NA)
#' obs(150)  # right:  list(t = 100, omega = "right", t_upper = NA)
#'
#' # No censoring (default)
#' obs_norc <- observe_scheme0()
#' obs_norc(999)  # exact: list(t = 999, omega = "exact", t_upper = NA)
observe_scheme0 <- function(tau = Inf) {
  force(tau)
  function(t_true) {
    if (t_true <= tau) {
      list(t = t_true, omega = "exact", t_upper = NA_real_)
    } else {
      list(t = tau, omega = "right", t_upper = NA_real_)
    }
  }
}


#' Scheme 1: Periodic inspection
#'
#' Creates an observation functor for Scheme 1, where the system is
#' inspected at regular intervals of width \code{delta}. The failure time
#' is known only to lie within the inspection interval containing it
#' (interval-censored). Systems surviving past \code{tau} are
#' right-censored.
#'
#' This corresponds to the periodic inspection design common in
#' reliability engineering, where continuous monitoring is infeasible
#' and inspections occur at scheduled times \code{delta, 2*delta, ...}.
#'
#' @param delta Inspection interval width (positive numeric).
#' @param tau Study end time / right-censoring time (positive numeric or
#'   \code{Inf} for no right-censoring).
#' @return Observation functor: a function with signature
#'   \code{function(t_true)} returning a list with components:
#'   \describe{
#'     \item{t}{lower bound of the inspection interval, or \code{tau} if
#'       right-censored}
#'     \item{omega}{\code{"interval"} or \code{"right"}}
#'     \item{t_upper}{upper bound of the inspection interval, or \code{NA}
#'       if right-censored}
#'   }
#' @export
#' @examples
#' obs <- observe_scheme1(delta = 10, tau = 100)
#' obs(25)   # interval: list(t = 20, omega = "interval", t_upper = 30)
#' obs(150)  # right:    list(t = 100, omega = "right", t_upper = NA)
#' obs(10)   # interval: list(t = 10, omega = "interval", t_upper = 20)
observe_scheme1 <- function(delta, tau = Inf) {
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


#' Scheme 2: Complete monitoring (exact observation)
#'
#' Creates an observation functor for Scheme 2, where all component
#' lifetimes are observed exactly. This is the trivial case where
#' complete data is available -- no censoring of any kind.
#'
#' Included for completeness in the scheme hierarchy and for use in
#' Fisher information comparisons across observation schemes.
#'
#' @return Observation functor: a function with signature
#'   \code{function(t_true)} returning a list with components:
#'   \describe{
#'     \item{t}{exact failure time}
#'     \item{omega}{\code{"exact"}}
#'     \item{t_upper}{\code{NA} (not used)}
#'   }
#' @export
#' @examples
#' obs <- observe_scheme2()
#' obs(42.5)  # exact: list(t = 42.5, omega = "exact", t_upper = NA)
observe_scheme2 <- function() {
  function(t_true) {
    list(t = t_true, omega = "exact", t_upper = NA_real_)
  }
}


#' Right-censoring observation scheme
#'
#' Creates an observation functor that applies right-censoring at time
#' \code{tau}. Systems that fail before \code{tau} are observed exactly;
#' systems surviving past \code{tau} are right-censored.
#'
#' This is a convenience alias: functionally identical to
#' \code{\link{observe_scheme0}}, but uses the traditional censoring
#' terminology rather than the scheme hierarchy.
#'
#' @param tau Censoring time (positive numeric).
#' @return Observation functor: a function with signature
#'   \code{function(t_true)} returning a list with components:
#'   \describe{
#'     \item{t}{observed time}
#'     \item{omega}{\code{"exact"} or \code{"right"}}
#'     \item{t_upper}{\code{NA} (not used)}
#'   }
#' @export
#' @examples
#' obs <- observe_right_censor(tau = 100)
#' obs(50)   # exact: list(t = 50, omega = "exact", t_upper = NA)
#' obs(150)  # right: list(t = 100, omega = "right", t_upper = NA)
observe_right_censor <- function(tau) {
  observe_scheme0(tau)
}


#' Periodic inspection observation scheme
#'
#' Creates an observation functor for periodic inspections. This is a
#' convenience alias for \code{\link{observe_scheme1}}, using the
#' traditional censoring terminology.
#'
#' @param delta Inspection interval width (positive numeric).
#' @param tau Study end time (positive numeric or \code{Inf}).
#' @return Observation functor (see \code{\link{observe_scheme1}} for
#'   details).
#' @export
#' @examples
#' obs <- observe_periodic(delta = 10, tau = 100)
#' obs(25)  # interval: list(t = 20, omega = "interval", t_upper = 30)
observe_periodic <- function(delta, tau = Inf) {
  observe_scheme1(delta, tau)
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
#' @return Observation functor: a function with signature
#'   \code{function(t_true)} that randomly selects a constituent scheme
#'   and returns its output.
#' @importFrom stats runif
#' @export
#' @examples
#' obs <- observe_mixture(
#'   observe_scheme0(tau = 100),
#'   observe_scheme1(delta = 10, tau = 100),
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
