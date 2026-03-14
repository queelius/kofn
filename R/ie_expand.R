# ===========================================================================
# Inclusion-Exclusion Expansion for Exponential Parallel Systems
# ===========================================================================
#
# The core computational trick for exponential parallel systems: the product
#   prod_{i != j} (1 - exp(-lambda_i * t))
# is expanded via inclusion-exclusion into a signed sum of exponentials,
# making all integrals (left/interval censoring) analytical.
#
# Number of terms = 2^m, so only practical for moderate component counts
# (m <= ~15).
# ===========================================================================


#' Inclusion-exclusion expansion of a product of CDFs
#'
#' For a set of exponential rates `lam`, computes the inclusion-exclusion
#' expansion of \eqn{\prod_{i} (1 - e^{-\lambda_i t})}.
#'
#' Returns sign and rate-sum vectors such that:
#' \deqn{\prod_i (1 - e^{-\lambda_i t}) = \sum_k \text{sign}_k \cdot
#'   e^{-\text{rate\_sum}_k \cdot t}}
#'
#' The number of terms is \eqn{2^m} where \eqn{m} is `length(lam)`.
#'
#' @param lam Numeric vector of positive rates (one per component).
#' @return A list with elements:
#'   \describe{
#'     \item{sign}{Integer vector of signs (\eqn{+1} or \eqn{-1}).}
#'     \item{rate_sum}{Numeric vector of summed rates for each term.}
#'   }
#' @examples
#' ie <- ie_expand(c(1, 2))
#' # Product: (1 - exp(-t)) * (1 - exp(-2t))
#' # = 1 - exp(-t) - exp(-2t) + exp(-3t)
#' ie$sign      # c(1, -1, -1, 1)
#' ie$rate_sum  # c(0, 1, 2, 3)
#'
#' @export
ie_expand <- function(lam) {
  m <- length(lam)
  if (m == 0L) {
    return(list(sign = 1, rate_sum = 0))
  }

  n_terms <- 2L^m
  signs <- integer(n_terms)
  rate_sums <- numeric(n_terms)

  for (k in seq_len(n_terms) - 1L) {
    bits <- as.integer(intToBits(k)[seq_len(m)])
    signs[k + 1L] <- (-1L)^sum(bits)
    rate_sums[k + 1L] <- sum(lam[bits == 1L])
  }

  list(sign = signs, rate_sum = rate_sums)
}


#' Compute w_j(t) = f_j(t) * prod_{i != j} F_i(t) for exponential components
#'
#' For component j with rate \eqn{\lambda_j} and other rates
#' \eqn{\lambda_{-j}}:
#' \deqn{w_j(t) = \lambda_j e^{-\lambda_j t}
#'   \prod_{i \neq j} (1 - e^{-\lambda_i t})}
#'
#' Uses the inclusion-exclusion expansion to express this as a finite sum
#' of exponentials.
#'
#' @param t Scalar time point (positive numeric).
#' @param par Numeric vector of rates (length m), one per component.
#' @param j Component index (integer, 1-based).
#' @return Scalar value of \eqn{w_j(t)}.
#' @seealso [w_j_integral()] for the closed-form integral of \eqn{w_j}.
#' @export
w_j_exact <- function(t, par, j) {
  lam_j <- par[j]
  lam_others <- par[-j]
  ie <- ie_expand(lam_others)
  lam_j * sum(ie$sign * exp(-(lam_j + ie$rate_sum) * t))
}


#' Vectorized w_j over multiple time points
#'
#' Same as [w_j_exact()] but vectorized over `t` for a single component j.
#'
#' @param t Numeric vector of time points.
#' @param par Numeric vector of rates (length m).
#' @param j Component index (integer, 1-based).
#' @return Numeric vector of \eqn{w_j(t)} values.
#' @keywords internal
w_j_exact_vec <- function(t, par, j) {
  lam_j <- par[j]
  lam_others <- par[-j]
  ie <- ie_expand(lam_others)
  # Each row = one IE term, each col = one time point
  lam_j * colSums(ie$sign * exp(-outer(lam_j + ie$rate_sum, t)))
}


#' Closed-form integral of w_j(t) over an interval
#'
#' Computes \eqn{\int_a^b w_j(t)\, dt} in closed form using the
#' inclusion-exclusion expansion. Each term integrates as:
#' \deqn{\int_a^b e^{-r t}\, dt = \frac{e^{-ra} - e^{-rb}}{r}}
#'
#' @param a Lower bound of integration (non-negative numeric scalar).
#' @param b Upper bound of integration (positive numeric scalar, `b > a`).
#' @param par Numeric vector of rates (length m).
#' @param j Component index (integer, 1-based).
#' @return Scalar integral value.
#' @seealso [w_j_exact()] for the pointwise evaluation.
#' @export
w_j_integral <- function(a, b, par, j) {
  lam_j <- par[j]
  lam_others <- par[-j]
  ie <- ie_expand(lam_others)
  rates <- lam_j + ie$rate_sum

  terms <- ifelse(
    rates > 0,
    (exp(-rates * a) - exp(-rates * b)) / rates,
    b - a
  )
  lam_j * sum(ie$sign * terms)
}


#' System CDF for exponential parallel systems
#'
#' Computes \eqn{F_{sys}(t) = \prod_j (1 - e^{-\lambda_j t})} using the
#' inclusion-exclusion expansion.
#'
#' @param t Scalar time point (non-negative numeric).
#' @param par Numeric vector of rates (length m).
#' @return Scalar CDF value \eqn{P(T_{sys} \le t)}.
#' @seealso [S_sys_exp()] for the survival function.
#' @export
F_sys_exp <- function(t, par) {
  ie <- ie_expand(par)
  sum(ie$sign * exp(-ie$rate_sum * t))
}


#' System survival function for exponential parallel systems
#'
#' Computes \eqn{S_{sys}(t) = 1 - F_{sys}(t)} for a parallel system with
#' exponential components.
#'
#' @param t Scalar time point (non-negative numeric).
#' @param par Numeric vector of rates (length m).
#' @return Scalar survival probability \eqn{P(T_{sys} > t)}.
#' @seealso [F_sys_exp()] for the CDF.
#' @export
S_sys_exp <- function(t, par) {
  1 - F_sys_exp(t, par)
}
