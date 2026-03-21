#' @importFrom stats pgamma
NULL

# ===========================================================================
# Truncated Weibull Moments
# ===========================================================================
#
# Truncated moments of the Weibull distribution, used by the EM algorithm
# for Weibull parallel system estimation.
#
# For T ~ Weibull(alpha, beta), the key identity is:
#   U = (T/beta)^alpha ~ Exp(1)
#
# This substitution reduces truncated Weibull moments to incomplete gamma
# functions evaluated at v_max = (t/beta)^alpha.
# ===========================================================================


#' Vectorized truncated power moment of the Weibull distribution
#'
#' Computes \eqn{E[T^k | T < t_i]} for \eqn{T \sim \text{Weibull}(\alpha, \beta)}
#' simultaneously for a vector of truncation points.
#'
#' Uses the substitution \eqn{U = (T/\beta)^\alpha \sim \text{Exp}(1)}, so that
#' \eqn{T^k = \beta^k U^{k/\alpha}} and the truncation \eqn{T < t} becomes
#' \eqn{U < v_{\max}}. The result is:
#' \deqn{E[T^k | T < t] = \beta^k \frac{\gamma(k/\alpha + 1, v_{\max})}{1 - e^{-v_{\max}}}}
#' where \eqn{\gamma(s, x)} is the lower incomplete gamma function.
#'
#' Computation is done in log-space for numerical stability.
#'
#' @param k Numeric scalar. Power parameter (can be non-integer).
#' @param v_max_vec Numeric vector. Precomputed values of \eqn{(t_i/\beta)^\alpha}.
#' @param alpha Numeric scalar. Weibull shape parameter.
#' @param beta Numeric scalar. Weibull scale parameter.
#' @return Numeric vector of \eqn{E[T^k | T < t_i]} values, same length as
#'   \code{v_max_vec}.
#'
#' @details
#' For very small \eqn{v_{\max}} (below \code{1e-12}), the asymptotic
#' approximation \eqn{E[T^k | T < t] \approx t^k / (k/\alpha + 1)} is used
#' to avoid numerical issues.
#'
#' @examples
#' # E[T^2 | T < t] for Weibull(shape=2, scale=1) at t = 0.5, 1.0, 2.0
#' v_max <- (c(0.5, 1.0, 2.0) / 1)^2
#' trunc_pow_moment_vec(k = 2, v_max_vec = v_max, alpha = 2, beta = 1)
#'
#' @export
trunc_pow_moment_vec <- function(k, v_max_vec, alpha, beta) {
  s <- k / alpha + 1
  n <- length(v_max_vec)
  result <- numeric(n)

  # Small v_max: use asymptotic E[T^k | T < t] -> t^k / (k/alpha + 1)
  small <- v_max_vec < 1e-12
  if (any(small)) {
    t_small <- beta * v_max_vec[small]^(1 / alpha)
    result[small] <- t_small^k / s
  }

  normal <- !small
  if (any(normal)) {
    vm <- v_max_vec[normal]
    # Log-space computation for stability:
    # log E[T^k|T<t] = k*log(beta) + lgamma(s) + log(pgamma(vm,s)) - log(1-exp(-vm))
    log_moment <- k * log(beta) + lgamma(s) +
      pgamma(vm, shape = s, lower.tail = TRUE, log.p = TRUE) -
      log1p(-exp(-vm))
    result[normal] <- exp(log_moment)
    # Guard: replace NaN/Inf with 0
    bad <- !is.finite(result[normal])
    if (any(bad)) {
      result[which(normal)[bad]] <- 0
    }
  }

  result
}


#' Scalar truncated power moment of the Weibull distribution
#'
#' Convenience wrapper around [trunc_pow_moment_vec()] for a single truncation
#' point. Computes \eqn{E[T^k | T < t]} for \eqn{T \sim \text{Weibull}(\alpha, \beta)}.
#'
#' @param k Numeric scalar. Power parameter (can be non-integer).
#' @param t Numeric scalar. Truncation point (positive).
#' @param alpha Numeric scalar. Weibull shape parameter.
#' @param beta Numeric scalar. Weibull scale parameter.
#' @return Numeric scalar \eqn{E[T^k | T < t]}.
#'
#' @examples
#' # E[T^1 | T < 2] for Weibull(shape=1.5, scale=1)
#' trunc_pow_moment(k = 1, t = 2, alpha = 1.5, beta = 1)
#'
#' @seealso [trunc_pow_moment_vec()] for the vectorized version.
#' @export
trunc_pow_moment <- function(k, t, alpha, beta) {
  trunc_pow_moment_vec(k, (t / beta)^alpha, alpha, beta)
}


#' Vectorized truncated log-moment of the Weibull distribution
#'
#' Computes \eqn{E[\log T | T < t_i]} for \eqn{T \sim \text{Weibull}(\alpha, \beta)}
#' simultaneously for a vector of truncation points.
#'
#' Uses the identity:
#' \deqn{E[\log T | T < t] = \log \beta + \frac{1}{\alpha} E[\log U | U < v_{\max}]}
#' where \eqn{U \sim \text{Exp}(1)} and \eqn{v_{\max} = (t/\beta)^\alpha}.
#'
#' The inner expectation \eqn{E[\log U | U < v_{\max}]} is computed via
#' central finite difference of the lower incomplete gamma function at \eqn{s = 1}.
#'
#' @param v_max_vec Numeric vector. Precomputed values of \eqn{(t_i/\beta)^\alpha}.
#' @param alpha Numeric scalar. Weibull shape parameter.
#' @param beta Numeric scalar. Weibull scale parameter.
#' @return Numeric vector of \eqn{E[\log T | T < t_i]} values, same length as
#'   \code{v_max_vec}.
#'
#' @details
#' For very small \eqn{v_{\max}} (below \code{1e-10}), the asymptotic
#' approximation \eqn{E[\log T | T < t] \approx \log t - 1/\alpha} is used.
#' For large \eqn{v_{\max}} where numerical issues arise, the untruncated
#' moment \eqn{E[\log T] = \log \beta + \psi(1)/\alpha} is used as fallback,
#' where \eqn{\psi} is the digamma function.
#'
#' @examples
#' # E[log T | T < t] for Weibull(shape=2, scale=1) at t = 0.5, 1.0, 2.0
#' v_max <- (c(0.5, 1.0, 2.0) / 1)^2
#' trunc_log_moment_vec(v_max_vec = v_max, alpha = 2, beta = 1)
#'
#' @export
trunc_log_moment_vec <- function(v_max_vec, alpha, beta) {
  n <- length(v_max_vec)
  result <- numeric(n)

  # Small v_max: asymptotic E[log T | T < t] -> log(t) - 1/alpha
  small <- v_max_vec < 1e-10
  if (any(small)) {
    t_small <- beta * v_max_vec[small]^(1 / alpha)
    result[small] <- log(t_small) - 1 / alpha
  }

  normal <- !small
  if (any(normal)) {
    vm <- v_max_vec[normal]

    # Central finite difference: d/ds[gamma(s, x)]|_{s=1}
    # where gamma(s, x) = Gamma(s) * pgamma(x, s)
    h <- 1e-6
    s_plus <- 1 + h
    s_minus <- 1 - h
    gamma_plus <- gamma(s_plus) * pgamma(vm, shape = s_plus)
    gamma_minus <- gamma(s_minus) * pgamma(vm, shape = s_minus)
    d_gamma <- (gamma_plus - gamma_minus) / (2 * h)

    denom <- 1 - exp(-vm)
    elog_u <- d_gamma / denom  # E[log U | U < v_max]

    result[normal] <- log(beta) + elog_u / alpha

    # Guard NaN/Inf: use untruncated moment as fallback
    bad <- !is.finite(result[normal])
    if (any(bad)) {
      # Untruncated: E[log T] = log(beta) + digamma(1)/alpha
      result[which(normal)[bad]] <- log(beta) + digamma(1) / alpha
    }
  }

  result
}
