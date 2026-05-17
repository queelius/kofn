# ===========================================================================
# Generalized Masked k-out-of-n Likelihood
# ===========================================================================
#
# Extends the masked cause-of-failure framework from series systems (k=1,
# handled by maskedcauses) to general k-out-of-n systems.
#
# At system failure T_sys = T_{(k)}, exactly k components have failed.
# A candidate failed set C (with |C| >= k) is a superset of the true
# failed set F (with |F| = k). The likelihood marginalizes over all
# possible k-subsets of C:
#
#   L_i(theta) = sum_{F in C(C_i, k)} sum_{j in F}
#     f_j(t_i) * prod_{l in F\{j}} F_l(t_i) * prod_{l not in F} S_l(t_i)
#
# For k=1 (series), this reduces to sum_{j in C} w_j(t) -- the standard
# maskedcauses formula.
#
# For k=m (parallel), F = {1,...,m} is forced, and the sum over j gives
# f_sys(t) -- the system density with no masking.
# ===========================================================================


#' Masked k-out-of-n log-likelihood
#'
#' Computes the log-likelihood for a k-out-of-n system with candidate
#' failed sets (generalized masking). The data frame must contain columns
#' for system lifetime, observation type, and Boolean candidate set
#' indicators (\code{c1, c2, ..., cm}).
#'
#' At system failure, k components have failed. The candidate set
#' \code{C_i} is a superset of the true failed set. The likelihood
#' marginalizes over all \eqn{\binom{|C_i|}{k}}{choose(|C_i|, k)}
#' possible k-subsets.
#'
#' @param model A \code{kofn} model object.
#' @param candset Prefix for candidate set columns (default \code{"c"}).
#' @param ... Additional arguments (currently unused).
#' @return A function \code{function(df, par)} returning a scalar
#'   log-likelihood.
#'
#' @details
#' For series systems (\code{k = 1}), this reduces to the standard
#' masked cause-of-failure likelihood \eqn{\sum_{j \in C_i} w_j(t_i)}.
#' For parallel systems (\code{k = m}), the candidate set must be
#' \code{\{1, \ldots, m\}} and the likelihood equals the system density.
#'
#' The function uses the exponential or Weibull distribution functions
#' directly (not the IE expansion), so it works for any family supported
#' by \code{\link{parse_params}}.
#'
#' @examples
#' # 2-out-of-4 system with candidate failed sets
#' model <- kofn(k = 2, m = 4, component = dfr_exponential())
#' ll <- loglik_masked(model)
#'
#' # Manually construct masked data: k=2 failed, C = {1,2,3}
#' df <- data.frame(
#'   t = 1.5, omega = "exact",
#'   c1 = TRUE, c2 = TRUE, c3 = TRUE, c4 = FALSE
#' )
#' ll(df, c(1, 0.8, 0.6, 0.4))
#'
#' @export
loglik_masked <- function(model, candset = "c", ...) {
  m <- model$m
  k <- model$k
  lt <- model$lifetime
  om <- model$omega
  component <- model$component

  if (is.na(k)) stop("Masked likelihood requires a k-out-of-n system (k must not be NA)")

  cand_cols <- paste0(candset, seq_len(m))

  function(df, par) {
    if (any(!is.finite(par)) || any(par <= 0)) return(-Inf)

    pp <- parse_params(par, m, component)
    shapes <- pp$shapes
    scales <- pp$scales

    t_obs <- df[[lt]]
    omega_vals <- as.character(df[[om]])
    n <- length(t_obs)

    # Extract candidate set matrix
    cmat <- as.matrix(df[cand_cols])

    # Lazily construct dist.structure DGP for non-exact contributions.
    dgp <- NULL
    surv_fn <- NULL

    ll <- 0
    for (i in seq_len(n)) {
      if (omega_vals[i] != "exact") {
        # For non-exact observations, fall back to marginal system survival
        # (masking info not applicable for censored obs)
        if (omega_vals[i] == "right") {
          if (is.null(dgp)) {
            dgp <- kofn_dgp(k, m, component, par)
            surv_fn <- algebraic.dist::surv(dgp)
          }
          S_sys <- surv_fn(t_obs[i])
          if (S_sys <= 0) return(-Inf)
          ll <- ll + log(S_sys)
        }
        next
      }

      ti <- t_obs[i]
      cand <- which(cmat[i, ])

      if (length(cand) < k) {
        warning(sprintf("Observation %d: candidate set size %d < k=%d", i, length(cand), k))
        return(-Inf)
      }

      # Enumerate all k-subsets of the candidate set
      if (length(cand) == k) {
        # No masking: F = C exactly
        subsets <- list(cand)
      } else {
        subsets <- utils::combn(cand, k, simplify = FALSE)
      }

      # Compute f_j, F_j, S_j at t_i for all components
      f_vals <- vapply(seq_len(m), function(j)
        stats::dweibull(ti, shape = shapes[j], scale = scales[j]), numeric(1))
      F_vals <- vapply(seq_len(m), function(j)
        stats::pweibull(ti, shape = shapes[j], scale = scales[j]), numeric(1))
      S_vals <- 1 - F_vals

      # Sum over all candidate failed sets
      total <- 0
      for (F_set in subsets) {
        not_F <- setdiff(seq_len(m), F_set)

        # For each j in F: f_j(t) * prod_{l in F\j} F_l(t) * prod_{l not in F} S_l(t)
        for (j in F_set) {
          if (f_vals[j] <= 0) next
          others_in_F <- setdiff(F_set, j)
          log_term <- log(f_vals[j])
          if (length(others_in_F) > 0) {
            if (any(F_vals[others_in_F] <= 0)) next
            log_term <- log_term + sum(log(F_vals[others_in_F]))
          }
          if (length(not_F) > 0) {
            if (any(S_vals[not_F] <= 0)) next
            log_term <- log_term + sum(log(S_vals[not_F]))
          }
          total <- total + exp(log_term)
        }
      }

      if (total <= 0) return(-Inf)
      ll <- ll + log(total)
    }

    ll
  }
}


#' Generate masked k-out-of-n data
#'
#' Returns a closure that generates system lifetime data with candidate
#' failed sets. The true failed set (the k components that failed by
#' T_sys) is always included in the candidate set. Additional non-failed
#' components are included independently with probability \code{p_mask}.
#'
#' @param model A \code{kofn} model object.
#' @param candset Prefix for candidate set columns (default \code{"c"}).
#' @param ... Additional arguments (currently unused).
#' @return A function \code{function(theta, n, p_mask = 0, observe = NULL)}
#'   returning a data frame with columns for system lifetime, observation
#'   type, and Boolean candidate set indicators.
#'
#' @examples
#' model <- kofn(k = 2, m = 4, component = dfr_exponential())
#' gen <- rdata_masked(model)
#' set.seed(42)
#' df <- gen(theta = c(1, 0.8, 0.6, 0.4), n = 10, p_mask = 0.3)
#' head(df)
#'
#' @export
rdata_masked <- function(model, candset = "c", ...) {
  m <- model$m
  k <- model$k
  lt <- model$lifetime
  om <- model$omega
  component <- model$component

  if (is.na(k)) stop("Masked data generation requires a k-out-of-n system")

  function(theta, n, p_mask = 0, observe = NULL) {
    stopifnot(length(theta) == n_par_kofn(m, component))

    pp <- parse_params(theta, m, component)
    shapes <- pp$shapes
    scales <- pp$scales
    stopifnot(all(shapes > 0), all(scales > 0))

    # Generate component lifetimes
    comp_times <- matrix(0, nrow = n, ncol = m)
    for (j in seq_len(m)) {
      comp_times[, j] <- stats::rweibull(n, shape = shapes[j],
                                         scale = scales[j])
    }

    # Compute system lifetime and identify failed components
    sys_times <- numeric(n)
    failed_sets <- vector("list", n)
    for (i in seq_len(n)) {
      cinfo <- kofn_censoring(k, comp_times[i, ])
      sys_times[i] <- cinfo$T_sys
      # Failed set: components with status "exact" or "left"
      failed_sets[[i]] <- which(cinfo$status != "right")
    }

    # Apply observation functor
    if (is.null(observe)) observe <- observe_exact()
    obs_t      <- numeric(n)
    omega_vals <- character(n)
    for (i in seq_len(n)) {
      obs <- observe(sys_times[i])
      obs_t[i]      <- obs$t
      omega_vals[i] <- obs$omega
    }

    # Generate candidate failed sets
    cand_matrix <- matrix(FALSE, nrow = n, ncol = m)
    for (i in seq_len(n)) {
      if (omega_vals[i] == "right") next  # no masking for censored obs
      F_true <- failed_sets[[i]]
      cand_matrix[i, F_true] <- TRUE  # true failed set always included
      if (p_mask > 0) {
        survivors <- setdiff(seq_len(m), F_true)
        if (length(survivors) > 0) {
          cand_matrix[i, survivors] <- stats::runif(length(survivors)) < p_mask
        }
      }
    }

    # Build data frame
    df <- data.frame(obs_t, omega_vals, stringsAsFactors = FALSE)
    names(df) <- c(lt, om)
    for (j in seq_len(m)) {
      df[[paste0(candset, j)]] <- cand_matrix[, j]
    }

    attr(df, "comp_times") <- comp_times
    attr(df, "failed_sets") <- failed_sets
    attr(df, "par") <- theta
    df
  }
}
