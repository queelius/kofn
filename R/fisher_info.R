# ===========================================================================
# Fisher Information Comparison Across Observation Schemes
# ===========================================================================
#
# Compares the observed Fisher information under:
#   - Scheme 0: System-level only (T = max T_j, no component data)
#   - Scheme 1: Periodic inspection (interval-censored component times)
#   - Scheme 2: Complete data (all component times observed)
#
# The comparison uses the determinant of the observed Fisher information
# matrix as a scalar summary of total information. Efficiency ratios
# (det ratios) quantify the information gain from richer observation.
# ===========================================================================


#' Determinant of a positive-definite Hessian, or NA
#'
#' Computes the numerical Hessian, checks finite-ness and positive
#' definiteness, and returns \code{det(H)} or \code{NA}.
#'
#' @param neg_ll Negative log-likelihood function.
#' @param par Parameter vector.
#' @return Scalar determinant, or \code{NA_real_}.
#' @keywords internal
safe_hessian_det <- function(neg_ll, par) {
  H <- tryCatch(
    numDeriv::hessian(neg_ll, par, method = "Richardson"),
    error = function(e) NULL
  )
  if (is.null(H) || !all(is.finite(H))) return(NA_real_)
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  if (all(ev > 0)) det(H) else NA_real_
}


#' Compare Fisher information across observation schemes
#'
#' For a given parameter configuration, computes the observed Fisher
#' information under Scheme 0 (system only), Scheme 1 (periodic inspection),
#' and Scheme 2 (complete component data) via simulation.
#'
#' The determinant of the observed information matrix is used as a scalar
#' summary. Efficiency ratios indicate relative information content:
#' a ratio less than 1 means the denominator scheme carries more information.
#'
#' @param shapes Numeric vector. Weibull shape parameters (for Weibull family)
#'   or \code{NULL} (for exponential).
#' @param scales Numeric vector. Weibull scale parameters (for Weibull family)
#'   or \code{NULL} (for exponential).
#' @param rates Numeric vector. Exponential rate parameters (for exponential
#'   family, alternative to \code{shapes}/\code{scales}).
#' @param n Integer. Sample size per replicate.
#' @param delta Numeric scalar. Inspection interval width for Scheme 1.
#' @param n_rep Integer. Number of Monte Carlo replicates.
#' @param family Character. Component distribution family:
#'   \code{"exponential"} or \code{"weibull"}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{scheme0_det}}{Numeric vector of Scheme 0 information
#'       determinants (length \code{n_rep}).}
#'     \item{\code{scheme1_det}}{Numeric vector of Scheme 1 information
#'       determinants.}
#'     \item{\code{scheme2_det}}{Numeric vector of Scheme 2 (complete data)
#'       information determinants.}
#'     \item{\code{median_det}}{Named numeric vector of median determinants
#'       across schemes.}
#'     \item{\code{efficiency_01}}{Ratio of median Scheme 0 to Scheme 1
#'       determinant (< 1 means Scheme 1 is more informative).}
#'     \item{\code{efficiency_02}}{Ratio of median Scheme 0 to Scheme 2
#'       determinant.}
#'     \item{\code{efficiency_12}}{Ratio of median Scheme 1 to Scheme 2
#'       determinant.}
#'     \item{\code{n_valid}}{Named integer vector of valid (non-NA) replicate
#'       counts per scheme.}
#'   }
#'
#' @details
#' For Scheme 2 (complete data), the Fisher information is computed
#' analytically:
#' \itemize{
#'   \item Exponential: \eqn{I_{jj} = n / \lambda_j^2} (diagonal).
#'   \item Weibull: Block-diagonal with per-component 2x2 Fisher information
#'     matrices using the standard Weibull FIM formulas.
#' }
#'
#' For Schemes 0 and 1, the observed information is computed numerically
#' via \code{numDeriv::hessian} of the negative log-likelihood evaluated
#' at the true parameter values.
#'
#' @importFrom numDeriv hessian
#' @export
compare_fisher_info <- function(shapes = NULL, scales = NULL, rates = NULL,
                                n = 200L, delta = 1.0, n_rep = 50L,
                                family = "exponential") {
  family <- match.arg(family, c("exponential", "weibull"))

  if (family == "exponential") {
    if (is.null(rates)) stop("rates required for exponential family")
    m <- length(rates)
    par_true <- rates
  } else {
    if (is.null(shapes) || is.null(scales)) {
      stop("shapes and scales required for weibull family")
    }
    m <- length(shapes)
    par_true <- as.numeric(rbind(shapes, scales))
  }
  pp <- parse_params(par_true, m, family)
  shapes <- pp$shapes
  scales <- pp$scales

  n_par <- length(par_true)

  # Storage
  det_s0 <- numeric(n_rep)
  det_s1 <- numeric(n_rep)
  det_s2 <- numeric(n_rep)

  for (rep in seq_len(n_rep)) {
    # Generate component times
    comp_times <- matrix(0, nrow = n, ncol = m)
    for (j in seq_len(m)) {
      comp_times[, j] <- stats::rweibull(n, shape = shapes[j],
                                         scale = scales[j])
    }
    sys_times <- apply(comp_times, 1, max)

    # --- Scheme 2: Complete data Fisher info ---
    if (family == "exponential") {
      I_s2 <- diag(n / rates^2)
    } else {
      I_s2 <- matrix(0, nrow = n_par, ncol = n_par)
      gamma_euler <- -digamma(1)
      for (j in seq_len(m)) {
        alpha_j <- shapes[j]
        beta_j <- scales[j]
        idx <- (2 * (j - 1) + 1):(2 * j)
        I_11 <- (1 + (1 - gamma_euler)^2 + pi^2 / 6) / alpha_j^2
        I_12 <- -(1 + digamma(1)) / (alpha_j * beta_j)
        I_22 <- alpha_j^2 / beta_j^2
        I_s2[idx, idx] <- n * matrix(c(I_11, I_12, I_12, I_22), 2, 2)
      }
    }
    det_s2[rep] <- det(I_s2)

    # --- Scheme 0: System-level only ---
    neg_ll_s0 <- function(p) {
      if (any(!is.finite(p)) || any(p <= 0)) return(.Machine$double.xmax / 2)
      pp_s0 <- parse_params(p, m, family)
      ll <- 0
      for (i in seq_along(sys_times)) {
        f_sys <- weibull_f_sys(sys_times[i], pp_s0$shapes, pp_s0$scales)
        if (f_sys <= 0) return(.Machine$double.xmax / 2)
        ll <- ll + log(f_sys)
      }
      -ll
    }

    det_s0[rep] <- safe_hessian_det(neg_ll_s0, par_true)

    # --- Scheme 1: Periodic inspection ---
    df_s1 <- data.frame(t = sys_times)
    for (j in seq_len(m)) {
      lower_j <- floor(comp_times[, j] / delta) * delta
      df_s1[[paste0("comp_lower_", j)]] <- lower_j
      df_s1[[paste0("comp_upper_", j)]] <- lower_j + delta
    }

    lower_cols <- paste0("comp_lower_", seq_len(m))
    upper_cols <- paste0("comp_upper_", seq_len(m))

    neg_ll_s1 <- function(p) {
      if (any(!is.finite(p)) || any(p <= 0)) return(.Machine$double.xmax / 2)
      pp_s1 <- parse_params(p, m, family)
      sh <- pp_s1$shapes
      sc <- pp_s1$scales
      ll <- 0
      for (i in seq_len(n)) {
        f_sys <- weibull_f_sys(df_s1$t[i], sh, sc)
        if (f_sys <= 0) return(.Machine$double.xmax / 2)
        ll <- ll + log(f_sys)

        for (j in seq_len(m)) {
          a_lower <- df_s1[[lower_cols[j]]][i]
          a_upper <- df_s1[[upper_cols[j]]][i]
          F_upper <- stats::pweibull(a_upper, shape = sh[j], scale = sc[j])
          F_lower <- stats::pweibull(a_lower, shape = sh[j], scale = sc[j])
          prob <- F_upper - F_lower
          if (prob <= 0) return(.Machine$double.xmax / 2)
          ll <- ll + log(prob)
        }
      }
      -ll
    }

    det_s1[rep] <- safe_hessian_det(neg_ll_s1, par_true)
  }

  # Compute efficiency ratios (median over replicates, ignoring NAs)
  med_s0 <- stats::median(det_s0, na.rm = TRUE)
  med_s1 <- stats::median(det_s1, na.rm = TRUE)
  med_s2 <- stats::median(det_s2, na.rm = TRUE)

  list(
    scheme0_det = det_s0,
    scheme1_det = det_s1,
    scheme2_det = det_s2,
    median_det  = c(scheme0 = med_s0, scheme1 = med_s1, scheme2 = med_s2),
    efficiency_01 = med_s0 / med_s1,
    efficiency_02 = med_s0 / med_s2,
    efficiency_12 = med_s1 / med_s2,
    n_valid = c(scheme0 = sum(!is.na(det_s0)),
                scheme1 = sum(!is.na(det_s1)),
                scheme2 = n_rep)
  )
}
