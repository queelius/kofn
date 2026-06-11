# ===========================================================================
# Exponential k-out-of-n Model: S3 Methods
# ===========================================================================
#
# Provides S3 methods for the `exp_kofn` class created by
# kofn(component = dfr_exponential()).
#
# Each generic returns a closure, following the likelihood.model /
# maskedcauses convention.
#
# Strategy:
#   - Parallel (k = m) loglik uses the IE expansion (closed form for all
#     four observation types: exact, right, left, interval). This is
#     kept locally because IE-based integral expressions are kofn's
#     distinctive contribution.
#   - General-k loglik delegates to dist.structure: density at exact
#     observations, survival at right-censored, CDF at left, and
#     CDF differences at interval-censored.
# ===========================================================================


# Internal: log-likelihood contribution helpers using dist.structure
# closures bound to a particular DGP. Returned closure is reused across
# observations to avoid reconstructing the dist object per row.
ll_via_dgp <- function(dgp, t_obs, omega_vals, t_upper, lt_up_avail) {
  surv_fn <- algebraic.dist::surv(dgp)
  cdf_fn  <- algebraic.dist::cdf(dgp)
  dens_fn <- density(dgp)

  ll <- 0
  for (i in seq_along(t_obs)) {
    omi <- omega_vals[i]
    ti  <- t_obs[i]
    if (omi == "exact") {
      val <- dens_fn(ti)
    } else if (omi == "right") {
      val <- surv_fn(ti)
    } else if (omi == "left") {
      val <- cdf_fn(ti)
    } else if (omi == "interval") {
      if (!lt_up_avail) {
        stop("interval-censored observations require an upper-bound column")
      }
      tu <- t_upper[i]
      val <- cdf_fn(tu) - cdf_fn(ti)
    } else {
      stop(sprintf("Unknown omega value: %s", omi))
    }
    if (!is.finite(val) || val <= 0) return(-Inf)
    ll <- ll + log(val)
  }
  ll
}


# ---------------------------------------------------------------------------
# loglik.exp_kofn
# ---------------------------------------------------------------------------

#' Log-likelihood for exponential k-out-of-n model
#'
#' Returns a closure `function(df, par)` that computes the log-likelihood of
#' exponential component rates given system-level data.
#'
#' For parallel systems (`k = m`), the closed-form IE expansion is used,
#' supporting all four observation types: exact, right, left, and interval
#' censored. For general k-out-of-n systems, the per-observation
#' contributions are computed via `dist.structure::exp_kofn(k, par)` and
#' the algebraic.dist generics `density`, `surv`, and `cdf`.
#'
#' @param model An `exp_kofn` object created by [kofn()].
#' @param ... Additional arguments (ignored).
#' @return A function `function(df, par)` where `df` is a data frame with
#'   columns for lifetime, observation type, and (for interval-censored)
#'   upper bounds, and `par` is a numeric vector of `m` component rates.
#'
#' @examples
#' model <- kofn(k = 3, m = 3, component = dfr_exponential())
#' ll <- loglik(model)
#' set.seed(1)
#' df <- rdata(model)(c(1, 2, 3), n = 30)
#' ll(df, c(1, 2, 3))
#'
#' @method loglik exp_kofn
#' @export
loglik.exp_kofn <- function(model, ...) {
  m     <- model$m
  k     <- model$k
  lt    <- model$lifetime
  om    <- model$omega
  lt_up <- model$lifetime_upper

  if (isTRUE(k == m)) {
    # Parallel fast path: IE-based closed-form likelihood
    function(df, par) {
      if (any(par <= 0)) return(-Inf)
      if (length(par) != m) {
        stop(sprintf("Expected %d parameters but got %d", m, length(par)))
      }

      t_obs <- df[[lt]]
      omega_vals <- as.character(df[[om]])
      comps <- seq_len(m)
      ll <- 0

      # Exact: log f_sys(t)
      for (i in which(omega_vals == "exact")) {
        val <- sum(vapply(comps,
          function(j) w_j_exact(t_obs[i], par, j), numeric(1)))
        if (val <= 0) return(-Inf)
        ll <- ll + log(val)
      }

      # Right-censored: log S_sys(t)
      for (i in which(omega_vals == "right")) {
        s <- S_sys_exp(t_obs[i], par)
        if (s <= 0) return(-Inf)
        ll <- ll + log(s)
      }

      # Left-censored: log F_sys(t)
      for (i in which(omega_vals == "left")) {
        val <- sum(vapply(comps,
          function(j) w_j_integral(0, t_obs[i], par, j), numeric(1)))
        if (val <= 0) return(-Inf)
        ll <- ll + log(val)
      }

      # Interval-censored: log(F_sys(b) - F_sys(a))
      int_idx <- which(omega_vals == "interval")
      if (length(int_idx) > 0 && !(lt_up %in% names(df))) {
        stop(sprintf(
          "interval-censored observations require column '%s'", lt_up))
      }
      for (i in int_idx) {
        a <- t_obs[i]
        b <- df[[lt_up]][i]
        val <- sum(vapply(comps,
          function(j) w_j_integral(a, b, par, j), numeric(1)))
        if (val <= 0) return(-Inf)
        ll <- ll + log(val)
      }

      ll
    }
  } else {
    # General k path: delegate per-observation contributions to dist.structure.
    function(df, par) {
      if (any(par <= 0)) return(-Inf)
      if (length(par) != m) {
        stop(sprintf("Expected %d parameters but got %d", m, length(par)))
      }
      t_obs <- df[[lt]]
      omega_vals <- as.character(df[[om]])
      lt_up_avail <- (lt_up %in% names(df))
      t_upper <- if (lt_up_avail) df[[lt_up]] else NULL
      dgp <- kofn_dgp(k, m, model$component, par)
      ll_via_dgp(dgp, t_obs, omega_vals, t_upper, lt_up_avail)
    }
  }
}


# ---------------------------------------------------------------------------
# score.exp_kofn
# ---------------------------------------------------------------------------

#' Score function for exponential k-out-of-n model
#'
#' Returns a closure `function(df, par)` that computes the gradient of
#' the log-likelihood with respect to the rate parameters.
#'
#' Uses numerical differentiation via [numDeriv::grad()] applied to the
#' log-likelihood closure.
#'
#' @param model An `exp_kofn` object created by [kofn()].
#' @param ... Additional arguments (ignored).
#' @return A function `function(df, par)` returning a numeric vector of
#'   length `m` (the score vector).
#'
#' @examples
#' model <- kofn(k = 2, m = 2, component = dfr_exponential())
#' sc <- score(model)
#' set.seed(1)
#' df <- rdata(model)(c(1, 2), n = 30)
#' sc(df, c(1, 2))  # gradient at true parameters
#'
#' @method score exp_kofn
#' @export
score.exp_kofn <- function(model, ...) {
  ll_fn <- loglik(model)

  function(df, par) {
    numDeriv::grad(
      func = function(p) ll_fn(df, p),
      x = par,
      method = "Richardson"
    )
  }
}


# ---------------------------------------------------------------------------
# hess_loglik.exp_kofn
# ---------------------------------------------------------------------------

#' Hessian of the log-likelihood for exponential k-out-of-n model
#'
#' Returns a closure `function(df, par)` that computes the Hessian matrix
#' of the log-likelihood with respect to the rate parameters.
#'
#' Uses numerical differentiation via [numDeriv::hessian()] applied to the
#' log-likelihood closure.
#'
#' @param model An `exp_kofn` object created by [kofn()].
#' @param ... Additional arguments (ignored).
#' @return A function `function(df, par)` returning an `m x m` Hessian matrix.
#'
#' @examples
#' model <- kofn(k = 2, m = 2, component = dfr_exponential())
#' H <- hess_loglik(model)
#' set.seed(1)
#' df <- rdata(model)(c(1, 2), n = 30)
#' H(df, c(1, 2))  # 2x2 Hessian matrix
#'
#' @method hess_loglik exp_kofn
#' @export
hess_loglik.exp_kofn <- function(model, ...) {
  ll_fn <- loglik(model)

  function(df, par) {
    numDeriv::hessian(
      func = function(p) ll_fn(df, p),
      x = par,
      method = "Richardson"
    )
  }
}


# ---------------------------------------------------------------------------
# fit.exp_kofn
# ---------------------------------------------------------------------------

#' MLE fitting for exponential k-out-of-n model
#'
#' Returns a closure `function(df, par0, n_starts)` that computes the
#' maximum likelihood estimate of the component rate parameters.
#'
#' Uses multi-start optimization with L-BFGS-B as primary solver and
#' Nelder-Mead on the log-scale as fallback. Standard errors are computed
#' from the observed Fisher information (negative Hessian) at the MLE.
#'
#' @param object An `exp_kofn` object created by [kofn()].
#' @param ... Additional arguments (ignored).
#' @return A function `function(df, par0 = NULL, n_starts = 5L)` returning
#'   a fisher_mle object (from likelihood.model).
#'
#' @examples
#' \donttest{
#' model <- kofn(k = 2, m = 2, component = dfr_exponential())
#' set.seed(42)
#' df <- rdata(model)(c(1, 2), n = 50)
#' result <- fit(model)(df)
#' coef(result)
#' }
#'
#' @method fit exp_kofn
#' @export
fit.exp_kofn <- function(object, ...) {
  ll_fn <- loglik(object)
  m     <- object$m
  lt    <- object$lifetime
  om    <- object$omega

  function(df, par0 = NULL, n_starts = 5L) {
    if (is.null(par0)) {
      par0 <- default_init_exp(df, m, lt, om, object$lifetime_upper)
    }
    if (length(par0) != m) {
      stop(sprintf("Initial parameter vector has length %d but model has %d components",
                   length(par0), m))
    }

    neg_ll <- function(par) {
      val <- -ll_fn(df, par)
      if (!is.finite(val)) return(.Machine$double.xmax / 2)
      val
    }

    solve_mle(neg_ll, par0, n_par = m, n_starts = n_starts,
              nobs = nrow(df))
  }
}


# ---------------------------------------------------------------------------
# rdata.exp_kofn
# ---------------------------------------------------------------------------

#' Random data generation for exponential k-out-of-n model
#'
#' Returns a closure that generates random system-level data from the
#' exponential k-out-of-n data-generating process.
#'
#' Workflow:
#' 1. Generate i.i.d. exponential component lifetimes.
#' 2. Compute system lifetime as the (m - k + 1)-th order statistic.
#' 3. Apply the observation functor (exact observation by default).
#'
#' @param model An `exp_kofn` object created by [kofn()].
#' @param ... Additional arguments (ignored).
#' @return A function `function(theta, n, observe = NULL)`
#'   returning a data frame with columns \code{t} (system lifetime) and
#'   \code{omega} (observation type). Latent component lifetimes, true
#'   critical component, and the true parameters are stored as attributes.
#'
#' @examples
#' model <- kofn(k = 3, m = 3, component = dfr_exponential())
#' gen <- rdata(model)
#' set.seed(1)
#' df <- gen(theta = c(1, 2, 3), n = 20)
#' head(df)
#'
#' # With right-censoring
#' df2 <- gen(c(1, 2, 3), n = 20, observe = observe_right_censor(tau = 2))
#' table(df2$omega)
#'
#' @method rdata exp_kofn
#' @export
rdata.exp_kofn <- function(model, ...) {
  m   <- model$m
  k   <- model$k
  lt  <- model$lifetime
  om  <- model$omega
  lt_up <- model$lifetime_upper

  function(theta, n, observe = NULL) {
    if (length(theta) != m) {
      stop(sprintf("theta has length %d but model has %d components",
                   length(theta), m))
    }
    if (any(theta <= 0)) stop("All rates must be positive")

    # 1. Generate component lifetimes
    comp_lifetimes <- matrix(nrow = n, ncol = m)
    for (j in seq_len(m)) {
      comp_lifetimes[, j] <- stats::rexp(n, rate = theta[j])
    }

    # 2. Compute system lifetime + critical component (k-of-n is the
    # (m - k + 1)-th order statistic; the critical component is the one
    # whose failure equals T_sys, deterministic for absolutely continuous
    # components).
    sys_lifetime_vec <- numeric(n)
    critical_comp    <- integer(n)
    for (i in seq_len(n)) {
      cinfo <- kofn_censoring(k, comp_lifetimes[i, ])
      sys_lifetime_vec[i] <- cinfo$T_sys
      critical_comp[i]    <- cinfo$critical
    }

    # 3. Apply observation functor (default: exact, no censoring)
    if (is.null(observe)) {
      observe <- observe_exact()
    }
    obs_t       <- numeric(n)
    omega_vals  <- character(n)
    t_upper_vals <- rep(NA_real_, n)
    for (i in seq_len(n)) {
      obs <- observe(sys_lifetime_vec[i])
      obs_t[i]        <- obs$t
      omega_vals[i]   <- obs$omega
      t_upper_vals[i] <- obs$t_upper
    }

    # Build data frame
    df <- data.frame(obs_t, omega_vals, stringsAsFactors = FALSE)
    names(df) <- c(lt, om)
    if (any(omega_vals == "interval")) {
      df[[lt_up]] <- t_upper_vals
    }

    attr(df, "comp_lifetimes") <- comp_lifetimes
    attr(df, "critical_comp")  <- critical_comp
    attr(df, "theta")          <- theta
    df
  }
}


# ---------------------------------------------------------------------------
# assumptions.exp_kofn
# ---------------------------------------------------------------------------

#' Assumptions for exponential k-out-of-n model
#'
#' Returns a character vector listing the assumptions made by the
#' exponential k-out-of-n likelihood model.
#'
#' @param model An `exp_kofn` object created by [kofn()].
#' @param ... Additional arguments (ignored).
#' @return Character vector of assumptions.
#'
#' @examples
#' model <- kofn(k = 3, m = 3, component = dfr_exponential())
#' assumptions(model)
#'
#' @method assumptions exp_kofn
#' @export
assumptions.exp_kofn <- function(model, ...) {
  c(
    "independent component lifetimes",
    "exponential component lifetime distributions",
    sprintf("%d-out-of-%d system structure", model$k, model$m),
    "i.i.d. system observations"
  )
}


# ===========================================================================
# Internal helpers
# ===========================================================================

#' Compute default initial values for exponential rates
#'
#' Uses the method-of-moments estimator based on the harmonic number
#' relationship for i.i.d. parallel exponentials: E[T] = H_m / lambda.
#' Returns a spread of initial values to break permutation symmetry.
#'
#' @param df Data frame with lifetime observations.
#' @param m Number of components.
#' @param lifetime Column name for lifetime.
#' @param omega Column name for observation type.
#' @param lifetime_upper Column name for interval upper bound.
#' @return Numeric vector of length `m` with initial rate estimates.
#' @keywords internal
default_init_exp <- function(df, m, lifetime, omega, lifetime_upper) {
  t_obs      <- df[[lifetime]]
  omega_vals <- as.character(df[[omega]])

  # Use interval midpoints where applicable
  t_mid <- t_obs
  if (!is.null(lifetime_upper) && lifetime_upper %in% names(df)) {
    int_idx <- omega_vals == "interval"
    if (any(int_idx)) {
      t_mid[int_idx] <- (t_obs[int_idx] + df[[lifetime_upper]][int_idx]) / 2
    }
  }

  # For parallel exponential with equal rates lambda:
  #   E[T_sys] = H_m / lambda, where H_m = sum(1/k, k=1..m)
  # So lambda ~ H_m / mean(T_sys).
  usable <- omega_vals != "right"
  mean_t <- if (sum(usable) > 0) mean(t_mid[usable]) else mean(t_mid)
  mean_t <- max(mean_t, 0.01)  # guard against degenerate data
  H_m  <- sum(1 / seq_len(m))
  lam0 <- H_m / mean_t

  # Spread around the moment estimate to break permutation symmetry.
  lam0 * seq(0.5, 1.5, length.out = m)
}
