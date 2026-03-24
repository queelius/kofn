# ===========================================================================
# Exponential k-out-of-n Model — S3 Methods
# ===========================================================================
#
# Provides S3 methods for the `exp_kofn` class created by
# kofn(family = "exponential").
#
# Each generic returns a closure, following the likelihood.model / maskedcauses
# convention.
#
# For parallel systems (k=m), the inclusion-exclusion (IE) expansion yields
# closed-form expressions for all integrals. For general k-out-of-n systems,
# the loglik falls back to the general system density engine
# (loglik_system from R/system_density.R).
# ===========================================================================


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
#' censored. For general k-out-of-n systems, the computation delegates to
#' [loglik_system()] from the system density engine.
#'
#' @param model An `exp_kofn` object created by [kofn()].
#' @param ... Additional arguments (ignored).
#' @return A function `function(df, par)` where `df` is a data frame with
#'   columns for lifetime, observation type, and candidate set indicators,
#'   and `par` is a numeric vector of `m` component rates.
#'
#' @examples
#' model <- kofn(k = 3, m = 3, family = "exponential")
#' ll <- loglik(model)
#' set.seed(1)
#' df <- rdata(model)(c(1, 2, 3), n = 30)
#' ll(df, c(1, 2, 3))
#'
#' @method loglik exp_kofn
#' @export
loglik.exp_kofn <- function(model, ...) {
  sys   <- model$system
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

      # --- Exact: log f_sys(t) ---
      for (i in which(omega_vals == "exact")) {
        val <- sum(vapply(comps,
          function(j) w_j_exact(t_obs[i], par, j), numeric(1)))
        if (val <= 0) return(-Inf)
        ll <- ll + log(val)
      }

      # --- Right-censored: log S_sys(t) ---
      for (i in which(omega_vals == "right")) {
        s <- S_sys_exp(t_obs[i], par)
        if (s <= 0) return(-Inf)
        ll <- ll + log(s)
      }

      # --- Left-censored: log F_sys(t) ---
      for (i in which(omega_vals == "left")) {
        val <- sum(vapply(comps,
          function(j) w_j_integral(0, t_obs[i], par, j), numeric(1)))
        if (val <= 0) return(-Inf)
        ll <- ll + log(val)
      }

      # --- Interval-censored: log(F_sys(b) - F_sys(a)) ---
      for (i in which(omega_vals == "interval")) {
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
    # General k-out-of-n: delegate to system density engine
    function(df, par) {
      if (any(par <= 0)) return(-Inf)
      t_obs <- df[[lt]]
      loglik_system(t_obs, sys, par, family = "exponential")
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
#' model <- kofn(k = 2, m = 2, family = "exponential")
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
#' model <- kofn(k = 2, m = 2, family = "exponential")
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
#'   a list with components:
#'   \describe{
#'     \item{par}{Numeric vector of MLE rate estimates.}
#'     \item{se}{Standard errors (from observed Fisher information).}
#'     \item{loglik}{Maximized log-likelihood value.}
#'     \item{convergence}{Integer convergence code (0 = success).}
#'     \item{fisher_info}{Observed Fisher information matrix (or NULL).}
#'   }
#'
#' @examples
#' \donttest{
#' model <- kofn(k = 2, m = 2, family = "exponential")
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

    multistart_mle(neg_ll, par0, n_par = m, n_starts = n_starts,
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
#' 2. Compute system lifetime via the coherent system structure function.
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
#' model <- kofn(k = 3, m = 3, family = "exponential")
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
  sys <- model$system
  m   <- model$m
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

    # 2. Compute system lifetime + critical component
    sys_lifetime_vec <- numeric(n)
    critical_comp    <- integer(n)
    for (i in seq_len(n)) {
      cinfo <- system_censoring(sys, comp_lifetimes[i, ])
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
#' model <- kofn(k = 3, m = 3, family = "exponential")
#' assumptions(model)
#'
#' @method assumptions exp_kofn
#' @export
assumptions.exp_kofn <- function(model, ...) {
  sys_desc <- if (is.na(model$k)) {
    "general coherent system structure"
  } else {
    sprintf("%d-out-of-%d system structure", model$k, model$m)
  }
  c(
    "independent component lifetimes",
    "exponential component lifetime distributions",
    sys_desc,
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
  #   E[T_sys] = H_m / lambda,  where H_m = sum(1/k, k=1..m)
  # So lambda ~ H_m / mean(T_sys)
  usable <- omega_vals != "right"
  mean_t <- if (sum(usable) > 0) mean(t_mid[usable]) else mean(t_mid)
  mean_t <- max(mean_t, 0.01)  # guard against degenerate data
  H_m  <- sum(1 / seq_len(m))
  lam0 <- H_m / mean_t

  # Spread around the moment estimate to break permutation symmetry
  lam0 * seq(0.5, 1.5, length.out = m)
}
