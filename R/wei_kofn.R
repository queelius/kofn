# ===========================================================================
# Weibull k-out-of-n System Model (S3 methods)
# ===========================================================================
#
# Provides S3 methods for the `wei_kofn` class, created by
# kofn(family = "weibull"). Implements the likelihood_model concept from
# the likelihood.model package.
#
# Parameter convention: par is a vector of length 2*m with interleaved
# (shape_1, scale_1, shape_2, scale_2, ..., shape_m, scale_m).
#
# For parallel systems (k=1), an EM algorithm is available as an alternative
# to direct MLE. The EM treats the identity of the last-failing component
# as latent data, using truncated Weibull moments for the E-step and
# profile optimization over shape with closed-form scale for the M-step.
# ===========================================================================


# ---------------------------------------------------------------------------
# Internal: Weibull parallel system density components
# ---------------------------------------------------------------------------

#' Component weight for Weibull parallel system (internal)
#'
#' Computes w_j(t) = f_j(t) * prod_{i != j} F_i(t), the contribution of
#' component j to the parallel system density at time t.
#'
#' @param t Numeric scalar. Time point (positive).
#' @param shapes Numeric vector. Weibull shape parameters (length m).
#' @param scales Numeric vector. Weibull scale parameters (length m).
#' @param j Integer scalar. Component index (1 to m).
#' @return Numeric scalar w_j(t).
#' @keywords internal
weibull_w_j <- function(t, shapes, scales, j) {
  if (t <= 0) return(0)
  m <- length(shapes)
  f_j <- stats::dweibull(t, shape = shapes[j], scale = scales[j])
  if (f_j <= 0) return(0)
  log_F_others <- 0
  for (i in seq_len(m)) {
    if (i == j) next
    F_i <- stats::pweibull(t, shape = shapes[i], scale = scales[i])
    if (F_i <= 0) return(0)
    log_F_others <- log_F_others + log(F_i)
  }
  f_j * exp(log_F_others)
}


#' System density for Weibull parallel system (internal)
#'
#' Computes f_sys(t) = sum_j f_j(t) * prod_{i != j} F_i(t).
#'
#' @param t Numeric scalar. Time point.
#' @param shapes Numeric vector. Weibull shape parameters.
#' @param scales Numeric vector. Weibull scale parameters.
#' @return Numeric scalar f_sys(t).
#' @keywords internal
weibull_f_sys <- function(t, shapes, scales) {
  m <- length(shapes)
  sum(vapply(seq_len(m), function(j) weibull_w_j(t, shapes, scales, j),
             numeric(1)))
}


# ---------------------------------------------------------------------------
# Internal: extract shapes and scales from interleaved par vector
# ---------------------------------------------------------------------------

#' @keywords internal
wei_par_to_shapes <- function(par, m) {
  par[seq(1L, 2L * m, by = 2L)]
}

#' @keywords internal
wei_par_to_scales <- function(par, m) {
  par[seq(2L, 2L * m, by = 2L)]
}


# ---------------------------------------------------------------------------
# loglik.wei_kofn
# ---------------------------------------------------------------------------

#' Log-likelihood for Weibull k-out-of-n system
#'
#' Returns a closure \code{function(df, par)} that computes the log-likelihood
#' for a Weibull k-out-of-n system given data and parameters.
#'
#' For parallel systems (k = 1), uses the direct Weibull system density
#' \eqn{f_{sys}(t) = \sum_j f_j(t) \prod_{i \ne j} F_i(t)}.
#'
#' For general k, delegates to [loglik_system()] with \code{family = "weibull"}.
#'
#' @param model A \code{wei_kofn} object.
#' @param ... Additional arguments (currently unused).
#' @return A function \code{function(df, par)} returning a scalar log-likelihood.
#'
#' @export
loglik.wei_kofn <- function(model, ...) {
  m <- model$m
  k <- model$k
  system <- model$system
  lt <- model$lifetime

  if (k == 1L) {
    # Parallel system: direct density computation
    function(df, par) {
      if (any(!is.finite(par)) || any(par <= 0)) return(-Inf)
      shapes <- wei_par_to_shapes(par, m)
      scales <- wei_par_to_scales(par, m)
      t_obs <- df[[lt]]
      vals <- vapply(t_obs, function(t) {
        f <- weibull_f_sys(t, shapes, scales)
        if (f <= 0) return(-1e20)
        log(f)
      }, numeric(1))
      sum(vals)
    }
  } else {
    # General k-out-of-n: delegate to system density engine
    function(df, par) {
      if (any(!is.finite(par)) || any(par <= 0)) return(-Inf)
      t_obs <- df[[lt]]
      loglik_system(t_obs, system, par, family = "weibull")
    }
  }
}


# ---------------------------------------------------------------------------
# score.wei_kofn
# ---------------------------------------------------------------------------

#' Score function for Weibull k-out-of-n system
#'
#' Returns a closure that computes the score (gradient of log-likelihood)
#' via numerical differentiation using \code{numDeriv::grad}.
#'
#' @param model A \code{wei_kofn} object.
#' @param ... Additional arguments (currently unused).
#' @return A function \code{function(df, par)} returning a numeric vector
#'   (the gradient of the log-likelihood).
#'
#' @importFrom numDeriv grad
#' @export
score.wei_kofn <- function(model, ...) {
  ll_fn <- loglik(model, ...)
  function(df, par) {
    numDeriv::grad(func = function(p) ll_fn(df, p), x = par,
                   method = "Richardson")
  }
}


# ---------------------------------------------------------------------------
# hess_loglik.wei_kofn
# ---------------------------------------------------------------------------

#' Hessian of log-likelihood for Weibull k-out-of-n system
#'
#' Returns a closure that computes the Hessian matrix of the log-likelihood
#' via numerical differentiation using \code{numDeriv::hessian}.
#'
#' @param model A \code{wei_kofn} object.
#' @param ... Additional arguments (currently unused).
#' @return A function \code{function(df, par)} returning a square matrix
#'   (the Hessian of the log-likelihood).
#'
#' @importFrom numDeriv hessian
#' @export
hess_loglik.wei_kofn <- function(model, ...) {
  ll_fn <- loglik(model, ...)
  function(df, par) {
    numDeriv::hessian(func = function(p) ll_fn(df, p), x = par,
                      method = "Richardson")
  }
}


# ---------------------------------------------------------------------------
# fit.wei_kofn
# ---------------------------------------------------------------------------

#' Fit Weibull k-out-of-n system model
#'
#' Returns a closure that fits the model to data. Two methods are available:
#'
#' \describe{
#'   \item{\code{"em"}}{EM algorithm for parallel systems (k = 1 only).
#'     Treats the identity of the last-failing component as latent data.
#'     Uses truncated Weibull moments for the E-step and profile
#'     optimization over shape with closed-form scale for the M-step.}
#'   \item{\code{"mle"}}{Direct MLE via L-BFGS-B with Nelder-Mead fallback.
#'     Works for any k-out-of-n structure.}
#' }
#'
#' @param model A \code{wei_kofn} object.
#' @param ... Additional arguments (currently unused).
#' @return A function \code{function(df, par0, ...)} that returns a
#'   \code{fisher_mle} object (from the likelihood.model package).
#'
#' @details
#' The returned solver accepts these additional arguments:
#' \describe{
#'   \item{\code{par0}}{Initial parameter vector (length 2*m, interleaved
#'     shape/scale). If \code{NULL}, method-of-moments initialization is used.}
#'   \item{\code{n_starts}}{Number of multi-start attempts (default: 5).}
#'   \item{\code{tol}}{Convergence tolerance for EM (default: 1e-8).}
#'   \item{\code{maxiter}}{Maximum EM iterations (default: 2000).}
#'   \item{\code{shape_bounds}}{Bounds for shape optimization in EM (default:
#'     \code{c(0.05, 15)}).}
#'   \item{\code{verbose}}{Print iteration progress (default: FALSE).}
#' }
#'
#' @export
fit.wei_kofn <- function(model, ...) {
  m <- model$m
  k <- model$k
  method <- model$method
  system <- model$system
  lt <- model$lifetime
  ll_fn <- loglik(model, ...)
  n_par <- 2L * m

  if (method == "em") {
    # EM algorithm (parallel systems only)
    if (k != 1L) {
      stop("EM algorithm is only available for parallel systems (k = 1)")
    }
    em_solver(model, ll_fn, ...)
  } else {
    # Direct MLE
    mle_solver(model, ll_fn, ...)
  }
}


# ---------------------------------------------------------------------------
# Internal: EM solver for parallel Weibull systems
# ---------------------------------------------------------------------------

#' @keywords internal
em_solver <- function(model, ll_fn, ...) {
  m <- model$m
  lt <- model$lifetime
  n_par <- 2L * m

  function(df, par0 = NULL,
           n_starts = 5L,
           tol = 1e-8,
           maxiter = 2000L,
           shape_bounds = c(0.05, 15),
           verbose = FALSE,
           ...) {
    t_obs <- df[[lt]]
    n <- length(t_obs)
    stopifnot(n > 0, all(t_obs > 0))

    # --- Default initialization ---
    if (is.null(par0)) {
      mean_t <- mean(t_obs)
      H_m <- sum(1 / seq_len(m))
      scale0 <- mean_t / H_m
      init_shapes <- rep(1.0, m)
      init_scales <- scale0 * seq(0.5, 1.5, length.out = m)
    } else {
      stopifnot(length(par0) == n_par)
      init_shapes <- wei_par_to_shapes(par0, m)
      init_scales <- wei_par_to_scales(par0, m)
    }
    stopifnot(all(init_shapes > 0), all(init_scales > 0))

    # --- Log-likelihood for Weibull parallel system ---
    loglik_wei <- function(shapes, scales) {
      par_vec <- as.numeric(rbind(shapes, scales))
      ll_fn(df, par_vec)
    }

    # --- Single EM run from given initialization ---
    em_run <- function(shapes, scales) {
      ll_prev <- loglik_wei(shapes, scales)
      trace <- data.frame(iter = 0L, loglik = ll_prev)
      log_t <- log(t_obs)

      for (iter in seq_len(maxiter)) {
        # ========== E-step ==========
        # Precompute v_max_ij = (t_i / scales[j])^shapes[j]
        v_max_mat <- matrix(0, nrow = n, ncol = m)
        for (j in seq_len(m)) {
          v_max_mat[, j] <- (t_obs / scales[j])^shapes[j]
        }

        # Posterior p_ij = P(J=j | T=t_i, theta_old)
        p_mat <- matrix(0, nrow = n, ncol = m)
        for (i in seq_len(n)) {
          wj <- vapply(seq_len(m), function(j) {
            weibull_w_j(t_obs[i], shapes, scales, j)
          }, numeric(1))
          f_sys <- sum(wj)
          if (f_sys > 0) {
            p_mat[i, ] <- wj / f_sys
          } else {
            p_mat[i, ] <- rep(1 / m, m)
          }
        }

        # E[log T_j | T_j < t_i] for each component j (vectorized over i)
        elog_mat <- matrix(0, nrow = n, ncol = m)
        for (j in seq_len(m)) {
          elog_mat[, j] <- trunc_log_moment_vec(v_max_mat[, j], shapes[j],
                                                scales[j])
        }

        # A_j = sum_i [p_ij * log(t_i) + (1 - p_ij) * E[log T_j | T_j < t_i]]
        A <- colSums(p_mat * log_t + (1 - p_mat) * elog_mat)

        # ========== M-step ==========
        new_shapes <- numeric(m)
        new_scales <- numeric(m)

        for (j in seq_len(m)) {
          vm_j <- v_max_mat[, j]
          p_j <- p_mat[, j]

          # Profile Q_j(alpha):
          # = n*log(alpha) + (alpha-1)*A_j - n - n*log(B_j(alpha)/n)
          profile_Q <- function(alpha) {
            if (alpha <= 0) return(-Inf)
            t_pow <- t_obs^alpha
            epow <- trunc_pow_moment_vec(alpha, vm_j, shapes[j], scales[j])
            B <- sum(p_j * t_pow + (1 - p_j) * epow)
            if (B <= 0 || !is.finite(B)) return(-Inf)
            val <- n * log(alpha) + (alpha - 1) * A[j] - n - n * log(B / n)
            if (!is.finite(val)) -Inf else val
          }

          opt <- stats::optimize(profile_Q, interval = shape_bounds,
                                 maximum = TRUE)
          new_shapes[j] <- opt$maximum

          # Recover scale: beta^alpha = B(alpha) / n
          alpha_opt <- new_shapes[j]
          t_pow <- t_obs^alpha_opt
          epow <- trunc_pow_moment_vec(alpha_opt, vm_j, shapes[j], scales[j])
          B_opt <- sum(p_j * t_pow + (1 - p_j) * epow)
          new_scales[j] <- (B_opt / n)^(1 / alpha_opt)

          # Guard degenerate
          new_shapes[j] <- max(new_shapes[j], shape_bounds[1])
          new_scales[j] <- max(new_scales[j], 1e-10)
        }

        # ========== Convergence ==========
        ll_new <- loglik_wei(new_shapes, new_scales)
        trace <- rbind(trace, data.frame(iter = iter, loglik = ll_new))

        if (verbose) {
          message(sprintf("  EM iter %d: loglik = %.6f (delta = %.2e)",
                          iter, ll_new, ll_new - ll_prev))
        }

        if (abs(ll_new - ll_prev) / max(abs(ll_prev), 1) < tol) {
          shapes <- new_shapes
          scales <- new_scales
          ll_prev <- ll_new
          break
        }

        shapes <- new_shapes
        scales <- new_scales
        ll_prev <- ll_new
      }

      list(shapes = shapes, scales = scales, loglik = ll_prev,
           converged = iter < maxiter, iterations = iter, trace = trace)
    }

    # --- Multi-start ---
    best <- em_run(init_shapes, init_scales)

    if (n_starts > 1L) {
      for (s in seq_len(n_starts - 1L)) {
        # Random perturbation of initialization
        sh_s <- init_shapes * exp(stats::rnorm(m, sd = 0.3))
        sc_s <- init_scales * exp(stats::rnorm(m, sd = 0.5))
        sh_s <- pmax(sh_s, shape_bounds[1])
        sc_s <- pmax(sc_s, 1e-10)

        res_s <- em_run(sh_s, sc_s)
        if (res_s$converged &&
            (!best$converged || res_s$loglik > best$loglik)) {
          best <- res_s
        }
      }
    }

    # --- Standard errors via numerical Hessian ---
    par_est <- as.numeric(rbind(best$shapes, best$scales))

    neg_ll_for_hess <- function(p) {
      val <- -ll_fn(df, p)
      if (!is.finite(val)) return(.Machine$double.xmax / 2)
      val
    }
    H <- tryCatch(
      numDeriv::hessian(neg_ll_for_hess, x = par_est, method = "Richardson"),
      error = function(e) NULL
    )

    hessian_mat <- if (!is.null(H) && all(is.finite(H))) H else {
      matrix(NA_real_, nrow = n_par, ncol = n_par)
    }
    score_val <- tryCatch(
      numDeriv::grad(func = function(p) ll_fn(df, p), x = par_est,
                     method = "Richardson"),
      error = function(e) rep(NA_real_, n_par)
    )

    make_mle_result(
      par        = par_est,
      loglik     = best$loglik,
      hessian    = hessian_mat,
      score      = score_val,
      nobs       = n,
      converged  = best$converged,
      iterations = best$iterations,
      shapes     = par_est[seq(1L, n_par, by = 2L)],
      scales     = par_est[seq(2L, n_par, by = 2L)]
    )
  }
}


# ---------------------------------------------------------------------------
# Internal: direct MLE solver
# ---------------------------------------------------------------------------

#' @keywords internal
mle_solver <- function(model, ll_fn, ...) {
  m <- model$m
  system <- model$system
  lt <- model$lifetime
  n_par <- 2L * m

  function(df, par0 = NULL,
           n_starts = 5L,
           ...) {
    t_obs <- df[[lt]]
    n <- length(t_obs)

    # Default initialization
    if (is.null(par0)) {
      mean_t <- max(mean(t_obs, na.rm = TRUE), 0.01)
      scale0 <- mean_t * seq(0.5, 1.5, length.out = m)
      par0 <- as.numeric(rbind(rep(1.0, m), scale0))
    }
    stopifnot(length(par0) == n_par)

    neg_ll <- function(p) {
      val <- -ll_fn(df, p)
      if (!is.finite(val)) return(.Machine$double.xmax / 2)
      val
    }

    neg_ll_log <- function(log_p) neg_ll(exp(log_p))

    fit_from_init <- function(init_par) {
      result <- tryCatch(
        stats::optim(par = init_par, fn = neg_ll, method = "L-BFGS-B",
                     lower = rep(1e-10, n_par), upper = rep(Inf, n_par)),
        error = function(e) {
          list(par = rep(NA_real_, n_par), value = Inf, convergence = 99L)
        }
      )
      if (result$convergence != 0) {
        result2 <- tryCatch(
          stats::optim(par = log(pmax(init_par, 1e-10)), fn = neg_ll_log,
                       method = "Nelder-Mead"),
          error = function(e) {
            list(par = rep(NA_real_, n_par), value = Inf, convergence = 99L)
          }
        )
        if (result2$convergence == 0 && result2$value < result$value) {
          result <- list(par = exp(result2$par), value = result2$value,
                         convergence = result2$convergence)
        }
      }
      result
    }

    # Multi-start
    best <- fit_from_init(par0)
    if (n_starts > 1L) {
      for (s in seq_len(n_starts - 1L)) {
        init_s <- par0 * exp(stats::rnorm(n_par, sd = 0.5))
        init_s <- pmax(init_s, 1e-10)
        res_s <- fit_from_init(init_s)
        if (res_s$convergence == 0 &&
            (best$convergence != 0 || res_s$value < best$value)) {
          best <- res_s
        }
      }
    }

    # Hessian of neg log-likelihood and score at MLE
    H <- tryCatch(
      numDeriv::hessian(neg_ll, x = best$par, method = "Richardson"),
      error = function(e) NULL
    )
    hessian_mat <- if (!is.null(H) && all(is.finite(H))) H else {
      matrix(NA_real_, nrow = n_par, ncol = n_par)
    }
    score_val <- tryCatch(
      numDeriv::grad(func = function(p) ll_fn(df, p), x = best$par,
                     method = "Richardson"),
      error = function(e) rep(NA_real_, n_par)
    )

    make_mle_result(
      par        = best$par,
      loglik     = -best$value,
      hessian    = hessian_mat,
      score      = score_val,
      nobs       = n,
      converged  = (best$convergence == 0),
      shapes     = best$par[seq(1L, n_par, by = 2L)],
      scales     = best$par[seq(2L, n_par, by = 2L)]
    )
  }
}


# ---------------------------------------------------------------------------
# rdata.wei_kofn
# ---------------------------------------------------------------------------

#' Generate random data from a Weibull k-out-of-n system
#'
#' Returns a closure \code{function(theta, n, observe = NULL)} that generates
#' system lifetime data. Component lifetimes are Weibull-distributed with
#' parameters specified by \code{theta} (interleaved shape/scale).
#'
#' @param model A \code{wei_kofn} object.
#' @param ... Additional arguments (currently unused).
#' @return A function \code{function(theta, n, observe = NULL)} returning a
#'   data frame with column \code{t} (system lifetimes). The data frame has
#'   attributes:
#'   \describe{
#'     \item{\code{comp_times}}{n x m matrix of component lifetimes.}
#'     \item{\code{par}}{The parameter vector used for generation.}
#'   }
#'   If \code{observe} is provided (a function), it is applied to the data
#'   frame to add observation-scheme-specific columns (e.g., censoring).
#'
#' @export
rdata.wei_kofn <- function(model, ...) {
  m <- model$m
  k <- model$k
  system <- model$system
  lt <- model$lifetime

  function(theta, n, observe = NULL) {
    stopifnot(length(theta) == 2L * m)
    shapes <- wei_par_to_shapes(theta, m)
    scales <- wei_par_to_scales(theta, m)
    stopifnot(all(shapes > 0), all(scales > 0))

    # Generate component lifetimes
    comp_times <- matrix(0, nrow = n, ncol = m)
    for (j in seq_len(m)) {
      comp_times[, j] <- stats::rweibull(n, shape = shapes[j],
                                         scale = scales[j])
    }

    # Compute system lifetimes
    if (k == 1L) {
      # Parallel: T = max
      sys_times <- apply(comp_times, 1, max)
    } else {
      # General k-out-of-n: use system_lifetime via cut sets
      sys_times <- vapply(seq_len(n), function(i) {
        system_lifetime(system, comp_times[i, ])
      }, numeric(1))
    }

    df <- data.frame(sys_times)
    names(df) <- lt
    attr(df, "comp_times") <- comp_times
    attr(df, "par") <- theta

    if (!is.null(observe)) {
      df <- observe(df)
    }

    df
  }
}


# ---------------------------------------------------------------------------
# assumptions.wei_kofn
# ---------------------------------------------------------------------------

#' Assumptions for Weibull k-out-of-n system model
#'
#' Returns a list of assumptions made by the model.
#'
#' @param model A \code{wei_kofn} object.
#' @param ... Additional arguments (currently unused).
#' @return A character vector of assumptions.
#'
#' @export
assumptions.wei_kofn <- function(model, ...) {
  m <- model$m
  k <- model$k
  method <- model$method

  assumptions <- c(
    sprintf("System is %d-out-of-%d (system fails when %d components fail)",
            k, m, k),
    sprintf("Component count: m = %d", m),
    "Component lifetimes are independent",
    "Each component lifetime follows a Weibull distribution",
    sprintf("Parameter vector: %d parameters (shape_j, scale_j for j=1,...,%d)",
            2L * m, m),
    "Observation: system lifetime only (Scheme 0)"
  )

  if (method == "em" && k == 1L) {
    assumptions <- c(assumptions,
      "Estimation: EM algorithm with latent last-failing component",
      "E-step uses truncated Weibull moments",
      "M-step uses profile optimization over shape, closed-form scale"
    )
  } else {
    assumptions <- c(assumptions,
      "Estimation: direct MLE via L-BFGS-B with Nelder-Mead fallback"
    )
  }

  assumptions
}
