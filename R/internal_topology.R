# ===========================================================================
# Internal topology helpers (k-out-of-n specialized)
# ===========================================================================
#
# kofn delegates DGP and topology queries to dist.structure, but several
# downstream functions (rdata.exp_kofn, rdata.wei_kofn, masked, scheme1,
# fisher_info) need:
#
#   1. A k-of-n system_censoring that returns kofn's historical
#      (T_sys, critical, status) shape, where `critical` is the index of
#      the "binding" component (the component whose failure caused the
#      system failure under the k-of-n rule). dist.structure's analog
#      returns (system_time, component_status), without identifying the
#      critical component.
#
#   2. A small wrapper to construct dist.structure exp_kofn / wei_kofn
#      objects from the kofn model's (k, family, par) without re-deriving
#      shapes or rates.
#
# These are package-private helpers (not exported).
# ===========================================================================


# k-of-m system censoring (kofn :F convention).
#
# In the :F convention, k is the number of component failures that
# trigger system failure: k=1 is series (first failure ends it), k=m is
# parallel (system survives until the m-th failure). System lifetime is
# therefore the k-th order statistic of the component lifetimes.
#
# Returns:
#   T_sys: scalar system lifetime (the k-th order statistic).
#   critical: index of the component whose failure equals T_sys
#     (deterministic for absolutely continuous components).
#   status: per-component label in c("exact", "left", "right").
kofn_censoring <- function(k, times) {
  m <- length(times)
  stopifnot(k >= 1L, k <= m)
  ord <- order(times)
  T_sys <- times[ord[k]]
  critical <- ord[k]
  status <- ifelse(times < T_sys, "left",
                   ifelse(times > T_sys, "right", "exact"))
  list(T_sys = T_sys, critical = critical, status = status)
}


# Construct a dist.structure k-of-n DGP from a kofn model + parameter vec.
#
# Convention conversion: kofn uses the :F convention (k_kofn = number of
# component failures that trigger system failure), while dist.structure
# uses the :G convention (k_dist = number of functioning components
# required for the system to function). They convert via
#   k_dist = m - k_kofn + 1
# Examples:
#   kofn series  (k_kofn = 1, fails on first failure) -> k_dist = m (all required)
#   kofn parallel (k_kofn = m, fails on last failure) -> k_dist = 1 (any one required)
kofn_dgp <- function(k, m, family, par) {
  k_dist <- m - k + 1L
  if (family == "exponential") {
    if (length(par) != m) {
      stop(sprintf("Expected %d rates but got %d", m, length(par)))
    }
    return(dist.structure::exp_kofn(k_dist, par))
  }
  # weibull: interleaved (shape_1, scale_1, shape_2, scale_2, ...)
  if (length(par) != 2L * m) {
    stop(sprintf("Expected %d Weibull parameters but got %d",
                 2L * m, length(par)))
  }
  shapes <- par[seq(1L, 2L * m, by = 2L)]
  scales <- par[seq(2L, 2L * m, by = 2L)]
  dist.structure::wei_kofn(k_dist, shapes, scales)
}


# Per-component distribution constructors at given parameters.
# Returns a list of m algebraic.dist `dist` objects with $surv and $cdf
# (closures, not the kofn dist.R lightweight $surv/$cdf interface).
# Used by masked.R and scheme1.R for component-level CDF/survival.
kofn_components <- function(k, m, family, par) {
  if (family == "exponential") {
    return(lapply(par, function(r) algebraic.dist::exponential(r)))
  }
  shapes <- par[seq(1L, 2L * m, by = 2L)]
  scales <- par[seq(2L, 2L * m, by = 2L)]
  lapply(seq_len(m), function(j) {
    algebraic.dist::weibull_dist(shape = shapes[j], scale = scales[j])
  })
}
