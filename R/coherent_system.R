# ===========================================================================
# Coherent System Engine
# ===========================================================================
#
# Representation and algorithms for coherent systems via minimal path sets.
#
# A coherent system is defined by its minimal path sets -- minimal sets of
# components whose simultaneous functioning guarantees system functioning.
#
# Dual to path sets are minimal cut sets -- minimal sets of components whose
# simultaneous failure causes system failure.  These are computed via the
# Berge transversal (hitting set) algorithm.
# ===========================================================================


# ---------------------------------------------------------------------------
# Helper: set operations on lists of integer vectors
# ---------------------------------------------------------------------------

#' Remove non-minimal sets from a list
#'
#' Given a list of integer vectors, removes any set that is a proper superset
#' of another set in the list. Used internally for minimizing path and cut
#' set collections.
#'
#' @param sets List of integer vectors.
#' @return List of integer vectors containing only minimal elements.
#' @export
minimize_sets <- function(sets) {
  n <- length(sets)
  if (n == 0L) return(sets)

  is_minimal <- rep(TRUE, n)
  for (i in seq_len(n)) {
    if (!is_minimal[i]) next
    for (j in seq_len(n)) {
      if (i == j || !is_minimal[j]) next
      # if sets[[i]] is a proper subset of sets[[j]], remove sets[[j]]
      if (all(sets[[i]] %in% sets[[j]]) && length(sets[[i]]) < length(sets[[j]])) {
        is_minimal[j] <- FALSE
      }
    }
  }
  sets[is_minimal]
}


# ---------------------------------------------------------------------------
# Minimal cut sets via Berge transversal algorithm
# ---------------------------------------------------------------------------

#' Compute minimal cut sets from minimal path sets
#'
#' A cut set must hit (intersect) every minimal path set. The minimal
#' transversals of the path set hypergraph are exactly the minimal cut sets.
#'
#' Algorithm: iterative construction -- start with the empty set and for each
#' path set P, replace each partial transversal T with
#' \eqn{\{T \cup \{p\} : p \in P\}}{T union {p} : p in P}, then minimize.
#'
#' @param min_paths List of integer vectors representing minimal path sets.
#' @param m Number of components (optional; not used but accepted for
#'   compatibility).
#' @return List of integer vectors representing minimal cut sets.
#' @export
min_cuts_from_paths <- function(min_paths, m = NULL) {
  if (length(min_paths) == 0L) return(list())

  # Initialize: single empty transversal
  transversals <- list(integer(0))

  for (path in min_paths) {
    new_transversals <- list()
    for (T in transversals) {
      # Extend T by adding each element of path
      for (p in path) {
        new_T <- sort(unique(c(T, p)))
        new_transversals <- c(new_transversals, list(new_T))
      }
    }
    # Remove duplicates then non-minimals
    # Deduplicate by canonical string key
    keys <- vapply(new_transversals, function(s) paste(sort(s), collapse = "-"),
                   character(1))
    new_transversals <- new_transversals[!duplicated(keys)]
    transversals <- minimize_sets(new_transversals)
  }

  transversals
}


# ---------------------------------------------------------------------------
# S3 constructor: coherent_system
# ---------------------------------------------------------------------------

#' Create a coherent system from minimal path sets
#'
#' A coherent system is defined by its minimal path sets. The constructor
#' precomputes minimal cut sets eagerly using the Berge transversal algorithm.
#'
#' @param min_paths List of integer vectors, where each vector contains
#'   1-based component indices forming a minimal path set.
#' @param m Number of components. If `NULL`, inferred as the maximum index
#'   appearing in `min_paths`.
#' @return An object of class `"coherent_system"` with elements:
#'   \describe{
#'     \item{min_paths}{List of integer vectors (minimal path sets).}
#'     \item{min_cuts}{List of integer vectors (minimal cut sets).}
#'     \item{m}{Integer number of components.}
#'   }
#' @examples
#' # A 2-out-of-3 system: paths are all pairs
#' sys <- coherent_system(list(c(1,2), c(1,3), c(2,3)), m = 3)
#' sys
#'
#' @export
coherent_system <- function(min_paths, m = NULL) {
  # Validate and normalise paths
  if (!is.list(min_paths) || length(min_paths) == 0L) {
    stop("min_paths must be a non-empty list of integer vectors")
  }
  min_paths <- lapply(min_paths, function(p) sort(as.integer(p)))

  # Infer m
  if (is.null(m)) {
    m <- max(unlist(min_paths))
  }
  m <- as.integer(m)

  # Precompute minimal cut sets
  cuts <- min_cuts_from_paths(min_paths, m)

  structure(
    list(
      min_paths = min_paths,
      min_cuts  = cuts,
      m         = m
    ),
    class = "coherent_system"
  )
}


#' Print method for coherent_system objects
#'
#' @param x A `coherent_system` object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns `x`.
#' @export
print.coherent_system <- function(x, ...) {
  cat(sprintf("Coherent system: m=%d components\n", x$m))
  cat(sprintf("  Minimal path sets (%d):\n", length(x$min_paths)))
  for (p in x$min_paths) {
    cat(sprintf("    {%s}\n", paste(p, collapse = ", ")))
  }
  cat(sprintf("  Minimal cut sets (%d):\n", length(x$min_cuts)))
  for (cs in x$min_cuts) {
    cat(sprintf("    {%s}\n", paste(cs, collapse = ", ")))
  }
  invisible(x)
}


# ---------------------------------------------------------------------------
# Standard system constructors
# ---------------------------------------------------------------------------

#' Create a series system
#'
#' A series system has a single minimal path set containing all components:
#' \eqn{\{\{1, 2, \ldots, m\}\}}{{1, 2, ..., m}}. The system functions only
#' when all components function. Equivalently, it is an m-out-of-m system.
#'
#' @param m Number of components (positive integer).
#' @return A `coherent_system` object.
#' @examples
#' sys <- series_system(3)
#' phi(sys, c(TRUE, TRUE, TRUE))   # TRUE
#' phi(sys, c(TRUE, FALSE, TRUE))  # FALSE
#'
#' @export
series_system <- function(m) {
  coherent_system(list(seq_len(m)), m = m)
}


#' Create a parallel system
#'
#' A parallel system has m singleton minimal path sets:
#' \eqn{\{\{1\}, \{2\}, \ldots, \{m\}\}}{{1}, {2}, ..., {m}}. The system
#' functions when at least one component functions. Equivalently, it is a
#' 1-out-of-m system.
#'
#' @param m Number of components (positive integer).
#' @return A `coherent_system` object.
#' @examples
#' sys <- parallel_system(3)
#' phi(sys, c(FALSE, FALSE, TRUE))  # TRUE
#' phi(sys, c(FALSE, FALSE, FALSE)) # FALSE
#'
#' @export
parallel_system <- function(m) {
  coherent_system(as.list(seq_len(m)), m = m)
}


#' Create a k-out-of-n system
#'
#' A k-out-of-n system functions when at least k of its n components
#' function. Minimal path sets are all subsets of size k.
#'
#' Special cases: k=1 gives a parallel system, k=n gives a series system.
#'
#' @param k Minimum number of functioning components (positive integer,
#'   \eqn{1 \le k \le n}{1 <= k <= n}).
#' @param n Total number of components (positive integer).
#' @return A `coherent_system` object.
#' @examples
#' sys <- kofn_system(2, 3)  # 2-out-of-3
#' phi(sys, c(TRUE, TRUE, FALSE))   # TRUE
#' phi(sys, c(TRUE, FALSE, FALSE))  # FALSE
#'
#' @export
kofn_system <- function(k, n) {
  stopifnot(k >= 1L, k <= n, n >= 1L)
  paths <- utils::combn(n, k, simplify = FALSE)
  coherent_system(paths, m = n)
}


#' Create a consecutive-k-out-of-n system (linear)
#'
#' A consecutive-k-out-of-n system fails when k consecutive components all
#' fail. Cut sets are all runs of k consecutive indices:
#' \eqn{\{i, i+1, \ldots, i+k-1\}}{i, i+1, ..., i+k-1}. Paths are computed
#' via cut-to-path duality using the Berge transversal algorithm.
#'
#' @param k Run length for failure (positive integer, \eqn{1 \le k \le n}).
#' @param n Total number of components (positive integer).
#' @return A `coherent_system` object.
#' @examples
#' sys <- consecutive_k_system(2, 4)
#' phi(sys, c(FALSE, FALSE, TRUE, TRUE))  # FALSE (components 1,2 consecutive fail)
#' phi(sys, c(FALSE, TRUE, FALSE, TRUE))  # TRUE  (no two consecutive failures)
#'
#' @export
consecutive_k_system <- function(k, n) {
  stopifnot(k >= 1L, k <= n)
  cuts <- lapply(seq_len(n - k + 1L), function(i) seq(i, i + k - 1L))
  # Convert cuts to paths via transversal (dual: transversals of cuts = paths)
  paths <- min_cuts_from_paths(cuts)
  coherent_system(paths, m = n)
}


#' Create a bridge system
#'
#' The standard bridge (reliability) system with 5 components. Minimal path
#' sets are \eqn{\{1,4\}, \{2,5\}, \{1,3,5\}, \{2,3,4\}}.
#'
#' @return A `coherent_system` object with `m = 5`.
#' @examples
#' sys <- bridge_system()
#' sys
#'
#' @export
bridge_system <- function() {
  paths <- list(
    c(1L, 4L),
    c(2L, 5L),
    c(1L, 3L, 5L),
    c(2L, 3L, 4L)
  )
  coherent_system(paths, m = 5L)
}


# ---------------------------------------------------------------------------
# Accessor: minimal cut sets
# ---------------------------------------------------------------------------

#' Get minimal cut sets from a coherent system
#'
#' @param system A `coherent_system` object.
#' @return List of integer vectors representing minimal cut sets.
#' @export
min_cuts <- function(system) {
  stopifnot(inherits(system, "coherent_system"))
  system$min_cuts
}


# ---------------------------------------------------------------------------
# Structure function phi
# ---------------------------------------------------------------------------

#' Evaluate the structure function of a coherent system
#'
#' Returns `TRUE` if the system functions given a logical component state
#' vector. Uses the path representation: the system works if and only if at
#' least one minimal path set is fully operational.
#'
#' @param system A `coherent_system` object.
#' @param x Logical (or coercible to logical) vector of length `m`, where
#'   `x[j] = TRUE` means component j is functioning.
#' @return Logical scalar: `TRUE` if the system functions, `FALSE` otherwise.
#' @examples
#' sys <- parallel_system(3)
#' phi(sys, c(TRUE, FALSE, FALSE))   # TRUE
#' phi(sys, c(FALSE, FALSE, FALSE))  # FALSE
#'
#' @export
phi <- function(system, x) {
  stopifnot(inherits(system, "coherent_system"))
  x <- as.logical(x)
  for (path in system$min_paths) {
    if (all(x[path])) return(TRUE)
  }
  FALSE
}


# ---------------------------------------------------------------------------
# System lifetime
# ---------------------------------------------------------------------------

#' Compute the system lifetime from component lifetimes
#'
#' For a coherent system with minimal cut sets \eqn{C_1, \ldots, C_r},
#' the system lifetime is:
#' \deqn{T_{sys} = \min_{C} \max_{j \in C} T_j}
#'
#' This is the series-of-parallel (min-of-max) representation dual to the
#' parallel-of-series path representation.
#'
#' @param system A `coherent_system` object.
#' @param times Numeric vector of component lifetimes (length `m`).
#' @return Scalar system lifetime.
#' @examples
#' sys <- parallel_system(3)
#' system_lifetime(sys, c(1.5, 2.3, 0.8))  # max = 2.3
#'
#' @export
system_lifetime <- function(system, times) {
  stopifnot(inherits(system, "coherent_system"))
  cuts <- system$min_cuts
  if (length(cuts) == 0L) stop("No minimal cut sets -- degenerate system?")
  min(vapply(cuts, function(C) max(times[C]), numeric(1)))
}


# ---------------------------------------------------------------------------
# Censoring pattern
# ---------------------------------------------------------------------------

#' Determine per-component censoring status from a coherent system
#'
#' Given component lifetimes and a coherent system, identifies:
#' \describe{
#'   \item{T_sys}{The system lifetime.}
#'   \item{critical}{Index of the critical component (the one in the binding
#'     cut set whose maximum equals `T_sys`; ties broken by lowest index).}
#'   \item{status}{Per-component censoring status:
#'     `"exact"` if \eqn{T_j = T_{sys}},
#'     `"left"` if \eqn{T_j < T_{sys}} (component failed before system),
#'     `"right"` if \eqn{T_j > T_{sys}} (component still alive at system
#'     failure).}
#' }
#'
#' @param system A `coherent_system` object.
#' @param times Numeric vector of component lifetimes (length `m`).
#' @return A list with elements `T_sys`, `critical`, and `status`.
#' @examples
#' sys <- parallel_system(2)
#' system_censoring(sys, c(1.5, 2.3))
#' # T_sys = 2.3, critical = 2, status = c("left", "exact")
#'
#' @export
system_censoring <- function(system, times) {
  stopifnot(inherits(system, "coherent_system"))
  m <- system$m
  T_sys <- system_lifetime(system, times)

  # Find the binding cut set: the one achieving T_sys
  cuts <- system$min_cuts
  cut_maxima <- vapply(cuts, function(C) max(times[C]), numeric(1))
  binding_idx <- which(cut_maxima == T_sys)[1L]
  binding_cut <- cuts[[binding_idx]]

  # Within the binding cut, the critical component is the one with T_j = T_sys
  crit_candidates <- binding_cut[times[binding_cut] == T_sys]
  critical <- crit_candidates[1L]  # break ties by taking first (lowest index)

  # Assign per-component status (vectorized)
  tol <- .Machine$double.eps^0.5
  status <- ifelse(
    abs(times - T_sys) < tol * max(1, T_sys), "exact",
    ifelse(times < T_sys, "left", "right")
  )

  list(T_sys = T_sys, critical = critical, status = status)
}


# ---------------------------------------------------------------------------
# Critical states
# ---------------------------------------------------------------------------

#' Find critical states for a component
#'
#' Returns a matrix where each row is a state vector \eqn{x_{-j}} (components
#' other than j) such that component j is critical:
#' \deqn{\phi(x_{-j}, x_j=1) = 1 \quad \text{and} \quad \phi(x_{-j}, x_j=0) = 0}
#'
#' That is, j is the pivotal component -- flipping its state changes the
#' system state.
#'
#' @param system A `coherent_system` object.
#' @param j Component index (integer, 1-based).
#' @return A logical matrix with \eqn{m-1} columns (one per component other
#'   than j). Each row is a state vector making j critical. Column names are
#'   `"comp1"`, `"comp2"`, etc. (excluding j).
#' @export
critical_states <- function(system, j) {
  stopifnot(inherits(system, "coherent_system"))
  m <- system$m
  others <- setdiff(seq_len(m), j)

  # Edge case: m=1, no other components -- j is trivially critical
  if (length(others) == 0L) {
    return(matrix(integer(0), nrow = 1L, ncol = 0L))
  }

  # Enumerate all 2^(m-1) states for components other than j
  n_other <- length(others)
  n_states <- 2L^n_other
  result_rows <- list()

  for (k in seq_len(n_states) - 1L) {
    bits <- as.logical(as.integer(intToBits(k))[seq_len(n_other)])
    x <- logical(m)
    x[others] <- bits

    # j is critical if flipping it changes the system state
    x[j] <- TRUE
    sys_on <- phi(system, x)
    x[j] <- FALSE
    if (sys_on && !phi(system, x)) {
      result_rows <- c(result_rows, list(bits))
    }
  }

  if (length(result_rows) == 0L) {
    mat <- matrix(logical(0), nrow = 0L, ncol = n_other)
  } else {
    mat <- do.call(rbind, result_rows)
  }
  colnames(mat) <- paste0("comp", others)
  mat
}


# ---------------------------------------------------------------------------
# System signature
# ---------------------------------------------------------------------------

#' Generate all permutations of 1:m
#'
#' Returns a matrix with \eqn{m!} rows and \eqn{m} columns, where each row
#' is a permutation of `1:m`.
#'
#' @param m Positive integer.
#' @return Integer matrix of dimension \eqn{m! \times m}{m! x m}.
#' @keywords internal
permutations <- function(m) {
  if (m == 1L) return(matrix(1L, nrow = 1L, ncol = 1L))

  sub_perms <- permutations(m - 1L)
  n_sub <- nrow(sub_perms)
  result <- matrix(0L, nrow = m * n_sub, ncol = m)
  row_idx <- 1L

  for (pos in seq_len(m)) {
    for (i in seq_len(n_sub)) {
      perm <- integer(m)
      perm[pos] <- m
      rest <- sub_perms[i, ]
      other_pos <- seq_len(m)[-pos]
      perm[other_pos] <- rest
      result[row_idx, ] <- perm
      row_idx <- row_idx + 1L
    }
  }
  result
}


#' Compute the system signature
#'
#' The system signature \eqn{s} is a probability vector of length \eqn{m}
#' where \eqn{s_k = P(T_{sys} = T_{(k)})}{s_k = P(T_sys = T_(k))} under the
#' assumption that component lifetimes are i.i.d. continuous.
#'
#' Computed by exhaustive enumeration over all \eqn{m!} permutations of
#' component lifetime orderings. Practical only for moderate \eqn{m}
#' (say \eqn{m \le 10}).
#'
#' @param system A `coherent_system` object.
#' @return Numeric vector of length `m` summing to 1.
#' @examples
#' system_signature(series_system(3))    # c(1, 0, 0)
#' system_signature(parallel_system(3))  # c(0, 0, 1)
#' system_signature(kofn_system(2, 3))   # c(0, 2/3, 1/3)
#'
#' @export
system_signature <- function(system) {
  stopifnot(inherits(system, "coherent_system"))
  m <- system$m
  perms <- permutations(m)
  n_perms <- nrow(perms)
  s <- numeric(m)

  for (i in seq_len(n_perms)) {
    perm <- perms[i, ]
    # Assign fictitious times by rank: component perm[k] has lifetime k
    times <- numeric(m)
    times[perm] <- seq_len(m)
    # Compute system lifetime order statistic index
    T_sys <- system_lifetime(system, times)
    k <- as.integer(round(T_sys))
    s[k] <- s[k] + 1L
  }

  s / n_perms
}
