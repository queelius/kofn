# ===========================================================================
# Re-exports from likelihood.model and generics
# ===========================================================================
#
# These generics follow the likelihood.model interface convention:
# each returns a closure that can be evaluated at data + parameters.
# ===========================================================================

#' @importFrom likelihood.model loglik
#' @export
likelihood.model::loglik

#' @importFrom likelihood.model score
#' @export
likelihood.model::score

#' @importFrom likelihood.model hess_loglik
#' @export
likelihood.model::hess_loglik

#' @importFrom likelihood.model assumptions
#' @export
likelihood.model::assumptions

#' @importFrom likelihood.model rdata
#' @export
likelihood.model::rdata

#' @importFrom generics fit
#' @export
generics::fit

# ===========================================================================
# Re-exports from flexhaz: component-distribution prototypes
# ===========================================================================
#
# Users construct kofn models via `component = dfr_exponential()` or
# `component = dfr_weibull()`. Re-exporting these makes them callable
# without an explicit `library(flexhaz)` in user code and vignettes.
# ===========================================================================

#' @importFrom flexhaz dfr_exponential
#' @export
flexhaz::dfr_exponential

#' @importFrom flexhaz dfr_weibull
#' @export
flexhaz::dfr_weibull

# ===========================================================================
# Re-exports from dist.structure
# ===========================================================================

#' @importFrom dist.structure ncomponents
#' @export
dist.structure::ncomponents
