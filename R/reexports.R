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
