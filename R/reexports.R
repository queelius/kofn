# ===========================================================================
# S3 Generics for kofn
# ===========================================================================
#
# These generics follow the likelihood.model interface convention:
# each returns a closure that can be evaluated at data + parameters.
#
# Once likelihood.model's dependency chain is fixed (algebraic.mle),
# these should be replaced with re-exports from likelihood.model.
# ===========================================================================

#' Compute log-likelihood function
#'
#' @param model A model object
#' @param ... Additional arguments
#' @return A closure `function(df, par)` returning the log-likelihood
#' @export
loglik <- function(model, ...) UseMethod("loglik")

#' Compute score (gradient of log-likelihood)
#'
#' @param model A model object
#' @param ... Additional arguments
#' @return A closure `function(df, par)` returning the score vector
#' @export
score <- function(model, ...) UseMethod("score")

#' Compute Hessian of log-likelihood
#'
#' @param model A model object
#' @param ... Additional arguments
#' @return A closure `function(df, par)` returning the Hessian matrix
#' @export
hess_loglik <- function(model, ...) UseMethod("hess_loglik")

#' Fit model via maximum likelihood
#'
#' @param model A model object
#' @param ... Additional arguments
#' @return A closure `function(df, par0, ...)` returning MLE results
#' @export
fit <- function(model, ...) UseMethod("fit")

#' Generate random data from model
#'
#' @param model A model object
#' @param ... Additional arguments
#' @return A closure `function(theta, n, ...)` returning a data frame
#' @export
rdata <- function(model, ...) UseMethod("rdata")

#' List model assumptions
#'
#' @param model A model object
#' @param ... Additional arguments
#' @return A list of assumption strings
#' @export
assumptions <- function(model, ...) UseMethod("assumptions")
