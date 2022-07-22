#' Efficient approximate leave-one-out cross-validation using the \pkg{loo} package
#'
#' This function uses the \pkg{loo} package to compute PSIS-LOO CV, efficient
#' approximate leave-one-out (LOO) cross-validation for Bayesian models using Pareto
#' smoothed importance sampling (PSIS).
#'
#' @param object The jsdmStanFit model object
#' @param ... Other arguments passed to the \code{\link[loo]{loo}} function
#' @importFrom loo loo
#' @export
#' @aliases loo
#'
#' @return A list with class \code{c("psis_loo","loo")}, as detailed in the
#'   \code{\link[loo]{loo}} documentation
#'
#'
loo.jsdmStanFit <- function(object, ...) {
  rstan::loo(object$fit, ...)
}
