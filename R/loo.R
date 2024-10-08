#' Efficient approximate leave-one-out cross-validation using the \pkg{loo} package
#'
#' This function uses the \pkg{loo} package to compute PSIS-LOO CV, efficient
#' approximate leave-one-out (LOO) cross-validation for Bayesian models using Pareto
#' smoothed importance sampling (PSIS). This requires that the model was fit using
#' \code{log_lik = TRUE}.
#'
#' @param x The jsdmStanFit model object
#' @param ... Other arguments passed to the \code{\link[loo]{loo}} function
#' @importFrom loo loo
#' @export
#' @aliases loo
#'
#' @return A list with class \code{c("psis_loo","loo")}, as detailed in the
#'   \code{\link[loo]{loo}} documentation
#'
#'
loo.jsdmStanFit <- function(x, ...) {
  rstan::loo(x$fit, ...)
}
