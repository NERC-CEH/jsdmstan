#' The 'jsdmstan' package.
#'
#' @description A package for fitting joint species distribution models (jSDMs) in
#'   Stan.
#'
#'   This package can be used to both simulate and fit data according to a
#'   multivariate generalised linear mixed model (MGLMM) or a generalised linear
#'   latent variable model (GLLVM). All models are fit using the \pkg{rstan} package.
#'   Summary functions are provided, as are interfaces to the \pkg{bayesplot}
#'   plotting functions
#'
#' @docType package
#' @name jsdmstan-package
#' @aliases jsdmstan
#' @import Rcpp
#'
#' @references Stan Development Team (NA). RStan: the R interface to Stan. R package
#'   version 2.26.1. https://mc-stan.org
#'
NULL

# Suppress R CMD check note
#' @importFrom RcppParallel CxxFlags
#' @importFrom rlang .data
#' @import BH
NULL
