#' The 'jsdmstan' package.
#'
#' @description A package for fitting joint species distribution models (jSDMs)
#'   in Stan.
#'
#'   This package can be used to both simulate and fit data according to a
#'   multivariate generalised linear mixed model (MGLMM) or a generalised linear
#'   latent variable model (GLLVM).
#'
#' @docType package
#' @name jsdmstan-package
#' @aliases jsdmstan
#' @useDynLib jsdmstan, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#' @importFrom RcppParallel RcppParallelLibs CxxFlags
#'
#' @references Stan Development Team (NA). RStan: the R interface to Stan. R
#'   package version 2.26.1. https://mc-stan.org
#'
NULL
