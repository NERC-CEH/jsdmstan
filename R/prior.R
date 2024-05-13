#' Create prior object for jsdmstan models and data simulation
#'
#' This function creates all of the potential priors used within a \code{jsdmstan}
#' model and can be used as the input to all [stan_jsdm()] family of functions and
#' the [jsdm_sim_data()] functions.
#'
#' Each prior has to be specified as a character string corresponding to the
#' appropriate stan command. The most common versions of these are supported by the
#' simulated data functions, however there are functions that can be fed to the stan
#' fitting procedure that will not be able to be used as input for [jsdm_sim_data()].
#' Parameters \code{sigmas_preds}, \code{sigma_a}, \code{sigmas_species},
#' \code{sigma_L}, \code{sigma}, and \code{kappa} are fixed to be positive only in
#' the stan code and this cannot be changed. Parameters \code{cor_preds} and
#' \code{cor_species} are assumed to be the Cholesky factor of a correlation matrix.
#' All other parameters are real numbers. For all parameters that represent vectors
#' or matrices the prior has to be the same across the entire vector or matrix (note
#' that for the species latent variable loadings in a GLLVM model the prior is set on
#' the non-zero matrix components \code{L} and not on the entire matrix).
#'
#' Prior distributions supported by [jsdm_sim_data()] are \code{"normal(mean, sd)"},
#' \code{"student_t(df, mu, sigma)"}, \code{"cauchy(location, scale)"},
#' \code{"gamma(shape, scale)"}, \code{"inv_gamma(shape, scale)"} and
#' \code{"lkj_corr_cholesky(eta)"}.
#'
#' @seealso [sim_helpers] for a description of the parameterisations used within the
#'   data simulation functions
#'
#'
#' @param sigmas_preds The standard deviation of the covariate effects, constrained
#'   to be positive (default standard normal)
#' @param z_preds The covariate effects (default standard normal)
#' @param cor_preds The correlation matrix on the covariate effects (npred by npred
#'   correlation matrix) (default
#'   \code{"lkj_corr(1)"})
#' @param betas If covariate effects are unstructured, the prior on the covariate
#'   effects
#' @param a The site level intercepts (default standard normal)
#' @param a_bar The mean site level intercept
#' @param sigma_a The standard deviation of the site level intercepts, constrained to
#'   be positive and default prior is half standard normal
#' @param sigmas_species For MGLMM method, the standard deviations of the species
#'   covariances, constrained to be positive (default half standard normal)
#' @param z_species For MGLMM method, the S by N matrix of species covariance by site
#'   (default standard normal)
#' @param cor_species For MGLMM method, the correlation between species represented
#'   as a nspecies by nspecies correlation matrix (default \code{"lkj_corr(1)"})
#' @param LV For GLLVM method, the per site latent variable loadings (default
#'   standard normal)
#' @param L For GLLVM method, the non-zero species latent variable loadings (default
#'   standard normal)
#' @param sigma_L For GLLVM method, the variance of the species loadings, constrained
#'   to be positive (default standard normal)
#' @param sigma For Gaussian response, the standard deviation parameter. Constrained
#'   to be positive (default standard normal)
#' @param kappa For negative binomial response, the negative binomial variance
#'   parameter. Constrained to be positive (default standard normal)
#' @param zi For zero-inflated poisson or negative binomial, the proportion of
#'   inflated zeros (default beta distribution with both alpha and beta parameters
#'   set to 1).
#'
#' @return An object of class \code{"jsdmprior"} taking the form of a named list
#' @export
#'
#' @examples
#' pr <- jsdm_prior(kappa = "gamma(0.01,0.01)")
#' pr
#'
jsdm_prior <- function(sigmas_preds = "normal(0,1)",
                       z_preds = "normal(0,1)",
                       cor_preds = "lkj_corr(1)",
                       betas = "normal(0,1)",
                       a = "normal(0,1)",
                       a_bar = "normal(0,1)",
                       sigma_a = "normal(0,1)",
                       sigmas_species = "normal(0,1)",
                       z_species = "normal(0,1)",
                       cor_species = "lkj_corr(1)",
                       LV = "normal(0,1)",
                       L = "normal(0,1)",
                       sigma_L = "normal(0,1)",
                       sigma = "normal(0,1)",
                       kappa = "normal(0,1)",
                       zi = "beta(1,1)") {
  res <- list(
    sigmas_preds = sigmas_preds, z_preds = z_preds, cor_preds = cor_preds,
    betas = betas,
    a = a, a_bar = a_bar, sigma_a = sigma_a,
    sigmas_species = sigmas_species, z_species = z_species,
    cor_species = cor_species,
    LV = LV, L = L, sigma_L = sigma_L,
    sigma = sigma, kappa = kappa, zi = zi
  )
  if (!(all(sapply(res, is.character)))) {
    stop("All arguments must be supplied as character vectors")
  }
  class(res) <- c("jsdmprior", "list")
  return(res)
}

#' @describeIn jsdm_prior Print method for object of class \code{jsdmprior}
#'
#' @param x Object of class \code{jsdmprior}
#' @param ... Currently unused
#'
#' @export
print.jsdmprior <- function(x, ...) {
  df <- data.frame(
    Parameter = names(x),
    Group = c(
      rep("covariate_effects", 4),
      rep("site_intercept", 3),
      rep("mglmm", 3),
      rep("gllvm", 3),
      "gaussian", "neg_binomial","zero_inflation"
    ),
    Constraint = c(
      "lower=0", rep("none", 5), rep("lower=0", 2),
      rep("none", 4), rep("lower=0", 3),"lower=0,upper=1"
    ),
    Prior = unlist(unname(x))
  )
  print(df)
}
