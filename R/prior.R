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
#' Parameters \code{sigmas_b}, \code{sigma_a}, \code{sigmas_u}, \code{sigma_L},
#' \code{sigma}, \code{kappa}, \code{etasq} and \code{rho} are fixed to be positive
#' only in the stan code and this cannot be changed. Parameters \code{L_Rho_preds}
#' and \code{L_Rho_species} are assumed to be the Cholesky factor of a correlation
#' matrix. All other parameters are real numbers. For all parameters that represent
#' vectors or matrices the prior has to be the same across the entire vector or
#' matrix (note that for the species latent variable loadings in a GLLVM model the
#' prior is set on the non-zero matrix components \code{L} and not on the entire
#' matrix).
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
#' @param sigmas_b The standard deviation of the covariate effects, constrained to be
#'   positive (default standard normal)
#' @param z_preds The covariate effects (default standard normal)
#' @param L_Rho_preds The correlation matrix on the covariate effects (npred by npred
#'   matrix represented as a Cholesky factor of a correlation matrix)
#' @param a The site level intercepts (default standard normal)
#' @param a_bar The mean site level intercept
#' @param sigma_a The standard deviation of the site level intercepts, constrained to
#'   be positive and default prior is half standard normal
#' @param sigmas_u For MGLMM method, the standard deviations of the species
#'   covariances, constrained to be positive (default half standard normal)
#' @param z_species For MGLMM method, the S by N matrix of species covariance by site
#'   (default standard normal)
#' @param L_Rho_species For MGLMM method, the correlation between species represented
#'   as a Cholesky factor correlation matrix, default \code{"lkj_corr_cholesky(1)"}
#' @param LV For GLLVM method, the per site latent variable loadings (default
#'   standard normal)
#' @param L For GLLVM method, the non-zero species latent variable loadings (default
#'   standard normal)
#' @param sigma_L For GLLVM method, the variance of the species loadings, constrained
#'   to be positive (default standard normal)
#' @param sigma For Gaussian response, the standard deviation parameter. Constrained
#'   to be positive, default standard normal
#' @param kappa For negative binomial response, the negative binomial variance
#'   parameter. Constrained to be positive, default standard normal
#' @param etasq For phylogenetic models, the variance parameter of the Matern kernel.
#'   See [cov_matern()] for details.
#' @param rho For phylogenetic models, the length scale parameter of the Matern
#'   kernel. See [cov_matern()] for details.
#'
#' @return An object of class \code{"jsdmprior"} taking the form of a named list
#' @export
#'
#' @examples
#' pr <- jsdm_prior(L_Rho_preds = "normal(0,10)")
#' pr
#'
jsdm_prior <- function(sigmas_b = "normal(0,1)",
                       z_preds = "normal(0,1)",
                       L_Rho_preds = "lkj_corr_cholesky(1)",
                       a = "normal(0,1)",
                       a_bar = "normal(0,1)",
                       sigma_a = "normal(0,1)",
                       sigmas_u = "normal(0,1)",
                       z_species = "normal(0,1)",
                       L_Rho_species = "lkj_corr_cholesky(1)",
                       LV = "normal(0,1)",
                       L = "normal(0,1)",
                       sigma_L = "normal(0,1)",
                       sigma = "normal(0,1)",
                       kappa = "normal(0,1)",
                       etasq = "inv_gamma(10,.1)",
                       rho = "inv_gamma(10,.1)"){
  res <- list(sigmas_b = sigmas_b, z_preds = z_preds, L_Rho_preds = L_Rho_preds,
              a = a, a_bar = a_bar, sigma_a = sigma_a,
              sigmas_u = sigmas_u, z_species = z_species,
              L_Rho_species = L_Rho_species,
              LV = LV, L= L, sigma_L = sigma_L,
              sigma = sigma, kappa = kappa,
              etasq = etasq, rho = rho)
  if(!(all(sapply(res, is.character)))){
    stop("All arguments must be supplied as character vectors")
  }
  class(res) <- c("jsdmprior","list")
  return(res)

}

#' @describeIn jsdm_prior Print method for object of class \code{jsdmprior}
#'
#' @param x Object of class \code{jsdmprior}
#' @param ... Currently unused
#'
#' @return
#' @export
print.jsdmprior <- function(x, ...){
  df <- data.frame(Parameter = names(x),
                   Group = c(rep("species_intercept",3),
                             rep("site_intercept",3),
                             rep("mglmm",3),
                             rep("gllvm",3),
                             "gaussian","neg_binomial","phylo","phylo"),
                   Constraint = c("lower=0",rep("none",4),rep("lower=0",2),
                                  rep("none",4),rep("lower=0",5)),
                   Prior = unlist(unname(x)))
  print(df)
}
