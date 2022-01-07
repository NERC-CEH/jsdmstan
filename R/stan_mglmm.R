#'Run mglmm models in Stan
#'
#'@details
#'
#'@param Y Matrix of species by sites. Rows are assumed to be sites, columns are
#'  assumed to be species
#'
#'@param X The covariates matrix, with rows being site and columns being covariates
#'
#'@param species_intercept Whether the model should be fit with an intercept by
#'  species, by default TRUE
#'
#'@param dat_list Alternatively, data can be given to the model as a list containing
#'  Y, X, N, S, K, and site_intercept. See output of mglmm_sim_data for an example of
#'  how this can be formatted.
#'
#'@param family The response family for the model, required to be one of "gaussian",
#'  "bernoulli", "poisson" or "neg_binomial"
#'
#'@param site_intercept Whether the model should be fit with a site intercept, by
#'  default FALSE
#'
#'@param phylo A distance matrix between species (not necessarily phylogenetic). The
#'  default FALSE does not incorporate phylogenetic information.
#'
#'@param covar The covariance function as a character string, options are Matérn
#'  kernel with \eqn{\nu} 1/2 ("matern_05"), 3/2 ("matern_15"), 5/2 ("matern_25"), or
#'  infinite ("matern_inf"). Matérn kernel with infinite nu is equivalent to the
#'  squared exponential kernel ("sq_exponential"), and with \eqn{\nu} = 1/2 the
#'  exponential kernel ("exponential").
#'
#'@param delta The constant added to the diagonal of the covariance matrix in the
#'  phylogenetic model to keep matrix semipositive definite, by default 1e-5.
#'
#'@param ... Arguments passed to rstan::sampling
#'
#'@export

stan_mglmm <- function(Y = NULL, X = NULL, species_intercept = TRUE,
                       dat_list = NULL, family, site_intercept = FALSE,
                       phylo = FALSE, covar = "matern_05", delta = 1e-5, ...){
  family <- match.arg(family, c("gaussian","bernoulli","poisson","neg_binomial"))

  if(!isFALSE(phylo)){
    if(!is.matrix(phylo)) stop("Phylo must be either FALSE or a matrix")
    if(!all(dim(phylo) == ncol(Y)))
      stop("Phylo must be a square matrix with dimensions equal to ncol(Y)")
    if(!isSymmetric(phylo) | any(eigen(phylo)$values<0))
      stop("Phylo must be symmetric and semipositive definite")

    covar <- match.arg(covar, c("matern_05","exponential","matern_15","matern_25",
                                "sq_exponential","matern_inf"))
    nu05 <- switch(covar, "matern_05" = 0L,
                   "exponential" = 0L,
                   "matern_15" = 1L,
                   "matern_25" = 2L,
                   "sq_exponential" = 3L,
                   "matern_inf" = 3L)
  }

  # do things if data not given as list:
  if(is.null(dat_list)){
    if(is.null(X) & !species_intercept)
      stop("Model requires a species intercept if there are no covariates")

    S <- ncol(Y)
    N <- nrow(Y)
    if(is.null(X)){
      K <- 1
      X <- matrix(1, nrow = N, ncol = 1)
    } else {
      K <- ncol(X) + 1*species_intercept
      if(species_intercept)
        X <- cbind(matrix(1, nrow = N, ncol = 1), X)
    }
    site_intercept <- as.integer(site_intercept)


    data_list <- list(Y = Y, S = S, K = K, X = X,
                      site_intercept = site_intercept)

    if(!isFALSE(phylo)){
      data_list$Dmat <- phylo
      data_list$nu05 <- nu05
      data_list$delta <- delta
    }

  } else {
    if(!all(c("Y","K","S","N","X","site_intercept") %in% names(dat_list)))
      stop("If supplying data as a list must have entries Y, K, S, N, X and site_intercept")

    if(!isFALSE(phylo)){
      if(!all(c("Dmat","nu05","delta") %in% names(dat_list)))
        stop("Phylo models require Dmat, nu05 and delta in dat_list")
    }

    data_list <- dat_list


  }

  if(!is.matrix(data_list$Y)){
    if(is.data.frame(data_list$Y)){
      if(all(sapply(data_list$Y, is.numeric))){
        data_list$Y <- as.matrix(data_list$Y)
      } else{
        stop("Non-numeric column(s) detected in Y")
      }
    } else{
      stop("Y must be either a numeric matrix or a dataframe consisting of only numeric columns")
    }
  }

  # Check if Y is appropriate given family
  if(identical(family, "bernoulli")){
    if(!isTRUE(all.equal(data_list$Y,
                         matrix(as.numeric(as.logical(data_list$Y)),
                                nrow=nrow(data_list$Y)) ) ))
      stop("Y matrix is not binary")
  } else if(family %in% c("poisson","neg_binomial")){
    if(any(apply(data_list$Y, 1:2, function(x) x%%1 != 0)))
      stop("Y matrix is not composed of integers")
  }


  # Need to specify pars to ignore
  # ignore_pars <- c("L_Rho_preds","L_Rho_species")

  # Fit models
  if(isFALSE(phylo)){
    model_fit <- rstan::sampling(switch(family,
                                        gaussian = stanmodels$mglmm_gaussian,
                                        bernoulli = stanmodels$mglmm_bernoulli,
                                        neg_binomial = stanmodels$mglmm_negbin,
                                        poisson = stanmodels$mglmm_poisson),
                                 data = data_list,# pars = ignore_pars, include = FALSE,
                                 ...)
  } else {
    model_fit <- rstan::sampling(switch(family,
                                        gaussian = stanmodels$mglmm_gaussian_phylo,
                                        bernoulli = stanmodels$mglmm_bernoulli_phylo,
                                        neg_binomial = stanmodels$mglmm_negbin_phylo,
                                        poisson = stanmodels$mglmm_poisson_phylo),
                                 data = data_list,# pars = ignore_pars, include = FALSE,
                                 ...)
  }
  # Turn into jsdmStanFit
  sites <- if(!is.null(rownames(data_list$Y))) rownames(data_list$Y) else
    as.character(1:nrow(data_list$Y))
  species <- if(!is.null(colnames(data_list$Y))) colnames(data_list$Y) else
    as.character(1:ncol(data_list$Y))
  preds <- if(!is.null(colnames(data_list$X))) colnames(data_list$X) else
    as.character(1:ncol(data_list$X))
  if(isTRUE(species_intercept)){
    preds <- c("Intercept", preds)
  }

  model_output <- list(fit = model_fit,
                       jsdm_type = "mglmm",
                       species = species,
                       sites = sites,
                       preds = preds)

  class(model_output) <- "jsdmStanFit"

  return(model_output)
}
