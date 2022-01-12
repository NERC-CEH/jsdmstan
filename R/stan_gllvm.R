#' Run gllvm models in Stan
#'
#' @details
#'
#' @param Y Matrix of species by sites. Rows are assumed to be sites, columns
#'   are assumed to be species
#'
#' @param D The number of latent variables
#'
#' @param X The covariates matrix, with rows being site and columns being
#'   covariates
#'
#' @param intercept Whether the model should be fit with an intercept
#'
#' @param family The response family for the model, required to be one of
#'   "gaussian", "bernoulli", "poisson" or "neg_binomial"
#'
#' @param save_data If the data used to fit the model should be saved in the
#'   model object, by default TRUE.
#'
#' @param ... Arguments passed to rstan::sampling
#'
#' @export

stan_gllvm <- function(Y = NULL, D = NULL, X = NULL, species_intercept = TRUE,
                       dat_list = NULL, family, site_intercept = FALSE,
                       save_data = TRUE, ...){
  family <- match.arg(family, c("gaussian","bernoulli","poisson","neg_binomial"))

  stopifnot(is.logical(species_intercept),
            is.logical(site_intercept),
            is.logical(save_data))

  # do things if data not given as list:
  if(is.null(dat_list)){
    if(is.null(X) & !species_intercept)
      stop("Model requires a species intercept if there are no covariates")

    if(is.null(D))
      stop("Must have at least one latent variable")
    if(D < 1)
      stop("Must have at least one latent variable")

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

    data_list <- list(Y = Y, D = D, S = S, K = K, X = X,
                      site_intercept = site_intercept)
  } else {
    if(!all(c("Y","D","K","S","N","X","site_intercept") %in% names(dat_list)))
      stop("If supplying data as a list must have entries Y, D, K, S, N, X and site_intercept")

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
  ignore_pars <- c("LV_uncor","Lambda_uncor", "L")
  # if(data_list$K == 1){
  #   ignore_pars <- c(ignore_pars, "L_Rho_preds", "L_Rho_preds_vect")
  # }

  # Fit models
  model_fit <- rstan::sampling(switch(family,
                                      gaussian = stanmodels$gllvm_gaussian,
                                      bernoulli = stanmodels$gllvm_bernoulli,
                                      neg_binomial = stanmodels$gllvm_negbin,
                                      poisson = stanmodels$gllvm_poisson),
                               data = data_list, pars = ignore_pars, include = FALSE,
                               ...)

  # upper triangle of factor loadings
  if(data_list$D>1){
    nr <- data_list$D
    upper_trinames <- character()
    while(nr>1){
      upper_trinames <- c(upper_trinames, paste0("[",1:(nr-1),",",nr,"]"))
      nr <- nr-1
    }

  }

  # Turn into jsdmStanFit
  sites <- if(!is.null(rownames(data_list$Y))) rownames(data_list$Y) else
    as.character(1:nrow(data_list$Y))
  species <- if(!is.null(colnames(data_list$Y))) colnames(data_list$Y) else
    as.character(1:ncol(data_list$Y))
  preds <- if(!is.null(colnames(data_list$X))) colnames(data_list$X) else
    as.character(1:ncol(data_list$X))
  if(isTRUE(species_intercept) & !("Intercept" %in% preds)){
    preds <- c("Intercept", preds)
  }

  if(isTRUE(save_data)){
    dat <- data_list
  } else{
    dat <- list()
  }


  model_output <- list(fit = model_fit,
                       jsdm_type = "gllvm",
                       family = family,
                       species = species,
                       sites = sites,
                       preds = preds,
                       data_list = dat,
                       n_latent = as.integer(round(data_list$D,0)),
                       phylo = NULL)

  class(model_output) <- "jsdmStanFit"

  return(model_output)
}

