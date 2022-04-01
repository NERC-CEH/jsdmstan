#' Fit jsdm models in Stan
#'
#'@param Y Matrix of species by sites. Rows are assumed to be sites, columns
#'   are assumed to be species
#'
#' @param X The covariates matrix, with rows being site and columns being
#'   covariates
#'
#'@param method Whether to fit a GLLVM or MGLMM model, details in description
#'
#' @param D The number of latent variables within a GLLVM model
#'
#' @param family The response family for the model, required to be one of
#'   \code{"gaussian"}, \code{"bernoulli"}, \code{"poisson"} or \code{"neg_binomial"}
#'
#'@param species_intercept Whether the model should be fit with an intercept by
#'  species, by default \code{TRUE}
#'
#'@param dat_list Alternatively, data can be given to the model as a list containing
#'  Y, X, N, S, K, and site_intercept. See output of mglmm_sim_data for an example of
#'  how this can be formatted.
#'
#'@param site_intercept Whether the model should be fit with a site intercept, by
#'  default \code{FALSE}
#'
#'@param phylo A distance matrix between species (not necessarily phylogenetic). The
#'  default \code{FALSE} does not incorporate phylogenetic information.
#'
#'@param covar The covariance function as a character string, options are Matérn
#'  kernel with \eqn{\nu} 1/2 (\code{"matern_05"}), 3/2 (\code{"matern_15"}), 5/2 (\code{"matern_25"}), or
#'  infinite (\code{"matern_inf"}). Matérn kernel with infinite nu is equivalent to the
#'  squared exponential kernel (\code{"sq_exponential"}), and with \eqn{\nu} = 1/2 the
#'  exponential kernel (\code{"exponential"}).
#'
#'@param delta The constant added to the diagonal of the covariance matrix in the
#'  phylogenetic model to keep matrix semipositive definite, by default 1e-5.
#'
#' @param save_data If the data used to fit the model should be saved in the
#'   model object, by default TRUE.

#' @param ... Arguments passed to \code{rstan::sampling}
#'
#' @return A \code{jsdmStanFit} object
#' @export
stan_jsdm <- function(Y = NULL, X = NULL, species_intercept = TRUE, method,
                      dat_list = NULL, family, site_intercept = FALSE, D = NULL,
                      phylo = FALSE, covar = "matern_05", delta = 1e-5,
                      save_data = TRUE, ...){
  family <- match.arg(family, c("gaussian","bernoulli","poisson","neg_binomial"))

  stopifnot(is.logical(species_intercept),
            is.logical(site_intercept),
            is.logical(save_data))

  if(!isFALSE(phylo)){
    if(method == "gllvm") stop("Phylo only supported with MGLMM")
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

  # validate data
  data_list <- validate_data(Y = Y, X = X, species_intercept = species_intercept,
                             D = D, site_intercept = site_intercept,
                             dat_list = dat_list, phylo = phylo,
                             family = family, method = method)

  # Create stancode
  model_code <- jsdm_stancode(data_list = data_list, family = family,
                              method = method, prior = NULL, phylo = phylo)

  # Compile model
  model_comp <- rstan::stan_model(model_code = model_code)

  # Fit model
  model_fit <- rstan::sampling(model_comp, data = data_list, ...)

  # upper triangle of factor loadings
  if(method == "gllvm"){
    if(data_list$D>1){
      nr <- data_list$D
      upper_trinames <- character()
      while(nr>1){
        upper_trinames <- c(upper_trinames, paste0("[",1:(nr-1),",",nr,"]"))
        nr <- nr-1
      }
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
  dat <- if(isTRUE(save_data)){
    data_list
  } else list()

  model_output <- list(fit = model_fit,
                       jsdm_type = method,
                       family = family,
                       species = species,
                       sites = sites,
                       preds = preds,
                       data_list = dat,
                       n_latent = ifelse(method == "gllvm",
                                         as.integer(round(data_list$D,0)),integer()))

  class(model_output) <- "jsdmStanFit"

  return(model_output)
}

#'@describeIn stan_jsdm Alias for stan_jsdm with method = "mglmm"
stan_mglmm <- function(Y = NULL, X = NULL, species_intercept = TRUE,
                       dat_list = NULL, family, site_intercept = FALSE,
                       phylo = FALSE, covar = "matern_05", delta = 1e-5,
                       save_data = TRUE, ...){
  stan_jsdm(Y = Y, X = X, species_intercept = species_intercept, method = "mglmm",
            dat_list = dat_list, family, site_intercept = site_intercept,
            phylo = phylo, covar = covar, delta = delta,
            save_data = save_data, ...)
}

#'@describeIn stan_jsdm Alias for stan_jsdm with method = "gllvm"
stan_gllvm <- function(Y = NULL, X = NULL, D = NULL, species_intercept = TRUE,
                       dat_list = NULL, family, site_intercept = FALSE,
                       phylo = FALSE, covar = "matern_05", delta = 1e-5,
                       save_data = TRUE, ...){
  stan_jsdm(Y = Y, X = X, D = D, species_intercept = species_intercept,
            method = "gllvm",
            dat_list = dat_list, family, site_intercept = site_intercept,
            phylo = phylo, covar = covar, delta = delta,
            save_data = save_data, ...)
}



validate_data <- function(Y, D, X, species_intercept,
                          dat_list, family, site_intercept,phylo,
                          method){
  method <- match.arg(method, c("gllvm","mglmm"))

  # do things if data not given as list:
  if(is.null(dat_list)){
    if(is.null(X) & !species_intercept)
      stop("Model requires a species intercept if there are no covariates")

    if(method == "gllvm"){
      if(is.null(D))
        stop("Must have at least one latent variable")
      if(D < 1)
        stop("Must have at least one latent variable")
    }

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

    if(method == "mglmm"){
      data_list <- list(Y = Y, S = S, K = K, X = X,
                        site_intercept = site_intercept)
    } else if(method == "gllvm"){
      data_list <- list(Y = Y, D = D, S = S, K = K, X = X,
                        site_intercept = site_intercept)
    }
    if(!isFALSE(phylo)){
      data_list$Dmat <- phylo
      data_list$nu05 <- nu05
      data_list$delta <- delta
    }

  } else {
    if(!all(c("Y","K","S","N","X","site_intercept") %in% names(dat_list)))
      stop("If supplying data as a list must have entries Y, K, S, N, X and site_intercept")

    if(!("D" %in% names(dat_list)) & method == "gllvm")
      stop("If supplying data as a list must have a D entry")

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
    if(!isTRUE(all.equal(unname(data_list$Y),
                         matrix(as.numeric(as.logical(data_list$Y)),
                                nrow=nrow(data_list$Y)) ) ))
      stop("Y matrix is not binary")
  } else if(family %in% c("poisson","neg_binomial")){
    if(!any(apply(data_list$Y, 1:2, is.wholenumber)))
      stop("Y matrix is not composed of integers")
  }

  return(data_list)
}
