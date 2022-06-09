#'Fit jsdm models in Stan
#'
#'This function fits joint Species Distribution models in Stan, using either a
#'generalised linear latent variable model approach (\code{method = "gllvm"}), or a
#'multivariate generalised linear mixed model approach (\code{method = "mglmm"}).
#'
#'@param formula The formula of covariates that the species means are modelled from
#'
#'@param data Dataframe or list of covariates.
#'
#'@param Y Matrix of species by sites. Rows are assumed to be sites, columns are
#'  assumed to be species
#'
#'@param X The covariates matrix, with rows being site and columns being covariates.
#'  Ignored in favour of data when formula approach is used to specify model.
#'
#'@param method Whether to fit a GLLVM or MGLMM model, details in description
#'
#'@param D The number of latent variables within a GLLVM model
#'
#'@param family The response family for the model, required to be one of
#'  \code{"gaussian"}, \code{"bernoulli"}, \code{"poisson"} or \code{"neg_binomial"}
#'
#'@param species_intercept Whether the model should be fit with an intercept by
#'  species, by default \code{TRUE}
#'
#'@param dat_list Alternatively, data can be given to the model as a list containing
#'  Y, X, N, S, K, and site_intercept. See output of jsdm_sim_data for an example of
#'  how this can be formatted.
#'
#'@param site_intercept Whether the model should be fit with a site intercept, by
#'  default \code{FALSE}
#'
#'@param phylo A distance matrix between species (not necessarily phylogenetic). The
#'  default \code{FALSE} does not incorporate phylogenetic information.
#'
#'@param covar The covariance function as a character string, options are Matérn
#'  kernel with \eqn{\nu} 1/2 (\code{"matern_05"}), 3/2 (\code{"matern_15"}), 5/2
#'  (\code{"matern_25"}), or infinite (\code{"matern_inf"}). Matérn kernel with
#'  infinite nu is equivalent to the squared exponential kernel
#'  (\code{"sq_exponential"}), and with \eqn{\nu} = 1/2 the exponential kernel
#'  (\code{"exponential"}).
#'
#'@param delta The constant added to the diagonal of the covariance matrix in the
#'  phylogenetic model to keep matrix semipositive definite, by default 1e-5.
#'
#'@param save_data If the data used to fit the model should be saved in the model
#'  object, by default TRUE.
#'
#'@param iter A positive integer specifying the number of iterations for each chain,
#'  default 4000.
#'
#'@param ... Arguments passed to \code{\link[rstan]{sampling}}
#'
#'@return A \code{jsdmStanFit} object, comprising a list including the StanFit
#'  object, the data used to fit the model plus a few other bits of information. See
#'  [jsdmStanFit] for details.
#'
#' @examples
#'
#' \donttest{
#'  # MGLMM - specified by using the mglmm aliases and with direct reference to Y and
#'  # X matrices:
#'  mglmm_data <- mglmm_sim_data(N = 100, S = 10, K = 3,
#'                               family = "gaussian")
#'  mglmm_fit <- stan_mglmm(Y = mglmm_data$Y, X = mglmm_data$X,
#'                          family = "gaussian")
#'  mglmm_fit
#'
#'  # You can also run a model by supplying the data as a list:
#'  gllvm_data <- jsdm_sim_data(method = "gllvm", N = 100, S = 6, D = 2,
#'                              family = "bernoulli")
#'  gllvm_fit <- stan_jsdm(dat_list = gllvm_data, method = "gllvm",
#'                         family = "bernoulli")
#'  gllvm_fit
#'
#'
#'}
#'@export
stan_jsdm <- function(X, ...) UseMethod("stan_jsdm")

#' @describeIn stan_jsdm this is the default way of doing things
#' @export
stan_jsdm.default <- function(X = NULL, Y = NULL, species_intercept = TRUE, method,
                      dat_list = NULL, family, site_intercept = FALSE, D = NULL,
                      phylo = FALSE, covar = "matern_05", delta = 1e-5,
                      save_data = TRUE, iter = 4000, ...){
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
  } else{
    nu05 <- NULL
  }

  # validate data
  data_list <- validate_data(Y = Y, X = X, species_intercept = species_intercept,
                             D = D, site_intercept = site_intercept,
                             dat_list = dat_list, phylo = phylo,
                             family = family, method = method, nu05 = nu05,
                             delta = delta)

  # Create stancode
  model_code <- jsdm_stancode(family = family,
                              method = method, prior = NULL, phylo = phylo)

  # Compile model
  model_comp <- rstan::stan_model(model_code = model_code)

  model_args <- list(...)
  if(any(c("pars","include") %in% names(model_args)))
    warning("Specified pars and include arguments are ignored")
  model_args$object <- model_comp
  model_args$data <- data_list
  model_args$iter <- iter
  model_args$pars <- if(method == "gllvm") c("L","LV_uncor","Lambda_uncor") else NA
  model_args$include <- ifelse(method == "gllvm",FALSE,TRUE)

  # Fit model
  model_fit <- do.call(rstan::sampling, model_args)

  model_output <- model_to_jsdmstanfit(model_fit = model_fit, method = method,
                                       data_list = data_list,
                                       species_intercept = species_intercept,
                                       family = family, save_data = save_data)

  return(model_output)
}

#' @describeIn stan_jsdm Formula interface
#' @export
stan_jsdm.formula <- function(formula, data = list(), ...){
  mf <- stats::model.frame(formula = formula, data = data)
  X <- stats::model.matrix(attr(mf, "terms"), data = mf)
  if(!is.null(stats::model.response(mf)))
    warning("Response variable in formula is ignored")

  est <- stan_jsdm.default(X, species_intercept = FALSE, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}

#'Alias for \code{stan_jsdm} with \code{method = "mglmm"}
#'@inheritParams stan_jsdm
#'@export
stan_mglmm <- function(X, ...) UseMethod("stan_mglmm")

#' @describeIn stan_mglmm Default
#' @export
stan_mglmm.default <- function(X = NULL, Y = NULL, species_intercept = TRUE,
                               dat_list = NULL, family, site_intercept = FALSE,
                               phylo = FALSE, covar = "matern_05", delta = 1e-5,
                               save_data = TRUE, iter = 4000, ...){
  stan_jsdm.default(X = X, Y = Y, species_intercept = species_intercept, method = "mglmm",
                    dat_list = dat_list, family, site_intercept = site_intercept,
                    phylo = phylo, covar = covar, delta = delta,
                    save_data = save_data, iter = iter, ...)
}
#' @describeIn stan_mglmm Formula interface
#' @export
stan_mglmm.formula <- function(formula, data = list(), ...){
  stan_jsdm.formula(formula, data = data, method = "mglmm", ...)
}


#'Alias for \code{stan_jsdm} with \code{method = "gllvm"}
#'@inheritParams stan_jsdm
#'@export
stan_gllvm <- function(X, ...) UseMethod("stan_gllvm")

#' @describeIn stan_gllvm Default
#' @export
stan_gllvm.default <- function(X = NULL, Y = NULL, D = NULL, species_intercept = TRUE,
                       dat_list = NULL, family, site_intercept = FALSE,
                       phylo = FALSE, covar = "matern_05", delta = 1e-5,
                       save_data = TRUE, iter = 4000, ...){
  stan_jsdm.default(X = X, Y = Y, D = D, species_intercept = species_intercept,
                    method = "gllvm",
                    dat_list = dat_list, family, site_intercept = site_intercept,
                    phylo = phylo, covar = covar, delta = delta,
                    save_data = save_data, iter = iter, ...)
}
#' @describeIn stan_gllvm Formula interface
#' @export
stan_gllvm.formula <- function(formula, data = list(), ...){
  stan_jsdm.formula(formula, data = data, method = "gllvm", ...)
}


validate_data <- function(Y, D, X, species_intercept,
                          dat_list, family, site_intercept,phylo,
                          method, nu05, delta){
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
      colnames(X) <- "(Intercept)"
    } else {
      K <- ncol(X) + 1*species_intercept
      if(isTRUE(species_intercept)){
        X <- cbind(matrix(1, nrow = N, ncol = 1), X)
        colnames(X)[1] <- "(Intercept)"
      }
    }
    site_intercept <- as.integer(site_intercept)

    if(method == "mglmm"){
      data_list <- list(Y = Y, S = S, K = K, X = X, N = N,
                        site_intercept = site_intercept)
    } else if(method == "gllvm"){
      data_list <- list(Y = Y, D = D, S = S, K = K, X = X, N = N,
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

    if(isTRUE(species_intercept)){
      dat_list$X <- cbind(matrix(1, nrow = dat_list$N, ncol = 1), dat_list$X)
      colnames(dat_list$X)[1] <- "(Intercept)"
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

model_to_jsdmstanfit <- function(model_fit, method, data_list, species_intercept,
                                 family, save_data){
  # upper triangle of factor loadings
  # if(method == "gllvm"){
  #   if(data_list$D>1){
  #     nr <- data_list$D
  #     upper_trinames <- character()
  #     while(nr>1){
  #       upper_trinames <- c(upper_trinames, paste0("[",1:(nr-1),",",nr,"]"))
  #       nr <- nr-1
  #     }
  #   }
  # }

  # Turn into jsdmStanFit
  sites <- if(!is.null(rownames(data_list$Y))) rownames(data_list$Y) else
    as.character(1:nrow(data_list$Y))
  species <- if(!is.null(colnames(data_list$Y))) colnames(data_list$Y) else
    as.character(1:ncol(data_list$Y))
  preds <- if(!is.null(colnames(data_list$X))) colnames(data_list$X) else
    as.character(1:ncol(data_list$X))
  if(isTRUE(species_intercept) & !("(Intercept)" %in% preds)){
    preds <- c("(Intercept)", preds)
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
