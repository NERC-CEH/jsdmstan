#' Update a jsdmStanFit model object with new data or Stan arguments
#'
#' This function allows you to update a jsdmStanFit model with new data or new
#' arguments to Stan. It does not support changes that require recompiling stancode -
#' for that you should use [stan_jsdm()]. Changes to the number of sites, species or
#' covariates do not require recompiling stancode and can therefore be done using this
#' function.
#'
#' @param object The jsdmStanFit model object
#' @param newY New Y data, by default \code{NULL}
#' @param newX New X data, by default \code{NULL}
#' @param newD New number of latent variables, by default \code{NULL}
#' @param newNtrials New number of trials (binomial model only), by default
#' \code{NULL}
#' @param newZi_X New predictor data for the zi parameter in zero-inflated models,
#' by default \code{NULL}. In cases where the model was originally fit with the
#' same X and zi_X data and only newX is supplied to update.jsdmStanFit the zi_X
#' data will also be set to newX.
#' @param newShp_X New predictor data for the family parameter in models where
#' the family parameter is modelled in response to data,
#' by default \code{NULL}. In cases where the model was originally fit with the
#' same X and shp_X data and only newX is supplied to update.jsdmStanFit the shp_X
#' data will also be set to newX.
#' @param save_data Whether to save the data in the jsdmStanFit object, by default
#'  \code{TRUE}
#' @param ... Arguments passed to [rstan::sampling()]
#'
#' @return An object of class \code{jsdmStanFit}
#' @export update
#' @export
#'
#' @examples
#' \dontrun{
#' # MGLMM - specified by using the mglmm aliases and with direct reference to Y and
#' # X matrices:
#' mglmm_data <- mglmm_sim_data(
#'   N = 100, S = 10, K = 3,
#'   family = "gaussian"
#' )
#' mglmm_fit <- stan_mglmm(
#'   Y = mglmm_data$Y, X = mglmm_data$X,
#'   family = "gaussian"
#' )
#' mglmm_fit2 <- update(mglmm_fit, iter = 4000)
#'
#' # You can also run a model by supplying the data as a list:
#' gllvm_data <- jsdm_sim_data(
#'   method = "gllvm", N = 100, S = 6, D = 2,
#'   family = "bernoulli"
#' )
#' gllvm_fit <- stan_jsdm(
#'   dat_list = gllvm_data, method = "gllvm",
#'   family = "bernoulli"
#' )
#' gllvm_fit
#' gllvm_data <- jsdm_sim_data(
#'   method = "gllvm", N = 500, S = 4, D = 2,
#'   family = "bernoulli"
#' )
#' gllvm_fit2 <- update(gllvm_fit, newY = gllvm_data$Y)
#' gllvm_fit2
#' }
update.jsdmStanFit <- function(object, newY = NULL, newX = NULL, newD = NULL,
                               newNtrials = NULL, newZi_X = NULL,
                               newShp_X = NULL,
                               save_data = TRUE, ...) {
  if (length(object$data_list) == 0) {
    stop("Update requires the original data to be saved in the model object")
  }
  # Use new options if specified, otherwise original options
  if (is.null(newX)) {
    X <- object$data_list$X
    if (ncol(X) == 1L) {
      if (colnames(X) == "(Intercept)") {
        X <- NULL
      }
    }
  } else {
    X <- newX
  }
  if (is.null(newY)) {
    Y <- object$data_list$Y
  } else {
    Y <- newY
  }
  family <- object$family$family
  method <- object$jsdm_type
  if(!is.null(newD)){
    D <- newD
  } else{
    D <- object$data_list$D
  }
  if(family == "binomial") {
    if(!is.null(newNtrials)){
      Ntrials <- ntrials_check(newNtrials,  nrow(Y))
    } else{
      Ntrials <- object$data_list$Ntrials
    }
  }
  if ("zi" %in% object$family$params_dataresp){
    if(is.null(newZi_X)) {
      if(isTRUE(all.equal(object$data_list$X, object$family$data_list$zi_X)) & !is.null(newX)){
        zi_X <- newX
      } else{
        zi_X <- object$family$data_list$zi_X
      }
    } else {
      zi_X <- newZi_X
    }
  } else{
    zi_X <- NULL
  }
  if (any(c("sigma","kappa") %in% object$family$params_dataresp)){
    if(is.null(newShp_X)) {
      if(isTRUE(all.equal(object$data_list$X, object$family$data_list$shp_X)) & !is.null(newX)){
        shp_X <- newX
      } else{
        shp_X <- object$family$data_list$shp_X
      }
    } else {
      shp_X <- newShp_X
    }
  } else{
    shp_X <- NULL
  }

  species_intercept <- "(Intercept)" %in% colnames(object$data_list$X)

  site_intercept <- ifelse("ngrp" %in% names(object$data_list), "grouped",
                           ifelse("a" %in% get_parnames(object), "ungrouped",
                                  "none"))
  site_groups <- if(site_intercept == "grouped"){
    object$data_list$grps} else{NULL}

  # validate data
  data_list <- validate_data(
    Y = Y, X = X, species_intercept = species_intercept,
    D = D, site_intercept = site_intercept, site_groups = site_groups,
    dat_list = NULL,
    family = family, method = method, Ntrials = Ntrials,
    zi_X = zi_X, shp_X = shp_X
  )

  # get original stan model
  stanmod <- rstan::get_stanmodel(object$fit)

  model_args <- list(...)
  if (any(c("pars", "include") %in% names(model_args))) {
    warning("Specified pars and include arguments are ignored")
  }
  model_args$object <- stanmod
  model_args$data <- data_list
  if (!"warmup" %in% names(model_args)) {
    if ("iter" %in% names(model_args)) {
      model_args$warmup <- 0.5 * model_args$iter
    } else {
      model_args$warmup <- object$fit@sim$warmup
    }
  }
  if (!"iter" %in% names(model_args)) {
    model_args$iter <- object$fit@sim$iter
  }
  if (!"chains" %in% names(model_args)) {
    model_args$chains <- object$fit@sim$chains
  }
  if (!"thin" %in% names(model_args)) {
    model_args$thin <- object$fit@sim$thin
  }
  model_args$pars <- if (method == "gllvm") c("L", "LV_uncor", "Lambda_uncor") else NA
  model_args$include <- ifelse(method == "gllvm", FALSE, TRUE)

  # Fit model
  model_fit <- do.call(rstan::sampling, model_args)

  model_output <- model_to_jsdmstanfit(
    model_fit = model_fit, method = method,
    data_list = data_list,
    species_intercept = species_intercept,
    family = family, save_data = save_data
  )

  return(model_output)
}
