#' Update a jsdmStanFit model object with new data or Stan arguments
#'
#' This function allows you to update a jsdmStanFit model with new data or new
#' arguments to Stan. It does not support changes that require recompiling stancode -
#' for that you should use [stan_jsdm()]. Changes to the number of sites, species or
#' covariates do not require recompiling stancode and can therefore be done using
#' this function.
#'
#' @param object The jsdmStanFit model object
#' @param newY New Y data, by default NULL
#' @param newX New X data, by default NULL
#' @param ... Arguments passed to [rstan::sampling()]
#'
#' @return An object of class jsdmStanFit
#' @export
#'
#' @examples
#' \donttest{
#'  # MGLMM - specified by using the mglmm aliases and with direct reference to Y and
#'  # X matrices:
#'  mglmm_data <- mglmm_sim_data(N = 100, S = 10, K = 3,
#'                               family = "gaussian")
#'  mglmm_fit <- stan_mglmm(Y = mglmm_data$Y, X = mglmm_data$X,
#'                          family = "gaussian")
#'  mglmm_fit2 <- update(mglmm_fit, iter = 4000)
#'
#'  # You can also run a model by supplying the data as a list:
#'  gllvm_data <- jsdm_sim_data(method = "gllvm", N = 100, S = 6, D = 2,
#'                              family = "bernoulli")
#'  gllvm_fit <- stan_jsdm(dat_list = gllvm_data, method = "gllvm",
#'                         family = "bernoulli")
#'  gllvm_fit
#'  gllvm_data <- jsdm_sim_data(method = "gllvm", N = 500, S = 4, D = 2,
#'                              family = "bernoulli")
#'  gllvm_fit2 <- update(gllvm_fit, newY = gllvm_data$Y)
#'  gllvm_fit2
#'
#'}
update.jsdmStanFit <- function(object, newY = NULL, newX = NULL, save_data = TRUE,
                               ...){
  if(length(object$data_list) == 0)
    stop("Update requires the original data to be saved in the model object")
  # Use new options if specified, otherwise original options
  if(is.null(newX)){
    X <- object$data_list$X
    if(ncol(X) == 1L){
      if(colnames(X) == "(Intercept)"){
        X <- NULL
      }
    }
  } else {
    X <- newX
  }
  if(is.null(newY)){
    Y <- object$data_list$Y
  } else{
    Y <- newY
  }
  family <- object$family
  method <- object$jsdm_type
  D <- object$data_list$D
  species_intercept <- "(Intercept)" %in% colnames(object$data_list$X)
  site_intercept <- object$data_list$site_intercept
  phylo <- object$data_list$phylo
  if(!isFALSE(phylo)){
    nu05 <- object$data_list$nu05
    delta <- object$data_list$delta
  } else{
    nu05 <- 0L
    delta <- 1e-5
  }

  # validate data
  data_list <- validate_data(Y = Y, X = X, species_intercept = species_intercept,
                             D = D, site_intercept = site_intercept,
                             dat_list = NULL, phylo = phylo,
                             family = family, method = method, nu05 = nu05,
                             delta = delta)

  # get original stan model
  stanmod <- rstan::get_stanmodel(object$fit)

  model_args <- list(...)
  if(any(c("pars","include") %in% names(model_args)))
    warning("Specified pars and include arguments are ignored")
  model_args$object <- stanmod
  model_args$data <- data_list
  if (!"warmup" %in% names(model_args)){
    if("iter" %in% names(model_args)){
      model_args$warmup <- 0.5*model_args$iter
    } else{
      model_args$warmup <- object$fit@sim$warmup
    }
  }
  if (!"iter" %in% names(model_args)){
    model_args$iter <- object$fit@sim$iter
  }
  if (!"chains" %in% names(model_args)){
    model_args$chains <- object$fit@sim$chains
  }
  if (!"thin" %in% names(model_args)){
    model_args$thin <- object$fit@sim$thin
  }
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
