#' Run gllvm models in Stan
#'
#' @details
#'
#' @param Y Matrix of species by sites. Rows are assumed to be sites, columns are
#'   assumed to be species
#'
#' @param D The number of latent variables
#'
#' @param X The covariates matrix, with rows being site and columns being covariates
#'
#' @param intercept Whether the model should be fit with an intercept
#'
#' @param family The response family for the model, required to be one of "gaussian",
#'   "bernoulli", "poisson" or "neg_binomial"
#'
#' @param ... Arguments passed to rstan::sampling
#'
#' @export

stan_gllvm <- function(Y = NULL, D = NULL, X = NULL, intercept = TRUE,
                       dat_list = NULL, family, ...){
  family <- match.arg(family, c("gaussian","bernoulli","poisson","neg_binomial"))

  # do things if data not given as list:
  if(is.null(dat_list)){
    if(is.null(X) & !intercept)
      stop("Model requires an intercept if there are no covariates")

    if(D<1)
      stop("Must have at least one latent variable")

    S <- ncol(Y)
    N <- nrow(Y)
    if(is.null(X)){
      K <- 1
      X <- matrix(1, nrow = N, ncol = 1)
    } else {
      K <- ncol(X) + 1*intercept
      if(intercept)
        X <- cbind(matrix(1, nrow = N, ncol = 1), X)
    }

    data_list <- list(Y = Y, D = D, S = S, K = K, X = X)
  } else {
    if(!all(c("Y","D","K","S","N","X") %in% names(dat_list)))
      stop("If supplying data as a list must have entries Y, D, K, S, N and X")

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

  return(model_fit)
}
