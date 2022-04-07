#' Access the posterior distribution of the linear predictor
#'
#' Extract the posterior draws of the linear predictor, possibly transformed by
#' the inverse-link function.
#'
#' @aliases posterior_linpred
#'
#' @param object The model object
#'
#' @param transform Should the linear predictor be transformed using the
#'   inverse-link function. The default is \code{FALSE}, in which case the
#'   untransformed linear predictor is returned.
#'
#' @param newdata New data, by default \code{NULL} and uses original data
#'
#' @param ndraws Number of draws, by default the number of samples in the
#'   posterior. Will be sampled randomly from the chains if fewer than the
#'   number of samples.
#'
#' @param draw_ids The IDs of the draws to be used, as a numeric vector
#'
#' @param newdata_type What form is the new data in, at the moment only
#'   supplying covariates is supported.
#'
#' @param list_index Whether to return the output list indexed by the number of
#'   draws (default), species, or site.
#' @param ... Currently unused
#'
#' @return A list of linear predictors. If list_index is \code{"draws"} (the default)
#'   the list will have length equal to the number of draws with each element of
#'   the list being a site x species matrix. If the list_index is \code{"species"} the
#'   list will have length equal to the number of species with each element of
#'   the list being a draws x sites matrix. If the list_index is \code{"sites"} the
#'   list will have length equal to the number of sites with each element of the
#'   list being a draws x species matrix.
#'
#' @seealso [posterior_predict.jsdmStanFit()]
#'
#'@importFrom rstantools posterior_linpred
#'@export posterior_linpred
#' @export
posterior_linpred.jsdmStanFit <- function(object, transform = FALSE,
                                          newdata = NULL, ndraws = NULL,
                                          draw_ids = NULL, newdata_type = "X",
                                          list_index = "draws", ...){
  if(newdata_type != "X")
    stop("Currently only data on covariates is supported.")
  stopifnot(is.logical(transform))
  if(isTRUE(transform) & object$family == "gaussian")
    warning("No inverse-link transform performed for Gaussian response models.")
  if(!is.null(ndraws) & !is.null(draw_ids))
    message("Both ndraws and draw_ids have been specified, ignoring ndraws")
  if(!is.null(draw_ids)){
    if(any(!is.wholenumber(draw_ids)))
      stop("draw_ids must be a vector of positive integers")
  }


  list_index <- match.arg(list_index, c("draws", "species", "sites"))

  n_sites <- length(object$sites)
  n_species <- length(object$species)
  n_preds <- length(object$preds)
  method <- object$jsdm_type
  orig_data_used <- is.null(newdata)

  if(is.null(newdata) & length(object$data_list) == 0)
    stop(paste("Original data must be included in model object if no new data",
               "is provided."))
  if(is.null(newdata))
    newdata <- object$data_list$X

  newdata <- validate_newdata(newdata, preds = object$preds,
                              newdata_type = newdata_type)

  model_pars <- "betas"
  if(method == "gllvm")
    model_pars <- c(model_pars, "LV","Lambda","sigma_L") else
      if(method == "mglmm")
        model_pars <- c(model_pars, "u")
  if(isTRUE(grep("a_bar", names(object$fit))))
    model_pars <- c(model_pars, "a","sigma_a", "a_bar")

  model_est <- extract(object, pars = model_pars)

  n_iter <- dim(model_est[[1]])[1]

  if(!is.null(draw_ids)){
    if(max(draw_ids)>n_iter)
      stop(paste("Maximum of draw_ids (",max(draw_ids),
                 ") is greater than number of iterations (",n_iter,")"))

    draw_id <- draw_ids
  } else{
    if(!is.null(ndraws)){
      if(n_iter < ndraws){
        warning(paste("There are fewer samples than ndraws specified, defaulting",
                      "to using all iterations"))
        ndraws <- n_iter
      }
      draw_id <- sample.int(n_iter, ndraws)
      model_est <- lapply(model_est, function(x){
        switch(length(dim(x)),
               `1` = x[draw_id,drop=FALSE],
               `2` = x[draw_id,,drop=FALSE],
               `3` =  x[draw_id,,,drop=FALSE])
      })

    } else{
      draw_id <- seq_len(n_iter)
    }
  }

  model_pred_list <- lapply(seq_along(draw_id), function(d){
    if(method == "gllvm"){
      if(orig_data_used){
        LV_sum <- t((model_est$Lambda[d,,] * model_est$sigma_L[d]) %*% model_est$LV[d,,])
      } else {
        LV_sum <- 0
      }
    } else if(method == "mglmm"){
      if(orig_data_used){
        u_ij <- model_est$u[d,,]
      } else {
        u_ij <- 0
      }
    }
    if(isTRUE(grep("a_bar", names(object$fit)))){
      alpha <- model_est$a_bar[d,] + model_est$a[d,] * model_est$sigma_a[d]
    } else{
      alpha<- 0
    }
    if(is.vector(newdata)){
      newdata <- matrix(newdata, ncol = 1)
    }
    mu <- switch(method,
                 "gllvm" = LV_sum + newdata %*% model_est$betas[d,,] + alpha,
                 "mglmm" = newdata %*% model_est$betas[d,,] + alpha + u_ij)
    if(isTRUE(transform)){
      mu <- apply(mu, 1:2, function(x) {
        switch(object$family,
               "gaussian" = x,
               "bernoulli" = inv_logit(x),
               "poisson" = exp(x),
               "neg_binomial" = exp(x))
      }
      )
    }

    return(mu)

  })

  if(list_index != "draws"){
    model_pred_list <- switch_indices(model_pred_list, list_index)
  }

  return(model_pred_list)

}

#' Draw from the posterior predictive distribution
#'
#' Draw from the posterior predictive distribution of the outcome.
#'
#' @aliases posterior_predict
#'
#' @inheritParams posterior_linpred.jsdmStanFit
#'
#' @return A list of linear predictors. If list_index is \code{"draws"} (the default)
#'   the list will have length equal to the number of draws with each element of
#'   the list being a site x species matrix. If the list_index is \code{"species"} the
#'   list will have length equal to the number of species with each element of
#'   the list being a draws x sites matrix. If the list_index is \code{"sites"} the
#'   list will have length equal to the number of sites with each element of
#'   the list being a draws x species matrix.
#'
#' @seealso [posterior_linpred.jsdmStanFit()]
#'
#'@importFrom rstantools posterior_predict
#'@export posterior_predict
#' @export
posterior_predict.jsdmStanFit <- function(object, newdata = NULL,
                                          newdata_type = "X", ndraws = NULL,
                                          draw_ids = NULL,
                                          list_index = "draws", ...){
  transform <- ifelse(object$family == "gaussian", FALSE, TRUE)
  post_linpred <- posterior_linpred(object, newdata = newdata, ndraws = ndraws,
                                    newdata_type = newdata_type, draw_ids = draw_ids,
                                    transform = transform, list_index = "draws")

  if(object$family == "gaussian")
    mod_sigma <- rstan::extract(object$fit, pars = "sigma", permuted = FALSE)
  if(object$family == "neg_binomial")
    mod_kappa <- rstan::extract(object$fit, pars = "kappa", permuted = FALSE)

  n_sites <- length(object$sites)
  n_species <- length(object$species)

  post_pred <- lapply(post_linpred, function(x, family = object$family){
    x2 <- x
    x2 <- apply(x2, 1:2, function(x) {
      switch(object$family,
             "gaussian" = stats::rnorm(1, x, mod_sigma),
             "bernoulli" = stats::rbinom(1, 1, x) ,
             "poisson" = stats::rpois(1, x),
             "neg_binomial" = rgampois(1, x, mod_kappa))
    }
    )
    x2
  })

  if(list_index != "draws"){
    post_pred <- switch_indices(post_pred, list_index)
  }

  return(post_pred)



}



# internal ~~~~

validate_newdata <- function(newdata, preds, newdata_type){

  preds_nointercept <- preds[preds != "Intercept"]

  if(!all(preds_nointercept %in% colnames(newdata)))
    stop(paste("New data does not have matching column names to model fit.\n",
               "Model has column names:",paste0(preds_nointercept),"\n"))

  newdata <- newdata[,preds_nointercept]
  if("Intercept" %in% preds){
    newdata <- cbind(Intercept = 1, newdata)
    newdata <- newdata[,preds]
  }
  return(newdata)

}

switch_indices <- function(res_list, list_index){
  if(list_index == "species"){
    res_list <- lapply(seq_len(dim(res_list[[1]])[2]), function(sp){
      t(sapply(res_list, "[", , sp))
    })
  } else if(list_index == "sites"){
    res_list <- lapply(seq_len(dim(res_list[[1]])[1]), function(st){
      t(sapply(res_list, "[", st, ))
    })
  } else stop("List index not valid")
}
