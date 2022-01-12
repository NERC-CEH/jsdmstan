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
#'   inverse-link function. The default is FALSE, in which case the
#'   untransformed linear predictor is returned.
#'
#' @param newdata New data, by default NULL and uses original data
#'
#' @param draws Number of draws, by default the number of samples in the
#'   posterior. Will be sampled randomly from the chains if fewer than the
#'   number of samples.
#'
#' @param newdata_type What form is the new data in, at the moment only
#'   supplying covariates is supported.
#'
#' @param list_index Whether to return the output list indexed by the number of
#'   draws (default), species, or site.
#'
#' @export
#' @return A list of linear predictors. If list_index is "draws" (the default)
#'   the list will have length equal to the number of draws with each element of
#'   the list being a site x species matrix. If the list_index is "species" the
#'   list will have length equal to the number of species with each element of
#'   the list being a draws x sites matrix. If the list_index is "sites" the
#'   list will have length equal to the number of sites with each element of the
#'   list being a draws x species matrix.
posterior_linpred.jsdmStanFit <- function(object, transform = FALSE,
                                          newdata = NULL, draws = NULL,
                                          newdata_type = "X",
                                          list_index = "draws", ...){
  if(newdata_type != "X")
    stop("Currently only data on covariates is supported.")
  stopifnot(is.logical(transform))
  if(isTRUE(transform) & object$family == "gaussian")
    warning("No inverse-link transform performed for Gaussian response models.")

  list_index <- match.arg(list_index, c("draws", "species", "sites"))

  dots <- list(...)

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

  model_est <- rstan::extract(object$fit, pars = model_pars)

  n_iter <- dim(model_est[[1]])[1]

  if(!is.null(draws)){
    if(n_iter < draws){
      warning(paste("There are fewer samples than draws specified, defaulting",
                    "to using all iterations"))
      draws <- n_iter
    } else {
      # draw_id <- seq_len(draws)*(floor(n_iter/draws))
      draw_id <- sample.int(n_iter, draws)
      model_est <- lapply(model_est, function(x){
        switch(length(dim(x)),
               `1` = x[draw_id],
               `2` = x[draw_id,],
               `3` =  x[draw_id,,])
      })
    }
  } else{
    draw_id <- seq_len(n_iter)
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
      # for(i in 1:nrow(mu)){
      #   for(j in 1:ncol(mu)){
      #     mu[i,j] <- switch(object$family,
      #                       "gaussian" = mu[i,j],
      #                       "bernoulli" = inv_logit(mu[i,j]),
      #                       "poisson" = exp(mu[i,j]),
      #                       "neg_binomial" = exp(mu[i,j]))
      #   }
      # }
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
#' @param object The model object
#'
#' @param newdata New data, by default NULL and uses original data
#'
#' @param draws Number of draws, by default the number of samples in the
#'   posterior
#'
#' @param newdata_type What form is the new data in, at the moment only
#'   supplying covariates is supported.
#'
#' @param list_index Whether to return the output list indexed by the number of
#'   draws (default), species, or site.
#'
#' @export
#' @return A list of linear predictors. If list_index is "draws" (the default)
#'   the list will have length equal to the number of draws with each element of
#'   the list being a site x species matrix. If the list_index is "species" the
#'   list will have length equal to the number of species with each element of
#'   the list being a draws x sites matrix. If the list_index is "sites" the
#'   list will have length equal to the number of sites with each element of
#'   the list being a draws x species matrix.
posterior_predict.jsdmStanFit <- function(object, newdata = NULL,
                                          newdata_type = "X", draws = NULL,
                                          list_index = "draws"){
  post_linpred <- posterior_linpred(object, newdata = newdata, draws = draws,
                                    newdata_type = newdata_type,
                                    transform = TRUE, list_index = "draws")

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
             "gaussian" = rnorm(1, x, mod_sigma),
             "bernoulli" = rbinom(1, 1, x) ,
             "poisson" = rpois(1, x),
             "neg_binomial" = rnbinom(1, x, mod_kappa))
    }
    )
    # for(i in 1:nrow(x)){
    #   for(j in 1:ncol(x)){
    #     x2[i,j] <- switch(family,
    #                       "gaussian" = rnorm(1,x[i,j],mod_sigma),
    #                       "neg_binomial" = rnbinom(1,x[i,j],mod_kappa),
    #                       "poisson" = rpois(1, x[i,j]),
    #                       "bernoulli" = rbinom(1,1, x[i,j]))
    #   }
    # }
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
