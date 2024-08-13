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
#'   list being a draws x species matrix. Note that in the zero-inflated case this is
#'   only the linear predictor of the non-zero-inflated part of the model.
#'
#' @seealso [posterior_predict.jsdmStanFit()]
#'
#' @importFrom rstantools posterior_linpred
#' @export posterior_linpred
#' @export
posterior_linpred.jsdmStanFit <- function(object, transform = FALSE,
                                          newdata = NULL, ndraws = NULL,
                                          draw_ids = NULL, newdata_type = "X",
                                          list_index = "draws", ...) {
  if (newdata_type != "X") {
    stop("Currently only data on covariates is supported.")
  }
  stopifnot(is.logical(transform))
  if (isTRUE(transform) & object$family == "gaussian") {
    warning("No inverse-link transform performed for Gaussian response models.")
  }
  if (!is.null(ndraws) & !is.null(draw_ids)) {
    message("Both ndraws and draw_ids have been specified, ignoring ndraws")
  }
  if (!is.null(draw_ids)) {
    if (any(!is.wholenumber(draw_ids))) {
      stop("draw_ids must be a vector of positive integers")
    }
  }


  list_index <- match.arg(list_index, c("draws", "species", "sites"))

  n_sites <- length(object$sites)
  n_species <- length(object$species)
  n_preds <- length(object$preds)
  method <- object$jsdm_type
  orig_data_used <- is.null(newdata)

  if (is.null(newdata) & length(object$data_list) == 0) {
    stop(paste(
      "Original data must be included in model object if no new data",
      "is provided."
    ))
  }
  if (is.null(newdata)) {
    newdata <- object$data_list$X
  }

  newdata <- validate_newdata(newdata,
    preds = object$preds,
    newdata_type = newdata_type
  )

  model_pars <- "betas"
  if (method == "gllvm") {
    model_pars <- c(model_pars, "LV", "Lambda", "sigma_L")
  } else
  if (method == "mglmm") {
    model_pars <- c(model_pars, "u")
  }
  if (isTRUE(grep("a_bar", names(object$fit)))) {
    model_pars <- c(model_pars, "a", "sigma_a", "a_bar")
  }

  model_est <- extract(object, pars = model_pars)
  n_iter <- dim(model_est[[1]])[1]

  draw_id <- draw_id_check(draw_ids = draw_ids, n_iter = n_iter, ndraws = ndraws)

  model_est <- lapply(model_est, function(x) {
    switch(length(dim(x)),
           `1` = x[draw_id, drop = FALSE],
           `2` = x[draw_id, , drop = FALSE],
           `3` =  x[draw_id, , , drop = FALSE]
    )
  })

  model_pred_list <- lapply(seq_along(draw_id), function(d) {
    if (method == "gllvm") {
      if (orig_data_used) {
        if(object$n_latent>1){
          LV_sum <- t((model_est$Lambda[d, , ] * model_est$sigma_L[d]) %*% model_est$LV[d, , ])
        } else{
          LV_sum <- t((matrix(model_est$Lambda[d, , ], ncol = 1) * model_est$sigma_L[d])
                      %*% matrix(model_est$LV[d, , ], nrow = 1))
        }
      } else {
        LV_sum <- 0
      }
    } else if (method == "mglmm") {
      if (orig_data_used) {
        u_ij <- model_est$u[d, , ]
      } else {
        u_ij <- 0
      }
    }
    if (isTRUE(grep("a_bar", names(object$fit)))) {
      alpha <- model_est$a_bar[d, ] + model_est$a[d, ] * model_est$sigma_a[d]
    } else {
      alpha <- 0
    }
    if (is.vector(newdata)) {
      newdata <- matrix(newdata, ncol = 1)
    }
    mu <- switch(method,
      "gllvm" = LV_sum + newdata %*% model_est$betas[d, , ] + alpha,
      "mglmm" = newdata %*% model_est$betas[d, , ] + alpha + u_ij
    )
    if (isTRUE(transform)) {
      mu <- apply(mu, 1:2, function(x) {
        switch(object$family,
          "gaussian" = x,
          "bernoulli" = inv_logit(x),
          "poisson" = exp(x),
          "neg_binomial" = exp(x),
          "binomial" = inv_logit(x),
          "zi_poisson" = exp(x),
          "zi_neg_binomial" = exp(x)
        )
      })
    }

    return(mu)
  })

  if (list_index != "draws") {
    model_pred_list <- switch_indices(model_pred_list, list_index)
  }

  return(model_pred_list)
}

#'Draw from the posterior predictive distribution
#'
#'Draw from the posterior predictive distribution of the outcome.
#'
#'@aliases posterior_predict
#'
#'@inheritParams posterior_linpred.jsdmStanFit
#'
#'@param Ntrials For the binomial distribution the number of trials, given as either
#'  a single integer which is assumed to be constant across sites or as a site-length
#'  vector of integers.
#'
#'@param include_zi For the zero-inflated poisson distribution, whether to include
#'  the zero-inflation in the prediction. Defaults to \code{TRUE}.
#'
#'@return A list of linear predictors. If list_index is \code{"draws"} (the default)
#'  the list will have length equal to the number of draws with each element of the
#'  list being a site x species matrix. If the list_index is \code{"species"} the
#'  list will have length equal to the number of species with each element of the
#'  list being a draws x sites matrix. If the list_index is \code{"sites"} the list
#'  will have length equal to the number of sites with each element of the list being
#'  a draws x species matrix.
#'
#'@seealso [posterior_linpred.jsdmStanFit()]
#'
#'@importFrom rstantools posterior_predict
#'@export posterior_predict
#'@export
posterior_predict.jsdmStanFit <- function(object, newdata = NULL,
                                          newdata_type = "X", ndraws = NULL,
                                          draw_ids = NULL,
                                          list_index = "draws",
                                          Ntrials = NULL,
                                          include_zi = TRUE, ...) {
  transform <- ifelse(object$family == "gaussian", FALSE, TRUE)
  if (!is.null(ndraws) & !is.null(draw_ids)) {
    message("Both ndraws and draw_ids have been specified, ignoring ndraws")
  }
  if (!is.null(draw_ids)) {
    if (any(!is.wholenumber(draw_ids))) {
      stop("draw_ids must be a vector of positive integers")
    }
  }
  n_iter <- length(object$fit@stan_args)*(object$fit@stan_args[[1]]$iter -
                                            object$fit@stan_args[[1]]$warmup)
  draw_id <- draw_id_check(draw_ids = draw_ids, n_iter = n_iter, ndraws = ndraws)


  post_linpred <- posterior_linpred(object,
    newdata = newdata,
    newdata_type = newdata_type, draw_ids = draw_id,
    transform = transform, list_index = "draws"
  )

  if (object$family == "gaussian") {
    mod_sigma <- extract(object, pars = "sigma")[[1]][draw_id,]
  } else  if (object$family == "neg_binomial") {
    mod_kappa <- extract(object, pars = "kappa")[[1]][draw_id,]
  } else  if(object$family == "binomial"){
    if(is.null(newdata)) {
      Ntrials <- object$data_list$Ntrials
    } else {
      Ntrials <- ntrials_check(Ntrials, nrow(newdata))
    }
  } else if(object$family == "zi_poisson"){
    mod_zi <- extract(object, pars = "zi")[[1]][draw_id,]
  } else if(object$family == "zi_neg_binomial"){
    mod_zi <- extract(object, pars = "zi")[[1]][draw_id,]
    mod_kappa <- extract(object, pars = "kappa")[[1]][draw_id,]
  }

  n_sites <- length(object$sites)
  n_species <- length(object$species)

  post_pred <- lapply(seq_along(post_linpred),
                      function(x, family = object$family) {
    x2 <- post_linpred[[x]]
    if(family == "binomial"){
      for(i in 1:nrow(x2)){
        for(j in 1:ncol(x2)){
          x2[i,j] <- stats::rbinom(1, Ntrials[i], x2[i,j])
        }
      }
    } else {
      for(i in seq_len(nrow(x2))){
        for(j in seq_len(ncol(x2))){
          x2[i,j] <- switch(
            object$family,
            "gaussian" = stats::rnorm(1, x2[i,j], mod_sigma[x,j]),
            "bernoulli" = stats::rbinom(1, 1, x2[i,j]),
            "poisson" = stats::rpois(1, x2[i,j]),
            "neg_binomial" = rgampois(1, x2[i,j], mod_kappa[x,j]),
            "zi_poisson" = if(isTRUE(include_zi)){
              (1-stats::rbinom(1, 1, mod_zi[x,j]))*stats::rpois(1, x2[i,j])
            } else {
              stats::rpois(1, x2[i,j])
            },
            "zi_neg_binomial" = if(isTRUE(include_zi)){
              (1-stats::rbinom(1, 1, mod_zi[x,j]))*rgampois(1, x2[i,j], mod_kappa[x,j])
            } else {
              rgampois(1, x2[i,j], mod_kappa[x,j])
            }
          )
      }
      }
    }
    x2
  })

  if (list_index != "draws") {
    post_pred <- switch_indices(post_pred, list_index)
  }

  return(post_pred)
}



# internal ~~~~

validate_newdata <- function(newdata, preds, newdata_type) {
  preds_nointercept <- preds[preds != "(Intercept)"]

  if (!all(preds_nointercept %in% colnames(newdata))) {
    stop(paste(
      "New data does not have matching column names to model fit.\n",
      "Model has column names:", paste0(preds_nointercept), "\n"
    ))
  }

  newdata <- newdata[, preds_nointercept]
  if ("(Intercept)" %in% preds) {
    newdata <- cbind(`(Intercept)` = 1, newdata)
    newdata <- newdata[, preds]
  }
  return(newdata)
}

switch_indices <- function(res_list, list_index) {
  if (list_index == "species") {
    res_list <- lapply(seq_len(dim(res_list[[1]])[2]), function(sp) {
      t(sapply(res_list, "[", , sp))
    })
  } else if (list_index == "sites") {
    res_list <- lapply(seq_len(dim(res_list[[1]])[1]), function(st) {
      t(sapply(res_list, "[", st, ))
    })
  } else {
    stop("List index not valid")
  }
}

draw_id_check <- function(draw_ids, n_iter, ndraws){
  if (!is.null(draw_ids)) {
    if (max(draw_ids) > n_iter) {
      stop(paste(
        "Maximum of draw_ids (", max(draw_ids),
        ") is greater than number of iterations (", n_iter, ")"
      ))
    }

    draw_id <- draw_ids
  } else {
    if (!is.null(ndraws)) {
      if (n_iter < ndraws) {
        warning(paste(
          "There are fewer samples than ndraws specified, defaulting",
          "to using all iterations"
        ))
        ndraws <- n_iter
      }
      draw_id <- sample.int(n_iter, ndraws)
    } else {
      draw_id <- seq_len(n_iter)
    }
  }
  return(draw_id)
}
