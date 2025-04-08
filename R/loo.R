#' Efficient approximate leave-one-out cross-validation using the \pkg{loo} package
#'
#' This function uses the \pkg{loo} package to compute PSIS-LOO CV, efficient
#' approximate leave-one-out (LOO) cross-validation for Bayesian models using Pareto
#' smoothed importance sampling (PSIS). This requires that the model was fit using
#' \code{log_lik = TRUE}.
#'
#' @param x The jsdmStanFit model object
#' @param calc_reff Whether to calculate the relative efficiencies for loo, by
#'   default \code{TRUE}. If set to \code{FALSE} then relative efficiency is
#'   assumed to be 1.
#' @param cores The number of cores the loo functions use, by default uses the
#'   mc.cores option (or 1, if unspecified).
#' @param ... Other arguments passed to the \code{\link[loo]{loo}} function
#' @importFrom loo loo
#' @export
#' @aliases loo
#'
#' @return A list with class \code{c("psis_loo","loo")}, as detailed in the
#'   \code{\link[loo]{loo}} documentation
#'
#'
loo.jsdmStanFit <- function(x, calc_reff = TRUE,
                              cores = getOption("mc.cores", 1), ...) {
  # setup
  jsdm_type <- x$jsdm_type
  # get family
  family <- x$family$family

  # get all parameters to extract
  extr_pars <- .par_extract_names(family = x$family,
                                  jsdm_type = x$jsdm_type)
  parameter_draws_1 <- extract(object = x,
                               pars = extr_pars)

  ll <- llfun(data_i = x$data_list, draws = parameter_draws_1,
              jsdm_type = jsdm_type,
              family = family, log = log)

  # calculate r_eff if necessary
  if(isTRUE(calc_reff)){
    nchains <- length(x$fit@stan_args)
    chain_id <- rep(seq_len(nchains),
                    each = x$fit@stan_args[[1]][["iter"]] -
                      x$fit@stan_args[[1]][["warmup"]])
    r_eff <- loo::relative_eff(exp(ll),
                               cores = cores,
                               chain_id = chain_id)

  } else{
    r_eff <- 1
  }
  loo::loo(ll, r_eff = r_eff, cores = cores, ...)
}

llfun <- function(data_i, draws, jsdm_type, family, log = TRUE){
  x <- data_i$X
  y <- data_i$Y

  if(family == "binomial"){
    Ntrials <- rep(data_i$Ntrials, nrow(y))
  }

  ll <- sapply(seq_len(dim(draws[[1]])[1]), function(di){
    if(family == "gaussian"){
      if("shp_betas" %in% names(draws)){
        sigma <- c(exp(data_i$shp_X %*% draws$shp_betas[di,,]))
      } else{
        sigma <- rep(draws$sigma[di,],each = ncol(y))
      }
    } else if(family %in% c("neg_binomial","zi_neg_binomial")){
      if("shp_betas" %in% names(draws)){
        kappa <- c(exp(data_i$shp_X %*% draws$shp_betas[di,,]))
      } else{
        kappa <- rep(draws$kappa[di,],each = ncol(y))
      }
    }
    if(family %in% c("zi_neg_binomial","zi_poisson")){
      if("zi_betas" %in% names(draws)){
        zi <- data_i$zi_X %*% draws$zi_betas[di,,]
      } else{
        zi <- rep(draws$zi[di,],each=ncol(y))
      }
    }
    linpred_i <- switch(
      jsdm_type,
      "gllvm" = (x %*% draws$betas[di,,]) + t((draws$Lambda[di,,] * draws$sigma_L[di]) %*% draws$LV[di,,]),
      "mglmm" = (x %*% draws$betas[di,,]) + draws$u[di,,]
    )
    ll_vals <- switch(
      family,
      "gaussian" = stats::dnorm(x = c(y), mean = c(linpred_i), sd = sigma,
                                log = log),
      "poisson" = stats::dpois(x = c(y), exp(c(linpred_i)), log = log),
      "neg_binomial" = stats::dnbinom(x = c(y), mu = exp(c(linpred_i)),
                                      size = kappa, log = log),
      "bernoulli" = stats::dbinom(x = c(y), size = 1, prob = inv_logit(c(linpred_i)),
                                  log = log),
      "binomial" = stats::dbinom(x = c(y), size = c(Ntrials), prob = inv_logit(c(linpred_i)),
                                 log = log),
      "zi_poisson" = dzipois(x = c(y), mu = c(linpred_i), zi = c(zi), log = log),
      "zi_neg_binomial" = dzinb(x = c(y), mu = c(linpred_i), kappa = c(kappa),
                                zi = c(zi), log = log)
    )
  })
  return(t(ll))

}
dzipois <- function(x, mu, zi, log = TRUE){
  if(isFALSE(all(is.wholenumber(x)))){
    warning("Not all values are integers")
  }
  ll_1 <- ifelse(x == 0,
                 log(exp(stats::dbinom(1,1,inv_logit(zi), log = TRUE)) +
                       exp(stats::dbinom(0,1,inv_logit(zi), log = TRUE) +
                             stats::dpois(x, exp(mu), log = TRUE))),
                 stats::dbinom(0,1,inv_logit(zi), log = TRUE) +
                   stats::dpois(x, exp(mu),log = TRUE)
  )
  if(isFALSE(log))
    ll_1 <- exp(ll_1)
  return(ll_1)
}

dzinb <- function(x, mu, kappa, zi, log = TRUE){
  if(isFALSE(all(is.wholenumber(x)))){
    warning("Not all values are integers")
  }
  ll_1 <- ifelse(x == 0,
                 log(exp(stats::dbinom(1,1,inv_logit(zi), log = TRUE)) +
                       exp(stats::dbinom(0,1,inv_logit(zi), log = TRUE) +
                             stats::dnbinom(x, mu = exp(mu),
                                            size = kappa, log = TRUE))),
                 stats::dbinom(0,1,inv_logit(zi), log = TRUE) +
                   stats::dnbinom(x, mu = exp(mu),
                                  size = kappa, log = TRUE))
  if(isFALSE(log))
    ll_1 <- exp(ll_1)
  return(ll_1)
}


# for calling moment match directly - currently doesn't work
# post_draws_jsdmstan <- function(x, ...){
#   as.matrix(x$fit, ...)
# }
# log_lik_i_jsdmstan <- function(x, i, ...){
#   extr_pars <- .par_extract_names(family = x$family,
#                                   jsdm_type = x$jsdm_type)
#   par_draws <- extract(x, extr_pars)
#   ll_i <- llfun(data_i = x$data_list,
#                 draws = par_draws,
#                 jsdm_type = x$jsdm_type,
#                 family = x$family$family,
#                 log = TRUE)[,i]
#   nchains <- length(x$fit@stan_args)
#   ll_im <- matrix(ll_i, ncol = nchains)
#   ll_im
#
# }
# unconstrain_pars_jsdmstan <- function(x, pars, ...){
#   skeleton <- .create_skeleton(x$fit@sim$pars_oi, x$fit@par_dims[x$fit@sim$pars_oi])
#   if(x$jsdm_type == "gllvm"){
#     skeleton$LV_uncor <- skeleton$LV
#     skeleton$Lambda_uncor <- skeleton$Lambda
#     skeleton$L <- c(skeleton$Lambda)[-1]
#     pars_uncor <- pars[,grepl(c("LV|Lambda"),colnames(pars))]
#     colnames(pars_uncor) <- gsub("LV","LV_uncor",colnames(pars_uncor))
#     colnames(pars_uncor) <- gsub("Lambda","Lambda_uncor",colnames(pars_uncor))
#     pars_L <- pars[,grepl("Lambda",colnames(pars))]
#     pars_L <- pars_L[,-1]
#     colnames(pars_L) <- paste0("L[",1:ncol(pars_L),"]")
#     pars_uncor <- cbind(pars_uncor, pars_L)
#     pars <- cbind(pars, pars_uncor)
#   }
#   upars <- apply(pars, 1, FUN = function(theta) {
#     rstan::unconstrain_pars(x$fit, .rstan_relist(theta, skeleton))
#   })
#   t(upars)
# }
# log_prob_upars_jsdmstan <- function(x, upars, ...){
#   apply(upars, 1, rstan::log_prob, object = x$fit,
#         adjust_transform = TRUE, gradient = FALSE)
# }
# log_lik_i_upars_jsdmstan <- function(x, upars, i, parameter_name = "log_lik",
#                                      ...){
#   S <- nrow(upars)
#   out <- numeric(S)
#   for (s in seq_len(S)) {
#     out[s] <- rstan::constrain_pars(x$fit, upars = upars[s, ])[[parameter_name]][i]
#   }
#   out
# }
# .create_skeleton <- function(pars, dims) {
#   out <- lapply(seq_along(pars), function(i) {
#     len_dims <- length(dims[[i]])
#     if (len_dims < 1) return(0)
#     return(array(0, dim = dims[[i]]))
#   })
#   names(out) <- pars
#   out
# }
# .rstan_relist <- function(x, skeleton) {
#   out <- utils::relist(x, skeleton)
#   for (i in seq_along(skeleton)) {
#     dim(out[[i]]) <- dim(skeleton[[i]])
#   }
#   out
# }

# helpers
.par_extract_names <- function(family,jsdm_type){
  extr_pars <- c("betas")
  if(family$family == "gaussian"){
    if(length(family$params_dataresp)>0){
      extr_pars <- c(extr_pars, "shp_betas")
    } else{
      extr_pars <- c(extr_pars, "sigma")
    }
  } else if(family$family %in% c("neg_binomial","zi_neg_binomial")){
    if("kappa" %in% family$params_dataresp){
      extr_pars <- c(extr_pars, "shp_betas")
    } else{
      extr_pars <- c(extr_pars, "kappa")
    }
    if(family$family == "zi_neg_binomial"){
      if("zi" %in% family$params_dataresp){
        extr_pars <- c(extr_pars, "zi_betas")
      } else{
        extr_pars <- c(extr_pars, "zi")
      }
    }
  } else if(family$family == "zi_poisson"){
    if("zi" %in% family$params_dataresp){
      extr_pars <- c(extr_pars, "zi_betas")
    } else{
      extr_pars <- c(extr_pars, "zi")
    }
  }
  if(jsdm_type == "gllvm"){
    extr_pars <- c(extr_pars, "Lambda","sigma_L","LV")
  } else if(jsdm_type == "mglmm"){
    extr_pars <- c(extr_pars, "u")
  }
  return(extr_pars)
}
