#' jsdmStanFit fitted jsdm_stan model
#'
#' @slot jsdm_type A length one character vector describing type of jSDM
#'
#' @slot species A character vector of the species names
#'
#' @slot sites A character vector of the site IDs
#'
#' @slot preds A character vector of the measured predictors included
#'
#' @slot n_latent A length one integer vector representing number of latent
#'   variables (in gllvm type fits) or NA in all other cases
setClass("jsdmStanFit",
         contains = "stanfit",
         slots = c(
           jsdm_type = "character",
           species = "character",
           sites = "character",
           preds = "character",
           n_latent = "integer"
         ),
         prototype = list(
           jsdm_type = NA_character_,
           species = NA_character_,
           sites = NA_character_,
           preds = NA_character_,
           n_latent = NA_integer_
         )
         )


# jsdmstanfit methods

#' @describeIn jsdmStanFit Print the default summary for the model
setMethod("show", "jsdmStanFit", function(object) {
  s <- summary(object, prob_pars_only = TRUE)

  cat("Model type: ", object@jsdm_type, if(object@jsdm_type == "gllvm"){
    paste0(" with ", object@n_latent, " latent variables")
  }, "\n",
  "  Number of species: ", length(object@species), "\n",
  "  Number of sites: ", length(object@sites), "\n",
  "  Number of predictors: ", length(object@preds), "\n",
  "\n",
  "Model run on ", length(object@stan_args), " chains with ",
  object@stan_args[[1]]$iter, " iterations per chain (",
  object@stan_args[[1]]$warmup," warmup).\n\n",
  sep = ""
  )
  if(nrow(s)>0){
    cat("Parameters with Rhat > 1.01, or ESS < 500:\n")

    print(s)
  }
})

#' @describeIn jsdmStanFit Summarise the model fit and data structure and give
#'   information on the parameter estimates
#'
#' @param object The model object
#'
#' @param prob_quantiles The quantiles to summarise the parameter estimates with, by
#'   default the 15% and 85% quantiles
#'
#' @param digit_summary The number of digits to round the results to
#'
#' @param prob_pars_only Whether to limit output to parameters with Rhat > 1.01 or
#'   effective sample size < 500, by default FALSE
#'
#' @param na_filter Whether to remove parameters with NAs in Rhat - this includes the
#'   parameters fixed to zero or one, such as the upper triangle of the cholesky
#'   factor of the correlation matrix. By default TRUE
setMethod("summary", "jsdmStanFit", function(object,
                                             prob_quantiles = c(0.15,0.85),
                                             digit_summary = 3,
                                             prob_pars_only = FALSE,
                                             na_filter = TRUE) {
  samps <- rstan::extract(object, permuted = FALSE)
  rhat <- apply(samps, 3, rstan::Rhat)
  bulk_ess <- round(apply(samps, 3, rstan::ess_bulk),0)
  tail_ess <- round(apply(samps, 3, rstan::ess_tail),0)

  mean <- apply(samps, 3, mean)
  sd <- apply(samps, 3, sd)
  quants <- t(apply(samps, 3, quantile, probs = prob_quantiles))

  s <- cbind(mean = mean, sd = sd, quants, Rhat = rhat,
             Bulk.ESS = bulk_ess, Tail.ESS = tail_ess)

  if(isTRUE(prob_pars_only)){
    prob_pars <- rhat>1.01 | bulk_ess < 500 | tail_ess < 500
    s <- s[prob_pars,,drop=FALSE]
  }
  if(isTRUE(na_filter))  s <- s[!is.na(s[,"Rhat"]),,drop=FALSE]

  s <- round(s, digit_summary)

  return(summary = s)
})
