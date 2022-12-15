#' jsdmStanFit class
#'
#' This is the jsdmStanFit class, which stan_gllvm and stan_mglmm both create.
#'
#' @name jsdmStanFit
#'
#' @section Elements for \code{jsdmStanFit} objects:
#' \describe{
#'  \item{\code{jsdm_type}}{
#'  A length one character vector describing type of jSDM
#'  }
#'  \item{\code{family}}{
#'  A character vector describing response family
#'  }
#'  \item{\code{species}}{
#'  A character vector of the species names
#'  }
#'  \item{\code{sites}}{
#'   character vector of the site IDs
#'  }
#'  \item{\code{preds}}{
#'  A character vector of the measured predictors included
#'  }
#'  \item{\code{data_list}}{
#'  A list containing the original data used to fit the model
#'   (empty when save_data is set to \code{FALSE})
#'  }
#'  \item{\code{n_latent}}{
#'  A length one integer vector representing number of latent
#'   variables (in gllvm type fits) or \code{NA} in all other cases
#'  }
#'  }
#'
#'
jsdmStanFit_empty <- function() {
  res <- list(
    jsdm_type = "None",
    family = character(),
    species = character(),
    sites = character(),
    preds = character(),
    data_list = list(),
    n_latent = integer()
  )
  class(res) <- "jsdmStanFit"
  res
}

# Check if model object is okay

# jsdmstanfit methods

#' Print the default summary for the model
#'
#' This prints out a summary for the models which includes the type of model fit, the
#' number of species, sites and predictors as well as a summary of any parameters
#' with Rhat > 1.01 or effective sample size to total number of samples ratio < 0.05
#'
#' @param x The \code{jsdmStanFit} model object
#' @param ... Other arguments passed to [summary.jsdmStanFit]
#'
#' @export
print.jsdmStanFit <- function(x, ...) {
  rhat_prob <- rhat(x)
  prob_rhat <- names(stats::na.omit(rhat_prob[rhat_prob > 1.01]))
  neff_prob <- neff_ratio(x)
  prob_neff <- names(stats::na.omit(neff_prob[neff_prob < 0.05]))
  prob_pars <- union(prob_rhat, prob_neff)

  if (length(prob_pars) > 0) {
    s <- summary(x, pars = prob_pars, ...)
  }

  cat("Model type: ", x$jsdm_type, if (x$jsdm_type == "gllvm") {
    paste0(" with ", x$n_latent, " latent variables")
  }, "\n",
  "  Number of species: ", length(x$species), "\n",
  "  Number of sites: ", length(x$sites), "\n",
  "  Number of predictors: ", length(x$preds), "\n",
  "\n",
  "Model run on ", length(x$fit@stan_args), " chains with ",
  x$fit@stan_args[[1]]$iter, " iterations per chain (",
  x$fit@stan_args[[1]]$warmup, " warmup).\n\n",
  sep = ""
  )
  if (length(prob_pars) > 0) {
    cat("Parameters with Rhat > 1.01, or Neff/N < 0.05:\n")

    print(s)
  } else {
    cat("No parameters with Rhat > 1.01 or Neff/N < 0.05")
  }
}

#' Summarise the model fit and data structure and give summaries for the parameters
#'
#' This returns a matrix of parameter summaries including a summary of the parameter
#' estimates, R-hat, bulk ESS and tail ESS. This can be limited to parameters with
#' Rhat > 1.01 or ESS < 500 by setting \code{prob_pars_only = TRUE}.
#'
#' @param object The model object
#'
#' @param prob_quantiles The quantiles to summarise the parameter estimates with, by
#'   default the 15% and 85% quantiles
#'
#' @param digit_summary The number of digits to round the results to
#'
#' @param prob_pars_only Whether to limit output to parameters with Rhat > 1.01 or
#'   effective sample size < 500, by default \code{FALSE}
#'
#' @param na_filter Whether to remove parameters with NAs in Rhat - this includes the
#'   parameters fixed to zero or one, such as the upper triangle of the cholesky
#'   factor of the correlation matrix. By default \code{TRUE}
#'
#' @param pars Parameters to compute the summary of, by default \code{NULL} i.e. all
#'   parameters included
#'
#' @param log_lik Whether the log_lik parameters should be included, default
#'   \code{FALSE}
#'
#' @param ... Arguments passed to [extract()]
#'
#' @examples
#' \dontrun{
#'
#' gllvm_data <- jsdm_sim_data(
#'   method = "gllvm", N = 100, S = 6, D = 2,
#'   family = "bernoulli"
#' )
#' gllvm_fit <- stan_jsdm(
#'   dat_list = gllvm_data, method = "gllvm",
#'   family = "bernoulli"
#' )
#' gllvm_summ <- summary(gllvm_fit)
#' head(gllvm_summ, 20)
#'
#' summary(gllvm_fit, prob_quantiles = c(0.05, 0.5, 0.95))
#' }
#'
#' @export
summary.jsdmStanFit <- function(object,
                                prob_quantiles = c(0.15, 0.85),
                                digit_summary = 3,
                                prob_pars_only = FALSE,
                                pars = NULL,
                                na_filter = TRUE, log_lik = FALSE, ...) {
  if (is.null(pars)) {
    full_pars <- get_parnames(object, log_lik = log_lik)
    samps <- extract(object, pars = full_pars, return_array = TRUE, ...)
  } else {
    samps <- extract(object, pars = pars, return_array = TRUE, ...)
  }
  rhat <- apply(samps, 3, rstan::Rhat)
  bulk_ess <- round(apply(samps, 3, rstan::ess_bulk), 0)
  tail_ess <- round(apply(samps, 3, rstan::ess_tail), 0)

  mean <- apply(samps, 3, mean)
  sd <- apply(samps, 3, sd)
  quants <- t(apply(samps, 3, stats::quantile, probs = prob_quantiles))

  s <- cbind(
    mean = mean, sd = sd, quants, Rhat = rhat,
    Bulk.ESS = bulk_ess, Tail.ESS = tail_ess
  )

  if (isTRUE(prob_pars_only)) {
    prob_pars <- rhat > 1.01 | bulk_ess < 500 | tail_ess < 500
    s <- s[prob_pars, , drop = FALSE]
  }
  if (isTRUE(na_filter)) s <- s[!is.na(s[, "Rhat"]), , drop = FALSE]

  s <- round(s, digit_summary)

  return(summary = s)
}

#' Extract samples from jsdmStanFit object
#'
#' This function extracts named parameters from a jsdmStanFit object, with optional
#' regular expression matching.
#'
#' @param object The jsdmStanFit model object
#' @param pars A character vector of parameter names
#' @param permuted Whether the draws should be randomly permuted, by default
#'   \code{FALSE}
#' @param inc_warmup Whether the warmup period should be included, by default
#'   \code{FALSE}
#' @param include Whether the parameters specified by \code{pars} should be included
#'   or excluded from the result, by default \code{TRUE} for inclusion
#' @param regexp Whether regular expression matching should be used to match the
#'   contents of \code{pars} to the parameter names, by default \code{FALSE}
#' @param return_array Whether to return the output as a 3 dimensional array
#'   (\code{TRUE}) or a named list (\code{FALSE}, the default)
#' @param ... Arguments passed to [get_parnames()]
#'
#' @return If \code{return_array = FALSE} returns a named list with each parameter
#'   group being an element of the list. Each list element is an array with the first
#'   dimension being the iteration (all chains are appended) and the other dimensions
#'   coming from the parameter dimensions. If \code{return_array = TRUE} then a 3
#'   dimensional array is returned with the first dimension being the iterations, the
#'   second the chains and the third the parameters.
#'
#' @export
extract.jsdmStanFit <- function(object, pars = NULL, permuted = FALSE,
                                inc_warmup = FALSE, include = TRUE, regexp = FALSE,
                                return_array = FALSE, ...) {
  if (isTRUE(return_array) & isTRUE(permuted)) {
    warning("If return_array is TRUE then permuted is set to FALSE")
    permuted <- FALSE
  }
  if (is.null(pars)) {
    pars <- get_parnames(object, ...)
  } else if (regexp) {
    full_pars <- get_parnames(object)
    pars <- grep(paste(pars, collapse = "|"), full_pars, value = TRUE)
  }
  pexp <- rstan::extract(object$fit, pars, permuted, inc_warmup, include)
  if (isFALSE(permuted) & isFALSE(return_array)) {
    pexp2 <- matrix(pexp, prod(dim(pexp)[1:2]), dim(pexp)[3])
    colnames(pexp2) <- dimnames(pexp)[[3]]
    rownames(pexp2) <- paste0(
      "Chain", rep(seq(1, dim(pexp)[2], 1), each = dim(pexp)[1]),
      "_Iter", rep(seq(1, dim(pexp)[1], 1), dim(pexp)[2])
    )

    pexp <- pars_indexes(pexp2)
  }

  pexp
}
#' @export
#' @describeIn extract.jsdmStanFit Generic method
extract <- function(object, ...) {
  UseMethod("extract")
}

pars_indexes <- function(x) {
  pars_unique <- unique(sapply(strsplit(colnames(x), "\\["), "[", 1))
  x_list <- lapply(pars_unique, function(y) x[, grepl(y, colnames(x))])
  x_list <- lapply(x_list, function(y) {
    if (class(y)[1] == "numeric") {
      val <- array(y, dim = length(y))
    } else {
      ynames <- colnames(y)
      ninit <- nrow(y)
      ndim <- sum(grepl(",", ynames[1])) + 1
      if (ndim == 1) {
        val <- y
      } else {
        dims <- apply(sapply(
          strsplit(
            unlist(regmatches(
              ynames,
              gregexpr("\\[\\K\\d+,\\d+(?=\\])",
                ynames,
                perl = TRUE
              )
            )),
            ","
          ), as.numeric
        ), 1, max)
        val <- array(y, dim = c(ninit, dims))
      }
    }
    val
  })
  names(x_list) <- pars_unique
  x_list
}

#' Extract quantities useful for model summaries
#'
#' These are methdos for extracting various useful summaries from models, including
#' the model parameter names, NUTS parameters, the log posterior, r-hat and n-eff
#' ratio.
#'
#' @name jsdmstan-extractors
#'
#' @return{
#' \code{get_parnames()} returns a character vector of model parameter names.
#'
#' \code{nuts_params()} returns a molten data frame (see [reshape2::melt()]). The
#' data frame should have columns "Parameter" (factor), "Iteration" (integer),
#' "Chain" (integer), and "Value" (numeric).
#'
#' \code{log_posterior()} returns a molten data frame (see [reshape2::melt()]). The
#' data frame should have columns "Chain" (integer), "Iteration" (integer),
#' and "Value" (numeric).
#'
#' \code{rhat()}, \code{neff_ratio()} both return named numeric vectors.
#' }
NULL

#' @describeIn jsdmstan-extractors Get the model parameter names
#' @param log_lik Whether the log_lik parameters should be included in the parameter
#'   list
#'
#' @export
get_parnames <- function(object, log_lik = FALSE) {
  if (!inherits(object, "jsdmStanFit")) {
    stop("This only works for jsdmStanFit objects")
  }
  parnames <- names(object$fit)
  if (isFALSE(log_lik)) {
    parnames <- parnames[!grepl("log_lik", parnames)]
  }
  parnames <- parnames[!grepl("_uncor", parnames)]
  return(parnames)
}

#' @describeIn jsdmstan-extractors Get the NUTS parameters
#' @aliases nuts_params
#' @param object The \code{jsdmStanFit} model object
#' @param ... Arguments passed on to the \pkg{bayesplot} equivalent for stanFit
#'   objects
#' @importFrom bayesplot nuts_params
#' @export nuts_params
#' @export
nuts_params.jsdmStanFit <- function(object, ...) {
  bayesplot::nuts_params(object$fit, ...)
}

#' @describeIn jsdmstan-extractors Get the log posterior
#' @aliases log_posterior
#' @importFrom bayesplot log_posterior
#' @export log_posterior
#' @export
log_posterior.jsdmStanFit <- function(object, ...) {
  bayesplot::log_posterior(object$fit, ...)
}

#' @describeIn jsdmstan-extractors Get the R-hat
#' @aliases rhat
#' @importFrom bayesplot rhat
#' @export rhat
#' @export
rhat.jsdmStanFit <- function(object, ...) {
  bayesplot::rhat(object$fit, ...)
}

#' @describeIn jsdmstan-extractors Get the n_eff ratio
#' @aliases neff_ratio
#' @importFrom bayesplot neff_ratio
#' @export neff_ratio
#' @export
neff_ratio.jsdmStanFit <- function(object, ...) {
  bayesplot::neff_ratio(object$fit, ...)
}
