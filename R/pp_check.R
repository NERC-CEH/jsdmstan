#' Posterior predictive checks for \code{jsdmStanFit} objects
#'
#' This function acts as an interface to [bayesplot::pp_check()], by
#' default calculating summary statistics for each site (i.e. row in the response
#' matrix) which are then plotted and compared to the data used to fit the model.
#'
#' This function takes a \code{jsdmStanFit} object and tries to extract statistics
#' that can give useful summaries within a posterior prediction (or retrodiction)
#' using the [bayesplot::pp_check()] function. It uses the [jsdm_statsummary()]
#' function to get summary statistics and then supplies them to the specified
#' \code{ppc_} function from the \pkg{bayesplot} package. For a full list of
#' available plotting functions call [bayesplot::available_ppc()].
#'
#' @aliases pp_check
#' @param object An object of class \code{jsdmStanFit}
#' @param plotfun The ppc plot function to use, given as a character string. The
#'   default is to call [ppc_dens_overlay][bayesplot::PPC-distributions]. Can be
#'   specified as either the entire name of function as a character string or without
#'   the ppc_ prefix
#' @inheritParams jsdm_statsummary
#' @param ... Additional arguments passed to [jsdm_statsummary()].
#'
#' @return A ggplot object that can be further customised using the \pkg{ggplot2}
#'   package.
#'
#' @examples
#' \dontrun{
#' # First simulate data and fit the jsdmStan model:
#' mglmm_data <- mglmm_sim_data(
#'   N = 100, S = 10, K = 3,
#'   family = "gaussian"
#' )
#' mglmm_fit <- stan_mglmm(
#'   Y = mglmm_data$Y, X = mglmm_data$X,
#'   family = "gaussian"
#' )
#'
#' # The default is to plot a density overlay:
#' pp_check(mglmm_fit)
#'
#' # Other plot functions can be called, such as a ribbon plot:
#' pp_check(mglmm_fit, plotfun = "ribbon")
#'
#' # Instead of calculating the sum over sites other statistics can be calculated,
#' # e.g. the mean of each species:
#' pp_check(mglmm_fit,
#'   plotfun = "ecdf_overlay", summary_stat = "mean",
#'   calc_over = "species", ndraws = 20
#' )
#' }
#'
#' @importFrom bayesplot pp_check
#' @export pp_check
#' @export
pp_check.jsdmStanFit <- function(object, plotfun = "dens_overlay", species = NULL,
                                 sites = NULL, summary_stat = "sum",
                                 calc_over = "site", ndraws = NULL, ...) {
  # check ppc plot type
  valid_types <- as.character(bayesplot::available_ppc())
  plotfun <- ifelse(grepl("^ppc_", plotfun), plotfun, paste0("ppc_", plotfun))
  if (!plotfun %in% valid_types) {
    stop(paste(
      "plotfun:", plotfun, "is not a valid ppc type. ",
      "Valid types are:\n", paste(valid_types, collapse = ", ")
    ))
  }
  ppc_fun <- get(plotfun, asNamespace("bayesplot"))
  dots <- list(...)

  # get draw IDs
  ndraws_given <- "ndraws" %in% names(match.call())
  nsamps <- dim(object$fit)[1] * dim(object$fit)[2]
  if (ndraws_given) {
    if (is.null(ndraws)) {
      draw_ids <- seq(1, nsamps, 1)
    } else {
      draw_ids <- sample.int(nsamps, ndraws)
    }
  } else {
    aps_plotfuns <- c(
      "ppc_error_scatter_avg", "ppc_error_scatter_avg_vs_x",
      "ppc_intervals", "ppc_intervals_grouped", "ppc_loo_pit",
      "ppc_loo_intervals", "ppc_loo_ribbon", "ppc_ribbon",
      "ppc_ribbon_grouped", "ppc_rootogram", "ppc_scatter_avg",
      "ppc_scatter_avg_grouped", "ppc_stat", "ppc_stat_2d",
      "ppc_stat_freqpoly_grouped", "ppc_stat_grouped",
      "ppc_violin_grouped"
    )
    if (plotfun %in% aps_plotfuns) {
      draw_ids <- seq(1, nsamps, 1)
      message(
        "Using all posterior draws for ppc plot type '",
        plotfun, "' by default."
      )
    } else {
      draw_ids <- sample.int(nsamps, 10)
      message(
        "Using 10 posterior draws for ppc plot type '",
        plotfun, "' by default."
      )
    }
  }

  if (is.character(summary_stat)) {
    stat_fun <- get(summary_stat)
  } else if (inherits(summary_stat, "function")) {
    stat_fun <- summary_stat
  }
  y <- apply(
    object$data_list$Y, switch(calc_over,
      "site" = 1,
      "species" = 2
    ),
    stat_fun
  )

  # Extract all data
  yrep <- jsdm_statsummary(object,
    species = species, sites = sites,
    summary_stat = summary_stat, calc_over = calc_over,
    draw_ids = draw_ids, post_type = "predict", ...
  )

  # prepare plotting arguments
  ppc_args <- list(y = y, yrep = yrep)

  for_pred <- union(
    names(dots) %in% names(formals(jsdm_statsummary)),
    names(dots) %in% names(formals(posterior_linpred.jsdmStanFit))
  )
  ppc_args <- c(ppc_args, dots[!for_pred])

  do.call(ppc_fun, ppc_args)
}



#' Extract summary statistics for a \code{jsdmStanFit} model
#'
#' This function extracts the predicted Y values for within the models and then
#' calculates summary statistics for each simulated community. The default is to sum
#' all the predicted scores for each site.
#'
#' @param object A \code{jsdmStanFit} model object
#' @param species Which species to include in the summary statistic, by default all
#' @param sites Which sites to include in the summary statistic, by default all
#' @param summary_stat The summary statistic to be used, by default \code{sum} but
#'   any function can be used.
#' @param post_type The type of posterior prediction to be used, either
#'   \code{"linpred"} for [posterior_linpred.jsdmStanFit()] or \code{"predict"} for
#'   [posterior_predict.jsdmStanFit()]
#' @param calc_over Whether to calculate the summary statistic by site or species, by
#'   default \code{species}
#' @param simplify Whether to simplify the output into a matrix, by default
#'   \code{TRUE}
#' @param ndraws Number of draws, by default the number of samples in the
#'   posterior. Will be sampled randomly from the chains if fewer than the
#'   number of samples.
#' @param draw_ids The IDs of the draws to be used, as a numeric vector
#' @param ... Arguments passed to the posterior prediction function
#'
#' @return If \code{simplify = TRUE} then a matrix where each row is a draw and each
#'   column is either a site (if \code{calc_over = "site"}) or a species (if
#'   \code{calc_over = "species"}).
#'
#' @seealso pp_check.jsdmStanFit
#'
#' @examples
#' \dontrun{
#' # First simulate data and fit the jsdmStan model:
#' gllvm_data <- gllvm_sim_data(
#'   N = 100, S = 9, D = 2,
#'   family = "bernoulli"
#' )
#' gllvm_fit <- stan_gllvm(dat_list = gllvm_data, family = "bernoulli")
#'
#' # The default is to return a matrix:
#' jsdm_statsummary(gllvm_fit)
#'
#' # The above returns the linear predictor, while we may want to get the posterior
#' # prediction instead:
#' jsdm_statsummary(gllvm_fit, post_type = "predict")
#'
#' # This can be limited to a specific set of species and/or sites:
#' jsdm_statsummary(gllvm_fit, species = 1:5, sites = seq(5, 95, 10))
#' }
#'
#' @export
jsdm_statsummary <- function(object, species = NULL, sites = NULL,
                             summary_stat = "sum", post_type = "linpred",
                             calc_over = "site", simplify = TRUE,
                             draw_ids = NULL, ndraws = NULL,
                             ...) {
  if (!inherits(object, "jsdmStanFit")) {
    stop("jsdm_summary only works for jsdmStanFit objects")
  }
  if (!is.null(species) & !is.character(species)) {
    if (any(!is.wholenumber(species))) {
      stop(paste(
        "Species must be either a character vector of species names or an",
        "integer vector of species positions in the input data columns"
      ))
    }
  }
  if (!is.null(sites) & !is.character(sites)) {
    if (any(!is.wholenumber(sites))) {
      stop(paste(
        "Sites must be either a character vector of site names or an",
        "integer vector of sites positions in the input data columns"
      ))
    }
  }
  calc_over <- match.arg(calc_over, c("site", "species"))
  post_type <- match.arg(post_type, c("linpred", "predict"))

  post_fun <- get(paste0("posterior_", post_type), asNamespace("jsdmstan"))
  post_args <- list(...)
  post_args$object <- object
  post_args$list_index <- "draws"
  post_args$draw_ids <- draw_ids
  post_args$ndraws <- ndraws

  post_res <- do.call(post_fun, post_args)

  if (is.character(summary_stat)) {
    stat_fun <- get(summary_stat)
  } else if (inherits(summary_stat, "function")) {
    stat_fun <- summary_stat
  }

  # Limit to species that have been selected:
  if (!is.null(species)) {
    if (is.character(species)) {
      species_names <- dimnames(post_res[[1]])[[2]]
      if (any(!(species %in% species_names))) {
        stop("Species specified are not found in the model fit object")
      }
      species <- match(species, species_names)
    }
    post_res <- lapply(post_res, "[", , species)
  }
  # Limit to sites that have been selected:
  if (!is.null(sites)) {
    if (is.character(sites)) {
      sites_names <- dimnames(post_res[[1]])[[1]]
      if (any(!(sites %in% sites_names))) {
        stop("Sites specified are not found in the model fit object")
      }
      sites <- match(sites, sites_names)
    }
    post_res <- lapply(post_res, "[", sites, )
  }

  # calculate summary statistic over sites:
  res <- lapply(post_res, function(x) {
    apply(x, switch(calc_over,
      "site" = 1,
      "species" = 2
    ), stat_fun)
  })

  if (simplify) {
    res <- do.call(rbind, res)
  }

  return(res)
}
