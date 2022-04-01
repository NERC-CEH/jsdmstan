#' Posterior predictive checks for \code{jsdmStanFit} objects
#'
#' @param object An object of class \code{jsdmStanFit}
#' @param plotfun The ppc plot function to use, given as a character string. The
#'   default is to call \code{\link[bayesplot:PPC-distributions]{ppc_dens_overlay}}. Can be specified as
#' @inheritParams jsdm_statsummary
#' @param ... Additional arguments passed to \link[jsdmstan]{jsdm_statsummary}.
#'
#' @return
#' @export
#'
pp_check.jsdmStanFit <- function(object, plotfun = "dens_overlay", species = NULL,
                                 sites = NULL, summary_stat = "sum",
                                 calc_over = "site", ndraws = NULL, ...){
  # check ppc plot type
  valid_types <- as.character(bayesplot::available_ppc())
  plotfun <- ifelse(grepl("^ppc_",plotfun),plotfun, paste0("ppc_",plotfun))
  if(!plotfun %in% valid_types)
    stop(paste("plotfun:",plotfun, "is not a valid ppc type. ",
               "Valid types are:\n", paste(valid_types, collapse = ", ")))
  ppc_fun <- get(plotfun, asNamespace("bayesplot"))
  dots <- list(...)

  # get draw IDs
  ndraws_given <- "ndraws" %in% names(match.call())
  nsamps <- dim(object$fit)[1] * dim(object$fit)[2]
  if(ndraws_given){
    if(is.null(ndraws)){
      draw_ids <- seq(1,nsamps,1)
    } else{
      draw_ids <- sample.int(nsamps, ndraws_given)
    }
  } else{
    aps_plotfuns <- c(
      "error_scatter_avg", "error_scatter_avg_vs_x",
      "intervals", "intervals_grouped", "loo_pit",
      "loo_intervals", "loo_ribbon", "ribbon",
      "ribbon_grouped", "rootogram", "scatter_avg",
      "scatter_avg_grouped", "stat", "stat_2d",
      "stat_freqpoly_grouped", "stat_grouped",
      "violin_grouped"
    )
    if (plotfun %in% aps_plotfuns) {
      draw_ids <- seq(1,nsamps,1)
      message("Using all posterior draws for ppc plot type '",
              plotfun, "' by default.")
    } else {
      draw_ids <- sample.int(nsamps, 10)
      message("Using 10 posterior draws for ppc plot type '",
              plotfun, "' by default.")
    }
  }

  if(is.character(summary_stat)){
    stat_fun <- get(summary_stat)
  } else if(class(summary_stat) == "function"){
    stat_fun <- summary_stat
  }
  y <- apply(object$data_list$Y, switch(calc_over,"site" = 1,
                                        "species" = 2),
             stat_fun)

  # Extract all data
  yrep <- jsdm_statsummary(object, species = species, sites = sites,
                           summary_stat = summary_stat, calc_over = calc_over,
                           draw_ids = draw_ids, post_type = "predict", ...)

  # prepare plotting arguments
  ppc_args <- list(y = y, yrep = yrep)

  for_pred <- union(names(dots) %in% names(formals(jsdm_statsummary)),
                    names(dots) %in% names(formals(posterior_linpred.jsdmStanFit)))
  ppc_args <- c(ppc_args, dots[!for_pred])

  do.call(ppc_fun, ppc_args)


}



#' Extract summary statistics for a \code{jsdmStanFit} model
#'
#' @param object A \code{jsdmStanFit} model object
#' @param species Which species to include in the summary statistic, by default all
#' @param sites Which sites to include in the summary statistic, by default all
#' @param summary_stat The summary statistic to be used, by default \code{sum} but any
#'   function can be used.
#' @param post_type The type of posterior prediction to be used, either
#'   \code{"linpred"} for \link[jsdmstan]{posterior_linpred} or \code{"predict"} for
#'   \link[jsdmstan]{posterior_predict}
#' @param calc_over Whether to calculate the summary statistic by site or species, by
#'   default \code{species}
#' @param simplify Whether to simplify the output into a matrix, by default
#'   \code{TRUE}
#' @param ... Arguments passed to the posterior prediction function
#'
#' @return
#' @export
jsdm_statsummary <- function(object, species = NULL, sites = NULL,
                             summary_stat = "sum", post_type = "linpred",
                             calc_over = "site", simplify = TRUE,
                             ...){
  if(class(object) != "jsdmStanFit")
    stop("jsdm_summary only works for jsdmStanFit objects")
  if(!is.null(species) & !is.character(species)){
    if(any(!is.wholenumber(species)))
      stop(paste("Species must be either a character vector of species names or an",
                 "integer vector of species positions in the input data columns"))
  }
  if(!is.null(sites) & !is.character(sites)){
    if(any(!is.wholenumber(sites)))
      stop(paste("Sites must be either a character vector of site names or an",
                 "integer vector of sites positions in the input data columns"))
  }
  calc_over <- match.arg(calc_over, c("site","species"))
  post_type <- match.arg(post_type, c("linpred","predict"))

  post_fun <- get(paste0("posterior_",post_type), asNamespace("jsdmstan"))
  post_args <- list(...)
  post_args$object <- object
  post_args$list_index <- "draws"

  post_res <- do.call(post_fun, post_args)

  if(is.character(summary_stat)){
    stat_fun <- get(summary_stat)
  } else if(class(stat_summary) == "function"){
    stat_fun <- stat_summary
  }

  # Limit to species that have been selected:
  if(!is.null(species)){
    if(is.character(species)){
      species_names <- dimnames(post_res[[1]])[[2]]
      if(any(!(species %in% species_names)))
        stop("Species specified are not found in the model fit object")
      species <- match(species, species_names)
    }
    post_res <- lapply(post_res, "[",,species)
  }
  # Limit to sites that have been selected:
  if(!is.null(sites)){
    if(is.character(sites)){
      sites_names <- dimnames(post_res[[1]])[[1]]
      if(any(!(sites %in% sites_names)))
        stop("Sites specified are not found in the model fit object")
      sites <- match(sites, sites_names)
    }
    post_res <- lapply(post_res, "[",sites,)
  }

  # calculate summary statistic over sites:
  res <- lapply(post_res,function(x) apply(x,switch(calc_over,"site" = 1,
                                                    "species" = 2),stat_fun))

  if(simplify){
    res <- do.call(rbind, res)
  }

  return(res)


}
