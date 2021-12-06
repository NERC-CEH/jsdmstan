#' jsdmStanFit fitted jsdm_stan model
#'
#' @slot jsdm_type A length one character vector describing type of jSDM
#'
#' @slot n_species A length one integer vector representing number of species
#'
#' @slot n_site A length one integer vector representing number of sites
#'
#' @slot n_preds A length two integer vector representing number of preds then
#'   whether intercept included
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
           site = NA_character_,
           preds = NA_character_,
           n_latent = NA_integer_
         )
         )
