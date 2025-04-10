% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_predict.R
\name{posterior_zipred}
\alias{posterior_zipred}
\title{Access the posterior distribution of the linear predictor for zero-inflation
parameter}
\usage{
posterior_zipred(
  object,
  transform = FALSE,
  newdata = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  list_index = "draws"
)
}
\arguments{
\item{object}{The model object}

\item{transform}{Should the linear predictor be transformed using the
inverse-link function. The default is \code{FALSE}, in which case the
untransformed linear predictor is returned.}

\item{newdata}{New data, by default \code{NULL} and uses original data}

\item{ndraws}{Number of draws, by default the number of samples in the
posterior. Will be sampled randomly from the chains if fewer than the
number of samples.}

\item{draw_ids}{The IDs of the draws to be used, as a numeric vector}

\item{list_index}{Whether to return the output list indexed by the number of
draws (default), species, or site.}
}
\value{
A list of linear predictors. If list_index is \code{"draws"} (the
default) the list will have length equal to the number of draws with each
element of the list being a site x species matrix. If the list_index is
\code{"species"} the list will have length equal to the number of species
with each element of the list being a draws x sites matrix. If the
list_index is \code{"sites"} the list will have length equal to the number
of sites with each element of the list being a draws x species matrix.
}
\description{
Extract the posterior draws of the linear predictor for the zero-inflation
parameter, possibly transformed by the inverse-link function.
}
\seealso{
\code{\link[=posterior_predict.jsdmStanFit]{posterior_predict.jsdmStanFit()}}
}
