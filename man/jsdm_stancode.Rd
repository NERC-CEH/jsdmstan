% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jsdm_stancode.R
\name{jsdm_stancode}
\alias{jsdm_stancode}
\alias{print.jsdmstan_model}
\title{Make stancode for the jsdm model}
\usage{
jsdm_stancode(
  method,
  family,
  prior = jsdm_prior(),
  site_intercept = "none",
  beta_param = "cor",
  zi_param = "constant",
  shp_param = "constant"
)

\method{print}{jsdmstan_model}(x, ...)
}
\arguments{
\item{method}{The method, one of \code{"gllvm"} or \code{"mglmm"}}

\item{family}{is the response family, must be one of \code{"gaussian"},
\code{"neg_binomial"}, \code{"poisson"}, \code{"binomial"},
\code{"bernoulli"}, \code{"zi_poisson"}, or
\code{"zi_neg_binomial"}. Regular expression
matching is supported.}

\item{prior}{The prior, given as the result of a call to \code{\link[=jsdm_prior]{jsdm_prior()}}}

\item{site_intercept}{Whether a site intercept should be included, potential
values \code{"none"} (no site intercept), \code{"grouped"} (a site intercept
with hierarchical grouping) or \code{"ungrouped"} (site intercept with no
grouping)}

\item{beta_param}{The parameterisation of the environmental covariate effects, by
default \code{"cor"}. See details for further information.}

\item{zi_param}{For the zero-inflated families, whether the zero-inflation parameter
is a species-specific constant (default, \code{"constant"}), or varies by
environmental covariates (\code{"covariate"}).}

\item{shp_param}{For the families with shape parameters, whether the shape parameter
is a species-specific constant (default, \code{"constant"}), or varies by
environmental covariates (\code{"covariate"}).}

\item{x}{The jsdm_stancode object}

\item{...}{Currently unused}
}
\value{
A character vector of Stan code, class "jsdmstan_model"
}
\description{
This function returns the Stan code used to fit the model as specified by the data
list, family and method.
}
\details{
Environmental covariate effects (\code{"betas"}) can be parameterised in two
ways. With the \code{"cor"} parameterisation all covariate effects are assumed
to be constrained by a correlation matrix between the covariates. With the
\code{"unstruct"} parameterisation all covariate effects are assumed to draw
from a simple distribution with no correlation structure. Both parameterisations
can be modified using the prior object.
}
\section{Functions}{
\itemize{
\item \code{print(jsdmstan_model)}: A printing function for jsdmstan_model objects

}}
\examples{
jsdm_stancode(family = "gaussian", method = "gllvm")
jsdm_stancode(family = "poisson", method = "mglmm")

}
