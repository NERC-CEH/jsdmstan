% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jsdmstanfit-class.R
\name{summary.jsdmStanFit}
\alias{summary.jsdmStanFit}
\title{Summarise the model fit and data structure and give summaries for the parameters}
\usage{
\method{summary}{jsdmStanFit}(
  object,
  prob_quantiles = c(0.15, 0.85),
  digit_summary = 3,
  prob_pars_only = FALSE,
  pars = NULL,
  na_filter = TRUE,
  log_lik = FALSE,
  ...
)
}
\arguments{
\item{object}{The model object}

\item{prob_quantiles}{The quantiles to summarise the parameter estimates with, by
default the 15\% and 85\% quantiles}

\item{digit_summary}{The number of digits to round the results to}

\item{prob_pars_only}{Whether to limit output to parameters with Rhat > 1.01 or
effective sample size < 500, by default \code{FALSE}}

\item{pars}{Parameters to compute the summary of, by default \code{NULL} i.e. all
parameters included}

\item{na_filter}{Whether to remove parameters with NAs in Rhat - this includes the
parameters fixed to zero or one, such as the upper triangle of the cholesky
factor of the correlation matrix. By default \code{TRUE}}

\item{log_lik}{Whether the log_lik parameters should be included, default
\code{FALSE}}

\item{...}{Arguments passed to \code{\link[=extract]{extract()}}}
}
\description{
This returns a matrix of parameter summaries including a summary of the parameter
estimates, R-hat, bulk ESS and tail ESS. This can be limited to parameters with
Rhat > 1.01 or ESS < 500 by setting \code{prob_pars_only = TRUE}.
}
\examples{
\dontrun{

gllvm_data <- jsdm_sim_data(
  method = "gllvm", N = 100, S = 6, D = 2,
  family = "bernoulli"
)
gllvm_fit <- stan_jsdm(
  dat_list = gllvm_data, method = "gllvm",
  family = "bernoulli"
)
gllvm_summ <- summary(gllvm_fit)
head(gllvm_summ, 20)

summary(gllvm_fit, prob_quantiles = c(0.05, 0.5, 0.95))
}

}
