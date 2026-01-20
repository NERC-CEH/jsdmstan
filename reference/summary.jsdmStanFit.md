# Summarise the model fit and data structure and give summaries for the parameters

This returns a matrix of parameter summaries including a summary of the
parameter estimates, R-hat, bulk ESS and tail ESS. This can be limited
to parameters with Rhat \> 1.01 or ESS \< 500 by setting
`prob_pars_only = TRUE`.

## Usage

``` r
# S3 method for class 'jsdmStanFit'
summary(
  object,
  prob_quantiles = c(0.15, 0.85),
  digit_summary = 3,
  prob_pars_only = FALSE,
  pars = NULL,
  na_filter = TRUE,
  log_lik = FALSE,
  ...
)
```

## Arguments

- object:

  The model object

- prob_quantiles:

  The quantiles to summarise the parameter estimates with, by default
  the 15% and 85% quantiles

- digit_summary:

  The number of digits to round the results to

- prob_pars_only:

  Whether to limit output to parameters with Rhat \> 1.01 or effective
  sample size \< 500, by default `FALSE`

- pars:

  Parameters to compute the summary of, by default `NULL` i.e. all
  parameters included

- na_filter:

  Whether to remove parameters with NAs in Rhat - this includes the
  parameters fixed to zero or one, such as the upper triangle of the
  cholesky factor of the correlation matrix. By default `TRUE`

- log_lik:

  Whether the log_lik parameters should be included, default `FALSE`

- ...:

  Arguments passed to
  [`extract()`](https://nerc-ceh.github.io/jsdmstan/reference/extract.jsdmStanFit.md)

## Examples

``` r
if (FALSE) { # \dontrun{

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
} # }
```
