# Extract site- or species-level summary statistics for predictions from a `jsdmStanFit` model

This function extracts the predicted Y values for within the models and
then calculates summary statistics for each simulated community. The
default is to sum all the predicted scores for each site.

## Usage

``` r
jsdm_statsummary(
  object,
  species = NULL,
  sites = NULL,
  summary_stat = "sum",
  post_type = "linpred",
  calc_over = "sites",
  simplify = TRUE,
  draw_ids = NULL,
  ndraws = NULL,
  ...
)
```

## Arguments

- object:

  A `jsdmStanFit` model object

- species:

  Which species to include in the summary statistic, by default all

- sites:

  Which sites to include in the summary statistic, by default all

- summary_stat:

  The summary statistic to be used, by default `sum` but any function
  can be used.

- post_type:

  The type of posterior prediction to be used, either `"linpred"` for
  [`posterior_linpred.jsdmStanFit()`](https://nerc-ceh.github.io/jsdmstan/reference/posterior_linpred.jsdmStanFit.md)
  or `"predict"` for
  [`posterior_predict.jsdmStanFit()`](https://nerc-ceh.github.io/jsdmstan/reference/posterior_predict.jsdmStanFit.md)

- calc_over:

  Whether to calculate the summary statistic by sites or species, by
  default `species`

- simplify:

  Whether to simplify the output into a matrix, by default `TRUE`

- draw_ids:

  The IDs of the draws to be used, as a numeric vector

- ndraws:

  Number of draws, by default the number of samples in the posterior.
  Will be sampled randomly from the chains if fewer than the number of
  samples.

- ...:

  Arguments passed to the posterior prediction function

## Value

If `simplify = TRUE` then a matrix where each row is a draw and each
column is either a site (if `calc_over = "sites"`) or a species (if
`calc_over = "species"`).

## See also

pp_check.jsdmStanFit

## Examples

``` r
if (FALSE) { # \dontrun{
# First simulate data and fit the jsdmStan model:
gllvm_data <- gllvm_sim_data(
  N = 100, S = 9, D = 2,
  family = "bernoulli"
)
gllvm_fit <- stan_gllvm(dat_list = gllvm_data, family = "bernoulli")

# The default is to return a matrix:
jsdm_statsummary(gllvm_fit)

# The above returns the linear predictor, while we may want to get the posterior
# prediction instead:
jsdm_statsummary(gllvm_fit, post_type = "predict")

# This can be limited to a specific set of species and/or sites:
jsdm_statsummary(gllvm_fit, species = 1:5, sites = seq(5, 95, 10))
} # }
```
