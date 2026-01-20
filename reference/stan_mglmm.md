# Alias for `stan_jsdm` with `method = "mglmm"`

Alias for `stan_jsdm` with `method = "mglmm"`

## Usage

``` r
stan_mglmm(X, ...)

# Default S3 method
stan_mglmm(
  X = NULL,
  Y = NULL,
  species_intercept = TRUE,
  dat_list = NULL,
  family,
  site_intercept = "none",
  prior = jsdm_prior(),
  save_data = TRUE,
  iter = 4000,
  ...
)

# S3 method for class 'formula'
stan_mglmm(formula, data = list(), ...)
```

## Arguments

- X:

  The covariates matrix, with rows being site and columns being
  covariates. Ignored in favour of data when formula approach is used to
  specify model.

- ...:

  Arguments passed to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)

- Y:

  Matrix of species by sites. Rows are assumed to be sites, columns are
  assumed to be species

- species_intercept:

  Whether the model should be fit with an intercept by species, by
  default `TRUE`

- dat_list:

  Alternatively, data can be given to the model as a list containing Y,
  X, N, S, K, and site_intercept. See output of
  [`jsdm_sim_data()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_sim_data.md)
  for an example of how this can be formatted.

- family:

  is the response family, must be one of `"gaussian"`, `"neg_binomial"`,
  `"poisson"`, `"binomial"`, `"bernoulli"`, or `"zi_poisson"`. Regular
  expression matching is supported.

- site_intercept:

  Whether a site intercept should be included, potential values `"none"`
  (no site intercept), `"grouped"` (a site intercept with hierarchical
  grouping) or `"ungrouped"` (site intercept with no grouping)

- prior:

  Set of prior specifications from call to
  [`jsdm_prior()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_prior.md)

- save_data:

  If the data used to fit the model should be saved in the model object,
  by default TRUE.

- iter:

  A positive integer specifying the number of iterations for each chain,
  default 4000.

- formula:

  The formula of covariates that the species means are modelled from

- data:

  Dataframe or list of covariates.

## Methods (by class)

- `stan_mglmm(default)`: Default

- `stan_mglmm(formula)`: Formula interface
