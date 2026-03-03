# Fit jsdm models in Stan

This function fits joint Species Distribution models in Stan, using
either a generalised linear latent variable model approach
(`method = "gllvm"`), or a multivariate generalised linear mixed model
approach (`method = "mglmm"`).

## Usage

``` r
stan_jsdm(X, ...)

# Default S3 method
stan_jsdm(
  X = NULL,
  Y = NULL,
  species_intercept = TRUE,
  method,
  dat_list = NULL,
  family,
  site_intercept = "none",
  D = NULL,
  prior = jsdm_prior(),
  site_groups = NULL,
  beta_param = "unstruct",
  Ntrials = NULL,
  censoring = "none",
  censor_points = NULL,
  cens_ID = NULL,
  zi_param = "constant",
  zi_X = NULL,
  shp_param = "constant",
  shp_X = NULL,
  spl_smooth = NULL,
  save_data = TRUE,
  iter = 4000,
  init = NULL,
  ...
)

# S3 method for class 'formula'
stan_jsdm(
  formula,
  data = list(),
  Y = NULL,
  zi_formula = NULL,
  shp_formula = NULL,
  zi_param = "constant",
  shp_param = "constant",
  ...
)
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

- method:

  Whether to fit a GLLVM or MGLMM model, details in description

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

- D:

  The number of latent variables within a GLLVM model

- prior:

  Set of prior specifications from call to
  [`jsdm_prior()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_prior.md)

- site_groups:

  If the site intercept is grouped, a vector of group identities per
  site

- beta_param:

  The parameterisation of the environmental covariate effects, by
  default `"unstruct"`. See details for further information.

- Ntrials:

  For the binomial distribution the number of trials, given as either a
  single integer which is assumed to be constant across sites or as a
  site-length vector of integers.

- censoring:

  If the response is left-censored (`"left"`) or not censored (default,
  `"none"`).

- censor_points:

  For left-censored models, the values at which censoring occurs, to be
  provided as a S-length vector.

- cens_ID:

  For left-censored models, an alternative way of declaring where
  censoring occurs where a N by S binary matrix is supplied with value 0
  where the measurement is uncensored and a value of 1 where the
  measurement is censored.

- zi_param:

  For the zero-inflated families, whether the zero-inflation parameter
  is a species-specific constant (default, `"constant"`), or varies by
  covariates (`"covariate"`). Set to `"covariate"` if `zi_formula` is
  specified.

- zi_X:

  If `zi_param = "covariate"`, the matrix of predictors that the
  zero-inflation is modelled in response to. If there is not already an
  intercept column (identified by all values being equal to one), one
  will be added to the front of the matrix. Overridden by `zi_formula`
  when formula approach is used.

- shp_param:

  For the families with shape parameters, whether the shape parameter is
  a species-specific constant (default, `"constant"`), or varies by
  covariates (`"covariate"`). Set to `"covariate"` if `shp_formula` is
  specified.

- shp_X:

  If `shp_param = "covariate"`, the matrix of predictors that the shape
  parameter is modelled in response to. If there is not already an
  intercept column (identified by all values being equal to one), one
  will be added to the front of the matrix. Overridden by `shp_formula`
  when formula approach is used.

- spl_smooth:

  The construction of the spline process, do not supply directly as user
  but instead include a spline within the formula interface as would be
  provided in [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html)

- save_data:

  If the data used to fit the model should be saved in the model object,
  by default TRUE.

- iter:

  A positive integer specifying the number of iterations for each chain,
  default 4000.

- init:

  Initialisation values for the sampling, see
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).
  If unspecified and `method = "mglmm"` then set to `"0"`.

- formula:

  The formula of covariates that the species means are modelled from

- data:

  Dataframe or list of covariates.

- zi_formula:

  For the zero-inflated families, the formula of any covariate effect
  upon the zi parameter. Covariates are sourced from the `data`
  argument. Only works if main effect is specified using formula
  argument.

- shp_formula:

  For the families with shape parameters, the formula of any covariate
  effect upon the shape parameter. Covariates are sourced from the
  `data` argument. Only works if main effect is specified using formula
  argument.

## Value

A `jsdmStanFit` object, comprising a list including the StanFit object,
the data used to fit the model plus a few other bits of information. See
[jsdmStanFit](https://nerc-ceh.github.io/jsdmstan/reference/jsdmStanFit.md)
for details.

## Details

Environmental covariate effects (`"betas"`) can be parameterised in two
ways. With the `"cor"` parameterisation all covariate effects are
assumed to be constrained by a correlation matrix between the
covariates. With the `"unstruct"` parameterisation all covariate effects
are assumed to draw from a simple distribution with no correlation
structure. Both parameterisations can be modified using the prior
object. Families supported are the Gaussian family, the negative
binomial family, the Poisson family, the binomial family (with number of
trials specificied using the `Ntrials` parameter), the Bernoulli family
(the special case of the binomial family where number of trials is equal
to one), the zero-inflated Poisson and the zero-inflated negative
binomial. For both zero-inflated families the zero-inflation is assumed
to be a species-specific constant.

## Methods (by class)

- `stan_jsdm(default)`: this is the default way of doing things

- `stan_jsdm(formula)`: Formula interface

## Examples

``` r
if (FALSE) { # \dontrun{
# MGLMM - specified by using the mglmm aliases and with direct reference to Y and
# X matrices:
mglmm_data <- mglmm_sim_data(
  N = 100, S = 10, K = 3,
  family = "gaussian"
)
mglmm_fit <- stan_mglmm(
  Y = mglmm_data$Y, X = mglmm_data$X,
  family = "gaussian"
)
mglmm_fit

# You can also run a model by supplying the data as a list:
gllvm_data <- jsdm_sim_data(
  method = "gllvm", N = 100, S = 6, D = 2,
  family = "bernoulli"
)
gllvm_fit <- stan_jsdm(
  dat_list = gllvm_data, method = "gllvm",
  family = "bernoulli"
)
gllvm_fit
} # }
```
