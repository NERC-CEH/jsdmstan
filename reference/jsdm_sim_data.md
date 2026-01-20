# Generate simulated data within a variety of jSDM methodologies

The `jsdm_sim_data` function can simulate data with either a
multivariate generalised mixed model (MGLMM) or a generalised linear
latent variable model (GLLVM). The `gllvm_sim_data` and `mglmm_sim_data`
are aliases for `jsdm_sim_data` that set `method` to `"gllvm"` and
`"mglmm"` respectively.

## Usage

``` r
jsdm_sim_data(
  S,
  N = NULL,
  D = NULL,
  K = 0L,
  family,
  method = c("gllvm", "mglmm"),
  species_intercept = TRUE,
  Ntrials = NULL,
  site_intercept = "none",
  beta_param = "unstruct",
  zi_param = "constant",
  zi_k = NULL,
  shp_param = "constant",
  shp_k = NULL,
  censoring = "none",
  censor_points = NULL,
  prior = jsdm_prior(),
  X = NULL,
  zi_X = NULL,
  shp_X = NULL
)

gllvm_sim_data(...)

mglmm_sim_data(...)
```

## Arguments

- S:

  is number of species

- N:

  is number of sites

- D:

  is number of latent variables, used within gllvm method

- K:

  is number of covariates, by default `0`

- family:

  is the response family, must be one of `"gaussian"`, `"neg_binomial"`,
  `"poisson"`, `"binomial"`, `"bernoulli"`, `"zi_poisson"`, or
  `"zi_neg_binomial"`. Regular expression matching is supported.

- method:

  is the jSDM method to use, currently either `"gllvm"` or `"mglmm"` -
  see details for more information.

- species_intercept:

  Whether to include an intercept in the predictors, must be `TRUE` if
  `K` is `0`. Defaults to `TRUE`.

- Ntrials:

  For the binomial distribution the number of trials, given as either a
  single integer which is assumed to be constant across sites or as a
  site-length vector of integers.

- site_intercept:

  Whether a site intercept should be included, potential values `"none"`
  (no site intercept) or `"ungrouped"` (site intercept with no
  grouping). Defaults to no site intercept, grouped is not supported
  currently.

- beta_param:

  The parameterisation of the environmental covariate effects, by
  default `"unstruct"`. See details for further information.

- zi_param:

  For the zero-inflated families, whether the zero-inflation parameter
  is a species-specific constant (default, `"constant"`), or varies by
  environmental covariates (`"covariate"`).

- zi_k:

  If `zi_param="covariate"`, the number of environmental covariates that
  the zero-inflation parameter responds to. The default (`NULL`) is that
  the zero-inflation parameter responds to exactly the same covariate
  matrix as the mean parameter. Otherwise, a different set of random
  environmental covariates are generated, plus an intercept (not
  included in zi_k) and used to predict zero-inflation. Will be ignored
  if zi_X is supplied.

- shp_param:

  For families with shape parameters, whether the shape parameter is a
  species-specific constant (default, `"constant"`), or varies by
  environmental covariates (`"covariate"`).

- shp_k:

  If `shp_param="covariate"`, the number of environmental covariates
  that the shape parameter responds to. The default (`NULL`) is that the
  shape parameter responds to exactly the same covariate matrix as the
  mean parameter. Otherwise, a different set of random environmental
  covariates are generated and used to predict the shape parameter. Will
  be ignored if shp_X is supplied.

- censoring:

  If the response is left-censored (`"left"`) or not censored (default,
  `"none"`).

- censor_points:

  The values at which censoring occurs, to be provided as either a
  S-length vector (if `censoring = "left" or "right"`) or a named list

- prior:

  Set of prior specifications from call to
  [`jsdm_prior()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_prior.md)

- X:

  The X matrix to be used to simulate data, by default `NULL` - i.e. the
  X matrix is simulated using a random draw from a standard normal
  distribution.

- zi_X:

  The zi_X matrix to be used to simulate data, by default `NULL` - i.e.
  the zi_X matrix is simulated using a random draw from a standard
  normal distribution.

- shp_X:

  The shp_X matrix to be used to simulate data, by default `NULL` - i.e.
  the shp_X matrix is simulated using a random draw from a standard
  normal distribution.

- ...:

  Arguments passed to jsdm_sim_data

## Details

This simulates data based on a joint species distribution model with
either a generalised linear latent variable model approach or a
multivariate generalised linear mixed model approach.

Models can be fit with or without "measured predictors", and if measured
predictors are included then the species have species-specific parameter
estimates. These can either be simulated completely independently, or
have information pooled across species. If information is pooled this
can be modelled as either a random draw from some mean and standard
deviation or species covariance can be modelled together (this will be
the covariance used in the overall model if the method used has
covariance).

Environmental covariate effects (`"betas"`) can be parameterised in two
ways. With the `"cor"` parameterisation all covariate effects are
assumed to be constrained by a correlation matrix between the
covariates. With the `"unstruct"` parameterisation all covariate effects
are assumed to draw from a simple distribution with no correlation
structure. Both parameterisations can be modified using the prior
object.

## Functions

- `gllvm_sim_data()`: Alias for `jsdm_sim_data` with `method = "gllvm"`

- `mglmm_sim_data()`: Alias for `jsdm_sim_data` with `method = "mglmm"`
