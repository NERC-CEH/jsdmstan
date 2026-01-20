# Update a jsdmStanFit model object with new data or Stan arguments

This function allows you to update a jsdmStanFit model with new data or
new arguments to Stan. It does not support changes that require
recompiling stancode - for that you should use
[`stan_jsdm()`](https://nerc-ceh.github.io/jsdmstan/reference/stan_jsdm.md).
Changes to the number of sites, species or covariates do not require
recompiling stancode and can therefore be done using this function.

## Usage

``` r
# S3 method for class 'jsdmStanFit'
update(
  object,
  newY = NULL,
  newX = NULL,
  newD = NULL,
  newNtrials = NULL,
  newZi_X = NULL,
  newShp_X = NULL,
  save_data = TRUE,
  ...
)
```

## Arguments

- object:

  The jsdmStanFit model object

- newY:

  New Y data, by default `NULL`

- newX:

  New X data, by default `NULL`

- newD:

  New number of latent variables, by default `NULL`

- newNtrials:

  New number of trials (binomial model only), by default `NULL`

- newZi_X:

  New predictor data for the zi parameter in zero-inflated models, by
  default `NULL`. In cases where the model was originally fit with the
  same X and zi_X data and only newX is supplied to update.jsdmStanFit
  the zi_X data will also be set to newX.

- newShp_X:

  New predictor data for the family parameter in models where the family
  parameter is modelled in response to data, by default `NULL`. In cases
  where the model was originally fit with the same X and shp_X data and
  only newX is supplied to update.jsdmStanFit the shp_X data will also
  be set to newX.

- save_data:

  Whether to save the data in the jsdmStanFit object, by default `TRUE`

- ...:

  Arguments passed to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)

## Value

An object of class `jsdmStanFit`

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
mglmm_fit2 <- update(mglmm_fit, iter = 4000)

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
gllvm_data <- jsdm_sim_data(
  method = "gllvm", N = 500, S = 4, D = 2,
  family = "bernoulli"
)
gllvm_fit2 <- update(gllvm_fit, newY = gllvm_data$Y)
gllvm_fit2
} # }
```
