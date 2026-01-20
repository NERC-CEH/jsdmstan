# Create prior object for jsdmstan models and data simulation

This function creates all of the potential priors used within a
`jsdmstan` model and can be used as the input to all
[`stan_jsdm()`](https://nerc-ceh.github.io/jsdmstan/reference/stan_jsdm.md)
family of functions and the
[`jsdm_sim_data()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_sim_data.md)
functions.

## Usage

``` r
jsdm_prior(
  sigmas_preds = "normal(0,1)",
  z_preds = "normal(0,1)",
  cor_preds = "lkj_corr(1)",
  betas = "normal(0,1)",
  sigma_a = "normal(0,1)",
  sigmas_species = "normal(0,1)",
  cor_species = "lkj_corr(1)",
  cor_species_chol = "lkj_corr_cholesky(1)",
  LV = "normal(0,1)",
  L = "normal(0,1)",
  sigma_L = "normal(0,1)",
  sigma = "normal(0,1)",
  kappa = "normal(0,1)",
  shape = "gamma(1,1)",
  zi = "beta(1,1)",
  zi_betas = "normal(0,1)",
  shp_betas = "normal(0,1)"
)

# S3 method for class 'jsdmprior'
print(x, ...)
```

## Arguments

- sigmas_preds:

  The standard deviation of the covariate effects, constrained to be
  positive (default standard normal)

- z_preds:

  The covariate effects (default standard normal)

- cor_preds:

  The correlation matrix on the covariate effects (npred by npred
  correlation matrix) (default `"lkj_corr(1)"`)

- betas:

  If covariate effects are unstructured, the prior on the covariate
  effects

- sigma_a:

  The standard deviation of the site level intercepts, constrained to be
  positive and default prior is half standard normal

- sigmas_species:

  For MGLMM method, the standard deviations of the species covariances,
  constrained to be positive (default half standard normal)

- cor_species:

  For MGLMM method (data simulation only), the correlation between
  species as the nspecies by nspecies correlation matrix (default
  `"lkj_corr(1)"`)

- cor_species_chol:

  For MGLMM method (model fit only), the correlation between species
  represented as a Cholesky decomposition of the nspecies by nspecies
  correlation matrix (default `"lkj_corr_cholesky(1)"`). Note that
  lkj_corr_cholesky(eta) is the cholesky decomposition of lkj_corr(eta)

- LV:

  For GLLVM method, the per site latent variable loadings (default
  standard normal)

- L:

  For GLLVM method, the non-zero species latent variable loadings
  (default standard normal)

- sigma_L:

  For GLLVM method, the variance of the species loadings, constrained to
  be positive (default standard normal)

- sigma:

  For Gaussian response, the standard deviation parameter. Constrained
  to be positive (default standard normal)

- kappa:

  For negative binomial response, the negative binomial variance
  parameter. Constrained to be positive (default standard normal)

- shape:

  For gamma response, the shape parameter. Constrained to be positive
  (default gamma(0.1,0.1))

- zi:

  For zero-inflated poisson or negative binomial with no environmental
  covariate effects upon the zero-inflation, the proportion of inflated
  zeros (default beta distribution with both alpha and beta parameters
  set to 1).

- zi_betas:

  For zero-inflated poisson or negative binomial with environmental
  effects upon the zero-inflation, the covariate effects on the
  zero-inflation on the logit scale

- shp_betas:

  For gaussian or negative binomial with environmental effects upon the
  family parameter (i.e. sigma or kappa), the covariate effects on the
  family parameter on the log scale

- x:

  Object of class `jsdmprior`

- ...:

  Currently unused

## Value

An object of class `"jsdmprior"` taking the form of a named list

## Details

Each prior has to be specified as a character string corresponding to
the appropriate stan command. The most common versions of these are
supported by the simulated data functions, however there are functions
that can be fed to the stan fitting procedure that will not be able to
be used as input for
[`jsdm_sim_data()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_sim_data.md).
Parameters `sigmas_preds`, `sigma_a`, `sigmas_species`, `sigma_L`,
`sigma`, and `kappa` are fixed to be positive only in the stan code and
this cannot be changed. Parameters `cor_preds` and `cor_species` are
assumed to be the Cholesky factor of a correlation matrix. All other
parameters are real numbers. For all parameters that represent vectors
or matrices the prior has to be the same across the entire vector or
matrix (note that for the species latent variable loadings in a GLLVM
model the prior is set on the non-zero matrix components `L` and not on
the entire matrix).

Prior distributions supported by
[`jsdm_sim_data()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_sim_data.md)
are `"normal(mean, sd)"`, `"student_t(df, mu, sigma)"`,
`"cauchy(location, scale)"`, `"gamma(shape, scale)"`,
`"inv_gamma(shape, scale)"` and `"lkj_corr_cholesky(eta)"`.

## Functions

- `print(jsdmprior)`: Print method for object of class `jsdmprior`

## See also

[sim_helpers](https://nerc-ceh.github.io/jsdmstan/reference/sim_helpers.md)
for a description of the parameterisations used within the data
simulation functions

## Examples

``` r
pr <- jsdm_prior(kappa = "gamma(0.01,0.01)")
pr
#>           Parameter             Group      Constraint                Prior
#> 1      sigmas_preds covariate_effects         lower=0          normal(0,1)
#> 2           z_preds covariate_effects            none          normal(0,1)
#> 3         cor_preds covariate_effects            none          lkj_corr(1)
#> 4             betas covariate_effects            none          normal(0,1)
#> 5           sigma_a    site_intercept            none          normal(0,1)
#> 6    sigmas_species             mglmm         lower=0          normal(0,1)
#> 7       cor_species             mglmm            none          lkj_corr(1)
#> 8  cor_species_chol             mglmm            none lkj_corr_cholesky(1)
#> 9                LV             gllvm            none          normal(0,1)
#> 10                L             gllvm            none          normal(0,1)
#> 11          sigma_L             gllvm         lower=0          normal(0,1)
#> 12            sigma          gaussian         lower=0          normal(0,1)
#> 13            kappa      neg_binomial         lower=0     gamma(0.01,0.01)
#> 14            shape             gamma         lower=0           gamma(1,1)
#> 15               zi    zero_inflation lower=0,upper=1            beta(1,1)
#> 16         zi_betas    zero_inflation            none          normal(0,1)
#> 17        shp_betas            family            none          normal(0,1)
```
