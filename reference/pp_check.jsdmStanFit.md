# Posterior predictive checks for `jsdmStanFit` objects

This function acts as an interface to
[`bayesplot::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html),
by default calculating summary statistics for each site (i.e. row in the
response matrix) which are then plotted and compared to the data used to
fit the model.

## Usage

``` r
# S3 method for class 'jsdmStanFit'
pp_check(
  object,
  plotfun = "dens_overlay",
  species = NULL,
  sites = NULL,
  summary_stat = "sum",
  calc_over = "sites",
  ndraws = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `jsdmStanFit`

- plotfun:

  The ppc plot function to use, given as a character string. The default
  is to call
  [ppc_dens_overlay](https://mc-stan.org/bayesplot/reference/PPC-distributions.html).
  Can be specified as either the entire name of function as a character
  string or without the ppc\_ prefix. If `plotfun == "pairs"` then a
  pairs plot is produced with a diagonal of density plots for the
  selected species, the upper triangle showing the recovery of the
  correlation between the selected species, and the lower triangle
  showing the plotted relationships between the species in the data and
  one draw of the posterior prediction.

- species:

  Which species to include in the summary statistic, by default all

- sites:

  Which sites to include in the summary statistic, by default all

- summary_stat:

  The summary statistic to be used, by default `sum` but any function
  can be used.

- calc_over:

  Whether to calculate the summary statistic by sites or species, by
  default `species`

- ndraws:

  Number of draws, by default the number of samples in the posterior.
  Will be sampled randomly from the chains if fewer than the number of
  samples.

- ...:

  Additional arguments passed to
  [`jsdm_statsummary()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_statsummary.md).

## Value

A ggplot object that can be further customised using the ggplot2
package.

## Details

This function takes a `jsdmStanFit` object and tries to extract
statistics that can give useful summaries within a posterior prediction
(or retrodiction) using the
[`bayesplot::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)
function. It uses the
[`jsdm_statsummary()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_statsummary.md)
function to get summary statistics and then supplies them to the
specified `ppc_` function from the bayesplot package. For a full list of
available plotting functions call
[`bayesplot::available_ppc()`](https://mc-stan.org/bayesplot/reference/available_ppc.html).

## Examples

``` r
if (FALSE) { # \dontrun{
# First simulate data and fit the jsdmStan model:
mglmm_data <- mglmm_sim_data(
  N = 100, S = 10, K = 3,
  family = "gaussian"
)
mglmm_fit <- stan_mglmm(
  Y = mglmm_data$Y, X = mglmm_data$X,
  family = "gaussian"
)

# The default is to plot a density overlay:
pp_check(mglmm_fit)

# Other plot functions can be called, such as a ribbon plot:
pp_check(mglmm_fit, plotfun = "ribbon")

# Instead of calculating the sum over sites other statistics can be calculated,
# e.g. the mean of each species:
pp_check(mglmm_fit,
  plotfun = "ecdf_overlay", summary_stat = "mean",
  calc_over = "species", ndraws = 20
)

# A pairs plot - limiting to only a subset of species for graphical simplicity
pp_check(mglmm_fit, plotfun = "pairs", species = 1:4)
} # }
```
