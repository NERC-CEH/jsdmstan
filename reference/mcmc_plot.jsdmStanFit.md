# MCMC plots implemented in bayesplot

An interface for calling the MCMC plotting functions implemented in the
bayesplot package.

## Usage

``` r
# S3 method for class 'jsdmStanFit'
mcmc_plot(
  x,
  pars = NULL,
  plotfun = "intervals",
  sample_n = 10,
  inc_warmup = FALSE,
  include = TRUE,
  regexp = FALSE,
  ...
)

mcmc_plot(x, ...)
```

## Arguments

- x:

  The jsdmStanFit model object

- pars:

  The parameters to plot, by default a random sample of twenty of the
  parameters fit within the model

- plotfun:

  The MCMC plot function from bayesplot to be used

- sample_n:

  If `pars = NULL` then the number of random non-sigma parameters to
  include (details in description)

- inc_warmup:

  Whether to include the warmup period in the plots, by default `FALSE`

- include:

  Whether to include or exclude the parameters specified by pars, by
  default `TRUE` (i.e. include)

- regexp:

  If pars should be treated as a regular expression for matching to
  parnames, by default `FALSE`

- ...:

  Other arguments to be passed to the MCMC plotting function

## Value

A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object
that can be customised using the ggplot2 package

## Details

This acts as an interface to the plotting functions prefaced with
`mcmc_` within the
[bayesplot](https://mc-stan.org/bayesplot/reference/bayesplot-package.html)
package. The default is to plot an interval plot for the parameters
specified, for a full list of available plot options run
[`bayesplot::available_mcmc()`](https://mc-stan.org/bayesplot/reference/available_ppc.html)
or the error message when you set plotfun to an unrecognised plot
function will give you a list of options. These jSDM models have a lot
of parameters so selecting a subset is recommended. If pars is set to
`NULL` (the default) all parameters with either sigma or kappa in their
name will be plotted along with a random selection of the other
parameters (total number of other parameters set by `sample_n`). To see
the name of the parameters within the model use
[`get_parnames()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdmstan-extractors.md) -
and if you want to plot all parameters (there will be hundreds in any
reasonably sized model) set `pars = get_parnames(x)`.

## See also

[`plot.jsdmStanFit()`](https://nerc-ceh.github.io/jsdmstan/reference/plot.jsdmStanFit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# First simulate data and fit the model:
gllvm_data <- jsdm_sim_data(
  method = "gllvm", N = 100, S = 6, D = 2,
  family = "bernoulli"
)
gllvm_fit <- stan_jsdm(
  dat_list = gllvm_data, method = "gllvm",
  family = "bernoulli"
)

# Default is to plot the intervals:
mcmc_plot(gllvm_fit)

# Other plot types include options to see parameter recovery (if the
# jsdm_sim_data functions are used the original parameters are saved within the
# data list)
mcmc_plot(gllvm_fit,
  plotfun = "recover_intervals",
  pars = c("LV[2,20]", "sigmas_preds[1]", "sigma_L"),
  true = c(
    gllvm_data$pars$LV[2, 20],
    gllvm_data$pars$beta_sds,
    gllvm_data$pars$L_sigma
  )
)
} # }
```
