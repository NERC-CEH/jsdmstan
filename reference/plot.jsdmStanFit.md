# Plot the traceplots and density plots for parameters within a jsdmStanFit object

This function takes parameters from the jsdmStanfit object and plots
them using the
[`bayesplot::mcmc_combo()`](https://mc-stan.org/bayesplot/reference/MCMC-combos.html)function
from the bayesplot package.

## Usage

``` r
# S3 method for class 'jsdmStanFit'
plot(
  x,
  pars = NULL,
  combo = c("dens", "trace"),
  N = 5L,
  ask = TRUE,
  inc_warmup = FALSE,
  include = TRUE,
  sample_n = 10,
  regexp = FALSE,
  plot = TRUE,
  newpage = TRUE,
  ...
)
```

## Arguments

- x:

  The `jsdmStanFit` model object

- pars:

  The parameters to plot, by default a random sample of twenty of the
  parameters fit within the model

- combo:

  Which combination of plot types within
  [`bayesplot::mcmc_combo()`](https://mc-stan.org/bayesplot/reference/MCMC-combos.html)
  to use, by default `c("dens", "trace")`

- N:

  The number of plots per page, default `5`

- ask:

  Whether to ask before plotting a new page, default `TRUE`

- inc_warmup:

  Whether to include the warmup period in the plots, by default `FALSE`

- include:

  Whether to include or exclude the parameters specified by pars, by
  default `TRUE` (i.e. include)

- sample_n:

  If `pars = NULL` then the number of random non-sigma parameters to
  include (details in description)

- regexp:

  If pars should be treated as a regular expression for matching to
  parnames, by default `FALSE`

- plot:

  Whether to plot the plots, default `TRUE`

- newpage:

  Whether the first plot should be plotted on a new page, by default
  `TRUE`

- ...:

  Arguments passed to
  [`bayesplot::mcmc_combo()`](https://mc-stan.org/bayesplot/reference/MCMC-combos.html)

## Value

An invisible list of the plots#

## Details

This acts as an interface to the
[`bayesplot::mcmc_combo()`](https://mc-stan.org/bayesplot/reference/MCMC-combos.html)
function, with the default being to plot a density plot and a trace plot
for each parameter specified, although this can be changed by setting
the `combo` argument. These jSDM models have a lot of parameters so
selecting a subset is recommended. If pars is set to `NULL` (the
default) all parameters with either sigma or kappa in their name will be
plotted along with a random selection of the other parameters (total
number of other parameters set by `sample_n`). To see the name of the
parameters within the model use
[`get_parnames()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdmstan-extractors.md) -
and if you want to plot all parameters (there will be hundreds in any
reasonably sized model) set `pars = get_parnames(x)`.

## See also

[`mcmc_plot.jsdmStanFit()`](https://nerc-ceh.github.io/jsdmstan/reference/mcmc_plot.jsdmStanFit.md)
for more plotting options.

## Examples

``` r
if (FALSE) { # \dontrun{
# First simulate data and get model fit:
mglmm_data <- mglmm_sim_data(
  N = 100, S = 10, K = 3,
  family = "gaussian"
)
mglmm_fit <- stan_mglmm(
  Y = mglmm_data$Y, X = mglmm_data$X,
  family = "gaussian"
)

# The default plot:
plot(mglmm_fit)

# Plotting specifically the cor_species parameters:
plot(mglmm_fit, pars = "cor_species", regexp = TRUE)

# Increasing the number of randomly sampled parameters to plot:
plot(mglmm_fit, sample_n = 20)
} # }
```
