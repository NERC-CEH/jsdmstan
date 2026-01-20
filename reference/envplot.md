# Plotting environmental effects on species

Plotting environmental effects on species

## Usage

``` r
envplot(
  object,
  include_intercept = FALSE,
  nrow = NULL,
  y_labels = NULL,
  plotfun = "intervals",
  widths = NULL,
  species = NULL,
  preds = NULL,
  ...
)
```

## Arguments

- object:

  The jsdmStanFit model object

- include_intercept:

  Whether to include the intercept in the plots

- nrow:

  The number of rows within the plot

- y_labels:

  Which plots should have annotated y axes. Needs to be given as an
  integer vector

- plotfun:

  Which plot function from mcmc_plot should be used, by default
  `"intervals"`

- widths:

  The widths of the plots

- species:

  The species to be included in the plot, by default `NULL` and all
  species are included. Should be specified as a vector of species names
  or numeric indices.

- preds:

  The predictors to be included in the plot, by default `NULL` and all
  predictors are included. Should be specified as a vector of predictor
  names or numeric indices.

- ...:

  Other arguments that are passed to the MCMC plotting function (see
  [`mcmc_plot.jsdmStanFit()`](https://nerc-ceh.github.io/jsdmstan/reference/mcmc_plot.jsdmStanFit.md))

## Value

An object of class `"bayesplot_grid"`, for more information see
[`bayesplot::bayesplot_grid()`](https://mc-stan.org/bayesplot/reference/bayesplot_grid.html)
