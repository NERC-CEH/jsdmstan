# Plot modelled correlations between species

Plot modelled correlations between species

## Usage

``` r
corrplot(
  object,
  species = NULL,
  plotfun = "intervals",
  nrow = NULL,
  widths = NULL,
  ...
)
```

## Arguments

- object:

  The jsdmStanFit model object

- species:

  Which species should be included - this plots the correlations between
  the specified species and all other species, so by default if all
  species are included this will result in a plot per species with all
  other species there and thus duplicate entries across all the plots

- plotfun:

  Which plotting function from bayesplot to use, by default
  `"intervals"`

- nrow:

  How many rows the grid of plots has

- widths:

  Whether the widths of the different rows of plots should vary

- ...:

  Other arguments passed to mcmc_plot

## Value

An object of class `"bayesplot_grid"`, for more information see
[`bayesplot::bayesplot_grid()`](https://mc-stan.org/bayesplot/reference/bayesplot_grid.html)
