# Plotting an ordination plot for a GLLVM model

This function takes a GLLVM model fit and plots an ordination plot with
a random (or specified) selection of draws

## Usage

``` r
ordiplot(
  object,
  choices = c(1, 2),
  type = "species",
  summary_stat = "median",
  ndraws = 0,
  draw_ids = NULL,
  size = c(2, 1),
  alpha = c(1, 0.5),
  shape = c(18, 16),
  geom = "point",
  errorbar_range = 0.75,
  errorbar_linewidth = 1
)
```

## Arguments

- object:

  The `jsdmStanFit` model object

- choices:

  Which latent variables to plot as dimensions, by default `c(1,2)`

- type:

  Whether to plot sites or species, default `"species"`.

- summary_stat:

  The summary statistic used to plot overall averages of the posterior
  sample. By default this is `"median"`, and `NULL` will result in no
  summary being included

- ndraws:

  How many individual draws to include in plot, by default `0`. Setting
  this to `0` will result in no individual draws being included

- draw_ids:

  Which draws to include in plot (overrides `ndraws`)

- size:

  The size of the points in the graph, specified as a two-element vector
  with the first being used for the summary points and the second the
  individual draws, default `c(2,1)`

- alpha:

  The transparency/alpha of the points in the graph, specified as a
  two-element vector with the first being used for the summary points
  and the second the individual draws, default `c(1,0.5)`

- shape:

  The shape of the points in the graph, specified as a two-element
  vector with the first being used for the summary points and the second
  the individual draws, default `c(18,16)`

- geom:

  Which geom from ggplot2 is used for the summary statistic by default
  `"point"`, or alternatively can be `"text"`

- errorbar_range:

  If specified, the central range of the data that is covered by the
  errorbar. Needs to be given as either `NULL` (i.e. no errorbar), or a
  number between 0 and 1, default `0.75`.

- errorbar_linewidth:

  The linewidth of the error bar, by default 1.

## Value

A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object
that can be customised using the ggplot2 package

## Examples

``` r
if (FALSE) { # \dontrun{
# First simulate data and fit the model:
gllvm_data <- jsdm_sim_data(
  method = "gllvm", N = 100, S = 6, D = 3,
  family = "bernoulli"
)
gllvm_fit <- stan_jsdm(
  dat_list = gllvm_data, method = "gllvm",
  family = "bernoulli"
)

ordiplot(gllvm_fit)
# now plot the 1st and 3rd latent variables against each other for the sites:
ordiplot(gllvm_fit, choices = c(1, 3), type = "sites")
} # }
```
