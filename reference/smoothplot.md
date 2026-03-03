# Plot the smooth effect of any spline within a model

Plot the smooth effect of any spline within a model

## Usage

``` r
smoothplot(
  object,
  ndraws = 20,
  alpha = 0.4,
  draw_ids = NULL,
  summarise = NULL,
  type = "link",
  facet_args = list(),
  select_smooths = NULL
)
```

## Arguments

- object:

  The jsdmStanFit model object

- ndraws:

  The number of draws of the predictor

- alpha:

  The alpha (transparency) value of the lines

- draw_ids:

  Which draws to include in plot (overrides `ndraws`)

- summarise:

  Whether to include a summary, by default `NULL` i.e. no summary. Can
  also be given as the name of a function (e.g. `"mean"`). Note that
  this summarisation only occurs over the number of draws specified.

- type:

  Whether to return plot on the link scale (default) or the response
  scale

- facet_args:

  Argument passed to
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html),
  specified as a list

- select_smooths:

  Which plots to include, given either as `NULL` (i.e. plot all smooths)
  or as a list containing two numeric vectors named "site" and "species"
  representing the indices of the site-specific smooths to plot and the
  species-specific smooths to plot

## Value

a list of ggplot objects
