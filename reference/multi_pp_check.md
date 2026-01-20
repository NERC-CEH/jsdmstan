# Multiple pp_check plots per species

Multiple pp_check plots per species

## Usage

``` r
multi_pp_check(
  object,
  plotfun = "dens_overlay",
  species = NULL,
  ndraws = NULL,
  xlim = NULL,
  ylim = NULL,
  titles = character(),
  subtitles = character(),
  legends = TRUE,
  save_gg_objects = TRUE,
  grid_args = list(...),
  ...
)
```

## Arguments

- object:

  The jsdmStanFit model object

- plotfun:

  The ppc plot function to use, given as a character string. The default
  is to call
  [ppc_dens_overlay](https://mc-stan.org/bayesplot/reference/PPC-distributions.html).
  Can be specified as either the entire name of function as a character
  string or without the ppc\_ prefix.

- species:

  Which species should be included, by default all

- ndraws:

  How many draws should be used within the plots

- xlim, ylim:

  Optionally, numeric vectors of length 2 specifying lower and upper
  limits for the axes that will be shared across all plots.

- titles, subtitles:

  Optional character vectors of plot titles and subtitles. If specified,
  `titles` and `subtitles` must must have length equal to the number of
  plots specified.

- legends:

  If any of the plots have legends should they be displayed? Defaults to
  `TRUE`.

- save_gg_objects:

  If `TRUE`, the default, then the ggplot objects specified in `...` or
  via the `plots` argument are saved in a list in the `"bayesplots"`
  component of the returned object. Setting this to `FALSE` will make
  the returned object smaller but these individual plot objects will not
  be available.

- grid_args:

  An optional named list of arguments to pass to
  [`gridExtra::arrangeGrob()`](https://rdrr.io/pkg/gridExtra/man/arrangeGrob.html)
  (`nrow`, `ncol`, `widths`, etc.).

- ...:

  Other options passed to the bayesplot PPC function or
  posterior_predict. Note that the same value will be used for every
  species included (e.g. binwidth for a histogram cannot be specified by
  species).

## Value

An object of class `"bayesplot_grid"`, for more information see
[`bayesplot::bayesplot_grid()`](https://mc-stan.org/bayesplot/reference/bayesplot_grid.html)
