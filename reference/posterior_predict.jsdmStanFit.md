# Draw from the posterior predictive distribution

Draw from the posterior predictive distribution of the outcome.

## Usage

``` r
# S3 method for class 'jsdmStanFit'
posterior_predict(
  object,
  newdata = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  list_index = "draws",
  Ntrials = NULL,
  include_zi = TRUE,
  include_shp = TRUE,
  zi_newdata = NULL,
  shp_newdata = NULL,
  ...
)
```

## Arguments

- object:

  The model object

- newdata:

  New data, by default `NULL` and uses original data

- ndraws:

  Number of draws, by default the number of samples in the posterior.
  Will be sampled randomly from the chains if fewer than the number of
  samples.

- draw_ids:

  The IDs of the draws to be used, as a numeric vector

- list_index:

  Whether to return the output list indexed by the number of draws
  (default), species, or site.

- Ntrials:

  For the binomial distribution the number of trials, given as either a
  single integer which is assumed to be constant across sites or as a
  site-length vector of integers.

- include_zi:

  For the zero-inflated distributions, whether to include the
  zero-inflation in the prediction. Defaults to `TRUE`.

- include_shp:

  For the distributions where the shape parameter is responsive to
  environmental predictors, whether to include this effect in the
  prediction. Defaults to `TRUE`.

- zi_newdata:

  For the zero-inflated distributions, the data used to fit any
  parameter effect upon the zi parameter. Defaults to `NULL` and uses
  original data.

- shp_newdata:

  For the distributions where the shape parameter is responsive to
  environmental predictors, any updated data for this effect in the
  prediction. Defaults to `NULL` and uses original data.

- ...:

  Currently unused

## Value

A list of linear predictors. If list_index is `"draws"` (the default)
the list will have length equal to the number of draws with each element
of the list being a site x species matrix. If the list_index is
`"species"` the list will have length equal to the number of species
with each element of the list being a draws x sites matrix. If the
list_index is `"sites"` the list will have length equal to the number of
sites with each element of the list being a draws x species matrix.

## See also

[`posterior_linpred.jsdmStanFit()`](https://nerc-ceh.github.io/jsdmstan/reference/posterior_linpred.jsdmStanFit.md)
