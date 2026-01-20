# Access the posterior distribution of the linear predictor

Extract the posterior draws of the linear predictor, possibly
transformed by the inverse-link function.

## Usage

``` r
# S3 method for class 'jsdmStanFit'
posterior_linpred(
  object,
  transform = FALSE,
  newdata = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  list_index = "draws",
  ...
)
```

## Arguments

- object:

  The model object

- transform:

  Should the linear predictor be transformed using the inverse-link
  function. The default is `FALSE`, in which case the untransformed
  linear predictor is returned.

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

- ...:

  Currently unused

## Value

A list of linear predictors. If list_index is `"draws"` (the default)
the list will have length equal to the number of draws with each element
of the list being a site x species matrix. If the list_index is
`"species"` the list will have length equal to the number of species
with each element of the list being a draws x sites matrix. If the
list_index is `"sites"` the list will have length equal to the number of
sites with each element of the list being a draws x species matrix. Note
that in the zero-inflated case this is only the linear predictor of the
non-zero-inflated part of the model.

## See also

[`posterior_predict.jsdmStanFit()`](https://nerc-ceh.github.io/jsdmstan/reference/posterior_predict.jsdmStanFit.md)
