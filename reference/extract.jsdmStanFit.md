# Extract samples from jsdmStanFit object

This function extracts named parameters from a jsdmStanFit object, with
optional regular expression matching.

## Usage

``` r
# S3 method for class 'jsdmStanFit'
extract(
  object,
  pars = NULL,
  permuted = FALSE,
  inc_warmup = FALSE,
  include = TRUE,
  regexp = FALSE,
  return_array = FALSE,
  ...
)

extract(object, ...)
```

## Arguments

- object:

  The jsdmStanFit model object

- pars:

  A character vector of parameter names

- permuted:

  Whether the draws should be randomly permuted, by default `FALSE`

- inc_warmup:

  Whether the warmup period should be included, by default `FALSE`

- include:

  Whether the parameters specified by `pars` should be included or
  excluded from the result, by default `TRUE` for inclusion

- regexp:

  Whether regular expression matching should be used to match the
  contents of `pars` to the parameter names, by default `FALSE`

- return_array:

  Whether to return the output as a 3 dimensional array (`TRUE`) or a
  named list (`FALSE`, the default)

- ...:

  Arguments passed to
  [`get_parnames()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdmstan-extractors.md)

## Value

If `return_array = FALSE` returns a named list with each parameter group
being an element of the list. Each list element is an array with the
first dimension being the iteration (all chains are appended) and the
other dimensions coming from the parameter dimensions. If
`return_array = TRUE` then a 3 dimensional array is returned with the
first dimension being the iterations, the second the chains and the
third the parameters.

## Functions

- `extract()`: Generic method
