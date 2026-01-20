# Extract quantities useful for model summaries

These are methdos for extracting various useful summaries from models,
including the model parameter names, NUTS parameters, the log posterior,
r-hat and n-eff ratio.

## Usage

``` r
get_parnames(object, log_lik = FALSE)

# S3 method for class 'jsdmStanFit'
nuts_params(object, ...)

# S3 method for class 'jsdmStanFit'
log_posterior(object, ...)

# S3 method for class 'jsdmStanFit'
rhat(object, ...)

# S3 method for class 'jsdmStanFit'
neff_ratio(object, ...)
```

## Arguments

- object:

  The `jsdmStanFit` model object

- log_lik:

  Whether the log_lik parameters should be included in the parameter
  list

- ...:

  Arguments passed on to the bayesplot equivalent for stanFit objects

## Value

`get_parnames()` returns a character vector of model parameter
names.`nuts_params()` returns a molten data frame (see
[`reshape2::melt()`](https://rdrr.io/pkg/reshape2/man/melt.html)). The
data frame should have columns "Parameter" (factor), "Iteration"
(integer), "Chain" (integer), and "Value" (numeric).`log_posterior()`
returns a molten data frame (see
[`reshape2::melt()`](https://rdrr.io/pkg/reshape2/man/melt.html)). The
data frame should have columns "Chain" (integer), "Iteration" (integer),
and "Value" (numeric).`rhat()`, `neff_ratio()` both return named numeric
vectors.

## Functions

- `get_parnames()`: Get the model parameter names

- `nuts_params(jsdmStanFit)`: Get the NUTS parameters

- `log_posterior(jsdmStanFit)`: Get the log posterior

- `rhat(jsdmStanFit)`: Get the R-hat

- `neff_ratio(jsdmStanFit)`: Get the n_eff ratio
