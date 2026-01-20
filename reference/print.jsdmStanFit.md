# Print the default summary for the model

This prints out a summary for the models which includes the type of
model fit, the number of species, sites and predictors as well as a
summary of any parameters with Rhat \> 1.01 or effective sample size to
total number of samples ratio \< 0.05

## Usage

``` r
# S3 method for class 'jsdmStanFit'
print(x, ...)
```

## Arguments

- x:

  The `jsdmStanFit` model object

- ...:

  Other arguments passed to
  [summary.jsdmStanFit](https://nerc-ceh.github.io/jsdmstan/reference/summary.jsdmStanFit.md)
