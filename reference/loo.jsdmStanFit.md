# Efficient approximate leave-one-out cross-validation using the loo package

This function uses the loo package to compute PSIS-LOO CV, efficient
approximate leave-one-out (LOO) cross-validation for Bayesian models
using Pareto smoothed importance sampling (PSIS). This requires that the
model was fit using `log_lik = TRUE`.

## Usage

``` r
# S3 method for class 'jsdmStanFit'
loo(x, calc_reff = TRUE, cores = getOption("mc.cores", 1), ...)
```

## Arguments

- x:

  The jsdmStanFit model object

- calc_reff:

  Whether to calculate the relative efficiencies for loo, by default
  `TRUE`. If set to `FALSE` then relative efficiency is assumed to be 1.

- cores:

  The number of cores the loo functions use, by default uses the
  mc.cores option (or 1, if unspecified).

- ...:

  Other arguments passed to the
  [`loo`](https://mc-stan.org/loo/reference/loo.html) function

## Value

A list with class `c("psis_loo","loo")`, as detailed in the
[`loo`](https://mc-stan.org/loo/reference/loo.html) documentation
