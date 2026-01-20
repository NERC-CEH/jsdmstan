# Helper functions for simulating data

The `rlkj` function is for generating random LKJ correlation matrices
and the `rgampois` function generates random draws from the Stan's
alternative parameterisation of the negative binomial distribution.

## Usage

``` r
rgampois(n, mu, scale)

rlkj(n, eta = 1, cholesky = FALSE)

rinvgamma(n, shape, scale)

rstudentt(n, df, mu, sigma)
```

## Arguments

- n:

  The number of samples to create/dimension of correlation matrix

- mu:

  The mean used within the negative binomial parameterisation and the
  Student T distribution

- scale:

  The phi parameter that controls overdispersion of the negative
  binomial distribution (see details for description), or the scale
  parameter used within the inverse gamma distribution (see
  [`stats::rgamma()`](https://rdrr.io/r/stats/GammaDist.html))

- eta:

  The shape parameter of the LKJ distribution

- cholesky:

  Whether the correlation matrix should be returned as the Cholesky
  decomposition, by default `FALSE`

- shape:

  The shape parameter of the inverse gamma distribution (see
  [`stats::rgamma()`](https://rdrr.io/r/stats/GammaDist.html))

- df:

  The degrees of freedom parameter within the Student T distribution
  (see details)

- sigma:

  The scale of the Student T distribution (see details)

## Details

The Lewandowski-Kurowicka-Joe (LKJ) distribution is a prior distribution
for correlation matrices, with the shape parameter eta. If eta is 1 then
the density is uniform over the correlation matrix, ith eta \> 1 then
the the probability concentrates around the identity matrix while is 0
\< eta \< 1 the probability concentrates away from the identity matrix.

The alternative parameterisation of the negative binomial distribution
is:

\$\$NegBinomial2(y \| mu, scale) = binom(y+scale-1,y) (mu/mu+scale)^y
(scale/mu + scale)^scale\$\$

Where the mean of the distribution is mu and the variance is \\mu +
(mu^2/scale)\\

The `rlkj` function are sourced from Ben Goodrich's response on the Stan
google mailing list. (see link
<https://groups.google.com/g/stan-users/c/3gDvAs_qwN8/m/Xpgi2rPlx68J)>).
The `rgampois` function is sourced from the rethinking package by
Richard McElreath.

The alternative parameterisation of the Student T distribution is by
mean (mu) and scale (sigma) to be consistent with the Stan
parameterisation rather than the parameterisation in
[`stats::rt()`](https://rdrr.io/r/stats/TDist.html).
