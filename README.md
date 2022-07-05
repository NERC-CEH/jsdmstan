# jsdmstan: joint Species Distribution Models in Stan

<!-- badges: start -->
[![R-CMD-check](https://github.com/NERC-CEH/jsdmstan/workflows/R-CMD-check/badge.svg)](https://github.com/NERC-CEH/jsdmstan/actions)
[![Codecov test coverage](https://codecov.io/gh/NERC-CEH/jsdmstan/branch/main/graph/badge.svg)](https://codecov.io/gh/NERC-CEH/jsdmstan?branch=main)
<!-- badges: end -->

This is an R package for running joint Species Distribution Models (jSDM) in [Stan](https://mc-stan.org/). jSDMs are models where multiple response variables (i.e. species) are fit at the same time, and the covariance between these species are used to inform the model results. For a review of jSDMs see Warton et al. (2015) So many variables: joint modelling in community ecology. *TREE*, 30:766-779 DOI: [10.1016/j.tree.2015.09.007](http://doi.org/10.1016/j.tree.2015.09.007).

This package can fit data to a Multivariate Generalised Linear Mixed Model (MGLMM) or a Generalised Linear Latent Variable Model (GLLVM), and also provides functionality for simulating data under these scenarios and an interface to the [bayesplot](https://mc-stan.org/bayesplot/) package for a wide variety of plotting options.

Example code:

```
library(jsdmstan)
```

Simulate data:
```
nsites <- 200
nspecies <- 9
ncovar <- 2
mglmm_test_data <- mglmm_sim_data(N = nsites, S = nspecies, 
                                  K = ncovar, family = "pois")
```

Fit Stan model:
```
dat <- as.data.frame(mglmm_test_data$X)
mglmm_fit <- stan_jsdm(~ V1 + V2, data = dat, Y = mglmm_test_data$Y, 
                       family = "pois", method = "mglmm")
```

Plot results:
```
plot(mglmm_fit)
mcmc_plot(mglmm_fit, plotfun = "rhat_hist")
```

***

This work was funded by the Natural Environment Research Council (part of UK Research and Innovation) under the UK-SCAPE Programme delivering National Capability (Grant Reference NE/R016429/1).

![UKCEH logo](./man/figures/UKCEH-Logo.png) ![UKSCAPE logo](./man/figures/uk_scape_logo.png)
