# jsdmstan: joint Species Distribution Models in Stan

This is an R package for running joint Species Distribution Models (jSDM) in [Stan](https://mc-stan.org/). jSDMs are models where multiple response variables (i.e. species) are fit at the same time, and the covariance between these species are used to inform the model results. For a review of jSDMs see Warton et al. (2015) So many variables: joint modelling in community ecology. TREE, 30:766-779 DOI: [10.1016/j.tree.2015.09.007](http://doi.org/10.1016/j.tree.2015.09.007).

This package can fit data to a Multivariate Generalised Linear Mixed Model (MGLMM) or a Generalised Linear Latent Variable Model (GLLVM), and also provides functionality for simulating data under these scenarios and an interface to the [bayesplot](https://mc-stan.org/bayesplot/) package for a wide variety of plotting options.
