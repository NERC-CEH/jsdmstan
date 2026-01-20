# Package index

## Simulate data

Simulate data according to the jsdm model structures

- [`jsdm_sim_data()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_sim_data.md)
  [`gllvm_sim_data()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_sim_data.md)
  [`mglmm_sim_data()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_sim_data.md)
  : Generate simulated data within a variety of jSDM methodologies
- [`rgampois()`](https://nerc-ceh.github.io/jsdmstan/reference/sim_helpers.md)
  [`rlkj()`](https://nerc-ceh.github.io/jsdmstan/reference/sim_helpers.md)
  [`rinvgamma()`](https://nerc-ceh.github.io/jsdmstan/reference/sim_helpers.md)
  [`rstudentt()`](https://nerc-ceh.github.io/jsdmstan/reference/sim_helpers.md)
  : Helper functions for simulating data

## Model

Run jsdmstan models

- [`stan_jsdm()`](https://nerc-ceh.github.io/jsdmstan/reference/stan_jsdm.md)
  : Fit jsdm models in Stan

- [`stan_gllvm()`](https://nerc-ceh.github.io/jsdmstan/reference/stan_gllvm.md)
  :

  Alias for `stan_jsdm` with `method = "gllvm"`

- [`stan_mglmm()`](https://nerc-ceh.github.io/jsdmstan/reference/stan_mglmm.md)
  :

  Alias for `stan_jsdm` with `method = "mglmm"`

- [`update(`*`<jsdmStanFit>`*`)`](https://nerc-ceh.github.io/jsdmstan/reference/update.jsdmStanFit.md)
  : Update a jsdmStanFit model object with new data or Stan arguments

## Set-up functions

Functions used to set up the prior structure (used in both data
simulation and model fitting) and the Stan code used within model
fitting

- [`jsdm_stancode()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_stancode.md)
  [`print(`*`<jsdmstan_model>`*`)`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_stancode.md)
  : Make stancode for the jsdm model
- [`jsdm_prior()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_prior.md)
  [`print(`*`<jsdmprior>`*`)`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_prior.md)
  : Create prior object for jsdmstan models and data simulation

## Summaries

Extract useful quantities from the model objects. This includes overall
summaries of model parameters and model diagnostics.

- [`extract()`](https://nerc-ceh.github.io/jsdmstan/reference/extract.jsdmStanFit.md)
  : Extract samples from jsdmStanFit object

- [`loo(`*`<jsdmStanFit>`*`)`](https://nerc-ceh.github.io/jsdmstan/reference/loo.jsdmStanFit.md)
  :

  Efficient approximate leave-one-out cross-validation using the loo
  package

- [`get_parnames()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdmstan-extractors.md)
  [`nuts_params(`*`<jsdmStanFit>`*`)`](https://nerc-ceh.github.io/jsdmstan/reference/jsdmstan-extractors.md)
  [`log_posterior(`*`<jsdmStanFit>`*`)`](https://nerc-ceh.github.io/jsdmstan/reference/jsdmstan-extractors.md)
  [`rhat(`*`<jsdmStanFit>`*`)`](https://nerc-ceh.github.io/jsdmstan/reference/jsdmstan-extractors.md)
  [`neff_ratio(`*`<jsdmStanFit>`*`)`](https://nerc-ceh.github.io/jsdmstan/reference/jsdmstan-extractors.md)
  : Extract quantities useful for model summaries

- [`print(`*`<jsdmStanFit>`*`)`](https://nerc-ceh.github.io/jsdmstan/reference/print.jsdmStanFit.md)
  : Print the default summary for the model

- [`summary(`*`<jsdmStanFit>`*`)`](https://nerc-ceh.github.io/jsdmstan/reference/summary.jsdmStanFit.md)
  : Summarise the model fit and data structure and give summaries for
  the parameters

- [`print(`*`<jsdmStanFamily>`*`)`](https://nerc-ceh.github.io/jsdmstan/reference/print.jsdmStanFamily.md)
  : Print jsdmStanFamily object

## Graphical summaries

Plots of various model parameters, including jSDM specific summary plots
as well as more general model summaries

- [`corrplot()`](https://nerc-ceh.github.io/jsdmstan/reference/corrplot.md)
  : Plot modelled correlations between species

- [`envplot()`](https://nerc-ceh.github.io/jsdmstan/reference/envplot.md)
  : Plotting environmental effects on species

- [`ordiplot()`](https://nerc-ceh.github.io/jsdmstan/reference/ordiplot.md)
  : Plotting an ordination plot for a GLLVM model

- [`mcmc_plot()`](https://nerc-ceh.github.io/jsdmstan/reference/mcmc_plot.jsdmStanFit.md)
  :

  MCMC plots implemented in bayesplot

- [`plot(`*`<jsdmStanFit>`*`)`](https://nerc-ceh.github.io/jsdmstan/reference/plot.jsdmStanFit.md)
  : Plot the traceplots and density plots for parameters within a
  jsdmStanFit object

## Posterior prediction

Extract posterior distributions and predict from them, including
graphical posterior checks

- [`posterior_linpred(`*`<jsdmStanFit>`*`)`](https://nerc-ceh.github.io/jsdmstan/reference/posterior_linpred.jsdmStanFit.md)
  : Access the posterior distribution of the linear predictor

- [`posterior_predict(`*`<jsdmStanFit>`*`)`](https://nerc-ceh.github.io/jsdmstan/reference/posterior_predict.jsdmStanFit.md)
  : Draw from the posterior predictive distribution

- [`posterior_shppred()`](https://nerc-ceh.github.io/jsdmstan/reference/posterior_shppred.md)
  : Access the posterior distribution of the linear predictor for the
  shape parameter

- [`posterior_zipred()`](https://nerc-ceh.github.io/jsdmstan/reference/posterior_zipred.md)
  : Access the posterior distribution of the linear predictor for
  zero-inflation parameter

- [`pp_check(`*`<jsdmStanFit>`*`)`](https://nerc-ceh.github.io/jsdmstan/reference/pp_check.jsdmStanFit.md)
  :

  Posterior predictive checks for `jsdmStanFit` objects

- [`multi_pp_check()`](https://nerc-ceh.github.io/jsdmstan/reference/multi_pp_check.md)
  : Multiple pp_check plots per species

- [`jsdm_statsummary()`](https://nerc-ceh.github.io/jsdmstan/reference/jsdm_statsummary.md)
  :

  Extract site- or species-level summary statistics for predictions from
  a `jsdmStanFit` model

## Other

Objects and classes used within jsdmstan

- [`jsdmStanFit`](https://nerc-ceh.github.io/jsdmstan/reference/jsdmStanFit.md)
  : jsdmStanFit class
- [`jsdmStanFamily`](https://nerc-ceh.github.io/jsdmstan/reference/jsdmStanFamily.md)
  : jsdmStanFamily class
- [`jsdmstan-package`](https://nerc-ceh.github.io/jsdmstan/reference/jsdmstan-package.md)
  [`jsdmstan`](https://nerc-ceh.github.io/jsdmstan/reference/jsdmstan-package.md)
  : The 'jsdmstan' package.
- [`bunce71`](https://nerc-ceh.github.io/jsdmstan/reference/bunce71.md)
  : Data from broadleaved woodlands across Great Britain, 1971
- [`pinewood`](https://nerc-ceh.github.io/jsdmstan/reference/pinewood.md)
  : Data from a survey of native pinewoods in Scotland, 1971
