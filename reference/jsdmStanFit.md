# jsdmStanFit class

This is the jsdmStanFit class, which stan_gllvm and stan_mglmm both
create.

## Elements for `jsdmStanFit` objects

- `jsdm_type`:

  A length one character vector describing type of jSDM

- `family`:

  A jsdmStanFamily object describing characteristics of family

- `species`:

  A character vector of the species names

- `sites`:

  character vector of the site IDs

- `preds`:

  A character vector of the measured predictors included

- `data_list`:

  A list containing the original data used to fit the model (empty when
  save_data is set to `FALSE`)

- `n_latent`:

  A length one integer vector representing number of latent variables
  (in gllvm type fits) or `NA` in all other cases
