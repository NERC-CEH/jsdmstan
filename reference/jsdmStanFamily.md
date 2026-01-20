# jsdmStanFamily class

This is the jsdmStanFamily class, which occupies a slot within any
jsdmStanFit object.

## Elements for `jsdmStanFamily` objects

- `family`:

  A length one character vector describing family used to fit object.
  Options are `"gaussian"`, `"poisson"`, `"bernoulli"`,
  `"neg_binomial"`, `"binomial"`, `"zi_poisson"`, `"zi_neg_binomial"`,
  or `"multiple"`.

- `params`:

  A character vector that includes all the names of the family-specific
  parameters.

- `params_dataresp`:

  A character vector that includes any named family-specific parameters
  that are modelled in response to data.

- `preds`:

  A character vector of the measured predictors included if family
  parameters are modelled in response to data. If family parameters are
  not modelled in response to data this is left empty.

- `data_list`:

  A list containing the original data used to fit the model (empty when
  save_data is set to `FALSE` or family parameters are not modelled in
  response to data).
