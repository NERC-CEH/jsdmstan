test_that("is whole number works", {
  expect_false(is.wholenumber(5.5))
  expect_true(is.wholenumber(1e4))
  expect_true(is.wholenumber(9L))
  expect_false(is.wholenumber("sdodhg"))
})

test_fit <- jsdmStanFit_empty()

test_that("plot errors when expected", {
  expect_error(
    plot(test_fit, N = 3.5),
    "N must be a positive integer"
  )
  expect_error(
    plot(test_fit, pars = c("fiona", "is", "bored")),
    "Please specify pars within the model, use get_parnames to find the names"
  )
  expect_error(
    plot(test_fit, sample_n = 5.6),
    "If pars is NULL then sample_n must be a positive integer"
  )
})

test_that("mcmc_plot errors when expected", {
  expect_error(
    mcmc_plot(test_fit, plotfun = "fiona"),
    "Invalid plotfun argument"
  )
  expect_error(
    mcmc_plot(test_fit, pars = NULL, sample_n = 8.6),
    "If pars is NULL then sample_n must be a positive integer"
  )
})

test_that("ordiplot errors when expected", {
  expect_error(
    ordiplot(list()),
    "Only objects of class jsdmStanFit are supported"
  )
  expect_error(
    ordiplot(test_fit),
    "Only gllvm models are supported"
  )
  test_fit$jsdm_type <- "gllvm"
  expect_error(
    ordiplot(test_fit, choices = 1:3),
    "Only two latent variables can be plotted at once"
  )
})
