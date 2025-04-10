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


test_that("envplot errors when expected", {
  expect_error(
    envplot(list()),
    "Only objects of class jsdmStanFit are supported"
  )
  expect_error(
    envplot(test_fit),
    "This beta parameterisation currently unsupported"
  )
  expect_error(
    envplot(test_fit, species = list()),
    "character or integer vector"
  )
  expect_error(
    envplot(test_fit, preds = matrix()),
    "character or integer vector"
  )
  expect_error(
    envplot(test_fit, species = c(1,5,8)),
    "corresponds to the species indices"
  )
  expect_error(
    envplot(test_fit, preds = c(1.2,0,8)),
    "corresponds to the predictor indices"
  )
  expect_error(
    envplot(test_fit, species = c("eight","bottles")),
    "must only include"
  )
  expect_error(
    envplot(test_fit, preds = c("eight","bottles")),
    "must only include"
  )
})

test_that("corrplot errors when expected", {
  expect_error(
    corrplot(list()),
    "Only objects of class jsdmStanFit with method mglmm are supported"
  )
  test_fit2 <- test_fit
  test_fit2$jsdm_type <- "gllvm"
  expect_error(
    corrplot(test_fit2),
    "Only objects of class jsdmStanFit with method mglmm are supported"
  )

  test_fit2$jsdm_type <- "mglmm"
  expect_error(
    corrplot(test_fit2, species = -1),
    "Species must be either a"
  )

  expect_error(
    corrplot(test_fit2, plotfun = "argh"),
    "not a valid ppc type"
  )

  test_fit2$species <- LETTERS[1:6]
  expect_error(
    corrplot(test_fit2, species = LETTERS[9:11]),
    "Species specified are not found in the model fit object"
  )

})
