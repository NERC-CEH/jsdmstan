test_that("is whole number works", {
  expect_false(is.wholenumber(5.5))
  expect_true(is.wholenumber(1e4))
  expect_true(is.wholenumber(9L))
  expect_false(is.wholenumber("sdodhg"))
})

test_fit <- jsdmStanFit_empty()

test_that("plot errors when expected", {
  expect_error(plot(test_fit, N = 3.5),
               "N must be a positive integer")
  expect_error(plot(test_fit, pars = c("fiona","is","bored")),
               "Please specify pars within the model, use get_parnames to find the names")
  expect_error(plot(test_fit, sample_n = 5.6),
               "If pars is NULL then sample_n must be a positive integer")
})

