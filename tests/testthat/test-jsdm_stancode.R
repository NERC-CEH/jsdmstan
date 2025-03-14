# jsdm_stancode checks
test_that("jsdm_stancode errors appropriately", {
  expect_error(jsdm_stancode(method = "mglmm", family = "nothing"))

  expect_error(jsdm_stancode(method = "mglmm", family = "poisson",
                             prior = list("betas" = "normal(0,1)")),
               "Prior must be given as a jsdmprior object")
})

jsdm_code <- jsdm_stancode(method = "mglmm", family = "zi_neg_binomial",
                           site_intercept = "grouped",
                           zi_param = "covariate")

test_that("jsdm_stancode print works", {
  expect_output(print(jsdm_code))
})


# jsdmStanFamily checks
test_that("jsdmStanFamily print works", {
  expect_output(print(jsdmStanFamily_empty()))
})
