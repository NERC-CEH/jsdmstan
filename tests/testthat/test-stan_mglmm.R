test_that("stan_mglmm fails with wrong inputs", {
  expect_error(stan_mglmm(Y = matrix(rnorm(100),
                                     nrow=20),
                          family = "bern"),
               "Y matrix is not binary")
  expect_error(stan_mglmm(Y = matrix(rnorm(100),
                                     nrow=20),
                          family = "pois"),
               "Y matrix is not composed of integers")
})

test_that("stan_mglmm returns right type of object", {
  set.seed(04012022)
  mglmm_data <- mglmm_sim_data(N = 100, S = 8, family = "bern",
                               site_intercept = TRUE)
  mglmm_fit <- stan_mglmm(Y = mglmm_data$Y, X = mglmm_data$X,
                          family = "bern",
                          refresh = 0, site_intercept = TRUE)

  expect_s3_class(mglmm_fit, "jsdmStanFit")

  # gaussian - now in test-stan_jsdm
  # mglmm_data <- mglmm_sim_data(N = 100, S = 5, family = "gauss")
  # mglmm_fit <- stan_mglmm(dat_list = mglmm_data, family = "gauss",
  #                         refresh = 0)
  #
  # expect_s3_class(mglmm_fit, "jsdmStanFit")

  # poisson
  mglmm_data <- jsdm_sim_data(method = "mglmm", N = 200, S = 8, K = 3,
                              family = "poisson")
  mglmm_fit <- stan_jsdm(X = NULL, dat_list = mglmm_data, family = "poisson",
                         refresh = 0, method = "mglmm",
                         control = list(adapt_delta = 0.95), iter = 4000)

  expect_s3_class(mglmm_fit, "jsdmStanFit")

  # neg bin
  mglmm_data <- mglmm_sim_data(N = 127, S = 12, K = 2, family = "neg_bin")
  mglmm_fit <- stan_mglmm(Y = mglmm_data$Y, X = mglmm_data$X, family = "neg_bin",
                          refresh = 0)

  expect_s3_class(mglmm_fit, "jsdmStanFit")
})
