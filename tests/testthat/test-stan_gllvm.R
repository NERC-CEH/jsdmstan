test_that("stan_gllvm fails with wrong inputs", {
  test_Y <- matrix(sample(0:1,size=100,replace=TRUE),
                                     nrow=20)
  expect_error(stan_gllvm(Y = test_Y, family = "bern"),
               "Must have at least one latent variable")

  expect_error(stan_gllvm(Y = matrix(rnorm(100),
                                     nrow=20),
                          D = 2, family = "bern"),
               "Y matrix is not binary")
  expect_error(stan_gllvm(Y = matrix(rnorm(100),
                                     nrow=20),
                          D = 2, family = "pois"),
               "Y matrix is not composed of integers")
})

test_that("stan_gllvm returns right type of object", {
  set.seed(04012022)
  gllvm_data <- gllvm_sim_data(N = 500, S = 16, D = 2, K = 5,
                               site_intercept = TRUE, family = "bern")
  gllvm_fit <- stan_gllvm(dat_list = gllvm_data, iter = 500, refresh = 0,
                          site_intercept = TRUE, family = "bern")

  expect_s3_class(gllvm_fit, "jsdmStanFit")

  # gaussian
  gllvm_data <- gllvm_sim_data(N = 100, S = 9, D = 2, K = 2, family = "gauss")
  gllvm_fit <- stan_gllvm(dat_list = gllvm_data, iter = 500, refresh = 0,
                          family = "gauss")

  expect_s3_class(gllvm_fit, "jsdmStanFit")

  # poisson
  gllvm_data <- gllvm_sim_data(N = 200, S = 12, D = 3, K = 3, family = "pois")
  gllvm_fit <- stan_gllvm(dat_list = gllvm_data, iter = 500, refresh = 0,
                          family = "pois")

  expect_s3_class(gllvm_fit, "jsdmStanFit")

  # neg binomial
  gllvm_data <- gllvm_sim_data(N = 100, S = 8, D = 2, family = "neg_bin")
  gllvm_fit <- stan_gllvm(dat_list = gllvm_data, iter = 500, refresh = 0,
                          family = "neg_b")

  expect_s3_class(gllvm_fit, "jsdmStanFit")

})
