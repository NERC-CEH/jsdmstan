set.seed(187750)
mglmm_data <- mglmm_sim_data(N = 30, S = 8, family = "gaussian", K = 3)
df <- as.data.frame(mglmm_data$X)

suppressWarnings(mglmm_fit <- stan_jsdm(~ V1 + V2 + V3,
                                        data = df, Y = mglmm_data$Y,
                                        family = "gaussian", method = "mglmm",
                                        refresh = 0, chains = 2, iter = 200
))
test_that("update works", {
  mglmm_data <- mglmm_sim_data(N = 20, S = 5, family = "gaussian", K = 2)
  suppressWarnings(mglmm_fit2 <- update(mglmm_fit,
                                        newY = mglmm_data$Y,
                                        newX = mglmm_data$X,
                                        refresh = 0,
                                        chains = 2, iter = 200
  ))
  suppressWarnings(mglmm_fit3 <- update(mglmm_fit,
                                        refresh = 0, iter = 100
  ))

  expect_s3_class(mglmm_fit2, "jsdmStanFit")
  expect_s3_class(mglmm_fit3, "jsdmStanFit")

  jsdm_empty <- jsdmStanFit_empty()
  expect_error(
    update(jsdm_empty),
    "Update requires the original data to be saved in the model object"
  )
})

gllvm_data <- gllvm_sim_data(N = 100, S = 8, family = "bern", D = 3,
                             site_intercept = "ungrouped")
gllvm_data$grps <- rep(1:20, each = 5)
gllvm_data$ngrp <- 20

test_that("site_intercept models update", {
  suppressWarnings(gllvm_fit <- stan_gllvm(X = NULL, dat_list = gllvm_data,
                                           site_intercept = "grouped",
                                           family = "bern",
                                           refresh = 0, chains = 1, iter = 200
  ))
  expect_s3_class(gllvm_fit, "jsdmStanFit")

  suppressWarnings(gllvm_fit2 <- update(gllvm_fit, newD = 2,
                                        refresh = 0, chains = 1, iter = 200
  ))
  expect_s3_class(gllvm_fit2, "jsdmStanFit")
})

gllvm_data <- gllvm_sim_data(N = 100, S = 8, family = "binomial", D = 3,
                             site_intercept = "ungrouped", Ntrials = 20)

test_that("binomial models update", {
  suppressWarnings(gllvm_fit <- stan_gllvm(dat_list = gllvm_data,
                                           family = "binomial",
                                           refresh = 0, chains = 1, iter = 200
  ))
  expect_s3_class(gllvm_fit, "jsdmStanFit")

  suppressWarnings(gllvm_fit2 <- update(gllvm_fit, newD = 2, newNtrials = 30,
                                        refresh = 0, chains = 1, iter = 200
  ))
  expect_s3_class(gllvm_fit2, "jsdmStanFit")
})

zip_data <- gllvm_sim_data(N = 97, S = 9,D = 2, family = "zi_poisson",K=1,
                           zi_param = "covariate")

test_that("zi models update", {
  suppressWarnings(gllvm_fit <- stan_gllvm(dat_list = zip_data,
                                           family = "zi_poisson",
                                           refresh = 0, chains = 1, iter = 200
  ))
  expect_s3_class(gllvm_fit, "jsdmStanFit")

  expect_error(gllvm_fit2 <- update(gllvm_fit, newD = 3,
                                        newZi_X = matrix(1:74, nrow = 74),
                                        refresh = 0, chains = 1, iter = 200
  ), "Number of rows of zi_X")


  expect_message(suppressWarnings(
    gllvm_fit2 <- update(gllvm_fit, newD = 3,
                         newZi_X = matrix(rnorm(194), nrow = 97),
                         refresh = 0, chains = 1, iter = 200
  )), "No column names")
  expect_equal(ncol(gllvm_fit2$data_list$zi_X),3)
  expect_equal(ncol(gllvm_fit3$data_list$X),2)
  expect_s3_class(gllvm_fit2, "jsdmStanFit")

  expect_message(suppressWarnings(
    gllvm_fit3 <- update(gllvm_fit, newD = 1,
                         newX = matrix(rnorm(194),nrow=97),
                         refresh = 0, chains = 1, iter = 200
    )), "No column names")
  expect_equal(ncol(gllvm_fit3$data_list$zi_X),3)
  expect_equal(ncol(gllvm_fit3$data_list$X),3)
  expect_s3_class(gllvm_fit3, "jsdmStanFit")
})


zinb_data <- gllvm_sim_data(N = 97, S = 9,D = 2, family = "zi_neg_binomial",
                           shp_param = "covariate", shp_k = 1)

test_that("shp models update", {
  suppressWarnings(gllvm_fit <- stan_gllvm(dat_list = zinb_data,
                                           family = "zi_neg_binomial",
                                           refresh = 0, chains = 1, iter = 200
  ))
  expect_s3_class(gllvm_fit, "jsdmStanFit")

  expect_error(gllvm_fit2 <- update(gllvm_fit, newD = 3,
                                    newShp_X = matrix(1:74, nrow = 74),
                                    refresh = 0, chains = 1, iter = 200
  ), "Number of rows of shp_X")


  expect_message(suppressWarnings(
    gllvm_fit2 <- update(gllvm_fit, newD = 3,
                         newShp_X = matrix(rnorm(194), nrow = 97),
                         refresh = 0, chains = 1, iter = 200
    )), "No column names")
  expect_equal(ncol(gllvm_fit2$data_list$shp_X),3)
  expect_s3_class(gllvm_fit2, "jsdmStanFit")
})
