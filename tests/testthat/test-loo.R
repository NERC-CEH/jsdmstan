set.seed(548935)
mglmm_data <- mglmm_sim_data(N = 30, S = 8, family = "neg_binomial", K = 3)
df <- as.data.frame(mglmm_data$X)

suppressWarnings(mglmm_fit <- stan_jsdm(~ V1 + V2 + V3,
                                        data = df, Y = mglmm_data$Y,
                                        family = "neg_binomial", method = "mglmm",
                                        refresh = 0, chains = 2, iter = 200
))

test_that("loo works for mglmm", {
  expect_warning(
    mglmm_loo <- loo(mglmm_fit),
    "Pareto k diagnostic"
  )
  expect_s3_class(mglmm_loo, "psis_loo")
})


gllvm_data <- gllvm_sim_data(N = 63, S = 5, D = 2, family = "zi_poisson", K = 1,
                             zi_param = "covariate")
df <- as.data.frame(gllvm_data$X)

suppressWarnings(gllvm_fit <- stan_jsdm(~ V1,
                                        data = df, Y = gllvm_data$Y,
                                        family = "zi_poisson", method = "gllvm",
                                        zi_param = "covariate", D = 2,
                                        refresh = 0, chains = 2, iter = 200
))

test_that("loo works for gllvm", {
  expect_warning(
    gllvm_loo <- loo(gllvm_fit),
    "Pareto k diagnostic"
  )
  expect_s3_class(gllvm_loo, "psis_loo")
})

# test that density functions return numbers
test_that("dzi functions work", {
  x <- dzipois(sample(0:5,20,replace=TRUE),rnorm(20),rnorm(20), log = TRUE)
  expect_length(x, 20)
  expect_true(all(x<0))


  x <- dzinb(sample(0:5,23,replace=TRUE),mu = rnorm(23),
             kappa = abs(rnorm(23)), zi = rnorm(23), log = TRUE)
  expect_length(x, 23)
  expect_true(all(x<0))

})


set.seed(6742678)
Xdat <- as.matrix(data.frame(V1 = sample(0:2, 87, replace = TRUE),
                             V2 = sample(0:1, 87, replace = TRUE),
                             V3 = rnorm(87)))
df <- as.data.frame(Xdat)

test_that("loo works for mglmm zi_nb", {
  expect_message(
    mglmm_data <- mglmm_sim_data(S = 13, X = Xdat, zi_X = Xdat[,1,drop=FALSE],
                                 shp_X = Xdat[,2,drop=FALSE],
                                 family = "zi_neg_binomial"),
    "without an intercept")
  expect_named(mglmm_data,
               c("Y","pars","N","S","D","K","X","zi_k","zi_X","shp_k","shp_X"))
  expect_named(mglmm_data$pars,
               c("betas","sigmas_species","cor_species","z_species","shp_betas",
                 "zi_betas","kappa"))

  suppressWarnings(mglmm_fit <- stan_jsdm(~ V1 + V2 + V3,
                                          data = df, Y = mglmm_data$Y,
                                          zi_formula = ~ V1, shp_formula = ~V2,
                                          family = "zi_neg_binomial", method = "mglmm",
                                          refresh = 0, chains = 2, iter = 200
  ))

  expect_warning(
    mglmm_loo <- loo(mglmm_fit),
    "Pareto k diagnostic"
  )
  expect_s3_class(mglmm_loo, "psis_loo")
})






test_that("loo works for gllvm gaussian", {
  set.seed(37687168)
  Xdat <- as.matrix(data.frame(V1 = sample(0:2, 63, replace = TRUE),
                               V2 = sample(0:1, 63, replace = TRUE),
                               V3 = rnorm(63)))

  df <- as.data.frame(Xdat)

  expect_message(
    gauss_data <- gllvm_sim_data(S = 18, X = Xdat, D = 2,
                                 shp_X = Xdat[,2,drop=FALSE],
                                 family = "gaussian"),
    "without an intercept")
  expect_named(gauss_data,
               c("Y","pars","N","S","D","K","X","shp_k","shp_X"))
  expect_named(gauss_data$pars,
               c("betas","L","LV","sigma_L","sigma","shp_betas"))
  suppressWarnings(gauss_fit <- stan_jsdm(~ V1 + V2 + V3,
                                          data = df, Y = gauss_data$Y,
                                          shp_formula = ~V2, D = 2,
                                          family = "gaussian", method = "gllvm",
                                          refresh = 0, chains = 2, iter = 200
  ))
  expect_warning(
    gauss_loo <- loo(gauss_fit),
    "Pareto k diagnostic"
  )
  expect_s3_class(gauss_loo, "psis_loo")
})
