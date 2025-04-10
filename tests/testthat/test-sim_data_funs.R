set.seed(23359)
test_that("gllvm_sim_data errors with bad inputs", {
  expect_error(
    gllvm_sim_data(
      N = 100, S = 10, D = 2,
      family = "bern", species_intercept = FALSE
    ),
    "If K is 0 then a species intercept is required"
  )

  expect_error(
    gllvm_sim_data(N = "a", S = 10, D = 2, family = "pois"),
    "N must be a positive integer"
  )

  expect_error(
    gllvm_sim_data(N = 200, S = 5, D = 3, family = "neg_bin",
                   site_intercept = "grouped"),
    "Grouped site intercept not supported"
  )

  expect_error(
    gllvm_sim_data(N = 200, S = 8, D = 2, family = "binomial", Ntrials = "1"),
    "Ntrials must be a positive integer"
  )

  expect_error(
    gllvm_sim_data(N = 200, S = 8, D = 2, family = "zi_poisson",
                   zi_param = "covariate", zi_k = -2),
    "zi_k must be either NULL or a positive integer"
  )

})

test_that("gllvm_sim_data returns a list of correct length", {
  gllvm_sim <- gllvm_sim_data(N = 100, S = 8, D = 2, family = "bern",
                              beta_param = "unstruct")
  expect_named(gllvm_sim, c(
    "Y", "pars", "N", "S", "D", "K", "X"
  ))
  gllvm_sim <- gllvm_sim_data(N = 100, S = 8, D = 2, family = "bern",
                              site_intercept = "ungrouped")
  expect_named(gllvm_sim, c(
    "Y", "pars", "N", "S", "D", "K", "X"
  ))
})

test_that("mglmm_sim_data errors with bad inputs", {
  expect_error(
    mglmm_sim_data(
      N = 100, S = 10, species_intercept = FALSE,
      family = "gauss"
    ),
    "If K is 0 then a species intercept is required"
  )

  expect_error(
    mglmm_sim_data(N = "a", S = 10, family = "neg_bin"),
    "N must be a positive integer"
  )
  expect_error(
    mglmm_sim_data(N = 200, S = 5, family = "gaussian",
                   site_intercept = "grouped"),
    "Grouped site intercept not supported"
  )

  expect_error(
    mglmm_sim_data(N = 100, S = 5, family = "binomial"),
    "Number of trials must be specified"
  )

  expect_error(
    mglmm_sim_data(N = 50, S = 8, family = "binomial", Ntrials = c(1,3)),
    "Ntrials must be of length"
  )
})


test_that("mglmm_sim_data returns a list of correct length", {
  mglmm_sim <- mglmm_sim_data(N = 100, S = 8, family = "bern",
                              K = 3, species_intercept = FALSE)
  expect_named(mglmm_sim, c(
    "Y", "pars", "N", "S", "D", "K", "X"
  ))
  mglmm_sim <- mglmm_sim_data(N = 100, S = 8, family = "pois",
                              beta_param = "unstruct",
                              site_intercept = "ungrouped")
  expect_named(mglmm_sim, c(
    "Y", "pars", "N", "S", "D", "K", "X"
  ))
  gllvm_sim <- jsdm_sim_data(N = 100, S = 12,D = 2,
                             family = "binomial", method = "gllvm",
                             Ntrials = 19)
  expect_named(gllvm_sim, c(
    "Y", "pars", "N", "S", "D", "K", "X", "Ntrials"
  ))
  expect_length(gllvm_sim$Ntrials, 100)
  gllvm_sim <- jsdm_sim_data(N = 100,S = 12,D=2,
                             family = "zi_neg_binomial", method = "gllvm",
                             zi_param = "covariate")
  expect_named(gllvm_sim, c(
    "Y", "pars", "N", "S", "D", "K", "X", "zi_k", "zi_X"
  ))
  expect_equal(dim(gllvm_sim$Y),c(100,12))
})

test_that("jsdm_sim_data returns all appropriate pars", {
  mglmm_sim <- jsdm_sim_data(100,12,family = "gaussian", method = "mglmm",
                             beta_param = "cor", K = 3)
  expect_named(mglmm_sim$pars, c(
    "betas","sigmas_preds","z_preds","cor_preds","sigmas_species",
    "cor_species","z_species","sigma"
  ))
  gllvm_sim <- jsdm_sim_data(100,12,D=2,family = "neg_bin", method = "gllvm",
                             beta_param = "unstruct",
                             site_intercept = "ungrouped")
  expect_named(gllvm_sim$pars, c(
    "betas","a_bar","sigma_a","a","L","LV","sigma_L","kappa"
  ))

  gllvm_sim2 <- jsdm_sim_data(100,12,D=2,family = "zi_poisson", method = "gllvm",
                              beta_param = "unstruct",
                              site_intercept = "ungrouped")
  expect_named(gllvm_sim2$pars, c(
    "betas","a_bar","sigma_a","a","L","LV","sigma_L","zi"
  ))


  gllvm_sim3 <- jsdm_sim_data(100,9,K=2,D=2,family = "zi_neg_bin", method = "gllvm",
                              beta_param = "unstruct",
                              site_intercept = "ungrouped", zi_param = "covariate",
                              zi_k = 1)
  expect_named(gllvm_sim3$pars, c(
    "betas","a_bar","sigma_a","a","L","LV","sigma_L","zi_betas","kappa"
  ))

})

test_that("prior specification works", {
  jsdm_sim <- jsdm_sim_data(
    N = 100, S = 8, K = 2, family = "gaus", method = "mglmm",
    prior = jsdm_prior(
      z_species = "student_t(3,0,2)",
      sigma = "gamma(1,1)"
    )
  )
  expect_named(jsdm_sim, c(
    "Y", "pars", "N", "S", "D", "K", "X"
  ))
  jsdm_sim <- jsdm_sim_data(
    N = 50, S = 8, D = 2,
    family = "neg_binomial", method = "gllvm",
    prior = jsdm_prior(
      sigma_L = "inv_gamma(10,0.1)",
      kappa = "cauchy(1,1)"
    )
  )

  expect_error(
    jsdm_sim_data(
      N = 50, S = 5, D = 2, family = "bern", method = "gllvm",
      prior = list(LV = "weibull(3,5)")
    ),
    "prior object must be of class jsdmprior"
  )

  expect_error(
    jsdm_sim_data(
      N = 50, S = 5, K = 2, D = 2,
      family = "gaus", method = "gllvm",
      prior = jsdm_prior(LV = "weibull(3,5)")
    ),
    "Not all prior distributions specified are supported."
  )
})

test_that("Can simulate data by providing X directly",{
  Xdat <- matrix(runif(100), ncol = 2)

  expect_error(jsdm_sim_data(S = 4, family = "gaussian",
                             method = "mglmm", X = "for"),
               "must be a matrix")

  expect_message(simdat <- jsdm_sim_data(S = 4, family = "gaussian",
                                         method = "mglmm", X = Xdat),
                 "no column names")
  expect_named(simdat, c("Y","pars","N","S","D","K","X"))
  expect_equal(Xdat, unname(simdat$X))

  Xdat2 <- matrix(runif(100), ncol = 1)

  expect_message(simdat2 <- jsdm_sim_data(S = 13, family = "neg_binomial",
                                         method = "gllvm", D = 2, X = Xdat2),
                 "no column names")
  expect_named(simdat2, c("Y","pars","N","S","D","K","X"))
  expect_equal(Xdat2, unname(simdat2$X))
})

test_that("Shape parameter as a function of data works",{
  expect_error(jsdm_sim_data(S = 9, N = 50, family = "gaussian",
                             method = "gllvm", D = 3, shp_param = "blank"),
               "should be one of")
  testdat <- jsdm_sim_data(S = 9, N = 50, K = 2, family = "gaussian",
                           method = "gllvm", D = 3, shp_param = "covariate")
  expect_named(testdat, c("Y","pars","N","S","D","K","X","shp_k", "shp_X"))
  expect_equal(testdat$X, testdat$shp_X[,-1])

  testdat <- jsdm_sim_data(S = 7, N = 74, K = 2, family = "neg_binomial",
                           method = "mglmm", shp_param = "covariate", shp_k = 3)
  expect_named(testdat, c("Y","pars","N","S","D","K","X","shp_k", "shp_X"))
  expect_false(identical(testdat$X, testdat$shp_X[,-1]))



})

test_that("Supplying matrix to shape parameter works",{
  expect_error(jsdm_sim_data(S = 7, N = 74, K = 2, family = "neg_binomial",
                             method = "mglmm", shp_param = "covariate",
                             shp_X = matrix(runif(10),ncol=1)),
               "shp_X must have N rows")
  shp_Xdat <- matrix(runif(81),ncol=1)
  colnames(shp_Xdat) <- "V1"
  expect_message(testdat <- jsdm_sim_data(S = 7, N = 81, K = 2,
                                          family = "neg_binomial",
                                          method = "mglmm", shp_param = "covariate",
                                          shp_X = shp_Xdat),
                 "without an intercept"
  )
  expect_named(testdat, c("Y","pars","N","S","D","K","X","shp_k", "shp_X"))
})


test_that("Supplying matrix to zi parameter works",{
  expect_error(jsdm_sim_data(S = 7, N = 74, K = 2, family = "zi_neg_binomial",
                             method = "mglmm", zi_param = "covariate",
                             zi_X = matrix(runif(10),ncol=1)),
               "zi_X must have N rows")
  zi_Xdat <- matrix(runif(81),ncol=1)
  colnames(zi_Xdat) <- "V1"
  expect_message(testdat <- jsdm_sim_data(S = 7, N = 81, K = 2,
                                          family = "zi_poisson",
                                          method = "gllvm",D = 1,
                                          zi_param = "covariate",
                                          zi_X = zi_Xdat),
                 "without an intercept"
  )
  expect_named(testdat, c("Y","pars","N","S","D","K","X","zi_k", "zi_X"))

  expect_message(testdat <- jsdm_sim_data(S = 7, N = 81, K = 2,
                                          family = "zi_neg_binomial",
                                          method = "mglmm", zi_param = "covariate",
                                          zi_X = zi_Xdat, shp_param = "covariate",
                                          shp_X = zi_Xdat),
                 "without an intercept"
  )
  expect_named(testdat, c("Y","pars","N","S","D","K","X","zi_k", "zi_X","shp_k", "shp_X"))
  expect_named(testdat$pars,c("betas", "sigmas_species","cor_species","z_species",
                              "shp_betas","zi_betas","kappa"))
})
