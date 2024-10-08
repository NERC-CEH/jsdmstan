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
    "N and S must be positive integers"
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
    "N and S must be positive integers"
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
  gllvm_sim <- jsdm_sim_data(100,12,D=2,family = "binomial", method = "gllvm",
                             Ntrials = 19)
  expect_named(gllvm_sim, c(
    "Y", "pars", "N", "S", "D", "K", "X", "Ntrials"
  ))
  expect_length(gllvm_sim$Ntrials, 100)
  gllvm_sim <- jsdm_sim_data(100,12,D=2,family = "zi_neg_binomial", method = "gllvm",
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
