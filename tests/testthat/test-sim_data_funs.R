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

  expect_error(gllvm_sim_data(
    N = 100, S = 8, D = 2, phylo = TRUE,
    family = "neg_bin"
  ))
})

test_that("gllvm_sim_data returns a list of correct length", {
  gllvm_sim <- gllvm_sim_data(N = 100, S = 8, D = 2, family = "bern")
  expect_named(gllvm_sim, c(
    "Y", "pars", "N", "S", "D", "K", "X",
    "site_intercept"
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
})


test_that("mglmm_sim_data returns a list of correct length", {
  mglmm_sim <- mglmm_sim_data(N = 100, S = 8, family = "bern")
  expect_named(mglmm_sim, c(
    "Y", "pars", "N", "S", "D", "K", "X",
    "site_intercept"
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
    "Y", "pars", "N", "S", "D", "K", "X",
    "site_intercept"
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
