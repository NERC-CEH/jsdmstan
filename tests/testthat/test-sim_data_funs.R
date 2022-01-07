test_that("gllvm_sim_data errors with bad inputs", {
  expect_error(gllvm_sim_data(N = 100, S = 10, D = 2,
                              family = "bern", species_intercept = FALSE),
               "If K is 0 then a species intercept is required")

  expect_error(gllvm_sim_data(N = "a", S = 10, D = 2, family = "pois"),
               "N and S must be positive integers")

  expect_error(gllvm_sim_data(N = 100, S = 8, D = 2, phylo  = TRUE,
                              family = "neg_bin"))
})

test_that("gllvm_sim_data returns a list of correct length", {
  gllvm_sim <- gllvm_sim_data(N = 100, S = 8, D = 2, family = "bern")
  expect_named(gllvm_sim, c("Y", "pars", "N","S","D","K","X",
                            "site_intercept"))
})

test_that("mglmm_sim_data errors with bad inputs", {
  expect_error(mglmm_sim_data(N = 100, S = 10, species_intercept = FALSE,
                              family = "gauss"),
               "If K is 0 then a species intercept is required")

  expect_error(mglmm_sim_data(N = "a", S = 10, family = "neg_bin"),
               "N and S must be positive integers")
})

test_that("mglmm_sim_data errors without phylo details", {
  expect_error(mglmm_sim_data(N = 100, S = 5, family = "bern", phylo = TRUE),
               "Need to specify delta and nu05 arguments for phylo")

  expect_error(mglmm_sim_data(N = 100, S = 5, family = "gauss", phylo = TRUE,
                              nu05 = 4, delta = 1e-10),
               "nu05 must be integer in range 0-3, given as:4")
})

test_that("mglmm_sim_data returns a list of correct length", {
  mglmm_sim <- mglmm_sim_data(N = 100, S = 8, family = "bern")
  expect_named(mglmm_sim, c("Y", "pars", "N","S","D","K","X",
                            "site_intercept"))
  mglmm_sim_phylo <- mglmm_sim_data(N = 100, S = 5, phylo = TRUE, delta = 1e-5,
                                    nu05= 1, family = "gauss")
  expect_named(mglmm_sim_phylo, c("Y", "pars", "N","S","D","K","X",
                                  "site_intercept","Dmat","delta","nu05"))
})
