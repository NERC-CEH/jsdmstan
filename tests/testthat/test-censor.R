test_that("sim_data errors appropriately", {
  expect_error(jsdm_sim_data(100,10,D=2,family = "gamma", method = "gllvm",
                             censoring = "left"),
               "please supply the censoring points")

  expect_error(jsdm_sim_data(7,89,D=2,family = "lognormal", method = "gllvm",
                             censoring = "left", censor_points = 1),
               "S-length vector")

  expect_error(jsdm_sim_data(7,89,D=2,family = "lognormal", method = "gllvm",
                             censoring = "left", censor_points = "a"),
               "S-length vector")

  expect_error(jsdm_sim_data(7,89,D=2,family = "lognormal", method = "gllvm",
                             censoring = "right", censor_points = c(1,2)),
               "S-length vector")
})

set.seed(676981)
gamma_dat <- jsdm_sim_data(9,75,D=2,family = "gamma", method = "gllvm",
                           censoring = "left",
                           censor_points = rep(c(1,1.5,2),each=3))
test_that("gamma sim_data returns correct censored data",{
  expect_named(gamma_dat, c("Y","pars","N","S","D","K","X",
                            "Y_uncensored","censor_points","cens_ID",
                            "censoring"))

  expect_false(anyNA(gamma_dat$Y))

  expect_false(anyNA(gamma_dat$Y_uncensored))

  expect_true(all(gamma_dat$Y_uncensored>0))
  expect_true(all(gamma_dat$Y>0))

  expect_true(all(sapply(1:9, function(i) all(gamma_dat$Y[,i] >= gamma_dat$censor_points[i]))))

  expect_length(gamma_dat$censor_points, 9)
  expect_type(gamma_dat$censor_points, "double")

  expect_equal(gamma_dat$censoring, "left")
})

set.seed(733148)
lnorm_dat <- jsdm_sim_data(8,53,K = 3,family = "lognorm", method = "mglmm",
                           censoring = "left",
                           prior = jsdm_prior(betas = "normal(0,0.25)",
                                              sigma = "normal(0,0.5)"),
                           censor_points = rep(0.5,8))
test_that("lnorm sim_data returns correct censored data",{
  expect_named(lnorm_dat, c("Y","pars","N","S","D","K","X",
                            "Y_uncensored","censor_points","cens_ID",
                            "censoring"))

  expect_false(anyNA(lnorm_dat$Y))

  expect_false(anyNA(lnorm_dat$Y_uncensored))

  expect_true(all(lnorm_dat$Y_uncensored>0))
  expect_true(all(lnorm_dat$Y>0))

  expect_true(all(sapply(1:8, function(i) all(lnorm_dat$Y[,i] >= lnorm_dat$censor_points[i], na.rm = TRUE))))

  expect_length(lnorm_dat$censor_points, 8)
  expect_type(lnorm_dat$censor_points, "double")

  expect_equal(lnorm_dat$censoring, "left")
})

suppressWarnings(gamma_mod <- stan_jsdm(dat_list = gamma_dat, family = "gamma",
                                        method = "gllvm",
                                        refresh = 0, iter = 500, chains = 2))
test_that("gamma censor model runs", {
  expect_s3_class(gamma_mod, "jsdmStanFit")
  expect_s3_class(gamma_mod$family, "jsdmStanFamily")
  expect_named(gamma_mod$family, c("family","params","params_dataresp",
                                   "preds","data_list",
                                   "censoring","censoring_data"))
})
test_that("print works", {
  expect_output(print(gamma_mod))
})

test_that("summary works", {
  gamma_cens_summ <- summary(gamma_mod)
  expect_equal(colnames(gamma_cens_summ), c(
    "mean", "sd", "15%", "85%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  gamma_cens_summ <- summary(gamma_mod,
                        pars = "beta", regexp = TRUE,
                        prob_quantiles = c(0.25, 0.5, 0.75)
  )
  expect_equal(colnames(gamma_cens_summ), c(
    "mean", "sd", "25%", "50%", "75%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  expect_match(rownames(gamma_cens_summ), "beta", all = TRUE)
})


test_that("nuts_params works", {
  expect_named(nuts_params(gamma_mod), c("Chain", "Iteration", "Parameter", "Value"))
})

test_that("log_posterior works", {
  expect_named(log_posterior(gamma_mod), c("Chain", "Iteration", "Value"))
})

test_that("rhat works", {
  expect_type(rhat(gamma_mod), "double")
  expect_named(rhat(gamma_mod))
})

test_that("neff_ratio works", {
  expect_type(neff_ratio(gamma_mod), "double")
  expect_named(neff_ratio(gamma_mod))
})

test_that("update works", {
  cens_matrix <- matrix(runif(9*75,2,3),nrow=75)
  cens_ID <- 1*(gamma_dat$Y < cens_matrix)
  suppressWarnings(gamma_mod2 <- update(gamma_mod,
                                        cens_ID = cens_ID,
                                        refresh = 0))
  expect_s3_class(gamma_mod2, "jsdmStanFit")
  expect_s3_class(gamma_mod2$family, "jsdmStanFamily")
  expect_named(gamma_mod2$family, c("family","params","params_dataresp",
                                    "preds","data_list",
                                    "censoring","censoring_data"))

})


# lognormal censoring ####
suppressWarnings(lnorm_mod <- stan_jsdm(dat_list = lnorm_dat, family = "lognormal",
                                        method = "mglmm",
                                        prior = jsdm_prior(betas = "normal(0,0.25)",
                                                           sigma = "normal(0,0.5)"),
                                        refresh = 0, iter = 500, chains = 2))
test_that("lnorm censor model runs", {
  expect_s3_class(lnorm_mod, "jsdmStanFit")
  expect_s3_class(lnorm_mod$family, "jsdmStanFamily")
  expect_named(lnorm_mod$family, c("family","params","params_dataresp",
                                   "preds","data_list",
                                   "censoring","censoring_data"))
})
test_that("print works", {
  expect_output(print(lnorm_mod))
})

test_that("summary works", {
  lnorm_cens_summ <- summary(lnorm_mod)
  expect_equal(colnames(lnorm_cens_summ), c(
    "mean", "sd", "15%", "85%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  lnorm_cens_summ <- summary(lnorm_mod,
                             pars = "beta", regexp = TRUE,
                             prob_quantiles = c(0.25, 0.5, 0.75)
  )
  expect_equal(colnames(lnorm_cens_summ), c(
    "mean", "sd", "25%", "50%", "75%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  expect_match(rownames(lnorm_cens_summ), "beta", all = TRUE)
})


test_that("nuts_params works", {
  expect_named(nuts_params(lnorm_mod), c("Chain", "Iteration", "Parameter", "Value"))
})

test_that("log_posterior works", {
  expect_named(log_posterior(lnorm_mod), c("Chain", "Iteration", "Value"))
})

test_that("rhat works", {
  expect_type(rhat(lnorm_mod), "double")
  expect_named(rhat(lnorm_mod))
})

test_that("neff_ratio works", {
  expect_type(neff_ratio(lnorm_mod), "double")
  expect_named(neff_ratio(lnorm_mod))
})

test_that("update works", {
  suppressWarnings(lnorm_mod2 <- update(lnorm_mod, censor_points = rep(0.7,8),
                                        refresh = 0))
  expect_s3_class(lnorm_mod2, "jsdmStanFit")
  expect_s3_class(lnorm_mod2$family, "jsdmStanFamily")
  expect_named(lnorm_mod2$family, c("family","params","params_dataresp",
                                   "preds","data_list",
                                   "censoring","censoring_data"))

})
