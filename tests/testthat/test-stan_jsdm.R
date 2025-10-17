# stan_jsdm formula tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set.seed(5488935)
mglmm_data <- mglmm_sim_data(N = 30, S = 8, family = "gaussian", K = 3)
df <- as.data.frame(mglmm_data$X)

suppressWarnings(mglmm_fit <- stan_jsdm(~ V1 + V2 + V3,
  data = df, Y = mglmm_data$Y,
  family = "gaussian", method = "mglmm",
  refresh = 0, chains = 2, iter = 200
))
test_that("formula structure works", {
  expect_s3_class(mglmm_fit, "jsdmStanFit")

  df$V3 <- as.factor(cut(df$V3, c(-Inf, -0.5, 0.5, Inf), c("a", "b", "c")))
  suppressWarnings(mglmm_fit <- stan_mglmm(~ V1 * V2 + V3,
    data = df, Y = mglmm_data$Y,
    family = "gaussian",
    refresh = 0, chains = 1, iter = 200
  ))
  expect_s3_class(mglmm_fit, "jsdmStanFit")
})

# jsdmstanfit object interaction tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("print works", {
  expect_output(print(mglmm_fit))
})

test_that("summary works", {
  mglmm_summ <- summary(mglmm_fit)
  expect_equal(colnames(mglmm_summ), c(
    "mean", "sd", "15%", "85%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  mglmm_summ <- summary(mglmm_fit,
    pars = "beta", regexp = TRUE,
    prob_quantiles = c(0.25, 0.5, 0.75)
  )
  expect_equal(colnames(mglmm_summ), c(
    "mean", "sd", "25%", "50%", "75%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  expect_match(rownames(mglmm_summ), "beta", all = TRUE)
})


test_that("nuts_params works", {
  expect_named(nuts_params(mglmm_fit), c("Chain", "Iteration", "Parameter", "Value"))
})

test_that("log_posterior works", {
  expect_named(log_posterior(mglmm_fit), c("Chain", "Iteration", "Value"))
})

test_that("rhat works", {
  expect_type(rhat(mglmm_fit), "double")
  expect_named(rhat(mglmm_fit))
})

test_that("neff_ratio works", {
  expect_type(neff_ratio(mglmm_fit), "double")
  expect_named(neff_ratio(mglmm_fit))
})


# GLLVM tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("stan_gllvm fails with wrong inputs", {
  test_Y <- matrix(sample(0:1, size = 100, replace = TRUE),
    nrow = 20
  )
  expect_error(
    stan_gllvm(Y = test_Y, family = "bern"),
    "Must have at least one latent variable"
  )

  expect_error(
    stan_gllvm(
      Y = matrix(rnorm(100),
        nrow = 20
      ),
      D = 2, family = "bern"
    ),
    "Y matrix is not binary"
  )
  expect_error(
    stan_gllvm(
      Y = matrix(rnorm(100),
        nrow = 20
      ),
      D = 2, family = "pois"
    ),
    "Y matrix is not composed of integers"
  )

  expect_error(
    stan_gllvm(Y = matrix(rnorm(100), nrow = 20),
               X = matrix(rnorm(60), nrow = 20),
               family = "bern", D = -1),
    "Must have at least one latent variable"
  )

  expect_error(
    suppressMessages(stan_gllvm(Y = matrix(sample.int(10,100,
                                                      replace = TRUE),
                                           nrow = 20),
                                X = matrix(rnorm(60), nrow = 20),
                                family = "binomial", D = 2, Ntrials = "a")),
    "Ntrials must be a positive integer"
  )

  expect_error(
    suppressMessages(stan_gllvm(Y = matrix(sample.int(10,100,
                                                      replace = TRUE), nrow = 20),
                                X = matrix(rnorm(60), nrow = 20),
                                family = "binomial", D = 2, Ntrials = c(1,3))),
    "Ntrials must be of length"
  )
})

test_that("stan_gllvm returns right type of object", {
  set.seed(04012022)
  gllvm_data <- gllvm_sim_data(
    N = 36, S = 16, D = 2, K = 5,
    site_intercept = "ungrouped", family = "bern"
  )
  suppressWarnings(gllvm_fit <- stan_jsdm(
    dat_list = gllvm_data, refresh = 0,
    method = "gllvm", site_intercept = "ungrouped",
    family = "bern",
    beta_param = "unstruct",
    chains = 2, iter = 200
  ))

  expect_s3_class(gllvm_fit, "jsdmStanFit")

  # poisson
  gllvm_data <- jsdm_sim_data(
    N = 55, S = 12, D = 3, K = 3, family = "pois",
    method = "gllvm"
  )
  suppressWarnings(gllvm_fit <- stan_gllvm(
    dat_list = gllvm_data, refresh = 0,
    chains = 2, iter = 200, family = "pois"
  ))

  expect_s3_class(gllvm_fit, "jsdmStanFit")

  # neg binomial
  gllvm_data <- gllvm_sim_data(N = 20, S = 8, D = 2, K = 2, family = "neg_bin")
  suppressWarnings(gllvm_fit <- stan_gllvm(
    Y = as.data.frame(gllvm_data$Y), X = as.data.frame(gllvm_data$X),
    D = gllvm_data$D, refresh = 0, chains = 2, iter = 200,
    family = "neg_b"
  ))

  expect_s3_class(gllvm_fit, "jsdmStanFit")

  # binomial
  gllvm_data <- gllvm_sim_data(N = 20, S = 8, D = 2, K = 2, family = "binomial",
                               Ntrials = 20)
  suppressWarnings(gllvm_fit <- stan_gllvm(
    Y = as.data.frame(gllvm_data$Y), X = as.data.frame(gllvm_data$X),
    D = gllvm_data$D, refresh = 0, chains = 2, iter = 200,
    family = "binomial", Ntrials = 20
  ))

  expect_s3_class(gllvm_fit, "jsdmStanFit")
})

# MGLMM tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("stan_mglmm fails with wrong inputs", {
  expect_error(
    stan_mglmm(
      Y = matrix(rnorm(100),
        nrow = 20
      ),
      family = "bern"
    ),
    "Y matrix is not binary"
  )
  expect_error(
    stan_mglmm(
      Y = matrix(rnorm(100),
        nrow = 20
      ),
      family = "pois"
    ),
    "Y matrix is not composed of integers"
  )
  expect_error(
    stan_mglmm(
      Y = matrix(sample.int(100,500,replace=TRUE), nrow = 100),
      X = data.frame(V1 = rnorm(100),
                     V2 = rnorm(100)),
      family = "binomial"
    ),
    "Number of trials must be specified"
  )
})

test_that("stan_mglmm returns right type of object", {
  set.seed(04012022)
  mglmm_data <- mglmm_sim_data(
    N = 32, S = 8, family = "bern",
    site_intercept = "ungrouped"
  )
  suppressWarnings(mglmm_fit <- stan_mglmm(
    Y = mglmm_data$Y, X = mglmm_data$X,
    family = "bern",
    beta_param = "unstruct",
    refresh = 0, site_intercept = "ungrouped", chains = 2, iter = 200
  ))

  expect_s3_class(mglmm_fit, "jsdmStanFit")

  # poisson
  mglmm_data <- jsdm_sim_data(
    method = "mglmm", N = 27, S = 8, K = 3,
    family = "poisson"
  )
  suppressWarnings(mglmm_fit <- stan_jsdm(
    X = NULL, dat_list = mglmm_data,
    family = "poisson",
    refresh = 0, method = "mglmm", iter = 200,
    chains = 2
  ))

  expect_error(stan_jsdm(dat_list = mglmm_data, family = "binomial",
                         method = "mglmm"),
               "Binomial models require Ntrials")

  expect_s3_class(mglmm_fit, "jsdmStanFit")

  # neg bin
  mglmm_data <- mglmm_sim_data(N = 87, S = 12, K = 2, family = "neg_bin")
  suppressWarnings(mglmm_fit <- stan_mglmm(
    Y = mglmm_data$Y, X = mglmm_data$X,
    family = "neg_bin",
    refresh = 0, chains = 2, iter = 200
  ))

  expect_s3_class(mglmm_fit, "jsdmStanFit")

  # binomial
  mglmm_data <- mglmm_sim_data(N = 51, S = 6, K = 2, family = "binomial",
                               Ntrials=sample.int(20,51,replace = TRUE))
  suppressWarnings(mglmm_fit <- stan_mglmm(
    dat_list = mglmm_data,
    family = "binomial",
    refresh = 0, chains = 2, iter = 200
  ))

  expect_s3_class(mglmm_fit, "jsdmStanFit")

  # binomial

})

# stan_jsdm site_intercept tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mglmm_data <- mglmm_sim_data(N = 100, S = 8, family = "gaussian", K = 3,
                             site_intercept = "ungrouped")
df <- as.data.frame(mglmm_data$X)
grps <- rep(1:20, each = 5)

test_that("site intercept errors correctly", {
  expect_error(stan_mglmm(~ V1 + V2, data = df, Y = mglmm_data$Y,
                          site_intercept = "fauh",family = "gaussian"))
  expect_error(stan_gllvm(~V1 + V2, data = df, D = 2, Y = mglmm_data$Y,
                          site_intercept = "grouped", family = "gaussian"),
               "If site_intercept is grouped then groups must be supplied to site_groups")
})

test_that("site intercept models run", {
  suppressWarnings(mglmm_fit <- stan_mglmm(~ V1 * V2,
                                           data = df, Y = mglmm_data$Y,
                                           site_intercept = "ungrouped",
                                           family = "gaussian",
                                           beta_param = "unstruct",
                                           refresh = 0, chains = 1, iter = 200
  ))
  expect_s3_class(mglmm_fit, "jsdmStanFit")


  suppressWarnings(mglmm_fit <- stan_mglmm(~ V1 * V2,
                                           data = df, Y = mglmm_data$Y,
                                           site_intercept = "grouped",
                                           site_groups = grps,
                                           family = "gaussian",
                                           beta_param = "cor",
                                           refresh = 0, chains = 1, iter = 200
  ))
  expect_s3_class(mglmm_fit, "jsdmStanFit")

})


set.seed(9598098)
zip_sim_data <- gllvm_sim_data(N = 100, S = 7, D = 2, K = 2, family = "zi_poisson",
                               zi_param = "covariate")

test_that("zi_poisson works okay", {
  suppressWarnings(zip_fit <- stan_gllvm(
    dat_list = zip_sim_data, family = "zi_poisson",zi_param="covariate",
    refresh = 0, chains = 2, iter = 200
  ))
  expect_s3_class(zip_fit, "jsdmStanFit")
  expect_s3_class(zip_fit$family, "jsdmStanFamily")
  expect_output(print(zip_fit$family),
                "is modelled in response to")
  expect_named(zip_fit$family,
               c("family" ,"params" ,"params_dataresp","preds","data_list"))
  expect_named(zip_fit$family$data_list, "zi_X")
})

set.seed(7684116)
zinb_sim_data <- mglmm_sim_data(N = 67, S = 6, K = 3, family = "zi_neg_binomial",
                                zi_param = "covariate", shp_param = "covariate",
                                shp_k = 1, zi_k = 1)

test_that("shp and zi models error properly",{

  expect_error(stan_mglmm(dat_list = zinb_sim_data, family = "zi_neg_binomial",
                          zi_param = "fail"))

  expect_error(stan_mglmm(dat_list = zinb_sim_data, family = "zi_neg_binomial",
                          shp_param = "fail"))

  expect_error(stan_mglmm(dat_list = zinb_sim_data, family = "poisson",
                          shp_param = "covariate"),
               "family parameter in response to data")

  expect_error(stan_jsdm(~V1 + V2, data = as.data.frame(zinb_sim_data$X),
                         method = "mglmm", family = "zi_neg_binomial",
                         Y = zinb_sim_data$Y,
                         zi_formula = ~ 1),
               "Intercept-only models")

  expect_error(stan_jsdm(~V1 + V2, data = as.data.frame(zinb_sim_data$X),
                         method = "gllvm", family = "gaussian",
                         Y = zinb_sim_data$Y,
                         shp_formula = ~ 1),
               "Intercept-only models")
})


zinb_data <- as.data.frame(cbind(zinb_sim_data$X,
                                 zinb_sim_data$shp_X[,-1],
                                 zinb_sim_data$zi_X[,-1]))
colnames(zinb_data) <- c("X1","X2","X3","S1","Z1")

test_that("shp models run okay", {
  suppressWarnings(zinb_fit <- stan_mglmm(~X1+X2+X3, data = zinb_data, Y = zinb_sim_data$Y,
                                          zi_formula = ~Z1, shp_formula = ~S1,
                                          family = "zi_neg_binomial",
                                          refresh = 0, chains = 2, iter = 200
  ))
  expect_s3_class(zinb_fit, "jsdmStanFit")
  expect_s3_class(zinb_fit$family, "jsdmStanFamily")
  expect_output(print(zinb_fit$family),
                "is modelled in response to")
  expect_named(zinb_fit$family,
               c("family" ,"params" ,"params_dataresp","preds","data_list"))
  expect_named(zinb_fit$family$data_list, c("zi_X","shp_X"))
})
