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

test_that("loo works", {
  expect_warning(
    mglmm_loo <- loo(mglmm_fit),
    "Pareto k diagnostic"
  )
  expect_s3_class(mglmm_loo, "psis_loo")
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

  expect_s3_class(mglmm_fit, "jsdmStanFit")

  # neg bin
  mglmm_data <- mglmm_sim_data(N = 87, S = 12, K = 2, family = "neg_bin")
  suppressWarnings(mglmm_fit <- stan_mglmm(
    Y = mglmm_data$Y, X = mglmm_data$X,
    family = "neg_bin",
    refresh = 0, chains = 2, iter = 200
  ))

  expect_s3_class(mglmm_fit, "jsdmStanFit")
})

# stan_jsdm site_intercept tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mglmm_data <- mglmm_sim_data(N = 100, S = 8, family = "gaussian", K = 3,
                             site_intercept = "ungrouped")
df <- as.data.frame(mglmm_data$X)
grps <- rep(1:20, each = 5)

gllvm_data <- gllvm_sim_data(N = 100, S = 8, family = "bern", D = 3,
                             site_intercept = "ungrouped")
gllvm_data$grps <- rep(1:20, each = 5)
gllvm_data$ngrp <- 20

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
