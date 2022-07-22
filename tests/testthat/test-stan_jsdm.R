# stan_jsdm formula tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set.seed(5488935)
mglmm_data <- mglmm_sim_data(N = 100, S = 8, family = "gaussian", K = 3)
df <- as.data.frame(mglmm_data$X)

suppressWarnings(mglmm_fit <- stan_jsdm(~V1 + V2 + V3, data = df, Y = mglmm_data$Y,
                                        family = "gaussian", method = "mglmm",
                                        refresh = 0, chains = 2, iter = 200))
test_that("formula structure works", {

  expect_s3_class(mglmm_fit, "jsdmStanFit")

  df$V3 <- as.factor(cut(df$V3, c(-Inf,-0.5,0.5,Inf), c("a","b","c")))
  suppressWarnings(mglmm_fit <- stan_mglmm(~V1*V2 + V3, data = df, Y = mglmm_data$Y,
                                           family = "gaussian",
                                           refresh = 0, chains = 1, iter = 200))
  expect_s3_class(mglmm_fit, "jsdmStanFit")
})

# jsdmstanfit object interaction tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("print works",{
  expect_output(print(mglmm_fit))
})

test_that("summary works", {
  mglmm_summ <- summary(mglmm_fit)
  expect_equal(colnames(mglmm_summ), c("mean", "sd", "15%", "85%", "Rhat",
                                       "Bulk.ESS", "Tail.ESS"))
  mglmm_summ <- summary(mglmm_fit, pars = "beta", regexp = TRUE,
                        prob_quantiles = c(0.25,0.5,0.75))
  expect_equal(colnames(mglmm_summ), c("mean", "sd", "25%", "50%", "75%", "Rhat",
                                       "Bulk.ESS", "Tail.ESS"))
  expect_match(rownames(mglmm_summ),"beta", all = TRUE)

})

test_that("update works", {
  mglmm_data <- mglmm_sim_data(N = 50, S = 5, family = "gaussian", K = 2)
  suppressWarnings(mglmm_fit2 <- update(mglmm_fit, newY = mglmm_data$Y,
                                        newX = mglmm_data$X,
                                        refresh = 0,
                                        chains = 2, iter = 200))
  suppressWarnings(mglmm_fit3 <- update(mglmm_fit,
                                        refresh = 0, iter = 100))

  expect_s3_class(mglmm_fit2, "jsdmStanFit")
  expect_s3_class(mglmm_fit3, "jsdmStanFit")

  jsdm_empty <- jsdmStanFit_empty()
  expect_error(update(jsdm_empty),
               "Update requires the original data to be saved in the model object")
})

test_that("nuts_params works", {
  expect_named(nuts_params(mglmm_fit), c("Chain","Iteration","Parameter","Value"))
})

test_that("log_posterior works", {
  expect_named(log_posterior(mglmm_fit), c("Chain","Iteration","Value"))
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
  expect_warning(mglmm_loo <- loo(mglmm_fit),
                 "Pareto k diagnostic")
  expect_s3_class(mglmm_loo, "psis_loo")
})

# GLLVM tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  suppressWarnings(gllvm_fit <- stan_jsdm(dat_list = gllvm_data, refresh = 0,
                                          method = "gllvm", site_intercept = TRUE,
                                          family = "bern",
                                          chains = 2, iter = 500))

  expect_s3_class(gllvm_fit, "jsdmStanFit")

  # poisson
  gllvm_data <- jsdm_sim_data(N = 200, S = 12, D = 3, K = 3, family = "pois",
                              method = "gllvm")
  suppressWarnings(gllvm_fit <- stan_gllvm(dat_list = gllvm_data, refresh = 0,
                                           chains = 2, iter = 500, family = "pois"))

  expect_s3_class(gllvm_fit, "jsdmStanFit")

  # neg binomial
  gllvm_data <- gllvm_sim_data(N = 100, S = 8, D = 2, family = "neg_bin")
  gllvm_fit <- stan_gllvm(Y = gllvm_data$Y, X = gllvm_data$X,
                          D = gllvm_data$D, refresh = 0, chains= 2,
                          family = "neg_b")

  expect_s3_class(gllvm_fit, "jsdmStanFit")

})

# MGLMM tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("stan_mglmm fails with wrong inputs", {
  expect_error(stan_mglmm(Y = matrix(rnorm(100),
                                     nrow=20),
                          family = "bern"),
               "Y matrix is not binary")
  expect_error(stan_mglmm(Y = matrix(rnorm(100),
                                     nrow=20),
                          family = "pois"),
               "Y matrix is not composed of integers")
})

test_that("stan_mglmm returns right type of object", {
  set.seed(04012022)
  mglmm_data <- mglmm_sim_data(N = 100, S = 8, family = "bern",
                               site_intercept = TRUE)
  mglmm_fit <- stan_mglmm(Y = mglmm_data$Y, X = mglmm_data$X,
                          family = "bern",
                          refresh = 0, site_intercept = TRUE, chains = 2)

  expect_s3_class(mglmm_fit, "jsdmStanFit")

  # poisson
  mglmm_data <- jsdm_sim_data(method = "mglmm", N = 200, S = 8, K = 3,
                              family = "poisson")
  suppressWarnings(mglmm_fit <- stan_jsdm(X = NULL, dat_list = mglmm_data,
                                          family = "poisson",
                                          refresh = 0, method = "mglmm", iter = 200,
                                          chains= 2))

  expect_s3_class(mglmm_fit, "jsdmStanFit")

  # neg bin
  mglmm_data <- mglmm_sim_data(N = 127, S = 12, K = 2, family = "neg_bin")
  suppressWarnings(mglmm_fit <- stan_mglmm(Y = mglmm_data$Y, X = mglmm_data$X,
                                           family = "neg_bin",
                                           refresh = 0, chains = 2, iter = 200))

  expect_s3_class(mglmm_fit, "jsdmStanFit")
})


# Phylo tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("phylo models error appropriately",{
  expect_error(mglmm_sim_data(N = 500, S = 12, phylo = TRUE, family = "poisson"),
               "Need to specify delta and covar arguments for phylo")
  phylo_data <- mglmm_sim_data(N = 500, S = 12, phylo = TRUE, family = "poisson",
                               delta = 1e-5, covar = "matern_05")
  expect_error(stan_jsdm(X = phylo_data$X, Y = phylo_data$Y, method = "gllvm",
                         family = "poisson", phylo = TRUE),
               "Phylo only supported with MGLMM")
  expect_error(stan_jsdm(X = phylo_data$X, Y = phylo_data$Y, method = "mglmm",
                         family = "poisson", phylo = TRUE),
               "Phylo must be either FALSE or a matrix")
  expect_error(stan_jsdm(X = phylo_data$X, Y = phylo_data$Y, method = "mglmm",
                         family = "poisson", phylo = matrix(1:100,nrow=20)),
               "Phylo must be a square matrix with dimensions equal")

})

test_that("phylo poisson matern 25 models run", {
  set.seed(46153123)
  phylo_data <- mglmm_sim_data(N = 123, S = 12, K = 3, phylo = TRUE,
                               family = "poisson",
                               delta = 1e-5, covar = "matern_25")
  expect_named(phylo_data, c("Y", "pars", "N","S","D","K","X",
                             "site_intercept","Dmat","delta","nu05"))
  suppressWarnings(phylo_fit <- stan_jsdm(~V1 + V2 + V3,
                                          data = as.data.frame(phylo_data$X),
                                          Y = phylo_data$Y,
                                          family = "poisson",
                                          method = "mglmm",
                                          phylo = phylo_data$Dmat,
                                          covar = "matern_25",
                                          iter = 200, chains = 2, refresh = 0,
                                          log_lik = FALSE))
  expect_s3_class(phylo_fit, "jsdmStanFit")
}
)
test_that("phylo sq exponential bern models run",{
  set.seed(46153123)
  phylo_data2 <- jsdm_sim_data(N = 100, S = 10, phylo = TRUE, family = "bern",
                               delta = 1e-4, covar = "sq_exponential",
                               method = "mglmm")
  expect_named(phylo_data2, c("Y", "pars", "N","S","D","K","X",
                              "site_intercept","Dmat","delta","nu05"))
  suppressWarnings(phylo_fit2 <- stan_jsdm(X = NULL, dat_list = phylo_data2,
                                           method = "mglmm",
                                           family = "bern",
                                           phylo = phylo_data2$Dmat,
                                           covar = "sq_exponential",
                                           iter = 200, chains = 2, refresh = 0))
  expect_s3_class(phylo_fit2, "jsdmStanFit")
})
test_that("phylo gaussian exponential models run", {
  phylo_data3 <- jsdm_sim_data(N = 100, S = 10, phylo = TRUE, family = "gaussian",
                               delta = 1e-5, covar = "exponential",
                               method = "mglmm")
  expect_named(phylo_data3, c("Y", "pars", "N","S","D","K","X",
                              "site_intercept","Dmat","delta","nu05"))
  suppressWarnings(phylo_fit3 <- stan_jsdm(X = NULL, dat_list = phylo_data3,
                                           method = "mglmm",
                                           family = "gaussian",
                                           phylo = phylo_data3$Dmat,
                                           covar = "exponential",
                                           iter = 200, chains = 2, refresh = 0,
                                           log_lik = FALSE))
  expect_s3_class(phylo_fit3, "jsdmStanFit")
})
test_that("phylo negbin matern 15 models run", {
  phylo_data4 <- jsdm_sim_data(N = 100, S = 10, phylo = TRUE,
                               family = "neg_binomial",
                               delta = 1e-5, covar = "matern_15",
                               method = "mglmm")
  expect_named(phylo_data4, c("Y", "pars", "N","S","D","K","X",
                              "site_intercept","Dmat","delta","nu05"))
  suppressWarnings(phylo_fit4 <- stan_jsdm(X = NULL, dat_list = phylo_data4,
                                           method = "mglmm",
                                           family = "neg_binomial",
                                           phylo = phylo_data4$Dmat,
                                           covar = "matern_15",
                                           iter = 200, chains = 2, refresh = 0,
                                           log_lik = FALSE))
  expect_s3_class(phylo_fit4, "jsdmStanFit")
})


