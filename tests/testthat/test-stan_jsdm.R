set.seed(5488935)
mglmm_data <- mglmm_sim_data(N = 100, S = 8, family = "gaussian", K = 3)
df <- as.data.frame(mglmm_data$X)

mglmm_fit <- stan_jsdm(~V1 + V2 + V3, data = df, Y = mglmm_data$Y,
                       family = "gaussian", method = "mglmm",
                       refresh = 0, control = list(adapt_delta = 0.99),
                       iter = 6000)
test_that("formula structure works", {

  expect_s3_class(mglmm_fit, "jsdmStanFit")

  df$V3 <- as.factor(cut(df$V3, c(-Inf,-0.5,0.5,Inf), c("a","b","c")))
  mglmm_fit <- stan_mglmm(~V1*V2 + V3, data = df, Y = mglmm_data$Y,
                          family = "gaussian",
                          refresh = 0, control = list(adapt_delta = 0.99),
                          iter = 6000)
  expect_s3_class(mglmm_fit, "jsdmStanFit")
})

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
  mglmm_fit2 <- update(mglmm_fit, newY = mglmm_data$Y, newX = mglmm_data$X,
                       refresh = 0, control = list(adapt_delta = 0.99),
                       iter = 4000)

  expect_s3_class(mglmm_fit2, "jsdmStanFit")
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

test_that("phylo models error appropriately",{
  expect_error(mglmm_sim_data(N = 500, S = 12, phylo = TRUE, family = "poisson"),
               "Need to specify delta and nu05 arguments for phylo")
  phylo_data <- mglmm_sim_data(N = 500, S = 12, phylo = TRUE, family = "poisson",
                               delta = 1e-5, nu05 = 1)
  expect_error(stan_jsdm(X = phylo_data$X, Y = phylo_data$Y, method = "gllvm",
                         family = "poisson", phylo = TRUE),
               "Phylo only supported with MGLMM")
  expect_error(stan_jsdm(X = phylo_data$X, Y = phylo_data$Y, method = "mglmm",
                         family = "poisson", phylo = TRUE),
               "Phylo must be either FALSE or a matrix")
  expect_error(stan_jsdm(X = phylo_data$X, Y = phylo_data$Y, method = "mglmm",
                         family = "poisson", phylo = matrix(1:100,nrow=20)),
               "Phylo must be a square matrix with dimensions equal")
  # phylo_fit <- stan_jsdm(X = phylo_data$X, Y = phylo_data$Y, method = "mglmm",
  #                        phylo = phylo_data$Dmat)
})
