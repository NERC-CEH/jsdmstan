set.seed(5497834)
bern_sim_data <- gllvm_sim_data(N = 300, S = 9, D = 2, K = 2, family = "bern")
bern_pred_data <- matrix(rnorm(100 * 2), nrow = 100)
colnames(bern_pred_data) <- c("V1", "V2")
suppressWarnings(bern_fit <- stan_gllvm(
  dat_list = bern_sim_data, family = "bern",
  refresh = 0, iter = 1000, chains = 2
))

test_that("posterior linpred errors appropriately", {
  expect_error(
    posterior_linpred(bern_fit, newdata_type = "F"),
    "Currently only data on covariates is supported."
  )


  expect_error(posterior_linpred(bern_fit, list_index = "bored"))

  bern_fit_emptydataslot <- bern_fit
  bern_fit_emptydataslot$data_list <- list()
  expect_error(posterior_linpred(bern_fit_emptydataslot))

  bad_newdata <- matrix(rnorm(100 * 2), nrow = 100)
  colnames(bad_newdata) <- c("NEWV1", "NEWV2")
  expect_error(
    posterior_linpred(bern_fit, newdata = bad_newdata),
    "New data does not have matching column names to model fit"
  )
})

test_that("posterior predictive errors appropriately", {
  expect_error(
    posterior_predict(bern_fit, newdata_type = "F"),
    "Currently only data on covariates is supported."
  )

  expect_error(
    posterior_predict(bern_fit, draw_ids = c(-5,-1,0)),
    "draw_ids must be a vector of positive integers"
  )

  expect_error(
    posterior_predict(bern_fit, draw_ids = c(1e5,6e4)),
    "is greater than number of iterations"
  )

  expect_error(posterior_predict(bern_fit, list_index = "bored"))
})

test_that("posterior_(lin)pred returns correct messages",{
  expect_message(
    posterior_predict(bern_fit, ndraws = 100, draw_ids = 50:80),
    "Both ndraws and draw_ids have been specified, ignoring ndraws")

  expect_warning(
    posterior_linpred(bern_fit, ndraws = 1e5),
    "There are fewer samples than ndraws specified")
})

test_that("posterior_(lin)pred works with gllvm", {
  bern_pred <- posterior_predict(bern_fit, ndraws = 100)

  expect_length(bern_pred, 100)
  expect_false(any(sapply(bern_pred, anyNA)))
  expect_false(any(sapply(bern_pred, function(x) x < 0)))

  bern_pred2 <- posterior_predict(bern_fit,
    newdata = bern_pred_data,
    ndraws = 50, list_index = "species"
  )

  expect_length(bern_pred2, 9)
  expect_false(any(sapply(bern_pred2, anyNA)))
  expect_false(any(sapply(bern_pred2, function(x) x < 0)))
})

pois_sim_data <- mglmm_sim_data(N = 300, S = 9, K = 2, family = "pois")
pois_pred_data <- matrix(rnorm(100 * 2), nrow = 100)
colnames(pois_pred_data) <- c("V1", "V2")
suppressWarnings(pois_fit <- stan_mglmm(
  dat_list = pois_sim_data, family = "pois",
  refresh = 0, chains = 2
))
test_that("posterior_(lin)pred works with mglmm", {
  pois_pred <- posterior_predict(pois_fit, ndraws = 100)

  expect_length(pois_pred, 100)
  expect_false(any(sapply(pois_pred, anyNA)))
  expect_false(any(sapply(pois_pred, function(x) x < 0)))

  pois_pred2 <- posterior_predict(pois_fit,
    newdata = pois_pred_data,
    ndraws = 50, list_index = "species"
  )

  expect_length(pois_pred2, 9)
  expect_false(any(sapply(pois_pred2, anyNA)))
  expect_false(any(sapply(pois_pred2, function(x) x < 0)))
})
