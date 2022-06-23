test_that("jsdm_statsummary errors correctly", {
  expect_error(jsdm_statsummary("eugsuf"),
               "jsdm_summary only works for jsdmStanFit objects")
  test_fit <- jsdmStanFit_empty()
  expect_error(jsdm_statsummary(test_fit, species = -3),
               paste("Species must be either a character vector of species names or an",
                     "integer vector of species positions in the input data columns"))
  expect_error(jsdm_statsummary(test_fit, sites = 4.5),
               paste("Sites must be either a character vector of site names or an",
                     "integer vector of sites positions in the input data columns"))
})

bern_sim_data <- mglmm_sim_data(N = 300, S = 9, K = 2, family = "bern")
colnames(bern_sim_data$Y) <- LETTERS[1:9]
rownames(bern_sim_data$Y) <- paste0("Site",1:300)
rownames(bern_sim_data$X) <- paste0("Site",1:300)
suppressWarnings(bern_fit <- stan_mglmm(dat_list = bern_sim_data, family = "bern",
                                        refresh = 0))

test_that("jsdm_statsummary errors with wrong names", {
  expect_error(jsdm_statsummary(bern_fit, species = LETTERS[7:12]),
               "Species specified are not found in the model fit object")

  expect_error(jsdm_statsummary(bern_fit, sites = as.character(7:12)),
               "Sites specified are not found in the model fit object")
})


gauss_sim_data <- gllvm_sim_data(N = 150, D = 2, S = 16, family = "gaussian")
colnames(gauss_sim_data$Y) <- LETTERS[1:16]
rownames(gauss_sim_data$Y) <- paste0("Site",1:150)
suppressWarnings(gauss_fit <- stan_gllvm(dat_list = gauss_sim_data, family = "gaussian",
                                         refresh = 0))


test_that("jsdm_statsummary returns correct form of output", {
  expect_type(jsdm_statsummary(bern_fit),
              "double")
  expect_equal(dim(jsdm_statsummary(bern_fit, summary_stat = "mean")), c(8000,300))
  expect_equal(dim(jsdm_statsummary(bern_fit, ndraws = 10, sites = 1:100)),
               c(10,100))
  test_stat <- jsdm_statsummary(bern_fit, post_type = "predict",
                                draw_ids = seq(1,1000,100))
  expect_false(anyNA(test_stat))
  expect_false(any(test_stat<0))


  expect_type(jsdm_statsummary(gauss_fit, post_type = "predict"),
              "double")
  expect_equal(dim(jsdm_statsummary(gauss_fit, species = 1:5,
                                    summary_stat = function(x) quantile(x, 0.1))),
               c(8000,150))
  expect_equal(dim(jsdm_statsummary(gauss_fit, ndraws = 30, calc_over = "species")),
               c(30,16))
  test_stat <- jsdm_statsummary(gauss_fit, post_type = "predict",
                                draw_ids = seq(1,1000,100))
  expect_false(anyNA(test_stat))
})

# pp_check
test_that("pp_check errors correctly", {
  test_fit <- jsdmStanFit_empty()
  expect_error(pp_check(test_fit, plotfun = "fiona"),
               "is not a valid ppc type")
})
