
# 1D nfs testing ####
set.seed(2)
dat <- mgcv::gamSim(1,n=400,dist="normal",scale=2)


# make a 1D example with less noise
dat$y2 <- dat$f2 + rnorm(400, sd=1)
dat$y3 <- dat$f2 + rnorm(400, sd=1)
dat$y4 <- dat$f2 + rnorm(400, sd=1)
dat$y5 <- dat$f2 + rnorm(400, sd=1)

Y <- dat[,c("y2","y3","y4","y5")]


suppressWarnings(
    nfs1_fit <- stan_jsdm(~s(x2), data = dat,
                          Y = Y, family ="gaussian", method = "mglmm",
                          chains = 2, iter = 200, refresh = 0))
test_that("model runs with one nfs spline", {
  expect_s3_class(nfs1_fit, "jsdmStanFit")

  expect_named(nfs1_fit$preds$spl_smooth$nfs,
               c("N","X","S1","nSr","nSc","nsp","Sr","Sn","nterms","ncoef","sm","sm_data"))

  expect_s3_class(nfs1_fit$preds$spl_smooth$nfs$sm[[1]],"mgcv.smooth")
})

test_that("print works", {
  expect_output(print(nfs1_fit))
})

test_that("summary works", {
  nfs1_fit_summ <- summary(nfs1_fit)
  expect_equal(colnames(nfs1_fit_summ), c(
    "mean", "sd", "15%", "85%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  nfs1_fit_summ <- summary(nfs1_fit,
                             pars = "beta", regexp = TRUE,
                             prob_quantiles = c(0.25, 0.5, 0.75)
  )
  expect_equal(colnames(nfs1_fit_summ), c(
    "mean", "sd", "25%", "50%", "75%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  expect_match(rownames(nfs1_fit_summ), "beta", all = TRUE)
})


test_that("nuts_params works", {
  expect_named(nuts_params(nfs1_fit), c("Chain", "Iteration", "Parameter", "Value"))
})

test_that("log_posterior works", {
  expect_named(log_posterior(nfs1_fit), c("Chain", "Iteration", "Value"))
})

test_that("rhat works", {
  expect_type(rhat(nfs1_fit), "double")
  expect_named(rhat(nfs1_fit))
})

test_that("neff_ratio works", {
  expect_type(neff_ratio(nfs1_fit), "double")
  expect_named(neff_ratio(nfs1_fit))
})


test_that("posterior_(lin)pred works", {
  nfs1_pred <- posterior_predict(nfs1_fit, ndraws = 100)

  expect_length(nfs1_pred, 100)
  expect_false(any(sapply(nfs1_pred, anyNA)))
})


test_that("pp_check works", {
  nfs1_pp <- pp_check(nfs1_fit, ndraws = 10)

  expect_s3_class(nfs1_pp, "gg")
})

test_that("smoothplot works", {
  nfs1_sp <- smoothplot(nfs1_fit, ndraws = 10)

  expect_s3_class(nfs1_sp[[1]], "gg")

  nfs1_sp <- smoothplot(nfs1_fit, ndraws = 10, summarise = "mean")

  expect_s3_class(nfs1_sp[[1]], "gg")
})


# 2 spline nfs testing ####
set.seed(2)
dat$y6 <- dat$f1 + dat$f2 + rnorm(400, sd=1)
dat$y7 <- dat$f1 + dat$f2 + rnorm(400, sd=1)
dat$y8 <- dat$f1 + dat$f2 + rnorm(400, sd=1)
dat$y9 <- dat$f1 + dat$f2 + rnorm(400, sd=1)

Y <- dat[,c("y6","y7","y8","y9")]

suppressWarnings(
  nfs2_fit <- stan_jsdm(~s(x2)+s(x1), data = dat,
                        Y = Y, family ="gaussian", method = "mglmm",
                        chains = 2, iter = 200, refresh = 0))
test_that("model runs with one nfs spline", {
  expect_s3_class(nfs2_fit, "jsdmStanFit")

  expect_named(nfs2_fit$preds$spl_smooth$nfs,
               c("N","X","S1","nSr","nSc","nsp","Sr","Sn","nterms","ncoef","sm","sm_data"))

  expect_s3_class(nfs2_fit$preds$spl_smooth$nfs$sm[[1]],"mgcv.smooth")
})

test_that("print works", {
  expect_output(print(nfs2_fit))
})

test_that("summary works", {
  nfs2_fit_summ <- summary(nfs2_fit)
  expect_equal(colnames(nfs2_fit_summ), c(
    "mean", "sd", "15%", "85%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  nfs2_fit_summ <- summary(nfs2_fit,
                           pars = "beta", regexp = TRUE,
                           prob_quantiles = c(0.25, 0.5, 0.75)
  )
  expect_equal(colnames(nfs2_fit_summ), c(
    "mean", "sd", "25%", "50%", "75%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  expect_match(rownames(nfs2_fit_summ), "beta", all = TRUE)
})


test_that("nuts_params works", {
  expect_named(nuts_params(nfs2_fit), c("Chain", "Iteration", "Parameter", "Value"))
})

test_that("log_posterior works", {
  expect_named(log_posterior(nfs2_fit), c("Chain", "Iteration", "Value"))
})

test_that("rhat works", {
  expect_type(rhat(nfs2_fit), "double")
  expect_named(rhat(nfs2_fit))
})

test_that("neff_ratio works", {
  expect_type(neff_ratio(nfs2_fit), "double")
  expect_named(neff_ratio(nfs2_fit))
})


test_that("posterior_(lin)pred works", {
  nfs2_pred <- posterior_predict(nfs2_fit, ndraws = 100)

  expect_length(nfs2_pred, 100)
  expect_false(any(sapply(nfs2_pred, anyNA)))
})


test_that("pp_check works", {
  nfs2_pp <- pp_check(nfs2_fit, ndraws = 10)

  expect_s3_class(nfs2_pp, "gg")
})

test_that("smoothplot works", {
  nfs2_sp <- smoothplot(nfs2_fit, ndraws = 10)

  expect_s3_class(nfs2_sp[[1]], "gg")
  expect_length(nfs2_sp, 2)


  nfs2_sp <- smoothplot(nfs2_fit, ndraws = 10, summarise = "median")

  expect_s3_class(nfs2_sp[[1]], "gg")
  expect_length(nfs2_sp, 2)
})


# 2D nfs testing ####
set.seed(2)
dat$y6 <- dat$f1 + dat$f2 + rnorm(400, sd=1)
dat$y7 <- dat$f1 + dat$f2 + rnorm(400, sd=1)
dat$y8 <- dat$f1 + dat$f2 + rnorm(400, sd=1)
dat$y9 <- dat$f1 + dat$f2 + rnorm(400, sd=1)

Y <- dat[,c("y6","y7","y8","y9")]

suppressWarnings(
  nfs3_fit <- stan_jsdm(~s(x2, x1), data = dat,
                        Y = Y, family ="gaussian", method = "mglmm",
                        chains = 2, iter = 200, refresh = 0))
test_that("model runs with one nfs spline", {
  expect_s3_class(nfs3_fit, "jsdmStanFit")

  expect_named(nfs3_fit$preds$spl_smooth$nfs,
               c("N","X","S1","nSr","nSc","nsp","Sr","Sn","nterms","ncoef","sm","sm_data"))

  expect_s3_class(nfs3_fit$preds$spl_smooth$nfs$sm[[1]],"mgcv.smooth")
})

test_that("print works", {
  expect_output(print(nfs3_fit))
})

test_that("summary works", {
  nfs3_fit_summ <- summary(nfs3_fit)
  expect_equal(colnames(nfs3_fit_summ), c(
    "mean", "sd", "15%", "85%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  nfs3_fit_summ <- summary(nfs3_fit,
                           pars = "beta", regexp = TRUE,
                           prob_quantiles = c(0.25, 0.5, 0.75)
  )
  expect_equal(colnames(nfs3_fit_summ), c(
    "mean", "sd", "25%", "50%", "75%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  expect_match(rownames(nfs3_fit_summ), "beta", all = TRUE)
})


test_that("nuts_params works", {
  expect_named(nuts_params(nfs3_fit), c("Chain", "Iteration", "Parameter", "Value"))
})

test_that("log_posterior works", {
  expect_named(log_posterior(nfs3_fit), c("Chain", "Iteration", "Value"))
})

test_that("rhat works", {
  expect_type(rhat(nfs3_fit), "double")
  expect_named(rhat(nfs3_fit))
})

test_that("neff_ratio works", {
  expect_type(neff_ratio(nfs3_fit), "double")
  expect_named(neff_ratio(nfs3_fit))
})


test_that("posterior_(lin)pred works", {
  nfs2_pred <- posterior_predict(nfs3_fit, ndraws = 100)

  expect_length(nfs2_pred, 100)
  expect_false(any(sapply(nfs2_pred, anyNA)))
})


test_that("pp_check works", {
  nfs2_pp <- pp_check(nfs3_fit, ndraws = 10)

  expect_s3_class(nfs2_pp, "gg")
})

test_that("smoothplot works", {
  nfs2_sp <- smoothplot(nfs3_fit, ndraws = 10)

  expect_s3_class(nfs2_sp[[1]], "gg")
  expect_length(nfs2_sp, 1)

  nfs2_sp <- smoothplot(nfs3_fit, ndraws = 10, summarise = "mean")

  expect_s3_class(nfs2_sp[[1]], "gg")
  expect_length(nfs2_sp, 1)
})

# 1D non-species fs testing ####
# simulated example from ?gam
set.seed(0)
## simulate data...
f1 <- function(x,a=2,b=-1) exp(a * x)+b
n <- 39;nf <- 13
fac <- as.factor(rep(LETTERS[1:nf],each=3))
x1 <- runif(39)
a <- rnorm(nf)*.2 + 2;b <- rnorm(nf)*.5
f <- f1(x1,a[fac],b[fac])
fac <- factor(fac)
y1 <- f + rnorm(n)
y2 <- f + rnorm(n)
y3 <- f + rnorm(n)
## so response depends on global smooths of x0 and
## x2, and a smooth of x1 for each level of fac.
df <- data.frame(x1 = x1, x2 = fac)
Y <- matrix(c(y1,y2,y3),nrow=39)
colnames(Y) <- LETTERS[20:22]
## fit model...
suppressWarnings(
  fs0_fit <- stan_jsdm(~s(x1,x2, bs = "fs"), data = df, Y = Y,
                       family = "gaussian", method = "mglmm",
                       chains = 2, iter = 200, refresh = 0))

test_that("model runs with one fs spline", {
  expect_s3_class(fs0_fit, "jsdmStanFit")

  expect_named(fs0_fit$preds$spl_smooth$fs,
               c("N","X","S1","nSr","nSc","nsp","Sr","Sn","nterms","ncoef","sm","sm_data"))

  expect_s3_class(fs0_fit$preds$spl_smooth$fs$sm[[1]],"mgcv.smooth")
})

test_that("print works", {
  expect_output(print(fs0_fit))
})

test_that("summary works", {
  fs0_fit_summ <- summary(fs0_fit)
  expect_equal(colnames(fs0_fit_summ), c(
    "mean", "sd", "15%", "85%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  fs0_fit_summ <- summary(fs0_fit,
                          pars = "beta", regexp = TRUE,
                          prob_quantiles = c(0.25, 0.5, 0.75)
  )
  expect_equal(colnames(fs0_fit_summ), c(
    "mean", "sd", "25%", "50%", "75%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  expect_match(rownames(fs0_fit_summ), "beta", all = TRUE)
})


test_that("nuts_params works", {
  expect_named(nuts_params(fs0_fit), c("Chain", "Iteration", "Parameter", "Value"))
})

test_that("log_posterior works", {
  expect_named(log_posterior(fs0_fit), c("Chain", "Iteration", "Value"))
})

test_that("rhat works", {
  expect_type(rhat(fs0_fit), "double")
  expect_named(rhat(fs0_fit))
})

test_that("neff_ratio works", {
  expect_type(neff_ratio(fs0_fit), "double")
  expect_named(neff_ratio(fs0_fit))
})


test_that("posterior_(lin)pred works", {
  fs0_pred <- posterior_predict(fs0_fit, ndraws = 100)

  expect_length(fs0_pred, 100)
  expect_false(any(sapply(fs0_pred, anyNA)))
})


test_that("pp_check works", {
  fs0_pp <- pp_check(fs0_fit, ndraws = 10)

  expect_s3_class(fs0_pp, "gg")
  expect_length(fs0_pp, 1)
})

test_that("smoothplot works", {
  fs0_sp <- smoothplot(fs0_fit, ndraws = 10)

  expect_s3_class(fs0_sp[[1]], "gg")

  fs0_sp <- smoothplot(fs0_fit, ndraws = 10, summarise = "mean")

  expect_s3_class(fs0_sp[[1]], "gg")
})


# 1D fs testing ####
# simulated example from ?gam
set.seed(0)
## simulate data...
f1 <- function(x,a=2,b=-1) exp(a * x)+b
n <- 500;nf <- 10
fac <- as.factor(rep(LETTERS[1:10],each=50))
x1 <- rep(runif(50),10)
a <- rnorm(nf)*.2 + 2;b <- rnorm(nf)*.5
f <- f1(x1,a[fac],b[fac])
fac <- factor(fac)
y <- f + rnorm(n)
## so response depends on global smooths of x0 and
## x2, and a smooth of x1 for each level of fac.
df <- data.frame(x1 = x1[1:50])
Y <- matrix(y, nrow = 50)
colnames(Y) <- LETTERS[1:10]
## fit model...
suppressWarnings(
  fs1_fit <- stan_jsdm(~s(x1,species, bs = "fs"), data = df, Y = Y,
                     family = "gaussian", method = "mglmm",
                     chains = 2, iter = 200, refresh = 0))

test_that("model runs with one fs spline", {
  expect_s3_class(fs1_fit, "jsdmStanFit")

  expect_named(fs1_fit$preds$spl_smooth$fs,
               c("N","X","S1","nSr","nSc","nsp","Sr","Sn","nterms","ncoef","sm","sm_data"))

  expect_s3_class(fs1_fit$preds$spl_smooth$fs$sm[[1]],"mgcv.smooth")
})

test_that("print works", {
  expect_output(print(fs1_fit))
})

test_that("summary works", {
  fs1_fit_summ <- summary(fs1_fit)
  expect_equal(colnames(fs1_fit_summ), c(
    "mean", "sd", "15%", "85%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  fs1_fit_summ <- summary(fs1_fit,
                           pars = "beta", regexp = TRUE,
                           prob_quantiles = c(0.25, 0.5, 0.75)
  )
  expect_equal(colnames(fs1_fit_summ), c(
    "mean", "sd", "25%", "50%", "75%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  expect_match(rownames(fs1_fit_summ), "beta", all = TRUE)
})


test_that("nuts_params works", {
  expect_named(nuts_params(fs1_fit), c("Chain", "Iteration", "Parameter", "Value"))
})

test_that("log_posterior works", {
  expect_named(log_posterior(fs1_fit), c("Chain", "Iteration", "Value"))
})

test_that("rhat works", {
  expect_type(rhat(fs1_fit), "double")
  expect_named(rhat(fs1_fit))
})

test_that("neff_ratio works", {
  expect_type(neff_ratio(fs1_fit), "double")
  expect_named(neff_ratio(fs1_fit))
})


test_that("posterior_(lin)pred works", {
  fs1_pred <- posterior_predict(fs1_fit, ndraws = 100)

  expect_length(fs1_pred, 100)
  expect_false(any(sapply(fs1_pred, anyNA)))
})


test_that("pp_check works", {
  fs1_pp <- pp_check(fs1_fit, ndraws = 10)

  expect_s3_class(fs1_pp, "gg")
  expect_length(fs1_pp, 1)
})

test_that("smoothplot works", {
  fs1_sp <- smoothplot(fs1_fit, ndraws = 10)

  expect_s3_class(fs1_sp[[1]], "gg")

  fs1_sp <- smoothplot(fs1_fit, ndraws = 10, summarise = "mean")

  expect_s3_class(fs1_sp[[1]], "gg")
})


# nfs plus fs testing ####
# simulated example from ?gam

## fit model...
suppressWarnings(
  nfsfs1_fit <- stan_jsdm(~s(x1) + s(x1,species, bs = "fs"), data = df, Y = Y,
                          family = "gaussian", method = "mglmm",
                          chains = 2, iter = 200, refresh = 0))

test_that("model runs with one fs spline", {
  expect_s3_class(nfsfs1_fit, "jsdmStanFit")

  expect_named(nfsfs1_fit$preds$spl_smooth$fs,
               c("N","X","S1","nSr","nSc","nsp","Sr","Sn","nterms","ncoef","sm","sm_data"))

  expect_s3_class(nfsfs1_fit$preds$spl_smooth$fs$sm[[1]],"mgcv.smooth")
  expect_s3_class(nfsfs1_fit$preds$spl_smooth$nfs$sm[[1]],"mgcv.smooth")
})

test_that("print works", {
  expect_output(print(nfsfs1_fit))
})

test_that("summary works", {
  nfsfs1_fit_summ <- summary(nfsfs1_fit)
  expect_equal(colnames(nfsfs1_fit_summ), c(
    "mean", "sd", "15%", "85%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  nfsfs1_fit_summ <- summary(nfsfs1_fit,
                          pars = "beta", regexp = TRUE,
                          prob_quantiles = c(0.25, 0.5, 0.75)
  )
  expect_equal(colnames(nfsfs1_fit_summ), c(
    "mean", "sd", "25%", "50%", "75%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  expect_match(rownames(nfsfs1_fit_summ), "beta", all = TRUE)
})


test_that("nuts_params works", {
  expect_named(nuts_params(nfsfs1_fit), c("Chain", "Iteration", "Parameter", "Value"))
})

test_that("log_posterior works", {
  expect_named(log_posterior(nfsfs1_fit), c("Chain", "Iteration", "Value"))
})

test_that("rhat works", {
  expect_type(rhat(nfsfs1_fit), "double")
  expect_named(rhat(nfsfs1_fit))
})

test_that("neff_ratio works", {
  expect_type(neff_ratio(nfsfs1_fit), "double")
  expect_named(neff_ratio(nfsfs1_fit))
})


test_that("posterior_(lin)pred works", {
  fs1_pred <- posterior_predict(nfsfs1_fit, ndraws = 100)

  expect_length(fs1_pred, 100)
  expect_false(any(sapply(fs1_pred, anyNA)))
})


test_that("pp_check works", {
  fs1_pp <- pp_check(nfsfs1_fit, ndraws = 10)

  expect_s3_class(fs1_pp, "gg")
})

test_that("smoothplot works", {
  fs1_sp <- smoothplot(nfsfs1_fit, ndraws = 10)

  expect_s3_class(fs1_sp[[1]], "gg")
  expect_length(fs1_sp, 2)

  fs1_sp <- smoothplot(nfsfs1_fit, ndraws = 10, summarise = "median")

  expect_s3_class(fs1_sp[[1]], "gg")
  expect_length(fs1_sp, 2)
})



# 2 spline fs testing ####

df$x2 <- rnorm(50)
suppressWarnings(
  fs2_fit <- stan_jsdm(~s(x2, species, bs = "fs")+s(x1, species, bs = "fs"), data = df,
                        Y = Y, family ="gaussian", method = "mglmm",
                        chains = 2, iter = 200, refresh = 0))
test_that("model runs with one fs spline", {
  expect_s3_class(fs2_fit, "jsdmStanFit")

  expect_named(fs2_fit$preds$spl_smooth$fs,
               c("N","X","S1","nSr","nSc","nsp","Sr","Sn","nterms","ncoef","sm","sm_data"))

  expect_s3_class(fs2_fit$preds$spl_smooth$fs$sm[[1]],"mgcv.smooth")
})

test_that("print works", {
  expect_output(print(fs2_fit))
})

test_that("summary works", {
  fs2_fit_summ <- summary(fs2_fit)
  expect_equal(colnames(fs2_fit_summ), c(
    "mean", "sd", "15%", "85%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  fs2_fit_summ <- summary(fs2_fit,
                           pars = "beta", regexp = TRUE,
                           prob_quantiles = c(0.25, 0.5, 0.75)
  )
  expect_equal(colnames(fs2_fit_summ), c(
    "mean", "sd", "25%", "50%", "75%", "Rhat",
    "Bulk.ESS", "Tail.ESS"
  ))
  expect_match(rownames(fs2_fit_summ), "beta", all = TRUE)
})


test_that("nuts_params works", {
  expect_named(nuts_params(fs2_fit), c("Chain", "Iteration", "Parameter", "Value"))
})

test_that("log_posterior works", {
  expect_named(log_posterior(fs2_fit), c("Chain", "Iteration", "Value"))
})

test_that("rhat works", {
  expect_type(rhat(fs2_fit), "double")
  expect_named(rhat(fs2_fit))
})

test_that("neff_ratio works", {
  expect_type(neff_ratio(fs2_fit), "double")
  expect_named(neff_ratio(fs2_fit))
})


test_that("posterior_(lin)pred works", {
  fs2_pred <- posterior_predict(fs2_fit, ndraws = 100)

  expect_length(fs2_pred, 100)
  expect_false(any(sapply(fs2_pred, anyNA)))
})


test_that("pp_check works", {
  fs2_pp <- pp_check(fs2_fit, ndraws = 10)

  expect_s3_class(fs2_pp, "gg")
})

test_that("smoothplot works", {
  fs2_sp <- smoothplot(fs2_fit, ndraws = 10)

  expect_s3_class(fs2_sp[[1]], "gg")
  expect_length(fs2_sp, 2)

  fs2_sp <- smoothplot(fs2_fit, ndraws = 10, summarise = "sd")

  expect_s3_class(fs2_sp[[1]], "gg")
  expect_length(fs2_sp, 2)
})


# 2D fs testing ####
# This model is insanely slow, so will not be included in casual tests
# suppressWarnings(
#   fs3_fit <- stan_jsdm(~s(x2, x1, species, bs  = "fs"), data = df,
#                         Y = Y, family ="gaussian", method = "mglmm",
#                         chains = 2, iter = 200, refresh = 0))
# test_that("model runs with one fs spline", {
#   expect_s3_class(fs3_fit, "jsdmStanFit")
#
#   expect_named(fs3_fit$preds$spl_smooth$fs,
#                c("N","X","S1","nSr","nSc","nsp","Sr","Sn","nterms","ncoef","sm","sm_data"))
#
#   expect_s3_class(fs3_fit$preds$spl_smooth$fs$sm[[1]],"mgcv.smooth")
# })
#
# test_that("print works", {
#   expect_output(print(fs3_fit))
# })
#
# test_that("summary works", {
#   fs3_fit_summ <- summary(fs3_fit)
#   expect_equal(colnames(fs3_fit_summ), c(
#     "mean", "sd", "15%", "85%", "Rhat",
#     "Bulk.ESS", "Tail.ESS"
#   ))
#   fs3_fit_summ <- summary(fs3_fit,
#                            pars = "beta", regexp = TRUE,
#                            prob_quantiles = c(0.25, 0.5, 0.75)
#   )
#   expect_equal(colnames(fs3_fit_summ), c(
#     "mean", "sd", "25%", "50%", "75%", "Rhat",
#     "Bulk.ESS", "Tail.ESS"
#   ))
#   expect_match(rownames(fs3_fit_summ), "beta", all = TRUE)
# })
#
#
# test_that("nuts_params works", {
#   expect_named(nuts_params(fs3_fit), c("Chain", "Iteration", "Parameter", "Value"))
# })
#
# test_that("log_posterior works", {
#   expect_named(log_posterior(fs3_fit), c("Chain", "Iteration", "Value"))
# })
#
# test_that("rhat works", {
#   expect_type(rhat(fs3_fit), "double")
#   expect_named(rhat(fs3_fit))
# })
#
# test_that("neff_ratio works", {
#   expect_type(neff_ratio(fs3_fit), "double")
#   expect_named(neff_ratio(fs3_fit))
# })
#
#
# test_that("posterior_(lin)pred works", {
#   fs2_pred <- posterior_predict(fs3_fit, ndraws = 100)
#
#   expect_length(fs2_pred, 100)
#   expect_false(any(sapply(fs2_pred, anyNA)))
# })
#
#
# test_that("pp_check works", {
#   fs2_pp <- pp_check(fs3_fit, ndraws = 10)
#
#   expect_s3_class(fs2_pp, "gg")
# })
#
# test_that("smoothplot works", {
#   fs2_sp <- smoothplot(fs3_fit, ndraws = 10)
#
#   expect_s3_class(fs2_sp[[1]], "gg")
#   expect_length(fs2_sp, 1)
#
#   fs2_sp <- smoothplot(fs3_fit, ndraws = 10, summarise = "mean")
#
#   expect_s3_class(fs2_sp[[1]], "gg")
#   expect_length(fs2_sp, 1)
# })
