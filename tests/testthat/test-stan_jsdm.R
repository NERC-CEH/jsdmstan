test_that("formula structure works", {
  set.seed(5488935)
  mglmm_data <- mglmm_sim_data(N = 100, S = 8, family = "gaussian", K = 3)
  df <- as.data.frame(mglmm_data$X)

  mglmm_fit <- stan_jsdm(~V1 + V2 + V3, data = df, Y = mglmm_data$Y,
                          family = "gaussian", method = "mglmm",
                          refresh = 0)

  df$V3 <- as.factor(cut(df$V3, c(-Inf,-0.5,0.5,Inf), c("a","b","c")))
  mglmm_fit <- stan_mglmm(~V1 + V2 + V3, data = df, Y = mglmm_data$Y,
                          family = "gaussian",
                          refresh = 0)
})
