
#' Make stancode for the jsdm model
#'
#' This function returns the Stan code used to fit the model as specified by the data
#' list, family and method.
#'
#' @details
#'   Environmental covariate effects (\code{"betas"}) can be parameterised in two
#'   ways. With the \code{"cor"} parameterisation all covariate effects are assumed
#'   to be constrained by a correlation matrix between the covariates. With the
#'   \code{"unstruct"} parameterisation all covariate effects are assumed to draw
#'   from a simple distribution with no correlation structure. Both parameterisations
#'   can be modified using the prior object.
#'
#' @param method The method, one of \code{"gllvm"} or \code{"mglmm"}
#' @param family The family, one of "\code{"gaussian"}, \code{"bernoulli"},
#'   \code{"poisson"} or \code{"neg_binomial"}
#' @param prior The prior, given as the result of a call to [jsdm_prior()]
#' @param log_lik Whether the log likelihood should be calculated in the generated
#'   quantities (by default \code{TRUE}), required for loo
#' @param site_intercept Whether a site intercept should be included, potential
#'   values \code{"none"} (no site intercept), \code{"grouped"} (a site intercept
#'   with hierarchical grouping) or \code{"ungrouped"} (site intercept with no
#'   grouping)
#' @param beta_param The parameterisation of the environmental covariate effects, by
#'   default \code{"cor"}. See details for further information.
#'
#' @return A character vector of Stan code, class "jsdmstan_model"
#' @export
#'
#' @examples
#' jsdm_stancode(family = "gaussian", method = "gllvm")
#' jsdm_stancode(family = "poisson", method = "mglmm")
#'
jsdm_stancode <- function(method, family, prior = jsdm_prior(),
                          log_lik = TRUE, site_intercept = "none",
                          beta_param = "cor") {
  # checks
  family <- match.arg(family, c("gaussian", "bernoulli", "poisson", "neg_binomial"))
  method <- match.arg(method, c("gllvm", "mglmm"))
  beta_param <- match.arg(beta_param, c("cor","unstruct"))
  site_intercept <- match.arg(site_intercept, c("none","grouped","ungrouped"))
  if (class(prior)[1] != "jsdmprior") {
    stop("Prior must be given as a jsdmprior object")
  }

  # data processing steps
  scode <- .modelcode(
    method = method, family = family,
    phylo = FALSE, prior = prior, log_lik = log_lik, site_intercept = site_intercept,
    beta_param = beta_param
  )
  class(scode) <- c("jsdmstan_model", "character")
  return(scode)
}


.modelcode <- function(method, family, phylo, prior, log_lik, site_intercept,
                       beta_param) {
  model_functions <- "
"
  data <- paste(
    " int<lower=1> N; // Number of sites
  int<lower=1> S; // Number of species
 ", ifelse(method == "gllvm",
      "int<lower=1> D; // Number of latent dimensions", ""
    ),
    "
  int<lower=0> K; // Number of predictor variables
  matrix[N, K] X; // Predictor matrix
", ifelse(site_intercept == "grouped",
            "
  int<lower=1> ngrp; // Number of groups in site intercept
  int<lower=0, upper = ngrp> grps[N]; // Vector matching sites to groups
  ",""),
  switch(family,
      "gaussian" = "real",
      "bernoulli" = "int<lower=0,upper=1>",
      "neg_binomial" = "int<lower=0>",
      "poisson" = "int<lower=0>"
    ), "Y[N,S]; //Species matrix"
  )
  transformed_data <- ifelse(method == "gllvm", "
  // Ensures identifiability of the model - no rotation of factors
  int<lower=1> M;
  M = D * (S - D) + choose(D, 2) + D;
", "")



  site_inter_par <- switch(site_intercept,
                           "ungrouped" = "
  // Site intercepts
  real a_bar;
  real<lower=0> sigma_a;
  vector[N] a;",
  "none" = "",
  "grouped" = "
  // Site intercepts
  real a_bar;
  real<lower=0> sigma_a;
  vector[ngrp] a;")
  species_pars <- switch(beta_param, "cor" = "
  //betas are hierarchical with covariance model
  vector<lower=0>[K] sigmas_preds;
  matrix[K, S] z_preds;
  // covariance matrix on betas by predictors
  cholesky_factor_corr[K] cor_preds;", "unstruct" = "
  matrix[K, S] betas;")
  mglmm_spcov_pars <- "
  // species covariances
  vector<lower=0>[S] sigmas_species;
  matrix[S, N] z_species;
  cholesky_factor_corr[S] cor_species;"
  gllvm_pars <- "
  // Factor parameters
  vector[M] L; // Non-zero factor loadings
  real<lower=0> sigma_L; // variance of species loadings
  // Latent variables
  matrix[D, N] LV_uncor; // Per-site latent variable"
  var_pars <- switch(family,
    "gaussian" = "
  real<lower=0> sigma; // Gaussian parameters",
    "bernoulli" = "",
    "neg_binomial" = "
  real<lower=0> kappa; // neg_binomial parameters",
    "poisson" = ""
  )

  pars <- paste(
    site_inter_par, species_pars,
    switch(method,
      "gllvm" = gllvm_pars,
      "mglmm" = mglmm_spcov_pars
    ),
    var_pars
  )
  transformed_pars <- paste(if(beta_param == "cor"){ "
  // covariance matrix on betas by preds
  matrix[K, S] betas;
  "} else {""}, switch(method,
    "gllvm" = "
  // Construct factor loading matrix
  matrix[S, D] Lambda_uncor;
  // Constraints to allow identifiability of loadings
  for (i in 1:(D-1)) {
    for (j in (i+1):(D)){
      Lambda_uncor[i,j] = 0;
    }
  }
  {
    int index;
    index = 0;
    for (j in 1:D) {
      for (i in j:S) {
        index = index + 1;
        Lambda_uncor[i, j] = L[index];
      }
    }
  }
  ",
    "mglmm" = "
  matrix[N, S] u;
  u = (diag_pre_multiply(sigmas_species, cor_species) * z_species)';
  "
  ), if(beta_param == "cor") {"
  betas = diag_pre_multiply(sigmas_preds, cor_preds) * z_preds;
"} else {""})



  gllvm_model <- switch(site_intercept, "ungrouped" = "
  // model
  matrix[N, S] LV_sum = ((Lambda_uncor * sigma_L) * LV_uncor)';
  matrix[N, S] alpha = rep_matrix(a_bar + a * sigma_a, S);
  mu = alpha + (X * betas) + LV_sum;
  ", "none" = "
  // model
  matrix[N, S] LV_sum = ((Lambda_uncor * sigma_L) * LV_uncor)';
  mu = (X * betas) + LV_sum;
  ", "grouped" = "
  matrix[N, S] alpha;
  matrix[N, S] LV_sum = ((Lambda_uncor * sigma_L) * LV_uncor)';
  {
    vector[ngrp] theta = a_bar + a * sigma_a;
    for (n in 1:N){
      alpha[n,] = rep_row_vector(theta[grps[n]],S);
    }
  }
  mu = alpha + (X * betas) + LV_sum;
  ")
  mglmm_model <- switch(site_intercept, "ungrouped" = "
  // model
  matrix[N, S] alpha = rep_matrix(a_bar + a * sigma_a, S);
  mu = alpha + (X * betas) + u;
  ", "none" = "
  // model
  mu = (X * betas) + u;
  ", "grouped" = "
  matrix[N, S] alpha;
  {
    vector[ngrp] theta = a_bar + a * sigma_a;
    for (n in 1:N){
      alpha[n,] = rep_row_vector(theta[grps[n]],S);
    }
  }
  mu = alpha + (X * betas) + u;
  ")
  model <- paste("
  matrix[N,S] mu;
  ", switch(method,
    "gllvm" = gllvm_model,
    "mglmm" = mglmm_model
  ))
  model_priors <- paste(
    ifelse(site_intercept %in% c("ungrouped","grouped"), paste("
  // Site-level intercept priors
  a ~ ", prior[["a"]], ";
  a_bar ~ ", prior[["a_bar"]], ";
  sigma_a ~ ", prior[["sigma_a"]], ";
  "), "
  "), switch(beta_param, "cor" = paste(
  "
  // Species parameter priors
  sigmas_preds ~ ", prior[["sigmas_preds"]], ";
  to_vector(z_preds) ~ ", prior[["z_preds"]], ";
  // covariance matrix priors
  cor_preds ~ ", prior[["cor_preds"]], ";
"), "unstruct" = paste(
  "
  // Species parameter priors
  to_vector(betas) ~ ", prior[["betas"]],";
  ")),
    switch(method,
      "gllvm" = paste("
  // Factor priors
  to_vector(LV_uncor) ~ ", prior[["LV"]], ";
  L ~ ", prior[["L"]], ";
  sigma_L ~ ", prior[["sigma_L"]], "; // Variance of factor loadings
"),
      "mglmm" = paste("
  // Species parameter priors
  sigmas_species ~ ", prior[["sigmas_species"]], ";
  to_vector(z_species) ~ ", prior[["z_species"]], ";
  cor_species ~ ", prior[["cor_species"]], ";
")
    ),
    switch(family,
      "gaussian" = paste("
  //Standard deviation parameters
  sigma ~ ", prior[["sigma"]], ";
"),
      "neg_binomial" = paste("
  //Scale parameter
  kappa ~ ", prior[["kappa"]], ";
"),
      "bern" = "",
      "poisson" = ""
    )
  )
  model_pt2 <- paste(
    "
  for(i in 1:N) Y[i,] ~ ",
    switch(family,
      "gaussian" = "normal(mu[i,], sigma);",
      "bernoulli" = "bernoulli_logit(mu[i,]);",
      "neg_binomial" = "neg_binomial_2_log(mu[i,], kappa);",
      "poisson" = "poisson_log(mu[i,]);"
    )
  )

  generated_quantities <- paste(
    ifelse(isTRUE(log_lik), "
  // Calculate linear predictor, y_rep, log likelihoods for LOO
  matrix[N, S] log_lik;
  ", ""),
    ifelse(method == "gllvm", "
  // Sign correct factor loadings and factors
  matrix[D, N] LV;
  matrix[S, D] Lambda;
  for(d in 1:D){
    if(Lambda_uncor[d,d] < 0){
      Lambda[,d] = -1 * Lambda_uncor[,d];
      LV[d,] = -1 * LV_uncor[d,];
    } else {
      Lambda[,d] = Lambda_uncor[,d];
      LV[d,] = LV_uncor[d,];
    }
  }", ""), ifelse(isTRUE(log_lik), paste(
      "
  {
    matrix[N, S] linpred;", switch(site_intercept, "ungrouped" = paste("
    linpred = rep_matrix(a_bar + a * sigma_a, S) + (X * betas) +",
    switch(method,
           "gllvm" = "((Lambda_uncor * sigma_L) * LV_uncor)'",
           "mglmm" = "u"
    ), ";
    "), "none" = paste("
    linpred = (X * betas) +",
    switch(method,
           "gllvm" = "((Lambda_uncor * sigma_L) * LV_uncor)'",
           "mglmm" = "u"
      ), ";
    "), "grouped" = paste("
  matrix[N, S] alpha;
  for (n in 1:N){
    alpha[n,] = rep_row_vector(a[grps[n]],S);
  }
  linpred = alpha + (X * betas) +",
  switch(method,
         "gllvm" = "((Lambda_uncor * sigma_L) * LV_uncor)'",
         "mglmm" = "u"
  ), ";
    ")),"
      for(i in 1:N) {
      for(j in 1:S) {
        log_lik[i, j] = ",
      switch(family,
        "gaussian" = "normal_lpdf",
        "bernoulli" = "bernoulli_logit_lpmf",
        "neg_binomial" = "neg_binomial_2_log_lpmf",
        "poisson" = "poisson_log_lpmf"
      ),
      "(Y[i, j] | linpred[i, j]",
      switch(family,
        "gaussian" = ", sigma)",
        "bernoulli" = ")",
        "neg_binomial" = ", kappa)",
        "poisson" = ")"
      ), ";
      }
    }
  }
  "
    ), "")
  )

  res <- paste(
    "//Generated by jsdmstan\n",
    "functions{\n",
    model_functions,
    "\n}\ndata{\n",
    data,
    "\n}\ntransformed data{\n",
    transformed_data,
    "\n}\nparameters{\n",
    pars,
    "\n}\ntransformed parameters{\n",
    transformed_pars,
    "\n}\nmodel{\n",
    model, "\n", model_priors, "\n", model_pt2, "\n",
    "\n}\ngenerated quantities{\n",
    generated_quantities,
    "\n}\n\n"
  )

  return(res)
}

#' @export
#' @describeIn jsdm_stancode A printing function for jsdmstan_model objects
#' @param x The jsdm_stancode object
#' @param ... Currently unused
print.jsdmstan_model <- function(x, ...) {
  cat(x)
}
