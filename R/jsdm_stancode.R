
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
#' @param family is the response family, must be one of \code{"gaussian"},
#'   \code{"neg_binomial"}, \code{"poisson"}, \code{"binomial"},
#'   \code{"bernoulli"}, \code{"zi_poisson"}, or
#'    \code{"zi_neg_binomial"}. Regular expression
#'   matching is supported.
#' @param prior The prior, given as the result of a call to [jsdm_prior()]
#' @param log_lik Whether the log likelihood should be calculated in the generated
#'   quantities (by default \code{TRUE}), required for loo
#' @param site_intercept Whether a site intercept should be included, potential
#'   values \code{"none"} (no site intercept), \code{"grouped"} (a site intercept
#'   with hierarchical grouping) or \code{"ungrouped"} (site intercept with no
#'   grouping)
#' @param beta_param The parameterisation of the environmental covariate effects, by
#'   default \code{"cor"}. See details for further information.
#' @param zi_param For the zero-inflated families, whether the zero-inflation parameter
#'   is a species-specific constant (default, \code{"constant"}), or varies by
#'   environmental covariates (\code{"covariate"}).
#' @param shp_param For the families with shape parameters, whether the shape parameter
#'   is a species-specific constant (default, \code{"constant"}), or varies by
#'   environmental covariates (\code{"covariate"}).
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
                          beta_param = "cor", zi_param = "constant",
                          shp_param = "constant") {
  # checks
  family <- match.arg(family, c("gaussian", "bernoulli", "poisson",
                                "neg_binomial","binomial","zi_poisson",
                                "zi_neg_binomial"))
  method <- match.arg(method, c("gllvm", "mglmm"))
  beta_param <- match.arg(beta_param, c("cor","unstruct"))
  site_intercept <- match.arg(site_intercept, c("none","grouped","ungrouped"))
  zi_param <- match.arg(zi_param, c("constant","covariate"))
  shp_param <- match.arg(shp_param, c("constant","covariate"))
  if(shp_param == "covariate" & family %in% c("poisson","bernoulli","binomial",
                                                 "zi_poisson"))
    stop(paste("Modelling the family parameter in response to data only works",
               "for Gaussian and negative binomial families"))
  if (class(prior)[1] != "jsdmprior") {
    stop("Prior must be given as a jsdmprior object")
  }

  # data processing steps
  scode <- .modelcode(
    method = method, family = family,
    phylo = FALSE, prior = prior, log_lik = log_lik, site_intercept = site_intercept,
    beta_param = beta_param, zi_param = zi_param, shp_param = shp_param
  )
  class(scode) <- c("jsdmstan_model", "character")
  return(scode)
}


.modelcode <- function(method, family, phylo, prior, log_lik, site_intercept,
                       beta_param, zi_param, shp_param) {
  model_functions <- "
"
  data <- paste(
    " int<lower=1> N; // Number of sites
  int<lower=1> S; // Number of species
",
ifelse(method == "gllvm",
      " int<lower=1> D; // Number of latent dimensions", ""
    ),
    "
  int<lower=0> K; // Number of predictor variables
  matrix[N, K] X; // Predictor matrix
",
ifelse(site_intercept == "grouped",
            "
  int<lower=1> ngrp; // Number of groups in site intercept
  int<lower=0, upper = ngrp> grps[N]; // Vector matching sites to groups
 ",""),
  switch(family,
      "gaussian" = "real",
      "bernoulli" = "int<lower=0,upper=1>",
      "neg_binomial" = "int<lower=0>",
      "poisson" = "int<lower=0>",
      "zi_poisson" = "int<lower=0>",
      "zi_neg_binomial" = "int<lower=0>",
      "binomial" = "int<lower=0>"
    ), "Y[N,S]; //Species matrix",
 ifelse(family == "binomial",
        "
  int<lower=0> Ntrials[N]; // Number of trials",""),
  ifelse(grepl("zi_", family),"
  int<lower=0> N_zero[S]; // number of zeros per species
  int<lower=0> N_nonzero[S]; //number of nonzeros per species
  int<lower=0> Sum_nonzero; //Total number of nonzeros across all species
  int<lower=0> Sum_zero; //Total number of zeros across all species
  int<lower=0> Y_nz[Sum_nonzero]; //Y values for nonzeros
  int<lower=0> ss[Sum_nonzero]; //species index for Y_nz
  int<lower=0> nn[Sum_nonzero]; //site index for Y_nz
  int<lower=0> sz[Sum_zero]; //species index for Y_z
  int<lower=0> nz[Sum_zero]; //site index for Y_z",""),
ifelse(grepl("zi_", family) & zi_param == "covariate","
  int<lower=1> zi_k; //number of covariates for env effects on zi
  matrix[N, zi_k] zi_X; //environmental covariate matrix for zi",""),
ifelse(shp_param == "covariate","
  int<lower=1> shp_k; //number of covariates for env effects on family parameter
  matrix[N, shp_k] shp_X; //environmental covariate matrix for family parameter","")
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
  corr_matrix[K] cor_preds;", "unstruct" = "
  matrix[K, S] betas;")
  mglmm_spcov_pars <- "
  // species covariances
  vector<lower=0>[S] sigmas_species;
  matrix[S, N] z_species;
  corr_matrix[S] cor_species;"
  gllvm_pars <- "
  // Factor parameters
  vector[M] L; // Non-zero factor loadings
  real<lower=0> sigma_L; // variance of species loadings
  // Latent variables
  matrix[D, N] LV_uncor; // Per-site latent variable"

  var_pars <- switch(family,
    "gaussian" = switch(shp_param, "constant" = "
  real<lower=0> sigma[S]; // Gaussian parameters", "covariate" = "
  matrix[shp_k,S] shp_betas; //environmental effects for family param"),
    "bernoulli" = "",
    "neg_binomial" = switch(shp_param, "constant" = "
  real<lower=0> kappa[S]; // neg_binomial parameters","covariate" = "
  matrix[shp_k,S] shp_betas; //environmental effects for family param"),
    "poisson" = "",
    "zi_poisson" = switch(zi_param,
    "constant" = "
  real<lower=0,upper=1> zi[S]; // zero-inflation parameter",
    "covariate" = "
  matrix[zi_k,S] zi_betas; //environmental effects for zi"),
  "zi_neg_binomial" = switch(zi_param,
                             "constant" = switch(shp_param, "constant" = "
  real<lower=0> kappa[S]; // neg_binomial parameters
  real<lower=0,upper=1> zi[S]; // zero-inflation parameter", "covariate" = "
  matrix[shp_k,S] shp_betas; //environmental effects for family param
  real<lower=0,upper=1> zi[S]; // zero-inflation parameter"),
                             "covariate" = switch(shp_param, "constant" = "
  real<lower=0> kappa[S]; // neg_binomial parameters
  matrix[zi_k,S] zi_betas; //environmental effects for zi", "covariate" = "
  matrix[shp_k,S] shp_betas; //environmental effects for family param
  matrix[zi_k,S] zi_betas; //environmental effects for zi"))
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
  ", ifelse(grepl("zi_",family),paste0("
  real mu_nz[Sum_nonzero];
  real mu_z[Sum_zero];
  int pos;
  int neg;",switch(zi_param,"constant" = "",
                   "covariate" = "
  real zi_nz[Sum_nonzero];
  real zi_z[Sum_zero];"),
                   ifelse(shp_param== "covariate" & family == "zi_neg_binomial", "
  real kappa_nz[Sum_nonzero];
  real kappa_z[Sum_zero];","")),""),
  switch(method,
    "gllvm" = gllvm_model,
    "mglmm" = mglmm_model
  ),ifelse(grepl("zi_",family),paste0(ifelse(zi_param == "covariate", "
  matrix[N,S] zi = zi_X * zi_betas;",""),
  ifelse(shp_param == "covariate", "
  matrix[N,S] kappa = exp(shp_X * shp_betas);",""),"
  for(i in 1:Sum_nonzero){
    mu_nz[i] = mu[nn[i],ss[i]];",
    switch(zi_param, "constant" = "",
           "covariate" = "
    zi_nz[i] = zi[nn[i],ss[i]];"),
    ifelse(shp_param == "covariate" & family == "zi_neg_binomial", "
    kappa_nz[i] = kappa[nn[i],ss[i]];",""),"
  }
  for(i in 1:Sum_zero){
    mu_z[i] = mu[nz[i],sz[i]];",
    switch(zi_param, "constant" = "",
           "covariate" = "
    zi_z[i] = zi[nz[i],sz[i]];"),
  ifelse(shp_param == "covariate" & family == "zi_neg_binomial", "
    kappa_z[i] = kappa[nz[i],sz[i]];",""),"
  }
  "),""),ifelse(shp_param == "covariate" & family != "zi_neg_binomial", paste("
  matrix[N,S]", switch(family, "gaussian" = "sigma",
                       "neg_binomial" = "kappa"),
                       "= exp(shp_X * shp_betas);",""),"")
  )
  # priors ####
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
      "gaussian" = switch(shp_param, "constant" = paste("
  //Standard deviation parameters
  sigma ~ ", prior[["sigma"]], ";
"), "covariate" = paste("
  //Standard deviation parameters
  to_vector(shp_betas) ~ ", prior[["shp_betas"]], ";
")),
      "neg_binomial" = switch(shp_param, "constant" = paste("
  //Scale parameter
  kappa ~ ", prior[["kappa"]], ";
"), "covariate" = paste("
  //Scale parameters
  to_vector(shp_betas) ~ ", prior[["shp_betas"]], ";
")),
      "bern" = "",
      "poisson" = "",
      "binomial" = "",
      "zi_poisson" = switch(zi_param,"constant" = paste("
  //zero-inflation parameter
  zi ~ ", prior[["zi"]], ";
"), "covariate" = paste("
  //zero-inflation parameter
  to_vector(zi_betas) ~ ", prior[["zi_betas"]], ";
")),
"zi_neg_binomial" = paste(switch(zi_param, "constant" = paste("
  //zero-inflation parameter
  zi ~ ", prior[["zi"]], ";
"), "covariate" = paste("
  //zero-inflation parameter
  to_vector(zi_betas) ~ ", prior[["zi_betas"]], ";
")), switch(shp_param, "constant" = paste("
  //Scale parameter
  kappa ~ ", prior[["kappa"]], ";
"), "covariate" = paste("
  //Scale parameter
  to_vector(shp_betas) ~ ", prior[["shp_betas"]], ";
"))
    )
  ))
  # model ####
  model_pt2 <- if(!grepl("zi_", family)){ paste(
    "
  for(i in 1:N) Y[i,] ~ ",
    switch(family,
      "gaussian" = switch(shp_param,"constant" = "normal(mu[i,], sigma);",
                          "covariate" = "normal(mu[i,],sigma[i,]);"),
      "bernoulli" = "bernoulli_logit(mu[i,]);",
      "neg_binomial" = switch(shp_param,
                              "constant" = "neg_binomial_2_log(mu[i,], kappa);",
                              "covariate" = "neg_binomial_2_log(mu[i,], kappa[i,]);"),
      "poisson" = "poisson_log(mu[i,]);",
      "binomial" = "binomial_logit(Ntrials[i], mu[i,]);"
    )
  )} else{paste("
  pos = 1;
  neg = 1;
  for(s in 1:S){
    target
      += N_zero[s]
           * log_sum_exp(",
           switch(zi_param,"constant" = "log(zi[s]),
                         log1m(zi[s])
                           +",
                  "covariate" = "bernoulli_logit_lpmf(1 | segment(zi_z, neg, N_zero[s])),
                         bernoulli_logit_lpmf(0 | segment(zi_z, neg, N_zero[s]))
                           +"),
  switch(family,
         "zi_poisson" = "poisson_log_lpmf(0 | segment(mu_z, neg, N_zero[s])));",
         "zi_neg_binomial" = paste0("neg_binomial_2_log_lpmf(0 | segment(mu_z, neg, N_zero[s]), ",
         switch(shp_param, "constant" = "kappa[s]));",
         "covariate" = "segment(kappa_z, neg, N_zero[s])));"))),"
    target += N_nonzero[s] * ",switch(zi_param,
    "constant" = "log1m(zi[s]);",
    "covariate" = "bernoulli_logit_lpmf(0 | segment(zi_nz, pos, N_nonzero[s]));"),"
    target +=",
    switch(family,
           "zi_poisson" = "poisson_log_lpmf(segment(Y_nz,pos,N_nonzero[s]) |
                                 segment(mu_nz, pos, N_nonzero[s]));",
           "zi_neg_binomial" = paste0("neg_binomial_2_log_lpmf(segment(Y_nz,pos,N_nonzero[s]) |
                                 segment(mu_nz, pos, N_nonzero[s]), ",
                                      switch(shp_param, "constant" = "kappa[s]);",
                                             "covariate" = "segment(kappa_nz, pos, N_nonzero[s]));"))),"
    pos = pos + N_nonzero[s];
    neg = neg + N_zero[s];
  }
")
  }
 # generated quantities ####
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
    matrix[N, S] linpred;",ifelse(grepl("zi", family) & zi_param == "covariate","
    matrix[N,S] zi = zi_X * zi_betas;",""),
      ifelse(shp_param == "covariate",paste("
    matrix[N,S]", switch(family, "gaussian" = "sigma",
    "neg_binomial" = "kappa",
    "zi_neg_binomial" = "kappa"),  "= exp(shp_X * shp_betas);"),""),
    switch(site_intercept, "ungrouped" = paste("
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
      for(j in 1:S) {",
      switch(family,
        "gaussian" = "log_lik[i, j] = normal_lpdf(Y[i, j] | linpred[i, j], sigma[j]);",
        "bernoulli" = "log_lik[i, j] = bernoulli_logit_lpmf(Y[i, j] | linpred[i, j]);",
        "neg_binomial" = "log_lik[i, j] = neg_binomial_2_log_lpmf(Y[i, j] | linpred[i, j], kappa[j]);",
        "poisson" = "log_lik[i, j] = poisson_log_lpmf(Y[i, j] | linpred[i, j]);",
        "binomial" = "log_lik[i, j] = binomial_logit_lpmf(Y[i, j] | Ntrials[i], linpred[i, j]);",
        "zi_poisson" = switch(zi_param,"constant" = "if (Y[i,j] == 0){
          log_lik[i, j] = log_sum_exp(bernoulli_lpmf(1 | zi[j]),
                          bernoulli_lpmf(0 |zi[j])
                          + poisson_log_lpmf(Y[i,j] | linpred[i,j]));
            } else {
              log_lik[i, j] = bernoulli_lpmf(0 | zi[j])
              + poisson_log_lpmf(Y[i,j] | linpred[i,j]);
            }",
                              "covariate" = "if (Y[i,j] == 0){
          log_lik[i, j] = log_sum_exp(bernoulli_logit_lpmf(1 | zi[i,j]),
                          bernoulli_logit_lpmf(0 |zi[i,j])
                          + poisson_log_lpmf(Y[i,j] | linpred[i,j]));
            } else {
              log_lik[i, j] = bernoulli_logit_lpmf(0 | zi[i,j])
              + poisson_log_lpmf(Y[i,j] | linpred[i,j]);
            }"),
        "zi_neg_binomial" = switch(zi_param,"constant" = paste0("if (Y[i,j] == 0){
          log_lik[i, j] = log_sum_exp(bernoulli_lpmf(1 | zi[j]),
                          bernoulli_lpmf(0 |zi[j])
                          + neg_binomial_2_log_lpmf(Y[i,j] | linpred[i,j],",
                          switch(shp_param, "constant" = "kappa[j]));",
                                 "covariate" = "kappa[,j]));"),"
            } else {
              log_lik[i, j] = bernoulli_lpmf(0 | zi[j])
              + neg_binomial_2_log_lpmf(Y[i,j] | linpred[i,j], kappa[j]);
            }"),
                                   "covariate" = paste0("if (Y[i,j] == 0){
          log_lik[i, j] = log_sum_exp(bernoulli_logit_lpmf(1 | zi[i,j]),
                          bernoulli_logit_lpmf(0 |zi[i,j])
                          + neg_binomial_2_log_lpmf(Y[i,j] | linpred[i,j],",
                          switch(shp_param, "constant" = "kappa[j]));",
                                 "covariate" = "kappa[,j]));"),"
            } else {
              log_lik[i, j] = bernoulli_logit_lpmf(0 | zi[i,j])
              + neg_binomial_2_log_lpmf(Y[i,j] | linpred[i,j],",
                          switch(shp_param, "constant" = "kappa[j]);",
                                 "covariate" = "kappa[,j]);"),"
            }"))
      ),"
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
