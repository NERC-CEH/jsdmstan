#' Make stancode for the jsdm model
#'
#' This function returns the Stan code used to fit the model as specified by the
#' data list, family and method.
#'
#' @details Environmental covariate effects (\code{"betas"}) can be
#'   parameterised in two ways. With the \code{"cor"} parameterisation all
#'   covariate effects are assumed to be constrained by a correlation matrix
#'   between the covariates. With the \code{"unstruct"} parameterisation all
#'   covariate effects are assumed to draw from a simple distribution with no
#'   correlation structure. Both parameterisations can be modified using the
#'   prior object.
#'
#' @param method The method, one of \code{"gllvm"} or \code{"mglmm"}
#' @param family is the response family, must be one of \code{"gaussian"},
#'   \code{"neg_binomial"}, \code{"poisson"}, \code{"binomial"},
#'   \code{"bernoulli"}, \code{"zi_poisson"}, or \code{"zi_neg_binomial"}.
#'   Regular expression matching is supported.
#' @param prior The prior, given as the result of a call to [jsdm_prior()]
#' @param site_intercept Whether a site intercept should be included, potential
#'   values \code{"none"} (no site intercept), \code{"grouped"} (a site
#'   intercept with hierarchical grouping) or \code{"ungrouped"} (site intercept
#'   with no grouping)
#' @param beta_param The parameterisation of the environmental covariate
#'   effects, by default \code{"cor"}. See details for further information.
#' @param zi_param For the zero-inflated families, whether the zero-inflation
#'   parameter is a species-specific constant (default, \code{"constant"}), or
#'   varies by environmental covariates (\code{"covariate"}).
#' @param shp_param For the families with shape parameters, whether the shape
#'   parameter is a species-specific constant (default, \code{"constant"}), or
#'   varies by environmental covariates (\code{"covariate"}).
#' @param censoring If the response is left-censored (\code{"left"}) or not
#'   censored (default, \code{"none"}).
#' @param site_smooth If there is a smooth over a predictor that is constant
#'   across all species, default \code{"none"}.
#' @param species_smooth If there is a factor smooth over a predictor that
#'   varies by species, default \code{"none"}.
#'
#' @return A character vector of Stan code, class "jsdmstan_model"
#' @export
#'
#' @examples
#' jsdm_stancode(family = "gaussian", method = "gllvm")
#' jsdm_stancode(family = "poisson", method = "mglmm")
#'
jsdm_stancode <- function(method, family, prior = jsdm_prior(),
                          site_intercept = "none",
                          beta_param = "cor", zi_param = "constant",
                          shp_param = "constant", censoring = "none",
                          site_smooth = "none", species_smooth = "none") {
  # checks
  family <- match.arg(family, c("gaussian","normal", "bernoulli", "poisson",
                                "neg_binomial","binomial","zi_poisson",
                                "zi_neg_binomial","lognormal","gamma","beta"))
  family <- ifelse(family == "normal", "gaussian", family)
  method <- match.arg(method, c("gllvm", "mglmm"))
  beta_param <- match.arg(beta_param, c("cor","unstruct"))
  site_intercept <- match.arg(site_intercept, c("none","grouped","ungrouped"))
  zi_param <- match.arg(zi_param, c("constant","covariate"))
  shp_param <- match.arg(shp_param, c("constant","covariate"))
  censoring <- match.arg(censoring, c("none","left"))
  site_smooth <- match.arg(site_smooth, c("none", "single", "multiple"))
  species_smooth <- match.arg(species_smooth, c("none", "single", "multiple"))
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
    phylo = FALSE, prior = prior, site_intercept = site_intercept,
    beta_param = beta_param, zi_param = zi_param, shp_param = shp_param,
    censoring = censoring, site_smooth = site_smooth, species_smooth = species_smooth
  )
  class(scode) <- c("jsdmstan_model", "character")
  return(scode)
}


.modelcode <- function(method, family, phylo, prior, site_intercept,
                       beta_param, zi_param, shp_param, censoring,
                       site_smooth, species_smooth) {
  model_functions <- "
"
  # data ####
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
  array[N] int<lower=0, upper = ngrp> grps; // Vector matching sites to groups
 ",""),"
  array[N,S]",
  switch(family,
      "gaussian" = "real",
      "lognormal" = "real<lower=0>",
      "bernoulli" = "int<lower=0,upper=1>",
      "neg_binomial" = "int<lower=0>",
      "poisson" = "int<lower=0>",
      "zi_poisson" = "int<lower=0>",
      "zi_neg_binomial" = "int<lower=0>",
      "binomial" = "int<lower=0>",
      "gamma" = "real<lower=0>",
      "beta" = "real<lower=0,upper=1>",
    ), "Y; //Species matrix",
 ifelse(family == "binomial",
        "
  array[N] int<lower=0> Ntrials; // Number of trials",""),
  ifelse(grepl("zi_", family) & zi_param == "constant" & shp_param == "constant","
  array[S] int<lower=0> N_zero; // number of zeros per species
  array[S] int<lower=0> N_nonzero; //number of nonzeros per species
  int<lower=0> Sum_nonzero; //Total number of nonzeros across all species
  int<lower=0> Sum_zero; //Total number of zeros across all species
  array[Sum_nonzero] int<lower=0> Y_nz; //Y values for nonzeros
  array[Sum_nonzero] int<lower=0> ss; //species index for Y_nz
  array[Sum_nonzero] int<lower=0> nn; //site index for Y_nz
  array[Sum_zero] int<lower=0> sz; //species index for Y_z
  array[Sum_zero] int<lower=0> nz; //site index for Y_z",""),
ifelse(grepl("zi_", family) & zi_param == "covariate","
  int<lower=1> zi_k; //number of covariates for env effects on zi
  matrix[N, zi_k] zi_X; //environmental covariate matrix for zi",""),
ifelse(shp_param == "covariate","
  int<lower=1> shp_k; //number of covariates for env effects on family parameter
  matrix[N, shp_k] shp_X; //environmental covariate matrix for family parameter",""),
ifelse(censoring == "left", "
  array[S] int<lower=0> N_cens; // number of cens per species
  array[S] int<lower=0> N_noncens; //number of noncens per species
  array[N,S] int<lower=0> J_cens; //Indices of cens across all species
  array[N,S] int<lower=0> J_noncens; //Indices of noncens across all species",""),
.smooth_datacode("nfs", site_smooth),
.smooth_datacode("fs", species_smooth)
)

  # transformed data ####
  transformed_data <- ifelse(method == "gllvm", "
  // Ensures identifiability of the model - no rotation of factors
  int<lower=1> M;
  M = D * (S - D) + choose(D, 2) + D;
", "")


  # parameters ####
  site_inter_par <- switch(site_intercept,
                           "ungrouped" = "
  // Site intercepts
  real<lower=0> sigma_a;
  vector[N] z_a;",
  "none" = "",
  "grouped" = "
  // Site intercepts
  real<lower=0> sigma_a;
  vector[ngrp] z_a;")
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
  cholesky_factor_corr[S] cor_species_chol;"
  gllvm_pars <- "
  // Factor parameters
  vector[M] L; // Non-zero factor loadings
  real<lower=0> sigma_L; // variance of species loadings
  // Latent variables
  matrix[D, N] LV_uncor; // Per-site latent variable"

  var_pars <- switch(family,
    "gaussian" = switch(shp_param, "constant" = "
  real<lower=0> sigma; // Gaussian parameters", "covariate" = "
  matrix[shp_k,S] shp_betas; //environmental effects for family param"),
    "lognormal" = switch(shp_param, "constant" = "
  real<lower=0> sigma; // Lognormal parameters", "covariate" = "
  matrix[shp_k,S] shp_betas; //environmental effects for family param"),
    "gamma" = switch(shp_param, "constant" = "
  vector<lower=0>[S] shape; // Gaussian parameters", "covariate" = "
  matrix[shp_k,S] shp_betas; //environmental effects for family param"),
    "bernoulli" = "",
    "neg_binomial" = switch(shp_param, "constant" = "
  array[S] real<lower=0> kappa; // neg_binomial parameters","covariate" = "
  matrix[shp_k,S] shp_betas; //environmental effects for family param"),
    "poisson" = "",
    "zi_poisson" = switch(zi_param,
    "constant" = "
  array[S] real<lower=0,upper=1> zi; // zero-inflation parameter",
    "covariate" = "
  matrix[zi_k,S] zi_betas; //environmental effects for zi"),
  "zi_neg_binomial" = switch(zi_param,
                             "constant" = switch(shp_param, "constant" = "
  array[S] real<lower=0> kappa; // neg_binomial parameters
  array[S] real<lower=0,upper=1> zi; // zero-inflation parameter", "covariate" = "
  matrix[shp_k,S] shp_betas; //environmental effects for family param
  array[S] real<lower=0,upper=1> zi; // zero-inflation parameter"),
                             "covariate" = switch(shp_param, "constant" = "
  array[S] real<lower=0> kappa; // neg_binomial parameters
  matrix[zi_k,S] zi_betas; //environmental effects for zi", "covariate" = "
  matrix[shp_k,S] shp_betas; //environmental effects for family param
  matrix[zi_k,S] zi_betas; //environmental effects for zi"))
  )
  sitesmooth_pars <- ifelse(site_smooth=="none","","
  vector[nfs_ncoef] nfs_b; // smooth coefs
  vector<lower=1e-16>[nfs_nsp] nfs_sp; // smoothing parameters")
  speciessmooth_pars <- ifelse(species_smooth=="none","","
  vector[fs_ncoef] fs_b; // smooth coefs
  vector<lower=1e-16>[fs_nsp] fs_sp; // smoothing parameters")

  pars <- paste(
    site_inter_par, species_pars,
    switch(method,
      "gllvm" = gllvm_pars,
      "mglmm" = mglmm_spcov_pars
    ),
    var_pars,sitesmooth_pars,speciessmooth_pars
  )

  # transformed parameters ####
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
  matrix[S, S] cor_species;
  matrix[N, S] u;
  cor_species = multiply_lower_tri_self_transpose(cor_species_chol);
  u = (diag_pre_multiply(sigmas_species, cor_species) * z_species)';
  "
  ), if(beta_param == "cor") {"
  betas = diag_pre_multiply(sigmas_preds, cor_preds) * z_preds;
"} else {""})


  model_pt1 <- paste(if(species_smooth == "none"){""} else{"
  matrix[N,S] fs_Xb = to_matrix(fs_X * fs_b, N, S);"
    },if(site_smooth == "none"){""} else{"
  matrix[N,S] nfs_Xb = rep_matrix(nfs_X * nfs_b, S);"
    },
    switch(site_intercept, "none" = "",
                            "ungrouped" = "
  matrix[N, S] alpha = rep_matrix(z_a * sigma_a, S);",
                            "grouped" = "
  matrix[N, S] alpha;"),
  switch(method,
         "gllvm" = "
  // model
  matrix[N, S] LV_sum = ((Lambda_uncor * sigma_L) * LV_uncor)';",
         "mglmm" = ""),
  if(site_intercept == "grouped"){"
  {
    vector[ngrp] theta = z_a * sigma_a;
    for (n in 1:N){
      alpha[n,] = rep_row_vector(theta[grps[n]],S);
    }
  }"
  } else {""},"
  mu =", ifelse(site_intercept == "none","","alpha +"),"(X * betas) +",
  ifelse(site_smooth == "none","", "nfs_Xb +"),
  ifelse(species_smooth == "none","", "fs_Xb +"),
  switch(method, "gllvm" = "LV_sum;",
         "mglmm" = "u;"))

  # model ####
  model <- paste("
  matrix[N,S] mu;
  ", ifelse(grepl("zi_",family) & shp_param == "constant" & zi_param == "constant",paste0("
  real mu_nz[Sum_nonzero];
  real mu_z[Sum_zero];
  int pos;
  int neg;",
                   ifelse(shp_param== "covariate" & family == "zi_neg_binomial", "
  array[Sum_nonzero] real kappa_nz;
  array[Sum_zero] real kappa_z;","")),""),
  model_pt1,
  ifelse(grepl("zi_",family),paste0(ifelse(zi_param == "covariate", "
  matrix[N,S] zi = zi_X * zi_betas;",""),
  ifelse(shp_param == "covariate", "
  matrix[N,S] kappa = exp(shp_X * shp_betas);",""),
  ifelse(zi_param =="constant" & shp_param == "constant","
  for(i in 1:Sum_nonzero){
    mu_nz[i] = mu[nn[i],ss[i]];
  }
  for(i in 1:Sum_zero){
    mu_z[i] = mu[nz[i],sz[i]];
  }
  ","")),""),
  ifelse(shp_param == "covariate" & family != "zi_neg_binomial", paste("
  matrix[N,S]", switch(family, "gaussian" = "sigma",
                       "lognormal" = "sigma",
                       "gamma" = "shape",
                       "neg_binomial" = "kappa"),
                       "= exp(shp_X * shp_betas);",""),""),
  ifelse(family == "gamma", "
  mu = exp(mu);",""),

  .smooth_codechunk("nfs",site_smooth),
  .smooth_codechunk("fs",species_smooth)
  )
  # priors
  model_priors <- .prior_codechunk(site_intercept, beta_param, method, family,
                                   shp_param, zi_param, prior)
  # model end
  model_pt2 <- if(!grepl("zi_", family) & censoring == "none"){ paste(
    "
  for(i in 1:N) Y[i,] ~ ",
    switch(family,
      "gaussian" = switch(shp_param,"constant" = "normal(mu[i,], sigma);",
                          "covariate" = "normal(mu[i,],sigma[i,]);"),
      "lognormal" = switch(shp_param,"constant" = "lognormal(mu[i,], sigma);",
                          "covariate" = "lognormal(mu[i,],sigma[i,]);"),
      "bernoulli" = "bernoulli_logit(mu[i,]);",
      "neg_binomial" = switch(shp_param,
                              "constant" = "neg_binomial_2_log(mu[i,], kappa);",
                              "covariate" = "neg_binomial_2_log(mu[i,], kappa[i,]);"),
      "poisson" = "poisson_log(mu[i,]);",
      "binomial" = "binomial_logit(Ntrials[i], mu[i,]);",
      "gamma" = switch(shp_param,"constant" = "gamma(shape, shape ./ to_vector(mu[i,]));",
                           "covariate" = "gamma(shape[i,],shape[i,] ./ to_vector(mu[i,]));")
    )
  )} else if(censoring == "left"){
    .lcens_modelcode(family, shp_param)
  } else {
    .zi_familycode(family, zi_param, shp_param)
  }
 # generated quantities ####
  generated_quantities <- paste(
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
  }", "")
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

.smooth_codechunk <- function(id,amount){
  paste0(
  switch(amount, "none" = "",
         "single" = paste0("
  matrix[",id,"_nSr, ",id,"_nSr] ",id,"_K1;
  ",id,"_K1 = rep_matrix(0, ",id,"_nSr, ",id,"_nSr);
  vector[",id,"_nSr] ",id,"_zero;
  ",id,"_zero = rep_vector(0, ",id,"_nSr);
  {
    int runi=0;

    // prior for smooth
    for (i in 1:",id,"_nsp) {
      ",id,"_K1 += ",id,"_S1[1:",id,"_nSr, (runi+1):(runi+",id,"_nSr)] * ",id,"_sp[i];
      runi += ",id,"_nSr;
    }
  }
  ",id,"_b ~ multi_normal_prec(",id,"_zero, ", id,"_K1);"),
         "multiple" = paste0("
  // smoother bits
  {
    int off = 1;
    int beta_off = 0;
    // iterate over smooth terms in the model
    for(s in 1:",id,"_nterms){
      // initialize penalty matrix for this term to zero
      matrix[",id,"_Sr[s], ",id,"_Sr[s]] ",id,"_K1;
      ",id,"_K1 = rep_matrix(0, ",id,"_Sr[s], ",id,"_Sr[s]);
      vector[",id,"_Sr[s]] ",id,"_zero;
      ",id,"_zero = rep_vector(0, ",id,"_Sr[s]);

      // iterate over penalties within term
      for(i in 1:",id,"_Sn[s]){
        ",id,"_K1 += ",id,"_S1[1:",id,"_Sr[s],
                 (off + ",id,"_Sr[s]*(i-1)) : (off - 1 +(",id,"_Sr[s]*(i-1))+",id,"_Sr[s])] * ",id,"_sp[i];
      } // end iteration over penalties within term

      // sample from prior of the spline coefficients for this term
      ",id,"_b[(beta_off+1):(beta_off + ",id,"_Sr[s])] ~ multi_normal_prec(",id,"_zero, ",id,"_K1);

      off += ",id,"_Sr[s]*",id,"_Sn[s];
      beta_off += ",id,"_Sr[s];

    } // end iteration over smooth terms
  } // end block for local variable definition
")))
}

.lcens_modelcode <- function(family, shp_param){
  paste("
  for (s in 1:S){
    target += ",
        switch(family,
               "gaussian" = switch(shp_param,"constant" = "normal_lpdf(Y[J_noncens[1:N_noncens[s],s],s] |
                                mu[J_noncens[1:N_noncens[s],s],s], sigma);",
                                   "covariate" = "normal_lpdf(Y[J_noncens[1:N_noncens[s],s],s] |
                                mu[J_noncens[1:N_noncens[s],s],s],sigma[,s]);"),
               "lognormal" = switch(shp_param,"constant" = "lognormal_lpdf(Y[J_noncens[1:N_noncens[s],s],s]  |
                                mu[J_noncens[1:N_noncens[s],s],s], sigma);",
                                    "covariate" = "lognormal_lpdf(Y[J_noncens[1:N_noncens[s],s],s]  |
                                mu[J_noncens[1:N_noncens[s],s],s],sigma[,s]);"),
               "gamma" = switch(shp_param,"constant" = "gamma_lpdf(Y[J_noncens[1:N_noncens[s],s],s] | shape[s],
                               shape[s] /  mu[J_noncens[1:N_noncens[s],s],s]);",
                                "covariate" = "gamma_lpdf(Y[J_noncens[1:N_noncens[s],s],s]  | shape[s],
                               shape[,s] / mu[J_noncens[1:N_noncens[s],s],s]);")
        ),"
    target += ",
        switch(family,
               "gaussian" = switch(shp_param,"constant" = "normal_lcdf(Y[J_cens[1:N_cens[s],s],s] | mu[J_cens[1:N_cens[s],s],s], sigma);",
                                   "covariate" = "normal_lcdf(Y[J_cens[1:N_cens[s],s],s] | mu[J_cens[1:N_cens[s],s],s],sigma[n,s]);"),
               "lognormal" = switch(shp_param,"constant" = "lognormal_lcdf(Y[J_cens[1:N_cens[s],s],s] | mu[J_cens[1:N_cens[s],s],s], sigma);",
                                    "covariate" = "lognormal_lcdf(Y[J_cens[1:N_cens[s],s],s]| mu[J_cens[1:N_cens[s],s],s],sigma[n,s]);"),
               "gamma" = switch(shp_param,"constant" = "gamma_lcdf(Y[J_cens[1:N_cens[s],s],s] | shape[s], shape[s] / mu[J_cens[1:N_cens[s],s],s]);",
                                "covariate" = "gamma_lcdf(Y[J_cens[1:N_cens[s],s],s] | shape[n,s], shape[n,s] / mu[J_cens[1:N_cens[s],s],s]);")
        ),"

  }"
  )
}

.zi_familycode <- function(family, zi_param, shp_param){
  if(zi_param == "constant" & shp_param == "constant"){paste("
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
  } else{paste(
    "
    for(s in 1:S){
    for (n in 1:N){
      if(Y[n,s] == 0){
        target += log_sum_exp(bernoulli_logit_lpmf(1 |",
    switch(zi_param, "covariate" =  "zi[n,s])",
           "constant" = "zi[s])"),",
                         bernoulli_logit_lpmf(0 | ",
    switch(zi_param, "covariate" =  "zi[n,s])",
           "constant" = "zi[s])"),"
                           +",
    switch(family,
           "zi_poisson" = "poisson_log_lpmf(0 | mu[n,s]))",
           "zi_neg_binomial" = paste("neg_binomial_2_log_lpmf(0 | mu[n,s],",
                                     switch(shp_param, "constant" = "kappa[s]))",
                                            "covariate" = "kappa[n,s]))"))),";
      } else {
        target += bernoulli_logit_lpmf(0 | ",
    switch(zi_param, "covariate" =  "zi[n,s])",
           "constant" = "zi[s])"),"
                           + ",
    switch(family,
           "zi_poisson" = "poisson_log_lpmf(Y[n,s] | mu[n,s])",
           "zi_neg_binomial" = paste("neg_binomial_2_log_lpmf(Y[n,s] | mu[n,s],",
                                     switch(shp_param, "constant" = "kappa[s])",
                                            "covariate" = "kappa[n,s])"))),";
      }
    }
  }"
  )}
}

.prior_codechunk <- function(site_intercept, beta_param, method, family,
                             shp_param, zi_param, prior){
  paste(
    ifelse(site_intercept %in% c("ungrouped","grouped"), paste("
  // Site-level intercept priors
  z_a ~ std_normal();
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
  to_vector(z_species) ~ std_normal();
  cor_species_chol ~ ", prior[["cor_species_chol"]], ";
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
         "lognormal" = switch(shp_param, "constant" = paste("
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
         "gamma" = switch(shp_param, "constant" = paste("
  //Standard deviation parameters
  shape ~ ", prior[["shape"]], ";
"), "covariate" = paste("
  //Standard deviation parameters
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
}

.smooth_datacode <- function(id, amount){
  paste(
    if(amount == "none"){""} else{paste0("
  int<lower=0> ",id,"_nsp; // number of smoothing parameters
  int<lower=0> ",id,"_ncoef; // number of non-hyperparameters
  int<lower=0> ",id,"_nSr; // penalty dimension for smooth
  int<lower=0> ",id,"_nSc; // penalty dimension for smooth
  matrix[",ifelse(id == "nfs","N","N*S"),", ",id,"_ncoef] ",id,"_X;  // design matrix for smooth
  matrix[",id,"_nSr, ",id,"_nSc] ",id,"_S1;  // penalties for smooth, as one long matrix",
     if(amount == "multiple"){paste0("
  int<lower=0> ",id,"_nterms;  // number of terms in the smooth model
  array[",id,"_nterms] int<lower=0> ",id,"_Sr;
  array[",id,"_nterms] int<lower=0> ",id,"_Sn;")
       } else {""})
    }
  )
}
