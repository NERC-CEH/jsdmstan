
#' Make stancode for the jsdm model
#'
#' This function returns the Stan code used to fit the model as specified by the data
#' list, family and method.
#'
#' @param data_list The data list
#' @param family The family, one of "\code{"gaussian"}, \code{"bernoulli"},
#'   \code{"poisson"} or \code{"neg_binomial"}
#' @param prior The prior, currently only default (i.e. \code{NULL}) supported
#' @param method The method, one of \code{"gllvm"} or \code{"mglmm"}
#' @param phylo Whether phylo should be included, only for MGLMM
#'
#' @return A character vector of Stan code
#' @export
#'
jsdm_stancode <- function(data_list, family, method, prior = NULL,
                          phylo = FALSE){
  # checks
  family <- match.arg(family, c("gaussian","bernoulli","poisson","neg_binomial"))
  method <- match.arg(method, c("gllvm","mglmm"))

  if(!all(c("Y","K","S","N","X","site_intercept") %in% names(data_list)))
    stop("Data list must have entries Y, K, S, N, X and site_intercept")


  # data processing steps

  scode <- .modelcode(method = method, family = family,
                      phylo = phylo, prior = prior)
  class(scode) <- c("jsdmstan_model","character")
  return(scode)

}


.modelcode <- function(method, family, phylo = NULL, prior = NULL){
  model_functions <- "
  matrix to_lower_tri(vector x, int nr, int nc){
    matrix[nr,nc] y;
    int pos = 1;
    for(i in 1:nr){
      for(j in 1:nc){
        if(i < j){
          y[i,j] = 0;
        } else{
          y[i,j] = x[pos];
          pos += 1;
        }
      }
    }
    return y;
  }

  vector lt_to_vector(matrix y){
    int nr = rows(y);
    int nc = cols(y);
    vector[nc*(nr-nc) + min(rows(y),cols(y)) * (min(rows(y),cols(y)) + 1) / 2] x;
    int pos = 1;
    // now fill vector with matrix elements
    for(i in 1:nr){
      for(j in 1:nc){
        if(i < j) continue;

        x[pos] = y[i,j];
        pos += 1;

      }
    }
    return x;
  }

"
  data <- paste("
  int<lower=1> N; // Number of samples
  int<lower=1> S; // Number of species
  ", ifelse(method == "gllvm",
            "int<lower=1> D; // Number of latent dimensions",""),
  "
  int<lower=0> K; // Number of predictor variables
  matrix[N, K] X; // Predictor matrix

  int<lower=0,upper=1> site_intercept; // whether to include a site intercept
  ",
  switch(family,"gaussian" = "real ",
         "bernoulli" = "int<lower=0,upper=1> ",
         "neg_binomial" = "int<lower=0> ",
         "poisson" = "int<lower=0> "), "Y[N,S]; //Species matrix")
  transformed_data <- ifelse(method == "gllvm","
  // Ensures identifiability of the model - no rotation of factors
  int<lower=1> M;
  M = D * (S - D) + D * (D - 1) / 2 + D;
","")



  site_inter_par <- "
  // Site intercepts
  real a_bar[site_intercept];
  real<lower=0> sigma_a[site_intercept];
  vector[N] a[site_intercept];"
  species_pars <- "
  //betas are hierarchical with covariance model
  vector<lower=0>[K] sigmas_b;
  matrix[K, S] z_preds;
  // covariance matrix on betas by predictors
  cholesky_factor_corr[K] L_Rho_preds;"
  mglmm_spcov_pars <- "
  // species covariances
  vector<lower=0>[S] sigmas_u;
  matrix[S, N] z_species;
  cholesky_factor_corr[S] L_Rho_species;"
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
  "poisson" = "")

  pars <- paste(site_inter_par, species_pars,
                switch(method, "gllvm" = gllvm_pars,
                       "mglmm" = mglmm_spcov_pars),
                       var_pars)
  transformed_pars <- paste("
  // covariance matrix on betas by preds
  matrix[K, S] betas;
  ", switch(method,
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
  u = (diag_pre_multiply(sigmas_u, L_Rho_species) * z_species)';
  "),"
  betas = diag_pre_multiply(sigmas_b, L_Rho_preds) * z_preds;
")



  gllvm_model <- "
  // model
  matrix[N, S] LV_sum = ((Lambda_uncor * sigma_L) * LV_uncor)';
  if(site_intercept == 1){
    matrix[N, S] alpha = rep_matrix(a_bar[1] + a[1,] * sigma_a[1], S);
    mu = alpha + (X * betas) + LV_sum;
  } else{
      mu = (X * betas) + LV_sum;
  }
  "
  mglmm_model <- "
  // model
  if(site_intercept == 1){
    matrix[N, S] alpha = rep_matrix(a_bar[1] + a[1,] * sigma_a[1], S);
    mu = alpha + (X * betas) + u;
  } else {
    mu = (X * betas) + u;
  }
  "
  model <- paste("
  matrix[N,S] mu;
  ", switch(method,
            "gllvm" = gllvm_model,
            "mglmm" = mglmm_model))
  if(is.null(prior)){
    model_priors <- paste("
  // Site-level intercept priors
  if(site_intercept == 1){
    a[1,] ~ std_normal();
    a_bar[1] ~ std_normal();
    sigma_a[1] ~ std_normal();
  }

  // Species parameter priors
  sigmas_b ~ std_normal();
  to_vector(z_preds) ~ std_normal();
  // covariance matrix priors
  L_Rho_preds ~ lkj_corr_cholesky(1);
",
switch(method,"gllvm" = "
  // Factor priors
  to_vector(LV_uncor) ~ std_normal();
  L ~ std_normal();
  sigma_L ~ std_normal(); // Variance of factor loadings
", "mglmm" = "
  // Species parameter priors
  sigmas_u ~ std_normal();
  to_vector(z_species) ~ std_normal();
  L_Rho_species ~ lkj_corr_cholesky(1);
"),
switch(family, "gaussian" = "
  //Standard deviation parameters
  sigma ~ std_normal();
", "neg_binomial" = "
  //Scale parameter
  kappa ~ std_normal();
", "bern" = "", "poisson" = ""))} else{
  stop("Specifying own priors currently unsupported")
}

  generated_quantities <- paste("
  // Calculate linear predictor, y_rep, log likelihoods for LOO
  matrix[N, S] log_lik;
  ",
  ifelse(method == "gllvm","
  // Sign correct factor loadings and factors
  matrix[D, N] LV;
  matrix[S, D] Lambda;
  //matrix[S, S] COV;
  //vector[M] Lambda_vect;
  for(d in 1:D){
    if(Lambda_uncor[d,d] < 0){
      Lambda[,d] = -1 * Lambda_uncor[,d];
      LV[d,] = -1 * LV_uncor[d,];
    } else {
      Lambda[,d] = Lambda_uncor[,d];
      LV[d,] = LV_uncor[d,];
    }
  }",""),"
  {
    matrix[N, S] linpred;
    if(site_intercept == 1){
      linpred = rep_matrix(a_bar[1] + a[1,] * sigma_a[1], S) + (X * betas) +",
  switch(method,"gllvm" = "((Lambda_uncor * sigma_L) * LV_uncor)'",
         "mglmm" = "u"),";
    } else{
      linpred = (X * betas) +",
  switch(method,"gllvm" = "((Lambda_uncor * sigma_L) * LV_uncor)'",
         "mglmm" = "u"),";
    }
        for(i in 1:N) {
      for(j in 1:S) {
        log_lik[i, j] = ",
  switch(family,
         "gaussian" = "normal_lpdf",
         "bernoulli" = "bernoulli_logit_lpmf",
         "neg_binomial" = "neg_binomial_2_log_lpmf",
         "poisson" = "poisson_log_lpmf"),
  "(Y[i, j] | linpred[i, j]",
  switch(family,
         "gaussian" = ", sigma)",
         "bernoulli" = ")",
         "neg_binomial" = ", kappa)",
         "poisson" = ")"), ";
      }
    }
  }
  ")


  if(isTRUE(phylo)){
    model_functions <- paste(model_functions,
                             '
/**
* Create covariance function based on matern
* @param x A distance matrix
* @param sq_eta
* @param rho
* @param delta Nugget added to diagonal
* @param nu05 Indicates what form the function should take
*/
matrix cov_matern(matrix x, real sq_eta, real rho, real delta, int nu05){
    int N = dims(x)[1];
    matrix[N,N] K;
    for (i in 1:(N-1)) {
      K[i,i] = sq_eta + delta;
      for (j in (i + 1):N) {
        if(nu05 == 0){
          K[i, j] = sq_eta * exp(- x[i,j]/rho );
        } else if(nu05 == 1){
          real dist_rho;
          dist_rho = x[i,j]/rho;
          K[i, j] = sq_eta*(1 + dist_rho)*exp(- dist_rho);
        } else if(nu05 == 2){
          real dist_rho;
          dist_rho = x[i,j]/rho;
          K[i, j] = sq_eta*(1 + dist_rho + pow(dist_rho, 2) / 3)*exp(- dist_rho);
        } else if(nu05 == 3){
          real dist_rho;
          K[i, j] = sq_eta * exp(- pow(dist_rho, 2) / 2);
          } else {
          reject("nu05 must be an integer in range 0-3, given as:",nu05);
        }
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_eta + delta;
    return K;
  }
'
    )
    data <- paste(data,"
  matrix[S, S] Dmat; // Distance matrix between species
  int<lower=0,upper=3> nu05; // whether matern cov has nu = 1/2, 3/2 or 5/2
  real<lower=0> delta; // constant added to diagonal of covariance
")
    pars <- paste(gsub("cholesky_factor_corr\\[S\\] L_Rho_species;","",pars),"
  // kernel parameters
  real<lower=0> etasq;
  real<lower=0> rho;
")
    transformed_pars <- gsub(
      "u = ",
      "matrix[S,S] L_Rho_species;
   L_Rho_species = cholesky_decompose(cov_matern(Dmat, etasq, rho, delta, nu05));
   u = ",
      transformed_pars)
  }

  res <- paste("//Generated by jsdmstan\n",
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
               model,"\n",model_priors,"\n",
               "\n}\ngenerated quantities{\n",
               generated_quantities,
               "\n}\n\n")

  return(res)

}

print.jsdmstan_model <- function(x){
  cat(x)
}
