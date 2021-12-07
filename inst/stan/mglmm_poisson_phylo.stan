functions{
#include /phylo/covariance_functions.stan
#include /include/tri_functions.stan
}
data {
  int<lower=1> N; // Number of samples
  int<lower=1> S; // Number of species
  int<lower=0> K; // Number of predictor variables
  matrix[N, K] X; // Predictor matrix

  matrix[S, S] Dmat; // Distance matrix between species
  int<lower=0,upper=3> nu05; // whether matern cov has nu = 1/2, 3/2 or 5/2
  real<lower=0> delta; // constant added to diagonal of covariance

  int<lower=0,upper=1> site_intercept; // whether to include a site intercept

  int<lower=0> Y[N, S]; //Species matrix
}
transformed data{
#include /mglmm/transformed_data.stan
}
parameters {
#include "/phylo/mglmm_pars.stan"

  // kernel parameters
  real<lower=0> etasq;
  real<lower=0> rho;


}
transformed parameters {
#include /phylo/mglmm_transformed_pars.stan
}
model {
  matrix[N, S] mu;
  // model
  if(site_intercept == 1){
    matrix[N, S] alpha = rep_matrix(a_bar[1] + a[1,] * sigma_a[1], S);
    mu = alpha + (X * betas) + u;
  } else{
    mu = (X * betas) + u;

  }

#include /phylo/mglmm_model_priors.stan


  for(i in 1:N) Y[i,] ~ poisson_log(mu[i,]);

}
generated quantities {
  // Calculate linear predictor, y_rep, log likelihoods for LOO
  matrix[N, S] log_lik;
  {
    matrix[N, S] linpred;
    if(site_intercept == 1){
      linpred = rep_matrix(a_bar[1] + a[1,] * sigma_a[1], S) + (X * betas) + u;
    } else{
      linpred = (X * betas) + u;
    }

    for(i in 1:N) {
      for(j in 1:S) {
        log_lik[i, j] = poisson_log_lpmf(Y[i, j] | linpred[i, j]);
      }
    }
  }

}
