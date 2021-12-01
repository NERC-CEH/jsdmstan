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

  int<lower=0,upper=1> Y[N, S]; //Species matrix
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

  // model
  matrix[N, S] alpha = rep_matrix(a_bar + a * sigma_a, S);
  matrix[N, S] mu = alpha + (X * betas) + u;

#include /phylo/mglmm_model_priors.stan

  for(i in 1:N) Y[i,] ~ bernoulli_logit(mu[i,]);

}
generated quantities {
  // Calculate linear predictor, y_rep, log likelihoods for LOO
  matrix[N, S] log_lik;
  {
    matrix[N, S] linpred = rep_matrix(a_bar + a * sigma_a, S) + (X * betas) + u;

    for(i in 1:N) {
      for(j in 1:S) {
        log_lik[i, j] = bernoulli_logit_lpmf(Y[i, j] | linpred[i, j]);
      }
    }
  }

}
