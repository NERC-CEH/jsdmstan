functions{
#include /include/tri_functions.stan
}
data {
  int<lower=1> N; // Number of samples
  int<lower=1> S; // Number of species
  int<lower=0> K; // Number of predictor variables
  matrix[N, K] X; // Predictor matrix

  int<lower=0> Y[N, S]; //Species matrix
}
transformed data{
#include /mglmm/transformed_data.stan
}
parameters {
#include "/mglmm/pars.stan"

  real<lower=0> kappa; // NegBin parameters
}
transformed parameters {
#include /mglmm/transformed_pars.stan
}
model {

  // model
  matrix[N, S] alpha = rep_matrix(a_bar + a * sigma_a, S);
  matrix[N, S] mu = alpha + (X * betas) + u;

#include /mglmm/model_priors.stan

  // NegBin scale parameter
  kappa ~ std_normal();

  for(i in 1:N) Y[i,] ~ neg_binomial_2_log(mu[i,], kappa);

}
generated quantities {
  // Calculate linear predictor, y_rep, log likelihoods for LOO
  matrix[N, S] log_lik;
  {
    matrix[N, S] linpred = rep_matrix(a_bar + a * sigma_a, S) + (X * betas) + u;

    for(i in 1:N) {
      for(j in 1:S) {
        log_lik[i, j] = neg_binomial_2_log_lpmf(Y[i, j] | linpred[i, j], kappa);
      }
    }
  }

}
