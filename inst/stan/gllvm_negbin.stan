functions{
#include /include/tri_functions.stan
}
data {
  int<lower=1> N; // Number of samples
  int<lower=1> S; // Number of species
  int<lower=1> D; // Number of latent dimensions
  int<lower=0> K; // Number of predictor variables
  matrix[N, K] X; // Predictor matrix

  int<lower=0> Y[N, S]; //Species matrix
}
transformed data{
#include "/gllvm/transformed_data.stan"
}
parameters {
#include /gllvm/pars.stan

  real<lower=0> kappa; // NegBin parameters
}
transformed parameters {
#include /gllvm/transformed_pars.stan
}
model {
  matrix[N, S] alpha = rep_matrix(a_bar + a * sigma_a, S);
  matrix[N, S] LV_sum = ((Lambda_uncor) * sigma_L * LV_uncor)';
  matrix[N, S] mu = alpha + (X * betas) + LV_sum;
#include /gllvm/model_priors.stan

  // NegBin scale parameter
  kappa ~ std_normal();

  // model

  for(i in 1:N) Y[i,] ~ neg_binomial_2_log(mu[i,], kappa);

}
generated quantities {
  // Calculate linear predictor, y_rep, log likelihoods for LOO
  matrix[N, S] log_lik;
#include /gllvm/sign_correct.stan
  {
    matrix[N, S] linpred = rep_matrix(a_bar + a * sigma_a, S) + (X * betas) +
                              ((Lambda_uncor) * sigma_L * LV_uncor)';

    for(i in 1:N) {
      for(j in 1:S) {
        log_lik[i, j] = neg_binomial_2_log_lpmf(Y[i, j] | linpred[i, j], kappa);
      }
    }
  }

}
