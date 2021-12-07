functions{
#include /include/tri_functions.stan
}
data {
  int<lower=1> N; // Number of samples
  int<lower=1> S; // Number of species
  int<lower=1> D; // Number of latent dimensions
  int<lower=0> K; // Number of predictor variables
  matrix[N, K] X; // Predictor matrix

  int<lower=0,upper=1> site_intercept; // whether to include a site intercept

  int<lower=0,upper=1> Y[N, S]; //Species matrix
}
transformed data{
#include "/gllvm/transformed_data.stan"
}
parameters {
#include "/gllvm/pars.stan"
}
transformed parameters {
#include /gllvm/transformed_pars.stan
}
model {
  matrix[N,S] mu;

  // model
  matrix[N, S] LV_sum = ((Lambda_uncor) * sigma_L * LV_uncor)';
  if(site_intercept == 1){
    matrix[N, S] alpha = rep_matrix(a_bar[1] + a[1,] * sigma_a[1], S);
    mu = alpha + (X * betas) + LV_sum;
  } else{
      mu = (X * betas) + LV_sum;
  }
#include /gllvm/model_priors.stan

  for(i in 1:N) Y[i,] ~ bernoulli_logit(mu[i,]);

}
generated quantities {

  // Calculate linear predictor, y_rep, log likelihoods for LOO
  matrix[N, S] log_lik;
#include /gllvm/sign_correct.stan
  {
    matrix[N, S] linpred;
    if(site_intercept == 1){
      linpred = rep_matrix(a_bar[1] + a[1,] * sigma_a[1], S) + (X * betas) +
                             ((Lambda_uncor) * sigma_L * LV_uncor)';
    } else{
      linpred = (X * betas) + ((Lambda_uncor) * sigma_L * LV_uncor)';
    }

    for(i in 1:N) {
      for(j in 1:S) {
        log_lik[i, j] = bernoulli_logit_lpmf(Y[i, j] | linpred[i, j]);
      }
    }
  }

}
