functions{
#include /include/tri_functions.stan
}
data {
  int<lower=1> N; // Number of samples
  int<lower=1> S; // Number of species
  int<lower=0> K; // Number of predictor variables
  matrix[N, K] X; // Predictor matrix

  int<lower=0,upper=1> site_intercept; // whether to include a site intercept

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
  matrix[N, S] mu;
  if(site_intercept == 1){
    matrix[N, S] alpha = rep_matrix(a_bar[1] + a[1,] * sigma_a[1], S);
    mu = alpha + (X * betas) + u;
  } else {
    mu = (X * betas) + u;
  }

#include /mglmm/model_priors.stan

  // NegBin scale parameter
  kappa ~ std_normal();

  for(i in 1:N) Y[i,] ~ neg_binomial_2_log(mu[i,], kappa);

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
        log_lik[i, j] = neg_binomial_2_log_lpmf(Y[i, j] | linpred[i, j], kappa);
      }
    }
  }

}
