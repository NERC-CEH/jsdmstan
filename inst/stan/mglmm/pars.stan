  // Site intercepts
  real a_bar[site_intercept];
  real<lower=0> sigma_a[site_intercept];
  vector[N] a[site_intercept];

  // Species parameters

  //betas are hierarchical with covariance model
  vector<lower=0>[K] sigmas_b;
  matrix[K, S] z_preds;
  // covariance matrix on betas by predictors
  cholesky_factor_corr[K] L_Rho_preds;

  // species covariances
  vector<lower=0>[S] sigmas_u;
  matrix[S, N] z_species;
  cholesky_factor_corr[S] L_Rho_species;
