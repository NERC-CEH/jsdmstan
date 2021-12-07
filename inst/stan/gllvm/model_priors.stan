  // Site-level intercept priors
  if(site_intercept == 1){
    a[1,] ~ std_normal();
    a_bar[1] ~ std_normal();
    sigma_a[1] ~ std_normal();
  }

  // Species parameter priors
  sigmas_b ~ std_normal();
  to_vector(z_species) ~ std_normal();
  // covariance matrix priors
  L_Rho_preds ~ lkj_corr_cholesky(1);


  // Factor priors
  to_vector(LV_uncor) ~ std_normal();
  L ~ std_normal();
  sigma_L ~ std_normal(); // Variance of factor loadings
