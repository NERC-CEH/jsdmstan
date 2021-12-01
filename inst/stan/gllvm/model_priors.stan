  // Site-level intercept priors
  a ~ std_normal();
  a_bar ~ std_normal();
  sigma_a ~ std_normal();

  // Species parameter priors
  sigmas_b ~ std_normal();
  to_vector(z_species) ~ std_normal();
  // covariance matrix priors
  L_Rho_preds ~ lkj_corr_cholesky(1);


  // Factor priors
  to_vector(LV_uncor) ~ std_normal();
  L ~ std_normal();
  sigma_L ~ std_normal(); // Variance of factor loadings
