  // Site-level intercept priors
  a ~ std_normal();
  a_bar ~ std_normal();
  sigma_a ~ std_normal();

  // Species parameter priors
  sigmas_b ~ std_normal();
  to_vector(z_preds) ~ std_normal();
  L_Rho_preds ~ lkj_corr_cholesky(1);

  // Species parameter priors
  sigmas_u ~ std_normal();
  to_vector(z_species) ~ std_normal();
  L_Rho_species ~ lkj_corr_cholesky(1);
