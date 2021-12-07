  // Site-level intercept priors
  if(site_intercept == 1){
    a[1,] ~ std_normal();
    a_bar[1] ~ std_normal();
    sigma_a[1] ~ std_normal();
  }
  // Species parameter priors
  sigmas_b ~ std_normal();
  to_vector(z_preds) ~ std_normal();
  L_Rho_preds ~ lkj_corr_cholesky(1);

  // Species parameter priors
  sigmas_u ~ std_normal();
  to_vector(z_species) ~ std_normal();

  // kernel priors
  etasq ~ inv_gamma(10,0.1);
  rho ~ inv_gamma(10,0.1);


