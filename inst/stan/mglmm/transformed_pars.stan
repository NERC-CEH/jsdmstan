
  matrix[K, S] betas;
  //vector[K_lower] L_Rho_preds_vector; //convert chol for post-processing
  matrix[N, S] u;
  //vector[S_lower] L_Rho_spec_vector; //convert chol for post-processing

  // covariance matrix on betas by preds
  betas = diag_pre_multiply(sigmas_b, L_Rho_preds) * z_preds;

  //L_Rho_preds_vector = lt_to_vector(L_Rho_preds);

  // species covariance matrix
  u = (diag_pre_multiply(sigmas_u, L_Rho_species) * z_species)';

  //L_Rho_spec_vector = lt_to_vector(L_Rho_species);

