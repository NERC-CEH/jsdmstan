  // covariance matrix on betas by preds
  matrix[K, S] betas;
  vector[K_lower] L_Rho_preds_vector;


  // Construct factor loading matrix
  matrix[S, D] Lambda_uncor;
  // Constraints to allow identifiability of loadings
  for (i in 1:(D-1)) {
    for (j in (i+1):(D)){
      Lambda_uncor[i,j] = 0;
    }
  }
  {
    int index;
    index = 0;
    for (j in 1:D) {
      for (i in j:S) {
        index = index + 1;
        Lambda_uncor[i, j] = L[index];
      }
    }
  }


  betas = diag_pre_multiply(sigmas_b, L_Rho_preds) * z_species;

  L_Rho_preds_vector = lt_to_vector(L_Rho_preds);

