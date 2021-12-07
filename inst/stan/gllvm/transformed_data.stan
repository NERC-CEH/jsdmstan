  //number of non-zero lower triangular chol factor entries (incl diagonal)
  //int<lower=1> K_lower;
  // Number of non-zero lower triangular factor loadings
  // Ensures identifiability of the model - no rotation of factors
  int<lower=1> M;
  M = D * (S - D) + D * (D - 1) %/% 2 + D;

  //K_lower = K * (K + 1) %/% 2;
