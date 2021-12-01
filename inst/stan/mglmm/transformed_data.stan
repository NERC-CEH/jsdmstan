  //number of non-zero lower triangular chol factor entries
  int<lower=1> S_lower;
  int<lower=1> K_lower;
  S_lower = S * (S + 1) / 2;
  K_lower = K * (K + 1) / 2;
