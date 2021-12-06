// Sign correct factor loadings and factors
  matrix[D, N] LV;
  matrix[S, D] Lambda;
  //matrix[S, S] COV;
  //vector[M] Lambda_vect;
  for(d in 1:D){
    if(Lambda_uncor[d,d] < 0){
      Lambda[,d] = -1 * Lambda_uncor[,d];
      LV[d,] = -1 * LV_uncor[d,];
    } else {
      Lambda[,d] = Lambda_uncor[,d];
      LV[d,] = LV_uncor[d,];
    }
  }
  //Lambda_vect = lt_to_vector(Lambda);

  // Calculate species covariance matrix
  //COV = multiply_lower_tri_self_transpose(Lambda);
