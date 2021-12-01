/**
* Create covariance function based on matern
* @param x A distance matrix
* @param sq_eta
* @param rho
* @param delta Nugget added to diagonal
* @param nu05 Indicates what form the function should take
*/
matrix cov_matern(matrix x, real sq_eta, real rho, real delta, int nu05){
    int N = dims(x)[1];
    matrix[N,N] K;
    for (i in 1:(N-1)) {
      K[i,i] = sq_eta + delta;
      for (j in (i + 1):N) {
        if(nu05 == 0){
          K[i, j] = sq_eta * exp(- x[i,j]/rho );
        } else if(nu05 == 1){
          real dist_rho;
          dist_rho = x[i,j]/rho;
          K[i, j] = sq_eta*(1 + dist_rho)*exp(- dist_rho);
        } else if(nu05 == 2){
          real dist_rho;
          dist_rho = x[i,j]/rho;
          K[i, j] = sq_eta*(1 + dist_rho + pow(dist_rho, 2) / 3)*exp(- dist_rho);
        } else if(nu05 == 3){
          real dist_rho;
          K[i, j] = sq_eta * exp(- pow(dist_rho, 2) / 2);
          } else {
          reject("nu05 must be an integer in range 0-3, given as:",nu05);
        }
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_eta + delta;
    return K;
  }
