#' Calculate covariance function based on a given data matrix
#'
#' Function for calculating the Matérn/ squared exponential covariance function given
#' a distance matrix
#'
#' @details This function takes a distance matrix and then calculates either the
#'   Matérn kernel with \eqn{\nu} 1/2, 3/2, 5/2, or infinite. Matérn kernel with
#'   infinite nu is equivalent to the squared exponential kernel, and with \eqn{\nu}
#'   = 1/2 the exponential kernel.
#'
#'   Matérn with \eqn{\nu = 1/2}, or exponential: \deqn{sq\_eta * \exp(-x_{ij}/rho)}
#'
#'   Matérn with \eqn{\nu = 3/2}: \deqn{sq\_eta * (1 + (\sqrt{3}x_{ij})/rho) * \exp(-\sqrt{3}x_{ij}/rho)}
#'
#'   Matérn with \eqn{\nu = 5/2}: \deqn{sq\_eta * (1 + \sqrt{5}x_{ij}/rho + 5(x_{ij}/rho)^2 /
#'   3) * \exp(-\sqrt{5}x_{ij}/rho)}
#'
#'   Matérn with \eqn{\nu = \infty}, or squared exponential: \deqn{sq\_eta *
#'   \exp(-(x_{ij}^2/2rho^2)}
#'
#'   For more information see the phylogenetic models vignette.
#'
#' @export
#'
#' @param x Distance matrix
#'
#' @param sq_eta Parameter for covariance function, see details
#'
#' @param rho Parameter for covariance function, see details
#'
#' @param delta = 1e-5 Nugget added to diagonal of resulting matrix to keep it
#'   positive definite
#'
#' @param nu05 Must be an integer in range 0-3. Indicates what type of covariance
#'   function is used. 0 is exponential, 1 is Matérn with nu = 1.5, 2 is Matérn with
#'   nu = 2.5 and 3 is the squared exponential.

cov_matern <- function(x, sq_eta, rho, delta = 1e-5, nu05) {
  N <- dim(x)[1]
  K <- matrix(ncol = N, nrow = N)
  for (i in 1:(N - 1)) {
    K[i, i] <- sq_eta + delta
    for (j in (i + 1):N) {
      if (nu05 == 0L) {
        K[i, j] <- sq_eta * exp(-x[i, j] / rho)
      } else if (nu05 == 1L) {
        K[i, j] <- sq_eta * (1 + (sqrt(3)*x[i, j]) / rho) * exp(-(sqrt(3)*x[i, j]) / rho)
      } else if (nu05 == 2L) {
        K[i, j] <- sq_eta * (1 + (sqrt(5)*x[i, j]) / rho + (5* (x[i, j]^2)) / (3 * (rho^2))) * exp(-(sqrt(5)*x[i, j]) / rho)
      } else if (nu05 == 3L) {
        K[i, j] <- sq_eta * exp(-((x[i, j] / rho)^2) / 2)
      } else {
        stop("nu05 must be integer in range 0-3, given as:", nu05)
      }
      K[j, i] <- K[i, j]
    }
  }
  K[N, N] <- sq_eta + delta
  return(K)
}

#' @describeIn cov_matern alias for fitting squared exponential
cov_sq_exp <- function(x, sq_eta, rho, delta = 1e-5) {
  cov_matern(x = x, sq_eta = sq_eta, rho = rho, delta = delta, nu05 = 3)
}
#' @describeIn cov_matern alias for fitting exponential
cov_exp <- function(x, sq_eta, rho, delta = 1e-5) {
  cov_matern(x = x, sq_eta = sq_eta, rho = rho, delta = delta, nu05 = 0)
}
