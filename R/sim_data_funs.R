#'Generate simulated data within a variety of jSDM methodologies
#'
#'The \code{jsdm_sim_data} function can simulate data with either a multivariate
#'generalised mixed model (MGLMM) or a generalised linear latent variable model
#'(GLLVM). The \code{gllvm_sim_data} and \code{mglmm_sim_data} are aliases for
#'\code{jsdm_sim_data} that set \code{method} to \code{"gllvm"} and \code{"mglmm"}
#'respectively.
#'
#'@details This simulates data based on a joint species distribution model with
#'  either a generalised linear latent variable model approach or a multivariate
#'  generalised linear mixed model approach.
#'
#'  Models can be fit with or without "measured predictors", and if measured
#'  predictors are included then the species have species-specific parameter
#'  estimates. These can either be simulated completely independently, or have
#'  information pooled across species. If information is pooled this can be modelled
#'  as either a random draw from some mean and standard deviation or species
#'  covariance can be modelled together (this will be the covariance used in the
#'  overall model if the method used has covariance).
#'
#'  It makes a whole load of assumptions, incl that all parameters and measured
#'  predictors are normally distributed with mean 0 and standard deviation 1, with
#'  variance parameters being the absolute of this distribution. The exceptions are
#'  that if there is a Matern (or exponential quadratic) kernel being fit the etasq
#'  and rho parameters are either simulated as an inverse gamma with shape 10 and
#'  scale 0.1 if the invgamma package is installed or as an absolute value of a
#'  normal distribution with mean 1 and standard deviation 0.2. If no phylogenetic
#'  matrix is supplied the correlations between species are simulated as a random
#'  draw from a LKJ correlation matrix with eta = 1.
#'
#'@export
#'
#'@param N is number of sites
#'
#'@param S is number of species
#'
#'@param D is number of latent variables, used within gllvm method
#'
#'@param K is number of covariates, by default 0
#'
#'@param family is the response family, must be one of "gaussian", "neg_binomial",
#'  "poisson" or "bernoulli"
#'
#'@param method is the jSDM method to use, currently either GLLVM or MGLMM - see
#'  details for more information.
#'
#'@param phylo is whether to randomly generate a phylogenetic tree that will be used
#'  to constrain beta estimates. Can be given as TRUE, FALSE or as a given distance
#'  matrix.
#'
#'@param species_intercept Whether to include an intercept in the predictors, must be
#'  TRUE if K is 0. Defaults to TRUE.
#'
#'@param site_intercept Whether to include a site intercept. Defaults to FALSE.
#'
#'@param delta Nugget added to diagonal of resulting matrix to keep it positive
#'  definite
#'
#'@param nu05 Must be an integer in range 0-3. Indicates what type of covariance
#'  function is used. 0 is exponential, 1 is Matérn with nu = 1.5, 2 is Matérn with
#'  nu = 2.5 and 3 is the squared exponential.
#'
#'@param eta Shape parameter in random generation of LKJ matrix, defaults to 1.

jsdm_sim_data <- function(N, S, D = NULL, K = 0L, family, method = c("gllvm","mglmm"),
                          phylo = FALSE, species_intercept = TRUE,
                          site_intercept = FALSE, delta = NULL, nu05 = NULL,
                          eta = 1){
  response <- match.arg(family, c("gaussian","neg_binomial","poisson","bernoulli"))

  if(any(!c(is.double(N), is.double(S))))
    stop("N and S must be positive integers")

  if(any(c(N,S) <= 0) | any(c(N,S)%%1 != 0))
    stop("N and S must be positive integers")

  if(K < 0 | K%%1 != 0)
    stop("K must be either 0 or a positive integer")
  if(K == 0 & isFALSE(species_intercept))
    stop("If K is 0 then a species intercept is required")

  if (isTRUE(phylo) & method == "gllvm") {
    stop("Phylogenetic sharing of information not currently supported for GLLVM")
  }
  if (isTRUE(phylo) & !requireNamespace("ape", quietly = TRUE)) {
    stop("The ape package is required for random generation of phylogenetic trees")
  }
  if (isTRUE(phylo) & any(c(is.null(delta), is.null(nu05)))) {
    stop("Need to specify delta and nu05 arguments for phylo")
  }

  if(method == "gllvm" & !is.double(D)){
    stop("gllvm method require D to be a positive integer")
  }


  # build phylogenetic tree - need to add check
  # if(phylo & )
  if(isTRUE(phylo)){

    if(requireNamespace("invgamma",quietly = TRUE)){
      sq_eta <- invgamma::rinvgamma(1,10,scale = 0.1)
      rho <- invgamma::rinvgamma(1,10,scale = 0.1)
    } else{
      sq_eta <- abs(stats::rnorm(1,1,0.2))
      rho <- abs(stats::rnorm(1,1,0.2))
    }

    if(isTRUE(phylo)){
      Dmat <- stats::cophenetic(ape::rtree(S))
    } else if(is.matrix(phylo)){
      if(all(dim(phylo)==S) & isSymmetric(phylo)){
        Dmat <- phylo
      } else stop("phylo matrix must be a symmetric square matrix of with row number equal to number of species")
    }


    # return matern covariance
    L_Rho_species <- t(chol(cov_matern(Dmat, sq_eta = sq_eta, rho = rho, delta = delta,
                                     nu05 = nu05)))
  } else if(isFALSE(phylo) & method=="mglmm"){
    L_Rho_species <- rlkj(S, cholesky = TRUE, eta = eta)
  }

  # now do covariates - if K = NULL then do intercept only
  if(K==0){
    x <- matrix(1, nrow = N, ncol = 1)
    colnames(x) <- "Intercept"
    J <- 1
  } else if(isTRUE(species_intercept)){
    x <- matrix(stats::rnorm(N*K), ncol = K, nrow = N)
    colnames(x) <- paste0("V",1:K)
    x <- cbind(Intercept = 1, x)
    J <- K + 1
  } else if(isFALSE(species_intercept)){
    x <- matrix(stats::rnorm(N*K), ncol = K, nrow = N)
    colnames(x) <- paste0("V",1:K)
    J <- K
  } else stop("K must be an integer value")

  # covariate parameters
  beta_sds <- abs(stats::rnorm(J))
  z_betas <- matrix(stats::rnorm(S*J), ncol = S, nrow = J)
  if(K == 0){
    beta_sim <- beta_sds %*% z_betas
    } else{
    L_Rho_preds <- rlkj(J, cholesky = TRUE, eta = eta)
    beta_sim <- (diag(beta_sds) %*% L_Rho_preds) %*% z_betas
  }

  mu_sim <- x %*% beta_sim


  ## site intercept
  if(isTRUE(site_intercept)){
    a_bar <- stats::rnorm(1)
    sigma_a <- abs(stats::rnorm(1))
    a <- stats::rnorm(N)
    a_i <- a_bar + a * sigma_a
  } else{
    a_i <- rep(0, N)
  }

  if(method == "mglmm"){
    # u covariance
    u_sds <- abs(stats::rnorm(S))
    u_ftilde <- matrix(stats::rnorm(S*N), nrow = S, ncol = N)
    u_ij <- t((diag(u_sds) %*% L_Rho_species) %*% u_ftilde)
  }

  if(method == "gllvm"){
    M <- D * (S - D) + D * (D - 1) / 2

    L_l <- stats::rnorm(M, 0, 1)
    L_d <- abs(stats::rnorm(D, 0, 1))
    L <- matrix(nrow = S, ncol = D)

    idx2 = 0;
    for (i in 1:(D-1)) { for (j in (i+1):(D)){ L[i,j] = 0 } }
    diag(L) <- L_d
    for (j in 1:D) {
      for (i in (j+1):S) {
        idx2 = idx2+1
        L[i,j] = L_l[idx2]
      }
    }

    L_sigma <- abs(stats::rnorm(1,0,1))

    LV <- matrix(stats::rnorm(N * D, 0, 1), nrow = D, ncol = N)

    LV_sum <- (L * L_sigma) %*% LV
  }

  # variance parameters
  if(response %in% c("neg_binomial", "gaussian"))
    sigma <- abs(stats::rnorm(1))

  Y <-  matrix(nrow = N, ncol = S)
  for(i in 1:N) {
    for(j in 1:S){

      mu_ij <- switch(method,
                      "gllvm" = a_i[i] + mu_sim[i,j] + LV_sum[j,i],
                      "mglmm" = a_i[i] + mu_sim[i,j] + u_ij[i,j])

      Y[i, j] <- switch(response,
                        "neg_binomial" = rgampois(1, mu = exp(mu_ij),
                                                  scale = sigma),
                        "gaussian" = stats::rnorm(1, mu_ij, sigma),
                        "poisson" = stats::rpois(1, exp(mu_ij)),
                        "bernoulli" = stats::rbinom(1, 1, inv_logit(mu_ij)))
    }
  }


  pars =  list(betas = beta_sim,
               beta_sds = beta_sds,
               z_betas = z_betas)
  if(isTRUE(phylo) | is.matrix(phylo)){
    pars$etasq <- sq_eta
    pars$rho <- rho
  }

  if(isTRUE(site_intercept)){
    pars$a_bar <- a_bar
    pars$sigma_a <- sigma_a
    pars$a <- a

  }
  if(method == "gllvm"){
    pars$L <- L
    pars$LV <- LV
    pars$L_sigma <- L_sigma
  }
  if(method == "mglmm"){
    pars$u_sds <- u_sds
    pars$L_Rho_species <- L_Rho_species
    pars$u_ftilde <- u_ftilde
  }
  if(isTRUE(species_intercept)){
    if(K > 0){
      x <- x[,2:ncol(x)]
    } else{
      x <- NULL
    }
  }
  output <- list(Y = Y, pars = pars, N = N, S = S, D = D, K = J, X = x,
                 site_intercept = as.integer(site_intercept))
  if(isTRUE(phylo)){
    output$Dmat <- Dmat
    output$delta <- delta
    output$nu05 <- nu05
  } else if(is.matrix(phylo)){
    output$Dmat <- phylo
    output$delta <- delta
    output$nu05 <- nu05
  }

  return(output)

}

#'
#'
#' @export
#' @describeIn jsdm_sim_data Alias for \code{jsdm_sim_data} with \code{method =
#'  "gllvm"}
#'
#'@param ... Arguments passed to jsdm_sim_data
gllvm_sim_data <- function(...){
  jsdm_sim_data(method = "gllvm", ...)
}

#'
#'@export
#'@describeIn jsdm_sim_data Alias for \code{jsdm_sim_data} with \code{method =
#'  "mglmm"}
#'
#'@param ... Arguments passed to jsdm_sim_data
mglmm_sim_data <- function(...){
  jsdm_sim_data(method = "mglmm", ...)
}

#' Helper functions for simulating data
#'
#' The \code{rlkj} function is for generating random LKJ correlation matrices and the
#' \code{rgampois} function generates random draws from the Stan's alternative
#' parameterisation of the negative binomial distribution.
#'
#' The Lewandowski-Kurowicka-Joe (LKJ) distribution is a prior distribution for
#' correlation matrices, with the shape parameter eta. If eta is 1 then the density
#' is uniform over the correlation matrix, ith eta > 1 then the the probability
#' concentrates around the identity matrix while is 0 < eta < 1 the probability
#' concentrates away from the identity matrix.
#'
#' The alternative parameterisation of the negative binomial distribution is:
#'
#' \deqn{NegBinomial2(y | mu, scale) = binom(y+scale-1,y) (mu/mu+scale)^y (scale/mu + scale)^scale}
#'
#' Where the mean of the distribution is mu and the variance is \eqn{mu + (mu^2/scale)}
#'
#' The \code{rlkj} function are sourced from Ben Goodrich's response on the Stan
#' google mailing list. (see link
#' \url{https://groups.google.com/g/stan-users/c/3gDvAs_qwN8/m/Xpgi2rPlx68J)}). The
#' `rgampois` function is sourced from the rethinking package by Richard McElreath.
#'
#' @param n The number of samples to create/dimension of correlation matrix
#' @param mu The mean of the negative binomial parameterisation
#' @param scale The phi parameter that controls overdispersion of the negative
#'   binomial distribution (see details for description)
#' @param eta The shape parameter of the LKJ distribution
#' @param cholesky Whether the correlation matrix should be returned as the Cholesky
#'   decomposition, by default \code{FALSE}
#'
#' @name sim_helpers
NULL

#' @rdname sim_helpers
#' @export
rgampois <- function(n, mu, scale) {
  # shape <- mu/scale
  prob <- scale/(mu + scale)
  stats::rnbinom(n, size = scale, prob = prob)
}

inv_logit <- function(x){
  1/(1 + exp(x))
}


# The following is code from Ben Goodrich's response on the stan google mailing list
# (see link https://groups.google.com/g/stan-users/c/3gDvAs_qwN8/m/Xpgi2rPlx68J)

rgbeta <-
  function(n, shape) {
    if(shape == Inf)     rep(0, n)
    else if(shape > 0)  -1 + 2 * stats::rbeta(n, shape, shape)
    else if(shape == 0) -1 + 2 * stats::rbinom(n, 1, 0.5)
    else stop("shape must be non-negative")
  }

#' @rdname sim_helpers
#' @export
rlkj <-
  function(n, eta = 1, cholesky = FALSE) {
    if (n < 2){
      stop("Dimension of correlation matrix must be >= 2")
    }
    if (eta < 1){
      stop("The value of eta must be >= 1")
    }
    alpha <- eta + (n - 2) / 2
    L <- matrix(0, n, n)
    L[1,1] <- 1
    L[-1,1] <- partials <- rgbeta(n - 1, alpha)
    if(n == 2) {
      L[2,2] <- sqrt(1 - L[2,1]^2)
      if(cholesky) return(L)
      Sigma <- tcrossprod(L)

      return(Sigma)
    }
    W <- log(1 - partials^2)
    for(i in 2:(n - 1)) {
      gap <- (i+1):n
      gap1 <- i:(n-1)
      alpha <- alpha - 0.5
      partials <- rgbeta(n - i, alpha)
      L[i,i] <- exp(0.5 * W[i-1])
      L[gap,i] <- partials * exp(0.5 * W[gap1])
      W[gap1] <- W[gap1] + log(1 - partials^2)
    }
    L[n,n] <- exp(0.5 * W[n-1])
    if(cholesky) return(L)
    Sigma <- tcrossprod(L)
    if(!cholesky) {
      ord <- sample(n)
      Sigma <- Sigma[ord,ord]
    }
    return(Sigma)
  }

