#' Generate simulated data within a variety of jSDM methodologies
#'
#' The \code{jsdm_sim_data} function can simulate data with either a
#' multivariate generalised mixed model (MGLMM) or a generalised linear latent
#' variable model (GLLVM). The \code{gllvm_sim_data} and \code{mglmm_sim_data}
#' are aliases for \code{jsdm_sim_data} that set \code{method} to \code{"gllvm"}
#' and \code{"mglmm"} respectively.
#'
#' @details This simulates data based on a joint species distribution model with
#'   either a generalised linear latent variable model approach or a
#'   multivariate generalised linear mixed model approach.
#'
#'   Models can be fit with or without "measured predictors", and if measured
#'   predictors are included then the species have species-specific parameter
#'   estimates. These can either be simulated completely independently, or have
#'   information pooled across species. If information is pooled this can be
#'   modelled as either a random draw from some mean and standard deviation or
#'   species covariance can be modelled together (this will be the covariance
#'   used in the overall model if the method used has covariance).
#'
#'   Environmental covariate effects (\code{"betas"}) can be parameterised in
#'   two ways. With the \code{"cor"} parameterisation all covariate effects are
#'   assumed to be constrained by a correlation matrix between the covariates.
#'   With the \code{"unstruct"} parameterisation all covariate effects are
#'   assumed to draw from a simple distribution with no correlation structure.
#'   Both parameterisations can be modified using the prior object.
#'
#' @export
#'
#' @param S is number of species
#'
#' @param N is number of sites
#'
#' @param D is number of latent variables, used within gllvm method
#'
#' @param K is number of covariates, by default \code{0}
#'
#' @param family is the response family, must be one of \code{"gaussian"},
#'   \code{"neg_binomial"}, \code{"poisson"}, \code{"binomial"},
#'   \code{"bernoulli"}, \code{"zi_poisson"}, or \code{"zi_neg_binomial"}.
#'   Regular expression matching is supported.
#'
#' @param method is the jSDM method to use, currently either \code{"gllvm"} or
#'   \code{"mglmm"} - see details for more information.
#'
#' @param species_intercept Whether to include an intercept in the predictors,
#'   must be \code{TRUE} if \code{K} is \code{0}. Defaults to \code{TRUE}.
#'
#' @param Ntrials For the binomial distribution the number of trials, given as
#'   either a single integer which is assumed to be constant across sites or as
#'   a site-length vector of integers.
#'
#' @param site_intercept Whether a site intercept should be included, potential
#'   values \code{"none"} (no site intercept) or \code{"ungrouped"} (site
#'   intercept with no grouping). Defaults to no site intercept, grouped is not
#'   supported currently.
#'
#' @param beta_param The parameterisation of the environmental covariate
#'   effects, by default \code{"unstruct"}. See details for further information.
#'
#' @param zi_param For the zero-inflated families, whether the zero-inflation
#'   parameter is a species-specific constant (default, \code{"constant"}), or
#'   varies by environmental covariates (\code{"covariate"}).
#'
#' @param zi_k If \code{zi_param="covariate"}, the number of environmental
#'   covariates that the zero-inflation parameter responds to. The default
#'   (\code{NULL}) is that the zero-inflation parameter responds to exactly the
#'   same covariate matrix as the mean parameter. Otherwise, a different set of
#'   random environmental covariates are generated, plus an intercept (not
#'   included in zi_k) and used to predict zero-inflation. Will be ignored if
#'   zi_X is supplied.
#'
#' @param shp_param For families with shape parameters, whether the shape
#'   parameter is a species-specific constant (default, \code{"constant"}), or
#'   varies by environmental covariates (\code{"covariate"}).
#'
#' @param shp_k If \code{shp_param="covariate"}, the number of environmental
#'   covariates that the shape parameter responds to. The default (\code{NULL})
#'   is that the shape parameter responds to exactly the same covariate matrix
#'   as the mean parameter. Otherwise, a different set of random environmental
#'   covariates are generated and used to predict the shape parameter. Will be
#'   ignored if shp_X is supplied.
#'
#' @param prior Set of prior specifications from call to [jsdm_prior()]
#'
#' @param X The X matrix to be used to simulate data, by default \code{NULL} -
#'   i.e. the X matrix is simulated using a random draw from a standard normal
#'   distribution.
#'
#' @param zi_X The zi_X matrix to be used to simulate data, by default \code{NULL} -
#'   i.e. the zi_X matrix is simulated using a random draw from a standard normal
#'   distribution.
#'
#' @param shp_X The shp_X matrix to be used to simulate data, by default \code{NULL} -
#'   i.e. the shp_X matrix is simulated using a random draw from a standard normal
#'   distribution.
jsdm_sim_data <- function(S, N = NULL, D = NULL, K = 0L, family,
                          method = c("gllvm", "mglmm"),
                          species_intercept = TRUE,
                          Ntrials = NULL,
                          site_intercept = "none",
                          beta_param = "unstruct",
                          zi_param = "constant", zi_k = NULL,
                          shp_param = "constant", shp_k = NULL,
                          prior = jsdm_prior(),
                          X = NULL, zi_X = NULL, shp_X = NULL) {
  response <- match.arg(family, c("gaussian", "neg_binomial", "poisson",
                                  "bernoulli", "binomial", "zi_poisson",
                                  "zi_neg_binomial"))
  site_intercept <- match.arg(site_intercept, c("none","ungrouped","grouped"))
  beta_param <- match.arg(beta_param, c("cor", "unstruct"))
  zi_param <- match.arg(zi_param, c("constant","covariate"))
  shp_param <- match.arg(shp_param, c("constant","covariate"))

  if(missing(method)){
    stop("method argument needs to be specified")
  }
  if(site_intercept == "grouped"){
    stop("Grouped site intercept not supported")
  }

  if (!is.double(S)) {
    stop("S must be a positive integer")
  }

  if (S <= 0| S %% 1 != 0) {
    stop("S must be a positive integer")
  }

  if(is.null(X)){
    if(!is.double(N))
      stop("N must be a positive integer")
    if(N <= 0 | N %% 1 != 0)
      stop("N must be a positive integer")
    if ((K < 0 | K %% 1 != 0)) {
      stop("K must be either 0 or a positive integer")
    }
    if ((K == 0 & isFALSE(species_intercept))) {
      stop("If K is 0 then a species intercept is required")
    }
  } else{
    if(!is.matrix(X)){
      stop("X must be a matrix")
    }
    K <- ncol(X)
    N <- nrow(X)
    species_intercept <- FALSE
    if(is.null(colnames(X))){
      message("X has been provided with no column names, assigning names")
      colnames(X) <- paste0("V", 1:K)
    }
  }



  if (method == "gllvm" & !is.double(D)) {
    stop("gllvm method require D to be a positive integer")
  }

  if (class(prior)[1] != "jsdmprior") {
    stop("prior object must be of class jsdmprior, produced by jsdm_prior()")
  }

  if(response == "binomial"){
    Ntrials <- ntrials_check(Ntrials = Ntrials, N = N)
  }

  if(is.null(zi_X)){
    if (!is.null(zi_k)) {
      if(zi_k < 1 | zi_k %% 1 != 0){
        stop("zi_k must be either NULL or a positive integer")
      }
      ZI_K <- zi_k
    } else {
      ZI_K <- K
    }
  } else{
    if(!is.matrix(zi_X)){
      stop("zi_X must be a matrix")
    }
    if(nrow(zi_X) != N){
      stop("zi_X must have N rows")
    }
    zi_k <- ncol(zi_X)
    ZI_K <- zi_k
    if(is.null(colnames(zi_X))){
      message("zi_X has been provided with no column names, assigning names")
      colnames(zi_X) <- paste0("V", 1:zi_k)
    }
    if(!any(apply(zi_X, 2, function(x) all(x == 1)))){
      message("zi_X has been provided without an intercept, so one has been added")
      zi_X <- cbind("(Intercept)" = 1, zi_X)
    }
  }

  if(is.null(shp_X)){
    if (!is.null(shp_k)) {
      if(shp_k < 1 | shp_k %% 1 != 0){
        stop("shp_k must be either NULL or a positive integer")
      }
      SHP_K <- shp_k
    } else {
      SHP_K <- K
    }
  } else{
    if(!is.matrix(shp_X)){
      stop("shp_X must be a matrix")
    }
    if(nrow(shp_X) != N){
      stop("shp_X must have N rows")
    }
    shp_k <- ncol(shp_X)
    SHP_K <- shp_k
    if(is.null(colnames(shp_X))){
      message("shp_X has been provided with no column names, assigning names")
      colnames(shp_X) <- paste0("V", 1:shp_k)
    }
    if(!any(apply(shp_X, 2, function(x) all(x == 1)))){
      message("shp_X has been provided without an intercept, so one has been added")
      shp_X <- cbind("(Intercept)" = 1, shp_X)
    }
  }




  # prior object breakdown
  prior_split <- lapply(prior, strsplit, split = "\\(|\\)|,")
  if (!all(sapply(prior_split, function(x) {
    x[[1]][1] %in% c(
      "normal",
      "inv_gamma",
      "lkj_corr_cholesky",
      "lkj_corr",
      "student_t",
      "cauchy",
      "gamma",
      "beta"
    )
  }))) {
    stop(paste(
      "Not all prior distributions specified are supported.",
      "Raise an issue at github.com/NERC-CEH/jsdmstan if you want a new",
      "distribution added to jsdm_sim_data"
    ))
  }
  D1 <- switch(method,
    "gllvm" = D,
    "mglmm" = 1
  )
  prior_func <- lapply(names(prior_split), function(x) {
    y <- prior_split[x]
    # print(str(y[[1]][[1]]))
    fun_name <- switch(y[[1]][[1]][1],
      "normal" = "rnorm",
      "inv_gamma" = "rinvgamma",
      "lkj_corr_cholesky" = "rlkj",
      "lkj_corr" = "rlkj",
      "student_t" = "rstudentt",
      "cauchy" = "rcauchy",
      "gamma" = "rgamma",
      "beta" = "rbeta"
    )
    fun_arg1 <- switch(x,
      "sigmas_preds" = K + 1 * species_intercept,
      "z_preds" = (K + 1 * species_intercept) * S,
      "cor_preds" = K + 1 * species_intercept,
      "betas" = (K + 1 * species_intercept) * S,
      "a" = N,
      "a_bar" = 1,
      "sigma_a" = 1,
      "sigmas_species" = S,
      "z_species" = S * N,
      "cor_species" = S,
      "LV" = D1 * N,
      "L" = D1 * (S - D1) + (D1 * (D1 - 1) / 2) + D1,
      "sigma_L" = 1,
      "sigma" = S,
      "kappa" = S,
      "zi" = S,
      "zi_betas" = S*(ZI_K+1),
      "shp_betas" = S*(SHP_K + 1)
    )
    fun_args <- as.list(c(fun_arg1, as.numeric(unlist(y[[1]][[1]])[-1])))

    if(x == "cor_preds" & grepl("lkj_corr_cholesky\\(",prior$cor_species))
      fun_args <- c(fun_args,1)
    if(x == "cor_species" & grepl("lkj_corr_cholesky\\(",prior$cor_species))
      fun_args <- c(fun_args,1)


    return(list(fun_name, fun_args))
  })
  names(prior_func) <- names(prior_split)


  # build species covariance matrix
  if (method == "mglmm") {
    cor_species <- do.call(
      match.fun(prior_func[["cor_species"]][[1]]),
      prior_func[["cor_species"]][[2]]
    )
  }

  # now do covariates - if K = NULL then do intercept only
  if(is.null(X)){
    if (K == 0) {
      x <- matrix(1, nrow = N, ncol = 1)
      colnames(x) <- "(Intercept)"
      J <- 1
    } else if (isTRUE(species_intercept)) {
      x <- matrix(stats::rnorm(N * K), ncol = K, nrow = N)
      colnames(x) <- paste0("V", 1:K)
      x <- cbind("(Intercept)" = 1, x)
      J <- K + 1
    } else if (isFALSE(species_intercept)) {
      x <- matrix(stats::rnorm(N * K), ncol = K, nrow = N)
      colnames(x) <- paste0("V", 1:K)
      J <- K
    } else {
      stop("K must be an integer value")
    }
  } else{
    x <- X
    J <- K
  }
  # print(str(x))

  # covariate parameters
  if(beta_param == "cor"){
    sigmas_preds <- abs(do.call(
      match.fun(prior_func[["sigmas_preds"]][[1]]),
      prior_func[["sigmas_preds"]][[2]]
    ))
    z_preds <- matrix(do.call(
      match.fun(prior_func[["z_preds"]][[1]]),
      prior_func[["z_preds"]][[2]]
    ), ncol = S, nrow = J)
    if (K == 0) {
      beta_sim <- sigmas_preds %*% z_preds
    } else {
      cor_preds <- do.call(
        match.fun(prior_func[["cor_preds"]][[1]]),
        prior_func[["cor_preds"]][[2]]
      )
      beta_sim <- (diag(sigmas_preds) %*% cor_preds) %*% z_preds
    }
  } else if (beta_param == "unstruct"){
    beta_sim <- matrix(do.call(
      match.fun(prior_func[["betas"]][[1]]),
      prior_func[["betas"]][[2]]
    ), ncol = S, nrow = J)
  }

  mu_sim <- x %*% beta_sim


  ## site intercept
  if (site_intercept %in% c("grouped","ungrouped")) {
    a_bar <- do.call(
      match.fun(prior_func[["a_bar"]][[1]]),
      prior_func[["a_bar"]][[2]]
    )
    sigma_a <- abs(do.call(
      match.fun(prior_func[["sigma_a"]][[1]]),
      prior_func[["sigma_a"]][[2]]
    ))
    a <- do.call(
      match.fun(prior_func[["a"]][[1]]),
      prior_func[["a"]][[2]]
    )
    a_i <- a_bar + a * sigma_a
  } else {
    a_i <- rep(0, N)
  }

  if (method == "mglmm") {
    # u covariance
    sigmas_species <- abs(do.call(
      match.fun(prior_func[["sigmas_species"]][[1]]),
      prior_func[["sigmas_species"]][[2]]
    ))
    z_species <- matrix(do.call(
      match.fun(prior_func[["z_species"]][[1]]),
      prior_func[["z_species"]][[2]]
    ), nrow = S, ncol = N)
    u_ij <- t((diag(sigmas_species) %*% cor_species) %*% z_species)
  }

  if (method == "gllvm") {
    L_l <- do.call(
      match.fun(prior_func[["L"]][[1]]),
      prior_func[["L"]][[2]]
    )
    L <- matrix(nrow = S, ncol = D)

    idx2 <- 0
    if(D > 1){
      for (i in 1:(D - 1)) {
        for (j in (i + 1):(D)) {
          L[i, j] <- 0
        }
      }
    }
    for (j in 1:D) {
      for (i in j:S) {
        idx2 <- idx2 + 1
        L[i, j] <- L_l[idx2]
      }
    }

    sigma_L <- abs(do.call(
      match.fun(prior_func[["sigma_L"]][[1]]),
      prior_func[["sigma_L"]][[2]]
    ))

    LV <- matrix(do.call(
      match.fun(prior_func[["LV"]][[1]]),
      prior_func[["LV"]][[2]]
    ), nrow = D, ncol = N)

    LV_sum <- (L * sigma_L) %*% LV
  }

  # variance parameters
  if (response == "gaussian") {
    if(shp_param == "constant"){
      sigma <- abs(do.call(
        match.fun(prior_func[["sigma"]][[1]]),
        prior_func[["sigma"]][[2]]
      ))
    } else if(shp_param == "covariate"){
      shp_betas <- matrix(do.call(
        match.fun(prior_func[["shp_betas"]][[1]]),
        prior_func[["shp_betas"]][[2]]
      ), ncol = S)
    }
  } else if (response == "neg_binomial") {
    if(shp_param == "constant"){
      kappa <- abs(do.call(
        match.fun(prior_func[["kappa"]][[1]]),
        prior_func[["kappa"]][[2]]
      ))
    } else if(shp_param == "covariate"){
      shp_betas <- matrix(do.call(
        match.fun(prior_func[["shp_betas"]][[1]]),
        prior_func[["shp_betas"]][[2]]
      ), ncol = S)
    }
  } else if (response  == "zi_poisson") {
    if(zi_param == "covariate"){
      zi_betas <- matrix(do.call(
        match.fun(prior_func[["zi_betas"]][[1]]),
        prior_func[["zi_betas"]][[2]]
      ), ncol = S)
    } else {
      zi <- do.call(
        match.fun(prior_func[["zi"]][[1]]),
        prior_func[["zi"]][[2]]
      )
    }
  } else if (response  == "zi_neg_binomial") {
    if(zi_param == "covariate"){
      zi_betas <- matrix(do.call(
        match.fun(prior_func[["zi_betas"]][[1]]),
        prior_func[["zi_betas"]][[2]]
      ), ncol = S)
    } else {
      zi <- do.call(
        match.fun(prior_func[["zi"]][[1]]),
        prior_func[["zi"]][[2]]
      )
    }
    if(shp_param == "constant"){
      kappa <- abs(do.call(
        match.fun(prior_func[["kappa"]][[1]]),
        prior_func[["kappa"]][[2]]
      ))
    } else if(shp_param == "covariate"){
      shp_betas <- matrix(do.call(
        match.fun(prior_func[["shp_betas"]][[1]]),
        prior_func[["shp_betas"]][[2]]
      ), ncol = S)
    }
  }
  # print(str(sigma))

  # zero-inflation in case of covariates
  if(grepl("zi_", response) & zi_param == "covariate"){
    if(is.null(zi_X)){
      if(is.null(zi_k)){
        zi_X <- x
      } else {
        zi_X <- matrix(stats::rnorm(N * zi_k), ncol = zi_k, nrow = N)
        colnames(zi_X) <- paste0("V", 1:zi_k)
        zi_X <- cbind("(Intercept)" = 1, zi_X)
      }
    }
    zi <- inv_logit(zi_X %*% zi_betas)
  }

  # zero-inflation in case of covariates
  if(shp_param == "covariate"){
    if(is.null(shp_X)){
      if(is.null(shp_k)){
        shp_X <- x
      } else {
        shp_X <- matrix(stats::rnorm(N * shp_k), ncol = shp_k, nrow = N)
        colnames(shp_X) <- paste0("V", 1:shp_k)
        shp_X <- cbind("(Intercept)" = 1, shp_X)

      }
    }
    if(family == "gaussian"){
     sigma <- exp(shp_X %*% shp_betas)
    } else if(family %in% c("neg_binomial","zi_neg_binomial")){
      kappa <- exp(shp_X %*% shp_betas)
    }
  }

  # generate Y
  Y <- matrix(nrow = N, ncol = S)
  for (i in 1:N) {
    for (j in 1:S) {
      mu_ij <- switch(method,
        "gllvm" = a_i[i] + mu_sim[i, j] + LV_sum[j, i],
        "mglmm" = a_i[i] + mu_sim[i, j] + u_ij[i, j]
      )
      # print(str(mu_ij))

      Y[i, j] <- switch(response,
        "neg_binomial" = rgampois(1,
          mu = exp(mu_ij),
          scale = switch(shp_param, "constant" = kappa[j],
                         "covariate" = kappa[i,j])
        ),
        "gaussian" = stats::rnorm(1, mu_ij, switch(shp_param, "constant" = sigma[j],
                                                   "covariate" = sigma[i,j])),
        "poisson" = stats::rpois(1, exp(mu_ij)),
        "bernoulli" = stats::rbinom(1, 1, inv_logit(mu_ij)),
        "binomial" = stats::rbinom(1, Ntrials[i], inv_logit(mu_ij)),
        "zi_poisson" = (1-stats::rbinom(
          1, 1,
          ifelse(zi_param == "covariate",
                 zi[i,j],zi[j])))*stats::rpois(1, exp(mu_ij)),
        "zi_neg_binomial" = (1-stats::rbinom(
          1, 1,
          ifelse(zi_param == "covariate",
                 zi[i,j],zi[j])))*rgampois(1, mu = exp(mu_ij), scale = kappa[j])
      )
    }
  }

  if(any(apply(Y, 2, function(x) all(x == 0)))){
    message(paste("Y contains an entirely empty column, which will not work for",
                  "jsdm fitting, it is recommended that the simulation is run again."))
  }



  pars <- list(
    betas = beta_sim
  )

  if(beta_param == "cor"){
    pars$sigmas_preds <- sigmas_preds
    pars$z_preds <- z_preds
    if(K != 0){
      pars$cor_preds <- cor_preds
    }
  }

  if (site_intercept == "ungrouped") {
    pars$a_bar <- a_bar
    pars$sigma_a <- sigma_a
    pars$a <- a
  }
  if (method == "gllvm") {
    pars$L <- L
    pars$LV <- LV
    pars$sigma_L <- sigma_L
  }
  if (method == "mglmm") {
    pars$sigmas_species <- sigmas_species
    pars$cor_species <- cor_species
    pars$z_species <- z_species
  }
  if (response == "gaussian") {
    pars$sigma <- sigma
  }
  if (response == "neg_binomial") {
    pars$kappa <- kappa
  }
  if(response == "gaussian" & shp_param == "constant"){
    pars$sigma <- sigma
  }
  if(response == "neg_binomial" & shp_param == "constant"){
    pars$kappa <- kappa
  }
  if(shp_param == "covariate"){
    pars$shp_betas <- shp_betas
  }
  if(grepl("zi_", response)){
    if(zi_param == "constant"){
      pars$zi <- zi
    } else if(zi_param == "covariate"){
      pars$zi_betas <- zi_betas
    }
  }
  if(response == "zi_neg_binomial"){
    pars$kappa <- kappa
  }
  if (isTRUE(species_intercept)) {
    if (K > 0) {
      x <- x[, 2:ncol(x), drop = FALSE]
    } else {
      x <- NULL
    }
  }
  output <- list(
    Y = Y, pars = pars, N = N, S = S, D = D, K = J, X = x
  )
  if(response == "binomial"){
    output$Ntrials <- Ntrials
  }
  if(grepl("zi_", response) & zi_param == "covariate"){
    output$zi_k <- ZI_K + 1
    output$zi_X <- zi_X
  }
  if(shp_param == "covariate"){
    output$shp_k <- SHP_K + 1
    output$shp_X <- shp_X
  }

  return(output)
}

#'
#'
#' @export
#' @describeIn jsdm_sim_data Alias for \code{jsdm_sim_data} with \code{method =
#'  "gllvm"}
#'
#' @param ... Arguments passed to jsdm_sim_data
gllvm_sim_data <- function(...) {
  jsdm_sim_data(method = "gllvm", ...)
}

#'
#' @export
#' @describeIn jsdm_sim_data Alias for \code{jsdm_sim_data} with \code{method =
#'  "mglmm"}
#'
#' @param ... Arguments passed to jsdm_sim_data
mglmm_sim_data <- function(...) {
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
#' The alternative parameterisation of the Student T distribution is by mean (mu) and
#' scale (sigma) to be consistent with the Stan parameterisation rather than the
#' parameterisation in [stats::rt()].
#'
#' @param n The number of samples to create/dimension of correlation matrix
#' @param mu The mean used within the negative binomial parameterisation and the
#'   Student T distribution
#' @param scale The phi parameter that controls overdispersion of the negative
#'   binomial distribution (see details for description), or the scale parameter used
#'   within the inverse gamma distribution (see [stats::rgamma()])
#' @param eta The shape parameter of the LKJ distribution
#' @param cholesky Whether the correlation matrix should be returned as the Cholesky
#'   decomposition, by default \code{FALSE}
#' @param shape The shape parameter of the inverse gamma distribution (see
#'   [stats::rgamma()])
#' @param df The degrees of freedom parameter within the Student T distribution (see
#'    details)
#' @param sigma The scale of the Student T distribution (see details)
#'
#' @name sim_helpers
NULL

#' @rdname sim_helpers
#' @export
rgampois <- function(n, mu, scale) {
  # shape <- mu/scale
  prob <- scale / (mu + scale)
  stats::rnbinom(n, size = scale, prob = prob)
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}


# The following is code from Ben Goodrich's response on the stan google mailing list
# (see link https://groups.google.com/g/stan-users/c/3gDvAs_qwN8/m/Xpgi2rPlx68J)

rgbeta <-
  function(n, shape) {
    if (shape == Inf) {
      rep(0, n)
    } else if (shape > 0) {
      -1 + 2 * stats::rbeta(n, shape, shape)
    } else if (shape == 0) {
      -1 + 2 * stats::rbinom(n, 1, 0.5)
    } else {
      stop("shape must be non-negative")
    }
  }

#' @rdname sim_helpers
#' @export
rlkj <-
  function(n, eta = 1, cholesky = FALSE) {
    if (n < 2) {
      stop("Dimension of correlation matrix must be >= 2")
    }
    if (eta < 1) {
      stop("The value of eta must be >= 1")
    }
    alpha <- eta + (n - 2) / 2
    L <- matrix(0, n, n)
    L[1, 1] <- 1
    L[-1, 1] <- partials <- rgbeta(n - 1, alpha)
    if (n == 2) {
      L[2, 2] <- sqrt(1 - L[2, 1]^2)
      if (cholesky) {
        return(L)
      }
      Sigma <- tcrossprod(L)

      return(Sigma)
    }
    W <- log(1 - partials^2)
    for (i in 2:(n - 1)) {
      gap <- (i + 1):n
      gap1 <- i:(n - 1)
      alpha <- alpha - 0.5
      partials <- rgbeta(n - i, alpha)
      L[i, i] <- exp(0.5 * W[i - 1])
      L[gap, i] <- partials * exp(0.5 * W[gap1])
      W[gap1] <- W[gap1] + log(1 - partials^2)
    }
    L[n, n] <- exp(0.5 * W[n - 1])
    if (cholesky) {
      return(L)
    }
    Sigma <- tcrossprod(L)
    if (!cholesky) {
      ord <- sample(n)
      Sigma <- Sigma[ord, ord]
    }
    return(Sigma)
  }

#' @rdname sim_helpers
#' @export
rinvgamma <- function(n, shape, scale) {
  1 / stats::rgamma(n, shape, scale = scale)
}
#' @rdname sim_helpers
#' @export
rstudentt <- function(n, df, mu, sigma) {
  mu + sigma * stats::rt(n, df = df)
}

ntrials_check <- function(Ntrials, N){
  if(is.null(Ntrials)){
    stop("Number of trials must be specified for the binomial distribution")
  }
  if(!is.double(Ntrials) & !is.integer(Ntrials)){
    stop("Ntrials must be a positive integer")
  }
  if(!(length(Ntrials) %in% c(1, N))){
    stop("Ntrials must be of length 1 or N")
  }
  if(length(Ntrials) == 1L){
    Ntrials <- rep(Ntrials, N)
  }
  Ntrials
}
