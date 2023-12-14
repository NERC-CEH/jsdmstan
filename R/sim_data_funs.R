#' Generate simulated data within a variety of jSDM methodologies
#'
#' The \code{jsdm_sim_data} function can simulate data with either a multivariate
#' generalised mixed model (MGLMM) or a generalised linear latent variable model
#' (GLLVM). The \code{gllvm_sim_data} and \code{mglmm_sim_data} are aliases for
#' \code{jsdm_sim_data} that set \code{method} to \code{"gllvm"} and \code{"mglmm"}
#' respectively.
#'
#' @details This simulates data based on a joint species distribution model with
#'   either a generalised linear latent variable model approach or a multivariate
#'   generalised linear mixed model approach.
#'
#'   Models can be fit with or without "measured predictors", and if measured
#'   predictors are included then the species have species-specific parameter
#'   estimates. These can either be simulated completely independently, or have
#'   information pooled across species. If information is pooled this can be modelled
#'   as either a random draw from some mean and standard deviation or species
#'   covariance can be modelled together (this will be the covariance used in the
#'   overall model if the method used has covariance).
#'
#'   Environmental covariate effects (\code{"betas"}) can be parameterised in two
#'   ways. With the \code{"cor"} parameterisation all covariate effects are assumed
#'   to be constrained by a correlation matrix between the covariates. With the
#'   \code{"unstruct"} parameterisation all covariate effects are assumed to draw
#'   from a simple distribution with no correlation structure. Both parameterisations
#'   can be modified using the prior object.
#'
#' @export
#'
#' @param N is number of sites
#'
#' @param S is number of species
#'
#' @param D is number of latent variables, used within gllvm method
#'
#' @param K is number of covariates, by default \code{0}
#'
#' @param family is the response family, must be one of \code{"gaussian"},
#'   \code{"neg_binomial"}, \code{"poisson"} or \code{"bernoulli"}. Regular
#'   expression matching is supported.
#'
#' @param method is the jSDM method to use, currently either \code{"gllvm"} or
#'   \code{"mglmm"} - see details for more information.
#'
#' @param species_intercept Whether to include an intercept in the predictors, must
#'   be \code{TRUE} if \code{K} is \code{0}. Defaults to \code{TRUE}.
#'
#' @param site_intercept Whether a site intercept should be included, potential
#'   values \code{"none"} (no site intercept) or \code{"ungrouped"} (site intercept
#'   with no grouping). Defaults to no site intercept, grouped is not supported
#'   currently.
#'
#' @param beta_param The parameterisation of the environmental covariate effects, by
#'   default \code{"unstruct"}. See details for further information.
#'
#' @param prior Set of prior specifications from call to [jsdm_prior()]
jsdm_sim_data <- function(N, S, D = NULL, K = 0L, family, method = c("gllvm", "mglmm"),
                          species_intercept = TRUE,
                          site_intercept = "none",
                          beta_param = "unstruct",
                          prior = jsdm_prior()) {
  response <- match.arg(family, c("gaussian", "neg_binomial", "poisson", "bernoulli"))
  site_intercept <- match.arg(site_intercept, c("none","ungrouped","grouped"))
  beta_param <- match.arg(beta_param, c("cor", "unstruct"))
  if(site_intercept == "grouped"){
    stop("Grouped site intercept not supported")
  }

  if (any(!c(is.double(N), is.double(S)))) {
    stop("N and S must be positive integers")
  }

  if (any(c(N, S) <= 0) | any(c(N, S) %% 1 != 0)) {
    stop("N and S must be positive integers")
  }

  if (K < 0 | K %% 1 != 0) {
    stop("K must be either 0 or a positive integer")
  }
  if (K == 0 & isFALSE(species_intercept)) {
    stop("If K is 0 then a species intercept is required")
  }

  if (method == "gllvm" & !is.double(D)) {
    stop("gllvm method require D to be a positive integer")
  }

  if (class(prior)[1] != "jsdmprior") {
    stop("prior object must be of class jsdmprior, produced by jsdm_prior()")
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
      "gamma"
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
      "gamma" = "rgamma"
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
      "sigma" = 1,
      "kappa" = 1
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
  if (K == 0) {
    x <- matrix(1, nrow = N, ncol = 1)
    colnames(x) <- "Intercept"
    J <- 1
  } else if (isTRUE(species_intercept)) {
    x <- matrix(stats::rnorm(N * K), ncol = K, nrow = N)
    colnames(x) <- paste0("V", 1:K)
    x <- cbind(Intercept = 1, x)
    J <- K + 1
  } else if (isFALSE(species_intercept)) {
    x <- matrix(stats::rnorm(N * K), ncol = K, nrow = N)
    colnames(x) <- paste0("V", 1:K)
    J <- K
  } else {
    stop("K must be an integer value")
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
    for (i in 1:(D - 1)) {
      for (j in (i + 1):(D)) {
        L[i, j] <- 0
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
    sigma <- abs(do.call(
      match.fun(prior_func[["sigma"]][[1]]),
      prior_func[["sigma"]][[2]]
    ))
  } else if (response == "neg_binomial") {
    kappa <- abs(do.call(
      match.fun(prior_func[["kappa"]][[1]]),
      prior_func[["kappa"]][[2]]
    ))
  }
  # print(str(sigma))

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
          scale = kappa
        ),
        "gaussian" = stats::rnorm(1, mu_ij, sigma),
        "poisson" = stats::rpois(1, exp(mu_ij)),
        "bernoulli" = stats::rbinom(1, 1, inv_logit(mu_ij))
      )
    }
  }
  # print(str(Y))


  pars <- list(
    betas = beta_sim
  )

  if(beta_param == "cor"){
    pars$sigmas_preds <- sigmas_preds
    pars$z_preds <- z_preds
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
  if (isTRUE(species_intercept)) {
    if (K > 0) {
      x <- x[, 2:ncol(x)]
    } else {
      x <- NULL
    }
  }
  output <- list(
    Y = Y, pars = pars, N = N, S = S, D = D, K = J, X = x
  )

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
  1 / (1 + exp(x))
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
