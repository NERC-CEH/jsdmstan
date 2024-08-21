#' Fit jsdm models in Stan
#'
#' This function fits joint Species Distribution models in Stan, using either a
#' generalised linear latent variable model approach (\code{method = "gllvm"}), or a
#' multivariate generalised linear mixed model approach (\code{method = "mglmm"}).
#'
#' @details
#'   Environmental covariate effects (\code{"betas"}) can be parameterised in two
#'   ways. With the \code{"cor"} parameterisation all covariate effects are assumed
#'   to be constrained by a correlation matrix between the covariates. With the
#'   \code{"unstruct"} parameterisation all covariate effects are assumed to draw
#'   from a simple distribution with no correlation structure. Both parameterisations
#'   can be modified using the prior object.
#'   Families supported are the Gaussian family, the negative binomial family,
#'   the Poisson family, the binomial family (with number of trials specificied
#'   using the \code{Ntrials} parameter), the Bernoulli family (the special case
#'   of the binomial family where number of trials is equal to one), the
#'   zero-inflated Poisson and the zero-inflated negative binomial. For both
#'   zero-inflated families the zero-inflation is assumed to be a species-specific
#'   constant.
#'
#' @param formula The formula of covariates that the species means are modelled from
#'
#' @param data Dataframe or list of covariates.
#'
#' @param Y Matrix of species by sites. Rows are assumed to be sites, columns are
#'   assumed to be species
#'
#' @param X The covariates matrix, with rows being site and columns being covariates.
#'   Ignored in favour of data when formula approach is used to specify model.
#'
#' @param method Whether to fit a GLLVM or MGLMM model, details in description
#'
#' @param D The number of latent variables within a GLLVM model
#'
#' @param family is the response family, must be one of \code{"gaussian"},
#'   \code{"neg_binomial"}, \code{"poisson"}, \code{"binomial"},
#'   \code{"bernoulli"}, or \code{"zi_poisson"}. Regular expression
#'   matching is supported.
#'
#' @param species_intercept Whether the model should be fit with an intercept by
#'   species, by default \code{TRUE}
#'
#' @param dat_list Alternatively, data can be given to the model as a list containing
#'   Y, X, N, S, K, and site_intercept. See output of [jsdm_sim_data()] for an
#'   example of how this can be formatted.
#'
#' @param site_intercept Whether a site intercept should be included, potential
#'   values \code{"none"} (no site intercept), \code{"grouped"} (a site intercept
#'   with hierarchical grouping) or \code{"ungrouped"} (site intercept with no
#'   grouping)
#'
#' @param site_groups If the site intercept is grouped, a vector of group identities
#'   per site
#'
#' @param Ntrials For the binomial distribution the number of trials, given as
#'   either a single integer which is assumed to be constant across sites or as
#'   a site-length vector of integers.
#'
#' @param prior Set of prior specifications from call to [jsdm_prior()]
#'
#' @param save_data If the data used to fit the model should be saved in the model
#'   object, by default TRUE.
#'
#' @param iter A positive integer specifying the number of iterations for each chain,
#'   default 4000.
#'
#' @param log_lik Whether the log likelihood should be calculated in the generated
#'   quantities (by default TRUE), required for loo
#'
#' @param beta_param The parameterisation of the environmental covariate effects, by
#'   default \code{"unstruct"}. See details for further information.
#'
#' @param zi_param For the zero-inflated families, whether the zero-inflation parameter
#'   is a species-specific constant (default, \code{"constant"}), or varies by
#'   environmental covariates (\code{"covariate"}).
#'
#' @param zi_X If \code{zi = "covariate"}, the matrix of environmental predictors
#'   that the zero-inflation is modelled in response to. If there is not already
#'   an intercept column (identified by all values being equal to one), one will
#'   be added to the front of the matrix.
#'
#' @param ... Arguments passed to [rstan::sampling()]
#'
#' @return A \code{jsdmStanFit} object, comprising a list including the StanFit
#'   object, the data used to fit the model plus a few other bits of information. See
#'   [jsdmStanFit] for details.
#'
#' @examples
#' \dontrun{
#' # MGLMM - specified by using the mglmm aliases and with direct reference to Y and
#' # X matrices:
#' mglmm_data <- mglmm_sim_data(
#'   N = 100, S = 10, K = 3,
#'   family = "gaussian"
#' )
#' mglmm_fit <- stan_mglmm(
#'   Y = mglmm_data$Y, X = mglmm_data$X,
#'   family = "gaussian"
#' )
#' mglmm_fit
#'
#' # You can also run a model by supplying the data as a list:
#' gllvm_data <- jsdm_sim_data(
#'   method = "gllvm", N = 100, S = 6, D = 2,
#'   family = "bernoulli"
#' )
#' gllvm_fit <- stan_jsdm(
#'   dat_list = gllvm_data, method = "gllvm",
#'   family = "bernoulli"
#' )
#' gllvm_fit
#' }
#' @export
stan_jsdm <- function(X, ...) UseMethod("stan_jsdm")

#' @describeIn stan_jsdm this is the default way of doing things
#' @export
stan_jsdm.default <- function(X = NULL, Y = NULL, species_intercept = TRUE, method,
                              dat_list = NULL, family, site_intercept = "none",
                              D = NULL, prior = jsdm_prior(), site_groups = NULL,
                              beta_param = "unstruct", Ntrials = NULL,
                              zi_param = "constant", zi_X = NULL,
                              save_data = TRUE, iter = 4000, log_lik = TRUE, ...) {
  family <- match.arg(family, c("gaussian", "bernoulli", "poisson",
                                "neg_binomial","binomial", "zi_poisson",
                                "zi_neg_binomial"))
  beta_param <- match.arg(beta_param, c("cor", "unstruct"))
  zi_param <- match.arg(zi_param, c("constant","covariate"))
  if(grepl("zi", family) & zi_param == "covariate"){
    if(is.null(zi_X) & is.null(dat_list)){
      message("If zi_param = 'covariate' and no zi_X matrix is supplied then the X matrix is used")
      zi_X <- X
    }
  } else{
    zi_X <- NULL
  }


  stopifnot(
    is.logical(species_intercept),
    is.logical(save_data)
  )
  if(site_intercept == "grouped" & is.null(site_groups) & is.null(dat_list))
    stop("If site_intercept is grouped then groups must be supplied to site_groups")

  # validate data
  data_list <- validate_data(
    Y = Y, X = X, species_intercept = species_intercept,
    D = D, site_intercept = site_intercept, site_groups = site_groups,
    dat_list = dat_list,
    family = family, method = method, Ntrials = Ntrials,
    zi_X = zi_X
  )

  # Create stancode
  model_code <- jsdm_stancode(
    family = family,
    method = method, prior = prior,
    log_lik = log_lik, site_intercept = site_intercept,
    beta_param = beta_param, zi_param = zi_param
  )

  # Compile model
  model_comp <- rstan::stan_model(model_code = model_code)

  model_args <- list(...)
  if (any(c("pars", "include") %in% names(model_args))) {
    warning("Specified pars and include arguments are ignored")
  }
  model_args$object <- model_comp
  model_args$data <- data_list
  model_args$iter <- iter
  model_args$pars <- if (method == "gllvm") c("L", "LV_uncor", "Lambda_uncor") else NA
  model_args$include <- ifelse(method == "gllvm", FALSE, TRUE)

  # Fit model
  model_fit <- do.call(rstan::sampling, model_args)

  model_output <- model_to_jsdmstanfit(
    model_fit = model_fit, method = method,
    data_list = data_list,
    species_intercept = species_intercept,
    family = family, save_data = save_data
  )

  return(model_output)
}

#' @describeIn stan_jsdm Formula interface
#' @export
stan_jsdm.formula <- function(formula, data = list(), ...) {
  mf <- stats::model.frame(formula = formula, data = data)
  X <- stats::model.matrix(attr(mf, "terms"), data = mf)
  if (!is.null(stats::model.response(mf))) {
    warning("Response variable in formula is ignored")
  }

  est <- stan_jsdm.default(X, species_intercept = FALSE, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}

#' Alias for \code{stan_jsdm} with \code{method = "mglmm"}
#' @inheritParams stan_jsdm
#' @export
stan_mglmm <- function(X, ...) UseMethod("stan_mglmm")

#' @describeIn stan_mglmm Default
#' @export
stan_mglmm.default <- function(X = NULL, Y = NULL, species_intercept = TRUE,
                               dat_list = NULL, family, site_intercept = "none",
                               prior = jsdm_prior(),
                               save_data = TRUE, iter = 4000, ...) {
  stan_jsdm.default(
    X = X, Y = Y, species_intercept = species_intercept, method = "mglmm",
    dat_list = dat_list, family, site_intercept = site_intercept,
    prior = prior,
    save_data = save_data, iter = iter, ...
  )
}
#' @describeIn stan_mglmm Formula interface
#' @export
stan_mglmm.formula <- function(formula, data = list(), ...) {
  stan_jsdm.formula(formula, data = data, method = "mglmm", ...)
}


#' Alias for \code{stan_jsdm} with \code{method = "gllvm"}
#' @inheritParams stan_jsdm
#' @export
stan_gllvm <- function(X, ...) UseMethod("stan_gllvm")

#' @describeIn stan_gllvm Default
#' @export
stan_gllvm.default <- function(X = NULL, Y = NULL, D = NULL, species_intercept = TRUE,
                               dat_list = NULL, family, site_intercept = "none",
                               prior = jsdm_prior(),
                               save_data = TRUE, iter = 4000, ...) {
  stan_jsdm.default(
    X = X, Y = Y, D = D, species_intercept = species_intercept,
    method = "gllvm", prior = prior,
    dat_list = dat_list, family, site_intercept = site_intercept,
    save_data = save_data, iter = iter, ...
  )
}
#' @describeIn stan_gllvm Formula interface
#' @export
stan_gllvm.formula <- function(formula, data = list(), ...) {
  stan_jsdm.formula(formula, data = data, method = "gllvm", ...)
}


# Internal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

validate_data <- function(Y, D, X, species_intercept,
                          dat_list, family, site_intercept,
                          method, site_groups, Ntrials,
                          zi_X) {
  method <- match.arg(method, c("gllvm", "mglmm"))

  # do things if data not given as list:
  if (is.null(dat_list)) {
    if (is.null(X) & !species_intercept) {
      stop("Model requires a species intercept if there are no covariates")
    }

    if (method == "gllvm") {
      if (is.null(D)) {
        stop("Must have at least one latent variable")
      }
      if (D < 1) {
        stop("Must have at least one latent variable")
      }
    }

    S <- ncol(Y)
    N <- nrow(Y)
    if (is.null(X)) {
      K <- 1
      X <- matrix(1, nrow = N, ncol = 1)
      colnames(X) <- "(Intercept)"
    } else {
      if(is.null(colnames(X))){
        message("No column names specified for X, assigning names")
        colnames(X) <- paste0("V",seq_len(ncol(X)))
      }
      K <- ncol(X) + 1 * species_intercept
      if(is.data.frame(X)){
        X <- as.matrix(X)
      }
      if (isTRUE(species_intercept)) {
        X <- cbind(matrix(1, nrow = N, ncol = 1), X)
        colnames(X)[1] <- "(Intercept)"
      }
    }
    if(grepl("zi", family) & !is.null(zi_X)){
      if(is.null(colnames(zi_X))){
        message("No column names specified for zi_X, assigning names")
        colnames(zi_X) <- paste0("V",seq_len(ncol(zi_X)))
      }
      if(!any(apply(zi_X, 2, function(x) all(x == 1)))){
        zi_X <- cbind("(Intercept)" = 1, zi_X)
      }
      zi_k <- ncol(zi_X)
    }

    if(site_intercept == "grouped"){
      if(length(site_groups) != N)
        stop("site_groups must be a vector of site identities of the same length as the number of sites")
      suppressWarnings(grps <- as.integer(as.factor(site_groups)))
      if(anyNA(grps))
        stop("site_groups must be coercible to a numeric vector with no NAs")
      ngrp <- max(grps)
    }

    if (method == "mglmm") {
      data_list <- list(
        Y = Y, S = S, K = K, X = X, N = N
      )
    } else if (method == "gllvm") {
      data_list <- list(
        Y = Y, D = D, S = S, K = K, X = X, N = N
      )
    }
    if(site_intercept == "grouped"){
      data_list$ngrp <- ngrp
      data_list$grps <- grps
    }
    if(family == "binomial"){
      data_list$Ntrials <- Ntrials
    }
    if(grepl("zi_",family) & !is.null(zi_X)){
      data_list$zi_k <- zi_k
      data_list$zi_X <- zi_X

      if(nrow(zi_X) != N){
        stop("Number of rows of zi_X must be equal to number of rows of X")
      }
    }
  } else {
    if (!all(c("Y", "K", "S", "N", "X") %in% names(dat_list))) {
      stop("If supplying data as a list must have entries Y, K, S, N, X")
    }

    if (!("D" %in% names(dat_list)) & method == "gllvm") {
      stop("If supplying data as a list must have a D entry")
    }

    if (identical(family, "binomial")) {
      if (!all(c("Ntrials") %in% names(dat_list))) {
        stop("Binomial models require Ntrials in dat_list")
      }
    }

    if (grepl("zi_", family) & !is.null(zi_X)) {
      if (!all(c("zi_X","zi_k") %in% names(dat_list))) {
        stop("Zero-inflated models with the covariate parameterisation of zi require zi_X and zi_k in dat_list")
      }
    }

    if (site_intercept == "grouped") {
      if (!all(c("ngrp","grps") %in% names(dat_list))) {
        stop("Grouped site intercept models require ngrp and grps in dat_list")
      }
    }

    if (isTRUE(species_intercept)) {
      dat_list$X <- cbind(matrix(1, nrow = dat_list$N, ncol = 1), dat_list$X)
      colnames(dat_list$X)[1] <- "(Intercept)"
    }

    data_list <- dat_list
  }

  if (!is.matrix(data_list$Y)) {
    if (is.data.frame(data_list$Y)) {
      if (all(sapply(data_list$Y, is.numeric))) {
        data_list$Y <- as.matrix(data_list$Y)
      } else {
        stop("Non-numeric column(s) detected in Y")
      }
    } else {
      stop("Y must be either a numeric matrix or a dataframe consisting of only numeric columns")
    }
  }

  # Check if Y is appropriate given family
  if (identical(family, "bernoulli")) {
    if (!isTRUE(all.equal(
      unname(data_list$Y),
      matrix(as.numeric(as.logical(data_list$Y)),
        nrow = nrow(data_list$Y)
      )
    ))) {
      stop("Y matrix is not binary")
    }
  } else if (family %in% c("poisson", "neg_binomial", "binomial",
                           "zi_poisson", "zi_neg_binomial")) {
    if (!any(apply(data_list$Y, 1:2, is.wholenumber))) {
      stop("Y matrix is not composed of integers")
    }
  }

  # check to make sure no completely blank columns in Y
  if(any(apply(data_list$Y, 2, function(x) all(x == 0)))){
    stop("Y contains an empty column, which cannot work for this model")
  }

  # Check if Ntrials is appropriate given
  if(identical(family, "binomial")) {
    data_list$Ntrials <- ntrials_check(data_list$Ntrials, data_list$N)
  }

  # create zeros/non-zeros for zero-inflated poisson
  if(grepl("zi_",family)) {
    if(any(apply(data_list$Y, 2, min)>0)){
      stop("Zero-inflated distributions require zeros to be present in all Y values.")
    }
    data_list$N_zero <- colSums(data_list$Y==0)
    data_list$N_nonzero <- colSums(data_list$Y>0)
    data_list$Sum_nonzero <- sum(data_list$N_nonzero)
    data_list$Sum_zero <- sum(data_list$N_zero)
    data_list$Y_nz <- c(as.matrix(data_list$Y))[c(as.matrix(data_list$Y))>0]
    data_list$nn <- rep(1:data_list$N,data_list$S)[c(data_list$Y>0)]
    data_list$ss <- rep(1:data_list$S,each=data_list$N)[c(data_list$Y>0)]
    data_list$nz <- rep(1:data_list$N,data_list$S)[c(data_list$Y==0)]
    data_list$sz <- rep(1:data_list$S,each=data_list$N)[c(data_list$Y==0)]
  }

  return(data_list)
}

model_to_jsdmstanfit <- function(model_fit, method, data_list, species_intercept,
                                 family, save_data) {
  # Turn into jsdmStanFit
  sites <- if (!is.null(rownames(data_list$Y))) {
    rownames(data_list$Y)
  } else {
    as.character(1:nrow(data_list$Y))
  }
  species <- if (!is.null(colnames(data_list$Y))) {
    colnames(data_list$Y)
  } else {
    as.character(1:ncol(data_list$Y))
  }
  preds <- if (!is.null(colnames(data_list$X))) {
    colnames(data_list$X)
  } else {
    as.character(1:ncol(data_list$X))
  }
  if (isTRUE(species_intercept) & !("(Intercept)" %in% preds)) {
    preds <- c("(Intercept)", preds)
  }
  dat <- if (isTRUE(save_data)) {
    data_list
  } else {
    list()
  }
  fam <- list(family = family,
              params = switch(family,
                              "gaussian" = "sigma",
                              "bernoulli" = character(),
                              "poisson" = character(),
                              "neg_binomial" = "kappa",
                              "binomial" = character(),
                              "zi_poisson" = "zi",
                              "zi_neg_binomial" = c("kappa","zi")),
              params_dataresp= character(),
              preds = character(),
              data_list = list())
  class(fam) <- "jsdmStanFamily"
  if(isTRUE(save_data) & ("zi_X" %in% names(data_list))){
    fam$params_dataresp <- "zi"
    fam$preds <- colnames(data_list$zi_X)
    fam$data_list <- list(zi_X = data_list$zi_X)
  }

  model_output <- list(
    fit = model_fit,
    jsdm_type = method,
    family = fam,
    species = species,
    sites = sites,
    preds = preds,
    data_list = dat,
    n_latent = ifelse(method == "gllvm",
      as.integer(round(data_list$D, 0)), integer()
    )
  )

  class(model_output) <- "jsdmStanFit"

  return(model_output)
}
