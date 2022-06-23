#' Plot the traceplots and density plots for parameters within a jsdmStanFit object
#'
#' This function takes parameters from the jsdmStanfit object and plots them using
#' the [bayesplot::mcmc_combo()]function from the \pkg{bayesplot} package.
#'
#' @details This acts as an interface to the [bayesplot::mcmc_combo()]
#'   function, with the default being to plot a density plot and a trace plot for
#'   each parameter specified, although this can be changed by setting the
#'   \code{combo} argument. These jSDM models have a lot of parameters so selecting a
#'   subset is recommended. If pars is set to \code{NULL} (the default) all
#'   parameters with either sigma or kappa in their name will be plotted along with a
#'   random selection of the other parameters (total number of other parameters set
#'   by \code{sample_n}). To see the name of the parameters within the model use
#'   [get_parnames()] - and if you want to plot all parameters (there will be
#'   hundreds in any reasonably sized model) set \code{pars = get_parnames(x)}.
#'
#' @param x The \code{jsdmStanFit} model object
#' @param pars The parameters to plot, by default a random sample of twenty of the
#'   parameters fit within the model
#' @param combo Which combination of plot types within [bayesplot::mcmc_combo()]
#'   to use, by default \code{c("dens", "trace")}
#' @param N The number of plots per page, default \code{5}
#' @param ask Whether to ask before plotting a new page, default \code{TRUE}
#' @param inc_warmup Whether to include the warmup period in the plots, by default
#'   \code{FALSE}
#' @param include Whether to include or exclude the parameters specified by pars, by
#'   default \code{TRUE} (i.e. include)
#' @param plot Whether to plot the plots, default \code{TRUE}
#' @param sample_n If \code{pars = NULL} then the number of random non-sigma
#'   parameters to include (details in description)
#' @param regexp If pars should be treated as a regular expression for matching to
#'   parnames, by default \code{FALSE}
#' @param newpage Whether the first plot should be plotted on a new page, by default
#'   \code{TRUE}
#' @param ... Arguments passed to [bayesplot::mcmc_combo()]
#'
#' @return An invisible list of the plots#
#'
#' @examples
#' \donttest{
#'  #First simulate data and get model fit:
#'  mglmm_data <- mglmm_sim_data(N = 100, S = 10, K = 3,
#'                               family = "gaussian")
#'  mglmm_fit <- stan_mglmm(Y = mglmm_data$Y, X = mglmm_data$X,
#'                          family = "gaussian")
#'
#'  # The default plot:
#'  plot(mglmm_fit)
#'
#'  # Plotting specifically the L_Rho_species parameters:
#'  plot(mglmm_fit, pars = "L_Rho_species", regexp = TRUE)
#'
#'  # Increasing the number of randomly sampled parameters to plot:
#'  plot(mglmm_fit, sample_n = 20)
#'
#' }
#'
#' @export
#'
#' @seealso [mcmc_plot.jsdmStanFit()] for more plotting options.
#'
plot.jsdmStanFit <- function(x, pars = NULL, combo = c("dens","trace"), N = 5L,
                             ask = TRUE, inc_warmup = FALSE, include = TRUE,
                             sample_n = 10,
                             regexp = FALSE, plot = TRUE, newpage = TRUE, ...){
  if(!is.wholenumber(N))
    stop("N must be a positive integer")
  if(!is.wholenumber(sample_n) & is.null(pars))
    stop("If pars is NULL then sample_n must be a positive integer")
  parnames <- get_parnames(x)
  model_pars <- par_sample(pars = pars, parnames = parnames, sample_n = sample_n,
                           regexp = regexp)

  model_est <- rstan::extract(x$fit, pars = model_pars, permuted = FALSE,
                              inc_warmup = inc_warmup, include = include)
  variables <- dimnames(model_est)[[3]]

  if (plot) {
    default_ask <- grDevices::devAskNewPage()
    on.exit(grDevices::devAskNewPage(default_ask))
    grDevices::devAskNewPage(ask = FALSE)
  }
  n_plots <- ceiling(length(variables) / N)
  plots <- vector(mode = "list", length = n_plots)
  for (i in seq_len(n_plots)) {
    sub_vars <- variables[((i - 1) * N + 1):min(i * N, length(variables))]
    sub_draws <- model_est[, , sub_vars, drop = FALSE]
    plots[[i]] <- bayesplot::mcmc_combo(
      sub_draws, combo = combo,  ...
    )
    if (plot) {
      plot(plots[[i]], newpage = newpage || i > 1)
      if (i == 1) {
        grDevices::devAskNewPage(ask = ask)
      }
    }
  }
  invisible(plots)
}


is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
  if(is.numeric(x)){
    (abs(x - round(x)) < tol) & x >=0
  } else(FALSE)
}
par_sample <- function(pars, parnames, sd = c("sigma","kappa"), sample_n, regexp){
  if(is.null(pars)){
    sd_pars <- grep(paste0(sd, collapse = "|"),parnames, value = TRUE)
    model_pars <- c(sd_pars,
                    sample(parnames[!parnames %in% sd_pars], sample_n))

  } else{
    if(isFALSE(regexp)){
      if(!all(pars %in% parnames)){
        stop("Please specify pars within the model, use get_parnames to find the names")
      }
      model_pars <- pars
    } else{
      model_pars <- grep(paste0(pars, collapse = "|"),
                         parnames, value = TRUE)
    }
  }
  model_pars
}

#' MCMC plots implemented in \pkg{bayesplot}
#'
#' An interface for calling the MCMC plotting functions implemented in the
#' \pkg{bayesplot} package.
#'
#' This acts as an interface to the plotting functions prefaced with \code{mcmc_}
#' within the \link[bayesplot:bayesplot-package]{bayesplot} package. The default is
#' to plot an interval plot for the parameters specified, for a full list of
#' available plot options run \code{bayesplot::available_mcmc()} or the error message
#' when you set plotfun to an unrecognised plot function will give you a list of options.
#' These jSDM models have a lot of parameters so selecting a subset is recommended.
#' If pars is set to \code{NULL} (the default) all parameters with either sigma or
#' kappa in their name will be plotted along with a random selection of the other
#' parameters (total number of other parameters set by \code{sample_n}). To see the
#' name of the parameters within the model use [get_parnames()] - and if you want to
#' plot all parameters (there will be hundreds in any reasonably sized model) set
#' \code{pars = get_parnames(x)}.
#'
#' @param x The jsdmStanFit model object
#' @param pars The parameters to plot, by default a random sample of twenty of the
#'   parameters fit within the model
#' @param plotfun The MCMC plot function from \pkg{bayesplot} to be used
#' @param sample_n If \code{pars = NULL} then the number of random non-sigma
#'   parameters to include (details in description)
#' @param inc_warmup Whether to include the warmup period in the plots, by default
#'   \code{FALSE}
#' @param include Whether to include or exclude the parameters specified by pars, by
#'   default \code{TRUE} (i.e. include)
#' @param regexp If pars should be treated as a regular expression for matching to
#'   parnames, by default \code{FALSE}
#' @param ... Other arguments to be passed to the MCMC plotting function
#'
#' @return A [ggplot][ggplot2::ggplot] object that can be customised using the
#'   \pkg{ggplot2} package
#'
#' @examples
#'
#' \donttest{
#'  # First simulate data and fit the model:
#'  gllvm_data <- jsdm_sim_data(method = "gllvm", N = 100, S = 6, D = 2,
#'                              family = "bernoulli")
#'  gllvm_fit <- stan_jsdm(dat_list = gllvm_data, method = "gllvm",
#'                         family = "bernoulli")
#'
#'  # Default is to plot the intervals:
#'  mcmc_plot(gllvm_fit)
#'
#'  # Other plot types include options to see parameter recovery (if the
#'  # jsdm_sim_data functions are used the original parameters are saved within the
#'  # data list)
#'  mcmc_plot(gllvm_fit, plotfun = "recover_intervals",
#'            pars = c("LV[2,20]","sigmas_b[1]","sigma_L"),
#'            true = c(gllvm_data$pars$LV[2,20],
#'                     gllvm_data$pars$beta_sds,
#'                     gllvm_data$pars$L_sigma))
#'
#' }
#'
#' @export
#'
#' @seealso [plot.jsdmStanFit()]
#'
mcmc_plot.jsdmStanFit <- function(x, pars = NULL, plotfun = "intervals",
                                  sample_n = 10, inc_warmup = FALSE, include = TRUE,
                                  regexp = FALSE, ...){
  valid_types <- gsub("^mcmc_","",as.character(bayesplot::available_mcmc()))
  plotfun <- ifelse(grepl("^mcmc_",plotfun),gsub("^mcmc_","",plotfun),plotfun)
  if(!plotfun %in% valid_types){
    stop(paste("Invalid plotfun argument. Valid plot functions are:",
               paste(valid_types, collapse = ", ")))
  }
  mcmc_fun <- get(paste0("mcmc_",plotfun), asNamespace("bayesplot"))

  if(!is.wholenumber(sample_n) & is.null(pars))
    stop("If pars is NULL then sample_n must be a positive integer")
  parnames <- get_parnames(x)
  model_pars <- par_sample(pars = pars, parnames = parnames, sample_n = sample_n,
                           regexp = regexp)

  mcmc_arg_names <- names(formals(mcmc_fun))
  mcmc_args <- list(...)
  if ("x" %in% mcmc_arg_names) {
    if (grepl("^nuts_", plotfun)) {
      # x refers to a molten data.frame of NUTS parameters
      mcmc_args$x <- nuts_params(x)
    } else {
      # x refers to a data.frame of draws
      draws <- rstan::extract(x$fit, pars = model_pars, permuted = FALSE,
                              inc_warmup = inc_warmup, include = include)
      sel_variables <- dimnames(draws)[[3]]
      if (plotfun %in% c("scatter", "hex") && length(sel_variables) != 2L) {
        stop(paste("Exactly 2 parameters must be selected for this plot function.",
                   "\nParameters selected: ", paste(sel_variables, collapse = ", ")))
      }
      mcmc_args$x <- draws
    }
  }
  if ("lp" %in% mcmc_arg_names) {
    mcmc_args$lp <- log_posterior(x, inc_warmup = inc_warmup)
  }
  if ("np" %in% mcmc_arg_names) {
    mcmc_args$np <- nuts_params(x, inc_warmup = inc_warmup)
  }
  interval_type <- plotfun %in% c("intervals", "areas")
  if ("rhat" %in% mcmc_arg_names && !interval_type) {
    mcmc_args$rhat <- rhat(x, pars = model_pars)
  }
  if ("ratio" %in% mcmc_arg_names) {
    mcmc_args$ratio <- neff_ratio(x, pars = model_pars)
  }
  do.call(mcmc_fun, mcmc_args)

}

#' @rdname mcmc_plot.jsdmStanFit
#' @export
mcmc_plot <- function(x, ...) {
  UseMethod("mcmc_plot")
}


#' Plotting an ordination plot for a GLLVM model
#'
#' This function takes a GLLVM model fit and plots an ordination plot with a random
#' (or specified) selection of draws
#'
#' @param object The \code{jsdmStanFit} model object
#' @param choices Which latent variables to plot as dimensions, by default
#'   \code{c(1,2)}
#' @param type Whether to plot sites or species, default \code{"species"}.
#' @param summary_stat The summary statistic used to plot overall averages of the
#'   posterior sample. By default this is \code{"mean"}, and \code{NULL} will result
#'   in no summary being included
#' @param ndraws How many individual draws to include in plot, by default \code{20}. Setting
#'   this to \code{0} will result in no individual draws being included
#' @param draw_ids Which draws to include in plot (overrides \code{ndraws})
#' @param size The size of the points in the graph, specified as a two-element vector
#'   with the first being used for the summary points and the second the individual
#'   draws, default \code{c(2,1)}
#' @param alpha The transparency/alpha of the points in the graph, specified as a
#'   two-element vector with the first being used for the summary points and the
#'   second the individual draws, default \code{c(1,0.5)}
#' @param shape The shape of the points in the graph, specified as a two-element
#'   vector with the first being used for the summary points and the second the
#'   individual draws, default \code{c(18,16)}
#'
#' @return A [ggplot][ggplot2::ggplot] object that can be customised using the
#'   \pkg{ggplot2} package
#' @export
#'
#' @examples
#'
#' \donttest{
#'  # First simulate data and fit the model:
#'  gllvm_data <- jsdm_sim_data(method = "gllvm", N = 100, S = 6, D = 3,
#'                              family = "bernoulli")
#'  gllvm_fit <- stan_jsdm(dat_list = gllvm_data, method = "gllvm",
#'                         family = "bernoulli")
#'
#' ordiplot(gllvm_fit)
#' # now plot the 1st and 3rd latent variables against each other for the sites:
#' ordiplot(gllvm_fit, choices = c(1,3), type = "sites")
#'
#' }
ordiplot <- function(object, choices = c(1,2), type = "species",
                     summary_stat = "mean", ndraws = 20, draw_ids = NULL,
                     size = c(2,1), alpha = c(1,0.5), shape = c(18,16)){
  if(class(object) != "jsdmStanFit")
    stop("Only objects of class jsdmStanFit are supported")
  if(object$jsdm_type != "gllvm")
    stop("Only gllvm models are supported")
  type <- match.arg(type, c("species","sites"))
  if(length(choices) != 2L)
    stop("Only two latent variables can be plotted at once")

  # Extract corrected latent variable scores for species OR sites
  ext_pars <- switch(type, "species" = "Lambda\\[",
                     "sites" = "LV\\[")
  model_est <- extract(object, pars = ext_pars, regexp = TRUE)

  n_iter <- dim(model_est[[1]])[1]
  varnames <- switch(type, "species"= c("draw","S","LV"),
                     "sites" = c("draw","LV","S"))


  if(!is.null(summary_stat)){
    model_est_copy <- model_est
  }

  if(ndraws > 0 | !is.null(draw_ids)){
    if(!is.null(draw_ids)){
      if(max(draw_ids)>n_iter)
        stop(paste("Maximum of draw_ids (",max(draw_ids),
                   ") is greater than number of iterations (",n_iter,")"))

      draw_id <- draw_ids
    } else{
      if(!is.null(ndraws)){
        if(n_iter < ndraws){
          warning(paste("There are fewer samples than ndraws specified, defaulting",
                        "to using all iterations"))
          ndraws <- n_iter
        }
        draw_id <- sample.int(n_iter, ndraws)


      } else{
        draw_id <- seq_len(n_iter)
      }
    }
    model_est <- lapply(model_est, function(x){
      switch(length(dim(x)),
             `1` = x[draw_id,drop=FALSE],
             `2` = x[draw_id,,drop=FALSE],
             `3` =  x[draw_id,,,drop=FALSE])
    })
    # Turn into long format
    ord_scores <- reshape2::melt(model_est[[1]],
                                 varnames = varnames)
    ord_scores <- subset(ord_scores, LV %in% choices)
    ord_scores[,"LV"] <- paste0("LV",ord_scores[,"LV"])
    ord_scores <- reshape2::dcast(ord_scores, S + draw ~ LV)
    ord_scores$S <- as.factor(ord_scores$S)
  }

  if(!is.null(summary_stat)){
    if(is.character(summary_stat)){
      stat_fun <- get(summary_stat)
    } else if(class(summary_stat) == "function"){
      stat_fun <- summary_stat
    }

    ord_scores_summary <- reshape2::melt(model_est_copy[[1]],
                                         varnames = varnames)
    ord_scores_summary <- subset(ord_scores_summary, LV %in% choices)
    ord_scores_summary[,"LV"] <- paste0("LV",ord_scores_summary[,"LV"])
    ord_scores_summary$S <- as.factor(ord_scores_summary$S)
    ord_scores_summary <- stats::aggregate(value ~ S + LV, ord_scores_summary, stat_fun)
    ord_scores_summary <- reshape2::dcast(ord_scores_summary, S ~ LV)
  }

  graph <- ggplot2::ggplot() +
    bayesplot::bayesplot_theme_get() +
    ggplot2::coord_equal()

  if(!is.null(summary_stat)){
    graph <- graph +
      ggplot2::geom_point(data = ord_scores_summary,
                          ggplot2::aes_string(paste0("LV",choices[1]),
                                              paste0("LV",choices[2]),
                                              colour = "S"),
                          size = size[1], alpha = alpha[1], shape = shape[1])
  }

  if(ndraws > 0 | !is.null(draw_ids)){
    graph <- graph +
      ggplot2::geom_point(data = ord_scores,
                          ggplot2::aes_string(paste0("LV",choices[1]),
                                              paste0("LV",choices[2]),
                                              colour = "S"),
                          size = size[2], alpha = alpha[2], shape = shape[2])
  }
  graph

}
