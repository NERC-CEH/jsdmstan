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
#' \dontrun{
#' # First simulate data and get model fit:
#' mglmm_data <- mglmm_sim_data(
#'   N = 100, S = 10, K = 3,
#'   family = "gaussian"
#' )
#' mglmm_fit <- stan_mglmm(
#'   Y = mglmm_data$Y, X = mglmm_data$X,
#'   family = "gaussian"
#' )
#'
#' # The default plot:
#' plot(mglmm_fit)
#'
#' # Plotting specifically the cor_species parameters:
#' plot(mglmm_fit, pars = "cor_species", regexp = TRUE)
#'
#' # Increasing the number of randomly sampled parameters to plot:
#' plot(mglmm_fit, sample_n = 20)
#' }
#'
#' @export
#'
#' @seealso [mcmc_plot.jsdmStanFit()] for more plotting options.
#'
plot.jsdmStanFit <- function(x, pars = NULL, combo = c("dens", "trace"), N = 5L,
                             ask = TRUE, inc_warmup = FALSE, include = TRUE,
                             sample_n = 10,
                             regexp = FALSE, plot = TRUE, newpage = TRUE, ...) {
  if (!is.wholenumber(N)) {
    stop("N must be a positive integer")
  }
  if (!is.wholenumber(sample_n) & is.null(pars)) {
    stop("If pars is NULL then sample_n must be a positive integer")
  }
  parnames <- get_parnames(x)
  model_pars <- par_sample(
    pars = pars, parnames = parnames, sample_n = sample_n,
    regexp = regexp
  )

  model_est <- rstan::extract(x$fit,
    pars = model_pars, permuted = FALSE,
    inc_warmup = inc_warmup, include = include
  )
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
      sub_draws,
      combo = combo, ...
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


is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (is.numeric(x)) {
    (abs(x - round(x)) < tol) & x >= 0
  } else {
    (FALSE)
  }
}
par_sample <- function(pars, parnames, sd = c("sigma", "kappa"), sample_n, regexp) {
  if (is.null(pars)) {
    sd_pars <- grep(paste0(sd, collapse = "|"), parnames, value = TRUE)
    model_pars <- c(
      sd_pars,
      sample(parnames[!parnames %in% sd_pars], sample_n)
    )
  } else {
    if (isFALSE(regexp)) {
      if (!all(pars %in% parnames)) {
        stop("Please specify pars within the model, use get_parnames to find the names")
      }
      model_pars <- pars
    } else {
      model_pars <- grep(paste0(pars, collapse = "|"),
        parnames,
        value = TRUE
      )
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
#' \dontrun{
#' # First simulate data and fit the model:
#' gllvm_data <- jsdm_sim_data(
#'   method = "gllvm", N = 100, S = 6, D = 2,
#'   family = "bernoulli"
#' )
#' gllvm_fit <- stan_jsdm(
#'   dat_list = gllvm_data, method = "gllvm",
#'   family = "bernoulli"
#' )
#'
#' # Default is to plot the intervals:
#' mcmc_plot(gllvm_fit)
#'
#' # Other plot types include options to see parameter recovery (if the
#' # jsdm_sim_data functions are used the original parameters are saved within the
#' # data list)
#' mcmc_plot(gllvm_fit,
#'   plotfun = "recover_intervals",
#'   pars = c("LV[2,20]", "sigmas_preds[1]", "sigma_L"),
#'   true = c(
#'     gllvm_data$pars$LV[2, 20],
#'     gllvm_data$pars$beta_sds,
#'     gllvm_data$pars$L_sigma
#'   )
#' )
#' }
#'
#' @export
#'
#' @seealso [plot.jsdmStanFit()]
#'
mcmc_plot.jsdmStanFit <- function(x, pars = NULL, plotfun = "intervals",
                                  sample_n = 10, inc_warmup = FALSE, include = TRUE,
                                  regexp = FALSE, ...) {
  valid_types <- gsub("^mcmc_", "", as.character(bayesplot::available_mcmc()))
  plotfun <- ifelse(grepl("^mcmc_", plotfun), gsub("^mcmc_", "", plotfun), plotfun)
  if (!plotfun %in% valid_types) {
    stop(paste(
      "Invalid plotfun argument. Valid plot functions are:",
      paste(valid_types, collapse = ", ")
    ))
  }
  mcmc_fun <- get(paste0("mcmc_", plotfun), asNamespace("bayesplot"))

  if (!is.wholenumber(sample_n) & is.null(pars)) {
    stop("If pars is NULL then sample_n must be a positive integer")
  }
  parnames <- get_parnames(x)
  model_pars <- par_sample(
    pars = pars, parnames = parnames, sample_n = sample_n,
    regexp = regexp
  )

  mcmc_arg_names <- names(formals(mcmc_fun))
  mcmc_args <- list(...)
  if ("x" %in% mcmc_arg_names) {
    if (grepl("^nuts_", plotfun)) {
      # x refers to a molten data.frame of NUTS parameters
      mcmc_args$x <- nuts_params(x)
    } else {
      # x refers to a data.frame of draws
      draws <- rstan::extract(x$fit,
        pars = model_pars, permuted = FALSE,
        inc_warmup = inc_warmup, include = include
      )
      sel_variables <- dimnames(draws)[[3]]
      if (plotfun %in% c("scatter", "hex") && length(sel_variables) != 2L) {
        stop(paste(
          "Exactly 2 parameters must be selected for this plot function.",
          "\nParameters selected: ", paste(sel_variables, collapse = ", ")
        ))
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
#' This function takes a GLLVM model fit and plots an ordination plot with a
#' random (or specified) selection of draws
#'
#' @param object The \code{jsdmStanFit} model object
#' @param choices Which latent variables to plot as dimensions, by default
#'   \code{c(1,2)}
#' @param type Whether to plot sites or species, default \code{"species"}.
#' @param summary_stat The summary statistic used to plot overall averages of
#'   the posterior sample. By default this is \code{"median"}, and \code{NULL}
#'   will result in no summary being included
#' @param ndraws How many individual draws to include in plot, by default
#'   \code{0}. Setting this to \code{0} will result in no individual draws
#'   being included
#' @param draw_ids Which draws to include in plot (overrides \code{ndraws})
#' @param size The size of the points in the graph, specified as a two-element
#'   vector with the first being used for the summary points and the second the
#'   individual draws, default \code{c(2,1)}
#' @param alpha The transparency/alpha of the points in the graph, specified as
#'   a two-element vector with the first being used for the summary points and
#'   the second the individual draws, default \code{c(1,0.5)}
#' @param shape The shape of the points in the graph, specified as a two-element
#'   vector with the first being used for the summary points and the second the
#'   individual draws, default \code{c(18,16)}
#' @param geom Which geom from \pkg{ggplot2} is used for the summary statistic
#'   by default \code{"point"}, or alternatively can be \code{"text"}
#' @param errorbar_range If specified, the central range of the data that is
#'   covered by the errorbar. Needs to be given as either \code{NULL} (i.e. no
#'   errorbar), or a number between 0 and 1, default \code{0.75}.
#' @param errorbar_linewidth The linewidth of the error bar, by default 1.
#'
#' @return A [ggplot][ggplot2::ggplot] object that can be customised using the
#'   \pkg{ggplot2} package
#' @export
#'
#' @examples
#' \dontrun{
#' # First simulate data and fit the model:
#' gllvm_data <- jsdm_sim_data(
#'   method = "gllvm", N = 100, S = 6, D = 3,
#'   family = "bernoulli"
#' )
#' gllvm_fit <- stan_jsdm(
#'   dat_list = gllvm_data, method = "gllvm",
#'   family = "bernoulli"
#' )
#'
#' ordiplot(gllvm_fit)
#' # now plot the 1st and 3rd latent variables against each other for the sites:
#' ordiplot(gllvm_fit, choices = c(1, 3), type = "sites")
#' }
ordiplot <- function(object, choices = c(1, 2), type = "species",
                     summary_stat = "median", ndraws = 0, draw_ids = NULL,
                     size = c(2, 1), alpha = c(1, 0.5), shape = c(18, 16),
                     geom = "point", errorbar_range = 0.75,
                     errorbar_linewidth = 1) {
  # try to remove CMD check note
  LV <- NULL
  if (!inherits(object, "jsdmStanFit")) {
    stop("Only objects of class jsdmStanFit are supported")
  }
  if (object$jsdm_type != "gllvm") {
    stop("Only gllvm models are supported")
  }
  type <- match.arg(type, c("species", "sites"))
  if (length(choices) != 2L) {
    stop("Only two latent variables can be plotted at once")
  }
  geom <- match.arg(geom, c("point","text"))
  if(!is.null(errorbar_range)){
    if(!is.numeric(errorbar_range)){
      stop("errorbar_range must be a number between 0 and 1")
    } else if(errorbar_range < 0 | errorbar_range > 1){
      stop("errorbar_range must be a number between 0 and 1")
    }
  }

  # Extract corrected latent variable scores for species OR sites
  ext_pars <- switch(type,
    "species" = "Lambda\\[",
    "sites" = "LV\\["
  )
  model_est <- extract(object, pars = ext_pars, regexp = TRUE)

  n_iter <- dim(model_est[[1]])[1]
  varnames <- switch(type,
    "species" = c("draw", "S", "LV"),
    "sites" = c("draw", "LV", "S")
  )


  if (!is.null(summary_stat) | !is.null(errorbar_range)) {
    model_est_copy <- model_est
  }

  if (ndraws > 0 | !is.null(draw_ids)) {
    if (!is.null(draw_ids)) {
      if (max(draw_ids) > n_iter) {
        stop(paste(
          "Maximum of draw_ids (", max(draw_ids),
          ") is greater than number of iterations (", n_iter, ")"
        ))
      }

      draw_id <- draw_ids
    } else {
      if (!is.null(ndraws)) {
        if (n_iter < ndraws) {
          warning(paste(
            "There are fewer samples than ndraws specified, defaulting",
            "to using all iterations"
          ))
          ndraws <- n_iter
        }
        draw_id <- sample.int(n_iter, ndraws)
      } else {
        draw_id <- seq_len(n_iter)
      }
    }
    model_est <- lapply(model_est, function(x) {
      switch(length(dim(x)),
        `1` = x[draw_id, drop = FALSE],
        `2` = x[draw_id, , drop = FALSE],
        `3` =  x[draw_id, , , drop = FALSE]
      )
    })
    # Turn into long format
    ord_scores <- reshape2::melt(model_est[[1]],
      varnames = varnames
    )
    ord_scores <- subset(ord_scores, LV %in% choices)
    ord_scores[, "LV"] <- paste0("LV", ord_scores[, "LV"])
    ord_scores <- reshape2::dcast(ord_scores, S + draw ~ LV)
    ord_scores$S <- as.factor(ord_scores$S)
  }

  if (!is.null(summary_stat) | !is.null(errorbar_range)) {
    if (is.character(summary_stat)) {
      stat_fun <- get(summary_stat)
    } else if (inherits(summary_stat, "function")) {
      stat_fun <- summary_stat
    } else if(is.null(summary_stat) & !is.null(errorbar_range)){
      stat_fun <- stats::median
    }

    ord_scores_summary <- reshape2::melt(model_est_copy[[1]],
      varnames = varnames
    )
    ord_scores_summary <- subset(ord_scores_summary, LV %in% choices)
    ord_scores_summary[, "LV"] <- paste0("LV", ord_scores_summary[, "LV"])
    ord_scores_summary$S <- as.factor(ord_scores_summary$S)
    ord_scores_summary <- stats::aggregate(value ~ S + LV, ord_scores_summary, stat_fun)
    ord_scores_summary <- reshape2::dcast(ord_scores_summary, S ~ LV)
  }

  if (!is.null(errorbar_range)) {
    errorbar_min = 0.5 - 0.5*errorbar_range
    errorbar_max = 0.5 + 0.5*errorbar_range
    ord_scores_summary_forerrorbar <- reshape2::melt(model_est_copy[[1]],
                                         varnames = varnames
    )
    ord_scores_summary_forerrorbar <- subset(ord_scores_summary_forerrorbar, LV %in% choices)
    ord_scores_summary_forerrorbar[, "LV"] <- paste0("LV", ord_scores_summary_forerrorbar[, "LV"])
    ord_scores_summary_forerrorbar$S <- as.factor(ord_scores_summary_forerrorbar$S)
    ord_scores_summary_forerrorbar <- stats::aggregate(value ~ S + LV, ord_scores_summary_forerrorbar,
                                                       function(x) stats::quantile(x, prob = c(errorbar_min, errorbar_max)))
    ord_scores_summary_forerrorbar$min <- ord_scores_summary_forerrorbar$value[,1]
    ord_scores_summary_forerrorbar$max <- ord_scores_summary_forerrorbar$value[,2]
    ord_scores_summary_forerrorbar <- ord_scores_summary_forerrorbar[,c("S","LV","min","max")]
    ord_scores_summary_forerrorbar <- reshape2::melt(ord_scores_summary_forerrorbar,
                                                     c("S","LV"))
    ord_scores_summary_forerrorbar$ID <- paste0(ord_scores_summary_forerrorbar$LV, "_",
                                                ord_scores_summary_forerrorbar$variable)
    ord_scores_summary_forerrorbar <- reshape2::dcast(ord_scores_summary_forerrorbar, S ~ ID)
    ord_scores_summary_forerrorbar <- cbind(ord_scores_summary_forerrorbar,
                                            ord_scores_summary[,-1])
  }

  graph <- ggplot2::ggplot() +
    bayesplot::bayesplot_theme_get() +
    ggplot2::coord_equal()

  if(!is.null(errorbar_range)){
    graph <- graph +
      ggplot2::geom_linerange(
        data = ord_scores_summary_forerrorbar,
        ggplot2::aes(x = .data[[paste0("LV", choices[1])]],
                     ymin = .data[[paste0("LV", choices[2], "_min")]],
                     ymax = .data[[paste0("LV", choices[2], "_max")]],
                     colour = .data[["S"]]),
        alpha = alpha[1], linewidth = errorbar_linewidth
      ) +
      ggplot2::geom_linerange(
        data = ord_scores_summary_forerrorbar,
        ggplot2::aes(y = .data[[paste0("LV", choices[2])]],
                     xmin = .data[[paste0("LV", choices[1], "_min")]],
                     xmax = .data[[paste0("LV", choices[1], "_max")]],
                     colour = .data[["S"]]),
        alpha = alpha[1], linewidth = errorbar_linewidth
      )
  }

  if (!is.null(summary_stat)) {
    if(geom == "point"){
      graph <- graph +
        ggplot2::geom_point(
          data = ord_scores_summary,
          ggplot2::aes(.data[[paste0("LV", choices[1])]],
                       .data[[paste0("LV", choices[2])]],
                       colour = .data[["S"]]
          ),
          size = size[1], alpha = alpha[1], shape = shape[1]
        )
    } else {
      graph <- graph +
        ggplot2::geom_text(
          data = ord_scores_summary,
          ggplot2::aes(.data[[paste0("LV", choices[1])]],
                       .data[[paste0("LV", choices[2])]],
                       label = switch(type,
                                      "species" = object$species,
                                      "sites" = object$sites)
          ),
          size = size[1], alpha = alpha[1]
        )
    }
  }

  if (ndraws > 0 | !is.null(draw_ids)) {
    graph <- graph +
      ggplot2::geom_point(
        data = ord_scores,
        ggplot2::aes(.data[[paste0("LV", choices[1])]],
                     .data[[paste0("LV", choices[2])]],
                     colour = .data[["S"]]
        ),
        size = size[2], alpha = alpha[2], shape = shape[2]
      )
  }
  graph
}

#' Plotting environmental effects on species
#'
#' @param object The jsdmStanFit model object
#' @param include_intercept Whether to include the intercept in the plots
#' @param plotfun Which plot function from mcmc_plot should be used, by default
#'   \code{"intervals"}
#' @param nrow The number of rows within the plot
#' @param y_labels Which plots should have annotated y axes. Needs to be given
#'   as an integer vector
#' @param widths The widths of the plots
#' @param species The species to be included in the plot, by default \code{NULL}
#'   and all species are included. Should be specified as a vector of species
#'   names or numeric indices.
#' @param preds The predictors to be included in the plot, by default \code{NULL}
#'   and all predictors are included. Should be specified as a vector of predictor
#'   names or numeric indices.
#' @param ... Other arguments that are passed to the MCMC plotting function (see
#'  [mcmc_plot.jsdmStanFit()])
#'
#' @return An object of class \code{"bayesplot_grid"}, for more information see
#'   [bayesplot::bayesplot_grid()]
#' @export
#'
envplot <- function(object, include_intercept = FALSE,
                    nrow = NULL, y_labels = NULL,
                    plotfun = "intervals", widths = NULL,
                    species = NULL, preds = NULL, ...){
  if (!inherits(object, "jsdmStanFit"))
    stop("Only objects of class jsdmStanFit are supported")
  if(is.null(preds)){
    preds <- object$preds$fixef
  } else if(is.numeric(preds)){
    if(min(preds)<1 | max(preds)>length(object$preds$fixef) |!all(is.wholenumber(preds)))
      stop("preds must be an integer vector that corresponds to the predictor indices")
    preds <- object$preds$fixef[preds]
  } else if(is.character(preds)){
    if(!all(preds %in% object$preds$fixef))
      stop("preds must only include predictors included in model fit")
  } else {
    stop("preds must be either a character or integer vector")
  }
  if(is.null(species)){
    species <- object$species
  } else if(is.numeric(species)){
    if(min(species)<1 | max(species)>length(object$species) |!all(is.wholenumber(species)))
      stop("species must be an integer vector that corresponds to the species indices")
    species <- object$species[species]
  } else if(is.character(species)){
    if(!all(species %in% object$species))
      stop("species must only include species included in model fit")
  } else{
    stop("species must be either a character or integer vector")
  }
  if(isFALSE(include_intercept)){
    if(all("preds" == "(Intercept)")){
      stop("include_intercept is FALSE but there are no other predictors")
    }
    preds <- preds[preds != "(Intercept)"]
  }
  pn <- get_parnames(object)
  # sort out species order
  beta_sind <- match(species, object$species)
  species <- species[order(beta_sind, species)]
  if(any(grepl("betas", pn))){
    pl_list <- lapply(preds, function(x){
      beta_pind <- match(x, object$preds$fixef)
      suppressMessages(
        pl <- mcmc_plot(object, plotfun = plotfun,
                pars = paste0("betas[",beta_pind,",",beta_sind,"]"), ...) +
          ggplot2::scale_y_discrete(labels = species) +
          ggplot2::labs(x = x)
      )
      return(pl)
    })
  } else
    stop("This beta parameterisation currently unsupported")

  if(is.null(nrow)){
    nrow <- ceiling(length(preds)/3)
  }
  if(is.null(y_labels)){
    ppr <- ceiling(length(preds)/nrow)
    y_labels <- 1 + ppr*seq(0,nrow-1)
  }
  y_nolabels <- seq_along(preds)[!seq_along(preds) %in% y_labels]
  for(i in y_nolabels){
    pl_list[[i]] <- pl_list[[i]] +
      ggplot2::theme(axis.text.y = ggplot2::element_blank())
  }

  bayesplot::bayesplot_grid(plots = pl_list,
                            grid_args = list(nrow = nrow,
                                             widths = widths))
}

#' Plot modelled correlations between species
#'
#' @param object The jsdmStanFit model object
#' @param species Which species should be included - this plots the correlations
#'   between the specified species and all other species, so by default if all
#'   species are included this will result in a plot per species with all other
#'   species there and thus duplicate entries across all the plots
#' @param plotfun Which plotting function from bayesplot to use, by default
#'   \code{"intervals"}
#' @param nrow How many rows the grid of plots has
#' @param widths Whether the widths of the different rows of plots should vary
#' @param ... Other arguments passed to mcmc_plot
#'
#' @return An object of class \code{"bayesplot_grid"}, for more information see
#'   [bayesplot::bayesplot_grid()]
#' @export
corrplot <- function(object, species = NULL,
                     plotfun = "intervals", nrow = NULL, widths = NULL, ...){
  if (!inherits(object, "jsdmStanFit"))
    stop("Only objects of class jsdmStanFit with method mglmm are supported")
  if(object$jsdm_type != "mglmm")
    stop("Only objects of class jsdmStanFit with method mglmm are supported")

  if (!is.null(species) & !is.character(species)) {
    if (any(!is.wholenumber(species))) {
      stop(paste(
        "Species must be either a character vector of species names or an",
        "integer vector of species positions in the input data columns"
      ))
    }
  }

  plotfun <- ifelse(grepl("^mcmc_", plotfun),
                    plotfun, paste0("mcmc_", plotfun))
  valid_types <- as.character(bayesplot::available_mcmc())
  if (!plotfun %in% valid_types) {
    stop(paste(
      "plotfun:", plotfun, "is not a valid ppc type. ",
      "Valid types are:\n", paste(valid_types, collapse = ", ")
    ))
  }
  ppc_fun <- get(plotfun, asNamespace("bayesplot"))

  if (!is.null(species)) {
    if (is.character(species)) {
      species_names <- object$species
      if (any(!(species %in% species_names))) {
        stop("Species specified are not found in the model fit object")
      }
      species <- match(species, species_names)
    }
    species_names <- object$species[species]
  } else{
    species_names <- object$species
    species <- seq_along(species_names)
  }

  compl_pars <- expand.grid(seq_along(object$species),
                            seq_along(object$species))
  compl_pars <- compl_pars[compl_pars[,1]>compl_pars[,2],]
  pl_list <- lapply(species_names, function(nm){
    x <- match(nm, object$species)
    pars <- compl_pars[compl_pars[,1] == x | compl_pars[,2] == x,]
    pars <- paste0("cor_species[",pars[,1],",",pars[,2],"]")
    suppressMessages(
      pl <- mcmc_plot(object, plotfun = plotfun,
                      pars = pars, ...) +
        ggplot2::scale_y_discrete(labels = object$species[object$species != nm]) +
        ggplot2::labs(x = nm)
    )
    return(pl)
  })

  if(is.null(nrow)){
    nrow <- ceiling(length(species)/3)
  }

  bayesplot::bayesplot_grid(plots = pl_list,
                            grid_args = list(nrow = nrow,
                                             widths = widths))
}

#' Plot the smooth effect of any spline within a model
#'
#' @param object The jsdmStanFit model object
#' @param ndraws The number of draws of the predictor
#' @param draw_ids Which draws to include in plot (overrides \code{ndraws})
#' @param alpha The alpha (transparency) value of the lines
#' @param summarise Whether to include a summary, by default \code{NULL} i.e. no
#'   summary. Can also be given as the name of a function (e.g. \code{"mean"}).
#'   Note that this summarisation only occurs over the number of draws
#'   specified.
#' @param type Whether to return plot on the link scale (default) or the
#'   response scale
#' @param facet_args Argument passed to [ggplot2::facet_wrap()], specified as a
#'   list
#' @param select_smooths Which plots to include, given either as \code{NULL}
#'   (i.e. plot all smooths) or as a list containing two numeric vectors named
#'   "site" and "species" representing the indices of the site-specific smooths
#'   to plot and the species-specific smooths to plot
#'
#' @returns a list of ggplot objects
#'
#' @export
smoothplot <- function(object, ndraws = 20, alpha = 0.4,
                       draw_ids = NULL, summarise = NULL,
                       type = "link", facet_args = list(), select_smooths = NULL){
  # get smooth terms and extract them
  if(!("spl_smooth" %in% names(object$preds)))
    stop("Model must have at least one smooth term")
  type <- match.arg(type, c("link","response"))
  pl_list <- list()
  if(length(object$preds$spl_smooth$nfs)>0){
    # preddata <- .smooth_dataprep(object$preds$spl_smooth$nfs$sm_data)

    betas <- extract(object, pars="nfs_b")[["nfs_b"]]
    n_iter <- dim(betas)[1]
    # take a subsample
    draw_id <- draw_id_check(draw_ids = draw_ids, n_iter = n_iter, ndraws = ndraws)
    betas <- betas[draw_id, ]

    off <- 0
    select_ind <- 0

    # looping over terms
    for(i in seq_along(object$preds$spl_smooth$nfs$sm)){
      select_ind <- select_ind + 1
      if(!is.null(select_smooths)){
        if(!select_ind %in% select_smooths$site){
          next
        }
      }

      ndf <- object$preds$spl_smooth$nfs$sm[[i]]$df

      if(inherits(object$preds$spl_smooth$nfs$sm[[i]],
                  "random.effect")){
        preddata <- unique(object$preds$spl_smooth$nfs$sm_data[,object$preds$spl_smooth$nfs$sm[[i]]$term,drop=FALSE])
        # make projection matrix
        Xp <- mgcv::PredictMat(object$preds$spl_smooth$nfs$sm[[i]],preddata)

        # make predictions (fiddly indexing)
        smoo <- Xp %*% t(betas[, (off+1):(off+ndf)])
        off <- off+ndf
        smoo <- as.data.frame(smoo)
        smoo$cov <- preddata[[object$preds$spl_smooth$nfs$sm[[i]]$term]]

        smoo <- reshape2::melt(smoo, id ="cov")

        if(!is.null(summarise)){
          smoo <- stats::aggregate(smoo[,c("value")], by = list(smoo[["cov"]]),
                                   FUN = summarise)
          smoo$cov <- smoo$Group.1
          smoo$value <- smoo$x
          smoo$variable <- "summary"
          alpha <- 1
        }
        if(type == "response"){
          if(object$family$family %in%
             c("poisson","lognormal","neg_binomial","zi_poisson",
               "zi_neg_binomial","gamma")){
            smoo$value <- exp(smoo$value)
          } else  if(object$family$family %in% c("binomial","bernoulli")){
            smoo$value <- inv_logit(smoo$value)
          }
        }

        term <- paste0(object$preds$spl_smooth$nfs$sm[[i]]$term, collapse = "__")

        # plot by term
        pl_list[[i]] <- ggplot2::ggplot(smoo) +
          ggplot2::geom_jitter(ggplot2::aes(x = .data[["cov"]],
                                          y = .data[["value"]]),
                             alpha = alpha, width = 0.1) +
          # ggplot2::stat_summary(ggplot2::aes(x = .data[["cov"]],
          #                                    y = .data[["value"]]),
          #                       fun.data = "mean_se") +
          bayesplot::bayesplot_theme_get() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0)) +
          ggplot2::labs(x = term)


      } else if("fac" %in% names(object$preds$spl_smooth$nfs$sm[[i]])){
        preddata <- .smooth_dataprep(object$preds$spl_smooth$nfs$sm_data[,object$preds$spl_smooth$nfs$sm[[i]]$term,drop=FALSE])

        ndf <- object$preds$spl_smooth$nfs$sm[[i]]$bs.dim
        nfs_fac <- object$preds$spl_smooth$nfs$sm[[i]]$fterm
        nfslevels <- levels(object$preds$spl_smooth$nfs$sm_data[[nfs_fac]])
        preddata_withlevels <- do.call(rbind,
                                       replicate(length(nfslevels), preddata,
                                                 simplify = FALSE))
        preddata_withlevels[[nfs_fac]] <- as.factor(rep(nfslevels,each=nrow(preddata)))

        # make projection matrix
        Xp <- mgcv::PredictMat(object$preds$spl_smooth$nfs$sm[[i]], preddata_withlevels)

        # make predictions (fiddly indexing)
        smoo <- Xp %*% t(betas[, (off+1):(off+ndf)])
        off <- off+ndf
        smoo <- as.data.frame(smoo)

        # put the covariate values into the data
        nfscov <- object$preds$spl_smooth$nfs$sm[[i]]$term[object$preds$spl_smooth$nfs$sm[[i]]$term != nfs_fac]

        pl_list[[i]] <- .smooth_factorplot(fscov = nfscov, smoo = smoo,
                                           preddata_withlevels = preddata_withlevels,
                                           summarise = summarise, fs_fac = nfs_fac,
                                           alpha= alpha, type = type,
                                           family =object$family$family,
                                           facet_args = facet_args)

      } else{
        preddata <- .smooth_dataprep(object$preds$spl_smooth$nfs$sm_data[,object$preds$spl_smooth$nfs$sm[[i]]$term,drop=FALSE])
        # make projection matrix
        Xp <- mgcv::PredictMat(object$preds$spl_smooth$nfs$sm[[i]], preddata)

        # make predictions (fiddly indexing)
        smoo <- Xp %*% t(betas[, (off+1):(off+ndf)])
        off <- off+ndf
        smoo <- as.data.frame(smoo)

        # put the covariate values into the data
        nfscov <- object$preds$spl_smooth$nfs$sm[[i]]$term
        if(length(nfscov) == 1){
          smoo$cov <- preddata[[object$preds$spl_smooth$nfs$sm[[i]]$term]]

          smoo <- reshape2::melt(smoo, id ="cov")

          if(!is.null(summarise)){
            smoo <- stats::aggregate(smoo[,c("value")], by = list(smoo[["cov"]]),
                                     FUN = summarise)
            smoo$cov <- smoo$Group.1
            smoo$value <- smoo$x
            smoo$variable <- "summary"
            alpha <- 1
          }
          if(type == "response"){
            if(object$family$family %in%
               c("poisson","lognormal","neg_binomial","zi_poisson",
                 "zi_neg_binomial","gamma")){
              smoo$value <- exp(smoo$value)
            } else  if(object$family$family %in% c("binomial","bernoulli")){
              smoo$value <- inv_logit(smoo$value)
            }
          }

          # add some other labels
          term <- paste0(object$preds$spl_smooth$nfs$sm[[i]]$term, collapse = "__")

          # allpred <- rbind(allpred, smoo)
          # plot by term
          pl_list[[i]] <- ggplot2::ggplot(smoo) +
            ggplot2::geom_line(ggplot2::aes(x = .data[["cov"]],
                                            y = .data[["value"]],
                                            group = .data[["variable"]]),
                               alpha = alpha) +
            bayesplot::bayesplot_theme_get() +
            ggplot2::labs(x = term)
        } else if (length(nfscov) == 2) {
          smoo$cov1 <- preddata[[nfscov[1]]]
          smoo$cov2 <- preddata[[nfscov[2]]]

          smoo <- reshape2::melt(smoo, id =c("cov1","cov2"))

          if(!is.null(summarise)){
            smoo <- stats::aggregate(smoo[,c("value")], by = list(smoo[["cov1"]],
                                                                  smoo[["cov2"]]),
                                     FUN = summarise)
            smoo$cov1 <- smoo$Group.1
            smoo$cov2 <- smoo$Group.2
            smoo$value <- smoo$x
            smoo$variable <- "summary"
          }

          # add some other labels
          term <- paste0(nfscov, collapse = "__")

          # allpred <- rbind(allpred, smoo)
          # plot by term
          pl_list[[i]] <- ggplot2::ggplot(smoo) +
            ggplot2::geom_tile(ggplot2::aes(x = .data[["cov1"]],
                                            y = .data[["cov2"]],
                                            fill = .data[["value"]])) +
            bayesplot::bayesplot_theme_get() +
            ggplot2::labs(x = nfscov[1],y = nfscov[2]) +
            ggplot2::scale_fill_distiller(palette = "RdBu", type = "div")
          if(is.null(summarise)){
            facet_args$facets <- ~.data[["variable"]]
            pl_list[[i]] <- pl_list[[i]] +
              ggplot2::facet_wrap(facet_args)
          }
        } else{
          stop("smoothplot can only work with one or two numeric predictors")
        }
      }
    }

  }
  fs_pl_list <- list()
  if(length(object$preds$spl_smooth$fs)>0){
    # preddata <- .smooth_dataprep(object$preds$spl_smooth$fs$sm_data)

    betas <- extract(object, pars="fs_b")[["fs_b"]]
    n_iter <- dim(betas)[1]
    # take a subsample
    draw_id <- draw_id_check(draw_ids = draw_ids, n_iter = n_iter, ndraws = ndraws)
    betas <- betas[draw_id, ]

    off <- 0
    select_ind <- 0

    # looping over terms
    for(i in seq_along(object$preds$spl_smooth$fs$sm)){
      select_ind <- select_ind + 1
      if(!is.null(select_smooths)){
        if(!select_ind %in% select_smooths$species){
          next
        }
      }

      preddata <- .smooth_dataprep(object$preds$spl_smooth$fs$sm_data[,object$preds$spl_smooth$fs$sm[[i]]$term,drop=FALSE])
      ndf <- object$preds$spl_smooth$fs$sm[[i]]$bs.dim

      fs_fac <- object$preds$spl_smooth$fs$sm[[i]]$fterm
      fslevels <- levels(object$preds$spl_smooth$fs$sm_data[[fs_fac]])
      preddata_withlevels <- do.call(rbind,
                                     replicate(length(fslevels), preddata,
                                               simplify = FALSE))
      preddata_withlevels[[fs_fac]] <- as.factor(rep(fslevels,each=nrow(preddata)))

      # make projection matrix
      Xp <- mgcv::PredictMat(object$preds$spl_smooth$fs$sm[[i]], preddata_withlevels)

      # make predictions (fiddly indexing)
      smoo <- Xp %*% t(betas[, (off+1):(off+ndf)])
      off <- off+ndf
      smoo <- as.data.frame(smoo)

      # put the covariate values into the data
      fscov <- object$preds$spl_smooth$fs$sm[[i]]$term[object$preds$spl_smooth$fs$sm[[i]]$term != fs_fac]

      fs_pl_list[[i]] <- .smooth_factorplot(fscov = fscov, smoo = smoo,
                                            preddata_withlevels = preddata_withlevels,
                                            summarise = summarise, fs_fac = fs_fac,
                                            alpha = alpha, type = type,
                                            family =object$family$family,
                                            facet_args = facet_args)
    }

  }
  return(c(pl_list,fs_pl_list))
}

.smooth_dataprep <- function(object_data){
  sm_data_num <- object_data[,sapply(object_data,is.numeric),
                             drop=FALSE]
  datamins <- sapply(sm_data_num, min)
  datamax <- sapply(sm_data_num, max)
  if (length(datamins) == 1){
    preddata <- matrix(nrow = 200, ncol = 1)
    preddata[,1] <- seq(datamins[1],datamax[1],length.out=200)
  } else if(length(datamins) == 2) {
    preddata1 <- seq(datamins[1],datamax[1],length.out=50)
    preddata2 <- seq(datamins[2],datamax[2],length.out=50)
    preddata <- expand.grid(preddata1,preddata2)
  } else {
    stop("smoothplot can only work with one or two numeric predictors in spline")
  }
  colnames(preddata) <- names(datamins)
  preddata <- as.data.frame(preddata)
  return(preddata)
}

.smooth_factorplot <- function(fscov, smoo, preddata_withlevels, summarise,
                               fs_fac, alpha, type, family, facet_args){
  if(length(fscov) == 1) {
    smoo$cov <- preddata_withlevels[[fscov]]
    smoo$group <- preddata_withlevels[[fs_fac]]

    smoo <- reshape2::melt(smoo, id =c("cov","group"))


    if(!is.null(summarise)){
      smoo <- stats::aggregate(smoo[,c("value")], by = list(smoo[["cov"]],smoo[["group"]]),
                               FUN = summarise)
      smoo$cov <- smoo$Group.1
      smoo$group <- smoo$Group.2
      smoo$value <- smoo$x
      smoo$variable <- "summary"
      alpha <- 1
    }
    if(type == "response"){
      if(family %in%
         c("poisson","lognormal","neg_binomial","zi_poisson",
           "zi_neg_binomial","gamma")){
        smoo$value <- exp(smoo$value)
      } else  if(family %in% c("binomial","bernoulli")){
        smoo$value <- inv_logit(smoo$value)
      }
    }


    # add some other labels
    term <- paste0(fscov, collapse = "__")
    facet_args$facets <- ~group

    # allpred <- rbind(allpred, smoo)
    # plot by term
    fs_pl <- ggplot2::ggplot(smoo) +
      ggplot2::geom_line(ggplot2::aes(x = .data[["cov"]],
                                      y = .data[["value"]],
                                      group = .data[["variable"]]),
                         alpha = alpha) +
      bayesplot::bayesplot_theme_get() +
      ggplot2::facet_wrap(facet_args) +
      ggplot2::labs(x = term)
  } else if (length(fscov) == 2) {
    smoo$cov1 <- preddata_withlevels[[fscov[1]]]
    smoo$cov2 <- preddata_withlevels[[fscov[2]]]
    smoo$group <- preddata_withlevels[[fs_fac]]

    smoo <- reshape2::melt(smoo, id =c("cov1","cov2","group"))

    if(!is.null(summarise)){
      smoo <- stats::aggregate(smoo[,c("value")], by = list(smoo[["cov1"]],
                                                            smoo[["cov2"]],
                                                            smoo[["group"]]),
                               FUN = summarise)
      smoo$cov1 <- smoo$Group.1
      smoo$cov2 <- smoo$Group.2
      smoo$group <- smoo$Group.3
      smoo$value <- smoo$x
      smoo$variable <- "summary"
    }


    # add some other labels
    term <- paste0(fscov, collapse = "__")

    # allpred <- rbind(allpred, smoo)
    # plot by term
    fs_pl <- ggplot2::ggplot(smoo) +
      ggplot2::geom_tile(ggplot2::aes(x = .data[["cov1"]],
                                      y = .data[["cov2"]],
                                      fill = .data[["value"]])) +
      bayesplot::bayesplot_theme_get() +
      ggplot2::labs(x = fscov[1],y=fscov[2]) +
      ggplot2::scale_fill_distiller(palette = "RdBu", type = "div")
    if(is.null(summarise)){
      facet_args$facets <- group~.data[["variable"]]
      fs_pl <- fs_pl +
        ggplot2::facet_grid(facet_args)
    } else {
      facet_args$facets <- ~group
      fs_pl <- fs_pl +
        ggplot2::facet_wrap(facet_args)
    }
  } else{
    stop("smoothplot can only work with one or two numeric predictors")
  }

  return(fs_pl)
}
