#' Plot the traceplots and density plots for parameters within a jsdmStanFit object
#'
#' This function takes parameters from the jsdmStanfit object and plots them using
#' the \code{mcmc_combo} function from the \code{bayesplot} package. These models
#' have a lot of parameters so selecting a subset is recommended. If pars is set to
#' \code{NULL} (the default) all parameters with either sigma or kappa in their name
#' will be plotted along with a random selection of the other parameters (total
#' number of other parameters set by \code{sample_n}).
#'
#' @param x The \code{jsdmStanFit} model object
#' @param pars The parameters to plot, by default a random sample of twenty of the
#'   parameters fit within the model
#' @param combo Which combination of plot types within \code{bayesplot::mcmc_combo}
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
#' @param ... Arguments passed to \code{bayesplot::mcmc_combo}
#'
#' @return An invisible list of the plots
#' @export
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
  model_pars <- par_sample(pars = pars, parnames = parnames, sample_n = sample_n)

  model_est <- rstan::extract(x$fit, pars = model_pars, permuted = FALSE,
                              inc_warmup = inc_warmup, include = include)
  variables <- dimnames(model_est)[[3]]

  if (plot) {
    default_ask <- devAskNewPage()
    on.exit(devAskNewPage(default_ask))
    devAskNewPage(ask = FALSE)
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
        devAskNewPage(ask = ask)
      }
    }
  }
  invisible(plots)
}


is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
  if(is.numeric(x)){
    abs(x - round(x)) < tol
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
#' \pkg{bayesplot} package
#'
#' @param x The jsdmStanFit model object
#' @param pars The parameters to plot, by default a random sample of twenty of the
#'   parameters fit within the model
#' @param type The MCMC plot type to be used
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
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object that can be customised using
#'   the \pkg{ggplot2} package
#' @export
#'
mcmc_plot.jsdmStanFit <- function(x, pars = NULL, type = "intervals",
                                  sample_n = 10, inc_warmup = FALSE, include = TRUE,
                                  regexp = FALSE, ...){
  valid_types <- gsub("^mcmc_","",as.character(bayesplot::available_mcmc()))
  if(!type %in% valid_types){
    stop(paste("Invalid plot type. Valid plot types are:",
               paste(valid_types, collapse = ", ")))
  }
  mcmc_fun <- get(paste0("mcmc_",type), asNamespace("bayesplot"))

  if(!is.wholenumber(sample_n) & is.null(pars))
    stop("If pars is NULL then sample_n must be a positive integer")
  parnames <- get_parnames(x)
  model_pars <- par_sample(pars = pars, parnames = parnames, sample_n = sample_n,
                           regexp = regexp)

  mcmc_arg_names <- names(formals(mcmc_fun))
  mcmc_args <- list(...)
  if ("x" %in% mcmc_arg_names) {
    if (grepl("^nuts_", type)) {
      # x refers to a molten data.frame of NUTS parameters
      mcmc_args$x <- nuts_params(x)
    } else {
      # x refers to a data.frame of draws
      draws <- extract(x, pars = model_pars, permuted = FALSE,
                       inc_warmup = inc_warmup, include = include)
      sel_variables <- dimnames(draws)[[3]]
      if (type %in% c("scatter", "hex") && length(sel_variables) != 2L) {
        stop(paste("Exactly 2 parameters must be selected for this type.",
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
  interval_type <- type %in% c("intervals", "areas")
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
