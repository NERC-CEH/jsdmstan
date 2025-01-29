#' jsdmStanFamily class
#'
#' This is the jsdmStanFamily class, which occupies a slot within any
#' jsdmStanFit object.
#'
#' @name jsdmStanFamily
#'
#' @section Elements for \code{jsdmStanFamily} objects:
#' \describe{
#'  \item{\code{family}}{
#'  A length one character vector describing family used to fit object. Options
#'  are \code{"gaussian"}, \code{"poisson"}, \code{"bernoulli"},
#'  \code{"neg_binomial"}, \code{"binomial"}, \code{"zi_poisson"},
#'  \code{"zi_neg_binomial"}, or \code{"multiple"}.
#'  }
#'  \item{\code{params}}{
#'  A character vector that includes all the names of the family-specific parameters.
#'  }
#'  \item{\code{params_dataresp}}{
#'  A character vector that includes any named family-specific parameters that are
#'  modelled in response to data.
#'  }
#'  \item{\code{preds}}{
#'  A character vector of the measured predictors included if family parameters
#'  are modelled in response to data. If family parameters are not modelled in
#'  response to data this is left empty.
#'  }
#'  \item{\code{data_list}}{
#'  A list containing the original data used to fit the model
#'   (empty when save_data is set to \code{FALSE} or family parameters are not
#'   modelled in response to data).
#'  }
#'  }
#'
NULL

jsdmStanFamily_empty <- function(){
  res <- list(family = character(),
              params = character(),
              params_dataresp= character(),
              preds = character(),
              data_list = list())
  class(res) <- "jsdmStanFamily"
  return(res)
}

# jsdmStanFamily methods

#' Print jsdmStanFamily object
#'
#' @param x A jsdmStanFamily object
#' @param ... Other arguments, not used at this stage.
#'
#' @export
print.jsdmStanFamily <- function(x, ...){
  cat(paste("Family:", x$family, "\n",
            ifelse(length(x$params)>0,
                   paste("With parameters:",
                         paste0(x$params, collapse = ", "),"\n"),
                   "")))
  if(length(x$params_dataresp)>0){
    pd_list <- lapply(seq_along(x$params_dataresp), function(i){
      paste("Family-specific parameter",
              paste0(x$params_dataresp[i],collapse=", "),
              "is modelled in response to", length(x$preds[[i]]),
              "predictors. These are named:",
              paste0(x$preds[[i]], collapse = ", "),"\n")
    }
    )
    cat(unlist(pd_list))

  }
}
