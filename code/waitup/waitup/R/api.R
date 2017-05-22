#' @importFrom Rcpp evalCpp
#' @useDynLib waitup, .registration = TRUE
NULL


## ---- density_api 

#' API function for weibull density.
#' @param x Numeric value to evaluate the density at.
#' @param ... Parameterization-specific parameter parameter values.
#' @param method Length one character vector for name of
#'               parameterization.
#' @return Weibull density at x.
#' @export weibull_pdf
weibull_pdf <- function(x, ..., method='naive') {
  if (method == 'naive') 
    return(weibull_pdf_1(x, ...))
  else if (method == 'naive-stan')
    return(weibull_pdf_1S(x, ...))
  else 
    stop("method not found.")
}

#' API function for weibull log density.
#' @param x Numeric value to evaluate the log density at.
#' @param ... Parameterization-specific parameter parameter values.
#' @param method Length one character vector for name of
#'               parameterization.
#' @return Weibull log density at x.
#' @export weibull_lpdf
weibull_lpdf <- function(x, ..., method='naive') {
  if (method == 'naive') 
    return(weibull_lpdf_1(x, ...))
  else if (method == 'naive-stan')
    return(weibull_lpdf_1S(x, ...))
  else
    stop("method not found.")
}

#' @describeIn weibull_pdf API function for gamma density.
#' @export gamma_pdf
gamma_pdf <- function(x, ..., method='naive') {
  if (method == 'naive') 
    return(gamma_pdf_1(x, ...))
  if (method == 'naive-stan') 
    return(gamma_pdf_1S(x, ...))
  else
    stop("method not found.")
}

#' @describeIn weibull_lpdf API function for gamma log density.
#' @export gamma_lpdf
gamma_lpdf <- function(x, ..., method='naive') {
  if (method == 'naive') 
    return(gamma_lpdf_1(x, ...))
  if (method == 'naive-stan') 
    return(gamma_lpdf_1S(x, ...))
  else
    stop("method not found.")
}

#' @describeIn weibull_pdf API function for Weibull-uniform mixture density.
#' @export weibull_uniform_mixture_pdf
weibull_uniform_mixture_pdf <- function(x, ..., method='naive') {
  if (method == 'naive') 
    return(weibull_uniform_mixture_pdf_1(x, ...))
  if (method == 'naive-stan') 
    return(weibull_uniform_mixture_pdf_1(x, ...))
  else
    stop("method not found.")
}

#' @describeIn weibull_lpdf API function for Weibull-uniform mixture log density.
#' @export weibull_uniform_mixture_lpdf
weibull_uniform_mixture_lpdf <- function(x, ..., method='naive') {
  if (method == 'naive') 
    return(weibull_uniform_mixture_lpdf_1(x, ...))
  if (method == 'naive-stan') 
    return(weibull_uniform_mixture_lpdf_1S(x, ...))
  else
    stop("method not found.")
}

#' @describeIn weibull_pdf API function for generalized gamma density.
#' @export generalized_gamma_pdf
generalized_gamma_pdf <- function(x, ..., method='Lawless-1992') {
  if (method == 'naive') 
    return(generalized_gamma_pdf_1(x, ...))
  if (method == 'naive-stan') { 
    d <- sapply(x, generalized_gamma_pdf_1S(x, ...))
    return(d)
  }
  if (method == 'Lawless-1992') 
    return(generalized_gamma_pdf_2(x, ...))
  if (method == 'Lawless-1992-stan') {
    d <- sapply(x, generalized_gamma_pdf_2S(x, ...))
    return(d)
  }
  if (method == 'Stacy-1962') 
    return(generalized_gamma_pdf_3(x, ...))
  if (method == 'Flexsurv') 
    return(generalized_gamma_pdf_4(x, ...))
  else
    stop("method not found.")
}

#' @describeIn weibull_lpdf API function for generalized gamma log density.
#' @export generalized_gamma_lpdf
generalized_gamma_lpdf <- function(x, ..., method='Lawless-1992') {
  if (method == 'naive') 
    return(generalized_gamma_lpdf_1(x, ...))
  if (method == 'naive-stan') 
    return(generalized_gamma_lpdf_1S(x, ...))
  if (method == 'Lawless-1992') 
    return(generalized_gamma_lpdf_2(x, ...))
  if (method == 'Lawless-1992-stan') 
    return(generalized_gamma_lpdf_2S(x, ...))
  else
    stop("method not found.")
}

