
## ---- weibull 

#' Naive implementation of the Weibull density.
#' @param x Value to evaluate density at.
#' @param alpha Shape parameter.
#' @param beta Scale parameter.
weibull_pdf_1 <- function(x, alpha, beta) {
  d <- alpha/beta*(x/beta)^(alpha-1) *
    exp(-(x/beta)^alpha)
  return(d)
}

#' @describeIn weibull_pdf_1 Naive implementation of the Weibull log-density.
weibull_lpdf_1 <- function(x, alpha, beta) {
  d <- log(alpha) - log(beta) +
    (alpha-1) * (log(x)-log(beta)) -
    (x/beta)^alpha
  return(d)
}

#' @describeIn weibull_pdf_1 Naive implementation of Weibull cumulative density.
weibull_cdf_1 <- function(x, alpha, beta) {
  d <- 1 - exp(-(x/beta)^(alpha))
  return(d)
}

#' @describeIn weibull_pdf_1 Naive implementation of Weibull log cumulative density.
weibull_lcdf_1 <- function(x, alpha, beta) {
  d <- log(weibull_cdf_1(x, alpha, beta))
  return(d)
}

