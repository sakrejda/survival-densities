
## ---- gamma 

#' Naive implementation of the Gamma density.
#' @param x Value to evaluate density at.
#' @param alpha Shape parameter.
#' @param beta Scale parameter.
gamma_pdf_1 <- function(x, alpha, beta) {
  d <- 1 / beta / gamma(alpha) *
    (x / beta)^(alpha - 1) * 
    exp(-x / beta)
  return(d)
}

#' @describeIn gamma_pdf_1 Naive implementation of the Gamma log
#' density.
gamma_lpdf_1 <- function(x, alpha, beta) {
  d <- -log(beta) - lgamma(alpha) +
    (alpha-1)*(log(x) - log(beta)) - x/beta
  return(d)
}

#' @describeIn gamma_pdf_1 Naive implementation of the Gamma
#' cumulative density.
gamma_cdf_1 <- function(x, alpha, beta) {
  d <- pgamma(x/beta,  alpha) 
  return(d)
}

#' @describeIn gamma_pdf_1 Naive implementation of the Gamma log
#' cumulative density.
gamma_lcdf_1 <- function(x, alpha, beta) {
  d <- log(gamma_cdf_1(x, alpha, beta))
  return(d)
}

