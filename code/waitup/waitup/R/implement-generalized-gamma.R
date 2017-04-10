
## ---- standard_generalized_gamma 

#' Naive implementation of the generalized Gamma density.
#' @param x Value to evaluate density at.
#' @param alpha Shape parameter.
#' @param beta Scale parameter.
#' @param nu Tail parameter.
generalized_gamma_pdf_1 <- function(x, alpha, beta, nu) {
  d <- nu/(beta*gamma(alpha/nu)) *
    (x/beta)^(alpha-1) * exp(-(x/beta)^nu)
  return(d)
}

#' @describeIn generalized_gamma_pdf_1 Naive implementation of the
#' generalized Gamma log density.
generalized_gamma_lpdf_1 <- function(x, alpha, beta, nu) {
  d <- log(nu) - log(beta) - lgamma(alpha/nu) +
    (alpha-1)*(log(x) - log(beta)) - (x/beta)^nu
  return(d) 
}

#' @describeIn generalized_gamma_pdf_1 Naive implementation of the
#' generalized Gamma cumulative density.
generalized_gamma_cdf_1 <- function(x, alpha, beta, nu) {
  d <- pgamma((x/beta)^nu, alpha/nu)
  return(d)
}

#' @describeIn generalized_gamma_pdf_1 Naive implementation of the
#' generalized Gamma log cumulative density.
generalized_gamma_lcdf_1 <- function(x, alpha, beta, nu) {
  d <- log(generalized_gamma_cdf_1(x, alpha, beta, nu))
  return(d)
}

## ---- lawless_generalized_gamma

#' Implementation of the generalized Gamma density from Lawless (1992)
#' @param x Value to evaluate density at.
#' @param k Shape parameter.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
generalized_gamma_pdf_2 <- function(x, k, mu, sigma) { 
  w <- (log(x)-mu)/sigma
  d <- (k^(k-.5))/sigma/gamma(k) * exp(
    sqrt(k)*w - k*exp(k^(-.5)*w)
  )/x
  return(d)
}

#' @describeIn generalized_gamma_pdf_2 Implementation of the generalized
#' Gamma log density from Lawless (1992)
generalized_gamma_lpdf_2 <- function(x, k, mu, sigma) {
  y <- log(x)
  w <- (y-mu)/sigma
  d <- (k-.5)*log(k) - log(sigma) - lgamma(k) +
    (sqrt(k)*w - k*exp(1/sqrt(k)*w)) - y
  return(d)
}

#' @describeIn generalized_gamma_pdf_2 Implementation of the generalized
#' Gamma cumulative density from Lawless (1992)
generalized_gamma_cdf_2 <- function(x, k, mu, sigma) {
   w <- (log(x) - mu)/sigma
   d <- pgamma(k*exp(1/sqrt(k)*w), k)
}

#' @describeIn generalized_gamma_pdf_2 Implementation of the generalized
#' Gamma cumulative log density from Lawless (1992)
generalized_gamma_lcdf_2 <- function(x, k, mu, sigma) {
   w <- (log(x) - mu)/sigma
   d <- pgamma(k*exp(1/sqrt(k)*w), k, log.p=TRUE)
}

## ---- stacy_generalized_gamma

#' Naive implementation of the generalized Gamma density from flexsurv
#' package (Stacy 1962)
#' @param x Value to evaluate density at.
#' @param alpha Shape parameter.
#' @param beta Scale parameter.
#' @param nu Tail parameter.
generalized_gamma_pdf_3 <- function(x, alpha, beta, nu) {
  d <- dgengamma.orig(x=x, shape=nu, scale=beta, k=alpha/nu)
  return(d)
}

#' @describeIn generalized_gamma_pdf_3 Naive implementation of the
#' generalized Gamma log density from flexsurv package (Stacy 1962)
generalized_gamma_lpdf_3 <- function(x, alpha, beta, nu) {
  d <- dgengamma.orig(x=x, shape=nu, scale=beta, k=alpha/nu, log=TRUE)
  return(d) 
}

#' @describeIn generalized_gamma_pdf_3 Naive implementation of the
#' generalized Gamma cumulative density from flexsurv package (Stacy 1962)
generalized_gamma_cdf_3 <- function(x, alpha, beta, nu) {
  d <- pgengamma.orig(q=x, shape=nu, scale=beta, k=alpha/nu)
  return(d) 
}

#' @describeIn generalized_gamma_pdf_3 Naive implementation of the
#' generalized Gamma log cumulative density from flexsurv package (Stacy 1962)
generalized_gamma_lcdf_3 <- function(x, alpha, beta, nu) {
  d <- pgengamma.orig(q=x, shape=nu, scale=beta, k=alpha/nu, log.p=TRUE)
  return(d) 
}

## ---- flexsurv_generalized_gamma

#' Implementation of the generalized Gamma density from R flexsurv
#' package 
#' @param x Value to evaluate density at.
#' @param k Shape parameter.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
generalized_gamma_pdf_4 <- function(x, k, mu, sigma) { 
  d <- dgengamma(x=x, mu=mu, sigma=sigma, Q=1/sqrt(k))
  return(d)
}

#' @describeIn generalized_gamma_pdf_4 Implementation of the generalized
#' Gamma log density from R flexsurv package.
generalized_gamma_lpdf_4 <- function(x, k, mu, sigma) {
  d <- dgengamma(x=x, mu=mu, sigma=sigma, Q=1/sqrt(k), log=TRUE)
  return(d)
}

#' @describeIn generalized_gamma_pdf_4 Implementation of the generalized
#' Gamma cumulative density from R flexsurv package.
generalized_gamma_cdf_4 <- function(x, k, mu, sigma) {
  d <- pgengamma(q=x, mu=mu, sigma=sigma, Q=1/sqrt(k))
  return(d)
}

#' @describeIn generalized_gamma_pdf_4 Implementation of the generalized
#' Gamma cumulative log density from R flexsurv package.
generalized_gamma_lcdf_4 <- function(x, k, mu, sigma) {
  d <- pgengamma(q=x, mu=mu, sigma=sigma, Q=1/sqrt(k), log.p=TRUE)
  return(d)
}
