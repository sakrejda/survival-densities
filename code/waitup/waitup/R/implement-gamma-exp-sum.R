
## ---- gamma-exp-sum

#' Naive implementation of the Gamma-Exponential sum density.
#' @param x Value to evaluate density at.
#' @param alpha Shape parameter for the gamma component.
#' @param beta Scale parameter for the gamma component.
#' @param delta Scale parameter for the exponential component.
gamma_exp_sum_pdf_1 <- function(x, alpha, beta, delta) {
  d <- delta^(alpha-1) / (delta-beta)^alpha * 
    exp(-x/delta) * pgamma(q=(1/beta-1/delta)*x, shape=alpha)
  return(d)
}

#' @describeIn gamma_exp_sum_pdf_1 Naive implementation of the
#' Gamma-Exponential sum log density.
gamma_exp_sum_lpdf_1 <- function(x, alpha, beta, delta) {
  d <- (alpha-1)*log(delta) - alpha*log(delta-beta) - x/delta +
    pgamma(q=(1/beta-1/delta)*x, log.p=TRUE, shape=alpha)
  return(d)
}

#' @describeIn gamma_exp_sum_pdf_1 Naive implementation of the
#' Gamma-Exponential sum cumulative density.
gamma_exp_sum_cdf_1 <- function(x, alpha, beta, delta) {
  d <- pgamma(q=x/beta, shape=alpha) - 
    (delta/(delta-beta))^alpha * exp(-x/delta) *
    pgamma(q=(1/beta-1/delta)*x, shape=alpha)
  return(d)
}

#' @describeIn gamma_exp_sum_pdf_1 Naive implementation of the
#' Gamma-Exponential sum log cumulative density.
gamma_exp_sum_lcdf_1 <- function(x, alpha, beta, delta) {
  d <- log(gamma_exp_sum_cdf_1(x, alpha, beta, delta))
  return(d)
}

