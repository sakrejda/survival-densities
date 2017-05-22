
## ---- gamma-exp-sum-gamma-mix

#' Naive implementation of the density of the mixture of a 
#' Gamma-Exponential sum with its Gamma component.
#' @param x Value to evaluate density at.
#' @param q Mixing parameter.
#' @param alpha Shape parameter.
#' @param beta Scale parameter.
#' @param delta Scale parameter.
gamma_exp_sum_gamma_mix_pdf_1 <- function(x, q, alpha, beta, delta) {
  dges <- gamma_exp_sum_pdf_1(x, alpha, beta, delta)
  dg <- dgamma(x=x, shape=alpha, scale=beta)
  d <- (1-q)*dg + q*dges
  return(d)
}

#' @describeIn gamma_exp_sum_gamma_mix_pdf_1
#' Naive implementation of the log density of the mixture of a 
#' Gamma-Exponential sum with its Gamma component.
gamma_exp_sum_gamma_mix_lpdf_1 <- function(x, q, alpha, beta, delta) {
  theta <- data.frame(x=x, q=q, alpha=alpha, beta=beta, delta=delta)
  N <- nrow(theta)
  dg <- dgamma(x=x, shape=alpha, scale=beta, log=TRUE)
  dges <- gamma_exp_sum_lpdf_1(theta[['x']], theta[['alpha']], theta[['beta']], theta[['delta']])
  d <- vector()
  for ( i in 1:N) {
    if (theta[i,'q'] != 1 && theta[i,'q'] != 0) {
      d[i] <- log_sum_exp(c(log(1-theta[i,'q']) + dg[i], 
                            log(  theta[i,'q']) + dges[i]))
    } else if (theta[i,'q'] == 0) {
      d[i] <- dg[i]
    } else if (1-theta[i,'q'] == 0) { 
      d[i] <- dges[i]
    }
  }
  return(d)
}

#' @describeIn gamma_exp_sum_gamma_mix_pdf_1
#' Naive implementation of the cumulative density of the mixture of a 
#' Gamma-Exponential sum with its Gamma component.
gamma_exp_sum_gamma_mix_cdf_1 <- function(x, q, alpha, beta, delta) {
  dg <- pgamma(q=x, shape=alpha, scale=beta)
  dges <- gamma_exp_sum_cdf_1(x, alpha, beta, delta)
  d <- (1-q)*dg + q*dges
  return(d)
}

#' @describeIn gamma_exp_sum_gamma_mix_pdf_1 Naive implementation of the
#' Weibull-uniform mixture log cumulative density.
gamma_exp_sum_gamma_mix_lcdf_1 <- function(x, q, alpha, beta, delta) {
  theta <- data.frame(x=x, q=q, alpha=alpha, beta=beta, delta=delta)
  N <- nrow(theta)
  dg <- pgamma(q=x, shape=alpha, scale=beta, log.p=TRUE)
  dges <- gamma_exp_sum_lcdf_1(x, alpha, beta)
  d <- vector()
  for ( i in 1:N ) {
    d[i] <- log_sum_exp(c(log(1-theta[i,'q'])+dges[i], 
                          log(  theta[i,'q'])-dg[i])) 
  }
  return(d)
}




