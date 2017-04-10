
## ---- weibull_uniform_mixture 

#' Naive implementation of the Weibull-uniform mixture density.
#' @param x Value to evaluate density at.
#' @param max Maximum value x could take.
#' @param q Mixing parameter.
#' @param alpha Shape parameter.
#' @param beta Scale parameter.
weibull_uniform_mixture_pdf_1 <- function(x, max, q, alpha, beta) {
  dw <- weibull_pdf_1(x, alpha, beta)/weibull_cdf_1(max, alpha, beta)
  du <- 1/max
  d <- (1-q)*dw + q*du
  return(d)
}

#' @describeIn weibull_uniform_mixture_pdf_1 Naive implementation of the Weibull-uniform mixture log density.
weibull_uniform_mixture_lpdf_1 <- function(x, max, q, alpha, beta) {
  theta <- data.frame(x=x, max=max, q=q, alpha=alpha, beta=beta)
  N <- nrow(theta)
  dw <- weibull_lpdf_1(theta[['x']], theta[['alpha']], theta[['beta']]) -
        weibull_lcdf_1(theta[['max']], theta[['alpha']], theta[['beta']])
  du <- -log(theta[['max']])
  d <- vector()
  for ( i in 1:N) {
    d[i] <- log_sum_exp(c(log(1-theta[i,'q']) + dw[i], 
                          log(  theta[i,'q']) + du[i]))
  }
  return(d)
}

#' @describeIn weibull_uniform_mixture_pdf_1 Naive implementation of the
#' Weibull-uniform mixture cumulative density.
weibull_uniform_mixture_cdf_1 <- function(x, max, q, alpha, beta) {
  dw <- weibull_cdf_1(x, alpha, beta)
  du <- x/max
  d <- (1-q)*dw + q*du
  return(d)
}

#' @describeIn weibull_uniform_mixture_pdf_1 Naive implementation of the
#' Weibull-uniform mixture log cumulative density.
weibull_uniform_mixture_lcdf_1 <- function(x, max, q, alpha, beta) {
  theta <- data.frame(x=x, max=max, q=q, alpha=alpha, beta=beta)
  N <- nrow(theta)
  N <- pmax(length(x), length(max), length(q), length(alpha), length(beta))
  dw <- weibull_lcdf_1(x, alpha, beta)
  d <- vector()
  for ( i in 1:N ) {
    d[i] <- log_sum_exp(c(log(1-theta[i,'q'])+dw[i], 
      log(  theta[i,'q'])-log(theta[i,'max'])+log(theta[i,'x']))) 
  }
  return(d)
}




