

## ---- log_gamma 

#' Density of the log generalized gamma directly from Lawless (1992)
#' @param x Value to evaluate density at.
#' @param k Shape parameter.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
generalized_log_gamma_pdf_1 <- function(x, k, mu, sigma) { 
  d <- (k^(k-.5))/sigma/gamma(k) * exp(
    sqrt(k)/sigma*(x-mu) - 
    k*exp((x-mu)/sigma/sqrt(k))
  )
  return(d)
}
