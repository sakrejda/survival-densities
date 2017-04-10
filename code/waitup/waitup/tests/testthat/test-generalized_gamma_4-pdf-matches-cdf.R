library(waitup)
context("generalized Gamma (4, Flexsurv) pdf-matches-cdf")
n_eps <- 5
n_points <- 10^2

test_that("generalized Gamma (Flexsurv) pdf matches generalized Gamma (Flexsurv) cdf.", {
  expect_pdf_matches_cdf(generalized_gamma_pdf_4, generalized_gamma_cdf_4,
    domain=c(.Machine$double.eps,100), arg_bounds=list(k=c(1.5, 5), 
      mu=c(1,10), sigma=c(1.5, 5)), n_arg_points=20)
})

test_that("generalized Gamma (Flexsurv) lpdf matches generalized Gamma (Flexsurv) lcdf.", {
  f <- function(x, ...) exp(generalized_gamma_lpdf_4(x, ...))
  g <- function(x, ...) exp(generalized_gamma_lcdf_4(x, ...))
  expect_pdf_matches_cdf(f, g,
    domain=c(.Machine$double.eps,100), arg_bounds=list(k=c(1.5, 5), 
      mu=c(1,10), sigma=c(1.5, 5)), n_arg_points=20)
})
