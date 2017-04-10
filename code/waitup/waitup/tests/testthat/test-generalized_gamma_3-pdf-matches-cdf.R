library(waitup)
context("generalized Gamma (3, Stacy 1962) pdf-matches-cdf")
n_eps <- 5
n_points <- 10^2


test_that("generalized Gamma (Stacy 1962) pdf matches generalized Gamma (Stacy 1962) cdf.", {
  expect_pdf_matches_cdf(generalized_gamma_pdf_3, generalized_gamma_cdf_3,
    domain=c(0,100), arg_bounds=list(alpha=c(1.5, 5), 
      beta=c(1,10), nu=c(1.5, 5)), n_arg_points=20)
})

test_that("generalized Gamma (Stacy 1962) lpdf matches generalized Gamma (Stacy 1962) lcdf.", {
  f <- function(x, ...) exp(generalized_gamma_lpdf_3(x, ...))
  g <- function(x, ...) exp(generalized_gamma_lcdf_3(x, ...))
  expect_pdf_matches_cdf(f, g,
    domain=c(0,100), arg_bounds=list(alpha=c(1.5, 5), 
      beta=c(1,10), nu=c(1.5, 5)), n_arg_points=20)
})

