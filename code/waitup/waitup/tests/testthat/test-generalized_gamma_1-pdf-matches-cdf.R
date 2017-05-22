library(waitup)
context("generalized Gamma (1, naive) pdf-matches-cdf")
n_eps <- 5
n_points <- 10^2

test_that("generalized Gamma pdf matches generalized Gamma cdf.", {
  expect_pdf_matches_cdf(generalized_gamma_pdf_1, generalized_gamma_cdf_1,
    domain=c(0,100), arg_bounds=list(alpha=c(1.5, 5), 
      beta=c(1,10), nu=c(1.5, 5)), n_arg_points=20)
})

test_that("generalized Gamma lpdf matches generalized Gamma lcdf.", {
  f <- function(x, ...) exp(generalized_gamma_lpdf_1(x, ...))
  g <- function(x, ...) exp(generalized_gamma_lcdf_1(x, ...))
  expect_pdf_matches_cdf(f, g,
    domain=c(0,100), arg_bounds=list(alpha=c(1.5, 5), 
      beta=c(1,10), nu=c(1.5, 5)), n_arg_points=20)
})

test_that("generalized Gamma (Stan) pdf matches generalized Gamma (Stan) cdf.", {
  f <- function(x, ...) sapply(x, generalized_gamma_pdf_1S, ...)
  g <- function(x, ...) sapply(x, generalized_gamma_cdf_1S, ...)
  expect_pdf_matches_cdf(f, g,
    domain=c(0,100), arg_bounds=list(alpha=c(1.5, 5), 
      beta=c(1,10), nu=c(1.5, 5)), n_arg_points=20)
})

test_that("generalized Gamma (Stan) lpdf matches generalized Gamma (Stan) lcdf.", {
  f <- function(x, ...) sapply(x, function(x) exp(generalized_gamma_lpdf_1S(x, ...)))
  g <- function(x, ...) sapply(x, function(x) exp(generalized_gamma_lcdf_1S(x, ...)))
  expect_pdf_matches_cdf(f, g,
    domain=c(0,100), arg_bounds=list(alpha=c(1.5, 5), 
      beta=c(1,10), nu=c(1.5, 5)), n_arg_points=20)
})




