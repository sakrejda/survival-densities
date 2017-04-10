library(waitup)
context("generalized Gamma (2, Lawless 1992) pdf-matches-cdf")
n_eps <- 5
n_points <- 10^2


test_that("generalized Gamma (Lawless 1992) pdf matches generalized Gamma (Lawless 1992) cdf.", {
  expect_pdf_matches_cdf(generalized_gamma_pdf_2, generalized_gamma_cdf_2,
    domain=c(.Machine$double.eps,100), arg_bounds=list(k=c(1.5, 5), 
      mu=c(1,10), sigma=c(1.5, 5)), n_arg_points=20)
})

test_that("generalized Gamma (Lawless 1992) lpdf matches generalized Gamma (Lawless 1992) lcdf.", {
  f <- function(x, ...) exp(generalized_gamma_lpdf_2(x, ...))
  g <- function(x, ...) exp(generalized_gamma_lcdf_2(x, ...))
  expect_pdf_matches_cdf(f, g,
    domain=c(.Machine$double.eps,100), arg_bounds=list(k=c(1.5, 5), 
      mu=c(1,10), sigma=c(1.5, 5)), n_arg_points=20)
})

test_that("generalized Gamma (Stan, Lawless 1992) pdf matches generalized Gamma (Lawless 1992) cdf.", {
  f <- function(x, ...) sapply(x, generalized_gamma_pdf_2S, ...)
  g <- function(x, ...) sapply(x, generalized_gamma_cdf_2S, ...)
  expect_pdf_matches_cdf(f, g,
    domain=c(.Machine$double.eps,100), arg_bounds=list(k=c(1.5, 5), 
      mu=c(1,10), sigma=c(1.5, 5)), n_arg_points=20)
})

test_that("generalized Gamma (Stan, Lawless 1992) lpdf matches generalized Gamma (Lawless 1992) lcdf.", {
  f <- function(x, ...) sapply(x, function(x) exp(generalized_gamma_lpdf_2S(x, ...)))
  g <- function(x, ...) sapply(x, function(x) exp(generalized_gamma_lcdf_2S(x, ...)))
  expect_pdf_matches_cdf(f, g,
    domain=c(.Machine$double.eps,100), arg_bounds=list(k=c(1.5, 5), 
      mu=c(1,10), sigma=c(1.5, 5)), n_arg_points=20)
})




