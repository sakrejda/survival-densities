library(waitup)
context("Weibull pdf matches cdf")
n_eps <- 5
n_points <- 10^2

test_that("Weibull pdf matches Weibull cdf.", {
  expect_pdf_matches_cdf(weibull_pdf_1, weibull_cdf_1,
    domain=c(0,100), arg_bounds=list(
      alpha=c(1.5, 5), beta=c(1,10)), n_arg_points=20)
})

test_that("Weibull lpdf matches Weibull lcdf.", {
  f <- function(x, ...) exp(weibull_lpdf_1(x, ...))
  g <- function(x, ...) exp(weibull_lcdf_1(x, ...))
  expect_pdf_matches_cdf(f, g,
    domain=c(0,100), arg_bounds=list(
      alpha=c(1.5, 5), beta=c(1,10)), n_arg_points=20)
})

test_that("Weibull (Stan) pdf matches Weibull (Stan) cdf.", {
  f <- function(x, ...) sapply(x, weibull_lpdf_1S, ...)
  g <- function(x, ...) sapply(x, weibull_lcdf_1S, ...)
  expect_pdf_matches_cdf(f, g,
    domain=c(0,100), arg_bounds=list(
      alpha=c(1.5, 5), beta=c(1,10)), n_arg_points=20)
})

test_that("Weibull (Stan) lpdf matches Weibull (Stan) lcdf.", {
  f <- function(x, ...) sapply(x, function(x) exp(weibull_lpdf_1S(x, ...)))
  g <- function(x, ...) sapply(x, function(x) exp(weibull_lcdf_1S(x, ...)))
  expect_pdf_matches_cdf(f, g,
    domain=c(0,100), arg_bounds=list(
      alpha=c(1.5, 5), beta=c(1,10)), n_arg_points=20)
})

