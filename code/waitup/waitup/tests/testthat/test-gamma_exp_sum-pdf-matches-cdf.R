library(waitup)
context("Gamma-Exponential sum (1, naive) pdf-matches-cdf")
n_eps <- 5
n_points <- 10^2

test_that("Gamma-Exponential sum pdf matches gamma_exp_sum cdf.", {
  expect_pdf_matches_cdf(gamma_exp_sum_pdf_1, gamma_exp_sum_cdf_1,
    domain=c(0,100), arg_bounds=list(
      alpha=c(1.5, 5), beta=c(1,10), delta=c(20, 30)), n_arg_points=20)
})

test_that("Gamma-Exponential sum lpdf matches gamma_exp_sum lcdf.", {
  f <- function(x, ...) exp(gamma_exp_sum_lpdf_1(x, ...))
  g <- function(x, ...) exp(gamma_exp_sum_lcdf_1(x, ...))
  expect_pdf_matches_cdf(f, g,
    domain=c(0,100), arg_bounds=list(
      alpha=c(1.5, 5), beta=c(1,10), delta=c(20, 30)), n_arg_points=20)
})

test_that("Gamma-Exponential sum (Stan) pdf matches gamma_exp_sum (Stan) cdf.", {
  f <- function(x, ...) sapply(x, gamma_exp_sum_pdf_1S, ...)
  g <- function(x, ...) sapply(x, gamma_exp_sum_cdf_1S, ...)
  expect_pdf_matches_cdf(f, g,
    domain=c(0,100), arg_bounds=list(
      alpha=c(1.5, 5), beta=c(1,10), delta=c(20, 30)), n_arg_points=20)
})

test_that("Gamma-Exponential sum (Stan) lpdf matches gamma_exp_sum (Stan) lcdf.", {
  f <- function(x, ...) sapply(x, function(x) exp(gamma_exp_sum_lpdf_1S(x, ...)))
  g <- function(x, ...) sapply(x, function(x) exp(gamma_exp_sum_lcdf_1S(x, ...)))
  expect_pdf_matches_cdf(f, g,
    domain=c(0,100), arg_bounds=list(
      alpha=c(1.5, 5), beta=c(1,10), delta=c(20, 30)), n_arg_points=20)
})


