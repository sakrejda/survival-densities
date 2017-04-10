library(waitup)
context("Gamma-Exponential sum (1, naive) cdf has cdf properties")
n_eps <- 5
n_points <- 10^2

test_that("Gamma-Exponential sum cumulative density is cdf.", {
  expect_is_cdf(gamma_exp_sum_cdf_1, domain=c(0, 200), 
    arg_bounds=list(alpha=c(1,5), beta=c(1,5), delta=c(6,20)),
    n_arg_points=20)
})


test_that("Gamma-Exponential sum exp(log cumulative) density is cdf.", {
  f <- function(x, ...) exp(gamma_exp_sum_lcdf_1(x, ...))
  expect_is_cdf(f, domain=c(0, 200), 
    arg_bounds=list(alpha=c(1,5), beta=c(1,5), delta=c(6,20)),
    n_arg_points=20)
})

test_that("Gamma-Exponential sum (Stan) cumulative density is cdf.", {
  f <- function(x, ...) sapply(x, gamma_exp_sum_cdf_1S, ...)
  expect_is_cdf(f, domain=c(0, 200), 
    arg_bounds=list(alpha=c(1,5), beta=c(1,5), delta=c(6,20)),
    n_arg_points=20)
})


test_that("Gamma-Exponential sum (Stan) exp(log cumulative) density is cdf.", {
  f <- function(x, ...) sapply(x, function(x) exp(gamma_exp_sum_lcdf_1S(x, ...)))
  expect_is_cdf(f, domain=c(0, 200), 
    arg_bounds=list(alpha=c(1,5), beta=c(1,5), delta=c(6,20)),
    n_arg_points=20)
})






