library(waitup)
context("Weibull (1, naive) cdf has cdf properties")
n_eps <- 5
n_points <- 10^2

test_that("Weibull cumulative density is cdf.", {
  expect_is_cdf(weibull_cdf_1, domain=c(0, 100), 
    arg_bounds=list(alpha=c(1,5), beta=c(1,5)),
    n_arg_points=20)
})

test_that("Weibull exp(log cumulative) density is cdf.", {
  f <- function(x, ...) exp(weibull_lcdf_1(x, ...))
  expect_is_cdf(f, domain=c(0, 100), 
    arg_bounds=list(alpha=c(1,5), beta=c(1,5)),
    n_arg_points=20)
})

test_that("Weibull (Stan) cumulative density is cdf.", {
  f <- function(x, ...) sapply(x, weibull_cdf_1S, ...)
  expect_is_cdf(f, domain=c(0, 100), 
    arg_bounds=list(alpha=c(1,5), beta=c(1,5)),
    n_arg_points=20)
})

test_that("Weibull (Stan) exp(log cumulative) density is cdf.", {
  f <- function(x, ...) sapply(x, function(x) exp(weibull_lcdf_1S(x, ...)))
  expect_is_cdf(f, domain=c(0, 100), 
    arg_bounds=list(alpha=c(1,5), beta=c(1,5)),
    n_arg_points=20)
})




