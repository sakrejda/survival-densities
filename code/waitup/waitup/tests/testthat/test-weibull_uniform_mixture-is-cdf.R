library(waitup)
context("Weibull-uniform mixture (1, naive) cdf has cdf properties")
n_eps <- 5
n_points <- 10^2


test_that("Weibull-uniform mixture cumulative density is cdf.", {
  expect_is_cdf(weibull_uniform_mixture_cdf_1, domain=c(0, 100), 
    arg_bounds=list(max=c(100,100), q=c(0.05, 0.95), alpha=c(1,5), beta=c(1,5)),
    n_arg_points=20)
})

test_that("Weibull-uniform mixture exp(log cumulative) density is cdf.", {
  f <- function(x, ...) exp(weibull_uniform_mixture_lcdf_1(x, ...))
  expect_is_cdf(f, domain=c(0, 100), arg_bounds=list(
    max=c(100,100), q=c(0.05,0.95), alpha=c(1,5), beta=c(1,5)),
    n_arg_points=20)
})

test_that("Weibull-uniform mixture (Stan) cumulative density is cdf.", {
  f <- function(x, ...) sapply(x, weibull_uniform_mixture_cdf_1S, ...)
  expect_is_cdf(f, domain=c(0, 100), 
    arg_bounds=list(max=c(100,100), q=c(0.05, 0.95), alpha=c(1,5), beta=c(1,5)),
    n_arg_points=20)
})

test_that("Weibull-uniform mixture (Stan) exp(log cumulative) density is cdf.", {
  f <- function(x, ...) sapply(x, function(x) exp(weibull_uniform_mixture_lcdf_1S(x, ...)))
  expect_is_cdf(f, domain=c(0, 100), arg_bounds=list(
    max=c(100,100), q=c(0.05,0.95), alpha=c(1,5), beta=c(1,5)),
    n_arg_points=20)
})


