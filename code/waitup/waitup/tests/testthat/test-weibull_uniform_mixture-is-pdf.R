library(waitup)
context("Weibull-uniform mixture (1, naive) pdf integrates")

test_that("Weibull-uniform mixture density integrates", {
  expect_is_pdf(f=weibull_uniform_mixture_pdf_1, domain=c(0,100), 
    arg_bounds=list(max=c(100, 100), q=c(0.05, 0.95), 
      alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})

test_that("Weibull-uniform mixture exp(log density) integrates", {
  f <- function(x, ...) exp(weibull_uniform_mixture_lpdf_1(x, ...))
  expect_is_pdf(f=f, domain=c(0,100), arg_bounds=list(
    max=c(100,100), q=c(0.05, 0.95), alpha=c(1.5, 5), beta=c(1,5)), 
    n_arg_points=20)
})

test_that("Weibull-uniform mixture (Stan) density integrates", {
  f <- function(x, ...) sapply(x, weibull_uniform_mixture_pdf_1S, ...)
  expect_is_pdf(f=f, domain=c(0,100), 
    arg_bounds=list(max=c(100, 100), q=c(0.05, 0.95), 
      alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})

test_that("Weibull-uniform mixture (Stan) exp(log density) integrates", {
  f <- function(x, ...) sapply(x, function(x) exp(weibull_uniform_mixture_lpdf_1S(x, ...)))
  expect_is_pdf(f=f, domain=c(0,100), arg_bounds=list(
    max=c(100,100), q=c(0.05, 0.95), alpha=c(1.5, 5), beta=c(1,5)), 
    n_arg_points=20)
})


