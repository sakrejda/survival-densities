library(waitup)
context("Weibull-uniform mixture pdf-matches-cdf")
n_eps <- 5
n_points <- 10^2

test_that("Weibull-uniform mixture pdf matches Weibull-uniform mixture cdf.", {
  expect_pdf_matches_cdf(
    weibull_uniform_mixture_pdf_1, weibull_uniform_mixture_cdf_1,
    domain=c(0,100), arg_bounds=list(max=c(100, 100), q=c(.05, .95),
      alpha=c(1.5, 5), beta=c(1,10)), n_arg_points=20)
})

test_that("Weibull-uniform mixture lpdf matches Weibull-uniform mixture lcdf.", {
  f <- function(x, ...) exp(weibull_uniform_mixture_lpdf_1(x, ...))
  g <- function(x, ...) exp(weibull_uniform_mixture_lcdf_1(x, ...))
  expect_pdf_matches_cdf(f, g,
    domain=c(0,100), arg_bounds=list(max=c(100,100), q=c(.05, .95),
      alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})

test_that("Weibull-uniform mixture (Stan) pdf matches Weibull-uniform mixture (Stan) cdf.", {
  f <- function(x, ...) sapply(x, weibull_uniform_mixture_pdf_1S, ...)
  g <- function(x, ...) sapply(x, weibull_uniform_mixture_cdf_1S, ...)
  expect_pdf_matches_cdf(f, g,
    domain=c(0,100), arg_bounds=list(max=c(100, 100), q=c(.05, .95),
      alpha=c(1.5, 5), beta=c(1,10)), n_arg_points=20)
})

test_that("Weibull-uniform mixture (Stan) lpdf matches Weibull-uniform mixture (Stan) lcdf.", {
  f <- function(x, ...) sapply(x, 
    function(x) exp(weibull_uniform_mixture_lpdf_1S(x, ...)))
  g <- function(x, ...) sapply(x, 
    function(x) exp(weibull_uniform_mixture_lcdf_1S(x, ...)))
  expect_pdf_matches_cdf(f, g,
    domain=c(0,100), arg_bounds=list(max=c(100,100), q=c(.05, .95),
      alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})


