library(waitup)
library(flexsurv)
context("Weibull pdf implementations match")

test_that("Weibull density matches R Weibull density", {
  g <- function(x, alpha, beta) dweibull(x=x, shape=alpha, scale=beta)
  expect_f_matches_g(f=weibull_pdf_1, g=g, domain=c(0,100), 
    arg_bounds=list(alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})

test_that("Weibull exp(log density) matches R Weibull density", {
  f <- function(x, ...) exp(weibull_lpdf_1(x, ...))
  g <- function(x, alpha, beta) dweibull(x=x, shape=alpha, scale=beta)
  expect_f_matches_g(f=f, g=g, domain=c(0,100), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})

test_that("Weibull (Stan) density matches R Weibull density", {
  f <- function(x, ...) sapply(x, weibull_pdf_1S, ...)
  g <- function(x, alpha, beta) dweibull(x=x, shape=alpha, scale=beta)
  expect_f_matches_g(f=f, g=g, domain=c(0,100), 
    arg_bounds=list(alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})

test_that("Weibull (Stan) exp(log density) matches R Weibull density", {
  f <- function(x, ...) sapply(x, function(x) exp(weibull_lpdf_1S(x, ...)))
  g <- function(x, alpha, beta) dweibull(x=x, shape=alpha, scale=beta)
  expect_f_matches_g(f=f, g=g, domain=c(0,100), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})







