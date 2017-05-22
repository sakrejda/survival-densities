library(waitup)
library(flexsurv)
context("Gamma pdf implementations match")

test_that("Gamma density matches R gamma density", {
  g <- function(x, alpha, beta) dgamma(x=x, shape=alpha, scale=beta)
  expect_f_matches_g(f=gamma_pdf_1, g=g, domain=c(0,100), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})

test_that("Gamma exp(log density) matches R gamma density", {
  f <- function(x, ...) exp(gamma_lpdf_1(x, ...))
  g <- function(x, alpha, beta) dgamma(x=x, shape=alpha, scale=beta)
  expect_f_matches_g(f=f, g=g, domain=c(0,100), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})

test_that("Gamma (Stan) density matches R gamma density", {
  f <- function(x, ...) sapply(x, gamma_pdf_1S, ...)
  g <- function(x, alpha, beta) dgamma(x=x, shape=alpha, scale=beta)
  expect_f_matches_g(f=f, g=g, domain=c(0,100), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})

test_that("Gamma (Stan) exp(log density) matches R gamma density", {
  f <- function(x, ...) sapply(x, function(x) exp(gamma_lpdf_1S(x, ...)))
  g <- function(x, alpha, beta) dgamma(x=x, shape=alpha, scale=beta)
  expect_f_matches_g(f=f, g=g, domain=c(0,100), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})



