library(waitup)
context("Gamma-Exponential sum (1, naive) pdf integrates")

test_that("Gamma-Exponential sum density integrates", {
  expect_is_pdf(f=gamma_exp_sum_pdf_1, domain=c(0,Inf), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5), delta=c(6, 10)), n_arg_points=20)
})

test_that("Gamma-Exponential sum exp(log density) integrates", {
  f <- function(x, ...) exp(gamma_exp_sum_lpdf_1(x, ...))
  expect_is_pdf(f=f, domain=c(0,Inf), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5), delta=c(6, 10)), n_arg_points=20)
})

test_that("Gamma-Exponential sum (Stan) density integrates", {
  f <- function(x, ...) sapply(x, gamma_exp_sum_pdf_1S, ...)
  expect_is_pdf(f=f, domain=c(0,Inf), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5), delta=c(6, 10)), n_arg_points=20)
})

test_that("Gamma-Exponential sum (Stan) exp(log density) integrates", {
  f <- function(x, ...) sapply(x, function(x) exp(gamma_exp_sum_lpdf_1S(x, ...)))
  expect_is_pdf(f=f, domain=c(0,Inf), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5), delta=c(6, 10)), n_arg_points=20)
})




