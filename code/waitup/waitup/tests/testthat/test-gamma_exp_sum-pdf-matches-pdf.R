library(waitup)
library(flexsurv)
context("Gamma-Exponential sum pdf implementations match")


test_that("Gamma-Exponential sum (Stan) density matches R implementation", {
  f <- function(x, ...) sapply(x, gamma_exp_sum_pdf_1S, ...)
  g <- gamma_exp_sum_pdf_1
  expect_f_matches_g(f=f, g=g, domain=c(0,100), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5), delta=c(6,10)), n_arg_points=20)
})

test_that("Gamma-Exponential sum (Stan) exp(log density) matches R implementation", {
  f <- function(x, ...) sapply(x, function(x) gamma_exp_sum_lpdf_1S(x, ...))
  g <- gamma_exp_sum_lpdf_1
  expect_f_matches_g(f=f, g=g, domain=c(0,100), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5), delta=c(6,10)), n_arg_points=20)
})



