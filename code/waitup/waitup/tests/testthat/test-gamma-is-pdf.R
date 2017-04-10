library(waitup)
context("Gamma (1, naive) pdf integrates")

test_that("Gamma density integrates", {
  expect_is_pdf(f=gamma_pdf_1, domain=c(0,Inf), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})

test_that("Gamma exp(log density) integrates", {
  f <- function(x, ...) exp(gamma_lpdf_1(x, ...))
  expect_is_pdf(f=f, domain=c(0,Inf), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})

test_that("Gamma (Stan) density integrates", {
  f <- function(x, ...) sapply(x, gamma_pdf_1S, ...)
  expect_is_pdf(f=f, domain=c(0,Inf), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})

test_that("Gamma (Stan) exp(log density) integrates", {
  f <- function(x, ...) sapply(x, function(x) exp(gamma_lpdf_1S(x, ...)))
  expect_is_pdf(f=f, domain=c(0,Inf), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5)), n_arg_points=20)
})




