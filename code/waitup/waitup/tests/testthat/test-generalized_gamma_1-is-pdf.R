library(waitup)
context("generalized Gamma (1, naive) pdf integrates")

test_that("generalized Gamma density integrates", {
  expect_is_pdf(f=generalized_gamma_pdf_1, domain=c(0,Inf), 
    arg_bounds=list(alpha=c(1.5, 5), beta=c(1,5), nu=c(1.5, 5)), 
    n_arg_points=20)
})

test_that("generalized Gamma exp(log density) integrates", {
  f <- function(x, ...) exp(generalized_gamma_lpdf_1(x, ...))
  expect_is_pdf(f=f, domain=c(0,Inf), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5), nu=c(1.5, 5)), n_arg_points=20)
})

test_that("generalized Gamma (Stan) density integrates", {
  f <- function(x, ...) sapply(x, generalized_gamma_pdf_1S, ...)
  expect_is_pdf(f=f, domain=c(0,Inf), 
    arg_bounds=list(alpha=c(1.5, 5), beta=c(1,5), nu=c(1.5, 5)), 
    n_arg_points=20)
})

test_that("generalized Gamma (Stan) exp(log density) integrates", {
  f <- function(x, ...) sapply(x, function(x) exp(generalized_gamma_lpdf_1S(x, ...)))
  expect_is_pdf(f=f, domain=c(0,Inf), arg_bounds=list(
    alpha=c(1.5, 5), beta=c(1,5), nu=c(1.5, 5)), n_arg_points=20)
})



