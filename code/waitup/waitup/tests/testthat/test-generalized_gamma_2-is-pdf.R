library(waitup)
context("generalized Gamma (2, Lawless 1992) pdf integrates")

test_that("generalized Gamma (Lawless 1992) density integrates", {
  expect_is_pdf(f=generalized_gamma_pdf_2, domain=c(0,Inf), 
    arg_bounds=list(k=c(1, 10), mu=c(3,5), sigma=c(1,2)), 
    n_arg_points=30)
})

test_that("generalized Gamma (Lawless 1992) exp(log density) integrates", {
  f <- function(x, ...) exp(generalized_gamma_lpdf_2(x, ...))
  expect_is_pdf(f=f, domain=c(0,Inf), arg_bounds=list(
    k=c(1, 10), mu=c(3, 5), sigma=c(1,2)), n_arg_points=30)
})

test_that("generalized Gamma (Stan, Lawless 1992) density integrates", {
  f <- function(x, ...) sapply(x, generalized_gamma_pdf_2S, ...)
  expect_is_pdf(f=f, domain=c(0,Inf), 
    arg_bounds=list(k=c(1, 10), mu=c(3,5), sigma=c(1,2)), 
    n_arg_points=30)
})

test_that("generalized Gamma (Stan, Lawless 1992) exp(log density) integrates", {
  f <- function(x, ...) sapply(x, function(x) exp(generalized_gamma_lpdf_2S(x, ...)))
  expect_is_pdf(f=f, domain=c(0,Inf), arg_bounds=list(
    k=c(1, 10), mu=c(3, 5), sigma=c(1,2)), n_arg_points=30)
})



