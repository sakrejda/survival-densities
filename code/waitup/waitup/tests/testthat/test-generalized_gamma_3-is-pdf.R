library(waitup)
context("generalized Gamma (3, Stacy 1962) pdf integrates")

test_that("generalized Gamma (Stacy 1962) density integrates", {
  expect_is_pdf(f=generalized_gamma_pdf_3, domain=c(0,Inf), 
    arg_bounds=list(alpha=c(1, 10), beta=c(3,5), nu=c(1,2)), 
    n_arg_points=30)
})

test_that("generalized Gamma (Stacy 1962) exp(log density) integrates", {
  f <- function(x, ...) exp(generalized_gamma_lpdf_3(x, ...))
  expect_is_pdf(f=f, domain=c(0,Inf), arg_bounds=list(
    alpha=c(1, 10), beta=c(3, 5), nu=c(1,2)), n_arg_points=30)
})



