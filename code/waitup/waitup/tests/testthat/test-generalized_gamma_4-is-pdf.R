library(waitup)
context("generalized gamma (4, Flexsurv) pdf integrates")


test_that("generalized Gamma (Flexsurv) density integrates", {
  expect_is_pdf(f=generalized_gamma_pdf_4, domain=c(0,Inf), 
    arg_bounds=list(k=c(1, 10), mu=c(3,5), sigma=c(1,2)), 
    n_arg_points=30)
})

test_that("generalized Gamma (Flexsurv) exp(log density) integrates", {
  f <- function(x, ...) exp(generalized_gamma_lpdf_4(x, ...))
  expect_is_pdf(f=f, domain=c(0,Inf), arg_bounds=list(
    k=c(1, 10), mu=c(3, 5), sigma=c(1,2)), n_arg_points=30)
})





