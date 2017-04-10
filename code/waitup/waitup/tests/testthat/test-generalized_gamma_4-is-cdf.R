library(waitup)
context("generalized Gamma (4, Flexsurv) cdf has cdf properties")
n_eps <- 5
n_points <- 10^2


test_that("generalized Gamma (Flexsurv) cumulative density is cdf.", {
  expect_is_cdf(generalized_gamma_cdf_4, domain=c(0, 15000), 
    arg_bounds=list(k=c(1,5), mu=c(1,3), sigma=c(1,2)),
    n_arg_points=20)
})

test_that("generalized Gamma (Flexsurv) exp(log cumulative) density is cdf.", {
  f <- function(x, ...) exp(generalized_gamma_lcdf_4(x, ...))
  expect_is_cdf(f, domain=c(0, 15000), 
    arg_bounds=list(k=c(1,5), mu=c(1,3), sigma=c(1,2)),
    n_arg_points=20)
})



