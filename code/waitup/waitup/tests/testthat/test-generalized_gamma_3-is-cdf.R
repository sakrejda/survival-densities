library(waitup)
context("generalized Gamma (3, Stacy 1962) cdf has cdf properties")
n_eps <- 5
n_points <- 10^2


test_that("generalized Gamma (Stacy 1962) cumulative density is cdf.", {
  expect_is_cdf(generalized_gamma_cdf_3, domain=c(0, 100), 
    arg_bounds=list(alpha=c(1,5), beta=c(1,5), nu=c(1,5)),
    n_arg_points=20)
})

test_that("generalized Gamma (Stacy 1962) exp(log cumulative) density is cdf.", {
  f <- function(x, ...) exp(generalized_gamma_lcdf_3(x, ...))
  expect_is_cdf(f, domain=c(0, 100), 
    arg_bounds=list(alpha=c(1,5), beta=c(1,5), nu=c(1,5)),
    n_arg_points=20)
})

