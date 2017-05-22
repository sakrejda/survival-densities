library(waitup)
context("pdf integrates")

test_that("generalized log Gamma density integrates", {
  expect_is_pdf(f=generalized_log_gamma_pdf_1, domain=c(-100,100),
    arg_bounds=list(k=c(1.5, 5), mu=c(-2,2), sigma=c(1.5,4)),
    n_arg_points=20)
})




