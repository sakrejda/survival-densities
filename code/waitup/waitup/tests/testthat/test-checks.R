library(waitup)
context("assertions assert as intended")

test_that("Zero length argument check.", {
  start <- tryCatch(
    expr=check_is_not_zero_length(numeric()),
    error=function(e) return(e)
  )
  expect_equal(class(start)[1], 'simpleError')
  expect_equal(start$message, 'Expecting a non-zero length argument.')
})

test_that("Numeric argument check.", {
  start <- tryCatch(
    expr=check_is_numeric("3"),
    error=function(e) return(e)
  )
  expect_equal(class(start)[1], 'simpleError')
  expect_equal(start$message, 'Expecting a numeric argument.')
  start <- tryCatch(
    expr=check_is_numeric(list(a=1)),
    error=function(e) return(e)
  )
  expect_equal(class(start)[1], 'simpleError')
  expect_equal(start$message, 'Expecting a numeric argument.')
})

test_that("R dnorm is a pdf.", {
  arg_bounds <- list(mean=c(-5, 5), sd=c(.1, 3))
  o <- expect_is_pdf(dnorm, c(-Inf, Inf), arg_bounds, n_arg_points=10)
})

test_that("R pnorm is a cdf.", {
  arg_bounds <- list(mean=c(-5, 5), sd=c(.1, 3))
  o <- expect_is_cdf(pnorm, c(-30, 20), arg_bounds, n_arg_points=10)
})

test_that("R dnorm (pdf) matches pnorm (cdf).", {
  arg_bounds <- list(mean=c(-5, 5), sd=c(.1, 3))
  o <- expect_pdf_matches_cdf(dnorm, pnorm, c(-30, 20), 
                              arg_bounds, n_arg_points=10)
})





