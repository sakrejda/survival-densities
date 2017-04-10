## ---- test

#' Create a set of points in parameter space.
#' @param arg_bounds, n_arg_points.
#' @return data frame of points uniformly distributed in parameter
#'         space.
#' @export generate_parameter_points
generate_parameter_points <- function(arg_bounds, n_points) {
  n_points_per_dim <- n_points^(1/length(arg_bounds)) %>% ceiling
  o <- list()
  for (nom in names(arg_bounds)) {
    o[[nom]] <- runif(n=n_points_per_dim, min=arg_bounds[[nom]][1], 
      max=arg_bounds[[nom]][[2]])
  }
  points <- do.call(what=expand.grid, args=o)
  return(points)
}

#' Verify requirements of a pdf.
#' @param f the putative probability density function.
#' @param domain c(lower, upper) indicating the pdf's domain.
#' @param arg_bounds a list with one named element per argument to f,
#'        containing a 2-vector with the lower and upper bound for
#'        parameter values of f to consider.
#' @param n_arg_points An integer indicating how many points to generate
#'        within the parameter bounds.
#' @export expect_is_pdf
expect_is_pdf <- function(f, domain, arg_bounds, n_arg_points) {
  points <- generate_parameter_points(arg_bounds, n_arg_points)
  lb <- domain[1]
  ub <- domain[2]
  for ( r in 1:nrow(points)) {
    o  <- do.call(what=integrate, args=c(
      list(f=f, lower=lb, upper=ub, subdivisions=10^4), 
      as.list(points[r,]))
    )
    mass <- o$value
    error <- o$abs.error
    expect_lt(abs(1-mass), error^.5)
  }
}

#' Verify requirements of a cdf.
#' @param f the putative cumulative density function.
#' @param domain c(lower, upper) indicating the pdf's domain.
#' @param arg_bounds a list with one named element per argument to f,
#'        containing a 2-vector with the lower and upper bound for
#'        parameter values of f to consider.
#' @param n_arg_points An integer indicating how many points to generate
#'        within the parameter bounds.
#' @export expect_is_cdf
expect_is_cdf <- function(f, domain, arg_bounds, n_arg_points) {
  points <- generate_parameter_points(arg_bounds, n_arg_points)
  lb <- domain[1]
  ub <- domain[2]
  for ( r in 1:nrow(points)) {
    start <- do.call(what=f, args=c(list(lb), as.list(points[r,])))
    expect_lte(start, .Machine$double.eps^.25)
    stop <- do.call(what=f, args=c(list(ub), as.list(points[r,])))
    expect_lte(abs(stop-1),.Machine$double.eps^.25)
      for (x in seq(from=lb, to=ub, length.out=n_arg_points)) {
        y <- x+.Machine$double.eps^.5
        a <- do.call(what=f, args=c(list(x), as.list(points[r,])))
        b <- do.call(what=f, args=c(list(y), as.list(points[r,])))
        expect_gte(b-a, 0)
    }
  }
}



#' Verify a pdf matches a cdf.
#' @param pdf the putative probability density function.
#' @param cdf the putative cumulative density function.
#' @param domain c(lower, upper) indicating the pdf's domain.
#' @param arg_bounds a list with one named element per argument to f,
#'        containing a 2-vector with the lower and upper bound for
#'        parameter values of f to consider.
#' @param n_arg_points An integer indicating how many points to generate
#'        within the parameter bounds.
#' @export expect_pdf_matches_cdf
expect_pdf_matches_cdf <- function(pdf, cdf, domain, arg_bounds, n_arg_points) {
  points <- generate_parameter_points(arg_bounds, n_arg_points)
  lb <- domain[1]
  ub <- domain[2]
  for ( r in 1:nrow(points)) {
    for (x in seq(from=lb, to=ub, length.out=n_arg_points)) {
      a <- do.call(what=cdf, args=c(list(x), as.list(points[r,])))
      b <- do.call(what=integrate, args=c(
        list(f=pdf, lower=lb, upper=x, subdivisions=10^4), 
        as.list(points[r,]))
      )
      expect_lte(abs(a-b$value), 10^-4)
    }
  }
}

#' Verify a pdf matches a cdf.
#' @param f first pdf. 
#' @param g second pdf.
#' @param domain c(lower, upper) indicating the pdf's domain.
#' @param arg_bounds a list with one named element per argument to f,
#'        containing a 2-vector with the lower and upper bound for
#'        parameter values of f to consider.
#' @param n_arg_points An integer indicating how many points to generate
#'        within the parameter bounds.
#' @export expect_f_matches_g
expect_f_matches_g <- function(f, g, domain, arg_bounds, n_arg_points) {
  points <- generate_parameter_points(arg_bounds, n_arg_points)
  lb <- domain[1]
  ub <- domain[2]
  for ( r in 1:nrow(points)) {
    for (x in seq(from=lb, to=ub, length.out=n_arg_points)) {
      a <- do.call(what=f, args=c(list(x), as.list(points[r,])))
      b <- do.call(what=g, args=c(list(x), as.list(points[r,])))
      #expect_lte(abs(a-b), .Machine$double.eps^.25)
      expect_equal(a, b)
    }
  }
}
