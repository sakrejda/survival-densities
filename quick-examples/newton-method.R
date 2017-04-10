
newton_approx <- function(y, f, x0, k) {
  x <- x0
  g <- function(x) f(x) - y
  dgdx <- function(x,eps=10^-18) (g(x+eps)-g(x-eps))/(2*eps)
  for ( i in 0:(k-1)) {
    x <- x - g(x)/dgdx(x)
  }
  return(x)
}




