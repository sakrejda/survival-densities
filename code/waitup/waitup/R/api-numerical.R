## ---- numerical_utilities 


#' log_sum_exp implementation as suggested by Bill Dunlap
#' and modified to deal with Inf.
#' @param x A vector of numbers to exp-sum-log
#' @return log-sum-exp of the x argument.
#' @export log_sum_exp
log_sum_exp <- function(x) { 
  check_is_not_zero_length(x)
  check_is_numeric(x)
  x <- x[x != -Inf]
  if (length(x) == 0)
    return(-Inf)
  x <- sort(x); n <- length(x)
  if (x[n] == Inf)
    return(Inf)
  o <- log1p(sum(exp(x[-n] - x[n]))) + x[n]
  return(o)
}




