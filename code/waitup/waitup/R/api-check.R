check_is_not_zero_length <- function(x) {
  if (length(x) == 0) 
    stop("Expecting a non-zero length argument.")
}

check_is_numeric <- function(x) {
  if (is.numeric(x))
    return
  else 
    stop("Expecting a numeric argument.")
}





