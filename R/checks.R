
# Returns true if all values in vector are in [0,1]
check_probs <- function(x){
  return(!any(x > 1 | x < 0))
}

# Returns true if there are no NaNs in vector
check_noNaNs <- function(x){
  return(!any(is.na(x)))
}

# Returns true if there are no Inf or -Inf in vector
check_noInfs <- function(x){
  return(!any(is.infinite(x)))
}
