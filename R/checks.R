# check vector of probs
check_probs <- function(x){
  return(!any(x > 1 | x < 0))
}
