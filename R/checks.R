
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


# Checks input probs are:
# - a vector (coerced to be if passed a list)
check_input_probs <- function(x, name="probs"){

  # check if p is list
  if(is.list(x)){
    warning("argument ", name, " is a list, coerced to vector")
    x <- unlist(x)
  }

  # check p is vector
  # if(!is.vector(x)) warning("argument ", name, " is ", class(x) ," type, not a
  #                           vector; unexpected behavior may occur")
  if(!is.vector(x)) stop("argument ", name, " is ", class(x) ," type, not a vector")

  # check p is numeric
  if(!is.numeric(x)) stop("argument ", name, " is not numeric type")

  # check p are probabilities in [0,1]
  if(!check_probs(x)) stop("argument ", name, " contains values outside of [0,1]")

  return(x)
}

check_input_outcomes <- function(y, name="outcomes"){

  # check if p is list
  if(is.list(y)){
    warning("argument ", name, " is a list, coerced to vector")
    y <- unlist(y)
  }

  # check p is vector
  # if(!is.vector(y)) warning("argument ", name, " is ", class(y) ," type, not a
  #                           vector; unexpected behavior may occur")
  if(!is.vector(y)) stop("argument ", name, " is ", class(y) ," type, not a vector")

  # check p is numeric
  if(!is.numeric(y)) stop("argument ", name, " is not numeric type")

  # check p are probabilities in [0,1]
  if(any(!(y %in% c(0,1)))) stop("argument ", name, " contains non 0 or 1 values")

  return(y)
}

check_input_params <- function(params, name="params", tau=FALSE){

  # print(paste0("params: delta = ", params[1], ", gamma = ", params[2], " in check_input_params"))

  # check if params is list
  if(is.list(params)){
    warning("argument params is a list, coerced to vector")
    params <- unlist(params)
  }

  # check p is vector
  if(!is.vector(params)) stop("argument ", name, " is ", class(params) ," type, not a vector")

  # check length
  if(length(params) != 2) stop("argument ", name, " must be length 2")

  if(tau){
    # check tau - use gamma function because same range
    check_input_gamma(params[2], name=paste0(name, "[1]"))
  } else {
    #check delta
    check_input_delta(params[1], name=paste0(name, "[1]"))
  }

  #check gamma
  check_input_gamma(params[2], name=paste0(name, "[2]"))

  return(params)
}

check_input_delta <- function(delta, name="delta"){
  # check delta > 0 & numeric & size 1
  if(length(delta) != 1) stop("argument ", name, " must be single value")
  if(!is.numeric(delta)) stop("argument ", name, " is not numeric type")
  if(delta <= 0) stop("argument ", name, " must be greater than 0")
}

check_input_gamma <- function(gamma, name="gamma"){
  # check gamma in Reals & numeric & size 1
  if(length(gamma) != 1) stop("argument ", name, " must be single value")
  if(!is.numeric(gamma)) stop("argument ", name, " is not numeric type")
}
