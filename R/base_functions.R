#' Linear Log Odds (LLO) Recalibration Function
#'
#' Description goes here.
#'
#' Details go here. NEED TO CITE TURNER? Consider removing functionality for
#' lists/matrices (any other types?).
#'
#' @param p a numeric vector of probabilities to be LLO-adjusted. Must only
#'          contain values in \[0,1\].
#' @param delta numeric, must be > 0, parameter \eqn{\delta} in LLO
#'              recalibration function.
#' @param gamma numeric, parameter \eqn{\gamma} in LLO recalibration function.
#'
#' @return The LLO-adjusted vector of probabilities (ADD NOTATION FROM PAPER? OR KEEP IN DETAILS?)
#' @export
#'
#' @examples
LLO <- function(p, delta, gamma){

  ##################
  #  Input Checks  #
  ##################

  # check if p is list
  if(is.list(p)){
    warning("argument p is a list, will be coerced to vector")
    p <- unlist(p)
  }

  # check p is vector
  if(!is.vector(p)) warning("argument p is ", class(p) ," type, not a vector")

  # check p is numeric
  if(!is.numeric(p)) stop("argument p is not numeric type")

  # check p are probabilities in [0,1]
  if(!check_probs(p)) stop("argument p contains values outside of [0,1]")

  # check delta > 0 & numeric & size 1
  if(length(delta) != 1) stop("argument delta must be single value")
  if(!is.numeric(delta)) stop("argument delta is not numeric type")
  if(delta <= 0) stop("argument delta must be greater than 0")

  # check gamma in Reals & numeric & size 1
  if(length(gamma) != 1) stop("argument gamma must be single value")
  if(!is.numeric(gamma)) stop("argument gamma is not numeric type")

  ###################
  #  Function Code  #
  ###################

  p_llo <- (delta * (p^gamma)) / ((delta * (p^gamma)) + ((1-p)^gamma))


  ###################
  #  Output Checks  #
  ###################

  # check if return vector contains nans
  if(!check_noNaNs(p_llo)) warning("return value contains NaNs")

  # check p are probabilities in [0,1]
  if(!check_probs(p_llo[!is.na(p_llo)])) warning("return value contains values outside of [0,1]")

  return(p_llo)
}

#' Prelec Two Parameter Recalibration Function
#'
#' @param p a numeric vector of probabilities to be Prelec-adjusted. Must only
#'          contain values in \[0,1\].
#' @param alpha numeric, must be > 0, \eqn{\alpha} in Prelec two parameter function.
#' @param beta numeric, must be > 0, \eqn{\beta} in Prelec two parameter function.
#'
#' @return
#' @export
#'
#' @examples
prelec <- function(p, alpha, beta){

  ##################
  #  Input Checks  #
  ##################

  # check if p is list
  if(is.list(p)){
    warning("argument p is a list, will be coerced to vector")
    p <- unlist(p)
  }

  # check p is vector
  if(!is.vector(p)) warning("argument p is ", class(p) ," type, not a vector")

  # check p is numeric
  if(!is.numeric(p)) stop("argument p is not numeric type")

  # check p are probabilities in [0,1]
  if(!check_probs(p)) stop("argument p contains values outside of [0,1]")

  # check alpha > 0 & numeric & size 1
  if(length(alpha) != 1) stop("argument alpha must be single value")
  if(!is.numeric(alpha)) stop("argument alpha is not numeric type")
  if(alpha <= 0) stop("argument alpha must be greater than 0")

  # check beta > 0 & numeric & size 1
  if(length(beta) != 1) stop("argument beta must be single value")
  if(!is.numeric(beta)) stop("argument beta is not numeric type")
  if(beta <= 0) stop("argument beta must be greater than 0")

  ###################
  #  Function Code  #
  ###################

  p_prelec <- exp(-beta * ((-log(p))^alpha))

  ###################
  #  Output Checks  #
  ###################

  # check if return vector contains nans
  if(!check_noNaNs(p_prelec)) warning("return value contains NaNs")

  # check if return vector contains nans
  if(!check_noInfs(p_prelec)) warning("return value contains +/-Inf")

  return(p_prelec)
}




# Converts probs to logit scale
to_logit <- logit <- function(p){

  # check p are probabilities in [0,1]

  # better way to handle the rounding here? - check literature
  p <- ifelse(p < (10^(-300)), (10^(-300)), p)
  p <- ifelse(p > 0.9999999999999999, 0.9999999999999999, p)

  return(log(p/(1-p)))
}

to_prob <- function(x){
  return(exp(x)/(1+exp(x)))
}

# Likelihood
llo_lik <- function(params, x, y, log = FALSE, neg = FALSE, tau = FALSE){

  # check params are of right length, right values

  # check x's are between 0,1

  # check y's are 0s or 1s

  # check log & neg are logical

  # rounding off x's that are too close to zero or one
  x <- ifelse(x < (10^(-300)), (10^(-300)), x)
  x <- ifelse(x > 0.9999999999999999, 0.9999999999999999, x)

  if(tau){
    llo <- LLO(p = x, delta = exp(params[1]), gamma = params[2])
  } else {
    llo <- LLO(p = x, delta = params[1], gamma = params[2])
  }



  llo <- ifelse(llo < (10^(-300)), (10^(-300)), llo)
  llo <- ifelse(llo > 0.9999999999999999, 0.9999999999999999, llo)


  if(log){
    result <- sum(y * log(llo) + (1 - y) * log(1 - llo))
  } else{
    result <- prod((llo^y) * (1 - llo)^(1 - y))
  }

  if(neg){
    result <- -result
  }
  return(result)
}

llo_optim_wrap <- function(params, x, y, log = FALSE, neg = FALSE){
  if(params[1] <= 0){
    result <- -9999
    if(neg){
      result <- -result
    }
  } else {
    result <- llo_lik(params=params, x=x, y=y, log=log, neg=neg)
  }
  return(result)
}

llo_optim <- function(x, y, lower=c(0.0001, -15), upper=c(4e+08, 150), start=c(0.5,0.5), tau=TRUE){
  if(tau){
    lower[1] <- log(lower[1])
    upper[1] <- log(upper[1])
    start[1] <- log(start[1])
  }

  gradient <- function(params, x, y, tau, log=TRUE, neg=TRUE){
    if(tau){
      params[1] <- exp(params[1])
    }
    ddelta <- -sum((y / params[1]) - ((x^params[2])/(params[1] * x^params[2] + (1-x)^params[2])))
    dgamma <- -sum(y * log(x) + log(1-x) - ((params[1] * log(x) * x^params[2] + log(1-x) * (1-x)^params[2])/(params[1] * x^params[2] + (1-x)^params[2])) - y * log(1-x))
    return(c(ddelta, dgamma))
  }

  # opt <- optim(par=start, fn=llo_lik, gr=gradient,
  #              x=x, y=y,
  #              method = "L-BFGS-B",
  #              lower = lower, upper = upper, tau=tau, log = TRUE, neg = TRUE)
  opt <- optim(start, llo_lik, x=x, y=y, method = "Nelder-Mead",
                                                 neg = TRUE, log = TRUE, tau=tau)

  if(tau){
    opt$par[1] <- exp(opt$par[1])
  }
  return(opt)
}

# Likelihood Ratio Test
LLO_LRT <- function(x, y, params = c(1,1), optim_details = FALSE, start = c(0.5,0.5), lower = c(0.001, -5), upper = c(10,30)){

  # check params are of right length, right values

  # check x's are between 0,1

  # check y's are 0s or 1s

  # check optim details are logical

  # check start, lower, upper


  top <- llo_lik(params, x, y, log = TRUE)
  # optLRT <- stats::optim(start, llo_lik, x=x, y=y, method = "L-BFGS-B",
  #                 lower = lower, upper = upper, neg = TRUE, log = TRUE)
  # optLRT <- stats::optim(start, llo_optim_wrap, x=x, y=y, method = "Nelder-Mead",
  #                        neg = TRUE, log = TRUE)
  optLRT <- llo_optim(x,y,lower,upper,start)

  bottom <- -optLRT$value
  est_params <- optLRT$par
  val <- 2*(bottom-top)
  pval <- 1-stats::pchisq(val, 2)

  if(optim_details){
    results <- list(test_stat = val,
                    pval = pval,
                    est_params = est_params,
                    opt_value = bottom,
                    opt_counts = optLRT$counts,
                    opt_convergence = optLRT$convergence,
                    opt_message = optLRT$message)
  } else {
    results <- list(test_stat = val,
                    pval = pval,
                    est_params = est_params)
  }
  return(results)
}

# Likelihood Ratio Test
LLO_LRT_dev <- function(x, y, params = c(1,1), optim_details = FALSE, start = c(0.5,0.5), lower = c(0.001, -5), upper = c(10,30)){

  # check params are of right length, right values

  # check x's are between 0,1

  # check y's are 0s or 1s

  # check optim details are logical

  # check start, lower, upper


  top <- llo_lik(params, x, y, log = TRUE)
  optLRT <- stats::optim(start, llo_optim_wrap, x=x, y=y, method = "Nelder-Mead",
                         neg = TRUE, log = TRUE)
  bottom <- -optLRT$value
  est_params <- optLRT$par
  val <- 2*(bottom-top)
  pval <- 1-stats::pchisq(val, 2)

  if(optim_details){
    results <- list(test_stat = val,
                    pval = pval,
                    est_params = est_params,
                    opt_value = bottom,
                    opt_counts = optLRT$counts,
                    opt_convergence = optLRT$convergence,
                    opt_message = optLRT$message)
  } else {
    results <- list(test_stat = val,
                    pval = pval,
                    est_params = est_params)
  }
  return(results)
}


