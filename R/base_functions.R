######################################################
#  External Functions                                #
######################################################

#' Linear Log Odds (LLO) Recalibration Function
#'
#' Description goes here.
#'
#' Details go here. NEED TO CITE TURNER/Gonzalez & Wu? Consider removing functionality for
#' lists/matrices (any other types?).
#'
#' @param x a numeric vector of probabilities to be LLO-adjusted. Must only
#'          contain values in \[0,1\].
#' @param delta numeric, must be > 0, parameter \eqn{\delta} in LLO
#'              recalibration function.
#' @param gamma numeric, parameter \eqn{\gamma} in LLO recalibration function.
#' @param ... Additional arguments for internal use only.
#' @return The LLO-adjusted vector of probabilities (ADD NOTATION FROM PAPER? OR KEEP IN DETAILS?)
#' @export
#'
#' @examples
LLO <- function(x, delta, gamma, ...){

  ##################
  #  Input Checks  #
  ##################
  if(!exists("input_checks_off")){ input_checks_off <- FALSE }
  if(!input_checks_off){
    # check input probs are valid
    x <- check_input_probs(x, "x")

    # check delta > 0 & numeric & size 1
    check_input_delta(delta)

    # check gamma in Reals & numeric & size 1
    check_input_gamma(gamma)
  }

  ###################
  #  Function Code  #
  ###################

  x_llo <- (delta * (x^gamma)) / ((delta * (x^gamma)) + ((1-x)^gamma))

  ###################
  #  Output Checks  #
  ###################
  if(!exists("output_checks_off")){ output_checks_off <- FALSE }
  if(!output_checks_off){
    # check if return vector contains nans
    if(!check_noNaNs(x_llo)) warning("LLO return value contains NaNs")

    # check if return vector contains +/- Inf - typically not possible
    if(!check_noInfs(x_llo)) warning("LLO return value contains +/-Inf")

    # check x are probabilities in [0,1] - typically not possible to break
    if(!check_probs(x_llo[!is.na(x_llo)])) warning("LLO return value contains values outside of [0,1]")
  }

  return(x_llo)
}


# Likelihood Ratio Test
llo_lrt <- LLO_LRT <- function(x, y, params = c(1,1), optim_details = FALSE,
                    start = c(0.5,0.5), lower = c(0.001, -5), upper = c(10,30),
                    ...){
  # print("LLO_LRT")
  # print(paste0(params, " LLO_LRT"))
  # check params, start, lower, upper are of right length, right values
  params <- check_input_params(params)
  start <- check_input_params(start, name="start")
  lower <- check_input_params(lower, name="lower")
  upper <- check_input_params(upper, name="upper")

  # check x is vector, values in [0,1]
  x <- check_input_probs(x, name="x")

  # check y is vector, values are 0s or 1s
  y <- check_input_outcomes(y, name="y")

  # check optim_details is logical
  if(!is.logical(optim_details) & !(optim_details %in% c(0,1))) stop("argument log must be logical")

  # check x and y are the same length
  if(length(x) != length(y)) stop("x and y length differ")

  top <- llo_lik(params, x, y, log = TRUE)

  optLRT <- llo_optim(x, y, lower, upper, start, ...)

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
  # print("LLO_LRT end")
  return(results)
}

mle_recal <- function(x, y, probs_only=TRUE, optim_details = FALSE,
                      start = c(0.5,0.5), ...) {
  optLRT <- llo_optim(x, y, lower, upper, start, ...)
  est_params <- optLRT$par
  new_probs <- LLO(x=x, est_params[1], est_params[2], input_checks=FALSE)
  if(probs_only){
    return(new_probs)
  } else {
    return(list(probs = new_probs,
                MLEs = est_params))
  }
}

######################################################
#  Internal Functions                                #
######################################################

to_logit <- logit <- function(p){

  # check input probs are valid
  p <- check_input_probs(p, "p")

  # better way to handle the rounding here? - check literature
  p <- ifelse(p < (10^(-300)), (10^(-300)), p)
  p <- ifelse(p > 0.9999999999999999, 0.9999999999999999, p)

  return(log(p/(1-p)))
}

to_prob <- function(x){
  return(exp(x) / (1 + exp(x)))
}

# Likelihood
# ADD OPTION TO HAVE ANY TWO LEVELS FOR y
llo_lik <- function(params, x, y, log = FALSE, neg = FALSE, tau = FALSE){

  ##################
  #  Input Checks  #
  ##################
  # print("llo_lik")
  # print(paste0(params, " llo_lik"))
  # print(paste0(tau, " llo_lik"))

  # check params are of right length, right values
  #params <- check_input_params(params, tau=tau)
  # DONT CHECK HERE< it gets run too many times

  # check x is vector, values in [0,1]
  x <- check_input_probs(x, name="x")

  # check y is vector, values are 0s or 1s
  y <- check_input_outcomes(y, name="y")

  # check x and y are the same length
  if(length(x) != length(y)) stop("x and y length differ")

  # check log, neg, and tau are logical - allow 0 or 1
  if(!is.logical(log) & !(log %in% c(0,1))) stop("argument log must be logical")
  if(!is.logical(neg) & !(neg %in% c(0,1))) stop("argument neg must be logical")
  if(!is.logical(tau) & !(tau %in% c(0,1))) stop("argument tau must be logical")

  ###################
  #  Function Code  #
  ###################

  # rounding off x's that are too close to zero or one
  x <- ifelse(x < (10^(-300)), (10^(-300)), x)
  x <- ifelse(x > 0.9999999999999999, 0.9999999999999999, x)

  if(tau){
    llo <- LLO(x = x, delta = exp(params[1]), gamma = params[2])
  } else {
    llo <- LLO(x = x, delta = params[1], gamma = params[2])
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
  # print("llo_lik end")
  return(result)
}


llo_optim <- function(x, y, start=c(0.5,0.5), tau=TRUE, gr=nll_gradient, ...){

  # convert delta to tau
  if(tau){
    start[1] <- log(start[1])
  }

   opt <- optim(par=start, fn=llo_lik,
               ...,
               x=x, y=y, neg = TRUE, log = TRUE, tau=tau)

  if(tau){
    opt$par[1] <- exp(opt$par[1])
  }

  return(opt)
}


prelec <- function(p, alpha, beta){

  ##################
  #  Input Checks  #
  ##################

  # check input probs are valid
  p <- check_input_probs(p, "p")

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

  # check if return vector contains Infs - typically not possible
  if(!check_noInfs(p_prelec)) warning("return value contains +/-Inf")

  return(p_prelec)
}


nll_gradient <- function(params, x, y, tau, log=TRUE, neg=TRUE){
  if(tau){
    params[1] <- exp(params[1])
    ddelta <- -sum(y - ((params[1]* x^params[2])/(params[1] * x^params[2] + (1-x)^params[2])))
  } else{
    ddelta <- -sum((y / params[1]) - ((x^params[2])/(params[1] * x^params[2] + (1-x)^params[2])))
  }
  dgamma <- -sum(y * log(x) + log(1-x) - ((params[1] * log(x) * x^params[2] + log(1-x) * (1-x)^params[2])/(params[1] * x^params[2] + (1-x)^params[2])) - y * log(1-x))
  return(c(ddelta, dgamma))
}
