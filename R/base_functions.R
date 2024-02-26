######################################################
#  External Functions                                #
######################################################

#' Linear Log Odds (LLO) Recalibration Function
#'
#' Description goes here.
#'
#' Details go here. NEED TO CITE TURNER/Gonzalez & Wu? Consider removing
#' functionality for lists/matrices (any other types?).
#'
#' @param x a numeric vector of probabilities to be LLO-adjusted. Must only
#'   contain values in \[0,1\].
#' @param delta numeric, must be > 0, parameter \eqn{\delta} in LLO
#'   recalibration function.
#' @param gamma numeric, parameter \eqn{\gamma} in LLO recalibration function.
#' @param ... Additional arguments for internal use only.
#' @return The LLO-adjusted vector of probabilities (ADD NOTATION FROM PAPER? OR
#'   KEEP IN DETAILS?)
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


#' Likelihood Ratio Test for Calibration \[NEED MORE\]
#'
#' @inheritParams LLO
#' @param y a numeric vector of outcomes corresponding to probabilities in `x`.
#'   Must only contain 0 or 1.
#' @param optim_details Logical.  If `TRUE`, the list returned by `optim()` when
#'   minimizing the negative log likelihood is also returned by `llo_lrt()`.
#' @param ... Additional arguments to be passed to `optim()`.
#'
#' @return A list with the following attributes: \item{\code{test_stat}}{The
#'   test statistic from the likelihood ratio test formed as FILL IN}
#'   \item{\code{pval}}{The p-value from the likelihood ratio test.}
#'   \item{\code{mles}}{Maximum likelihood estimates for \eqn{\delta} and \eqn{\gamma}.}
#'   \item{\code{optim_details}}{If `optim_details = TRUE`, the list returned by
#'   `optim()` when minimizing the negative log likelihood, includes convergence
#'   information, number of iterations, and achieved negative log likelihood
#'   value and MLEs.}
#' @export
#'
#' @examples
llo_lrt <- function(x, y, optim_details = TRUE, ...){

  ##################
  #  Input Checks  #
  ##################
  if(!exists("input_checks_off")){ input_checks_off <- FALSE }
  if(!input_checks_off){

    # # check params, start, lower, upper are of right length, right values
    # params <- check_input_params(params)

    # check x is vector, values in [0,1]
    x <- check_input_probs(x, name="x")

    # check y is vector, values are 0s or 1s - Relax this?
    y <- check_input_outcomes(y, name="y")

    # check optim_details is logical
    if(!is.logical(optim_details) & !(optim_details %in% c(0,1))) stop("argument log must be logical")

    # check x and y are the same length
    if(length(x) != length(y)) stop("x and y length differ")
  }

  ###################
  #  Function Code  #
  ###################

  # Numerator of test statistic
  top <- llo_lik(c(1,1), x, y, log = TRUE)

  # Minimize log likelihood
  optLRT <- llo_optim(x, y, ...)

  # Denominator of test statistic
  bottom <- -optLRT$value

  # Extract MLEs
  est_params <- optLRT$par

  # Calc test statistic
  test_stat <- 2*(bottom-top)

  # Calc p-value
  pval <- 1-stats::pchisq(test_stat, 2)

  if(optim_details){
    results <- list(test_stat = test_stat,
                    pval = pval,
                    mles = est_params,
                    optim_details = optLRT)
  } else {
    results <- list(test_stat = test_stat,
                    pval = pval,
                    est_params = est_params)
  }

  return(results)
}

#' Title
#'
#' @inheritParams llo_lrt
#' @param probs_only blah
#'
#' @return
#' @export
#'
#' @examples
mle_recal <- function(x, y, probs_only=TRUE, optim_details = TRUE, ...) {
  optLRT <- llo_optim(x, y, lower, upper, ...)
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

# NEED TO MAKE NOT IN DOCUMENTATION ABOUT ROUNDING OFF
to_logit <- logit <- function(p){

  # check input probs are valid
  # p <- check_input_probs(p, "p")

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


llo_optim <- function(x, y, par=c(0.5,0.5), tau=TRUE, gr=nll_gradient, ...){

  # convert delta to tau
  # NEED HANDELING FOR TAU = FALSE BC BOUND ON DELTA!
  if(tau){
    par[1] <- log(par[1])
    if(exists("lower")){ lower[1] <- log(lower[1]) }
    if(exists("upper")){ upper[1] <- log(upper[1]) }
  }


  opt <- optim(par=par, fn=llo_lik,
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

# Gradient function for negative log likelihood
# add options for non log or non neg? why?
nll_gradient <- function(params, x, y, tau){
  if(tau){
    params[1] <- exp(params[1])
    ddelta <- -sum(y - ((params[1]* x^params[2])/
                          (params[1] * x^params[2] + (1-x)^params[2])))
  } else{
    ddelta <- -sum((y / params[1]) -
                     ((x^params[2])/(params[1] * x^params[2] + (1-x)^params[2])))
  }
  dgamma <- -sum(y * log(x) + log(1-x) -
                   ((params[1] * log(x) * x^params[2] + log(1-x) * (1-x)^params[2])/
                      (params[1] * x^params[2] + (1-x)^params[2])) - y * log(1-x))
  return(c(ddelta, dgamma))
}
