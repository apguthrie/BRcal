######################################################
#  External Functions                                #
######################################################

#' Linear Log Odds (LLO) Recalibration Function
#'
#' LLO-adjust predicted probabilities based on specified \eqn{\delta} and
#' \eqn{\gamma}.
#'
#' The Linear Log Odds (LLO) recalibration function can be written as
#' \deqn{c(x_i;\delta, \gamma) = \frac{\delta x_i^\gamma}{\delta x_i^\gamma +
#' (1-x_i)^\gamma}} where \eqn{x_i} is a predicted probability,
#' \eqn{\delta > 0} and \eqn{\gamma \in \mathbb{R}}.  Then \eqn{c(x_i;\delta,
#' \gamma)} is the corresponding LLO-adjusted probability that has been shifted
#' by \eqn{\delta} and scaled by \eqn{\gamma} on the log odds scale.  When
#' \eqn{\delta = \gamma = 1}, there is no shifting or scaling imposed on `x`.
#'
#' @param x a numeric vector of predicted probabilities of an event. Must only
#'   contain values in \[0,1\].
#' @param delta numeric, must be > 0, parameter \eqn{\delta} in LLO
#'   recalibration function.
#' @param gamma numeric, parameter \eqn{\gamma} in LLO recalibration function.
#' @return Vector of LLO-adjusted probabilities via specified \eqn{\delta} and
#'   \eqn{\gamma}.
#' @export
#'
#' @references Turner, B., Steyvers, M., Merkle, E., Budescu, D., and Wallsten,
#'   T. (2014) Forecast aggregation via recalibration, \emph{Machine Learning}
#'   95, 261–289.
#'
#'   Gonzalez, R., and Wu, G. (1999), On the shape of probability weighting
#'   function, \emph{Cognitive Psychology} 38, 129–66.
#'
#' @examples
#'
#' # Vector of probability predictions from 0 to 1
#' x1 <- seq(0, 1, by=0.1)
#' x1
#'
#' # LLO-adjusted predictions via delta = 2, gamma = 3
#' x1_llo23 <- LLO(x1, 2, 3)
#' x1_llo23  
#'
#' # LLO-adjusted predictions via delta = 1, gamma = 1
#' x1_llo11 <- LLO(x1, 1, 1)
#' x1_llo11  # no change
#'
#' # Create vector of 100 probability predictions
#' x2 <- runif(100)
#'
#' # LLO-adjust via delta = 2, gamma = 3
#' x2_llo23 <- LLO(x2, 2, 3)
#'
#' plot(x2, x2_llo23)
LLO <- function(x, delta, gamma){
  
  ##################
  #  Input Checks  #
  ##################
  
  # check input probs are valid
  x <- check_input_probs(x, "x")
  
  # check delta > 0 & numeric & size 1
  check_input_delta(delta)
  
  # check gamma in Reals & numeric & size 1
  check_input_gamma(gamma)
  
  ###################
  #  Function Code  #
  ###################
  
  x_llo <- LLO_internal(x=x, delta=delta, gamma=gamma)
  
  ###################
  #  Output Checks  #
  ###################
  
  # check if return vector contains nans
  if(!check_noNaNs(x_llo)) warning("LLO return value contains NaNs")
  
  # check if return vector contains +/- Inf - typically not possible
  if(!check_noInfs(x_llo)) warning("LLO return value contains +/-Inf")
  
  # check x are probabilities in [0,1] - typically not possible to break
  if(!check_probs(x_llo[!is.na(x_llo)])) warning("LLO return value contains values outside of [0,1]")
  
  
  return(x_llo)
}


#' Likelihood Ratio Test for Calibration
#'
#' Perform a likelihood ratio test for if calibration a set of probability
#' predictions, `x`, are well-calibrated given a corresponding set of binary
#' event outcomes, `y`. See Guthrie and Franck (2024).
#'
#' This likelihood ratio test is based on the following likelihood
#' \deqn{\pi(\mathbf{x}, \mathbf{y} | \delta, \gamma) = \prod_{i=1}^n
#' c(x_i;\delta, \gamma)^{y_i} \left[1-c(x_i;\delta, \gamma)\right]^{1-y_i}}
#' where \eqn{c(x_i; \delta, \gamma)} is the Linear in Log Odds
#' (\link[BRcal]{LLO}) function, \eqn{\delta>0} is the shift parameter on the
#' logs odds scale, and \eqn{\gamma \in \mathbb{R}} is the scale parameter on
#' the log odds scale.
#'
#' As \eqn{\delta = \gamma = 1} corresponds to no shift or scaling of
#' probabilities, i.e. `x` is well calibrated given corresponding outcomes `y`.
#' Thus the hypotheses for this test are as follows: \deqn{H_0: \delta = 1,
#' \gamma = 1 \text{  (Probabilities are well calibrated)}} \deqn{H_1: \delta
#' \neq 1 \text{ and/or } \gamma \neq 1 \text{  (Probabilities are poorly
#' calibrated)}.}
#'
#' The likelihood ratio test statistics for \eqn{H_0} is
#' \deqn{\lambda_{LR} = -2 log\left[\frac{\pi(\delta =1, \gamma=1|\mathbf{x},
#' \mathbf{y})}{\pi(\delta = \hat\delta_{MLE}, \gamma = \hat\gamma_{MLE}|
#' \mathbf{x},\mathbf{y})}\right]} where \eqn{\lambda_{LR}\stackrel{H_0}{\sim}{\chi^2_2}}
#' asymptotically under the null hypothesis \eqn{H_0}, and
#' \eqn{\hat{\delta}_{MLE}} and \eqn{\hat{\gamma}_{MLE}} are the maximum
#' likelihood estimates for \eqn{\delta} and \eqn{\gamma}.
#'
#' @inheritParams LLO
#' @param y a vector of outcomes corresponding to probabilities in `x`. Must
#'   only contain two unique values (one for "events" and one for "non-events").
#'   By default, this function expects a vector of 0s (non-events) and 1s
#'   (events).
#' @param optim_details Logical.  If `TRUE`, the list returned by \link[stats]{optim} when
#'   minimizing the negative log likelihood is also returned by this function.
#' @param event Value in `y` that represents an "event".  Default value is 1.
#' @param epsilon Amount by which probabilities are pushed away from 0 or 1
#'   boundary for numerical stability. If a value in `x` < `epsilon`, it will be
#'   replaced with `epsilon`.  If a value in `x` > `1-epsilon`, that value will
#'   be replaced with `1-epsilon`.
#' @param ... Additional arguments to be passed to \link[stats]{optim}.
#'
#' @return A list with the following attributes: \item{\code{test_stat}}{The
#'   test statistic \eqn{\lambda_{LR}} from the likelihood ratio test.}
#'   \item{\code{pval}}{The p-value from the likelihood ratio test.}
#'   \item{\code{MLEs}}{Maximum likelihood estimates for \eqn{\delta} and
#'   \eqn{\gamma}.}
#'   \item{\code{optim_details}}{If `optim_details = TRUE`, the list returned by
#'   \link[stats]{optim} when minimizing the negative log likelihood, includes convergence
#'   information, number of iterations, and achieved negative log likelihood
#'   value and MLEs.}
#'
#' @importFrom stats optim
#'
#' @export
#'
#' @references Guthrie, A. P., and Franck, C. T. (2024) Boldness-Recalibration
#'   for Binary Event Predictions, \emph{The American Statistician} 1-17.
#'
#' @examples
#' # Simulate 100 predicted probabilities
#' x <- runif(100)
#' # Simulated 100 binary event outcomes using `x`
#' y <- rbinom(100, 1, x)  # By construction, `x` is well calibrated.
#'
#' # Run the likelihood ratio test on `x` and `y`
#' llo_lrt(x, y, optim_details=FALSE)
#'
#' # Use optim_details = TRUE to see returned info from call to optim(),
#' # details useful for checking convergence
#' llo_lrt(x, y, optim_details=TRUE)  # no convergence problems in this example
#'
#' # Use different start value in `optim()` call, start at delta = 5, gamma = 5
#' llo_lrt(x, y, optim_details=TRUE, par=c(5,5))
#'
#' # Use `L-BFGS-B` instead of `Nelder-Mead`
#' llo_lrt(x, y, optim_details=TRUE, method = "L-BFGS-B")  # same result
#'
#' # What if events are defined by text instead of 0 or 1?
#' y2 <- ifelse(y==0, "Loss", "Win")
#' llo_lrt(x, y2, event="Win", optim_details=FALSE)  # same result
#'
#' # What if we're interested in the probability of loss instead of win?
#' x2 <- 1 - x
#' llo_lrt(x2, y2, event="Loss", optim_details=FALSE)
#'
#' # Push probabilities away from bounds by 0.000001
#' x3 <- c(runif(50, 0, 0.0001), runif(50, .9999, 1))
#' y3 <- rbinom(100, 1, 0.5)
#' llo_lrt(x3, y3, epsilon=0.000001)
#' 
llo_lrt <- function(x, y, event=1, optim_details=TRUE, 
                    epsilon=.Machine$double.eps,  ...){
  
  ##################
  #  Input Checks  #
  ##################
  
  # check x is vector, values in [0,1]
  x <- check_input_probs(x, name="x")
  
  # check y is vector, values are 0s or 1s or can be converted using event
  y <- check_input_outcomes(y, name="y", event=event)
  
  # check optim_details is logical
  if(!is.logical(optim_details) & !(optim_details %in% c(0,1))) stop("argument optim_details must be logical")
  
  # check x and y are the same length
  if(length(x) != length(y)) stop("x and y length differ")
  
  # check epsilon
  epsilon <- check_value01(epsilon, name="epsilon")
  
  ###################
  #  Function Code  #
  ###################
  
  # Numerator of test statistic
  top <- llo_lik(c(1,1), x, y, log = TRUE, epsilon=epsilon)
  
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
                    MLEs = est_params,
                    optim_details = optLRT)
  } else {
    results <- list(test_stat = test_stat,
                    pval = pval,
                    MLEs = est_params)
  }
  
  return(results)
}

#' Recalibration via Maximum Likelihood Estimates (MLEs)
#'
#' MLE recalibrate (i.e. LLO-adjust via \eqn{\hat{\delta}_{MLE}} and
#' \eqn{\hat{\gamma}_{MLE}} as specified in Guthrie and Franck (2024).
#'
#' Given a set of probability predictions `x`, the corresponding MLE
#' recalibrated set is \eqn{c(x; \hat{\delta}_{MLE}, \hat{\gamma}_{MLE})} (see
#' \link[BRcal]{LLO}).
#'
#' @inheritParams llo_lrt
#' @param probs_only Logical.  If `TRUE`, `mle_recal()` returns only the vector
#'   of MLE recalibrated probabilities.
#'
#' @return If `probs_only==TRUE`, `mle_recal()`returns a vector of MLE
#'   recalibrated probabilities.  Otherwise, `mle_recal()` returns a list with
#'   the following attributes:
#'   \item{\code{probs}}{The vector of MLE
#'   recalibrated probabilities.}
#'   \item{\code{MLEs}}{Maximum likelihood estimates for \eqn{\delta} and
#'   \eqn{\gamma}.}
#'   \item{\code{optim_details}}{If `optim_details = TRUE`, the list returned by
#'   \link[stats]{optim} when minimizing the negative log likelihood, includes convergence
#'   information, number of iterations, and achieved negative log likelihood
#'   value and MLEs.}
#'
#' @importFrom stats optim
#'
#' @export
#'
#' @references Guthrie, A. P., and Franck, C. T. (2024) Boldness-Recalibration
#'   for Binary Event Predictions, \emph{The American Statistician} 1-17.
#'
#' @examples
#' # Simulate 100 predicted probabilities
#' x <- runif(100)
#' # Simulated 100 binary event outcomes using `x` 
#' y <- rbinom(100, 1, x)
#' 
#' # MLE recalibrate `x`
#' mle_recal(x, y, optim_details=FALSE)  
#' 
#' # Just return the vector of MLE recalibrated probabilities
#' x_mle <- mle_recal(x, y, optim_details=FALSE, probs_only=TRUE)
#' x_mle
#' 
#' plot(x, x_mle)
#' 
#' # Use optim_details = TRUE to see returned info from call to optim(),
#' # details useful for checking convergence
#' mle_recal(x, y, optim_details=TRUE)  # no convergence problems in this example
#' 
#' # Use different start value in `optim()` call, start at delta = 5, gamma = 5
#' mle_recal(x, y, optim_details=TRUE, par=c(5,5))
#' 
#' # Use `L-BFGS-B` instead of `Nelder-Mead` 
#' mle_recal(x, y, optim_details=TRUE, method = "L-BFGS-B")  # same result
#' 
#' # What if events are defined by text instead of 0 or 1? 
#' y2 <- ifelse(y==0, "Loss", "Win")
#' mle_recal(x, y2, event="Win", optim_details=FALSE)  # same result
#' 
#' # What if we're interested in the probability of loss instead of win?
#' x2 <- 1 - x
#' mle_recal(x2, y2, event="Loss", optim_details=FALSE)
#' 
mle_recal <- function(x, y, probs_only=FALSE, event=1, optim_details=TRUE, ...){
  
  ##################
  #  Input Checks  #
  ##################
  
  # check x is vector, values in [0,1]
  x <- check_input_probs(x, name="x")
  
  # check y is vector, values are 0s or 1s or can be converted using event
  y <- check_input_outcomes(y, name="y", event=event)
  
  # check optim_details is logical
  if(!is.logical(optim_details) & !(optim_details %in% c(0,1))) stop("argument optim_details must be logical")
  
  # check probs_only is logical
  if(!is.logical(probs_only) & !(probs_only %in% c(0,1))) stop("argument probs_only must be logical")
  
  # check probs_only & optim_details are NOT both true
  if(probs_only & optim_details) warning("optim_details cannot be returned when probs_only is TRUE, print details instead")
  
  # check x and y are the same length
  if(length(x) != length(y)) stop("x and y length differ")
  
  ###################
  #  Function Code  #
  ###################
  
  val <- mle_recal_internal(x=x, y=y, probs_only=probs_only, optim_details=optim_details, ...)
  return(val)
}

######################################################
#  Internal Functions                                #
######################################################

LLO_internal <- function(x, delta, gamma){
  (delta * (x^gamma)) / ((delta * (x^gamma)) + ((1-x)^gamma))
}


mle_recal_internal <- function(x, y, probs_only=TRUE, optim_details = TRUE, ...) {
  
  
  ###################
  #  Function Code  #
  ###################
  
  optLRT <- llo_optim(x, y, ...)
  est_params <- optLRT$par
  new_probs <- LLO_internal(x=x, est_params[1], est_params[2])
  if(probs_only){
    if(optim_details) print(optLRT)
    return(new_probs)
  } else if(optim_details){
    return(list(probs = new_probs,
                MLEs = est_params,
                optim_details = optLRT))
  } else{
    return(list(probs = new_probs,
                MLEs = est_params))
  }
}

# NEED TO MAKE NOTE IN DOCUMENTATION ABOUT ROUNDING OFF
to_logit <- logit <- function(p, epsilon=.Machine$double.eps){
  
  # rounding off x's that are too close to zero or one
  p <- ifelse(p < 0 + epsilon, 0 + epsilon, p) 
  p <- ifelse(p > 1 - epsilon, 1 - epsilon, p)
  
  return(log(p/(1-p)))
}

to_prob <- function(x){
  return(exp(x) / (1 + exp(x)))
}

# Likelihood
llo_lik <- function(params, x, y, log = FALSE, neg = FALSE, tau = FALSE, epsilon=.Machine$double.eps){
  
  if(tau){
    llo <- LLO_internal(x = x, delta = exp(params[1]), gamma = params[2])
  } else {
    llo <- LLO_internal(x = x, delta = params[1], gamma = params[2])
  }
  
  # rounding off x's that are too close to zero or one
  llo <- ifelse(llo < 0 + epsilon, 0 + epsilon, llo) 
  llo <- ifelse(llo > 1 - epsilon, 1 - epsilon, llo)
  
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


llo_optim <- function(x, y, par=c(0.5,0.5), tau=TRUE, gr=nll_gradient, lower = 0, upper = Inf, ...){
  
  # print(lower)
  # print(par)
  # convert delta to tau
  # NEED HANDELING FOR TAU = FALSE BC BOUND ON DELTA!
  if(tau){
    par[1] <- log(par[1])
    lower[1] <- log(lower[1]) 
    upper[1] <- log(upper[1]) 
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
