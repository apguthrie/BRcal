
# also return br probs (add flag to toggle this?)
# add option to see nlopt printing vs our own printing
# NOTE THAT EVERYTHING IS PRINTED IN TERMS OF TAU (use a message?)

# suppress nloprt output by default, add in details that if choose to print, will be in terms of tau (if that's what the bakeoff says)
# otherwise report everything in terms of delta

#' Boldness Recalibration for Binary Events
#'
#' Perform Bayesian boldness-recalibration as specified in Guthrie and Franck
#' (2024).
#'
#'
#' NEED TO CITE NLOPTR & ALGS, mention default algs and such,
#'
#' When `tau=TRUE`, the optimization routine operates relative to \eqn{\tau =
#' log(\delta)} instead of \eqn{\delta}.  Specification of start location `x0`
#' and bounds `lb`, `ub` should still be specified in terms of \eqn{\delta}. The
#' `brcal` function will automatically convert from \eqn{\delta} to \eqn{\tau}.
#' The returned values in ...
#'
#' For more control over the optimization routine conducted by `nloptr()`, the
#' user may specify their own options via the `opts` argument. See the
#' documentation for `nloptr()` for full details. Note that any objective,
#' constraint, or gradient functions specified by the user will be overwritten
#' by those specified in this package.
#'
#' @inheritParams bayes_ms
#' @param t Desired level of calibration in \[0,1\].
#' @param tau Logical.  If `TRUE`, the optimization operates on \eqn{\tau =
#'   log(\delta)} instead of \eqn{\delta}. See details.
#' @param start_at_MLEs Logical. If `TRUE`, the optimizer will start at
#'   \eqn{x_0} = the maximum likelihood estimates for \eqn{\delta} and
#'   \eqn{\gamma}. Otherwise, the user must specify \eqn{x_0}.
#' @param x0 Vector with starting locations for \eqn{\delta} and \eqn{\gamma}.
#'   This argument is ignored when start_at_MLEs = TRUE.
#' @param lb Vector with lower bounds for \eqn{\delta} and \eqn{\gamma}. Use
#'   `-Inf` to indicate no lower bound.
#' @param ub Vector with upper bounds for \eqn{\delta} and \eqn{\gamma}. Use
#'   `Inf` to indicate no upper bound.
#' @param maxeval Value passed to `nloptr()` to stop optimization when the
#'   number of function evaluations exceeds `maxeval`.
#' @param maxtime Value passed to `nloptr()` to stop optimization when
#'   evaluation time (in seconds) exceeds `maxtime`.
#' @param xtol_rel_inner Value passed to `nloptr()` to stop the inner
#'   optimization routine when the parameter estimates for \eqn{\delta} and
#'   \eqn{\gamma} change by less than `xtol_rel_inner`.
#' @param xtol_rel_outer Value passed to `nloptr()` to stop the outer
#'   optimization routine when the parameter estimates for \eqn{\delta} and
#'   \eqn{\gamma} change by less than `xtol_rel_inner`.
#' @param print_level Value passed to `nloptr()` to control how much output is
#'   printed during optimization. Default is to print the most information
#'   allowable by `nloptr()`. Specify `0` to suppress all output.
#' @param opts List with options to be passed to `nloptr`.  See details.
#'
#' @return A list with the following attributes:
#'   \item{\code{nloptr}}{The list returned by `nloptr()` including convergence
#'   information, number of iterations, and more.}
#'   \item{\code{Pmc}}{The prior model probability for the calibrated model
#'   \eqn{M_c} specified in function call.}
#'   \item{\code{t}}{Desired level of
#'   calibration in \[0,1\] specified in function call.}
#'   \item{\code{BR_params}}{(100\eqn{*}t)% Boldness-recalibration estimates for \eqn{\delta} and
#'   \eqn{\gamma}.}
#'   \item{\code{sb}}{The Bayesian Information Criteria (BIC) for the
#'   calibrated model \eqn{M_c}.}
#'   \item{\code{probs}}{Vector of (100\eqn{*}t)% boldness-recalibrated
#'   probabilities.}
#' @export
#'
#' @references Guthrie, A. P., and Franck, C. T. (2024) Boldness-Recalibration
#'   for Binary Event Predictions. \emph{arxiv}.
#'
#' @examples
brcal <- function(x, y, t=0.95, Pmc=0.5, tau=FALSE, event=1,
                  start_at_MLEs=TRUE, x0=NULL,
                  lb=c(0.00001, -Inf), ub=c(Inf, Inf),
                  maxeval=500, maxtime=NULL,
                  xtol_rel_inner=1.0e-6,
                  xtol_rel_outer=1.0e-6,
                  print_level=3,
                  opts=NULL){


  ##################
  #  Input Checks  #
  ##################

  # check x is vector, values in [0,1]
  x <- check_input_probs(x, name="x")

  # check y is vector, values are 0s or 1s
  y <- check_input_outcomes(y, name="y", event=event)

  # check x and y are the same length
  if(length(x) != length(y)) stop("x and y length differ")

  # check t is valid calibration prob
  t <- check_input_probs(t, name="t")

  # check Pmc is valid prior model prob
  Pmc <- check_input_probs(Pmc, name="Pmc")

  # check tau is logical
  if(!is.logical(tau) & !(tau %in% c(0,1))){
    stop("argument tau must be logical")
  }

  # check start_at_MLEs is logical
  if(!is.logical(start_at_MLEs) & !(start_at_MLEs %in% c(0,1))){
    stop("argument start_at_MLEs must be logical")
  }

  # check x0 is specified if start_at_MLEs is not
  if(!start_at_MLEs){
    if(is.null(x0)) stop("must specify x0 when start_at_MLEs=FALSE")

    # check x0
    x0 <- check_input_params(x0, name="x0")
  }

  # check upper and lower bounds
  lb <- check_input_params(lb, name="lb")
  ub <- check_input_params(ub, name="ub")

  ###################
  #  Function Code  #
  ###################

  if(start_at_MLEs){
    bt <- bayes_ms_internal(x,y)
    x0 <- bt$MLEs
  }

  if(tau){
    x0[1] <- log(x0[1])
  }

  # print(missing(opts))

  if(missing(opts)){
    # print("inside opts ifelse")
    res <- nloptr::nloptr(x0 = x0,
                          eval_f = obj_f,
                          eval_grad_f = obj_grad_f,
                          eval_g_ineq = constr_g,
                          eval_jac_g_ineq = constr_grad_g,
                          lb = lb,
                          ub = ub,
                          opts = list(algorithm = "NLOPT_LD_AUGLAG",
                                      maxeval = maxeval,
                                      maxtime = maxtime,
                                      xtol_rel = xtol_rel_outer,
                                      print_level = print_level,
                                      local_opts = list(
                                        algorithm = "NLOPT_LD_SLSQP",
                                        lb = lb,
                                        ub = ub,
                                        eval_grad_f = obj_grad_f,
                                        eval_jac_g_ineq = constr_grad_g,
                                        xtol_rel = xtol_rel_inner)),
                          probs = x,
                          outs = y,
                          t = t,
                          tau = tau,
                          Pmc = Pmc)
  } else {

    if("eval_f" %in% names(opts$local_opts)){
      opts$local_opts$eval_f <- obj_f
    }
    if("eval_grad_f" %in% names(opts$local_opts)){
      opts$local_opts$eval_grad_f <- obj_grad_f
    }
    if("eval_g_ineq" %in% names(opts$local_opts)){
      opts$local_opts$eval_g_ineq <- constr_g
    }
    if("eval_jac_g_ineq" %in% names(opts$local_opts)){
      opts$local_opts$eval_jac_g_ineq <- constr_grad_g
    }

    res <- nloptr::nloptr(x0 = x0,
                          eval_f = obj_f,
                          eval_grad_f = obj_grad_f,
                          eval_g_ineq = constr_g,
                          eval_jac_g_ineq  = constr_grad_g,
                          lb = lb,
                          ub = ub,
                          opts = opts,
                          probs = x,
                          outs = y,
                          t = t,
                          tau = tau,
                          Pmc = Pmc)
  }


  if(tau){
    res$solution[1] <- exp(res$solution[1])
  }

  l <- list(nloptr = res,
            Pmc = Pmc,
            t=t,
            BR_params = c(res$solution[1], res$solution[2]),
            sb = -res$objective,
            probs = LLO_internal(x=x, res$solution[1], res$solution[2]))

  return(l)
}


######################################################
#  Internal Functions                                #
######################################################

obj_f <- function(x, probs, outs, t, tau, Pmc){
  if(tau)(
    x[1] <- exp(x[1])
  )
  # this assumes x is vector?
  probs_new <- LLO_internal(x=probs, x[1], x[2])
  return(-stats::sd(probs_new))
}

obj_grad_f <- function(x, probs, outs, t, tau, Pmc){
  if(tau)(
    x[1] <- exp(x[1])
  )

  n <- length(probs)
  probs_new <- LLO_internal(x=probs, x[1], x[2])
  sdp <- stats::sd(probs_new)
  meanp <- mean(probs_new)
  xmxg <- (probs * (1-probs))^x[2]
  denom <- (x[1] * probs^x[2] + (1-probs)^x[2])^2
  firstd <- xmxg/denom
  meanxmxg <- mean(firstd)
  pmmp <- (probs_new - meanp)
  dxmxgl <- x[1] * xmxg * log(probs / (1-probs))
  firstg <- dxmxgl/denom
  meang <- mean(firstg)
  innerd <- pmmp * (firstd - meanxmxg)
  innerg <- pmmp * (firstg - meang)

  grad_obj <- c((-1/((n-1) * sdp)) * sum(innerd),
                (-1/((n-1) * sdp)) * sum(innerg))

  return(grad_obj)
}

constr_g <- function(x, probs, outs, t, tau, Pmc){
  if(tau)(
    x[1] <- exp(x[1])
  )
  probs_new <- LLO_internal(x=probs, x[1], x[2])
  c1 <- bayes_ms_internal(x=probs_new, y=outs, Pmc=Pmc)$posterior_model_prob * -1 + t
  return(c1)
}

constr_grad_g <- function(x, probs, outs, t, tau, Pmc){
  if(tau)(
    x[1] <- exp(x[1])
  )

  n <- length(probs)
  probs_new <- LLO_internal(x=probs, x[1], x[2])
  bt <- bayes_ms_internal(x=probs_new, y=outs, Pmc = Pmc)
  pmp <- bt$posterior_model_prob
  pmp2 <- pmp^2
  dhat <- bt$MLEs[1]
  ghat <- bt$MLEs[2]

  pmps <- (pmp - pmp2)

  firstd <- (ghat - 1) * outs * (1/x[1])
  secondd <- (dhat * probs^(ghat * x[2])*ghat*x[1]^(ghat-1)) / (dhat * x[1]^ghat * probs^(ghat * x[2]) + (1-probs)^(ghat * x[2]))
  thirdd <- (probs^x[2]) / (x[1] * probs^(x[2]) + (1-probs)^x[2])

  firstg <- (ghat - 1) * (outs * log(probs) + (1-outs) * log(1-probs))
  secondg <- (ghat*(dhat * x[1]^ghat * probs^(x[2] * ghat) * log(probs) + (1-probs)^(x[2] * ghat)*log(1-probs))) / (dhat * x[1]^ghat * probs^(x[2]*ghat) + (1-probs)^(x[2]*ghat))
  thirdg <- (x[1] * probs^x[2] * log(probs) + (1-probs)^x[2] * log(1-probs)) / (x[1] * probs^x[2] + (1-probs)^x[2])

  grad_obj <- c(pmps * sum(firstd - secondd + thirdd),
                pmps * sum(firstg - secondg + thirdg))
  return(grad_obj)
}
