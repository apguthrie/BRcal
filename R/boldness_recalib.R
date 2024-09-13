######################################################
#  External Functions                                #
######################################################

#' Boldness-Recalibration for Binary Events
#'
#' Perform Bayesian boldness-recalibration as specified in Guthrie and Franck
#' (2024). Boldness-recalibration maximizes the spread in predictions (`x`)
#' subject to a constraint on the minimum tolerable posterior probability of
#' calibration (`t`).
#'
#' The objective function in boldness-recalibration is \deqn{ f(\delta, \gamma)
#' = -sd(\mathbf{x}')} and the constraint is \deqn{g(\delta, \gamma) =
#' -(P(M_c|\mathbf{y}, \mathbf{x}')-t) \leq 0.}  As both the objective and
#' constraint functions are non-linear with respect to \eqn{\delta} and
#' \eqn{\gamma}, we use \link[nloptr]{nloptr} for this optimization rather than
#' \link[stats]{optim}. Note that we use `x` to denote a vector of predicted
#' probabilities, `nloptr()` uses `x` to denote the parameters being optimized.
#' Thus, starting values for \eqn{\delta} and \eqn{\gamma} are passed via
#' argument `x0` and all output refers to the objective and constraint as `f(x)`
#' and `g(x)`.
#'
#' By default, this function uses the Augmented Lagrangian Algorithm (AUGLAG)
#' (Conn et. al. 1991, Birgin and Martinez 2008) as the outer optimization
#' routine and Sequential Least-Squares Quadratic Programming (SLSQP) (Dieter
#' 1988, Dieter 1994) as the inner optimization routine.
#'
#' @section Adjusting call to `nloptr()`:
#'
#' For more control over the optimization routine conducted by `nloptr()`, the
#' user may specify their own options via the `opts` argument.  Note that any
#' objective, constraint, or gradient functions specified by the user will be
#' overwritten by those specified in this package. See the documentation for
#' `nloptr()` and the NLopt website for full details
#' (<https://nlopt.readthedocs.io/en/latest/>).
#'
#' @section Adjusting call to `optim()`:
#'
#' While `optim()` is not used for the non-linear constrained optimization for
#' finding he boldness-recalibration parameters, it is used in the constraint
#' function as it involves the posterior model posterior.  Because of this,
#' we do allow users to pass additional arguments to `optim` to be used in this
#' calculation.  However, rather than use the `...`,
#' users should pass these arguments to `optim_options` via a list.
#'
#' @section Optimizing over \eqn{\tau}:
#'
#' When `tau=TRUE`, the optimization routine operates relative to \eqn{\tau =
#' log(\delta)} instead of \eqn{\delta}.  Specification of start location `x0`
#' and bounds `lb`, `ub` should still be specified in terms of \eqn{\delta}. The
#' `brcal` function will automatically convert from \eqn{\delta} to \eqn{\tau}.
#' In the returned list, `BR_params` will always report in terms of
#' \eqn{\delta}. However, the results returned in `nloptr` will reflect
#' whichever scale `nloptr()` optimized on.
#'
#' @inheritParams bayes_ms
#' @param t Minimum tolerable level of calibration in \[0,1\].
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
#' @param opts List with options to be passed to `nloptr`.
#' @param optim_options List with options to be passed to `optim`.
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
#' @importFrom stats sd
#' @importFrom stats optim
#' @importFrom nloptr nloptr
#'
#'
#'
#'
#' @references Birgin, E. G., and Martínez, J. M. (2008) Improving ultimate
#'   convergence of an augmented Lagrangian method, \emph{Optimization Methods
#'   and Software} vol. 23, no. 2, p. 177-195.
#'
#'   Conn, A. R., Gould, N. I. M., and Toint, P. L. (1991) A globally convergent
#'   augmented Lagrangian algorithm for optimization with general constraints
#'   and simple bounds, \emph{SIAM Journal of Numerical Analysis} vol. 28, no.
#'   2, p. 545-572.
#'
#'   Guthrie, A. P., and Franck, C. T. (2024) Boldness-Recalibration for Binary
#'   Event Predictions, \emph{The American Statistician} 1-17.
#'
#'   Johnson, S. G., The NLopt nonlinear-optimization package,
#'   <https://nlopt.readthedocs.io/en/latest/>.
#'
#'   Kraft, D. (1988) A software package for sequential quadratic programming",
#'   \emph{Technical Report} DFVLR-FB 88-28, Institut für Dynamik der
#'   Flugsysteme, Oberpfaffenhofen.
#'
#'   Kraft, D. (1994) Algorithm 733: TOMP-Fortran modules for optimal control
#'   calculations, \emph{ACM Transactions on Mathematical Software}, vol. 20,
#'   no. 3, pp. 262-281.
#'
#' @examples
#' 
#' # Simulate 50 predicted probabilities
#' x <- runif(50)
#' # Simulated 50 binary event outcomes using x
#' y <- rbinom(50, 1, x)  # By construction, x is well calibrated.
#'
#' # Perform 90% boldness-recalibration by setting t=0.9
#' # To suppress all output from nloptr() for each iteration use print_level=0
#' # (For reduced output at each iteration used print_level=1 or 2)
#' # To specify different starting values, use x0 and set start_at_MLEs=FALSE
#' brcal(x, y, t=0.9, x0=c(1,1), start_at_MLEs=FALSE, print_level=0)
#' 
#' # Adjust stopping criteria set max number of evaluations to 50 (maxeval) OR 
#' # stop after 0.5 second (maxtime)
#' # and set optimization bounds using lb and ub
#' brcal(x, y, maxeval = 50, maxtime = 0.5, lb=c(0.001, 0), ub=c(10, 10), print_level=0)
#' 
#' # Specify different options for nloptr & optim
#' brcal(x, y, opts=list(xtol_abs=0.01,
#'                      local_opts=list(algorithm="NLOPT_LD_MMA")), 
#'                      optim_options=list(method = "L-BFGS-B", lower = c(0, -1), 
#'                               upper = c(10, 25), control=list(factr=0.01)),
#'                               print_level=0)
#'                               
#' # Push probabilities away from bounds by 0.000001 and
#' # Stop outer optimization when parameters change by less than .001
#' x3 <- c(runif(25, 0, 0.0001), runif(25, .9999, 1))
#' y3 <- rbinom(50, 1, 0.5)
#' brcal(x3, y3, epsilon=0.000001, xtol_rel_outer = .01, print_level=0)
#' 
#' # See vignette for more examples
#' 
brcal <- function(x, y, t=0.95, Pmc=0.5, tau=FALSE, event=1,
                  start_at_MLEs=TRUE, x0=NULL,
                  lb=c(0.00001, -Inf), ub=c(Inf, Inf),
                  maxeval=500, maxtime=NULL,
                  xtol_rel_inner=1.0e-6,
                  xtol_rel_outer=1.0e-6,
                  print_level=3,
                  epsilon=.Machine$double.eps,
                  opts=NULL,
                  optim_options=NULL){
  
  
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
  t <- check_value01(t, name="t")
  
  # check Pmc is valid prior model prob
  Pmc <- check_value01(Pmc, name="Pmc")
  
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
  
  # check epsilon
  epsilon <- check_value01(epsilon, name="epsilon")
  
  ###################
  #  Function Code  #
  ###################
  
  if(start_at_MLEs){
    bt <- bayes_ms_internal(x,y, epsilon=epsilon)
    x0 <- bt$MLEs
  }
  
  if(tau){
    x0[1] <- log(x0[1])
    lb[1] <- log(lb[1])
    ub[1] <- log(ub[1])
  }
  
  if(missing(opts)){
    
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
                                        eval_grad_f = obj_grad_f,
                                        eval_jac_g_ineq = constr_grad_g,
                                        xtol_rel = xtol_rel_inner)),
                          probs = x,
                          outs = y,
                          t = t,
                          tau = tau,
                          Pmc = Pmc, 
                          epsilon=epsilon,
                          optim_options=optim_options)
  } else {
    
    
    
    # Check if algorithm, maxeval, maxtime, xtol_rel, or print_level are specified
    # use defaults if not
    if(!("algorithm" %in% names(opts))){
      opts$algorithm <- "NLOPT_LD_AUGLAG"
    }
    if(!("maxeval" %in% names(opts))){
      opts$maxeval <- maxeval
    }
    if(!("maxtime" %in% names(opts))){
      opts$maxtime <- maxtime
    }
    if(!("xtol_rel" %in% names(opts))){
      opts$xtol_rel <- xtol_rel_outer
    }
    if(!("print_level" %in% names(opts))){
      opts$print_level <- print_level
    }
    
    # Setting local opts
    if("local_opts" %in% names(opts)){
      
      # Overwrite any changes to the objective, constraint, or their jacobians
      opts$local_opts$eval_f <- obj_f
      opts$local_opts$eval_grad_f <- obj_grad_f
      opts$local_opts$eval_g_ineq <- constr_g
      opts$local_opts$eval_jac_g_ineq <- constr_grad_g
      
      # append default args if not specified by user
      if(!("algorithm" %in% names(opts$local_opts))){
        opts$local_opts$algorithm <- "NLOPT_LD_SLSQP"
      }
      if(!("xtol_rel" %in% names(opts$local_opts))){
        opts$local_opts$xtol_rel <- xtol_rel_inner
      }
      
    } else { # use default local opts
      opts$local_opts <- list(
        algorithm = "NLOPT_LD_SLSQP",
        eval_grad_f = obj_grad_f,
        eval_jac_g_ineq = constr_grad_g,
        xtol_rel = xtol_rel_inner)
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
                          Pmc = Pmc, 
                          epsilon=epsilon,
                          optim_options=optim_options)
    
    res$call$opts <- opts
  }
  
  
  
  # extract BR params
  br_params <- c(res$solution[1], res$solution[2])
  
  # # Make Call more useful when printed
  res$call$x0 <- round(res$x0, 6)
  res$call$lb <- res$lb
  res$call$ub <- res$ub
  
  nms <- names(res$call)[!(names(res$call) %in% c("", "eval_f",  "eval_grad_f", 
                                                  "eval_g_ineq", "eval_jac_g_ineq", 
                                                  "x0", "lb", "ub", "Pmc", "t", 
                                                  "tau", "epsilon"))]
  # print(nms)
  for(i in 1:length(nms)) {
    
    current <- nms[i]
    
    if(current == "opts"){
      if(is.null(opts)){
        
        nms2 <- names(res$call$opts)[!(names(res$call$opts) %in% c("", "eval_f",  "eval_grad_f", "eval_g_ineq", "eval_jac_g_ineq"))]
        
        for(j in 1:length(nms2)){
          
          current <- nms2[j]
          if(current == "local_opts"){
            nms3 <- names(res$call$opts$local_opts)[!(names(res$call$opts$local_opts) %in% c("", "eval_f",  "eval_grad_f", "eval_g_ineq", "eval_jac_g_ineq"))]
            
            for(k in 1:length(nms3)){
              current <- nms3[k]
              res$call$opts$local_opts[[current]] <- res$local_options[[current]]
            }
            
          }else{
            res$call$opts[[current]] <- res$options[[current]]
          }
        }
      } else{
        res$call$opts$local_opts$eval_f <- as.name("obj_f")
        res$call$opts$local_opts$eval_grad_f <- as.name("obj_grad_f")
        res$call$opts$local_opts$eval_g_ineq <- as.name("constr_g")
        res$call$opts$local_opts$eval_jac_g_ineq <- as.name("constr_grad_g")
      }
    } else{
      res$call[[current]] <- res[[current]]
    }
  }
  
  res$call$Pmc <- Pmc
  res$call$t <- t
  res$call$tau <- tau
  res$call$epsilon <- epsilon
  
  # Convert to delta scale if optimized on tau
  if(tau){
    br_params[1] <- exp(br_params[1])
  }
  
  # set up return list 
  l <- list(nloptr = res,
            Pmc = Pmc,
            t=t,
            BR_params = br_params,
            sb = -res$objective,
            probs = LLO_internal(x=x, br_params[1], br_params[2]))
  
  return(l)
}


######################################################
#  Internal Functions                                #
######################################################

obj_f <- function(x, probs, outs, t, tau, Pmc, epsilon, optim_options){
  if(tau)(
    x[1] <- exp(x[1])
  )
  # this assumes x is vector?
  probs_new <- LLO_internal(x=probs, x[1], x[2])
  return(-stats::sd(probs_new))
}

obj_grad_f <- function(x, probs, outs, t, tau, Pmc, epsilon, optim_options){
  if(tau)(
    x[1] <- exp(x[1])
  )
  n <- length(probs)
  probs_new <- LLO_internal(x=probs, x[1], x[2]) # x'
  sdp <- stats::sd(probs_new)
  meanp <- mean(probs_new)    # xbar'
  
  out <- (probs_new - meanp) 
  in1 <- (probs_new * (1 - probs_new))
  in2 <- log(probs / (1-probs))
  denom <- (n-1) * sdp
  
  d_grad <- -sum(out * (in1 - mean(in1))) / (denom *x[1])
  g_grad <- -sum(out * (in1*in2 - mean(in1*in2))) / denom
  
  grad_obj <- c(d_grad, g_grad)
  return(grad_obj)
}

constr_g <- function(x, probs, outs, t, tau, Pmc, epsilon, optim_options){
  if(tau)(
    x[1] <- exp(x[1])
  )
  probs_new <- LLO_internal(x=probs, x[1], x[2])
  bt <- do.call(bayes_ms_internal, c(list(x=probs_new, y=outs, Pmc = Pmc, epsilon=epsilon), optim_options))
  c1 <- bt$posterior_model_prob * -1 + t
  return(c1)
}

constr_grad_g <- function(x, probs, outs, t, tau, Pmc, epsilon, optim_options){
  if(tau)(
    x[1] <- exp(x[1])
  )
  n <- length(probs)
  probs_new <- LLO_internal(x=probs, x[1], x[2])
  bt <- do.call(bayes_ms_internal, c(list(x=probs_new, y=outs, Pmc = Pmc, epsilon=epsilon), optim_options))
  pmp <- bt$posterior_model_prob
  dhat <- bt$MLEs[1]
  ghat <- bt$MLEs[2]
  probs_MLE <- LLO_internal(x=probs_new, dhat, ghat)
  
  pmps <- pmp * (1-pmp)
  logp <- log(probs / (1-probs))
  
  d_grad <- pmps * (1/x[1]) * sum((ghat-1)*outs - ghat*probs_MLE + probs_new)
  g_grad <- pmps * sum(logp * ((ghat-1)*outs - ghat*probs_MLE + probs_new))
  
  grad_obj <- c(d_grad, g_grad)
  return(grad_obj)
}

