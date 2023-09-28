# add option to pass args to nloptr
# also return br probs (add flag to toggle this)

brcal <- function(x,y,t_level, x0=c(0.5,0.5), print_level=3, maxeval=300,
                  xtol_rel_outer=1.0e-6, start_at_MLEs=TRUE,
                  xtol_rel_inner=1.0e-6,
                  # algorithm="NLOPT_LD_SLSQP",
                  check_derivatives=FALSE){

  if(start_at_MLEs){
    bt <- bayes_testing(x,y)
    x0 <- bt$est_params
  }

  eval_f <- function(x, probs, outs, t_level){
    # this assumes x is vector?
    probs_new <- LLO(probs, x[1], x[2])
    return(-stats::sd(probs_new))
  }

  eval_grad_f <- function(x, probs, outs, t_level){


    n <- length(probs)
    probs_new <- LLO(probs, x[1], x[2])
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

  eval_g <- function(x, probs, outs, t_level){

    probs_new <- LLO(probs, x[1], x[2])
    c1 <- bayes_testing(probs_new, outs)$posterior_model_prob * -1 + t_level
    return(c1)
  }

  eval_grad_g <- function(x, probs, outs, t_level){


    n <- length(probs)
    probs_new <- LLO(probs, x[1], x[2])
    bt <- bayes_testing(probs_new, outs)
    pmp <- bt$posterior_model_prob
    pmp2 <- pmp^2
    dhat <- bt$est_params[1]
    ghat <- bt$est_params[2]

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

  # res <- nloptr::nloptr(x0 = x0,
  #                       eval_f = eval_f,
  #                       eval_grad_f = eval_grad_f,
  #                       lb = c(0.0001, -15),
  #                       ub = c(Inf, 1074),
  #                       eval_g_ineq = eval_g,
  #                       eval_jac_g_ineq  = eval_grad_g,
  #                       opts = list("algorithm"=algorithm,
  #                                   "xtol_rel"=xtol_rel,
  #                                   "print_level"=print_level,
  #                                   "maxeval"=maxeval, "check_derivatives"=check_derivatives),
  #                       probs = x,
  #                       outs = y,
  #                       t_level = t_level)


  res <- nloptr::nloptr(x0 = x0,
                        eval_f = eval_f,
                        eval_grad_f = eval_grad_f,
                        eval_g_ineq = eval_g,
                        eval_jac_g_ineq  = eval_grad_g,
                        opts = list(algorithm = "NLOPT_LD_AUGLAG",
                                    maxeval = maxeval,
                                    xtol_rel = xtol_rel_outer,
                                    print_level = print_level,
                                    #xtol_abs1 = 1e-4,
                                    local_opts = list(
                                      algorithm = "NLOPT_LD_SLSQP",
                                      eval_grad_f = eval_grad_f,
                                      eval_jac_g_ineq  = eval_grad_g,
                                      #xtol_abs1 = 1e-4,
                                      xtol_rel = xtol_rel_inner
                                    )),
                        probs = x,
                        outs = y,
                        t_level = t_level)

  return(res)
}


