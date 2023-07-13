# add option to pass args to nloptr
# also return br probs (add flag to toggle this)

brcal <- function(x,y,t_level, x0=c(0.5,0.5), print_level=3){

  eval_f <- function(x, probs, outs, t_level){
    # this assumes x is vector?
    probs_new <- LLO(probs, x[1], x[2])
    return(-stats::sd(probs_new) )
  }

  eval_g <- function(x, probs, outs, t_level){

    probs_new <- LLO(probs, x[1], x[2])
    c1 <- bayes_testing(probs_new, outs)$posterior_model_prob * -1 + t_level
    return(c1)
  }

  res <- nloptr::nloptr(x0 = x0,
                  eval_f = eval_f,
                  lb = c(0.0001, -Inf),
                  ub = c(Inf, Inf),
                  eval_g_ineq = eval_g,
                  opts = list("algorithm"="NLOPT_LN_COBYLA",
                              "xtol_rel"=1.0e-8,
                              "print_level"=print_level),
                  probs = x,
                  outs = y,
                  t_level = t_level)
  return(res)
}


brcal_dev <- function(x,y,t_level, x0=c(0.5,0.5), print_level=3){

  eval_f <- function(x, probs, outs, t_level){
    # this assumes x is vector?
    probs_new <- LLO(probs, x[1], x[2])
    return(-stats::sd(probs_new) )
  }

  eval_g <- function(x, probs, outs, t_level){

    probs_new <- LLO(probs, x[1], x[2])
    c1 <- bayes_testing_dev(probs_new, outs)$posterior_model_prob * -1 + t_level
    return(c1)
  }

  res <- nloptr::nloptr(x0 = x0,
                        eval_f = eval_f,
                        lb = c(0.0001, -Inf),
                        ub = c(Inf, Inf),
                        eval_g_ineq = eval_g,
                        opts = list("algorithm"="NLOPT_LN_COBYLA",
                                    "xtol_rel"=1.0e-8,
                                    "print_level"=print_level),
                        probs = x,
                        outs = y,
                        t_level = t_level)
  return(res)
}
