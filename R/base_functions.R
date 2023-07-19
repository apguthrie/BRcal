# Linear Log Odds Recalibration Function
LLO <- function(p, delta, gamma){
  # check p are probabilities in [0,1]

  # check delta > 0 & numeric & size 1

  # check gamma in Reals & numeric & size 1

  return((delta * (p^gamma)) / ((delta * (p^gamma)) + ((1-p)^gamma)))
}

# Converts probs to logit scale
logit <- function(p){

  # check p are probabilities in [0,1]

  # better way to handle the rounding here? - check literature
  p <- ifelse(p < (10^(-300)), (10^(-300)), p)
  p <- ifelse(p > 0.9999999999999999, 0.9999999999999999, p)

  return(log(p/(1-p)))
}

# Likelihood
llo_lik <- function(params, x, y, log = FALSE, neg = FALSE){

  # check params are of right length, right values

  # check x's are between 0,1

  # check y's are 0s or 1s

  # check log & neg are logical

  # rounding off x's that are too close to zero or one
  x <- ifelse(x < (10^(-300)), (10^(-300)), x)
  x <- ifelse(x > 0.9999999999999999, 0.9999999999999999, x)

  llo <- LLO(p = x, delta = params[1], gamma = params[2])

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
