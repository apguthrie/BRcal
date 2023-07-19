
# BIC for this likelihood
BIC_llo <- function(x, y, k, params = NA, lower = c(0.001, -5), upper = c(10,30)){
  n <- length(x)
  if(k == 0){
    #suppressWarnings(if(is.na(params)) stop("must specify null params when k = 0"))
    suppressWarnings(if(anyNA(params)) stop("must specify null params when k = 0"))

    result <- -2*llo_lik(params = params, x = x, y = y, log = TRUE)
  } else if(anyNA(params)){
    # optBayes <- stats::optim(c(0.5, 0.5), llo_lik, x=x, y=y, log = TRUE, neg = TRUE, method = "L-BFGS-B",
    #                   lower = lower, upper = upper)
    optBayes <- stats::optim(c(0.5, 0.5), llo_optim_wrap, x=x, y=y, method = "Nelder-Mead",
                           neg = TRUE, log = TRUE)
    max_lik <- -optBayes$value
    MLEs <- optBayes$par
    result <- list(BIC = k * log(n) - (2 * max_lik),
                   est_params = MLEs)
  } else {
    result <- list(BIC = k * log(n) - (2 * llo_lik(params = params, x = x, y = y, log = TRUE)),
                   params = params)
  }
  return(result)
}

BIC_llo_dev <- function(x, y, k, params = NA, lower = c(0.001, -5), upper = c(10,30)){
  n <- length(x)
  if(k == 0){
    #suppressWarnings(if(is.na(params)) stop("must specify null params when k = 0"))
    suppressWarnings(if(anyNA(params)) stop("must specify null params when k = 0"))

    result <- -2*llo_lik(params = params, x = x, y = y, log = TRUE)
  } else if(anyNA(params)){
    optBayes <- stats::optim(c(0.5, 0.5), llo_optim_wrap, x=x, y=y, log = TRUE, neg = TRUE, method = "Nelder-Mead")
                             #lower = lower, upper = upper)
    max_lik <- -optBayes$value
    MLEs <- optBayes$par
    result <- list(BIC = k * log(n) - (2 * max_lik),
                   est_params = MLEs)
  } else {
    result <- list(BIC = k * log(n) - (2 * llo_lik(params = params, x = x, y = y, log = TRUE)),
                   params = params)
  }
  return(result)
}

# Bayes factor - only approx BIC version for now
bayes_factor <- function(BIC1, BIC2, approx = TRUE){
  return(exp(-(1/2) * (BIC1 - BIC2)))
}

# Posterior model probability
post_mod_prob <- function(BF){
  return(1/(1+BF))
}

# Bayesian Calibration assessment function
bayes_testing <- function(x, y, k = 2, params_null = c(1,1), params = NA, optim_details = FALSE, lower = c(0.001, -5), upper = c(10,30)){
  # BIC under null
  BIC1 <- BIC_llo(x = x, y = y, k = 0, params = params_null, lower = lower, upper = upper)

  # BIC under alternative
  temp <- BIC_llo(x = x, y = y, k = k, params = params, lower = lower, upper = upper)
  BIC2 <- temp$BIC


  # Bayes factors
  ## Likelihood of h0/likelihood of h1
  BF12 <- bayes_factor(BIC1 = BIC1, BIC2 = BIC2)
  ## Likelihood of h1/likelihood of h0
  BF21 <- 1/BF12

  # Posterior Model Probabilities
  ## P(cal|data) = P(H0|data)
  post1 <- post_mod_prob(BF = BF21)

  ## P(not cal|data) = P(H1|data)
  post2 <- 1-post1

  if(anyNA(params)){
    est_params <- temp$est_params
    results <- list(BIC_H0 = BIC1,
                    BIC_H1 = BIC2,
                    BF = BF21,
                    posterior_model_prob = post1,
                    est_params = est_params)
  } else{
    results <- list(BIC_H0 = BIC1,
                    BIC_H1 = BIC2,
                    BF = BF21,
                    posterior_model_prob = post1,
                    params = params)
  }

  return(results)
}

bayes_testing_dev <- function(x, y, k = 2, params_null = c(1,1), params = NA, optim_details = FALSE, lower = c(0.001, -5), upper = c(10,30)){
  # BIC under null
  BIC1 <- BIC_llo_dev(x = x, y = y, k = 0, params = params_null, lower = lower, upper = upper)

  # BIC under alternative
  temp <- BIC_llo_dev(x = x, y = y, k = k, params = params, lower = lower, upper = upper)
  BIC2 <- temp$BIC


  # Bayes factors
  ## Likelihood of h0/likelihood of h1
  BF12 <- bayes_factor(BIC1 = BIC1, BIC2 = BIC2)
  ## Likelihood of h1/likelihood of h0
  BF21 <- 1/BF12

  # Posterior Model Probabilities
  ## P(cal|data) = P(H0|data)
  post1 <- post_mod_prob(BF = BF21)

  ## P(not cal|data) = P(H1|data)
  post2 <- 1-post1

  if(anyNA(params)){
    est_params <- temp$est_params
    results <- list(BIC_H0 = BIC1,
                    BIC_H1 = BIC2,
                    BF = BF21,
                    posterior_model_prob = post1,
                    est_params = est_params)
  } else{
    results <- list(BIC_H0 = BIC1,
                    BIC_H1 = BIC2,
                    BF = BF21,
                    posterior_model_prob = post1,
                    params = params)
  }

  return(results)
}
