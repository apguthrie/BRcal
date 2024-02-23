######################################################
#  External Functions                                #
######################################################


# Bayesian Calibration assessment function
bayes_ms <- bayes_testing <- function(x, y, Pmc = 0.5, params_null = c(1,1), optim_details = FALSE,
                                      event=1, ...){
  # check y only has two values
  y <- ifelse(y == event, 1, 0)

  # BIC under null (well calibrated model Mc)
  BIC_Mc <- BIC_llo(x = x, y = y, k = 0, params = params_null)

  # BIC under alternative (uncalibrated model Mu)
  temp <- BIC_llo(x = x, y = y, k = 2, params = NA, ...)
  BIC_Mu <- temp$BIC


  # Bayes factors
  ## Likelihood of h1/likelihood of h0
  BF_uc <- 1/bayes_factor(BIC1 = BIC1, BIC2 = BIC2)

  # Posterior Model Probabilities
  ## P(cal|data) = P(H0|data) = P(Mc|data)
  post <- post_mod_prob(BF = BF_uc, Pmc = Pmc)

  if(anyNA(params)){
    est_params <- temp$est_params
    results <- list(BIC_Mc = BIC_Mc,
                    BIC_Mu = BIC_Mu,
                    BF = BF_uc,
                    posterior_model_prob = post,
                    MLEs = est_params)
  } else{
    results <- list(BIC_Mc = BIC_Mc,
                    BIC_Mu = BIC_Mu,
                    BF = BF_uc,
                    posterior_model_prob = post,
                    MLEs = params)
  }

  return(results)
}

######################################################
#  Internal Functions                                #
######################################################

# BIC for this likelihood
BIC_llo <- function(x, y, k, params = NA, ...){
  n <- length(x)
  if(k == 0){
    #suppressWarnings(if(is.na(params)) stop("must specify null params when k = 0"))
    suppressWarnings(if(anyNA(params)) stop("must specify null params when k = 0"))
    result <- -2*llo_lik(params = params, x = x, y = y, log = TRUE)
  } else if(anyNA(params)){
    optBayes <- llo_optim(x,y, tau=TRUE, ...)
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
bayes_factor <- function(BIC1, BIC2,){
  return(exp(-(1/2) * (BIC1 - BIC2)))
}

# Posterior model probability
post_mod_prob <- function(BF, Pmc){
  Pmu <- 1 - Pmc
  return(1/(1+(BF*(Pmu/Pmc))))
}


