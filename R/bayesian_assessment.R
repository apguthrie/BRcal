######################################################
#  External Functions                                #
######################################################


# Bayesian Calibration assessment function
# IF ADDING OPT FOR DIFF NULL PARAMS, NEED TO RETURN THOSE AS WELL
#' Title
#'
#' @inheritParams llo_lrt
#' @param Pmc
#'
#' @return
#' @export
#'
#' @examples
bayes_ms <- function(x, y, Pmc = 0.5, optim_details = TRUE, event=1, ...){

  ##################
  #  Input Checks  #
  ##################

  # check x is vector, values in [0,1]
  x <- check_input_probs(x, name="x")

  # check y is vector, values are 0s or 1s - Relax this?
  y <- check_input_outcomes(y, name="y", event=event)

  # check Pmc is valid prior model prob
  Pmc <- check_input_probs(Pmc, name="Pmc")

  # check optim_details is logical
  if(!is.logical(optim_details) & !(optim_details %in% c(0,1))) stop("argument optim_details must be logical")

  # check x and y are the same length
  if(length(x) != length(y)) stop("x and y length differ")

  ###################
  #  Function Code  #
  ###################

  n <- length(x)
  params_null <- c(1,1)

  # BIC under null (well calibrated model Mc)
  # BIC_Mc <- BIC_llo(x = x, y = y, k = 0, params = params_null)
  BIC_Mc <- -2*llo_lik(params = params_null, x = x, y = y, log = TRUE)

  # Maximize likelihood
  optimlik <- llo_optim(x,y, tau=TRUE, ...)
  max_lik <- -optimlik$value
  MLEs <- optimlik$par

  # BIC under alternative (uncalibrated model Mu)
  BIC_Mu <- 2 * log(n) - (2 * max_lik)
  # temp <- BIC_llo(x = x, y = y, k = 2, params = NA, ...)
  # BIC_Mu <- temp$BIC

  # Bayes factors
  ## Likelihood of h1/likelihood of h0
  # BF_uc <- 1/bayes_factor(BIC1 = BIC1, BIC2 = BIC2)
  BF_uc <- exp(-(1/2) * (BIC_Mu - BIC_Mc))


  # Posterior Model Probabilities
  ## P(cal|data) = P(H0|data) = P(Mc|data)
  # post <- post_mod_prob(BF = BF_uc, Pmc = Pmc)
  Pmu <- 1 - Pmc
  post <- 1/(1+(BF_uc*(Pmu/Pmc)))

  # Return Value
  if(optim_details){
    results <- list(Pmc = Pmc,
                    BIC_Mc = BIC_Mc,
                    BIC_Mu = BIC_Mu,
                    BF = BF_uc,
                    posterior_model_prob = post,
                    MLEs = MLEs,
                    optim_details = optimlik)
  } else {
    results <- list(Pmc = Pmc,
                    BIC_Mc = BIC_Mc,
                    BIC_Mu = BIC_Mu,
                    BF = BF_uc,
                    posterior_model_prob = post,
                    MLEs = MLEs)
  }


  return(results)
}

######################################################
#  Internal Functions                                #
######################################################

# BIC for this likelihood
# BIC_llo <- function(x, y, k, params = NA, ...){
#   n <- length(x)
#   if(k == 0){
#     #suppressWarnings(if(is.na(params)) stop("must specify null params when k = 0"))
#     suppressWarnings(if(anyNA(params)) stop("must specify null params when k = 0"))
#     result <- -2*llo_lik(params = params, x = x, y = y, log = TRUE)
#   } else if(anyNA(params)){
#     optBayes <- llo_optim(x,y, tau=TRUE, ...)
#     max_lik <- -optBayes$value
#     MLEs <- optBayes$par
#     result <- list(BIC = k * log(n) - (2 * max_lik),
#                    est_params = MLEs)
#   } else {
#     result <- list(BIC = k * log(n) - (2 * llo_lik(params = params, x = x, y = y, log = TRUE)),
#                    params = params)
#   }
#   return(result)
# }


# # Bayes factor - only approx BIC version for now
# bayes_factor <- function(BIC1, BIC2,){
#   return(exp(-(1/2) * (BIC1 - BIC2)))
# }

# # Posterior model probability
# post_mod_prob <- function(BF, Pmc){
#   Pmu <- 1 - Pmc
#   return(1/(1+(BF*(Pmu/Pmc))))
# }


