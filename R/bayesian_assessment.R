######################################################
#  External Functions                                #
######################################################


#' Bayesian Model Selection-Based Calibration Assessment
#'
#' Perform Bayesian model selection-based calibration assessment as specified in
#' Guthrie and Franck (2024).
#'
#' Compares a well calibrated model, \eqn{M_c} where \eqn{\delta = \gamma = 1}
#' to an uncalibrated model, \eqn{M_u} where \eqn{\delta>0, \gamma \in
#' \mathbb{R}}.
#'
#'
#' @inheritParams llo_lrt
#' @param Pmc The prior model probability for the calibrated model \eqn{M_c}.
#'
#' @return A list with the following attributes:
#'   \item{\code{Pmc}}{The prior
#'   model probability for the calibrated model \eqn{M_c}.}
#'   \item{\code{BIC_Mc}}{The Bayesian Information Criteria (BIC) for the
#'   calibrated model \eqn{M_c}.}
#'   \item{\code{BIC_Mu}}{The Bayesian Information Criteria
#'   (BIC) for the uncalibrated model \eqn{M_u}.}
#'   \item{\code{BF}}{The Bayes Factor of uncalibrated model over calibrated
#'   model.}
#'   \item{\code{posterior_model_prob}}{The posterior model probability of the
#'   calibrated model \eqn{M_c} given the observed outcomes `y`, i.e. \eqn{P(M_c|y)}.}
#'   \item{\code{MLEs}}{Maximum likelihood estimates for \eqn{\delta} and
#'   \eqn{\gamma}.}
#'   \item{\code{optim_details}}{If `optim_details = TRUE`, the list returned by
#'   `optim()` when minimizing the negative log likelihood, includes convergence
#'   information, number of iterations, and achieved negative log likelihood
#'   value and MLEs.}
#' @export
#'
#' @importFrom stats optim
#'
#' @references Guthrie, A. P., and Franck, C. T. (2024) Boldness-Recalibration
#'   for Binary Event Predictions. \emph{arxiv}.
#'
#' @examples
bayes_ms <- function(x, y, Pmc = 0.5, event=1, optim_details = TRUE,  ...){

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
  if(!is.logical(optim_details) & !(optim_details %in% c(0,1))){
    stop("argument optim_details must be logical")
  }

  # check x and y are the same length
  if(length(x) != length(y)) stop("x and y length differ")

  ###################
  #  Function Code  #
  ###################

  results <- bayes_ms_internal(x=x, y=y, Pmc=Pmc, optim_details=optim_details, ...)
  return(results)
}

######################################################
#  Internal Functions                                #
######################################################

bayes_ms_internal <- function(x, y, Pmc = 0.5, optim_details = TRUE,  ...){
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

  # Bayes factors
  ## Likelihood of h1/likelihood of h0
  BF_uc <- exp(-(1/2) * (BIC_Mu - BIC_Mc))


  # Posterior Model Probabilities
  ## P(cal|data) = P(H0|data) = P(Mc|data)
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



