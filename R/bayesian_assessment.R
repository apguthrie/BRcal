######################################################
#  External Functions                                #
######################################################

#' Bayesian Model Selection-Based Calibration Assessment
#'
#' Perform Bayesian model selection-based approach to determine if a set of
#' predicted probabilities `x` is well calibrated given the corresponding set of
#' binary event outcomes `y` as described in Guthrie and Franck (2024).
#'
#' This function compares a well calibrated model, \eqn{M_c} where \eqn{\delta =
#' \gamma = 1} to an uncalibrated model, \eqn{M_u} where \eqn{\delta>0, \gamma \in
#' \mathbb{R}}.
#'
#' The posterior model probability of \eqn{M_c} given the observed
#' outcomes `y` (returned as `posterior_model_prob`) is expressed as \deqn{P(M_c|\mathbf{y})
#' = \frac{P(\mathbf{y}|M_c) P(M_c)}{P(\mathbf{y}|M_c) P(M_c) + P(\mathbf{y}|M_{u}) P(M_{u})}}
#' where \eqn{P(\mathbf{y}|M_i)} is the integrated likelihoof of `y` given
#' \eqn{M_i} and \eqn{P(M_i)} is the prior probability of model i, \eqn{i \in
#' \{c,u\}}. By default, this function uses \eqn{P(M_c) = P(M_u) = 0.5}. To set a
#' different prior for \eqn{P(M_c)}, use `Pmc`, and \eqn{P(M_u)} will be set to
#' `1 - Pmc`.
#'
#' The Bayes factor (returned as `BF`) compares \eqn{M_u} to \eqn{M_c}.  This
#' value is approximated via the following large sample Bayesian Information
#' Criteria (BIC) approximation (see Kass & Raftery 1995, Kass & Wasserman 1995) \deqn{BF =
#' \frac{P(\mathbf{y}|M_{u})}{P(\mathbf{y}|M_c)} = \approx exp\left\{
#' -\frac{1}{2}(BIC_u - BIC_c) \right\}} where the BIC for the calibrated model
#' (returned as `BIC_mc`) is \deqn{BIC_c = - 2 \times log(\pi(\delta = 1, \gamma =1|\mathbf{x},\mathbf{y}))}
#' and the BIC for the uncalibrated model (returned as `BIC_mu`) is \deqn{BIC_u =
#' 2\times log(n) - 2\times log(\pi(\hat\delta_{MLE}, \hat\gamma_{MLE}|\mathbf{x},\mathbf{y})).}
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
#'   Kass, R. E., and Raftery, A. E. (1995) Bayes factors. \emph{Journal of the 
#'   American Statistical Association}
#'
#'   Kass, R. E., and Wassermann, L. (1995) A reference bayesian test for nested 
#'   hypotheses and its relationship to the schwarz criterion. \emph{Journal of 
#'   the American Statistical Association}
#'   
#' @examples
#' # Simulate 100 predicted probabilities
#' x <- runif(100)
#' # Simulated 100 binary event outcomes using x
#' y <- rbinom(100, 1, x)  # By construction, x is well calibrated.
#'
#' # Use bayesian model selection approach to check calibration of x given outcomes y
#' bayes_ms(x, y, optim_details=FALSE)
#'
#' # To specify different prior model probability of calibration, use Pmc
#' # Prior model prob of 0.7:
#' bayes_ms(x, y, Pmc=0.7)
#' # Prior model prob of 0.2
#' bayes_ms(x, y, Pmc=0.2)
#'
#' # Use optim_details = TRUE to see returned info from call to optim(),
#' # details useful for checking convergence
#' bayes_ms(x, y, optim_details=TRUE)  # no convergence problems in this example
#'
#' # Pass additional arguments to optim() via ... (see optim() for details)
#' # Specify different start values via par in optim() call, start at delta = 5, gamma = 5:
#' bayes_ms(x, y, optim_details=TRUE, par=c(5,5))
#' # Specify different optimization algorithm via method, L-BFGS-B instead of Nelder-Mead:
#' bayes_ms(x, y, optim_details=TRUE, method = "L-BFGS-B")  # same result
#'
#' # What if events are defined by text instead of 0 or 1?
#' y2 <- ifelse(y==0, "Loss", "Win")
#' bayes_ms(x, y2, event="Win", optim_details=FALSE)  # same result
#'
#' # What if we're interested in the probability of loss instead of win?
#' x2 <- 1 - x
#' bayes_ms(x2, y2, event="Loss", optim_details=FALSE)
#'
#' # Push probabilities away from bounds by 0.000001
#' x3 <- c(runif(50, 0, 0.0001), runif(50, .9999, 1))
#' y3 <- rbinom(100, 1, 0.5)
#' bayes_ms(x3, y3, epsilon=0.000001)
#' 
bayes_ms <- function(x, y, Pmc = 0.5, event = 1, optim_details = TRUE, 
                     epsilon = .Machine$double.eps, ...){

  ##################
  #  Input Checks  #
  ##################

  # check x is vector, values in [0,1]
  x <- check_input_probs(x, name="x")

  # check y is vector, values are 0s or 1s - Relax this?
  y <- check_input_outcomes(y, name="y", event=event)

  # check Pmc is valid prior model prob
  Pmc <- check_value01(Pmc, name="Pmc")

  # check optim_details is logical
  if(!is.logical(optim_details) & !(optim_details %in% c(0,1))){
    stop("argument optim_details must be logical")
  }

  # check x and y are the same length
  if(length(x) != length(y)) stop("x and y length differ")

  # check epsilon
  epsilon <- check_value01(epsilon, name="epsilon")
  
  ###################
  #  Function Code  #
  ###################

  results <- bayes_ms_internal(x=x, y=y, Pmc=Pmc, optim_details=optim_details, epsilon=epsilon, ...)
  return(results)
}

######################################################
#  Internal Functions                                #
######################################################

bayes_ms_internal <- function(x, y, Pmc = 0.5, optim_details = TRUE, epsilon=.Machine$double.eps,  ...){
  n <- length(x)
  params_null <- c(1,1)

  # BIC under null (well calibrated model Mc)
  # BIC_Mc <- BIC_llo(x = x, y = y, k = 0, params = params_null)
  BIC_Mc <- -2*llo_lik(params = params_null, x = x, y = y, log = TRUE, epsilon=epsilon)

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



