#' Hockey Home Team Win Predictions data
#'
#' Home team win probability predictions and outcomes pertaining to the 2020-21
#' National Hockey League (NHL) Season. Probability predictions x were obtained
#' from FiveThirtyEight via downloadable spreadsheet on their website (see below
#' for link).  The win/loss game results were obtained by web-scraping from
#' NHL.com using the NHL API.
#'
#' @format ## `hockey`
#' A data frame with 868 rows and 4 columns:
#' \describe{
#'   \item{y}{game result, 1 = home team win, 0 = home team loss}
#'   \item{x}{predicted probabilities of a home team win from FiveThirtyEight}
#'   \item{rand}{uniformly random generated predicted probability of a home team from range \[0.26, 0.78\] }
#'   \item{winner}{game result (string), "home" = home team win, "away" = home team loss}
#' }
#' @source <https://data.fivethirtyeight.com/>
"hockey"

#' Foreclosure Monitoring Predictions data
#'
#' Foreclosure monitoring probability predictions and the true foreclosure 
#' status pertaining of 5,000 housing transactions in 2010 
#' from Wayne County, Michigan. These data were a randomly selected subset from
#' data from presented in Keefe et al. (2017). 
#'
#' @format ## `foreclosure`
#' A data frame with 5,000 rows and 3 columns:
#' \describe{
#'   \item{y}{sale type, 1 = foreclosure, 0 = regular sale}
#'   \item{x}{predicted probabilities of foreclosure}
#'   \item{year}{year of observed foreclosure or regular sale}
#' }
#' @source  Keefe, M.J., Franck, C.T., Woodall, W.H. (2017): Monitoring foreclosure 
#' rates with a spatially risk-adjusted bernoulli cusum chart for concurrent 
#' observations. \emph{Journal of Applied Statistics} 44(2), 325–341  
#' <https://doi.org/10.1080/02664763.2016.1169257>

"foreclosure"