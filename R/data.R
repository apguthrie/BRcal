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
