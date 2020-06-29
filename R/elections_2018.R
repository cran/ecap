#' FiveThirtyEight probabilities from the 2018 election cycle.
#'
#' A dataset containing FiveThirtyEight's probability estimates for the 2018 US elections. The races covered seats
#' in the House, Senate, and for Governor. They also used three different methods to compute these estimates. Classic, Deluxe,
#' and Lite. All of these have been included. The probabilities can be viewed from the perspective of a Democrat winning
#' or from the perspective Republican winning any given race. The final results of each race are also included.
#'
#' @format A data frame with 1518 rows and 1 variables:
#' \describe{
#' \item{cycle}{the year the race took place}
#' \item{branch}{what branch of the government the race was for}
#' \item{race}{the specific position the election was for}
#' \item{forecastdate}{date of forecast}
#' \item{version}{type of prediction method FiveThirtyEight used}
#' \item{Democrat_WinProbability}{probability of a democrat winning the race}
#' \item{Republican_WinProbability}{probability of a republican winning the race}
#' \item{category}{category}
#' \item{Democrat_Won}{binary variable indicating if a democrat won or not}
#' \item{Republican_Won}{binary variable indicating if a republican won or not}
#' \item{uncalled}{binary variable indicating if the race was uncalled or not}
#' }
#' @source \url{https://github.com/fivethirtyeight/data/tree/master/forecast-review}
"elections_2018"
