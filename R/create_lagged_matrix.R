#' Create Lagged Matrix
#'
#' This function generates a matrix where each column is a lagged version of the input time series.
#' The number of lags (p) is specified by the user. NA values are used to fill in where data
#' is not available due to lagging, but the output will automatically remove rows are not complete.
#'
#' @param time_series A numeric vector representing the time series data.
#' @param p An integer specifying the number of lags to create.
#'
#' @return A matrix with p columns, each representing a lagged version of the input time series.
#'         The matrix will have the same number of rows as the input time series, but the
#'         first p rows will contain NA values to account for the lag.
#'
#' @examples
#' ts <- 1:10
#' create_lagged_matrix(ts, 2)
#'
#' @export

create_lagged_matrix <- function(time_series, p) {
  # Get the length of the time series
  n <- length(time_series)

  # Create a matrix where each column is a lagged version of the original series
  lagged_matrix <- sapply(1:p, function(i) c(rep(NA, i), time_series[1:(n-i)]))

  return(lagged_matrix[-c(1:p),])
}
