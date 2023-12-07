#' Select Best Autoregressive Model
#'
#' This function selects the best autoregressive (AR) model based on a specified criterion. It iteratively estimates AR models for different lag orders and chooses the model with the best criterion value.
#'
#' @param time_series A numeric vector representing the time series data.
#' @param max_p An integer specifying the maximum lag order to consider for the AR model.
#' @param criterion A string specifying the criterion for model selection. Valid options are "aic", "bic", "gaic", "gbic", and "gbicp".
#'
#' @return A numeric vector containing two elements: the lag order of the best AR model and the corresponding criterion value.
#'
#' @examples
#' ts_data <- rnorm(100) # Example time series data
#' best_ar_model <- select_best_ar_model(ts_data, 5, "aic") # Selecting the best AR model using AIC
#'
#' @export
#'

select_best_ar_model <- function(time_series, max_p, criterion) {
  best_model <- NULL
  best_criterion_value <- Inf

  for (p in 1:max_p) {
    current_model <- estimate_ar_model(time_series, p, criterion)
    if (current_model < best_criterion_value) {
      best_model <- p
      best_criterion_value <- current_model
    }
  }

  return(c(best_model,best_criterion_value))
}
