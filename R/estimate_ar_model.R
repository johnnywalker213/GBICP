#' Estimate Autoregressive Model
#'
#' This function estimates an autoregressive (AR) model using four different model selection criteria: AIC, BIC, generalized BIC (GBIC), and generalized BIC with prior information. The method for GBIC with prior information is based on the reference: https://doi.org/10.1111/rssb.12023.
#'
#' @param time_series A numeric vector representing the time series data.
#' @param p An integer specifying the maximum lag order of the AR model.
#' @param criterion A string specifying the criterion for model selection. The valid options are "aic", "bic","gaic", "gbic", and "gbicp".
#'
#' @return The value of the specified criterion for the estimated AR model.
#'
#' @examples
#' ts_data <- rnorm(100) # Example time series data
#' ar_model_aic <- estimate_ar_model(ts_data, 1, "aic") # Estimate using AIC
#' ar_model_bic <- estimate_ar_model(ts_data, 1, "bic") # Estimate using BIC
#' ar_model_gbic <- estimate_ar_model(ts_data, 1, "gbic") # Estimate using GBIC
#' ar_model_gbicp <- estimate_ar_model(ts_data, 1, "gbicp") # Estimate using GBICP
#'
#' @references
#' Zhao, Z., & Shao, X. (2013). Model selection for autoregressive processes with fractional time trends. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 75(2), 323-341. https://doi.org/10.1111/rssb.12023
#'
#' @export
#'

estimate_ar_model <- function(time_series, p, criterion) {
  lagged_matrix <- create_lagged_matrix(time_series, p)
  #print(dim(lagged_matrix))
  n = length(time_series)
  Y <- time_series[-c(1:p)]
  X <- lagged_matrix

  #ar_model <- lm(Y ~ X-1)
  ar_model = forecast::Arima(time_series, order=c(p,0,0))
  n <- length(time_series)

  #beta = ar_model$coefficients
  beta = ar_model$coef[1:p]
  if (length(beta) == 1) {
    theta <- X * beta
  } else {
    theta <- X %*% beta
  }

  #residual = ar_model$residuals
  residual = time_series - ar_model$fitted
  residual = residual[-c(1:p)]

  mu_theta = 1/(theta)
  sigma_theta = diag(as.vector(-1/(theta)^2))
  An = var(residual)*t(X)%*%X
  Bn = t(X)%*%(diag(as.vector((residual)*(residual))))%*%X
  Hn = solve(An)%*%Bn

  if (criterion == "aic") {
    cr <- AIC(ar_model, k = 2)
    #return(list(order = p, model = ar_model, criterion = aic))
  } else if (criterion == "bic") {
    cr <- BIC(ar_model)
    #return(list(order = p, model = ar_model, criterion = bic))
  } else if (criterion == "gbicp") {
    cr = -2*logLik(ar_model)[1]+log(n)*p+sum(diag(Hn))-log(det(Hn))
  } else if (criterion == "gaic"){
    cr = -2*logLik(ar_model)[1]+2*sum(diag(Hn))
  } else if (criterion == "gbic")
    cr = -2*logLik(ar_model)[1]+log(n)*p-log(det(Hn))
  return(cr)
}
