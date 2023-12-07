#' Estimate Phi Function
#'
#' This function estimates the phi parameter of a given time series. It uses a
#' local polynomial regression (loess) to estimate the variance of the residuals.
#'
#' @param time_series A numeric vector representing the time series data.
#' @param p An integer specifying the number of lags to be used.
#' @param b A numeric value specifying the bandwidth for the local polynomial regression.
#' @return A list containing three elements:
#'   - `phi_hat`: the estimated phi parameter,
#'   - `sigma_hat_2`: the estimated variance of the residuals,
#'   - `sigma_t_2`: the estimated variances at each time point.
#' @examples
#' ts_data <- rnorm(100)
#' estimate_phi(ts_data, p = 5, b = 0.5)
#' @export
#' @importFrom stats var
#' @importFrom KernSmooth locpoly
estimate_phi <- function(time_series, p, b) {

  lagged_X = create_lagged_matrix(time_series, p)
  phi_hat = solve(t(lagged_X) %*% lagged_X) %*% t(lagged_X) %*% time_series[-c(1:p)]

  err_est = lagged_X %*% phi_hat

  res = time_series[-c(1:p)] - err_est
  res_sq = res^2

  est = KernSmooth::locpoly(x = seq_along(res_sq)/length(res_sq),
                y= res_sq,
                gridsize = length(res_sq),
                kernel = "normal",
                bandwidth = b)

  W = diag(est$y)
  phi_hat_new = solve(t(lagged_X) %*% W %*% lagged_X) %*% t(lagged_X) %*% W %*% time_series[-c(1:p)]
  err_est_new = lagged_X %*% phi_hat_new
  res_new = time_series[-c(1:p)] - err_est_new
  sigma_hat_2 = var(res_new)

  return(list('phi_hat' = phi_hat_new, 'sigma_hat_2' = sigma_hat_2, 'sigma_t_2' = est$y))
}

