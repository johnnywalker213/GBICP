#' Time Series Simulation Function
#'
#' Generates a simulated time series based on specified autoregressive (AR), moving average (MA) coefficients, and other parameters.
#'
#' @param ar_coef Vector of autoregressive coefficients. Default is \code{NULL}.
#' @param ma_coef Vector of moving average coefficients. Default is \code{NULL}.
#' @param n Length of the time series to generate.
#' @param sigma2 Variance of the noise in the time series. Default is 1.
#' @return A numeric vector representing the generated time series.
#' @examples
#' ts1 <- generate_ts(n = 100)
#' ts2 <- generate_ts(ar_coef = c(0.5, -0.2), n = 100)
#' ts3 <- generate_ts(ma_coef = c(0.5, 0.3), n = 100)
#' ts4 <- generate_ts(ar_coef = c(0.5, -0.2), ma_coef = c(0.5, 0.3), n = 100)
#' @export



generate_ts <- function(ar_coef = NULL, ma_coef = NULL, n, sigma2 = 1) {
  ts <- numeric(n)
  # generate time values from 0 to 10
  t <- seq(1, n, length.out = n)

  # generate normally distributed noise
  noise <- rnorm(n, mean = 0, sd = sqrt(sigma2))

  # create bell shape with time-varying scale factor
  #scale_factor <- 1*exp(-((t-midpoint)/n )^2/0.05)+0.1
  #noise_bell <- noise * scale_factor
  noise_bell = noise

  # Generate random noise if both AR and MA coefficients are NULL
  if (is.null(ar_coef) & is.null(ma_coef)) {
    ts <- noise_bell

    # Generate MA process if AR coefficients are NULL
  } else if (is.null(ar_coef)) {
    e <- noise
    for (i in 1:n) {
      if (i <= length(ma_coef)) {
        ma_term <- 0
        for (j in 1:i) {
          ma_term <- ma_term + ma_coef[j]*e[i-j+1]
        }
        ts[i] <- ma_term+noise_bell[i]
      } else {
        ma_term <- 0
        for (j in 1:length(ma_coef)) {
          ma_term <- ma_term + ma_coef[j]*e[i-j+1]
        }
        ts[i] <- ma_term+noise_bell[i]
      }
    }

    # Generate AR process if MA coefficients are NULL
  } else if (is.null(ma_coef)) {
    ts[1:length(ar_coef)] <- rnorm(length(ar_coef))
    for (i in (length(ar_coef)+1):n) {
      ar_term <- 0
      for (j in 1:length(ar_coef)) {
        ar_term <- ar_term + ar_coef[j]*ts[i-j]
      }
      ts[i] <- ar_term + +noise_bell[i]
    }

    # Generate ARMA process if both AR and MA coefficients are given
  } else {
    e <- noise
    ts[1:length(ar_coef)] <- rnorm(length(ar_coef))
    for (i in (length(ar_coef)+1):n) {
      ar_term <- 0
      for (j in 1:length(ar_coef)) {
        ar_term <- ar_term + ar_coef[j]*ts[i-j]
      }
      if (i <= length(ma_coef)) {
        ma_term <- 0
        for (j in 1:i) {
          ma_term <- ma_term + ma_coef[j]*e[i-j+1]
        }
        ts[i] <- ar_term + ma_term
      } else {
        ma_term <- 0
        for (j in 1:length(ma_coef)) {
          ma_term <- ma_term + ma_coef[j]*e[i-j+1]
        }
        ts[i] <- ar_term + ma_term
      }
      ts[i] <- ts[i] +noise_bell[i]
    }
  }

  return(ts)
}
