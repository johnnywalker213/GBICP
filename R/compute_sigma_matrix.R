#' Compute Sigma Matrix
#'
#' This function calculates the inverse covariance matrix for a given time series with a specified lag order (p).
#' It employs the function of 'estimate_phi' and the variance at each level to construct
#' the covariance matrix in a band limited manner.
#'
#' @param time_series A numeric vector representing the time series data.
#' @param p An integer indicating the lag order to be used in the calculation.
#'
#' @return An inverse covariance matrix (precision matrix) of the given time series.
#'
#' @examples
#' ts_data <- rnorm(100) # Example time series data
#' inv_cov_matrix <- compute_sigma_matrix(ts_data, 1) # Computing the inverse covariance matrix
#'
#' @export
#'
compute_sigma_matrix <- function(time_series, p) {
  n = length(time_series)
  # Initialize gamma matrix
  L_k <- diag(rep(1, n))

  # lambda inverse matrix diagonal
  lambda_inv = rep(0, n)

  # special case for (1,1)
  lambda_inv[1] = 1/var(time_series)

  if (p == 1) {
    lvl_p_result = estimate_phi(time_series, p=p, b=0.04)
    phi = lvl_p_result$phi_hat
    sigma_p_sq = lvl_p_result$sigma_hat_2[1]
    lambda_inv[(p+1):n] = 1/lvl_p_result$sigma_t_2
    #print(lvl_p_result$sigma_t_2)
    #print(lambda_inv)

    for (i in (p+1):n) {
      L_k[i, (i-p):(i-1)] <- -phi
    }
  } else if (p >= 2) {
    for(i in 2:p) {
      lvl_i_result = estimate_phi(time_series, p=i-1, b=0.04)
      betas = lvl_i_result$phi_hat
      sigma_i_sq = lvl_i_result$sigma_hat_2[1]

      # Compute gamma values
      L_k[i, 1:(i-1)] <- -betas
      lambda_inv[i] = 1/sigma_i_sq
    }
    lvl_p_result = estimate_phi(time_series, p=p, b=0.04)
    phi = lvl_p_result$phi_hat
    sigma_p_sq = lvl_p_result$sigma_hat_2[1]
    lambda_inv[(p+1):n] = 1/lvl_p_result$sigma_t_2

    for (i in (p+1):n) {
      L_k[i, (i-p):(i-1)] <- -phi
    }
  }

  Lambda_inv = diag(lambda_inv)
  sigma_matrix = L_k %*% Lambda_inv %*% t(L_k)
  return(sigma_matrix)
}
