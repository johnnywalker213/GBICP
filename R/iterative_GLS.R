#' Iterative Generalized Least Squares Estimation for Time Series
#'
#' This function performs an iterative generalized least squares (GLS) estimation on a given time series.
#' It starts by detrending the time series using a B-spline and ordinary least squares (OLS).
#' An autoregressive (AR) model is then fitted to the detrended series using a specified model selection criterion.
#' The covariance matrix of the detrended series is used to update the spline coefficients using GLS.
#' The process iteratively re-detrends the original time series until the L2 norm of the spline coefficients
#' between consecutive iterations is below a specified threshold or the maximum number of iterations is reached.
#'
#' @param Yt A numeric vector representing the time series data.
#' @param epsilon A numeric value specifying the convergence threshold for the L2 norm of spline coefficients.
#' @param b bandwidth for kernel smoothing when estimating time varying variance of the time series.
#' @param max_iter An integer specifying the maximum number of iterations for the GLS estimation process.
#' @param criterion A string specifying the criterion for AR model selection. Valid options include "aic", "bic","gaic", "gbic", and "gbicp".
#'
#' @return A list containing various components of the GLS estimation process, including the number of iterations (`iter`),
#'         updated spline coefficients (`beta`), AR model coefficients (`phi`), detrended series (`mu`), inverse covariance matrix (`Sigma_inv`),
#'         and other AR model parameters.
#'
#' @examples
#' ts_data <- rnorm(100) # Example time series data
#' gls_result <- iterative_GLS(ts_data) # Perform iterative GLS estimation
#'
#' @export
#'@importFrom splines2 bSpline
#'

iterative_GLS <- function(Yt, epsilon = 1e-6, b = 0.04, max_iter = 1000, criterion = 'gbicp') {

  # Number of data points
  n <- length(Yt)
  x <- (1:n)/n
  # Step 1: Fit a spline for initial guess
  knots <- quantile(x, probs = seq(0, 1, length.out = 7))[-c(1, 7)]
  X <- splines2::bSpline(x, knots = knots, degree = 4, intercept = TRUE)

  # Initialize beta estimates
  beta_old <- solve(t(X)%*%X)%*%t(X)%*%Yt
  mu= X %*% beta_old
  noise_x = Yt - mu

  ## select best AR(p)
  p = select_best_ar_model(time_series = noise_x, max_p = 6, criterion = criterion)[1]

  for (iter in 1:max_iter) {
    #print(iter)

    # Calculate mu

    mu= X %*% beta_old


    # Step 2: Calculate initial error
    noise_x = Yt - mu

    # Step 3: Fit AR(2) model
    lagged_X = create_lagged_matrix(time_series = noise_x, p = p)

    ar_result = estimate_phi(time_series = noise_x, p = p, b = b)
    phi_hat = ar_result$phi_hat
    sigma_hat_2 = ar_result$sigma_hat_2[1]

    noise_est = lagged_X%*%phi_hat
    # Step 4: Create Sigma_inv
    Sigma_inverse = compute_sigma_matrix(time_series = noise_x, p=p)




    # Step 5: Compute new beta estimates
    beta_new = solve(t(X)%*%Sigma_inverse%*%X, t(X)%*%Sigma_inverse%*%Yt)

    # Check for convergence
    if (sqrt(sum((beta_new - beta_old)^2)) < epsilon) {
      break
    }
    #print(sqrt(sum((beta_new - beta_old)^2)))

    # Update beta_old for next iteration
    beta_old <- beta_new
  }
  return(list('iter' = iter, "beta" = beta_new, "phi" = phi_hat, "mu" = X %*% beta_new, 'Sigma_inv' = Sigma_inverse, 'sigma_hat_2' = sigma_hat_2, 'p' = p))
}
