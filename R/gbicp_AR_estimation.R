library(splines)
library(stats)
library(forecast)
library(KernSmooth)


generate_ts <- function(ar_coef = NULL, ma_coef = NULL, n, sigma2 = 1, midpoint = 250) {
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


# Function for creating lagged matrix
create_lagged_matrix <- function(time_series, p) {
  # Get the length of the time series
  n <- length(time_series)

  # Create a matrix where each column is a lagged version of the original series
  lagged_matrix <- sapply(1:p, function(i) c(rep(NA, i), time_series[1:(n-i)]))

  return(lagged_matrix[-c(1:p),])
}


estimate_phi <- function(time_series, p, b) {

  create_lagged_matrix <- function(time_series, p) {
    n <- length(time_series)
    lagged_matrix <- sapply(1:p, function(i) c(rep(NA, i), time_series[1:(n-i)]))
    return(lagged_matrix[-c(1:p),])
  }

  lagged_X = create_lagged_matrix(time_series, p)
  phi_hat = solve(t(lagged_X) %*% lagged_X) %*% t(lagged_X) %*% time_series[-c(1:p)]

  err_est = lagged_X %*% phi_hat

  res = time_series[-c(1:p)] - err_est
  res_sq = res^2

  est = locpoly(x = seq_along(res_sq)/length(res_sq),
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

# Function to estimate AR(p) coefficients and calculate AIC or BIC
estimate_ar_model <- function(time_series, p, criterion, ARMA=T) {
  lagged_matrix <- create_lagged_matrix(time_series, p)
  #print(dim(lagged_matrix))
  n = length(time_series)
  Y <- time_series[-c(1:p)]
  X <- lagged_matrix

  #ar_model <- lm(Y ~ X-1)
  ar_model = Arima(time_series, order=c(p,0,0))
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

# Function to select the best AR(p) model
select_best_ar_model <- function(time_series, max_p, criterion) {
  best_model <- NULL
  best_criterion_value <- 99999

  for (p in 1:max_p) {
    current_model <- estimate_ar_model(time_series, p, criterion)
    if (current_model < best_criterion_value) {
      best_model <- p
      best_criterion_value <- current_model
    }
  }

  return(c(best_model,best_criterion_value))
}


# The iterative procedure function
iterative_GLS <- function(Yt, x, epsilon = 1e-6, b, max_iter = 1000) {

  # Number of data points
  n <- length(x)

  # Step 1: Fit a spline for initial guess
  knots <- quantile(x, probs = seq(0, 1, length.out = 7))[-c(1, 7)]
  X <- bSpline(x, knots = knots, degree = 4, intercept = TRUE)

  # Initialize beta estimates
  beta_old <- solve(t(X)%*%X)%*%t(X)%*%Yt
  mu= X %*% beta_old
  noise_x = Yt - mu

  ## select best AR(p)
  p = select_best_ar_model(time_series = noise_x, max_p = 6, criterion = 'gbicp')[1]

  for (iter in 1:max_iter) {
    #print(iter)

    # Calculate mu

    mu= X %*% beta_old


    # Step 2: Calculate initial error
    noise_x = Yt - mu

    # Step 3: Fit AR(2) model
    lagged_X = create_lagged_matrix(time_series = noise_x, p = p)

    ar_result = estimate_phi(time_series = noise_x, p = p, b = 0.04)
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
  return(list('iter' = iter, "beta" = beta_new, "phi" = phi_hat, "mu" = X %*% beta_new, 'Sigma_inv' = Sigma_inverse, 'phi_1' = phi_1, 'phi_2'=phi_2, 'sigma_hat_2' = sigma_hat_2, 'p' = p))
}

