library(ggplot2)
library(dplyr)
library(tibble)

nn_erf_kernel <- function(x, y, sigma, sigma_offset, sigma_x) {
  
  
  # see pages 91 eq. 4.29 from Gaussian Processes for Machine Learning MIT Press 2006 
  # C. E. Rasmussen & C. K. I. Williams for the half-integer solutions
  
  #print(paste0("length x is ", length(x)))
  #print(paste0("length y is ", length(y)))
  #print(paste0("length sigma_x is ", length(sigma_x)))
  
  if(length(x) != length(y)) stop(paste0("nn_erf_kernel error: mismatched input vector dimensions."))
  
  if(length(sigma_x) != length(x)) {
    if(length(sigma_x) == 1) {
      warning("sigma_x has length 1 and is being repeated to have same dimensions as x and y which are not length 1.")
      sigma_x <- rep(sigma_x, length(x))
    } else {
      stop("nn_erf_kernel error: The length of sigma_x is neither 1 nor equal to the length of x & y.")
    }
  }
  
  #if x or y comes in as a 1 dimensional data.frame or tibble() unlisting will correct things
  #if x or y is a matrix or vector unlist will do nothing.
  x <- matrix(c(1,x) %>% unlist(), ncol = 1) #pre-appending a 1 allows offset to work
  y <- matrix(c(1,y) %>% unlist(), ncol = 1) 
  
  # sigma_offset determines variability of offset
  # sigma_x determines variability in each x-direction
  S <- diag(c(sigma_offset^2, sigma_x^2))
  
  
  xSx <- t(x) %*% S %*% x
  xSy <- t(x) %*% S %*% y
  ySy <- t(y) %*% S %*% y
  
  #k_NN(x,y) =....
  sigma^2 * (2/pi) * asin(2*xSy*((1 + 2*xSx)*(1 + 2*ySy))^(-1/2))
}


rbf_kernel <- function(x, y, sigma, length_scale) {
  
  # x and y are vectors in the feature space
  if(length(x) != length(y)) stop(paste0("rbf_kernel: mismatched input vector dimensions."))
  
  r2 <- t(x - y) %*% (x-y)
  
  
  
  sigma^2 * exp(-r2 / (2*length_scale^2))
}

matern <- function(x, y, nu, length_scale) {
  # see pages 84-86 from Gaussian Processes for Machine Learning MIT Press 2006 
  # C. E. Rasmussen & C. K. I. Williams for the half-integer solutions
  
  if( ! nu %in% c(1/2, 3/2, 5/2)) stop("only nu = 1/2, 3/2, 5/2 developed so far.")
  
  # x and y are vectors in the feature space
  if(length(x) != length(y)) stop(paste0("rbf_kernel: mismatched input vector dimensions."))
  
  r <- sqrt(t(x - y) %*% (x - y))
  
  if(nu == 1/2) {
    #exponential covariance function
    k_r <- exp(-r/length_scale)
  } else if(nu == 3/2) {
    k_r <- (1 + sqrt(3)*r/length_scale)*exp(-sqrt(3)*r/length_scale)
  } else if(nu == 5/2) {
    k_r <- (1 + (sqrt(5)*r/length_scale) + (5*r^2/(3*length_scale^2)))*exp(-sqrt(5)*r/length_scale)
  }
  
  return(k_r)
}

exponential_covariance <- function(x,y, length_scale) {
  return(matern(x = x, y = y, nu = 1/2, length_scale = length_scale))
}

periodic_random <- function(x,y, sigma, length_scale) {
  # see pages 92 eq. 4.31 from Gaussian Processes for Machine Learning MIT Press 2006 
  # C. E. Rasmussen & C. K. I. Williams for the half-integer solutions
  
  # x and y are vectors in the feature space
  if(length(x) != length(y)) stop(paste0("rbf_kernel: mismatched input vector dimensions."))
  
  # get a [periodic] non-stationary kernel by running x -> u(x) and running u through
  # a stationary kernel (in this case, RBF with u(x) = (cos(x),sin(x)) for periodicity)
  
  sigma^2 * exp(-2 * sin(.5 * (x - y))^2 / length_scale^2)
}

compute_cov_matrix <- function(x1, x2= x1, kernel = "rbf", parameters) {
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(nrow = n1, ncol = n2)
  for (i in 1:n1) {
    for (j in 1:n2) {
      if(kernel == "rbf") K[i, j] <- rbf_kernel(x1[i], x2[j], parameters$sigma, parameters$length_scale)
      if(kernel == "matern") K[i, j] <- matern(x1[i], x2[j], parameters$nu, parameters$length_scale)
      if(kernel == "nn") K[i, j] <- nn_erf_kernel(x1[i], x2[j], parameters$sigma, parameters$sigma_offset, parameters$sigma_x)
      if(kernel == "periodic") K[i, j] <- periodic_random(x1[i], x2[j], parameters$sigma, parameters$length_scale)
    }
  }
  return(K)
}

gaussian_process <- function(x_train, y_train, x_test, sigma_n = 0, kernel = "rbf", parameters) {
  #sigma <- 1
  #length_scale <- 1
  #sigma_n <- 0.1
  
  K <- compute_cov_matrix(x_train, kernel = kernel, parameters = parameters)
  #K__s = K(X,X*), K_s_ = K(X*,X) = K__s^T
  K__s <- compute_cov_matrix(x_train, x_test, kernel = kernel, parameters = parameters)
  K_s_ <- compute_cov_matrix(x_test, x_train, kernel = kernel, parameters = parameters)
  K_ss <- compute_cov_matrix(x_test, kernel = kernel, parameters = parameters)
  
  #It's extremely likely that R's solve recognizes this as positive semi-definite
  #and uses cholesky decomposition anyway
  K_y <- K + sigma_n^2 * diag(1, nrow = nrow(K))
  K_inv <- solve(K_y) 
  
  # P(y|D,x) ~ N(k_star^T * K^(-1) * y, K_starstar * K^(-1) * K_star)
  # so mu = k_star^T * K^(-1) * y
  # and cov =  K_starstar * K^(-1) * K_star
  
  mu <- K_s_ %*% K_inv %*% y_train
  
  cov <- K_ss - K_s_ %*% K_inv %*% K__s
  
  
  log_p_y_given_x <- -.5*t(y_train) %*% K_inv %*% y_train - 0.5*log(det(K_y)) - 0.5*length(x_train)*log(2*pi)
  
  list(mean = mu, covariance = cov, log_marginal_likelihood = log_p_y_given_x)
}

gaussian_process_improved <- function(x_train, y_train, x_test, sigma_n = 0, kernel = "rbf", parameters) {
  #Cholesky decomposition is more numerically stable.
  
  # page 19 (section 2.3 Varying the Hyperparameters) 
  # Gaussian Processes for Machine Learning
  # Carl Edward Rasmussen and Christopher K. I. Williams
  # MIT Press, 2006. ISBN-10 0-262-18253-X, ISBN-13 978-0-262-18253-9.
  # https://gaussianprocess.org/gpml/chapters/RW2.pdf
  K <- compute_cov_matrix(x_train, kernel = kernel, parameters = parameters)
  #K__s = K(X,X*), K_s_ = K(X*,X) = K__s^T
  K__s <- compute_cov_matrix(x_train, x_test, kernel = kernel, parameters = parameters)
  K_s_ <- compute_cov_matrix(x_test, x_train, kernel = kernel, parameters = parameters)
  K_ss <- compute_cov_matrix(x_test, kernel = kernel, parameters = parameters)
  
  # R gives the upper triangular matrix, algo expects the lower triangular matrix
  #also possibly could just use chol2inv( ) function on K + sigma_n^2 * diag(1, nrow = nrow(K))
  L <- t(chol(K + sigma_n^2 * diag(1, nrow = nrow(K))))
  Linv = solve(L)
  #alpha is K^(-1) %*% y_train
  alpha <- t(Linv)%*%(Linv%*%y_train)
  fstarmean <- t(K__s) %*% alpha
  v <-  Linv %*% K__s
  var_fstar <- K_ss - t(v)%*%v
  
  # K = LL^T , det(K) = det(L)^2 = product(diag(L))^2
  # -(1/2)*log(det(K)) = -log(sqrt(det(K))) = -log(product(diag(L))) = - sum(log(diag(L)))
  log_p_y_given_x <- -.5*t(y_train) %*%  alpha - sum(log(diag(L))) - 0.5*length(x_train)*log(2*pi)
  
  list(mean = fstarmean, covariance = var_fstar, log_marginal_likelihood = log_p_y_given_x)
}
