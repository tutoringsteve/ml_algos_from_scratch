rbf_kernel <- function(x, y, sigma, length_scale) {
  # x and y are vectors of inputs
  sigma^2 * exp(- (x - y)^2 / (2 * length_scale^2))
}

matern <- function(x, y, nu, length_scale) {
  if( ! nu %in% c(1/2, 3/2, 5/2)) stop("only nu = 1/2, 3/2, 5/2 developed so far.")
  
  
  r <- abs(x - y)
  
  # see pages 84-86 from Gaussian Processes for Machine Learning MIT Press 2006 
  # C. E. Rasmussen & C. K. I. Williams for the half-integer solutions
  
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


compute_cov_matrix <- function(x1, x2= x1, kernel = "rbf", parameters) {
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(nrow = n1, ncol = n2)
  for (i in 1:n1) {
    for (j in 1:n2) {
      if(kernel == "rbf") K[i, j] <- rbf_kernel(x1[i], x2[j], parameters$sigma, parameters$length_scale)
      if(kernel == "matern") K[i, j] <- matern(x1[i], x2[j], parameters$nu, parameters$length_scale)
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
  K_inv <- solve(K + sigma_n^2 * diag(1, nrow = nrow(K)))
  
  # P(y|D,x) ~ N(k_star^T * K^(-1) * y, K_starstar * K^(-1) * K_star)
  # so mu = k_star^T * K^(-1) * y
  # and cov =  K_starstar * K^(-1) * K_star
  
  mu <- K_s_ %*% K_inv %*% y_train
  
  cov <- K_ss - K_s_ %*% K_inv %*% K__s
  
  list(mean = mu, covariance = cov)
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
  L <- t(chol(K + sigma_n^2 * diag(1, nrow = nrow(K))))
  Linv = solve(L)
  alpha <- t(Linv)%*%(Linv%*%y_train)
  fstarmean <- t(K__s) %*% alpha
  v <-  Linv %*% K__s
  var_fstar <- K_ss - t(v)%*%v
  log_p_y_given_x <- -.5*t(y_train) %*%  alpha - sum(log(get_diagonal(L))) - 0.5*length(x_train)*log(2*pi)
  
  list(mean = fstarmean, covariance = var_fstar, log_marginal_likelihood = log_p_y_given_x)
}


get_diagonal <- function(M) {
  if(dim(M)[1] != dim(M)[2]) stop("Not a square matrix")
  
  n <- nrow(M)
  diagonal <- numeric(n)
  
  for(i in 1:n) {
    diagonal[i] <- M[i,i]
  }
  
  return(diagonal)
}

# Example usage
x_train <- seq(-5, 5, length.out = 10)
noise_sd <- 0.1
y_train <- sin(x_train) + rnorm(length(x_train), sd = noise_sd)  # noisy observations
x_test <- seq(-6, 6, length.out = 1200)

results <- gaussian_process(x_train, y_train, x_test, sigma_n = noise_sd, kernel = "rbf", parameters = list(sigma = 1, length_scale = 1))
results2 <- gaussian_process_improved(x_train, y_train, x_test, sigma_n = noise_sd, kernel = "rbf", parameters = list(sigma = 1, length_scale = 1))
results_matern_3halves <- gaussian_process_improved(x_train, y_train, x_test, sigma_n = noise_sd, kernel = "matern", parameters = list(nu = 3/2, length_scale = 1))
results_matern_5halves <- gaussian_process_improved(x_train, y_train, x_test, sigma_n = noise_sd, kernel = "matern", parameters = list(nu = 5/2, length_scale = 1))

results <- results_matern_5halves


sds <- sqrt(get_diagonal(results$covariance))
sds2 <- sqrt(get_diagonal(results2$covariance))

all.equal(sds,sds2)

uppers <- results$mean + 1.96 * sds
lowers <- results$mean - 1.96 * sds
data_test <- tibble(x = as.vector(x_test), 
                    y = as.vector(results$mean), 
                    ymin = as.vector(lowers), 
                    ymax = as.vector(uppers))
data_train <-  tibble(x = x_train, y = y_train, ymin = y_train - 1.96*noise_sd, ymax = y_train + 1.96 * noise_sd)


# Plotting results
ggplot() +
  geom_line(data = data_test, aes(x,y), color = "red", size = 1) +
  geom_point(data = data_train, aes(x,y), color = "black", size = 3) +
  # plot the ribbon chart for the noise used to generate the train distribution
  #geom_ribbon(data = data_train, aes(x = x, ymin = ymin, ymax = ymax), fill = "green", alpha = 0.2)+
  # Confidence band for the trained GP
  geom_ribbon(data = data_test, aes(x = x, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4) +
  scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
  labs(title = "Gaussian Process Regression for y = sin(x) + N(0,0.1^2) for 10 train points\nThe blue indicates 95% confidence interval at each point.") +
  theme_classic()
  