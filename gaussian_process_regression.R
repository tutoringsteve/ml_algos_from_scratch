rbf_kernel <- function(x, y, sigma, length_scale) {
  # x and y are vectors of inputs
  sigma^2 * exp(- (x - y)^2 / (2 * length_scale^2))
}


compute_cov_matrix <- function(x1, x2= x1, sigma, length_scale) {
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(nrow = n1, ncol = n2)
  for (i in 1:n1) {
    for (j in 1:n2) {
      K[i, j] <- rbf_kernel(x1[i], x2[j], sigma, length_scale)
    }
  }
  return(K)
}

gaussian_process <- function(x_train, y_train, x_test, sigma, length_scale, sigma_n) {
  #sigma <- 1
  #length_scale <- 1
  #sigma_n <- 0.1
  
  K <- compute_cov_matrix(x_train, sigma = sigma, length_scale = length_scale)
  K__s <- compute_cov_matrix(x_train, x_test, sigma = sigma, length_scale = length_scale)
  K_s_ <- compute_cov_matrix(x_test, x_train, sigma = sigma, length_scale = length_scale)
  K_ss <- compute_cov_matrix(x_test, sigma = sigma, length_scale = length_scale)
  K_inv <- solve(K)
  
  # P(y|D,x) ~ N(k_star^T * K^(-1) * y, K_starstar * K^(-1) * K_star)
  # so mu = k_star^T * K^(-1) * y
  # and cov =  K_starstar * K^(-1) * K_star
  
  mu <- K_s_ %*% K_inv %*% y_train
  
  cov <- K_ss - K_s_ %*% K_inv %*% K__s
  
  list(mean = mu, covariance = cov)
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

results <- gaussian_process(x_train, y_train, x_test, sigma = 1, length_scale = 1, sigma_n = 0.1)



sds <- sqrt(get_diagonal(results$covariance))

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
  