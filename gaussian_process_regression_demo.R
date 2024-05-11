library(pracma) #for the error function erf(x)

source("gaussian_process_regression.R")

parameters_to_string <- function(parameters_list) {
  paste0(sapply(names(parameters_list), function(param_name) {
    paste(param_name, "=", round(parameters_list[[param_name]],4), sep = " ")
  }, simplify = TRUE), collapse = "; ")
}

# Example usage
x_train <- seq(-5, 5, length.out = 10)
noise_sd <- 0.1
y_train <- tibble(
  rbf = sin(x_train) + rnorm(length(x_train), sd = noise_sd),  # noisy observations
  rbf2 = sin(x_train) + rnorm(length(x_train), sd = noise_sd),
  matern_1half = sin(x_train) + rnorm(length(x_train), sd = noise_sd),
  matern_3halves = sin(x_train) + rnorm(length(x_train), sd = noise_sd),
  matern_5halves = sin(x_train) + rnorm(length(x_train), sd = noise_sd),
  nn_erf = erf(x_train) + rnorm(length(x_train), sd = noise_sd), #make it similar to the underlying transfer function
  periodic = sin(x_train) + rnorm(length(x_train), sd = noise_sd),
)

y_train$rbf2 <- y_train$rbf # so they have the same output
  
erf_data <- tibble(x = seq(-5, 5, 0.01), y = erf(seq(-5, 5, 0.01)))
sin_data <- tibble(x = seq(-5, 5, 0.01), y = sin(seq(-5, 5, 0.01)))

x_test <- seq(-6, 6, length.out = 600)

parameters <- list( 
  rbf = list(sigma = 1, length_scale = 1),
  rbf2 = list(sigma = 1, length_scale = 1),
  matern_1half = list(nu = 1/2, length_scale = 1),
  matern_3halves = list(nu = 3/2, length_scale = 1),
  matern_5halves = list(nu = 5/2, length_scale = 1),
  nn_erf = list(sigma = 1, sigma_offset = 1, sigma_x = 10),
  periodic = list(sigma = 1, length_scale = 1)
)

noise_sd_gp <- 0.1
sigma_n <- list(rbf = noise_sd_gp,
                rbf2 = noise_sd_gp,
                matern_1half = noise_sd_gp,
                matern_3halves = noise_sd_gp,
                matern_5halves = noise_sd_gp,
                nn_erf = noise_sd_gp,
                periodic = noise_sd_gp
)

results <- gaussian_process(x_train, y_train$rbf, x_test, sigma_n = sigma_n$rbf, kernel = "rbf", parameters = parameters$rbf)
results2 <- gaussian_process_improved(x_train, y_train$rbf2, x_test, sigma_n = sigma_n$rbf2, kernel = "rbf", parameters = parameters$rbf2)
results_matern_1half <- gaussian_process_improved(x_train, y_train$matern_1half, x_test, sigma_n = sigma_n$matern_1half, kernel = "matern", parameters = parameters$matern_1half)
#results_exponential_covariance <- results_matern_1half
results_matern_3halves <- gaussian_process_improved(x_train, y_train$matern_3halves, x_test, sigma_n = sigma_n$matern_3halves, kernel = "matern", parameters = parameters$matern_3halves)
results_matern_5halves <- gaussian_process_improved(x_train, y_train$matern_5halves, x_test, sigma_n = sigma_n$matern_5halves, kernel = "matern", parameters = parameters$matern_5halves)
results_nn_erf <- gaussian_process_improved(x_train, y_train$nn_erf, x_test, sigma_n = sigma_n$nn_erf, kernel = "nn", parameters = parameters$nn_erf)
results_periodic <- gaussian_process_improved(x_train, y_train$periodic, x_test, sigma_n = sigma_n$periodic, kernel = "periodic", parameters = parameters$periodic)

means <- list(
  rbf = results$mean,
  rbf2 = results2$mean,
  matern_1half = results_matern_1half$mean,
  matern_3halves = results_matern_3halves$mean,
  matern_5halves = results_matern_5halves$mean,
  nn_erf = results_nn_erf$mean,
  periodic = results_periodic$mean
)

sds <- list(
  rbf = sqrt(diag(results$covariance)),
  rbf2 = sqrt(diag(results2$covariance)),
  matern_1half = sqrt(diag(results_matern_1half$covariance)),
  matern_3halves = sqrt(diag(results_matern_3halves$covariance)),
  matern_5halves = sqrt(diag(results_matern_5halves$covariance)),
  nn_erf = sqrt(diag(results_nn_erf$covariance)),
  periodic = sqrt(diag(results_periodic$covariance))
)

lowers <- list(
  rbf = means$rbf - 1.96 * sds$rbf,
  rbf2 = means$rbf2 - 1.96 * sds$rbf2,
  matern_1half = means$matern_1half - 1.96 * sds$matern_1half,
  matern_3halves = means$matern_3halves - 1.96 * sds$matern_3halves,
  matern_5halves = means$matern_5halves - 1.96 * sds$matern_5halves,
  nn_erf = means$nn_erf - 1.96 * sds$nn_erf,
  periodic = means$periodic - 1.96 * sds$periodic
)

uppers <- list(
  rbf = means$rbf + 1.96 * sds$rbf,
  rbf2 = means$rbf2 + 1.96 * sds$rbf2,
  matern_1half = means$matern_1half + 1.96 * sds$matern_1half,
  matern_3halves = means$matern_3halves + 1.96 * sds$matern_3halves,
  matern_5halves = means$matern_5halves + 1.96 * sds$matern_5halves,
  nn_erf = means$nn_erf + 1.96 * sds$nn_erf,
  periodic = means$periodic + 1.96 * sds$periodic
)

all.equal(sds$rbf,sds$rbf2)


data_test <- list(  
  rbf            = tibble(x = as.vector(x_test), y = as.vector(results$mean),                ymin = as.vector(lowers$rbf),            ymax = as.vector(uppers$rbf)),
  rbf2           = tibble(x = as.vector(x_test), y = as.vector(results2$mean),               ymin = as.vector(lowers$rbf2),           ymax = as.vector(uppers$rbf2)),
  matern_1half   = tibble(x = as.vector(x_test), y = as.vector(results_matern_1half$mean),   ymin = as.vector(lowers$matern_1half),   ymax = as.vector(uppers$matern_1half)),
  matern_3halves = tibble(x = as.vector(x_test), y = as.vector(results_matern_3halves$mean), ymin = as.vector(lowers$matern_3halves), ymax = as.vector(uppers$matern_3halves)),
  matern_5halves = tibble(x = as.vector(x_test), y = as.vector(results_matern_5halves$mean), ymin = as.vector(lowers$matern_5halves), ymax = as.vector(uppers$matern_5halves)),
  nn_erf         = tibble(x = as.vector(x_test), y = as.vector(results_nn_erf$mean),         ymin = as.vector(lowers$nn_erf),         ymax = as.vector(uppers$nn_erf)),
  periodic       = tibble(x = as.vector(x_test), y = as.vector(results_periodic$mean),       ymin = as.vector(lowers$periodic),       ymax = as.vector(uppers$periodic))
)

data_train <-  list(
  rbf = tibble(x = x_train, y = y_train$rbf, ymin = y_train$rbf - 1.96*sigma_n$rbf, ymax = y_train$rbf + 1.96 * sigma_n$rbf),
  rbf2 = tibble(x = x_train, y = y_train$rbf2, ymin = y_train$rbf2 - 1.96*sigma_n$rbf2, ymax = y_train$rbf2 + 1.96 * sigma_n$rbf2),
  matern_1half = tibble(x = x_train, y = y_train$matern_1half, ymin = y_train$matern_1half - 1.96*sigma_n$matern_1half, ymax = y_train$matern_1half + 1.96 * sigma_n$matern_1half),
  matern_3halves = tibble(x = x_train, y = y_train$matern_3halves, ymin = y_train$matern_3halves - 1.96*sigma_n$matern_3halves, ymax = y_train$matern_3halves + 1.96 * sigma_n$matern_3halves),
  matern_5halves = tibble(x = x_train, y = y_train$matern_5halves, ymin = y_train$matern_5halves - 1.96*sigma_n$matern_5halves, ymax = y_train$matern_5halves + 1.96 * sigma_n$matern_5halves),
  nn_erf = tibble(x = x_train, y = y_train$nn_erf, ymin = y_train$nn_erf - 1.96*sigma_n$nn_erf, ymax = y_train$nn_erf + 1.96 * sigma_n$nn_erf),
  periodic = tibble(x = x_train, y = y_train$periodic, ymin = y_train$periodic - 1.96*sigma_n$periodic, ymax = y_train$periodic + 1.96 * sigma_n$periodic)
)

title_text <- list(
  rbf = paste0("Gaussian Process Regression for y = sin(x) + N(0,0.1^2) for ", length(x_train)," train points\nThe blue indicates 95% confidence interval at each point.",
               "\nRBF Kernel with parameters: ", parameters_to_string(parameters$rbf), "; noise = ", sigma_n$rbf, "."),
  rbf2 = paste0("Gaussian Process Regression for y = sin(x) + N(0,0.1^2) for ", length(x_train)," train points\nThe blue indicates 95% confidence interval at each point.",
                "\nRBF Kernel with parameters: ", parameters_to_string(parameters$rbf2), "; noise = ", sigma_n$rbf2, "."),
  matern_1half = paste0("Gaussian Process Regression for y = sin(x) + N(0,0.1^2) for ", length(x_train)," train points\nThe blue indicates 95% confidence interval at each point.",
                        "\nMatern Kernel with parameters: ", parameters_to_string(parameters$matern_1half), "; noise = ", sigma_n$matern_1half, "." ,"\nAKA exponetial covariance kernel."),
  matern_3halves = paste0("Gaussian Process Regression for y = sin(x) + N(0,0.1^2) for ", length(x_train)," train points\nThe blue indicates 95% confidence interval at each point.",
                          "\nMatern Kernel with parameters: ", parameters_to_string(parameters$matern_3halves), "; noise = ", sigma_n$matern_3halves, "."),
  matern_5halves = paste0("Gaussian Process Regression for y = sin(x) + N(0,0.1^2) for ", length(x_train)," train points\nThe blue indicates 95% confidence interval at each point.",
                          "\nMatern Kernel with parameters: ", parameters_to_string(parameters$matern_3halves), "; noise = ", sigma_n$matern_5halves, "."),
  nn_erf = paste0("Gaussian Process Regression for y = erf(x) + N(0,0.1^2) for ", length(x_train)," train points\nThe blue indicates 95% confidence interval at each point.",
                  "\nNeural Net with erf transer function h Kernel with parameters: ", parameters_to_string(parameters$nn_erf), "; noise = ", sigma_n$nn_erf, "."),
  periodic = paste0("Gaussian Process Regression for y = sin(x) + N(0,0.1^2) for ", length(x_train)," train points\nThe blue indicates 95% confidence interval at each point.",
                  "\nperiodic Kernel with parameters: ", parameters_to_string(parameters$periodic), "; noise = ", sigma_n$periodic, ".")
)


# Plotting results

graphs <- list(
  rbf = ggplot() +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test$rbf, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$rbf, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    # Confidence band for the trained GP
    geom_ribbon(data = data_test$rbf, aes(x = x, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4) +
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text$rbf) +
    theme_classic(),
  rbf2 = ggplot() +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test$rbf2, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$rbf2, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    # Confidence band for the trained GP
    geom_ribbon(data = data_test$rbf2, aes(x = x, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4) +
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text$rbf2) +
    theme_classic(),
  matern_1half = ggplot() +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test$matern_1half, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$matern_1half, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    # Confidence band for the trained GP
    geom_ribbon(data = data_test$matern_1half, aes(x = x, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4) +
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text$matern_1half) +
    theme_classic(),
  matern_3halves = ggplot() +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test$matern_3halves, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$matern_3halves, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    # Confidence band for the trained GP
    geom_ribbon(data = data_test$matern_3halves, aes(x = x, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4) +
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text$matern_3halves) +
    theme_classic(),
  matern_5halves = ggplot() +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test$matern_5halves, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$matern_5halves, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    # Confidence band for the trained GP
    geom_ribbon(data = data_test$matern_5halves, aes(x = x, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4) +
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text$matern_5halves) +
    theme_classic(),
  nn_erf = ggplot() +
    geom_line(data = erf_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test$nn_erf, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$nn_erf, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    # Confidence band for the trained GP
    geom_ribbon(data = data_test$nn_erf, aes(x = x, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4) +
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text$nn_erf) +
    theme_classic(),
  periodic = ggplot() +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test$periodic, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$periodic, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    # Confidence band for the trained GP
    geom_ribbon(data = data_test$periodic, aes(x = x, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4) +
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text$periodic) +
    theme_classic()
)

graphs$nn_erf

#############optimization experimentation #############

objectives <- list(
  
  rbf = function(hypers) {
    hypers <- as.list(hypers)
    -1*gaussian_process(x_train, y_train$rbf, x_test, sigma_n = sigma_n$rbf, kernel = "rbf", parameters = list(sigma = hypers$sigma, length_scale = hypers$length_scale))$log_marginal_likelihood
  },

  rbf2= function(hypers) {
    hypers <- as.list(hypers)
    -1*gaussian_process_improved(x_train, y_train$rbf2, x_test, sigma_n = sigma_n$rbf2, kernel = "rbf", parameters = list(sigma = hypers$sigma, length_scale = hypers$length_scale))$log_marginal_likelihood
  },
  matern_1half = function(hypers) {
    hypers <- as.list(hypers)
    -1*gaussian_process_improved(x_train, y_train$matern_1half, x_test, sigma_n = sigma_n$matern_1half, kernel = "matern", parameters = list(nu=3/2, length_scale = hypers$length_scale))$log_marginal_likelihood
  },
  matern_3halves = function(hypers) {
    hypers <- as.list(hypers)
    -1*gaussian_process_improved(x_train, y_train$matern_3halves, x_test, sigma_n = sigma_n$matern_3halves, kernel = "matern", parameters = list(nu=3/2, length_scale = hypers$length_scale))$log_marginal_likelihood
  },
  matern_5halves = function(hypers) {
    hypers <- as.list(hypers)
    -1*gaussian_process_improved(x_train, y_train$matern_5halves, x_test, sigma_n = sigma_n$matern_5halves, kernel = "matern", parameters = list(nu=5/2, length_scale = hypers$length_scale))$log_marginal_likelihood
  },
  nn_erf = function(hypers) {
    hypers <- as.list(hypers)
    -1*gaussian_process_improved(x_train, y_train$nn_erf, x_test, sigma_n = sigma_n$nn_erf, kernel = "nn", parameters = list(sigma =hypers$sigma, sigma_offset = hypers$sigma_offset, sigma_x = hypers$sigma_x))$log_marginal_likelihood
  },
  periodic = function(hypers) {
    hypers <- as.list(hypers)
    -1*gaussian_process_improved(x_train, y_train$periodic, x_test, sigma_n = sigma_n$periodic, kernel = "periodic", parameters = list(sigma =hypers$sigma, length_scale = hypers$length_scale))$log_marginal_likelihood
  }
)


initial_hypers <- list(
  rbf = list(sigma = parameters$rbf$sigma, length_scale = parameters$rbf$length_scale),
  rbf2 = list(sigma = parameters$rbf2$sigma, length_scale = parameters$rbf2$length_scale),
  matern_1half = list(length_scale = parameters$matern_1half$length_scale),
  matern_3halves = list(length_scale = parameters$matern_3halves$length_scale),
  matern_5halves = list(length_scale = parameters$matern_5halves$length_scale),
  nn_erf = list(sigma = parameters$nn_erf$sigma, sigma_offset = parameters$nn_erf$sigma_offset, sigma_x = parameters$nn_erf$sigma_x),
  periodic = list(sigma = parameters$periodic$sigma, length_scale = parameters$periodic$length_scale)
)


method <- "Nelder-Mead"

time_start <- Sys.time()
optim_nelder_mead_rbf <- optim(par = initial_hypers$rbf, fn = objectives$rbf, method = "Nelder-Mead")
(time_length_nelder_mead_rbf <- Sys.time() - time_start)

time_start <- Sys.time()
optim_nelder_mead_rbf2 <- optim(par = initial_hypers$rbf2, fn = objectives$rbf2, method = "Nelder-Mead")
(time_length_nelder_mead_rbf2 <- Sys.time() - time_start)

time_start <- Sys.time()
optim_nelder_mead_matern_1half <- optim(par = initial_hypers$matern_1half, fn = objectives$matern_1half, method = "Nelder-Mead")
(time_length_nelder_mead_matern_1half <- Sys.time() - time_start)

time_start <- Sys.time()
optim_nelder_mead_matern_3halves <- optim(par = initial_hypers$matern_3halves, fn = objectives$matern_3halves, method = "Nelder-Mead")
(time_length_nelder_mead_matern_3halves <- Sys.time() - time_start)

time_start <- Sys.time()
optim_nelder_mead_matern_5halves <- optim(par = initial_hypers$matern_5halves, fn = objectives$matern_5halves, method = "Nelder-Mead")
(time_length_nelder_mead_matern_5halves <- Sys.time() - time_start)

#time_start <- Sys.time()
#optim_nelder_mead_nn_erf <- optim(par = initial_hypers$nn_erf, fn = objectives$nn_erf, method = "Nelder-Mead")
#(time_length_nelder_mead_nn_erf <- Sys.time() - time_start)

time_start <- Sys.time()
optim_nelder_mead_periodic <- optim(par = initial_hypers$periodic, fn = objectives$periodic, method = "Nelder-Mead")
(time_length_nelder_mead_periodic <- Sys.time() - time_start)



#time_start <- Sys.time()
#optim_CG_rbf <- optim(par = initial_hypers$rbf, fn = objectives$rbf, method = "CG")
#(time_length_cgm_rbf <- Sys.time() - time_start)
#
#time_start <- Sys.time()
#optim_CG_rbf2 <- optim(par = initial_hypers$rbf2, fn = objectives$rbf2, method = "CG")
#(time_length_cgm_rbf2 <- Sys.time() - time_start)


results$log_marginal_likelihood
optim_nelder_mead_rbf$value
parameters$rbf
optim_nelder_mead_rbf$par

results2$log_marginal_likelihood
optim_nelder_mead_rbf2$value
parameters$rbf2
optim_nelder_mead_rbf2$par

results_matern_1half$log_marginal_likelihood
optim_nelder_mead_matern_1half$value
parameters$matern_1half
optim_nelder_mead_matern_1half$par

results_matern_3halves$log_marginal_likelihood
optim_nelder_mead_matern_3halves$value
parameters$matern_3halves
optim_nelder_mead_matern_3halves$par

results_matern_5halves$log_marginal_likelihood
optim_nelder_mead_matern_5halves$value
parameters$matern_5halves
optim_nelder_mead_matern_5halves$par

#results_nn_erf$log_marginal_likelihood
#optim_nelder_mead_nn_erf$value
#parameters$nn_erf
#optim_nelder_mead_nn_erf$par

results_periodic$log_marginal_likelihood
optim_nelder_mead_periodic$value
parameters$periodic
optim_nelder_mead_periodic$par

results_optimized <- list(
  rbf = gaussian_process_improved(x_train, y_train$rbf, x_test, sigma_n = sigma_n$rbf, kernel = "rbf", parameters = as.list(optim_nelder_mead_rbf$par)),
  rbf2 = gaussian_process_improved(x_train, y_train$rbf2, x_test, sigma_n = sigma_n$rbf2, kernel = "rbf", parameters = as.list(optim_nelder_mead_rbf2$par)),
  matern_1half = gaussian_process_improved(x_train, y_train$matern_1half, x_test, sigma_n = sigma_n$matern_1half, kernel = "matern", parameters = c(list(nu = 1/2), as.list(optim_nelder_mead_matern_1half$par))),
  matern_3halves = gaussian_process_improved(x_train, y_train$matern_3halves, x_test, sigma_n = sigma_n$matern_3halves, kernel = "matern", parameters = c(list(nu = 3/2), as.list(optim_nelder_mead_matern_3halves$par))),
  matern_5halves = gaussian_process_improved(x_train, y_train$matern_5halves, x_test, sigma_n = sigma_n$matern_5halves, kernel = "matern", parameters = c(list(nu = 5/2), as.list(optim_nelder_mead_matern_5halves$par))),
  #nn_erf = gaussian_process_improved(x_train, y_train$nn_erf, x_test, sigma_n = sigma_n$nn_erf, kernel = "nn", parameters = as.list(optim_nelder_mead_nn_erf$par)),
  periodic = gaussian_process_improved(x_train, y_train$periodic, x_test, sigma_n = sigma_n$periodic, kernel = "periodic", parameters = as.list(optim_nelder_mead_periodic$par))
)

sds_optimized <- list(
  rbf = sqrt(diag(results_optimized$rbf$covariance)),
  rbf2 = sqrt(diag(results_optimized$rbf2$covariance)),
  matern_1half = sqrt(diag(results_optimized$matern_1half$covariance)),
  matern_3halves = sqrt(diag(results_optimized$matern_3halves$covariance)),
  matern_5halves = sqrt(diag(results_optimized$matern_5halves$covariance)),
  nn_erf = sqrt(diag(results_optimized$nn_erf$covariance)),
  periodic = sqrt(diag(results_optimized$periodic$covariance))
)

lowers_optimized <- list(
  rbf = results_optimized$rbf$mean - 1.96 * sds_optimized$rbf,
  rbf2 = results_optimized$rbf2$mean - 1.96 * sds_optimized$rbf2,
  matern_1half = results_optimized$matern_1half$mean - 1.96 * sds$matern_1half,
  matern_3halves = results_optimized$matern_3halves$mean - 1.96 * sds$matern_3halves,
  matern_5halves = results_optimized$matern_5halves$mean - 1.96 * sds$matern_5halves,
  nn_erf = results_optimized$nn_erf$mean - 1.96 * sds$nn_erf,
  periodic = results_optimized$periodic$mean - 1.96 * sds$periodic
)

uppers_optimized <- list(
  rbf = results_optimized$rbf$mean + 1.96 * sds_optimized$rbf,
  rbf2 = results_optimized$rbf2$mean + 1.96 * sds_optimized$rbf2,
  matern_1half = results_optimized$matern_1half$mean + 1.96 * sds$matern_1half,
  matern_3halves = results_optimized$matern_3halves$mean + 1.96 * sds$matern_3halves,
  matern_5halves = results_optimized$matern_5halves$mean + 1.96 * sds$matern_5halves,
  nn_erf = results_optimized$nn_erf$mean + 1.96 * sds$nn_erf,
  periodic = results_optimized$periodic$mean + 1.96 * sds$periodic
)

all.equal(sds$rbf,sds$rbf2)


data_test_optimized <- list(  
  rbf            = tibble(x = as.vector(x_test), y = as.vector(results_optimized$rbf$mean),            ymin = as.vector(lowers_optimized$rbf),            ymax = as.vector(uppers_optimized$rbf)),
  rbf2           = tibble(x = as.vector(x_test), y = as.vector(results_optimized$rbf2$mean),           ymin = as.vector(lowers_optimized$rbf2),           ymax = as.vector(uppers_optimized$rbf2)),
  matern_1half   = tibble(x = as.vector(x_test), y = as.vector(results_optimized$matern_1half$mean),   ymin = as.vector(lowers_optimized$matern_1half),   ymax = as.vector(uppers_optimized$matern_1half)),
  matern_3halves = tibble(x = as.vector(x_test), y = as.vector(results_optimized$matern_3halves$mean), ymin = as.vector(lowers_optimized$matern_3halves), ymax = as.vector(uppers_optimized$matern_3halves)),
  matern_5halves = tibble(x = as.vector(x_test), y = as.vector(results_optimized$matern_5halves$mean), ymin = as.vector(lowers_optimized$matern_5halves), ymax = as.vector(uppers_optimized$matern_5halves)),
  nn_erf         = tibble(x = as.vector(x_test), y = as.vector(results_optimized$nn_erf$mean),         ymin = as.vector(lowers_optimized$nn_erf),         ymax = as.vector(uppers_optimized$nn_erf)),
  periodic       = tibble(x = as.vector(x_test), y = as.vector(results_optimized$periodic$mean),       ymin = as.vector(lowers_optimized$periodic),       ymax = as.vector(uppers_optimized$periodic))
)

title_text_optimized <- list(
  rbf = paste0("Gaussian Process Regression for y = sin(x) + N(0,0.1^2) for ", length(x_train)," train points\nThe blue indicates 95% confidence interval at each point.",
               "\nRBF Kernel with OPTIMIZED parameters: ", parameters_to_string(as.list(optim_nelder_mead_rbf$par)), "; noise = ", sigma_n$rbf, "."),
  rbf2 = paste0("Gaussian Process Regression for y = sin(x) + N(0,0.1^2) for ", length(x_train)," train points\nThe blue indicates 95% confidence interval at each point.",
                "\nRBF Kernel with OPTIMIZED parameters [nelder mead]: ", parameters_to_string(as.list(optim_nelder_mead_rbf2$par)), "; noise = ", sigma_n$rbf2, "."),
  matern_1half = paste0("Gaussian Process Regression for y = sin(x) + N(0,0.1^2) for ", length(x_train)," train points\nThe blue indicates 95% confidence interval at each point.",
                        "\nMatern Kernel with OPTIMIZED parameters: ", parameters_to_string(c(list(nu = 1/2), as.list(optim_nelder_mead_matern_1half$par))), "; noise = ", sigma_n$matern_1half, "." ,"\nAKA exponetial covariance kernel."),
  matern_3halves = paste0("Gaussian Process Regression for y = sin(x) + N(0,0.1^2) for ", length(x_train)," train points\nThe blue indicates 95% confidence interval at each point.",
                          "\nMatern Kernel with OPTIMIZED parameters: ", parameters_to_string(c(list(nu = 3/2), as.list(optim_nelder_mead_matern_3halves$par))), "; noise = ", sigma_n$matern_3halves, "."),
  matern_5halves = paste0("Gaussian Process Regression for y = sin(x) + N(0,0.1^2) for ", length(x_train)," train points\nThe blue indicates 95% confidence interval at each point.",
                          "\nMatern Kernel with OPTIMIZED parameters: ", parameters_to_string(c(list(nu = 5/2), as.list(optim_nelder_mead_matern_5halves$par))), "; noise = ", sigma_n$matern_5halves, "."),
  nn_erf = paste0("Gaussian Process Regression for y = erf(x) + N(0,0.1^2) for ", length(x_train)," train points\nThe blue indicates 95% confidence interval at each point.",
                  "\nNeural Net with erf transer function h Kernel.\nOPTIMIZED parameters: ", parameters_to_string(as.list(optim_nelder_mead_nn_erf$par)), "; noise = ", sigma_n$nn_erf, "."),
  periodic = paste0("Gaussian Process Regression for y = sin(x) + N(0,0.1^2) for ", length(x_train)," train points\nThe blue indicates 95% confidence interval at each point.",
                    "\nperiodic Kernel with OPTIMIZED parameters: ", parameters_to_string(as.list(optim_nelder_mead_periodic$par)), "; noise = ", sigma_n$periodic, ".")
)



graphs_optimized <- list(
  rbf = ggplot() +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test_optimized$rbf, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$rbf, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    # Confidence band for the trained GP
    geom_ribbon(data = data_test_optimized$rbf, aes(x = x, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4) +
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text_optimized$rbf) +
    theme_classic(),
  rbf2 = ggplot() +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test_optimized$rbf2, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$rbf2, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    # Confidence band for the trained GP
    geom_ribbon(data = data_test_optimized$rbf2, aes(x = x, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4) +
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text_optimized$rbf2) +
    theme_classic(),
  matern_1half = ggplot() +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test_optimized$matern_1half, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$matern_1half, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    # Confidence band for the trained GP
    geom_ribbon(data = data_test_optimized$matern_1half, aes(x = x, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4) +
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text_optimized$matern_1half) +
    theme_classic(),
  matern_3halves = ggplot() +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test_optimized$matern_3halves, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$matern_3halves, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    # Confidence band for the trained GP
    geom_ribbon(data = data_test_optimized$matern_3halves, aes(x = x, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4) +
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text_optimized$matern_3halves) +
    theme_classic(),
  matern_5halves = ggplot() +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test_optimized$matern_5halves, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$matern_5halves, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    # Confidence band for the trained GP
    geom_ribbon(data = data_test_optimized$matern_5halves, aes(x = x, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4) +
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text_optimized$matern_5halves) +
    theme_classic(),
  nn_erf = ggplot() +
    geom_line(data = erf_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test_optimized$nn_erf, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$nn_erf, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    # Confidence band for the trained GP
    geom_ribbon(data = data_test_optimized$nn_erf, aes(x = x, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4) +
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text_optimized$nn_erf) +
    theme_classic(),
  periodic = ggplot() +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test$periodic, aes(x,y), color = "purple", size = 1, linetype = "dashed") +
    geom_line(data = data_test_optimized$periodic, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$periodic, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    # Confidence band for the trained GP
    geom_ribbon(data = data_test$periodic, aes(x = x, ymin = ymin, ymax = ymax), fill = "purple", alpha = 0.2) +
    geom_ribbon(data = data_test_optimized$periodic, aes(x = x, ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4) +
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text_optimized$periodic) +
    theme_classic()
)
title_text_optimized$rbf2

graphs$rbf
graphs_optimized$rbf

graphs$rbf2
graphs_optimized$rbf2  

graphs$matern_1half
graphs_optimized$matern_1half

graphs$matern_3halves
graphs_optimized$matern_3halves  

graphs$matern_5halves
graphs_optimized$matern_5halves

graphs$nn_erf
graphs_optimized$nn_erf

graphs$periodic
graphs_optimized$periodic

alpha1 = 0.8
alpha2 = 0.8
ribbon_color_1 <- "red"
ribbon_color_2 <- "green"

graphs_compare <- list(
  rbf = ggplot() +
    # Confidence band for the trained GP
    geom_ribbon(data = data_test$rbf, aes(x = x, ymin = ymin, ymax = ymax), fill = ribbon_color_1, alpha = alpha1) +
    geom_ribbon(data = data_test_optimized$rbf, aes(x = x, ymin = ymin, ymax = ymax), fill = ribbon_color_2, alpha = alpha2) +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test$rbf, aes(x,y), color = "purple", size = 1, linetype = "dashed") +
    geom_line(data = data_test_optimized$rbf, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$rbf, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text_optimized$rbf) +
    theme_classic(),
  rbf2 = ggplot() +
    # Confidence band for the trained GP
    geom_ribbon(data = data_test$rbf2, aes(x = x, ymin = ymin, ymax = ymax), fill = ribbon_color_1, alpha = alpha1) +
    geom_ribbon(data = data_test_optimized$rbf2, aes(x = x, ymin = ymin, ymax = ymax), fill = ribbon_color_2, alpha = alpha2) +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test$rbf2, aes(x,y), color = "purple", size = 1, linetype = "dashed") +
    geom_line(data = data_test_optimized$rbf2, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$rbf2, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text_optimized$rbf2) +
    theme_classic(),
  matern_1half = ggplot() +
    # Confidence band for the trained GP
    geom_ribbon(data = data_test$matern_1half, aes(x = x, ymin = ymin, ymax = ymax), fill = ribbon_color_1, alpha = alpha1) +
    geom_ribbon(data = data_test_optimized$matern_1half, aes(x = x, ymin = ymin, ymax = ymax), fill = ribbon_color_2, alpha = alpha2) +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test$matern_1half, aes(x,y), color = "purple", size = 1, linetype = "dashed") +
    geom_line(data = data_test_optimized$matern_1half, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$matern_1half, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text_optimized$matern_1half) +
    theme_classic(),
  matern_3halves = ggplot() +
    # Confidence band for the trained GP
    geom_ribbon(data = data_test$matern_3halves, aes(x = x, ymin = ymin, ymax = ymax), fill = ribbon_color_1, alpha = alpha1) +
    geom_ribbon(data = data_test_optimized$matern_3halves, aes(x = x, ymin = ymin, ymax = ymax), fill = ribbon_color_2, alpha = alpha2) +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test$matern_3halves, aes(x,y), color = "purple", size = 1, linetype = "dashed") +
    geom_line(data = data_test_optimized$matern_3halves, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$matern_3halves, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text_optimized$matern_3halves) +
    theme_classic(),
  matern_5halves = ggplot() +
    # Confidence band for the trained GP
    geom_ribbon(data = data_test$matern_5halves, aes(x = x, ymin = ymin, ymax = ymax), fill = ribbon_color_1, alpha = alpha1) +
    geom_ribbon(data = data_test_optimized$matern_5halves, aes(x = x, ymin = ymin, ymax = ymax), fill = ribbon_color_2, alpha = alpha2) +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test$matern_5halves, aes(x,y), color = "purple", size = 1, linetype = "dashed") +
    geom_line(data = data_test_optimized$matern_5halves, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$matern_5halves, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text_optimized$matern_5halves) +
    theme_classic(),
  nn_erf = ggplot() +
    # Confidence band for the trained GP
    geom_ribbon(data = data_test$nn_erf, aes(x = x, ymin = ymin, ymax = ymax), fill = ribbon_color_1, alpha = alpha1) +
    geom_ribbon(data = data_test_optimized$nn_erf, aes(x = x, ymin = ymin, ymax = ymax), fill = ribbon_color_2, alpha = alpha2) +
    geom_line(data = erf_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test$nn_erf, aes(x,y), color = "purple", size = 1, linetype = "dashed") +
    geom_line(data = data_test_optimized$nn_erf, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$nn_erf, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text_optimized$nn_erf) +
    theme_classic(),
  periodic = ggplot() +
    # Confidence band for the trained GP
    geom_ribbon(data = data_test$periodic, aes(x = x, ymin = ymin, ymax = ymax), fill = ribbon_color_1, alpha = alpha1) +
    geom_ribbon(data = data_test_optimized$periodic, aes(x = x, ymin = ymin, ymax = ymax), fill = ribbon_color_2, alpha = alpha2) +
    geom_line(data = sin_data, aes(x,y), color = "black", size = 1) +
    geom_line(data = data_test$periodic, aes(x,y), color = "purple", size = 1, linetype = "dashed") +
    geom_line(data = data_test_optimized$periodic, aes(x,y), color = "red", size = 1) +
    geom_point(data = data_train$periodic, aes(x,y), color = "black", size = 3) +
    # plot the ribbon chart for the noise used to generate the train distribution
    scale_x_continuous(limits = c(-6,6), expand= c(0,0)) +
    labs(title = title_text_optimized$periodic) +
    theme_classic()
)

graphs_compare
