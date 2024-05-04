library(pracma) #for the error function erf(x)

source("gaussian_process_regression.R")

parameters_to_string <- function(parameters_list) {
  paste0(sapply(names(parameters$rbf), function(param_name) {
    paste(param_name, "=", parameters$rbf[[param_name]], sep = " ")
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

noise_sd_gp <- 0
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

sds <- list(
  rbf = sqrt(get_diagonal(results$covariance)),
  rbf2 = sqrt(get_diagonal(results2$covariance)),
  matern_1half = sqrt(get_diagonal(results_matern_1half$covariance)),
  matern_3halves = sqrt(get_diagonal(results_matern_3halves$covariance)),
  matern_5halves = sqrt(get_diagonal(results_matern_5halves$covariance)),
  nn_erf = sqrt(get_diagonal(results_nn_erf$covariance)),
  periodic = sqrt(get_diagonal(results_periodic$covariance))
)

lowers <- list(
  rbf = results$mean - 1.96 * sds$rbf,
  rbf2 = results2$mean - 1.96 * sds$rbf2,
  matern_1half = results_matern_1half$mean - 1.96 * sds$matern_1half,
  matern_3halves = results_matern_3halves$mean - 1.96 * sds$matern_3halves,
  matern_5halves = results_matern_5halves$mean - 1.96 * sds$matern_5halves,
  nn_erf = results_nn_erf$mean - 1.96 * sds$nn_erf,
  periodic = results_periodic$mean - 1.96 * sds$periodic
)

uppers <- list(
  rbf = results$mean + 1.96 * sds$rbf,
  rbf2 = results2$mean + 1.96 * sds$rbf2,
  matern_1half = results_matern_1half$mean + 1.96 * sds$matern_1half,
  matern_3halves = results_matern_3halves$mean + 1.96 * sds$matern_3halves,
  matern_5halves = results_matern_5halves$mean + 1.96 * sds$matern_5halves,
  nn_erf = results_nn_erf$mean + 1.96 * sds$nn_erf,
  periodic = results_periodic$mean + 1.96 * sds$periodic
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
  rbf = tibble(x = x_train, y = y_train$rbf, ymin = y_train$rbf - 1.96*noise_sd, ymax = y_train$rbf + 1.96 * noise_sd),
  rbf2 = tibble(x = x_train, y = y_train$rbf2, ymin = y_train$rbf2 - 1.96*noise_sd, ymax = y_train$rbf2 + 1.96 * noise_sd),
  matern_1half = tibble(x = x_train, y = y_train$matern_1half, ymin = y_train$matern_1half - 1.96*noise_sd, ymax = y_train$matern_1half + 1.96 * noise_sd),
  matern_3halves = tibble(x = x_train, y = y_train$matern_3halves, ymin = y_train$matern_3halves - 1.96*noise_sd, ymax = y_train$matern_3halves + 1.96 * noise_sd),
  matern_5halves = tibble(x = x_train, y = y_train$matern_5halves, ymin = y_train$matern_5halves - 1.96*noise_sd, ymax = y_train$matern_5halves + 1.96 * noise_sd),
  nn_erf = tibble(x = x_train, y = y_train$nn_erf, ymin = y_train$nn_erf - 1.96*noise_sd, ymax = y_train$nn_erf + 1.96 * noise_sd),
  periodic = tibble(x = x_train, y = y_train$periodic, ymin = y_train$periodic - 1.96*noise_sd, ymax = y_train$periodic + 1.96 * noise_sd)
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
