library(stringr)
library(tidymodels)
library(scales)
tidymodels_prefer()

#Defines the Lp norm, variable p, where p = 2 is Euclidean distance.
Lp <- function(x,y, p=2) {
  
  #tests for Lp
  if(!is.numeric(x)) { 
    stop("x is not numeric!")
  }
  
  if(!is.numeric(y)) { 
    stop("y is not numeric!")
  }
  
  if(length(x) != length(y)) {
    stop("x and y have different sizes!")
  }
  
  if(length(x) < 1) {
    stop("x length less than 1")
  }
  
  if(length(y) < 1) {
    stop("y length less than 1")
  }
     
  if(!is.numeric(p)) {
    stop("p is non-numeric!")
  } else {
    if(p <= 0) {
      stop("p is not positive!")
    }
  }
  
  sum(abs(x-y)^p)^(1/p)
}

# return object is a list with the scaled data such that mean = 0 and SD = 1, 
# along with the original column means and column standard deviations
# Works on matrices, tibbles, and data.frame's but knn algo is yet to be this robust
normalize_and_center <- function(M) {
  if(is.matrix(M)) {
    column_means <- colMeans(M, na.rm = T)
    column_sds <- apply(M, MARGIN = 2, FUN = sd, na.rm = T)
  } else if(is_tibble(M) || is.data.frame(M)) {
    column_means <- sapply(M, function(col) mean(col, na.rm = T))
    column_sds <- sapply(M, function(col) sd(col, na.rm = T))
  }
  
  normalized <- M
  for(i in 1:ncol(M)) {
    normalized[,i] <- (normalized[,i] - column_means[i])/column_sds[i]
  }

  return(list(normalized = normalized, centers = column_means, deviations = column_sds))
}

revert_normalization <- function(normalized, centers, deviations) {
  M <- normalized
  for(i in 1:ncol(M)) {
    M[,i] <- normalized[,i]*deviations[i] + centers[i]
  }
  
  return(M)
}

# applies knn algorithm to a single test data point
# only returns the indices in train of the k nearest neighbors to test
# 
# k determines number of nearest neighbors to look for
# normalize determines whether or not train and test are mean centered on 0/ and SD scaled to 1 when applying knn
# p determines the type of Lp metric
# returns original test data with nn and ypred columns appended.
knn_single_point <- function(train, test, k = 3, normalize = T, p = 2) {
  if(!is.numeric(k)) {
    stop(paste0("k has to be an integer value! Received k of type ", typeof(k), "."))
  }
  
  if(ceiling(k) != floor(k)) {
    stop("k, the number of neighbors, has to be an integer value!\nReceived k has a non-zero decimal component.")
  }
  
  if(k < 1) {
    stop("k, the number of neighbors, has to be greater than  or equal to 1!")
  }
  
  if(k > nrow(train)) {
    warning(paste0("Value of k set greater than number of training points.\n",
                   "k is being truncated to number of rows, and all indices are being returned."))
    return(1:nrow(train))
  }
  
  if(is.matrix(train)) {
    stop("matrix train not handled yet! come back soon :)")
  }
  # train = xy
  # test = test
  # p = 2
  # k = 3
  
  if(normalize) {
    normalized_and_centered <- normalize_and_center(train)
    train <- normalized_and_centered$normalized
    centers <- normalized_and_centered$centers
    deviations <- normalized_and_centered$deviations
    test <- (test - centers)/deviations
  } 
  
  distances <- c()
  for(i in 1:nrow(train)) {
    distances <- c(distances, Lp(train[i,] %>% as.numeric(), test %>% as.numeric(), p))
  }
  
  #returns indices of rows corresponding to the nearest neighbors.
  tibble(
    row_number = 1:nrow(train),
    distances = distances
  ) %>% 
    arrange(distances) %>%
    slice_head(n = k) %>%
    .$row_number
}

# applies knn algorithm to a single test data point
# returns the indices in train of the k nearest neighbors to test and the 
# non-normalized distances
# 
# k determines number of nearest neighbors to look for
# normalize determines whether or not train and test are mean centered on 0/ and SD scaled to 1 when applying knn
# p determines the type of Lp metric
# returns original test data with nn and ypred columns appended.
knn_single_point_nn_indices_and_distances <- function(train, test, k = 3, p = 2) {
  if(!is.numeric(k)) {
    stop(paste0("k has to be an integer value! Received k of type ", typeof(k), "."))
  }
  
  if(ceiling(k) != floor(k)) {
    stop("k, the number of neighbors, has to be an integer value!\nReceived k has a non-zero decimal component.")
  }
  
  if(k < 1) {
    stop("k, the number of neighbors, has to be greater than  or equal to 1!")
  }
  
  if(k > nrow(train)) {
    warning(paste0("Value of k set greater than number of training points.\n",
                   "k is being truncated to number of rows, and all indices are being returned."))
    return(1:nrow(train))
  }
  
  if(is.matrix(train)) {
    stop("matrix train not handled yet! come back soon :)")
  }
  
  distances <- c()
  for(i in 1:nrow(train)) {
    distances <- c(distances, Lp(train[i,] %>% as.numeric(), test %>% as.numeric(), p))
  }
  
  #returns indices of rows corresponding to the nearest neighbors.
  tibble(
    row_number = 1:nrow(train),
    distances = distances
  ) %>% 
    arrange(distances) %>%
    slice_head(n = k) %>%
    select(row_number, distances)
}

# return the single mode of a categorical vector, or return a random value among the ties
vote_best <- function(y) {
  y <- y
  freq_table_y <- table(y)
  best <- freq_table_y %>%
    as_tibble() %>%
    filter(n == max(freq_table_y))
    
  return(sample(best[,1], 1)[[1,1]])
}

# convert a single numeric vector into a tibble for the code to work
# We know train is a tibble with the same column meanings as test
# so get any row of train to form a tibble with the right columns
# then for each column i, populate the slice with the ith index of the
# vector test.
# 
# Does nothing if test is a list/data.frame/tibble/etc
if_necessary_convert_test_from_vector_to_tibble <- function(train, features, test) {
  if(typeof(test) != "list") {
    tibble_storage_for_test <- train %>%
      # make sure test has right feature set.
      select(all_of(features)) %>%
      slice(1)
    for(i in 1:length(test)) {
      tibble_storage_for_test[1,i] <- test[i]
    }
    test <- tibble_storage_for_test
  }
  
  return(test)
}

# applies knn algorithm to a test data set that can contain multiple points
# knn is applied to each point in test, against all the points in train
# normalize determines whether or not train and test are mean centered on 0/ and SD scaled to 1 when applying knn
# p determines the type of Lp metric
# returns original test data with nn and ypred columns appended.
knn <- function(train, test, features = names(train), target, categorical = T, k = 3, normalize = T, p = 2) {
 
  if(sum(!(names(test) %in% names(train))) > 0) {
    stop("The test set needs to have column names present in the training data!")
  }
  
  test <- if_necessary_convert_test_from_vector_to_tibble(train, features, test)
  
  nn <- list() #initialize empty for loop
  distances <- list() #initialize empty for loop
  ypred <- c() #initialize empty for loop
  for(i in 1:nrow(test)) {
    
    test_point <- test %>% slice(i)
    
    nn_i <- knn_single_point(train = train[features], test = test_point, k = k, normalize = normalize, p = p)
    nn <- c(nn, list(nn_i))
    
    distances_i <- c()
    for(j in 1:length(nn_i)) {
      distances_i <- c(distances_i, Lp(as.numeric(train[features][nn_i[j],]), as.numeric(test_point), p))
    }
    
    distances <- c(distances, list(distances_i))
    
    y <- train[nn[[i]],][target]
    if(categorical) {
      ypred <- c(ypred, vote_best(y))
    } else {
      ypred <- c(ypred, mean(y, na.rm = T))
    }
  }
  
  # return original test data with nn and ypred columns appended.
  test %>%
    bind_cols(
      tibble(
        ypred = ypred,
        nn = nn,
        distances = distances
      )
    )
}

################# code for decision boundary calculation #############
#expands a grid with num_points along each axis, within the described rectangular region
#returns a tibble with columns x and y
get_grid <- function(num_points, min_x, max_x, min_y, max_y) {
  range_x <- max_x - min_x
  range_y <- max_y - min_y
  xrange <- 0:num_points
  xrange <- xrange/(num_points/range_x)
  xrange <- xrange + min_x
  
  yrange <- 0:num_points
  yrange <- yrange/(num_points/range_y)
  yrange <- yrange + min_y
  
  expand_grid(
    x = xrange,
    y = yrange
  )
}

#runs the knn algorithm on the training data using a grid of points fine enough to plot a decision boundary without any gaps in the graph
# returns a 
get_knn_decision_boundary <- function(num_points = 50, min_x = -5.0, max_x = 5.0, min_y = min_x, max_y = max_x, 
                                      train, features = names(train), target, 
                                      categorical = T, k = 3, normalize = T, p = 2) {
  test_grid <- get_grid(num_points, min_x, max_x, min_x, max_x)
  
  # appends nn and ypred columns for nearest neighbors and predicted y-values to the training data which in this case is the 
  # num_points by num_points grid covering the rectangular region defined by [min_x,max_x] X [min_y,max_y] 
  knn(train=train, features = features, test = test_grid, target = target, categorical = categorical, k = k, normalize = normalize, p = p)
}




graph_knn_decision_boundary <- function(decision_boundary) {
  # decision_boundary <- get_knn_decision_boundary(num_points = num_points, min_x = min_x, max_x = max_x, min_y = min_y, max_y = max_y, 
  #                           train = train, target = target, categorical = categorical, k = k, normalize = normalize, p = p)
  
  
  p <- ggplot(data = decision_boundary) +
        geom_point(aes(x = x, y = y, color = ypred), shape = 15, size = 12) +
        labs(title = "Decision Boundary") +
        theme(legend.position="bottom", legend.text = element_text(size = 12 + size_label_text)) +
        guides(color=guide_legend(title="Class"))
  
  list(decision_boundary = decision_boundary, graph = p) %>%
    return()

    }
