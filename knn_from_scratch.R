library(stringr)
library(tidymodels)
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


#45-45-90 triangle check
# x<- c(0,1)
# y<- c(1,0)
# Lp(x,y, 2)

#30-60-90 triangle check
# Lp(c(0,0), c(1/2, sqrt(3)/2), 2)

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




###########tests of knn_single_point###############
{
  #double/real valued, but lies on integer grid points from -5 to 5 in x and y directions
  single_point_train_grid <- expand_grid(
    x= (-5:5)*1.0,
    y= (-5:5)*1.0
  )
  
  #randomly assigned real number betwixt the grid points
  single_test_point <- c(runif(2,-3, 3))
  indices_from_single_random_point_tested_on_grid <- knn_single_point(
                                                        train = single_point_train_grid,
                                                        test = single_test_point,
                                                        k = 3, 
                                                        normalize = T, 
                                                        p = 2)
  
  non_normalized_indices_from_single_random_point_tested_on_grid <- knn_single_point(
                                                                       train = single_point_train_grid,
                                                                       test = single_test_point,
                                                                       k = 3, 
                                                                       normalize = F, 
                                                                       p = 2)
  
  indices_and_distances_from_single_random_point_tested_on_grid <- knn_single_point_nn_indices_and_distances(
                                                                      train = single_point_train_grid,
                                                                      test = single_test_point,
                                                                      k = 3, 
                                                                      p = 2)
  
  # check that the original knn method, with normalize off, gives neighbors whose euclidean distance is same as 
  # the distances returned by the new knn_single_point_nn_indices_and_distances 
  all.equal(
    indices_and_distances_from_single_random_point_tested_on_grid$distances - 
      single_point_train_grid[indices_from_single_random_point_tested_on_grid,] %>%
      mutate(testx = single_test_point[1],
             testy = single_test_point[2]) %>%
      mutate(dist = sqrt((x - testx)^2 + (y - testy)^2)) %>%
      .$dist, 
    c(0,0,0))
  
  
  
  #test character "k"
  knn_single_point(train = single_point_train_grid,
                   test = c(runif(2,-3, 3)),
                   k = "character", 
                   normalize = T, 
                   p = 2)
  
  #test non-integer k
  knn_single_point(train = single_point_train_grid,
                   test = c(runif(2,-3, 3)),
                   k = 3.14, 
                   normalize = T, 
                   p = 2)
  
  #test k less than 1
  knn_single_point(train = single_point_train_grid,
                   test = c(runif(2,-3, 3)),
                   k = 0, 
                   normalize = T, 
                   p = 2)
  
  #test k greater than number of points in train
  #test non-integer k
  all.equal(knn_single_point(train = single_point_train_grid,
                   test = c(runif(2,-3, 3)),
                   k = nrow(single_point_train_grid) + 1, 
                   normalize = T, 
                   p = 2),
            1:nrow(single_point_train_grid) #should return all indices if k is greater than number of rows
            )
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

# applies knn algorithm to a test data set that can contain multiple points
# knn is applied to each point in test, against all the points in train
# normalize determines whether or not train and test are mean centered on 0/ and SD scaled to 1 when applying knn
# p determines the type of Lp metric
# returns original test data with nn and ypred columns appended.
knn <- function(train, test, target, categorical = T, k = 3, normalize = T, p = 2) {
 
  if(sum(!(names(test) %in% names(train))) > 0) {
    stop("The test set needs to have column names present in the training data!")
  }
  
  #convert a single numeric vector into a tibble for the code to work
  if(typeof(test) != "list") {
    tibble_storage_for_test <- train %>% slice(1)
    for(i in 1:length(test)) {
      tibble_storage_for_test[1,i] <- test[i]
    }
    test <- tibble_storage_for_test
  }
  
  
  nn <- list() #initialize empty for loop
  ypred <- c() #initialize empty for loop
  for(i in 1:nrow(test)) {
    nn <- c(nn, list(knn_single_point(train = train[names(test)], test = test %>% slice(i), k = k, normalize = normalize, p = p)))
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
        nn = nn,
        ypred = ypred
      )
    )
}


##########DEMO OF KNN ON 2D DATA GENERATED AS FOLLOWS:#############
# Set classes for points in a grid such that the class is based on distance to one of 
# two centers. The methods vary in randomness and how quickly the influence of the centers
# decays with distance.
{
(center1 <- sample(xy$x, 2, replace = T)*1.0)
(center2 <- sample(xy$x, 2, replace = T)*1.0)

centers <- tibble(
  x = c(center1[1], center2[1]),
  y = c(center1[2], center2[2])
  )

symbol1 <- "-"
symbol2 <- "+"

color1 <- "black"
color2 <- "black"

(xy2 <- xy %>%
  rowwise() %>% 
  mutate(dist_from_center1 = Lp(c(x,y), center1),
         dist_from_center2 = Lp(c(x,y), center2)) %>%
  mutate(definite_class_with_random_ties = case_when(dist_from_center1 < dist_from_center2 ~ symbol1,
                                                      dist_from_center1 > dist_from_center2 ~ symbol2,
                                                      TRUE ~ sample(c(symbol1, symbol2), 1))) %>%  #tie assigns randomly
  rowwise() %>%
  mutate(
    random_class_based_on_guassian = sample(x = c(symbol1, symbol2), size = 1,
      prob = case_when(dist_from_center1 == 0 ~ c(1,0), #guarantee class of center1
                       dist_from_center2 == 0 ~ c(0,1), #guarantee class of center2
            TRUE ~ c(exp(-1.0*(dist_from_center1)^2), exp(-1.0*(dist_from_center2)^2))/sum(exp(-1.0*(dist_from_center1)^2), exp(-1.0*(dist_from_center2)^2)))),
    random_class_based_on_dist = sample(x = c(symbol1, symbol2), size = 1,
        prob = case_when(dist_from_center1 == 0 ~ c(1,0), #guarantee class of center1 & handle div by 0 when on center1
                         dist_from_center2 == 0 ~ c(0,1), #guarantee class of center2 & handle div by 0 when on center2
                         TRUE ~ c(1/abs(dist_from_center1), 1/abs(dist_from_center2))/sum(1/abs(dist_from_center1), 1/abs(dist_from_center2)))
      )
  )
)


size_label_text <- 8
size_circle <- 6
library(scales)
test <- c(runif(2,-3, 3))
if(typeof(test) != "list") {
  tibble_storage_for_test <- train %>% slice(1)
  for(i in 1:length(test)) {
    tibble_storage_for_test[1,i] <- test[i]
  }
}
test <- tibble_storage_for_test

test_knn_based_on_gaussian <- knn(train=train, test = test, target = "random_class_based_on_guassian", k = 5, p = 2)
knn_on_test_predictions <- knn_on_test$ypred

knn_on_test_colors <- case_when(knn_on_test_predictions == symbol1 ~ hue_pal()(2)[1],
                                TRUE ~ hue_pal()(2)[2])
predictions_tibble <- tibble(
  x = test$x,
  y = test$y,
  predictions = knn_on_test_predictions,
  colors = knn_on_test_colors
)

ggplot(xy2, aes(x = x, y = y)) + 
  geom_point(data = centers %>% slice_head(n = 1), aes(x, y), color = hue_pal()(2)[1], size = size_circle) + 
  geom_point(data = centers %>% slice_tail(n = 1), aes(x, y), color = hue_pal()(2)[2], size = size_circle) +
  geom_text(aes(label = random_class_based_on_guassian, color = random_class_based_on_guassian), size = size_label_text) +
  #geom_point(data = test, aes(x, y), size = size_circle, shape = 21, color = "purple") +
  geom_point(data = xy[knn_single_point(train=xy, test = test, k = 5, p = 2),], aes(x = x, y = y), shape = 21, color = "blue", size = size_circle) +
  geom_point(data = knn(train=train, test = test, target = "random_class_based_on_guassian", k = 5, p = 2), aes(x = x, y = y), color = "purple", size = size_circle) +
  geom_text(data = knn(train=train, test = test, target = "random_class_based_on_guassian", k = 5, p = 2), aes(x = x, y = y, label = ypred, color = ypred), size = size_label_text*2, hjust = .5, vjust = 0.35) +
  labs(title = "Labels determined by distance from center (probabilistic, exponential decay)") +
  theme(legend.position="bottom", legend.text = element_text(size = 12 + size_label_text)) +
  theme(legend.position = "none") +
  guides(color=guide_legend(title="Class Color"))

ggplot(xy2, aes(x = x, y = y)) + 
  geom_point(data = centers %>% slice_head(n = 1), aes(x, y), color = hue_pal()(2)[1], size = size_circle) + 
  geom_point(data = centers %>% slice_tail(n = 1), aes(x, y), color = hue_pal()(2)[2], size = size_circle) +
  geom_text(aes(label = random_class_based_on_dist, color = random_class_based_on_dist), size = size_label_text) +
  geom_point(data = xy[knn_single_point(train=xy, test = test, k = 5, p = 2),], aes(x = x, y = y), shape = 21, color = "blue", size = size_circle) +
  geom_point(data = knn(train=train, test = test, target = "random_class_based_on_dist", k = 5, p = 2), aes(x = x, y = y), color = "purple", size = size_circle) +
  geom_text(data = knn(train=train, test = test, target = "random_class_based_on_dist", k = 5, p = 2), aes(x = x, y = y, label = ypred, color = ypred), size = size_label_text*2, hjust = .5, vjust = 0.35) +
  labs(title = "Labels determined by distance from center (probabilistic, linear decay)") +
  theme(legend.position="bottom", legend.text = element_text(size = 12 + size_label_text)) +
  theme(legend.position = "none") +
  guides(color=guide_legend(title="Class Color"))

ggplot(xy2, aes(x = x, y = y)) + 
  geom_point(data = centers %>% slice_head(n = 1), aes(x, y), color = hue_pal()(2)[1], size = size_circle) + 
  geom_point(data = centers %>% slice_tail(n = 1), aes(x, y), color = hue_pal()(2)[2], size = size_circle) +
  geom_text(aes(label = definite_class_with_random_ties, color = definite_class_with_random_ties), size = size_label_text) +
  geom_point(data = xy[knn_single_point(train=xy, test = test, k = 5, p = 2),], aes(x = x, y = y), shape = 21, color = "blue", size = size_circle) +
    geom_point(data = knn(train=train, test = test, target = "definite_class_with_random_ties", k = 5, p = 2), aes(x = x, y = y), color = "purple", size = size_circle) +
  geom_text(data = knn(train=train, test = test, target = "definite_class_with_random_ties", k = 5, p = 2), aes(x = x, y = y, label = ypred, color = ypred), size = size_label_text*2, hjust = .5, vjust = 0.35) +
  labs(title = "Labels determined by distance from center (determinstic except ties)") +
  theme(legend.position="bottom", legend.text = element_text(size = 12 + size_label_text)) +
  theme(legend.position = "none") +
  guides(color=guide_legend(title="Class Color"))
  #      text = element_text(size = size_label_text))

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
get_knn_decision_boundary <- function(num_points = 50, min_x = -5.0, max_x = 5.0, min_y = min_x, max_y = max_x, train, target, categorical = T, k = 3, normalize = T, p = 2) {
  test_grid <- get_grid(num_points, min_x, max_x, min_x, max_x)
  
  # appends nn and ypred columns for nearest neighbors and predicted y-values to the training data which in this case is the 
  # num_points by num_points grid covering the rectangular region defined by [min_x,max_x] X [min_y,max_y] 
  knn(train=train, test = test_grid, target = target, categorical = categorical, k = k, normalize = normalize, p = p)
}

# warning this takes a long time!
decision_boundary <- get_knn_decision_boundary(num_points = 50, min_x = -5.0, max_x = 5.0 , min_y = -5.0, max_y = 5.0, 
                                               train=train, target = "definite_class_with_random_ties", k = 5) 


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


graph_knn_decision_boundary(decision_boundary)
