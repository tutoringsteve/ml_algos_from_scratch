source("knn_from_scratch.R")

#45-45-90 triangle check
x<- c(0,1)
y<- c(1,0)
if(all.equal(sqrt(2),Lp(x,y, 2))) print("Lp passed 45-45-90 triangle check")

#30-60-90 triangle check
x <- c(0,0)
y <- c(1/2, sqrt(3)/2)
if(all.equal(1,Lp(x,y, 2))) print("Lp passed 30-60-90 triangle check")




# double/real valued, but lies on integer grid points from -5 to 5 in x and y directions
# will be used throughout testing code
train_grid <- expand_grid(
  x= (-5:5)*1.0,
  y= (-5:5)*1.0
)

# alias of train_grid
xy <- train_grid


###########tests of knn_single_point###############

#randomly assigned real number betwixt the grid points
single_test_point <- c(runif(2,-3, 3))
indices_from_single_random_point_tested_on_grid <- knn_single_point(
  train = train_grid,
  test = single_test_point,
  k = 3, 
  normalize = T, 
  p = 2)

non_normalized_indices_from_single_random_point_tested_on_grid <- knn_single_point(
  train = train_grid,
  test = single_test_point,
  k = 3, 
  normalize = F, 
  p = 2)

indices_and_distances_from_single_random_point_tested_on_grid <- knn_single_point_nn_indices_and_distances(
  train = train_grid,
  test = single_test_point,
  k = 3, 
  p = 2)

# check that the original knn method, with normalize off, gives neighbors whose euclidean distance is same as 
# the distances returned by the new knn_single_point_nn_indices_and_distances 
all.equal(
  indices_and_distances_from_single_random_point_tested_on_grid$distances - 
    train_grid[indices_from_single_random_point_tested_on_grid,] %>%
    mutate(testx = single_test_point[1],
           testy = single_test_point[2]) %>%
    mutate(dist = sqrt((x - testx)^2 + (y - testy)^2)) %>%
    .$dist, 
  c(0,0,0))



#test character "k" : should produce error.
knn_single_point(train = train_grid,
                 test = c(runif(2,-3, 3)),
                 k = "character", 
                 normalize = T, 
                 p = 2)

#test non-integer k : should produce error.
knn_single_point(train = train_grid,
                 test = c(runif(2,-3, 3)),
                 k = 3.14, 
                 normalize = T, 
                 p = 2)

#test k less than 1 : should produce error.
knn_single_point(train = train_grid,
                 test = c(runif(2,-3, 3)),
                 k = 0, 
                 normalize = T, 
                 p = 2)

#test k greater than number of points in train
#test non-integer k
all.equal(knn_single_point(train = train_grid,
                           test = c(runif(2,-3, 3)),
                           k = nrow(train_grid) + 1, 
                           normalize = T, 
                           p = 2),
          1:nrow(train_grid) #should return all indices if k is greater than number of rows
)

##########DEMO OF KNN ON 2D DATA GENERATED AS FOLLOWS:#############
# Set classes for points in a grid such that the class is based on distance to one of 
# two centers. The methods vary in randomness and how quickly the influence of the centers
# decays with distance.
{
  (center1 <- sample(train_grid$x, 2, replace = T)*1.0)
  (center2 <- sample(train_grid$x, 2, replace = T)*1.0)
  
  centers <- tibble(
    x = c(center1[1], center2[1]),
    y = c(center1[2], center2[2])
  )
  
  symbol1 <- "-"
  symbol2 <- "+"
  
  color1 <- "black"
  color2 <- "black"
  
  (train <- train_grid %>%
      rowwise() %>% 
      mutate(dist_from_center1 = Lp(c(x,y), center1),
             dist_from_center2 = Lp(c(x,y), center2)) %>%
      ungroup() %>% # Removes rowwise status
      mutate(
        definite_class_with_random_ties = case_when(dist_from_center1 < dist_from_center2 ~ symbol1,
                                                    dist_from_center1 > dist_from_center2 ~ symbol2,
                                                    #tie assigns randomly
                                                    TRUE ~ sample(c(symbol1, symbol2), 1)) 
      ) %>%
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
      ) %>%
      #removes the rowwise status
      ungroup()
  )
  
  features <- names(train_grid)
  
  test <- c(runif(2,-3, 3))
  test <- if_necessary_convert_test_from_vector_to_tibble(train, features, test)
  
  test_knn_based_on_gaussian <- knn(train=train, features = features, test = test, target = "random_class_based_on_guassian", k = 5, p = 2)
  knn_on_test_predictions <- test_knn_based_on_gaussian$ypred
  
  knn_on_test_colors <- case_when(knn_on_test_predictions == symbol1 ~ hue_pal()(2)[1],
                                  TRUE ~ hue_pal()(2)[2])
  predictions_tibble <- tibble(
    x = test$x,
    y = test$y,
    predictions = knn_on_test_predictions,
    colors = knn_on_test_colors
  )
  
  size_label_text <- 8
  size_circle <- 6
  
  ggplot(train, aes(x = x, y = y)) + 
    geom_point(data = centers %>% slice_head(n = 1), aes(x, y), color = hue_pal()(2)[1], size = size_circle) + 
    geom_point(data = centers %>% slice_tail(n = 1), aes(x, y), color = hue_pal()(2)[2], size = size_circle) +
    geom_text(aes(label = random_class_based_on_guassian, color = random_class_based_on_guassian), size = size_label_text) +
    #geom_point(data = test, aes(x, y), size = size_circle, shape = 21, color = "purple") +
    geom_point(data = xy[knn_single_point(train=xy, test = test, k = 5, p = 2),], 
               aes(x = x, y = y), shape = 21, color = "blue", size = size_circle) +
    geom_point(data = knn(train=train, features = features, test = test, target = "random_class_based_on_guassian", k = 5, p = 2), 
               aes(x = x, y = y), color = "purple", size = size_circle) +
    geom_text(data = knn(train=train, features = features, test = test, target = "random_class_based_on_guassian", k = 5, p = 2), 
              aes(x = x, y = y, label = ypred, color = ypred), size = size_label_text*2, hjust = .5, vjust = 0.35) +
    labs(title = "Labels determined by distance from center (probabilistic, exponential decay)") +
    theme(legend.position="bottom", legend.text = element_text(size = 12 + size_label_text)) +
    theme(legend.position = "none") +
    guides(color=guide_legend(title="Class Color"))
  
  
  ggplot(train, aes(x = x, y = y)) + 
    geom_point(data = centers %>% slice_head(n = 1), aes(x, y), color = hue_pal()(2)[1], size = size_circle) + 
    geom_point(data = centers %>% slice_tail(n = 1), aes(x, y), color = hue_pal()(2)[2], size = size_circle) +
    geom_text(aes(label = random_class_based_on_dist, color = random_class_based_on_dist), size = size_label_text) +
    geom_point(data = xy[knn_single_point(train=xy, test = test, k = 5, p = 2),], 
               aes(x = x, y = y), shape = 21, color = "blue", size = size_circle) +
    geom_point(data = knn(train=train, features = features, test = test, target = "random_class_based_on_dist", k = 5, p = 2), 
               aes(x = x, y = y), color = "purple", size = size_circle) +
    geom_text(data = knn(train=train, features = features, test = test, target = "random_class_based_on_dist", k = 5, p = 2), 
              aes(x = x, y = y, label = ypred, color = ypred), size = size_label_text*2, hjust = .5, vjust = 0.35) +
    labs(title = "Labels determined by distance from center (probabilistic, linear decay)") +
    theme(legend.position="bottom", legend.text = element_text(size = 12 + size_label_text)) +
    theme(legend.position = "none") +
    guides(color=guide_legend(title="Class Color"))
  
  ggplot(train, aes(x = x, y = y)) + 
    geom_point(data = centers %>% slice_head(n = 1), aes(x, y), color = hue_pal()(2)[1], size = size_circle) + 
    geom_point(data = centers %>% slice_tail(n = 1), aes(x, y), color = hue_pal()(2)[2], size = size_circle) +
    geom_text(aes(label = definite_class_with_random_ties, color = definite_class_with_random_ties), size = size_label_text) +
    geom_point(data = xy[knn_single_point(train=xy, test = test, k = 5, p = 2),], 
               aes(x = x, y = y), shape = 21, color = "blue", size = size_circle) +
    geom_point(data = knn(train=train, features = features, test = test, target = "definite_class_with_random_ties", k = 5, p = 2), 
               aes(x = x, y = y), color = "purple", size = size_circle) +
    geom_text(data = knn(train=train, features = features, test = test, target = "definite_class_with_random_ties", k = 5, p = 2), 
              aes(x = x, y = y, label = ypred, color = ypred), size = size_label_text*2, hjust = .5, vjust = 0.35) +
    labs(title = "Labels determined by distance from center (determinstic except ties)") +
    theme(legend.position="bottom", legend.text = element_text(size = 12 + size_label_text)) +
    theme(legend.position = "none") +
    guides(color=guide_legend(title="Class Color"))
  #      text = element_text(size = size_label_text))
  
}

########### test decision boundary code
{
  # warning this takes a long time!
  decision_boundary <- get_knn_decision_boundary(num_points = 25, min_x = -5.0, max_x = 5.0 , min_y = -5.0, max_y = 5.0, 
                                                 train=train, features=features, target = "definite_class_with_random_ties", k = 5) 
  
  graph_knn_decision_boundary(decision_boundary)
}
