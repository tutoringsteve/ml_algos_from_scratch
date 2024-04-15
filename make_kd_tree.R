source("kd_trees.R")

max_depth <- 3
root_node <- BuildKDTree(ames,"Sale_Price", max_depth = max_depth, min_data_in_leaf = 4)
root_node$print_short()


#find the directions to the test node!
test <- ames[sample.int(nrow(ames),1),]
direction <- c()
cuts <- list()
node <- root_node
leaf_data_where_test_is <- ames
while(!is.null(node)) {
  if(test[[node$cut_variable]] < node$cut_value) {
    leaf_data_where_test_is <- leaf_data_where_test_is %>% filter(get(node$cut_variable) < node$cut_value)
    cuts <- c(cuts, list(cut_variable = node$cut_variable, cut_value = node$cut_value))
    direction <- c(direction,"L")
    node <- node$left
  } else {
    leaf_data_where_test_is <- leaf_data_where_test_is %>% filter(get(node$cut_variable) >= node$cut_value)
    cuts <- c(cuts, list(cut_variable = node$cut_variable, cut_value = node$cut_value))
    direction <- c(direction,"R")
    node <- node$right
  }
}

direction
cuts
leaf_data_where_test_is



#create the xy - coordinates of points such that a nearest neighbor search will have to look across most leaves in the KD-Tree

xcuts <- c(0.45,0.50, 0.55)
ycuts <- c(0.50,0.25)

dens = 5

x11 <- runif(3*dens, 0, xcuts[1])
x12 <- runif(1*dens, xcuts[1], xcuts[2])
y1 <- runif(2*dens, ycuts[1], 1)
y2 <- runif(2*dens, 0, ycuts[1])
x3 <- runif(2*dens, xcuts[2], 1)
y21 <- runif(2*dens, 0, ycuts[2])
y22 <- runif(2*dens, ycuts[2], 1)
x21 <- runif(1*dens, xcuts[2], xcuts[3])
x22 <- runif(3*dens, xcuts[3], 1)

kd_test_grid_2d <- bind_rows(
  expand.grid(c(x11, x12),y1),
  expand.grid(x11,y2),
  expand.grid(x12,y2),
  expand.grid(x21,y21),
  expand.grid(x22,y21),
  expand.grid(c(x21, x22),y22)
) %>%
  rename(x = Var1, y = Var2) %>%
  mutate(cell = factor(case_when(
    x <= xcuts[2] & y > ycuts[1] ~ 1,
    x <= xcuts[1] & y <= ycuts[1] ~ 2,
    x > xcuts[1] & x <= xcuts[2] & y <= ycuts[1] ~ 3,
    x > xcuts[2] & x <= xcuts[3]  & y <= ycuts[2] ~ 4,
    x > xcuts[3] & y <= ycuts[2] ~ 5,
    TRUE ~ 6
  )))

filter(kd_test_grid_2d, x > xcuts[1], x <= xcuts[2] , y <= ycuts[1])
filter(kd_test_grid_2d, x > xcuts[2], y <= ycuts[2])

#create best possible test point
(test_best <- tibble(x = xcuts[1]/2, y = mean(c(ycuts[1], 1))))
#create worst possible test point
(test_worst <- tibble(x = mean(c(xcuts[1],xcuts[2])), y = ycuts[2]))

#plot the grid
(gg_2d_test_grid <- ggplot(kd_test_grid_2d) +
  geom_point(aes(x, y, color = cell)) +
  geom_segment(aes(x = 0, xend = xcuts[2], y = ycuts[1], yend = ycuts[1]), color = "black") +
  geom_segment(aes(x = xcuts[1], xend = xcuts[1], y = 0, yend = ycuts[1]), color = "black") +
  geom_segment(aes(x = xcuts[2], xend = xcuts[2], y = 0, yend = 1), color = "black") +
  geom_segment(aes(x = xcuts[2], xend = 1, y = ycuts[2], yend = ycuts[2]), color = "black") +
  geom_segment(aes(x = xcuts[3], xend = xcuts[3], y = 0, yend = ycuts[2]), color = "black") +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
  geom_segment(aes(x = 0, xend = 1, y = 1, yend = 1), color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1), color = "black") +
  geom_segment(aes(x = 1, xend = 1, y = 0, yend = 1), color = "black") +
  geom_point(data = test_best, aes(x, y), color = "black") +
  geom_label(data = test_best, aes(x, y, label = "best"), vjust = -0.25) +
  geom_point(data = test_worst, aes(x, y), color = "black") +
  geom_label(data = test_worst, aes(x, y, label = "worst"), vjust = -0.25) +
  theme_minimal())


max_depth <- 3
root_node_2d_test <- BuildKDTree(kd_test_grid_2d %>% mutate(cell = as.numeric(cell)),"cell", max_depth = max_depth, min_data_in_leaf = 4)
root_node_2d_test$print_short()


leaf_data_where_test_is <- kd_test_grid_2d
test <- leaf_data_where_test_is[sample.int(nrow(leaf_data_where_test_is),1),]

direction <- c()
cuts <- list()
node <- root_node_2d_test

while(!is.null(node)) {
  if(test[[node$cut_variable]] < node$cut_value) {
    leaf_data_where_test_is <- leaf_data_where_test_is %>% filter(get(node$cut_variable) < node$cut_value)
    cuts <- c(cuts, list(cut_variable = node$cut_variable, cut_value = node$cut_value))
    direction <- c(direction,"L")
    node <- node$left
  } else {
    leaf_data_where_test_is <- leaf_data_where_test_is %>% filter(get(node$cut_variable) >= node$cut_value)
    cuts <- c(cuts, list(cut_variable = node$cut_variable, cut_value = node$cut_value))
    direction <- c(direction,"R")
    node <- node$right
  }
}


#all.equal(leaf_data_where_test_is, get_leaf_data(root_node_2d_test, kd_test_grid_2d, test))
get_cuts_to(root_node_2d_test, kd_test_grid_2d, test)
get_direction_to(root_node_2d_test, kd_test_grid_2d, test)


gg_2d_test_grid +
  geom_point(data = leaf_data_where_test_is, aes(x = x, y = y), color = "red", size= 3) +
  geom_point(data = test, aes(x = x, y = y), color = "black", size= 4, shape = 21) +
  geom_point(data = test, aes(x = x, y = y), color = "black", size= 1)
  
dat <- circleFun(c(test$x,test$y),,npoints = 100)
