source("kd_trees.R")

max_depth <- 1
root <- BuildKDTree(ames,"Sale_Price", max_depth = max_depth, min_data_in_leaf = 4)
rootSimple <- BuildKDTreeSimple(ames,"Sale_Price", max_depth = max_depth, min_data_in_leaf = 4)
root$print_short()
rootSimple$print_short()


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
root_node_2d_test_simple <- BuildKDTreeSimple(kd_test_grid_2d %>% mutate(cell = as.numeric(cell)),"cell", max_depth = max_depth, min_data_in_leaf = 4)
root_node_2d_test$print_short()
root_node_2d_test_simple$print_short()


leaf_data_where_test_is <- kd_test_grid_2d
test <- leaf_data_where_test_is[sample.int(nrow(leaf_data_where_test_is),1),]

direction <- c()
cuts <- list()
node <- root_node_2d_test_simple

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


all.equal(leaf_data_where_test_is, get_leaf_data_at(root_node_2d_test_simple, kd_test_grid_2d, test))
get_cuts_to(root_node_2d_test_simple, kd_test_grid_2d, test)
get_direction_to(root_node_2d_test_simple, kd_test_grid_2d, test)
root_node_2d_test_simple$print_short()

gg_2d_test_grid +
  geom_point(data = leaf_data_where_test_is, aes(x = x, y = y), color = "red", size= 3) +
  geom_point(data = test, aes(x = x, y = y), color = "black", size= 4, shape = 21) +
  geom_point(data = test, aes(x = x, y = y), color = "black", size= 1)
  
#dat <- circleFun(c(test$x,test$y),,npoints = 100)


#TODO 1 graph all the splits.
#TODO 2 find the distance from the test point to each of the splits
#TODO 3 find the knn with the single point in the cell.

pre_order_traversal_2d_tree_get_cut_points <- function(root, xrange = c(0,1), yrange = c(0,1), 
                                                       points = tibble(x1 = c(), y1 = c(), x2 = c(), y2 = c(), var = c(), cut = c())) {
  if(!is.null(root)) {
    #current node
    
    var <- root$cut_variable
    cut <- root$cut_value
    
    if(var == "x") {
      points <- points %>% 
        bind_rows(tibble(x1 = cut, y1 = yrange[1], 
                         x2 = cut, y2 = yrange[2],
                         var = var, cut = cut))
      
      xrange_left <- c(xrange[1], cut)
      xrange_right <- c(cut, xrange[2])
      
      yrange_left <- yrange
      yrange_right <- yrange
    } else {
      points <- points %>% 
        bind_rows(tibble(x1 = xrange[1], y1 = cut, 
                         x2 = xrange[2], y2 = cut,
                         var = var, cut = cut))
      
      yrange_left <- c(yrange[1], cut)
      yrange_right <- c(cut, yrange[2])
      
      xrange_left <- xrange
      xrange_right <- xrange
    }
    
    points <- pre_order_traversal_2d_tree_get_cut_points(root$left, xrange_left, yrange_left, points)
    points <- pre_order_traversal_2d_tree_get_cut_points(root$right, xrange_right, yrange_right, points)
    
  }
  return(points) 
}


max_depth <- 2

root <- BuildKDTreeSimple(kd_test_grid_2d %>% 
                      mutate(cell = as.numeric(cell)),
                    "cell", 
                    max_depth = max_depth, 
                    min_data_in_leaf = 4)
root$print_short()
cut_points <- pre_order_traversal_2d_tree_get_cut_points(root)

cut_points$line <- factor(1:nrow(cut_points))

p <- ggplot()

for(i in 1:nrow(cut_points)) {
  point <- cut_points[i,]
  p <- p +
    geom_segment(data = point, aes(x = x1, xend = x2, y = y1, yend = y2, color = line, size = line)) +
    geom_text(data = point, aes(x = mean(x1,x2), y = mean(y1,y2), label = paste0(var," = ",round(cut,2))), color = "black")
}

p


cellx <- runif(100)
celly <- runif(100)

test_grid2 <- 
expand.grid(cellx - 2, celly*2)
expand.grid(cellx - 1, -1*celly*2)
expand.grid(cellx, celly*2)
expand.grid(cellx + 1, -1*celly*2)


kd_test_grid_2d_2 <- bind_rows(
  expand.grid(cellx - 2, celly*2) %>% mutate(cell = 1),
  expand.grid(cellx - 1, -1*celly*2)  %>% mutate(cell = 2),
  expand.grid(cellx + 1, (0.5-1*(celly*2.5)))  %>% mutate(cell = 3),
  expand.grid(cellx, 0.5+3*celly/2) %>% mutate(cell = 4)
) %>%
  rename(x = Var1, y = Var2)


#plot the grid
(gg_2d_test_grid_2 <- ggplot(kd_test_grid_2d_2) +
    geom_point(aes(x, y, color = factor(cell))) +
    #geom_segment(aes(x = 0, xend = xcuts[2], y = ycuts[1], yend = ycuts[1]), color = "black") +
    #geom_segment(aes(x = xcuts[1], xend = xcuts[1], y = 0, yend = ycuts[1]), color = "black") +
    #geom_segment(aes(x = xcuts[2], xend = xcuts[2], y = 0, yend = 1), color = "black") +
    #geom_segment(aes(x = xcuts[2], xend = 1, y = ycuts[2], yend = ycuts[2]), color = "black") +
    #geom_segment(aes(x = xcuts[3], xend = xcuts[3], y = 0, yend = ycuts[2]), color = "black") +
    #geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
    #geom_segment(aes(x = 0, xend = 1, y = 1, yend = 1), color = "black") +
    #geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1), color = "black") +
    #geom_segment(aes(x = 1, xend = 1, y = 0, yend = 1), color = "black") +
    #geom_point(data = test_best, aes(x, y), color = "black") +
    #geom_label(data = test_best, aes(x, y, label = "best"), vjust = -0.25) +
    #geom_point(data = test_worst, aes(x, y), color = "black") +
    #geom_label(data = test_worst, aes(x, y, label = "worst"), vjust = -0.25) +
    theme_minimal())


max_depth <- 2

root_2 <- BuildKDTreeSimple(kd_test_grid_2d_2 %>% 
                            mutate(cell = as.numeric(cell)),
                          "cell", 
                          max_depth = max_depth, 
                          min_data_in_leaf = 4)
root_2$print_short()



integer_range <- function(data, var) {
  if(var %in% names(data)) {
    if(!is.numeric(data[[var]])) {
      stop(paste0("Variable ", var, " in data is non-numeric."))
    } else {
      return(c(floor(min(data[[var]])), ceiling(max(data[[var]]))))
    }
  } else {
    stop(paste0("Variable ", var, " not in data!"))
  }
}


cut_points2 <- pre_order_traversal_2d_tree_get_cut_points(root_2, xrange = integer_range(kd_test_grid_2d_2, "x"), yrange = integer_range(kd_test_grid_2d_2, "y"))

cut_points2$line <-c(1,2,4,5,3,6,7)
cut_points2 <- cut_points2 %>%
  mutate(xmove = case_when(var == "y" ~ -0.15, TRUE ~ 0)) %>%
  mutate(ymove = case_when(var == "x" ~ -0.15, TRUE ~ 0))

library(viridis)
library(RColorBrewer)
display.brewer.pal(name = "RdPu")
p2 <- gg_2d_test_grid_2
#p2 <- ggplot()
for(i in 1:nrow(cut_points2)) {
  point <- cut_points2[i,]
  
  p2 <- p2 +
    geom_segment(data = point, aes(x = x1, xend = x2, y = y1, yend = y2), color = "black", size = 3) +
    geom_text(data = point, aes(x = mean(x1,x2) + xmove, y = mean(y1,y2) + ymove, 
                                label = paste0("cut #", line, ":\n", var," = ",round(cut,2))), color = "red", size = 4)
}

p2

