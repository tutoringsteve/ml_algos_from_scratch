source("kd_trees.R")


########### utility functions #######################

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

########### test data creation #######################

cellx <- runif(100)
celly <- runif(100)

kd_test_grid_2d <- bind_rows(
  expand.grid(cellx - 2, celly*2) %>% mutate(cell = 1),
  expand.grid(cellx - 1, -1*celly*2)  %>% mutate(cell = 2),
  expand.grid(cellx + 1, (0.5-1*(celly*2.5)))  %>% mutate(cell = 3),
  expand.grid(cellx, 0.5+3*celly/2) %>% mutate(cell = 4)
) %>%
  rename(x = Var1, y = Var2)


#plot the grid
(gg_2d_test_grid <- ggplot(kd_test_grid_2d) +
    geom_point(aes(x, y, color = factor(cell))) +
  theme_minimal())

########### KD-Tree creation #######################

max_depth <- 3
root <- BuildKDTree(kd_test_grid_2d %>% 
                              mutate(cell = as.numeric(cell)),
                            "cell", 
                            max_depth = max_depth, 
                            min_data_in_leaf = 4)
root$print_short()



########### cut line result graphing #######################


cut_points2 <- pre_order_traversal_2d_tree_get_cut_points(root, xrange = integer_range(kd_test_grid_2d, "x"), yrange = integer_range(kd_test_grid_2d, "y"))

# nodes order isn't same as node number so manually renumber
cut_points2$line <-c(1,2,4,5,3,6,7)
cut_points2 <- cut_points2 %>%
  # for graphing the cuts centered on but away from the cuts
  mutate(xmove = case_when(var == "y" ~ -0.15, TRUE ~ 0)) %>%
  mutate(ymove = case_when(var == "x" ~ -0.15, TRUE ~ 0))


p <- gg_2d_test_grid
# loop through the cuts to the build the plot
for(i in 1:nrow(cut_points2)) {
  point <- cut_points2[i,]
  
  p <- p +
    geom_segment(data = point, aes(x = x1, xend = x2, y = y1, yend = y2), color = "black", size = 3) +
    geom_text(data = point, aes(x = mean(x1,x2) + xmove, y = mean(y1,y2) + ymove, 
                                label = paste0("cut #", line, ":\n", var," = ",round(cut,2))), color = "red", size = 4)
}

# finally display
p
