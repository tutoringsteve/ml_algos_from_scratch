library(tidyverse)
library(tidymodels)
tidymodels_prefer()
if (!requireNamespace("R6", quietly = TRUE)) install.packages("R6")


KDNode <- R6Class(
  classname = "KDNode",
  public = list(
    data = NULL,
    cut_variable = NULL,
    cut_value = NULL,
    left = NULL,
    right = NULL,
    
    # Constructor
    initialize = function(data = NULL, cut_variable = NULL, cut_value = NULL, left = NULL, right = NULL) {
      self$data <- data
      self$cut_variable <- cut_variable
      self$cut_value <- cut_value
      self$left <- left
      self$right <- right
    },
    
    # Special print function
    print = function() {
      cat("KDNode:\n")
      cat("Data:\n")
      cat("dimension: ", dim(self$data), "\n")
      cat("names: ", names(self$data), "\n")
      cat(paste0(c("Cut Variable: ", self$cut_variable), collapse = ""), "\n")
      cat(paste0(c("Cut Value: ", self$cut_value), collapse = ""), "\n")
      if (!is.null(self$left)) {
        cat("Left Child: Exists\n")
        print(self$left)
      } else {
        cat("Left Child: NULL\n")
      }
      if (!is.null(self$right)) {
        cat("Right Child: Exists\n")
        print(self$right)
      } else {
        cat("Right Child: NULL\n")
      }
    },
    
    print_short = function(dir = NULL, layer = "root", depth = 0, tabs = 0) {
      if(!is.null(dir)) {
        cat(rep("\t",depth),paste(dir,paste0(layer,":")), nrow(self$data), "\n")
      } else {
        cat(rep("\t",depth),paste0(layer,":"), nrow(self$data), "\n")
      }
      if (!is.null(self$left)) {
        if(depth == 0) layer <- "child"
        if(depth == 1) layer <- "grand child"
        if(depth > 1) layer <- paste0(paste0(rep("great", depth - 1), collapse=" ")," grand child")
        self$left$print_short(dir = "left", layer = layer, depth = depth + 1)
      } 
      if (!is.null(self$right)) {
        if(depth == 0) layer <- "child"
        if(depth == 1) layer <- "grand child"
        if(depth > 1) layer <- paste0(paste0(rep("great", depth - 1), collapse=" ")," grand child")
        self$right$print_short(dir = "right", layer, depth = depth + 1)
      }
    }
  )
)

BuildKDTree <- function(data, target, last_cut_indices = c(), target_is_numeric = is.numeric(data[[target]]), depth = 0, max_depth = NULL, min_data_in_leaf = NULL) {
  
  if(!is.null(max_depth)) {
    if(depth > max_depth) {
      return(NULL)    
    }
  }
  
  if(nrow(data) <= 1) return(NULL)
  if(!is.null(min_data_in_leaf)) {
    if(nrow(data) <= min_data_in_leaf) {
      return(NULL)    
    }
  }
  
  #checks if a column for row numbers can be added
  if(!"row_num" %in% names(data)) {
    data$row_num <- 1:nrow(data)
  } else {
    if(all.equal(data$row_num, 1:nrow(data))) {
      stop("data has row_num column already and it's not equal to 1:nrow(data)!\nInvestigate!")
    }
  }
  
  #checks to see if parameter is a lie
  #TODO refactor this out!
  if(target_is_numeric && !is.numeric(data[[target]])) {
    stop("target isn't actually numeric!")
  }
  
  medians <- data[,-which(names(data) %in% c(target, "row_num"))] %>%
                summarise_if(is.numeric, median, na.rm = T)
  
  index_range <- 1:ncol(medians)
  #if last_cut_indices is non-empty filter them out from the candidate indices
  if(length(last_cut_indices)){
    index_range <- index_range[!index_range %in% last_cut_indices]
  }
  
  #obtains the indices in left/right nodes for data and values of target for each possible split
  target_in_each <- lapply(index_range, function(i) 
    list(var = names(medians[i]), 
         L = data %>% filter(data[names(medians)[i]] < medians[[i]]) %>% select(row_num, all_of(target)),
         R = data %>% filter(data[names(medians)[i]] >= medians[[i]]) %>% select(row_num, all_of(target))
    )
  )
  
  size_of_L <- sapply(target_in_each, function(x) nrow(x$L))
  (size_of_R <- sapply(target_in_each, function(x) nrow(x$R)))
  ideal_target_in_each <- target_in_each[(size_of_L > nrow(data)/4) & (size_of_R > nrow(data)/4)]
  #first variable should minimize the variance in sum of each node if target is numeric
  # if target is categorical it should cut such that 
  # dist(distribution(target in L), uniform_distribution(target)) is maximized!
  # we want to choose leaf nodes to avoid each leaf having equal numbers of the possible target classes
  if(target_is_numeric) {
    cut_var_index <- which.min(sapply(ideal_target_in_each, function(x) {
            left_var <- 0
            right_var <- 0
            if(nrow(x$L) > 0) left_var <- var(x$L[target])
            if(nrow(x$R) > 0) right_var <- var(x$R[target])
            
            return(left_var + right_var)
          }
        )
      )
    
    cut_var_index <- which(sapply(target_in_each, function(x) x$var) == ideal_target_in_each[[cut_var_index]]$var)
  } else{
    #TODO handle the categorical case!
    stop(paste0(c("non-numeric target not coded yet.",
                  "need to measure distance from uniform distribution in each!",
                  "review notes :)"), collapse = "\n"))
    #cut_var_index <- which.max(sapply(target_in_each, function(x) sum dist L from uniform , dist R from uniform
  }
  
  # add the cut variable index to the list of already cut indices
  # if the new list indexes over all of medians, then empty
  # last_cut_indices to start recycling variables again
  last_cut_indices <- c(last_cut_indices, cut_var_index)
  if(length(last_cut_indices) == ncol(medians)) last_cut_indices <- c()
  
  
  #build the tree depth first from left side through to right
  node <- KDNode$new(data = data, cut_variable = names(medians)[[cut_var_index]], cut_value = medians[[cut_var_index]])
  
  left_row_nums <- target_in_each[[cut_var_index]]$L$row_num
  right_row_nums <- target_in_each[[cut_var_index]]$R$row_num
  
  node$left <- BuildKDTree(data = data[left_row_nums,] %>% select(-row_num), target = target, target_is_numeric = target_is_numeric, depth = depth + 1, max_depth = max_depth, min_data_in_leaf = min_data_in_leaf)
  node$right <- BuildKDTree(data = data[right_row_nums,] %>% select(-row_num), target = target, target_is_numeric = target_is_numeric, depth = depth + 1, max_depth = max_depth, min_data_in_leaf = min_data_in_leaf)
  
  return(node)
}

  BuildKDTree(ames,"Sale_Price", max_depth = 10, min_data_in_leaf = 4)$print_short()
