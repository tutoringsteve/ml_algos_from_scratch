library(R6)
library(tidyverse)
library(tidymodels)
tidymodels_prefer()
if (!requireNamespace("R6", quietly = TRUE)) install.packages("R6")


KDNode <- R6Class(
  classname = "KDNode",
  public = list(
    rows = NULL,
    cut_variable = NULL,
    cut_value = NULL,
    left = NULL,
    right = NULL,
    
    # Constructor
    initialize = function(rows = NULL, cut_variable = NULL, cut_value = NULL, left = NULL, right = NULL) {
      self$rows <- rows
      self$cut_variable <- cut_variable
      self$cut_value <- cut_value
      self$left <- left
      self$right <- right
    },
    
    # Special print function
    print = function() {
      cat("KDNode:\n")
      cat("rows: ", self$rows, "\n")
      cat("num points: ", length(self$rows), "\n")
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
        cat(rep("\t",depth),paste(dir,paste0(layer," node #",":")), self$rows, "\n")
        cat(rep("\t",depth),paste(dir,paste0(layer," node #",":")), self$cut_variable, " ", self$cut_value, "\n")
      } else {
        cat(rep("\t",depth),paste0(layer," node #",":"), self$rows, "\n")
        cat(rep("\t",depth),paste0(layer," node #",":"), self$cut_variable, " ", self$cut_value, "\n")
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

# Builds a KD Tree by choosing among the numerical variables in data
# the variable that both keeps at least a quarter of data in both children
# and among those the one that minimizes the variance in each leaf.
# If no variable choice leads to both leaves having enough data, the criteria is
# softened by cutting it in half until one such candidate is found
# a list of already used variables is maintained and those variables are avoided
# until each variable is used once then the do not use list is rest.
BuildKDTree <- function(data, target, last_cut_indices = c(), target_is_numeric = is.numeric(data[[target]]), 
                        depth = 1, max_depth = NULL, min_data_in_leaf = NULL, verbose = FALSE) {
  
  if(verbose) print(paste0("depth: ", depth))
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
  
  if(verbose) cat("index range before filter: ", index_range, "\n")
  if(verbose) cat("last_cut_indeces: ", last_cut_indices, "\n")
  #if last_cut_indices is non-empty filter them out from the candidate indices
  if(length(last_cut_indices) > 0){
    index_range <- index_range[!index_range %in% last_cut_indices]
    if(verbose) cat("index range after filter: ", index_range, "\n")
  }
  
  #obtains the indices in left/right nodes for data and values of target for each possible split
  target_in_each <- lapply(index_range, function(i) 
    list(var = names(medians[i]), 
         L = data %>% filter(data[names(medians)[i]] < medians[[i]]) %>% select(row_num, all_of(target)),
         R = data %>% filter(data[names(medians)[i]] >= medians[[i]]) %>% select(row_num, all_of(target))
    )
  )
  
  sizes_of_L <- sapply(target_in_each, function(x) nrow(x$L))
  (sizes_of_R <- sapply(target_in_each, function(x) nrow(x$R)))
  # we want each sub branch to contain at least a quarter of the data.
  # so for each possible size of L and each possible size of R check condition, and then use 
  # vector boolean filtering to trim the list.
  
  
  ideal_fraction <- nrow(data)/4
  not_too_small <- (sizes_of_L > ideal_fraction) & (sizes_of_R > ideal_fraction)
  ideal_target_in_each <- target_in_each[not_too_small]
  
  # keep reducing ideal size until ideal_fraction rounds to zero or 
  # at least one of the targets can produce the ideal size.
  while(sum(not_too_small) == 0 & floor(ideal_fraction) > 0) {
    ideal_fraction <- ideal_fraction/2
    not_too_small <- (sizes_of_L > ideal_fraction) & (sizes_of_R > ideal_fraction)
    ideal_target_in_each <- target_in_each[not_too_small]
  }
  
  if(floor(ideal_fraction) == 0) warning(paste0("accepted a split with at most 1 element! At depth ", depth))
  
  # first variable should minimize the variance in sum of each node if target is numeric
  # if target is categorical it should cut such that 
  # dist(distribution(target in L), uniform_distribution(target)) is maximized!
  # we want to choose leaf nodes to avoid each leaf having equal numbers of the possible target classes
  if(target_is_numeric) {
    which_index_in_ideal_list <- which.min(sapply(ideal_target_in_each, function(x) {
            left_variance <- 0
            right_variance <- 0
            if(nrow(x$L) > 0) left_var <- var(x$L[target])
            if(nrow(x$R) > 0) right_var <- var(x$R[target])
            
            return(left_variance + right_variance)
          }
        )
      )
    
    # find which index of target_in_each had the variance of the preferred ideal target (which may have a reduced indices list
    # because of the filtering to ideal size)
    target_in_each_index <- which(sapply(target_in_each, function(x) x$var) == ideal_target_in_each[[which_index_in_ideal_list]]$var)
    cut_var_index <- index_range[[target_in_each_index]]
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
  
  if(verbose) print(paste0("depth = ", depth))
  if(verbose) print(paste0("cut_var_index = ", cut_var_index))
  if(verbose) print(paste0("last_cut_indices before = ", last_cut_indices))
  
  last_cut_indices <- c(last_cut_indices, cut_var_index)
  if(length(last_cut_indices) == ncol(medians)) {
      last_cut_indices <- c()
  }
  
  if(verbose) print(paste0("last_cut_indices after = ", last_cut_indices))

  #build the tree depth first from left side through to right
  #node <- KDNode$new(rows = nrow(data), cut_variable = names(medians)[[cut_var_index]], cut_value = medians[[cut_var_index]])
  
  left_row_nums <- target_in_each[[target_in_each_index]]$L$row_num
  right_row_nums <- target_in_each[[target_in_each_index]]$R$row_num
  
  left_child <- BuildKDTree(data = data[left_row_nums,] %>% select(-row_num), target = target, last_cut_indices = last_cut_indices, 
                            target_is_numeric = target_is_numeric, 
                            depth = depth + 1, max_depth = max_depth, min_data_in_leaf = min_data_in_leaf)

  right_child <- BuildKDTree(data = data[right_row_nums,] %>% select(-row_num), target = target, last_cut_indices = last_cut_indices, 
                             target_is_numeric = target_is_numeric, 
                             depth = depth + 1, max_depth = max_depth, min_data_in_leaf = min_data_in_leaf)

  return(KDNode$new(rows = nrow(data), cut_variable = names(medians)[[cut_var_index]], cut_value = medians[[cut_var_index]], left = left_child, right = right_child))
}



BuildKDTreeSimple <- function(data, target, target_is_numeric = is.numeric(data[[target]]), 
                        depth = 1, max_depth = NULL, min_data_in_leaf = NULL) {
  
  print(paste0("depth: ", depth))
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
  
  cut_var_index <- ((depth - 1) %% ncol(medians)) + 1
  cut_var <- names(medians)[[cut_var_index]]
  cut_val <- medians[[cut_var_index]]
  
  left_row_nums <- data %>%
    filter(get(cut_var) < cut_val) %>%
    .$row_num
  right_row_nums <- data %>%
    filter(get(cut_var) >= cut_val) %>%
    .$row_num
  
  left_child <- BuildKDTreeSimple(data = data[left_row_nums,] %>% select(-row_num), target = target, target_is_numeric = target_is_numeric, 
                            depth = depth + 1, max_depth = max_depth, min_data_in_leaf = min_data_in_leaf)
  
  right_child <- BuildKDTreeSimple(data = data[right_row_nums,] %>% select(-row_num), target = target, target_is_numeric = target_is_numeric, 
                             depth = depth + 1, max_depth = max_depth, min_data_in_leaf = min_data_in_leaf)
  
  return(KDNode$new(rows = nrow(data), cut_variable = cut_var, cut_value = cut_val, left = left_child, right = right_child))
}

list_KD_tree <- function(data, depth = 0) {
  data <- data %>% select_if(is.numeric)
  dim <- ncol(data)
  
  if(!dim) return(NULL)
  
  split_var <- (depth %% dim) + 1
  
  
  
}


############ helper functions ############
# get the data from the leaf of the tree defined by node, 
# get the cuts at each split leading to the leaf of the tree defined by node,
# and get the split directions at each split leading to the leaf of the tree 
# defined by node, with original data set data, where the point is located
# return everything in list
get_data_cuts_and_direction_to <- function(node, data, point) {
  direction <- c()
  cuts <- list()
  
  while(!is.null(node)) {
    if(point[[node$cut_variable]] < node$cut_value) {
      data <- data %>% filter(get(node$cut_variable) < node$cut_value)
      cuts <- c(cuts, list(cut_variable = node$cut_variable, cut_value = node$cut_value))
      direction <- c(direction,"L")
      node <- node$left
    } else {
      data <- data %>% filter(get(node$cut_variable) >= node$cut_value)
      cuts <- c(cuts, list(cut_variable = node$cut_variable, cut_value = node$cut_value))
      direction <- c(direction,"R")
      node <- node$right
    }
  }
  
  return(list(data_in_leaf = data, directions_to_leaf = direction, cuts_to_leaf = cuts))
}

# get the data from the leaf of the tree defined by node, 
# with original data set data, where the point is located
get_leaf_data_at <- function(node, data, point) {
  direction <- c()
  cuts <- list()
  
  while(!is.null(node)) {
    if(point[[node$cut_variable]] < node$cut_value) {
      data <- data %>% filter(get(node$cut_variable) < node$cut_value)
      cuts <- c(cuts, list(cut_variable = node$cut_variable, cut_value = node$cut_value))
      direction <- c(direction,"L")
      node <- node$left
    } else {
      data <- data %>% filter(get(node$cut_variable) >= node$cut_value)
      cuts <- c(cuts, list(cut_variable = node$cut_variable, cut_value = node$cut_value))
      direction <- c(direction,"R")
      node <- node$right
    }
  }
  
  return(data)
}

# get the cuts at each split leading to the leaf of the tree defined by node, 
# with original data set data, where the point is located
get_cuts_to <- function(node, data, point) {
  cuts <- list()
  
  while(!is.null(node)) {
    if(point[[node$cut_variable]] < node$cut_value) {
      cuts <- c(cuts, list(cut_variable = node$cut_variable, cut_value = node$cut_value))
      node <- node$left
    } else {
      cuts <- c(cuts, list(cut_variable = node$cut_variable, cut_value = node$cut_value))
      node <- node$right
    }
  }
  
  return(cuts)
}

# get the split directions at each split leading to the leaf of the tree defined by node, 
# with original data set data, where the point is located
get_direction_to <- function(node, data, point) {
  direction <- c()
  
  while(!is.null(node)) {
    if(point[[node$cut_variable]] < node$cut_value) {
      direction <- c(direction,"L")
      node <- node$left
    } else {
      direction <- c(direction,"R")
      node <- node$right
    }
  }
  
  return(direction)
}