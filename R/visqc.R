#' summarize data
#' 
#' summarizes a matrix or data.frame, where rows are samples
#' and columns are features
#' 
#' @param in_data matrix or data.frame
#' @param sample_classes which samples are in which class
#' @param avg_function which function to use for summary
#' 
#' @return data.frame
#' @export
summarize_data <- function(in_data, sample_classes=NULL, avg_function = mean){
  if (is.null(sample_classes)){
    sample_classes <- rep("A", nrow(in_data))
  }
  
  n_feature <- ncol(in_data)
  split_indices <- split(seq(1, nrow(in_data)), sample_classes)
  
  split_values <- lapply(split_indices, function(use_index){
    tmp_mean <- apply(in_data[use_index, , drop = FALSE], 2, avg_function)
    tmp_sd <- apply(in_data[use_index, , drop = FALSE], 2, sd)
    data.frame(mean = tmp_mean, var = c(tmp_sd, tmp_sd / tmp_mean), type = rep(c("sd", "rsd"), each = n_feature))
  })
  
  out_data <- do.call(rbind, split_values)
  out_data$type <- factor(out_data$type, ordered = TRUE, levels = c("sd", "rsd"))
  out_data
}


#' pairwise correlation
#' 
#' given a data matrix, calculate inter-row correlation values
#' 
#' @param data_matrix the data
#' @param use what data to use, "pairwise" or "complete"
#' @param exclude_na should NA values be excluded (default TRUE)
#' @param exclude_inf should Inf values be excluded (default TRUE)
#' @param exclude_0 should 0 values be excluded (default FALSE)
#' @param method which method of correlation to use
#' 
#' @return matrix
#' @export
pairwise_correlation <- function(data_matrix, use = "pairwise", exclude_na = TRUE, exclude_inf = TRUE, exclude_0 = FALSE, method = "pearson"){
  n_entry <- nrow(data_matrix)
  out_cor <- matrix(0, nrow = nrow(data_matrix), ncol = nrow(data_matrix))
  rownames(out_cor) <- colnames(out_cor) <- rownames(data_matrix)
  diag(out_cor) <- 1
  
  na_loc <- matrix(FALSE, nrow = n_entry, ncol = ncol(data_matrix))
  inf_loc <- na_loc
  zero_loc <- na_loc
  
  if (exclude_na){
    na_loc <- is.na(data_matrix)
  }
  
  if (exclude_inf){
    inf_loc <- is.infinite(data_matrix)
  }
  
  if (exclude_0){
    zero_loc <- data_matrix == 0
  }
  
  exclude_loc <- na_loc | zero_loc | inf_loc
  
  if (use == "complete"){
    keep_vals <- apply(exclude_loc, 2, function(x){sum(!x) == n_entry})
    keep_vals <- matrix(keep_vals, nrow = n_entry, ncol = ncol(data_matrix), byrow = FALSE)
  } else {
    keep_vals <- !exclude_loc
  }
  
  for (i in seq(1, n_entry)){
    for (j in seq(i, n_entry)){
      use_vals <- keep_vals[i, ] & keep_vals[j, ]
      
      if (sum(use_vals) > 1){
        #print(c(i, j))
        out_cor[i, j] <- out_cor[j, i] <- cor(data_matrix[i, use_vals], data_matrix[j, use_vals], method = method)
      }
    }
  }
  out_cor
}

#' pairwise correlation multicore
#' 
#' given a data matrix, calculate inter-row correlation values
#' 
#' @param data_matrix the data
#' @param use what data to use, "pairwise" or "complete"
#' @param exclude_na should NA values be excluded (default TRUE)
#' @param exclude_inf should Inf values be excluded (default TRUE)
#' @param exclude_0 should 0 values be excluded (default FALSE)
#' @param method which method of correlation to use
#' 
#' @return matrix
#' @export
pairwise_correlation_multicore <- function(data_matrix, use = "pairwise", exclude_na = TRUE, exclude_inf = TRUE, exclude_0 = FALSE, method = "pearson"){
  n_entry <- nrow(data_matrix)
  out_cor <- matrix(0, nrow = nrow(data_matrix), ncol = nrow(data_matrix))
  rownames(out_cor) <- colnames(out_cor) <- rownames(data_matrix)
  diag(out_cor) <- 1
  
  na_loc <- matrix(FALSE, nrow = n_entry, ncol = ncol(data_matrix))
  inf_loc <- na_loc
  zero_loc <- na_loc
  
  if (exclude_na){
    na_loc <- is.na(data_matrix)
  }
  
  if (exclude_inf){
    inf_loc <- is.infinite(data_matrix)
  }
  
  if (exclude_0){
    zero_loc <- data_matrix == 0
  }
  
  exclude_loc <- na_loc | zero_loc | inf_loc
  
  if (use == "complete"){
    keep_vals <- apply(exclude_loc, 2, function(x){sum(!x) == n_entry})
    keep_vals <- matrix(keep_vals, nrow = n_entry, ncol = ncol(data_matrix), byrow = FALSE)
  } else {
    keep_vals <- !exclude_loc
  }
  
  all_comparisons <- expand.grid(seq(1, n_entry), seq(1, n_entry))
  all_comparisons <- all_comparisons[(all_comparisons[,2] > all_comparisons[,1]), ]
  all_comparisons <- as.matrix(all_comparisons)
  
  correlation_function <- function(x){
    n_vals <- length(x)
    out_value <- rep(0, n_vals)
    do_comparison <- all_comparisons[x, ]
    
    for (i_val in seq_len(n_vals)){
      i <- do_comparison[1]
      j <- do_comparison[2]
      
      use_vals <- keep_vals[i, ] & keep_vals[j, ]
      
      if (sum(use_vals) > 1){
        #print(c(i, j))
        out_value[i_val] <- cor(data_matrix[i, use_vals], data_matrix[j, use_vals], method = method)
      }
    }
    
    out_value
  }
  
  all_cors <- parallel::pvec(seq(1, nrow(all_comparisons)), correlation_function)
  
  out_cor[all_comparisons[,1], all_comparisons[,2]] <- all_cors
  out_cor[all_comparisons[,2], all_comparisons[,1]] <- all_cors
  out_cor
}

#' pairwise distance
#' 
#' given a data matrix, calculate inter-row distances
#' 
#' @param data_matrix the data
#' @param use what data to use, "pairwise" or "complete"
#' @param exclude_na should NA values be excluded (default TRUE)
#' @param exclude_inf should Inf values be excluded (default TRUE)
#' @param exclude_0 should 0 values be excluded (default FALSE)
#' @param method which method of correlation to use
#' 
#' @return matrix
#' @export
pairwise_distance <- function(data_matrix, use = "pairwise", exclude_na = TRUE, exclude_inf = TRUE, exclude_0 = FALSE){
  n_entry <- nrow(data_matrix)
  out_cor <- matrix(0, nrow = nrow(data_matrix), ncol = nrow(data_matrix))
  rownames(out_cor) <- colnames(out_cor) <- rownames(data_matrix)
  diag(out_cor) <- 1
  
  na_loc <- matrix(FALSE, nrow = n_entry, ncol = ncol(data_matrix))
  inf_loc <- na_loc
  zero_loc <- na_loc
  
  if (exclude_na){
    na_loc <- is.na(data_matrix)
  }
  
  if (exclude_inf){
    inf_loc <- is.infinite(data_matrix)
  }
  
  if (exclude_0){
    zero_loc <- data_matrix == 0
  }
  
  exclude_loc <- na_loc | zero_loc | inf_loc
  
  if (use == "complete"){
    keep_vals <- apply(exclude_loc, 2, function(x){sum(!x) == n_entry})
    keep_vals <- matrix(keep_vals, nrow = n_entry, ncol = ncol(data_matrix), byrow = FALSE)
  } else {
    keep_vals <- !exclude_loc
  }
  
  for (i in seq(1, n_entry)){
    for (j in seq(i, n_entry)){
      use_vals <- keep_vals[i, ] & keep_vals[j, ]
      
      if (sum(use_vals) > 1){
        #print(c(i, j))
        out_cor[j, i] <- out_cor[i, j] <- sqrt(sum((data_matrix[i, use_vals] - data_matrix[j, use_vals]) ^ 2))
      }
    }
  }
  out_cor
}

#' pairwise non-zero
#' 
#' given a data matrix, how many entries are non-zero on a pairwise basis
#' 
#' @param data_matrix the data
#' @param exclude_na should NA values be excluded (default TRUE)
#' @param exclude_inf should Inf values be excluded (default TRUE)
#' @param method which method of correlation to use
#' 
#' @return matrix
#' @export
pairwise_nonzero <- function(data_matrix, use = "pairwise", exclude_na = TRUE, exclude_inf = TRUE){
  n_entry <- nrow(data_matrix)
  out_cor <- matrix(0, nrow = nrow(data_matrix), ncol = nrow(data_matrix))
  rownames(out_cor) <- colnames(out_cor) <- rownames(data_matrix)
  diag(out_cor) <- 1
  
  na_loc <- matrix(FALSE, nrow = n_entry, ncol = ncol(data_matrix))
  inf_loc <- na_loc
  
  if (exclude_na){
    na_loc <- is.na(data_matrix)
  }
  
  if (exclude_inf){
    inf_loc <- is.infinite(data_matrix)
  }
  
  zero_loc <- data_matrix == 0
  
  exclude_loc <- na_loc | inf_loc
  
  keep_vals <- !exclude_loc

  for (i in seq(1, n_entry)){
    for (j in seq(i, n_entry)){
      use_vals <- keep_vals[i, ] & keep_vals[j, ]
      
      zero_both <- zero_loc[i, ] & zero_loc[j, ]
      
      sum_non_zero <- sum(use_vals) - sum(zero_both)
      out_cor[j, i] <- out_cor[i, j] <- sum_non_zero
    }
  }
  out_cor
}

#' calculate F-ratio
#' 
#' given a data matrix of samples (rows) and features (columns), and a vector of classes (character or factor),
#' calculate an F-ratio for each feature.
#' 
#' @param data the data matrix, with samples (rows) and features (columns)
#' @param data_classes what are the classes of the rows
#' 
#' @return vector
#' @export
calculate_fratio <- function(data, data_classes){
  if(is.character(data_classes)){
    data_classes <- factor(data_classes)
  }
  
  all_means <- colMeans(data)
  
  split_indices <- split(seq(1, nrow(data)), data_classes)
  n_sample <- nrow(data)
  n_group <- length(split_indices)
  split_data <- lapply(split_indices, function(in_index){data[in_index, , drop = FALSE]})
  
  group_means <- lapply(split_data, colMeans)
  group_var <- lapply(split_data, function(in_data){apply(in_data, 2, var)})
  group_count <- lapply(split_data, nrow)
  
  weight_var <- function(count, var, sub1 = TRUE){
    if (sub1){
      count <- count - 1
    }
    count * var
  }
  within_var <- Map(weight_var, group_count, group_var)
  within_var <- do.call(rbind, within_var)
  within_var <- colSums(within_var) / (n_sample - n_group)
  
  between_var <- lapply(group_means, function(in_mean){(in_mean - all_means)^2})
  between_var <- Map(weight_var, group_count, between_var, FALSE)
  between_var <- do.call(rbind, between_var)
  between_var <- colSums(between_var) / (n_group - 1)
  
  f_ratio <- between_var / within_var
  f_ratio
}

#' calculate values from summaries
#' 
#' given a data.frame of means and variances, calculate mean sd at low end and
#' mean rsd at high end.
#' 
#' @param data data.frame of means and variances
#' @param low_cut means <= this value used for average sd
#' @param hi_cut means >= this value used for average rsd
#' 
#' @import dplyr
#' @export
#' @return vector
calc_sd_rsd <- function(data, low_cut, hi_cut = NULL){
  if (is.null(hi_cut)){
    hi_cut <- min(data[, "mean"]) > low_cut
  }
  
  sd_mn <- filter(data, mean <= low_cut, type == "sd") %>% summarise(., mean(var))
  rsd_mn <- filter(data, mean >= hi_cut, type == "rsd") %>% summarise(., mean(var))
  
  return(c(sd = sd_mn, rsd = rsd_mn))
}

#' calculate values from summaries v2
#' 
#' given a data.frame of means and variances, use a two step non-linear least squares.
#' The first step is done on the mean vs sd, then the estimates are used in a second
#' that estimates them using the mean vs rsd.
#' 
#' @param data data.frame of means and variances
#' @param ... other nls parameters
#' 
#' @import dplyr
#' @export
#' @return vector
calc_sd_rsd_nls <- function(data, ...){
  nl_sd <- filter(data, type == "sd") %>% nls(var ~ B + A*mean, data = ., start = list(A = 0, B = 0), ...)
  nl_rsd <- filter(data, type == "rsd") %>% nls(var ~ ((A * mean) + B) / mean, data = ., start = list(A = coef(nl_sd)["A"], B = coef(nl_sd)["B"]), ...)
  
  return(c(additive = coef(nl_rsd)["B"], proportional = coef(nl_rsd)["A"]))
}
