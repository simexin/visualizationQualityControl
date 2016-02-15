#' pairwise correlation
#'
#' given a data matrix, calculate inter-row correlation values.
#'
#' @param data_matrix the data
#' @param use what data to use, "pairwise" or "complete"
#' @param exclude_na should NA values be excluded (default TRUE)
#' @param exclude_inf should Inf values be excluded (default TRUE)
#' @param exclude_0 should 0 values be excluded (default FALSE)
#' @param method which method of correlation to use
#' 
#' @details The function returns a named list with:
#'   \describe{
#'     \item{cor}{the correlation matrix}
#'     \item{count}{how many points were used in the correlation}
#'     \item{keep}{a logical matrix indicating which points passed filtering}
#'     }
#' 
#'
#' @return list
#' @export
pairwise_correlation <- function(data_matrix, use = "pairwise", exclude_na = TRUE, exclude_inf = TRUE, exclude_0 = FALSE, method = "pearson"){
  n_entry <- nrow(data_matrix)
  out_cor <- matrix(0, nrow = nrow(data_matrix), ncol = nrow(data_matrix))
  rownames(out_cor) <- colnames(out_cor) <- rownames(data_matrix)
  
  out_count <- out_cor
  
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
        out_count[i, j] <- sum(use_vals)
      }
    }
  }
  return(list(cor = out_cor, count = out_count, keep = keep_vals))
}

#' pairwise correlation keep both zero
#'
#' given a data matrix, calculate inter-row correlation values, but if both
#' entries have zero, keep it
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
pairwise_correlation_both0 <- function(data_matrix, use = "pairwise", exclude_na = TRUE, exclude_inf = TRUE, exclude_0 = FALSE, method = "pearson"){
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
      both_zero <- zero_loc[i, ] & zero_loc[j, ]
      use_vals <- (keep_vals[i, ] & keep_vals[j, ]) | both_zero
      
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
#' @param use how to work with the data
#' @param exclude_na should NA values be excluded (default TRUE)
#' @param exclude_inf should Inf values be excluded (default TRUE)
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

# calculate median correlations
