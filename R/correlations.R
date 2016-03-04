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

#' calculate median correlations
#' 
#' Given a correlation matrix and optionally the sample class information,
#' calculates the median correlations of each sample to all other samples in
#' the same class. May be useful for determining outliers.
#' 
#' @param cor_matrix the sample - sample correlations
#' @param sample_classes the sample classes as a character or factor
#' 
#' @return data.frame
#' @export
#' 
#' @details The data.frame returned has two columns:
#' \describe{
#'   \item{med_cor}{the median correlation with other samples}
#'   \item{sample_class}{the class of the sample. If not provided, set to "C1"}
#' }
#' 
median_correlations <- function(cor_matrix, sample_classes = NULL){
  stopifnot(nrow(cor_matrix) == ncol(cor_matrix))
  n_sample <- nrow(cor_matrix)
  
  if (is.null(sample_classes)) {
    use_classes <- factor(rep("C1", n_sample))
  } else if (is.factor(sample_classes)) {
    use_classes <- sample_classes
  } else {
    sample_classes <- factor(sample_classes)
    use_classes <- sample_classes
  }
  
  split_classes <- split(seq(1, n_sample), use_classes)
  
  median_values_class <- lapply(split_classes, function(class_index){
    class_matrix <- cor_matrix[class_index, class_index]
    median_values <- vapply(seq(1, nrow(class_matrix)), function(in_row){
      median(class_matrix[in_row, -in_row])
    }, numeric(1))
    
    if (!is.null(rownames(class_matrix))) {
      names(median_values) <- rownames(class_matrix)
    }
    median_values
  })
  
  if (is.null(rownames(cor_matrix))) {
    sample_id <- seq(1, n_sample)
  } else {
    sample_id <- rownames(cor_matrix)
  }
  
  out_values <- data.frame(med_cor = unlist(median_values_class, use.names = FALSE),
                           sample_id = sample_id,
                           sample_class = use_classes)
  
  out_values
}

#' fraction of outliers
#' 
#' Calculates the fraction of entries in each sample that are more than \code{X}
#' standard deviations from the trimmed mean. See Details.
#' 
#' @param data the data matrix (samples are rows, columns are features)
#' @param sample_classes the sample classes
#' @param n_trim how many features to trim at each end (default is 3)
#' @param n_sd how many SD before treated as outlier (default is 5)
#' @param remove_0 should zeros be removed before calculating? (default is TRUE)
#' 
#' @details Based on the Gerlinski paper \href{https://dx.doi.org/10.1093/bioinformatics/btv425}{link}
#' for each feature (in a sample class), take the range across all the samples,
#' remove the \code{n_trim} lowest and highest values, and calculate the \code{mean}
#' and \code{sd}, and the actual upper and lower ranges of \code{n_sd} from the
#' \code{mean}. For each sample and feature, determine if \emph{within} or \emph{outside}
#' that limit. Fraction is reported as the number of features outside the range.
#' 
#' @export
#' @return data.frame
outlier_fraction <- function(data, sample_classes = NULL, n_trim = 3,
                             n_sd = 5, remove_0 = FALSE){
  n_sample <- nrow(data)
  
  if (is.null(sample_classes)) {
    use_classes <- factor(rep("C1", n_sample))
  } else if (is.factor(sample_classes)) {
    use_classes <- sample_classes
  } else {
    sample_classes <- factor(sample_classes)
    use_classes <- sample_classes
  }
  
  split_classes <- split(seq(1, n_sample), use_classes)
  
  calc_outlier <- function(data, n_trim, n_sd, remove_0){
    outlier_data <- apply(data, 2, function(x){
      if (remove_0) {
        y <- x[x != 0]
      } else {
        y <- x
      }
      n_y <- length(y)
      y <- sort(y)
      y_start <- n_trim + 1
      y_end <- n_y - (n_trim + 1)
      if ((y_end <= y_start) || (y_end <= 0)) {
        y_start <- 1
        y_end <- n_y
      }
      # need to fix cases with negatives
      y <- y[y_start:y_end]
      y_mean <- mean(y)
      y_sd <- sd(y)
      y_lo <- y_mean - (n_sd * y_sd)
      y_hi <- y_mean + (n_sd * y_sd)
      
      x_out <- !((x >= y_lo) & (x <= y_hi))
      
      if (remove_0) {
        x_out[x == 0] <- FALSE
      }
      x_out
    })
    return(outlier_data)
  }
  
  frac_outlier_class <- lapply(names(split_classes), function(class_name){
    class_index <- split_classes[[class_name]]
    is_outlier <- calc_outlier(data[class_index, , drop = FALSE], n_trim, n_sd, remove_0)
    frac_outlier <- rowSums(is_outlier) / ncol(data)
    data.frame(sample = class_index, class = class_name, frac = frac_outlier)
  })
  
  frac_outlier <- do.call(rbind, frac_outlier_class)
  frac_outlier
}

#' grp_cor_data
#' 
#' Example data used for demonstrating \code{median correlation}. A list with
#' 2 named entries:
#' 
#' \describe{
#'   \item{data}{a data matrix with 100 rows and 20 columns}
#'   \item{class}{a character vector of 20 entries denoting classes}
#' }
#' 
#' The data comes from two groups of samples, where there is ~0.85 correlation
#' within each group, and ~0.53 correlation between groups.
#'  
#' @format List with 2 entries, data and class
#' @source Robert M Flight
"grp_cor_data"

#' grp_exp_data
#' 
#' Example data that requires log-transformation before doing PCA or other QC.
#' A list with 2 named entries:
#' 
#' \describe{
#'   \item{data}{a data matrix with 1000 rows and 20 columns}
#'   \item{class}{a character vector of 20 entries denoting classes}
#' }
#' 
#' The data comes from two groups of samples, where there is ~0.80 correlation
#' within each group, and ~0.38 correlation between groups.
#'  
#' @format List with 2 entries, data and class
#' @source Robert M Flight
"grp_exp_data"
