#' find non-zero by percentage
#' 
#' Given a value matrix (features are rows, samples are columns), and sample classes, 
#' find those things that are not zero in at least a certain 
#' number of one of the classes
#' 
#' @param data_matrix the matrix of values to work with 
#' @param sample_classes the classes of each sample
#' @param keep_perc what percentage of samples in each class need a non-zero value
#' 
#' @export
#' @return matrix
filter_non_zero_percentage <- function(data_matrix, sample_classes, keep_perc = 0.75){
  uniq_classes <- unique(sample_classes)
  class_index <- lapply(uniq_classes, function(x){sample_classes %in% x})
  names(class_index) <- uniq_classes
  
  min_notzero <- sapply(class_index, function(x){round(sum(x) * keep_perc)})
  
  has_min <- apply(data_matrix, 2, function(in_col){
    n_pass <- sapply(class_index, function(index){sum(in_col[index] > 0)})
    keep <- sum(n_pass >= min_notzero) > 0
    keep
  })
  data_matrix[, has_min]
}
