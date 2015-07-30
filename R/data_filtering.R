#' find non-zero by percentage
#' 
#' Given a value matrix (features are columns, samples are rows), and sample classes, 
#' find those things that are not zero in at least a certain 
#' number of one of the classes
#' 
#' @param data_matrix the matrix of values to work with 
#' @param sample_classes the classes of each sample
#' @param keep_num what number of samples in each class need a non-zero value (see Details)
#' 
#' @details The number of samples that must be non-zero can be expressed either as a whole
#'   number (that is greater than one), or as a fraction that will be be multiplied 
#'   by the number of samples in each class to get the lower limits for each of the classes.
#' 
#' @export
#' @return matrix
filter_non_zero_percentage <- function(data_matrix, sample_classes = NULL, keep_num = 0.75){
  if (is.null(sample_classes)) {
    sample_classes <- rep("A", nrow(data_matrix))
  }
  uniq_classes <- unique(sample_classes)
  class_index <- lapply(uniq_classes, function(x){sample_classes %in% x})
  names(class_index) <- uniq_classes
  
  if (keep_num <= 1) {
    min_notzero <- sapply(class_index, function(x){round(sum(x) * keep_perc)})
  } else {
    min_notzero <- sapply(class_index, function(x){round(keep_num)})
  }
  
  
  has_min <- apply(data_matrix, 2, function(in_col){
    n_pass <- sapply(class_index, function(index){sum(in_col[index] > 0)})
    keep <- sum(n_pass >= min_notzero) > 0
    keep
  })
  data_matrix[, has_min]
}
