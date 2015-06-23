#' create set of disjoint colors
#' 
#' When multiple sample classes need to be visualized on a heatmap, it is useful
#' to be able to distinguish them by color. This function generates a set of colors
#' for sample classes
#' 
#' @param n_group how many groups should there be colors for
#' 
#' @importFrom colorspace rainbow_hcl
#' @export
generate_group_colors <- function(n_group){
  end_color <- 360 * (n_group - 1) / n_group
  group_color <- rainbow_hcl(n_group, c = 100, start = 0, end = end_color)
  group_color <- sample(group_color)
  group_color
}

#' cluster and reorder
#' 
#' given a matrix (maybe a distance matrix), cluster and then re-order using
#' dendsort.
#' 
#' @param similarity_matrix matrix of similarities
#' @param matrix_indices indices to reorder
#' @param transform should a transformation be applied to the data first
#' 
#' @import dendsort
#' @importFrom stats as.dist hclust as.dendrogram
#' @export
similarity_reorder <- function(similarity_matrix, matrix_indices=NULL, transform = "none"){
  if (is.null(matrix_indices)){
    matrix_indices <- seq(1, nrow(similarity_matrix))
  }
  
  if (class(similarity_matrix) != "dist"){
    similarity_matrix <- as.dist(similarity_matrix)
  }
  
  transform_data <- switch(none = similarity_matrix,
                           inverse = 1 / similarity_matrix,
                           sub_1 = 1 - similarity_matrix,
                           log = log(similarity_matrix))
  
  tmp_clust <- as.dendrogram(hclust(transform_data))
  new_sort <- dendsort(tmp_clust)
  matrix_indices[new_sort]
}


#' reorder by sample class
#' 
#' to avoid spurious visualization problems, it is useful in a heatmap visualization
#' to reorder the samples within each sample class. This function uses 
#' hierarchical clustering and \link{dendsort} to sort entries in a heatmap.
#' 
#' @param similarity_matrix matrix of similarities between objects
#' @param sample_classes list of indices of sample classes
#' @param transform a transformation to apply to the data
#' 
#' @import dendsort
#' @export
similarity_reorderbyclass <- function(similarity_matrix, sample_classes, transform="none"){
  new_order <- lapply(sample_classes, function(x){
    reorder_cluster(similarity_matrix[x, x], x, transform = transform)
  })
  
  out_order <- unlist(new_order)
}
