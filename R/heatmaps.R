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
#' 
#' @examples 
#' set.seed(1234)
#' mat <- matrix(rnorm(100, sd = 0.5), 10, 10)
#' rownames(mat) <- colnames(mat) <- letters[1:10]
#' 
#'
similarity_reorderbyclass <- function(similarity_matrix, sample_classes=NULL, transform="none"){
  if (is.null(sample_classes)){
    sample_classes <- list(none = seq(1, nrow(similarity_matrix)))
  }
  new_order <- lapply(sample_classes, function(x){
    reorder_cluster(similarity_matrix[x, x], x, transform = transform)
  })
  
  out_order <- unlist(new_order)
}


#' easier heatmaps
#' 
#' rolls some of the common \code{Heatmap} options into a single function call
#' to make life easier when creating lots of heatmaps. \strong{Note:} this function
#' does \emph{reorder} the rows or columns of the matrix before visualization,
#' the order going in is what will be presented.
#' 
#' @param matrix_data the matrix you want to plot as a heatmap
#' @param color_values the color mapping of values to colors (see Details)
#' @param row_color_data data for row annotations
#' @param row_color_list list for row annotations
#' @param col_color_data data for column annotations
#' @param col_color_list list for column annotations
#' 
#' @details This function uses the \code{ComplexHeatmap} package to produce
#' heatmaps with complex row- and column-color annotations. Both \code{row_color_data}
#' and \code{col_color_data} should be \code{data.frame}'s where each column describes
#' meta-data about the rows or columns of the matrix. The \code{row_color_list} and 
#' \code{col_color_list} provide the mapping of color to annotation, where each
#' \code{list} entry should be a named vector of colors, with the list entry
#' corresponding to a column entry in the data.frame, and the names of the colors
#' corresponding to annotations in that column.
#' 
#' @examples  
#' \dontrun{
#' library(circlize)
#' colormap <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
#' mat <- matrix(rnorm(100, sd = 0.5), 10, 10)
#' rownames(mat) <- colnames(mat) <- letters[1:10]
#' 
#' meta_data <- data.frame(grp = rep(c("grp1", "grp2"), each = 5),
#'                         set = c(rep("set1", 3), rep("set2", 7)))
#' annotation_color <- c(grp1 = "green", grp2 = "red", set1 = "blue",
#'                       set2 = "yellow")
#' 
#' row_data <- meta_data[, "grp", drop = FALSE]
#' col_data <- meta_data[, "set", drop = FALSE]
#' row_annotation = list(grp = annotation_color[1:2])
#' col_annotation = list(set = annotation_color[3:4])
#' 
#' generate_heatmap(mat, colormap, row_color_data = row_data, row_color_list = row_annotation,
#'                  col_color_data = col_data, col_color_list = col_annotation)
#' }
#' 
#' @import ComplexHeatmap
#' @export
generate_heatmap <- function(matrix_data, color_values, title = "", row_color_data = NULL, row_color_list, col_color_data = NULL, col_color_list){
  if (!is.null(row_color_data)){
    row_annot <- rowAnnotation(df = row_color_data, col = row_color_list)
  }
  if (!is.null(col_color_data)){
    col_annot <- HeatmapAnnotation(df = col_color_data, col = col_color_list)
  }
  
  heat_out <- Heatmap(matrix_data, col = color_values, cluster_rows = FALSE, cluster_columns = FALSE,
                      top_annotation = col_annot, column_title = title) + row_annot
  heat_out
}
