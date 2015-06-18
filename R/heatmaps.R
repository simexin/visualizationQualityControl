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
