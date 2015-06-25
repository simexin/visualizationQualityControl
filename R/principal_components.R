#' sensible principal components
#' 
#' given the results of \code{\link{prcomp}}, generate decent looking 
#' plot of selected components using \code{ggbiplot}.
#' 
#' @param pca_decomp the results of \code{prcomp}
#' @param princomps which principal components to plot
#' @param obj_class factor defining the classes of each object
#' @param dot_size size of the points (default is 1)
#'
#' @import ggbiplot
#' @import ggplot2
#' @export
visqc_pca <- function(pca_decomp, princomps = c(1, 2), obj_class = NULL, dot_size = 1){
  x_data <- pca_decomp$x
  n_obj <- nrow(x_data)
  
  if (is.null(obj_class)){
    obj_class <- data.frame(groups = as.factor(rep("", n_obj)))
  }
  
  
  ggbiplot(pca_decomp, obs.scale = 1, groups = obj_class$groups, var.axes = FALSE, choices = princomps) + geom_point(aes(color = groups), size = dot_size)
}
