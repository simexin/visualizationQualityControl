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

#' heatmap
#' 
#' generates a heatmap using ggplot2
#' 
#' @param dataMatrix a data matrix with values for the heatmap
#' @param cols the colors to use
#' @param limits the limits to use for the data
#' @return ggplot2 object
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
gg_heatmap <- function(dataMatrix, cols = NULL, limits = NULL){
  if (is.null(limits)){
    limits <- range(dataMatrix)
  }
  
  if (is.null(cols)){
    cols = grey.colors(10, start = 0, end = 1)
  }
  
  nEntry <- nrow(dataMatrix)
  dataIndices <- seq(1, nEntry)
  
  meltedData <- melt(dataMatrix)
  
  # checking if we have numeric values as indices, do they correspond to the
  # actual indices, and if not, change them to character
  if (is.integer(meltedData$Var1)){
    if (sum(dataIndices %in% meltedData$Var1) != nEntry){
      meltedData$Var1 <- as.character(meltedData$Var1)
      meltedData$Var2 <- as.character(meltedData$Var2)
    }
  }
  
  outPlot <- ggplot(meltedData, aes(x = Var1, y = Var2, fill = value)) + geom_tile()
  outPlot <- outPlot + scale_fill_gradientn(colours = cols, limits = limits)
  # outPlot <- outPlot + scale_x_continuous(expand = c(0, 0)) + 
  #  scale_y_continuous(expand = c(0, 0))
  outPlot <- outPlot + coord_equal()
  outPlot <- outPlot + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  return(outPlot)
}
