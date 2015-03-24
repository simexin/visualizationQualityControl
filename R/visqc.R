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
  
  meltedData <- melt(dataMatrix, as.is = TRUE)
  
  if (is.character(meltedData[, "Var1"])){
    meltedData$Var1 <- factor(meltedData$Var1, levels = rownames(dataMatrix), ordered = TRUE)
    meltedData$Var2 <- factor(meltedData$Var2, levels = colnames(dataMatrix), ordered = TRUE)
  }
  
  # checking if we have numeric values as indices, do they correspond to the
  # actual indices, and if not, change them to character
  if (is.integer(meltedData$Var1)){
    if (sum(dataIndices %in% meltedData$Var1) != nEntry){
      meltedData$Var1 <- factor(as.character(meltedData$Var1), levels = rownames(dataMatrix), ordered = TRUE)
      meltedData$Var2 <- factor(as.character(meltedData$Var2), levels = colnames(dataMatrix), ordered = TRUE)
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
