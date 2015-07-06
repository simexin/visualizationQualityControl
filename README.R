## ----load_data-----------------------------------------------------------
library(visualizationQualityControl)
data(all_intensity)
data(all_info)

## ----dopca---------------------------------------------------------------
int_zero <- all_intensity
int_zero[is.na(int_zero)] <- 0
has_3 <- apply(int_zero, 1, function(x){sum(x != 0) >= 2})
int_zero <- int_zero[has_3,]
all_pca <- prcomp(t(int_zero), center = TRUE, scale. = TRUE)

## ----plotpca-------------------------------------------------------------
visqc_pca(all_pca, groups = all_info$disease)

## ----calc_correlation----------------------------------------------------
all_cor <- pairwise_correlation(t(all_intensity))
dim(all_cor)

## ----plot_heatmap--------------------------------------------------------
library(circlize)
colormap <- colorRamp2(c(0.95, 1), c("black", "white"))
visqc_heatmap(all_cor, colormap)

## ----reorder-------------------------------------------------------------
all_order <- similarity_reorderbyclass(all_cor, all_info[, "disease", drop = FALSE], transform = "sub_1")
visqc_heatmap(all_cor, colormap, row_order = all_order$indices, column_order = all_order$indices)

## ----add_legend----------------------------------------------------------
color_legend <- generate_group_colors(2)
names(color_legend) <- c("normal", "cancer")

row_data <- all_info[, "disease", drop = FALSE]
row_annotation <- list(disease = color_legend)
visqc_heatmap(all_cor, colormap, row_color_data = row_data, row_color_list = row_annotation, col_color_data = row_data, col_color_list = row_annotation, row_order = all_order$indices, column_order = all_order$indices)

## ----add_legend2---------------------------------------------------------
all_info$random
color_legend <- generate_group_colors(6)
names(color_legend) <- c(unique(all_info$random), unique(all_info$disease))

col_data <- all_info
col_annotation <- list(random = color_legend[1:4], disease = color_legend[5:6])

row_data <- all_info[, "disease", drop = FALSE]
row_annotation <- list(disease = color_legend[5:6])
visqc_heatmap(all_cor, colormap, row_color_data = row_data, row_color_list = row_annotation, col_color_data = col_data, col_color_list = col_annotation, row_order = all_order$indices, column_order = all_order$indices)

