#' Make a Heatmap of the pairwise Jaccard distance between cluster ident of two
#' Seurat object
#'
#'
#' @param mat The pairwise Jaccard distance matrix returned by
#' \code{\link{PairWiseJaccardSets}}
#' @param col_low Color for low Jaccard index.
#' @param col_high Color for high Jaccard index.
#' @param title  The title of the heatmap
#' @param show_row_dend Whether or not show row dendrogram
#' @param show_column_dent Whether or not show column dendrogram
#'
#' @return A Heatmap representing the pair-wise Jaccard correlation
#' @export
#'
#' @examples
#'
PairWiseJaccardSetsHeatmap<- function(mat, title = NULL, col_low = "white", col_high= "red",
                                      show_row_dend = T, show_column_dent = T){
        cell_fun = function(j, i, x, y, width, height, fill) {
                grid.rect(x = x, y = y, width = width *0.99, height = height *0.99,
                          gp = gpar(col = "grey", fill = fill, lty = 1, lwd = 0.5))
        }

        col_fun<- circlize::colorRamp2(c(0, 1), c(col_low, col_high))
        ComplexHeatmap::Heatmap(mat, cluster_rows = T, cluster_columns = T,
                                show_row_names = T, show_column_names = T,
                                show_row_dend = show_row_dend,
                                show_column_dend = show_column_dend,
                                col = col_fun, rect_gp = gpar(type = "none"),
                                cell_fun = cell_fun,
                                name = "Jaccard distance",
                                column_title = title,
                                heatmap_legend_param = list(color_bar = "discrete"))


}


