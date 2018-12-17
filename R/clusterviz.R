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
                                      show_row_dend = T, show_column_dend = T){
        cell_fun = function(j, i, x, y, width, height, fill) {
                grid::grid.rect(x = x, y = y, width = width *0.99, height = height *0.99,
                          gp = grid::gpar(col = "grey", fill = fill, lty = 1, lwd = 0.5))
        }

        col_fun<- circlize::colorRamp2(c(0, 1), c(col_low, col_high))
        ComplexHeatmap::Heatmap(mat, cluster_rows = T, cluster_columns = T,
                                show_row_names = T, show_column_names = T,
                                show_row_dend = show_row_dend,
                                show_column_dend = show_column_dend,
                                col = col_fun, rect_gp = grid::gpar(type = "none"),
                                cell_fun = cell_fun,
                                name = "Jaccard distance",
                                column_title = title,
                                heatmap_legend_param = list(color_bar = "discrete"))


}


#' Plot the Jaccard index distribution using raincloud plot
#'
#' @param ident1 The cluster identity from the original full data set.
#' @param idents A list of cluster identity from the subsampled data sets.
#' @param title Title of the plot
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#'\dontrun{
#'data(idents)
#'## the pbmc here need to be fully processed.
#'JaccardRainCloudPlot(pbmc@@ident, idents)
#'}
#'
JaccardRainCloudPlot<- function(ident1, idents, title= NULL){
        mat_list<- purrr::map(idents, ~PairWiseJaccardSets(ident1 = ident1, ident2 = .x))
        mat_max<- purrr::map(mat_list, SelectHighestJaccard)
        mats<- purrr::reduce(mat_max, dplyr::bind_rows)

        g<- mats %>% tibble::as_tibble() %>% tibble::rownames_to_column(var = "bootstrap")  %>%
                tidyr::gather(-bootstrap, key= "cluster", value = "jaccard") %>%
                ggplot(aes(x = cluster, y = jaccard, fill = cluster)) +
                geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
                geom_point(aes(y = jaccard, color = cluster), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
                geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
                theme(legend.position="none") +
                ggtitle(title)
        return(g)
}


#' Scatter plot of number of stable clusters across different k
#' and percentage of cells in stable clusters across different k.
#'
#' @param stable_cluster A list of list returned by AssignStableCluster
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' ## see README.md
#' AssignStableCluster(ks_idents_original$k20, ks_idents$`20`)
#'
#'ks_stable<- purrr::map2(ks_idents_original, ks_idents, ~AssignStableCluster(ident1= .x, idents = .y))
#'}
BootParameterScatterPlot<- function(stable_cluster){
        percent<- purrr::map_dbl(stable_cluster, "percent_cell_in_stable")
        stable<- purrr::map_dbl(stable_cluster, "number_of_stable_cluster")
        g<- dplyr::bind_rows(percent, stable) %>%
                tibble::add_column(parameter = c("percent", "stable_cluster_num")) %>%
                tidyr::gather(-parameter, key = "bootpar", value= "value" ) %>%
                ggplot(aes(x = bootpar, y = value)) +
                geom_line(aes(group = parameter, col = parameter)) +
                geom_point() +
                facet_wrap(~parameter, scale = "free")
        return(g)

}
