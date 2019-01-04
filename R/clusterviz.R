#' Make a Heatmap of the pairwise Jaccard distance between cluster ident of two
#' Seurat object
#'
#'
#' @param ident1 a named factor vector. names are the cell names, the values are
#' the cluster id.
#' @param ident2 a named factor vector. names are the cell names, the values are
#' the cluster id.
#' @param col_low Color for low Jaccard index.
#' @param col_high Color for high Jaccard index.
#' @param title  The title of the heatmap
#' @param cluster_rows  cluster row or not, default FALSE
#' @param cluster_columns cluster columns or not, default FASLE
#' @param show_column_dend Whether or not show column dendrogram
#' @param show_row_dend Whether or not show row dendrogram
#' @param best_match Whether or not only show the best match of ident1 from ident2.
#' if set to TRUE, the Jaccard index matrix will be subsetted using the ident2 column
#' from the output of \code{\link{MatchClusters}}, the row order will be in order from cluster
#' 0 to the total number of clusters, the columns will be the best match of ident1 from ident2,
#' and the columns idents could be duplicated. e.g. single cluster from ident2 matches multiple
#' clusters in ident1.
#' @param ... other parameters pass to \code{\link[ComplexHeatmap]{Heatmap}}
#'
#' @return A Heatmap representing the pair-wise Jaccard correlation, rows are ident1,
#' columns are ident2
#' @export
#'
#' @examples
#'
PairWiseJaccardSetsHeatmap<- function(ident1, ident2, best_match = FALSE,
                                      title = NULL, col_low = "white", col_high= "red",
                                      cluster_rows = F, cluster_columns =F,
                                      show_row_dend = F, show_column_dend = F, ...){
        cell_fun = function(j, i, x, y, width, height, fill) {
                grid::grid.rect(x = x, y = y, width = width *0.99, height = height *0.99,
                          gp = grid::gpar(col = "grey", fill = fill, lty = 1, lwd = 0.5))
        }
        mat<- PairWiseJaccardSets(ident1, ident2)
        col_fun<- circlize::colorRamp2(c(0, 1), c(col_low, col_high))
        if (best_match){
                cluster_rows = F
                cluster_columns =F
                show_row_dend = F
                show_column_dend = F
                match_idx<- MatchClusters(ident1, ident2)
                ComplexHeatmap::Heatmap(mat[, match_idx$ident2],
                                        cluster_rows = cluster_rows, cluster_columns = cluster_columns,
                                        show_row_names = T, show_column_names = T,
                                        show_row_dend = show_row_dend,
                                        show_column_dend = show_column_dend,
                                        col = col_fun, rect_gp = grid::gpar(type = "none"),
                                        cell_fun = cell_fun,
                                        name = "Jaccard distance",
                                        column_title = title,
                                        heatmap_legend_param = list(color_bar = "discrete"),
                                        ...)
        }
        else{
                ComplexHeatmap::Heatmap(mat,
                                        cluster_rows = cluster_rows, cluster_columns = cluster_columns,
                                        show_row_names = T, show_column_names = T,
                                        show_row_dend = show_row_dend,
                                        show_column_dend = show_column_dend,
                                        col = col_fun, rect_gp = grid::gpar(type = "none"),
                                        cell_fun = cell_fun,
                                        name = "Jaccard distance",
                                        column_title = title,
                                        heatmap_legend_param = list(color_bar = "discrete"),
                                        ...)

        }

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
        SelectHighestJaccard<- function(mat){
                apply(mat, 1, max)

        }
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

## see https://jokergoo.github.io/circlize_book/book/the-chorddiagram-function.html

#' Plot ChordDiagram of cell identity changes between two runs of clusters.
#'
#' @param ident1 a named factor vector. names are the cell names, the values are
#' the cluster id.
#' @param ident2 a named factor vector. names are the cell names, the values are
#' the cluster id.
#' @param clusters_to_show_ident1 A character vector of cluster ids to show for ident1.
#' default is NULL, all clusters will be shown.
#' @param big.gap Gap between sectors of two cluster runs.
#' @param transparency Transparency of link colors, 0 means no transparency and 1 means full transparency.
#' see \code{\link[circlize]{chordDiagramFromMatrix}}
#' @param grid.col Grid colors which correspond to matrix rows/columns (or sectors).
#' The length of the vector should be either 1 or length(union(rownames(mat), colnames(mat))).
#' It's preferred that grid.col is a named vector of which names correspond to sectors.
#' If it is not a named vector, the order of grid.col corresponds to order of sectors.
#' see \code{\link[circlize]{chordDiagramFromMatrix}}
#' @param link.sort whether sort links on every sector based on the width of the links on it.
#' If it is set to "overall", all links are sorted regardless whether they are from rows or columns.
#' see \code{\link[circlize]{chordDiagramFromMatrix}}
#' @param link.decreasing for link.sort
#' @param directional Whether links have directions. 1 means the direction is from the first column
#' in df to the second column, -1 is the reverse, 0 is no direction, and 2 for two directional.
#' see \code{\link[circlize]{chordDiagramFromMatrix}}
#'
#' @return A data frame which contains positions of links. see \code{\link[circlize]{chordDiagramFromMatrix}}
#' @export
#'
#' @examples
ClusterIdentityChordPlot<- function(ident1, ident2,
                                    clusters_to_show_ident1 = NULL,
                                    big.gap = 10, transparency = 0.5,
                                    grid.col = NULL,
                                    link.sort = TRUE, link.decreasing = TRUE,
                                    directional = -1){
        mat<- PairWiseOverlappingIdents(ident1, ident2)
        if (!is.null(clusters_to_show_ident1)){
                mat<- mat[clusters_to_show_ident1, ]
        }
        rownames(mat)<- paste0("1_", rownames(mat))
        colnames(mat)<- paste0("2_", colnames(mat))
        circlize::circos.par(start.degree = 90, clock.wise = FALSE)
        circlize::chordDiagram(mat, big.gap = big.gap, transparency = transparency,
                               grid.col = grid.col,
                               link.sort = link.sort, link.decreasing = link.decreasing,
                               directional = directional)
        circlize::circos.clear()
}
