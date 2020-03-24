
#' Make a Barplot for cluster size
#'
#' @param ident a named factor vector. names are the cell names, the values are
#' the cluster id.
#' @param bar_col color for the bar. Default is blue.
#' @param label_number whether or not put cell number in each cluster on top of the bar
#'
#' @return a ggplot2 bar graph object
#' @export
#'
#' @examples
#' data(pbmc_small)
#' CLusterSizeBarplot(pbmc_small@@ident)
#'
ClusterSizeBarplot<- function(ident, bar_col = "blue", label_number = TRUE){
        g<- as.data.frame(table(ident)) %>%
                dplyr::rename(cluster = ident, size = Freq) %>%
                ggplot2::ggplot(aes(x = cluster, y = size)) +
                ggplot2::geom_bar(stat = "identity", fill = bar_col) +
                ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
        if (!label_number){
                return (g)

        } else {
                g<- g + ggplot2::geom_text(aes(label=size), vjust= -1.5, angle = 45)
                return (g)

        }
}

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
                                        name = "Jaccard index",
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
                                        name = "Jaccard index",
                                        column_title = title,
                                        heatmap_legend_param = list(color_bar = "discrete"),
                                        ...)

        }

}


#' Plot the Jaccard index distribution using raincloud plot
#'
#' @param idents1 A list of cluster identity from the subsampled data set
#' before reclustering. (cluster id copied from the original full data set)
#' @param idents2 A list of cluster identity from the subsampled data sets after
#' reclustering.
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
#'JaccardRainCloudPlot(idents, idents)
#'}
#'
JaccardRainCloudPlot<- function(idents1, idents2, title= NULL){
        mats<- AssignHighestJaccard(idents1, idents2)
        g<- mats %>% tibble::as_tibble() %>% tibble::rownames_to_column(var = "bootstrap")  %>%
                tidyr::gather(-bootstrap, key= "cluster", value = "jaccard") %>%
                dplyr::mutate(cluster = as.factor(as.numeric(.$cluster))) %>%
                ggplot2::ggplot(aes(x = cluster, y = jaccard, fill = cluster)) +
                geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
                ggplot2::geom_point(aes(y = jaccard, color = cluster), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
                ggplot2::geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
                ggplot2::theme_classic() +
                ggplot2::theme(legend.position="none") +
                ggplot2::ggtitle(title)
        return(g)
}


#' Plot a scatter plot for different clustering parameters
#'
#' x-axis is the parameters tested (e.g. many different k.param)
#' y-axis is the total number of clusters and total number of stable clusters based
#' on the jaccard cutoff as determined by AssignStableClusters, or precentage of cells
#' in stable clusters.
#'
#' @param stable_clusters a dataframe with list-columns for data, stable_cluster determined by
#' \code{\link{AssignStableCluster}} and the rest of the columns are pc, resolution and k_param.
#' @param fullsample_idents a dataframe with the list-column contain the original ident for
#' the full dataset. This is the direct output from the Snakemake workflow.
#' @param x_var one of "pc", "resolution" and "k_param".
#' @param y_var one of "number" or "percentage". If it is "number",
#' y-axis si the total number of clusters and total number of stable clusters.
#' @param facet_rows one of "pc", "resolution" and "k_param" for ggplot2 to facet.
#' @param facet_cols one of "pc", "resolution" and "k_param" for ggplot2 to facet.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
ParameterSetScatterPlot<- function(stable_clusters,
                                   fullsample_idents,
                                   x_var,
                                   y_var,
                                   facet_rows,
                                   facet_cols ) {

        df<- dplyr::left_join(stable_clusters, fullsample_idents) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(total = map_dbl(stable_cluster, ~ length(.x$stable_cluster))) %>%
                dplyr::mutate(stable = map_dbl(stable_cluster, ~ .x$number_of_stable_cluster)) %>%
                dplyr::mutate(percentage = map2_dbl(original_ident_full, stable_cluster,
                                                    function(x, y) CalculatePercentCellInStable(x,                                                                                      y$stable_cluster))) %>%
                dplyr::select(-data, - stable_cluster, -original_ident_full) %>%
                dplyr::mutate_if(is.character, function(x) as.factor(as.numeric(x))) %>%
                tidyr::gather(total:stable , key = "category", value = "number")
        ## plotting

        if (!all(c(x_var, y_var, facet_rows, facet_cols) %in% colnames(df))) {
                stop("x_var, faect_rows and facet_cols must be one of the parameter columns in the dataframe,\n
         y_var must be 'number' or 'percentage'.")
        }

        if (y_var == "percentage") {
                p<- ggplot2::ggplot(df, aes(x=.data[[x_var]], y = .data[[y_var]])) +
                        ggplot2::geom_point(color = "blue") +
                        ggplot2::geom_line(aes(group = 1), color = "red") +
                        ggplot2::scale_y_continuous(labels = scales::percent) +
                        ggplot2::facet_grid(rows = vars(.data[[facet_rows]]), cols = vars(.data[[facet_cols]])) +
                        ggplot2::xlab(x_var) +
                        ggplot2::ylab(y_var)
        }
        if (y_var == "number"){
                p<- ggplot2::ggplot(df, aes(x=.data[[x_var]], y = .data[[y_var]])) +
                        ggplot2::geom_point() +
                        ggplot2::geom_line(aes(group = category, color = category )) +
                        ggplot2::facet_grid(rows = vars(.data[[facet_rows]]), cols = vars(.data[[facet_cols]])) +
                        ggplot2::xlab(x_var) +
                        ggplot2::ylab(y_var)
        }

        return(p)

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


#' Plot raincloud plot for silhouette score
#'
#' @param silhouette_score a dataframe returned by \code{link[CalculateSilhouette]}
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#'
#' SilhouetteRainCloudPlot(CalculateSilhouette(pbmc_small, dims = 1:15))
SilhouetteRainCloudPlot<- function(silhouette_score){
                g<- ggplot2::ggplot(silhouette_score, aes(x = cluster, y = width, fill = cluster)) +
                geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
                ggplot2::geom_point(aes(y = width, color = cluster), position = position_jitter(width = .15), size = .5,
                           alpha = 0.8) +
                ggplot2::geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
                ggplot2::ylab("silhouette width") +
                ggplot2::theme_classic(base_size = 14) +
                ggplot2::theme(legend.position="none")
                return(g)
}
