
#' Calculate Silhouette width from PCA space for each cell after clustering
#' This is calculated from Seurat object
#' @param object A Seurat object with Idents set to cluster ids (factors)
#' @param dims default 1:50  dimension to use in the PCA space to calculate
#' eucledian distance
#'
#' @return a dataframe with silhouette width for each cell. see also \code{\link[cluster]{silhouette}}
#' @export
#'
#' @examples
#' CalculateSilhouette(pbmc_small, dims = 1:15)
#'
CalculateSilhouette<- function(object, dims = 1:50){
        if (length(dims) > ncol(object@reductions$pca@cell.embeddings)) {
                stop("please specify PCA dims smaller than calculated")
        }
        cell_distance<- dist(object@reductions$pca@cell.embeddings[, dims])
        # or as.integer
        cell_cluster<- as.numeric(as.character(Idents(object)))
        silhouette_score<- cluster::silhouette(cell_cluster, cell_distance)
        silhouette_score<- tibble::tibble(cluster = silhouette_score[,1],
                                          width = silhouette_score[,3],
                                          cell = colnames(object)) %>%
                dplyr::mutate(cluster = as.factor(cluster))
        return(silhouette_score)
}
