

#' Calculate jaccard distance for two sets of character vectors
#'
#' @param set1 character vector 1
#' @param set2 character vector 2
#'
#' @return jaccard distance
#' @export
#'
#' @examples
#' JaccardSets(sample(LETTERS, 10), sample(LETTERS, 10))
JaccardSets<- function(set1, set2){
        length(intersect(set1, set2))/length(unique(c(set1, set2)))
}


#' Calculate pair-wise Jaccard distance for @@ident slots from two Seurat objects
#'
#' Calculate pair-wise Jaccard distance for two named factor vector. e.g.
#' seurat_obj1@ident and seurat_obj2@ident
#'
#' @param ident1 a named factor vector. names are the cell names, the values are
#' the cluster id.
#' @param ident2  a named factor vector. names are the cell names, the values are
#' the cluster id.
#'
#' @return a matrix of pair-wise Jaccard distance. Rows are clusters from ident1,
#' columns are clusters from ident2
#' @export
#'
#' @examples
#' PairWiseJaccardSets(pbmc@ident, pbmc_small@ident)
#'
#'
PairWiseJaccardSets<- function(ident1, ident2){
        ident1.list<- split(names(ident1), ident1)
        ident2.list<- split(names(ident2), ident2)
        res<- c()
        for (i in seq_along(ident1.list)){
                ind<- purrr::map_dbl(ident2.list, ~JaccardSets(ident1.list[[i]], .x))
                res<- rbind(res, ind)
        }
        rownames(res)<- names(ident1.list)
        return(res)
}
#' Bootstrap for a fully processed Seurat object
#'
#' @param object A fully processed Seurat object.
#' @param n  Number of times you want to bootstrap.
#' @param rate A number between 0 and 1 for subsampling the cells.
#' @param ... Other parameters passed to \code{\link{PreprocessSubsetData}}
#'
#' @return A list of lists containing the ident from the subsetted reclustered
#' seurat objects.
#' @export
#'
#' @examples
#'

# # see https://github.com/satijalab/seurat/issues/457
# # parallelize Seurat functions. The authors decided to go with the future framework.
# scClusterBoot<- function(object, n = 4, workers = 4, rate = 0.8, ...){
#         multicoreParam <- BiocParallel::MulticoreParam(workers = workers)
#         BiocParallel::register(multicoreParam)
#         # the parameter n is not used inside the function
#         GetProcessedSubsetDataCluster<- function(n, ...){
#                 object<- RandomSubsetData(object, rate = rate)
#                 object<- PreprocessSubsetData(object, ...)
#                 return(list(ident = object@ident, pc.sig = object@meta.data$pc.sig))
#         }
#         boot_clusters<- BiocParallel::bplapply(1:n, GetProcessedSubsetDataCluster)
#         return(boot_clusters)
# }



# scClusterBoot<- function(object, n = 4, workers = 4, rate = 0.8, ...){
#         future::plan(multiprocess)
#         # the parameter n is not used inside the function
#         GetProcessedSubsetDataCluster<- function(n, ...){
#                 object<- RandomSubsetData(object, rate = rate)
#                 object<- PreprocessSubsetData(object, ...)
#                 return(list(ident = object@ident, pc.sig = object@meta.data$pc.sig))
#         }
#         boot_clusters<- future.apply::future_lapply(1:n, GetProcessedSubsetDataCluster)
#         return(boot_clusters)
# }

scClusterBoot<- function(object, n = 4, rate = 0.8, ...){
        # the parameter n is not used inside the function
        GetProcessedSubsetDataCluster<- function(n, ...){
                object<- RandomSubsetData(object, rate = rate)
                object<- PreprocessSubsetData(object, ...)
                return(list(ident = object@ident, pc.sig = object@meta.data$pc.sig))
        }
        boot_clusters<- lapply(1:n, GetProcessedSubsetDataCluster)
        return(boot_clusters)
}




