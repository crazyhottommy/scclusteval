

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
#' PairWiseJaccardSets(pbmc@@ident, pbmc_small@@ident)
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


#' Calculate pair-wise overlapping cluster identities for @@ident slots from two Seurat objects
#'
#'Calculate pair-wise overlapping cluster identities for two named factor vector. e.g.
#' seurat_obj1@ident and seurat_obj2@ident
#' @param ident1 a named factor vector. names are the cell names, the values are
#' the cluster id.
#' @param ident2 a named factor vector. names are the cell names, the values are
#' the cluster id.
#'
#' @return A matrix of pairwise number of common cell identities for each cluster.
#' @export
#'
#' @examples
#' PairWiseOverlappingIdents(pbmc@@ident, pbmc_small@@ident)
PairWiseOverlappingIdents<- function(ident1, ident2){
        ident1.list<- split(names(ident1), ident1)
        ident2.list<- split(names(ident2), ident2)
        res<- c()
        for (i in seq_along(ident1.list)){
                ind<- purrr::map_dbl(ident2.list, ~length(intersect(ident1.list[[i]], .x)))
                res<- rbind(res, ind)
        }
        rownames(res)<- names(ident1.list)
        return(res)

}


#' Match two run of cluster ids with highest Jaccard index
#'
#' @param ident1 a named factor vector. names are the cell names, the values are
#' the cluster id.
#' @param ident2 a named factor vector. names are the cell names, the values are
#' the cluster id.
#'
#' @return A tibble with two columns, column 1 is the cluster ids from ident1, column2
#' is the cluster ids from ident2.
#' @export
#'
#' @examples
#'  MatchClusters(pbmc@@ident, pbmc_small@@ident)
MatchClusters<- function(ident1, ident2){
        jaccard_mat<- PairWiseJaccardSets(ident1, ident2)

        get_corresponding_cluster<- function(x){
                id<- which.max(x)
                return(colnames(jaccard_mat)[id])
        }
        matching_ids<- apply(jaccard_mat, 1, get_corresponding_cluster)
        return(tibble::tibble(ident1 = names(matching_ids), ident2 = matching_ids))
}


#' Assign stable cluster
#'
#' @param ident1 The cluster identity from the original full data set.
#' @param idents A list of cluster identity from the subsampled data sets.
#' @param method what function to summarize the jaccard index across all simulations.
#' default is median
#' @param cutoff Cutoff for assigning stable clusters. A jaccard of 0.4 is used
#' for default
#'
#' @return A list containing the raw data for jaccard index for all simulations,
#' TRUE or FALSE of stable cluster for each cluster and a number of stable clusters
#' and percentage of cells in stable clusters.
#' @export
#'
#' @examples
#'
#' data(idents)
#' AssignStableCluster(pbmc@@ident, idents)
AssignStableCluster<- function(ident1, idents, method = median, cutoff = 0.4){
        mat_list<- purrr::map(idents, ~PairWiseJaccardSets(ident1 = ident1, ident2 = .x))
        SelectHighestJaccard<- function(mat){
                apply(mat, 1, max)

        }
        mat_max<- purrr::map(mat_list, SelectHighestJaccard)
        mats<- purrr::reduce(mat_max, dplyr::bind_rows)

        stable_cluster<- mats %>% dplyr::summarise_all(method) %>%
                dplyr::mutate_all(~ifelse(.x > cutoff, T, F)) %>% unlist()
        number_of_stable_cluster<- sum(stable_cluster)

        #calculate number/percentage of cells in stable clusters
        ident1.list<- split(names(ident1), ident1)
        number_of_cells_each_cluster<- purrr::map_int(ident1.list, length)
        percent_cell_in_stable<- sum(number_of_cells_each_cluster[stable_cluster])/sum(number_of_cells_each_cluster)
        return(list(jaccardIndex = mats, stable_cluster = stable_cluster,
                    percent_cell_in_stable = percent_cell_in_stable,
                    number_of_stable_cluster = number_of_stable_cluster))
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




