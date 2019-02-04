

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



#' Assign highest Jaccard index for each cluster of the subsampled data set before
#' reclustering with the cluster identites of subsampled data set after reclustering
#'
#' @param idents1 A list of cluster identity copied from the orginal data sets.
#' idents1 is a list of the cluster identity from the subsampled data sets before reclustering.
#' @param idents2 A list of cluster identity from the subsampled data sets.
#' idents2 is a list of the cluster identity from the subsampled data sets after reclustering.
#' The order of identities in idents1 and idents2 should correspond to each other.
#'
#' @return A matrix with dimention of #number of subsampling * #number of clusters in the
#' original data set.
#' @export
#'
#' @examples
AssignHighestJaccard<- function(idents1, idents2){
        mat_list<- purrr::map2(idents1, idents2,  ~PairWiseJaccardSets(ident1 = .x, ident2 = .y))
        SelectHighestJaccard<- function(mat){
                apply(mat, 1, max)

        }
        # or use the anonymous function
        mat_max<- purrr::map(mat_list, SelectHighestJaccard)
        mats<- purrr::reduce(mat_max, dplyr::bind_rows)
        return(mats)
}

#' Assign stable cluster
#'
#' @param idents1 A list of cluster identity copied from the orginal data sets.
#' idents1 is a list of the cluster identity from the subsampled data sets before reclustering.
#' @param idents2 A list of cluster identity from the subsampled data sets.
#' idents2 is a list of the cluster identity from the subsampled data sets after reclustering.
#' The order of identities in idents1 and idents2 should correspond to each other.
#' @param method what function to summarize the jaccard index across all simulations.
#' default is median
#' @param cutoff Cutoff for assigning stable clusters. A jaccard of 0.4 is used
#' for default
#'
#' @return A list containing the raw data for jaccard index for all simulations,
#' TRUE or FALSE of stable cluster for each cluster and a number of stable clusters.
#'
#' @export
#'
#' @examples
#'
#' data(idents)
#' dummy example, all clusters are stable.
#' AssignStableCluster(idents, idents)
#'
AssignStableCluster<- function(idents1, idents2, method = median, cutoff = 0.6){
        mats<- AssignHighestJaccard(idents1, idents2)
        stable_cluster<- mats %>% dplyr::summarise_all(method) %>%
                dplyr::mutate_all(~ifelse(.x > cutoff, T, F)) %>% unlist()
        number_of_stable_cluster<- sum(stable_cluster)
        return(list(jaccardIndex = mats, stable_cluster = stable_cluster,
                    number_of_stable_cluster = number_of_stable_cluster))
}


#' Calculate the percentage of cells in stable clusters in the full data set
#'
#' @param ident. A named factor vector. names are the cell names, the values are
#' the cluster id from the full data set.
#' @param stable_cluster. A logical vector for each of the original cluster indicating
#' it is stable or not, calculated from \code{\link{AssignStableCluster}}
#'
#' @return A percentage of cells in stable cluster
#' @export
#'
#' @examples

CalculatePercentCellInStable<- function(ident, stable_cluster){
        ident.list<- split(names(ident), ident)
        number_of_cells_each_cluster<- purrr::map_int(ident.list, length)
        percent_cell_in_stable<- sum(number_of_cells_each_cluster[stable_cluster])/sum(number_of_cells_each_cluster)
        return(percent_cell_in_stable)

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




