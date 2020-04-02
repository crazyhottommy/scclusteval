

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
        res<- matrix(nrow = length(ident1.list), ncol = length(ident2.list),
                     dimnames = list(names(ident1.list), names(ident2.list)))
        for (i in seq_along(ident1.list)){
                res[i, ]<- purrr::map_dbl(ident2.list, ~JaccardSets(ident1.list[[i]], .x))
        }
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
#' @param method what way to summarize the jaccard index across all simulations.
#' to determine a cluster is stable or not. options are "Jaccard_mean", "Jaccard_median" and "Jaccard_percent"
#' @param jaccard_cutoff Cutoff of the jaccard index to determin a cluster is stable or not.
#' it is the mean or median cutoff when the method is "jaccard_mean" or "jaccard_median" and it is
#' the cutoff for every subsampling when the method is "jaccard_percent"
#' @param percent_cutoff The percentage of jaccard index greater than jaccard_cutoff. Used
#' when method is "jaccard_percent". specify 0.6 when you mean 60%.
#'
#' @return A list containing the raw data for jaccard index for all simulations,
#' TRUE or FALSE of stable cluster for each cluster and a number of stable clusters.
#' A cluster is deemed as stable if the median (or mean) jaccard index is > cutoff.
#' in addtion, a stable_index is calculated, which is the pecentage of jaccard index >
#' cutoff for all the subsampling. e.g. for 100 times subsampling, 0.8 means 80% of the
#' time, the jaccard index is > cutoff. Sometimes, we see bimodal distrbution of the
#' 100 jaccard index, the percentage is a better measurement than the mean or median of the
#' 100 jaccard index.
#'
#' @export
#'
#' @examples
#'
#' data(idents)
#'
#' AssignStableCluster(idents, idents)
#'
AssignStableCluster<- function(idents1, idents2,
                               method = "jaccard_median",
                               jaccard_cutoff = 0.6,
                               percent_cutoff = 0.6){
        mats<- AssignHighestJaccard(idents1, idents2)

        stable_index<- (mats > jaccard_cutoff) %>%
                as.data.frame() %>%
                dplyr::summarise_all(mean) %>%
                unlist()

        if (method == "jaccard_mean"){
                stable_cluster<- mats %>%
                        dplyr::summarise_all(mean) %>%
                        dplyr::mutate_all(~ifelse(.x > jaccard_cutoff, TRUE, FALSE)) %>%
                        unlist()
                number_of_stable_cluster<- sum(stable_cluster)

        } else if (method == "jaccard_median"){
                stable_cluster<- mats %>%
                        dplyr::summarise_all(median) %>%
                        dplyr::mutate_all(~ifelse(.x > jaccard_cutoff, TRUE, FALSE)) %>%
                        unlist()
                number_of_stable_cluster<- sum(stable_cluster)
        } else if (method == "jaccard_percent"){
                number_of_stable_cluster<- sum(stable_index > percent_cutoff)
                stable_cluster<- stable_index > percent_cutoff

        } else {
                stop("please specify jaccard_mean, jaccard_median or jaccard_percent
                     for method")
        }

        return(list(jaccardIndex = mats, stable_cluster = stable_cluster,
                    number_of_stable_cluster = number_of_stable_cluster,
                    stable_index = stable_index))
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




