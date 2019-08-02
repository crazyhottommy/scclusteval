#' Randomly subset (cells) seurat object by a rate
#'
#' @param object Seurat object
#' @param rate a number betwee 0-1 for subsetting
#' @param random.subset.seed set a random seed for sampling, default is NULL.
#' @param ... any other parameters to \code{\link[Seurat]{subset}}
#'
#' @return Returns a randomly subsetted seurat object
#' @export
#'
#' @examples
#' pbmc_small
#' pbmc_small_subset<- RandomSubsetData(pbmc_small, 0.8)
#' dim(pbmc_small_subset@@meta.data)
#'
#'
# read this issue https://github.com/satijalab/seurat/issues/243
# Seurat V3 does not have do.clean =T any more
# see https://github.com/satijalab/seurat/issues/1792 use DietSeurat
RandomSubsetData<- function(object, rate, random.subset.seed = NULL, ...){
        ncells<- nrow(object@meta.data)
        ncells.subsample<- round(ncells * rate)

        set.seed(random.subset.seed)

        selected.cells<- sample(colnames(object), ncells.subsample)
        object<- subset(object, cells =  selected.cells,
                            ...)
        return(object)
}



