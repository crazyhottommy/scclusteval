#' Randomly subset (cells) seurat object by a rate
#'
#' @param object Seurat object
#' @param rate a number betwee 0-1 for subsetting
#' @param random.subset.seed set a random seed for sampling, default is NULL.
#' @param do.clean subset the seurat object in a clean way, include the raw.data
#' @param ... any other parameters to SubsetData from Seurat package
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
# I have to do set.seed(NULL) to make it truly random.
RandomSubsetData<- function(object, rate, random.subset.seed = NULL, do.clean = TRUE, ...){
        ncells<- nrow(object@meta.data)
        ncells.subsample<- round(ncells * rate)

        set.seed(random.subset.seed)

        selected.cells<- sample(object@cell.names, ncells.subsample)
        object<- SubsetData(object, cells.use =  selected.cells, do.clean = do.clean,
                            ...)
        return(object)
}



