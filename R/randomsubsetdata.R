#' Randomly subset (cells) seurat object by a rate
#'
#' @param object Seurat object
#' @param rate a number betwee 0-1 for subsetting
#' @param random.seed set a random seed for sampling, default is 1.
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
RandomSubsetData<- function(object, rate, random.seed = 1, do.clean = TRUE, ...){
        ncells<- nrow(object@meta.data)
        ncells_subsample<- round(ncells * rate)
        set.seed(random.seed)
        idx<- sample(1:ncells, ncells_subsample, replace = FALSE)
        selected_cells<- rownames(object@meta.data)[idx]
        object<- SubsetData(object, cells.use =  selected_cells, do.clean = do.clean, ...)
        return(object)
}



