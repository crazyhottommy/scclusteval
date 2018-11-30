
#' A wrapper for preprocessing subsetted Seurat object
#'
#' The wrapper does FindVeriableGenes, ScaleData, RunPCA, JackStraw to
#' determine how many PCs to use, ProjectPCA and FindClusters and retrun
#' a fully processed Seurat object. The input subsetted seurat object is
#' supposed to be fully processed as well. So the NormalizeData step is not
#' necessary.
#'
#' @param object A subsetted Seurat object created by RandomSubsetData
#' @param mean.function Function to compute x-axis value (average expression).
#' Default is to take the mean of the detected (i.e. non-zero) values
#' @param dispersion.function Function to compute y-axis value (dispersion).
#' Default is to take the standard deviation of all values
#' @param x.low.cutoff Bottom cutoff on x-axis for identifying variable genes
#' @param x.high.cutoff Top cutoff on x-axis for identifying variable genes
#' @param y.cutoff Bottom cutoff on y-axis for identifying variable genes
#' @param ... any other parameters for FindVariableGenes
#'
#' @return a fully processed Seurat object
#' @export
#'
#' @examples
#' pbmc_small_subset<- RandomSubsetData(pbmc_small, 0.8)
#' pbmc_small_subset_processed<- PreprocessSubsetData(pbmc_small_subset)
#' pbmc_small_subset_processed@meta.data


PreprocessSubsetData<- function(object, mean.function = ExpMean,
                                dispersion.function = logVMR,
                                x.low.cutoff = 0.05,
                                x.high.cutoff = 10,
                                y.cutoff = 0.5,
                                ...){
        object<- FindVariableGenes(object = object, mean.function = mean.function,
                                   dispersion.function = dispersion.function,
                                   x.low.cutoff = x.low.cutoff,
                                   x.high.cutoff = x.high.cutoff,
                                   y.cutoff = y.cutoff)

        object<- ScaleData(object = object, genes.use = object@var.genes,
                           vars.to.regress = c("percent.mito","nUMI"), block.size = 400,
                           min.cells.to.block=3000,
                           display.progress = TRUE, do.par = TRUE, num.cores = 6)

        object<- RunPCA(object = object, pc.genes = object@var.genes,
                        pcs.compute = 100, do.print = TRUE, pcs.print = 1:5,
                        genes.print = 5)

        object<- JackStraw( object = object, num.replicate = 100, num.pc = 85,
                do.par = T)

        object<- ProjectPCA(object = object, do.print = T,pcs.print = 1:5,
                              genes.print = 30,
                              do.center= T,pcs.store=85)

        object <- FindClusters(object = object, reduction.type = "pca",
                                dims.use = 1:50,
                                n.start = 10,
                                nn.eps = 0.5, resolution = 2, print.output = 0,
                                save.SNN = TRUE, force.recalc = TRUE)
        return(object)
}
