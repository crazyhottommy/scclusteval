
#' A wrapper for preprocessing subsetted Seurat object
#'
#' The wrapper does FindVeriableGenes, ScaleData, RunPCA, JackStraw to
#' determine how many PCs to use, ProjectPCA and FindClusters and retrun
#' a fully processed Seurat object. The input subsetted seurat object is
#' supposed to be fully processed as well. So the NormalizeData step is not
#' necessary.
#'
#' @param object A subsetted Seurat object created by RandomSubsetData
#' @param x.low.cutoff Bottom cutoff on x-axis for identifying variable genes
#' @param x.high.cutoff Top cutoff on x-axis for identifying variable genes
#' @param y.cutoff Bottom cutoff on y-axis for identifying variable genes
#' @param num.pc number of PCs to calculate in RunPCA, JackStraw and JackStrawPlot
#' step. The optimal PCs for FindClusters will be determined by only significant PCs
#' from JackStrawPlot or if pc.use is set, JackStraw step will be skipped and use pc.use
#' for FindClusters.
#' @param pc.use number of PCs used for FindClusters. if pc.use is set, JackStraw step
#' will be skipped and use pc.use for FindClusters. score.thresh and sig.pc.thresh will be ignored.
#' @param do.par use parallel processing for regressing out variables faster.
#' If set to TRUE, will use half of the machines available cores. see JackStraw.
#' @param num.cores If do.par = TRUE, specify the number of cores to use.
#' Note that for higher number of cores, larger free memory is needed.
#' If num.cores = 1 and do.par = TRUE, num.cores will be set to half of all
#' available cores on the machine. see JackStraw.
#' @param n.start Number of random start.
#' @param nn.eps Error bound when performing nearest neighbor seach using RANN;
#' default of 0.0 implies exact nearest neighbor search. See FindClusters.
#' @param resolution Value of the resolution parameter, use a value above (below)
#' 1.0 if you want to obtain a larger (smaller) number of communities. see FIndClusters.
#' @param k.param Defines k for the k-nearest neighbor algorithm.
#' @param score.thresh Threshold to use for the proportion test of PC significance.
#' @param sig.pc.thresh Threshold for the significance of a particular PC.
#' @param ... any other parameters for FindVariableGenes
#'
#' @return a fully processed Seurat object
#' @export
#'
#' @examples
#' pbmc_small_subset<- RandomSubsetData(pbmc_small, 0.8)
#' pbmc_small_subset_processed<- PreprocessSubsetData(pbmc_small_subset)
#' pbmc_small_subset_processed@meta.data


PreprocessSubsetData<- function(object,
                                x.low.cutoff = 0.05,
                                x.high.cutoff = 10,
                                y.cutoff = 0.5,
                                num.pc = 20,
                                pc.use = NULL,
                                do.par =TRUE,
                                num.cores = 2,
                                score.thresh = 1e-5,
                                sig.pc.thresh = 0.05,
                                n.start = 100,
                                nn.eps = 0,
                                resolution = 0.8,
                                k.param = 30,
                                ...){
        object<- FindVariableGenes(object = object,
                                   x.low.cutoff = x.low.cutoff,
                                   x.high.cutoff = x.high.cutoff,
                                   y.cutoff = y.cutoff)
        meta.data.colnames<- object@meta.data %>% colnames()
        vars.to.regress<- c("percent.mito","nUMI")
        # in case the seurat object does not have percent.mito in metadata
        vars.to.regress<- vars.to.regress[vars.to.regress %in% meta.data.colnames]
        # default is on variable features only, omit the features argument
        object<- ScaleData(object = object,
                           vars.to.regress = vars.to.regress, block.size = 1000,
                           min.cells.to.block=3000,
                           display.progress = TRUE, do.par = TRUE, num.cores = num.cores)

        object<- RunPCA(object = object, features = VariableFeatures(object = object),
                        npcs = num.pc, do.print = FALSE)

        if (is.null(pc.use)){
                object<- JackStraw( object = object, num.replicate = 100, num.cores = num.cores,
                                    do.par = T, dims = num.pc)

                object <- JackStrawPlot(object = object, dims = 1:num.pc, score.thresh = score.thresh)

                PC_pvalues<- object@dr$pca@jackstraw@overall.p.values

                ## determin how many PCs to use.
                pc.use<- min(which(PC_pvalues[,"Score"] > sig.pc.thresh)) -1

        }

        # add significant pc number to metadata, need to have names same as the cells
        pc.use.meta<- rep(pc.use, length(object@cell.names))
        names(pc.use.meta)<- object@cell.names
        object<- AddMetaData(object = object, metadata = pc.use.meta, col.name = "pc.use")
        object <- FindClusters(object = object, reduction.type = "pca",
                                dims.use = 1:pc.use,
                                n.start = n.start,
                                k.param = k.param,
                                nn.eps = nn.eps, resolution = resolution,
                                print.output = FALSE,
                                save.SNN = TRUE, force.recalc = TRUE)
        return(object)
}
