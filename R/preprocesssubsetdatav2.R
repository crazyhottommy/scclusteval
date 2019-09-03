#' A wrapper for preprocessing subsetted Seurat object using ScaleData
#'
#' The wrapper does FindVeriableGenes, ScaleData, RunPCA, JackStraw to
#' determine how many PCs to use, ProjectPCA and FindClusters and retrun
#' a fully processed Seurat object. The input subsetted seurat object is
#' supposed to be fully processed as well. So the NormalizeData step is not
#' necessary.
#'
#' @param object A subsetted Seurat object created by RandomSubsetData
#' @param num.pc number of PCs to calculate in RunPCA, JackStraw and JackStrawPlot
#' step. The optimal PCs for FindClusters will be determined by only significant PCs
#' from JackStrawPlot or if pc.use is set, JackStraw step will be skipped and use pc.use
#' for FindClusters.
#' @param pc.use number of PCs used for FindClusters. if pc.use is set, JackStraw step
#' will be skipped and use pc.use for FindClusters. score.thresh and sig.pc.thresh will be ignored.
#' @param n.start Number of random start.
#' @param nn.eps Error bound when performing nearest neighbor seach using RANN;
#' default of 0.0 implies exact nearest neighbor search. See FindClusters.
#' @param resolution Value of the resolution parameter, use a value above (below)
#' 1.0 if you want to obtain a larger (smaller) number of communities. see FIndClusters.
#' @param k.param Defines k for the k-nearest neighbor algorithm.
#' @param score.thresh Threshold to use for the proportion test of PC significance.
#' @param sig.pc.thresh Threshold for the significance of a particular PC.
#' @param ... any other parameters
#' @param variable.features.n number of variable features for \code{\link[Seurat]{SCTransform}}
#' @param workers number of CPUs to use for \code{\link[future]{plan}} parallel processing
#'
#' @return a fully processed Seurat object
#' @export
#'
#' @examples
#' pbmc_small_subset<- RandomSubsetData(pbmc_small, 0.8)
#' pbmc_small_subset_processed<- PreprocessSubsetData(pbmc_small_subset)
#' pbmc_small_subset_processed@meta.data


PreprocessSubsetDataV2<- function(object,
                                nfeatures = 2000,
                                num.pc = 20,
                                pc.use = NULL,
                                workers = 2,
                                score.thresh = 1e-5,
                                sig.pc.thresh = 0.05,
                                n.start = 100,
                                nn.eps = 0,
                                resolution = 0.8,
                                k.param = 30,
                                ...){
        meta.data.colnames<- object@meta.data %>% colnames()
        vars.to.regress<- c("percent.mt","nFeature_RNA")
        # in case the seurat object does not have percent.mito in metadata
        vars.to.regress<- vars.to.regress[vars.to.regress %in% meta.data.colnames]
        object<- FindVariableFeatures(object, selection.method = "vst", nfeatures = nfeatures)
        object<- ScaleData(object)

        object<- RunPCA(object = object, features = VariableFeatures(object = object),
                        npcs = num.pc)

        if (is.null(pc.use)){
                object<- JackStraw( object = object, num.replicate = 100, dims = num.pc)

                object <- ScoreJackStraw(object = object, dims = 1:num.pc, score.thresh = score.thresh)

                PC_pvalues<- object@reductions$pca@jackstraw@overall.p.values

                ## determin how many PCs to use.
                pc.use<- min(which(PC_pvalues[,"Score"] > sig.pc.thresh)) -1

        }

        # add significant pc number to metadata, need to have names same as the cells
        pc.use.meta<- rep(pc.use, length(colnames(object)))
        names(pc.use.meta)<- colnames(object)
        object<- AddMetaData(object = object, metadata = pc.use.meta, col.name = "pc.use")
        object<- FindNeighbors(object, dims = 1:pc.use, k.param = k.param, nn.eps = nn.eps,
                               verbose = FALSE, reduction = "pca", force.recalc = TRUE)
        object <- FindClusters(object = object, reduction.type = "pca",
                               n.start = n.start,
                               resolution = resolution,
                               verbose = FALSE)
        return(object)
}
