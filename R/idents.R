#' cluster identity of subsetted pbmc data
#'
#' The 2700 cell pbmc data were subsetted to 80 percent of the cells for 100 times.
#' Each time, we fully re-processed the subsetted data from FindVaraiableGenes to
#' FindClusters using k =30 and resolution = 0.6, and record the cluster identity
#' from the processed seurat@@ident
#' slot and saved in to a list of factor.
#'
#' @docType data
#'
#' @usage data(idents)
#'
#' @format A list of factors
#' @source \url{https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz}
#'
"idents"
