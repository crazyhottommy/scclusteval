

#' Read multiple 10x run into Seurat objects and merge into a single Seurat object
#'
#' Read multiple 10x run into Seurat objects and merge into a single Seurat object.
#' The names of the list of paths will be prepended to the cell name.
#'
#' @param input_folders A named list of folder path for each run.
#' @param do.normalize Whether or not normalize the data after mergeing, default is FALSE
#' @param ... Other parameters for CreatSeuratObject in the Seurat package
#'
#' @return A single merged Seurat object from mulitple 10x runs.
#' @export
#'
#' @examples
#' \dontrun{
#' library(fs)
#' library(here)
#' library(stringr)
#' input_folders<- dir_ls( path = here("data"), recursive = T) %>% path_dir() %>%
#' unique() %>% str_subset("mm10-1.2.0_premrna")
#' merged_seurat<- MergeMultipleSeuratObjects(input_folders)
#' }


MergeMultipleSeuratObjects<- function(input_folders, do.normalize = FALSE, ...){
        seurat_data<- purrr::map(input_folders, Read10X)
        #prefix the sample name to the cell name, otherwise merge seurat objects gives error
        add_sample_name_to_cell<- function(x, y){
                colnames(x)<- paste(y, colnames(x), sep = "_")
                return(x)
        }
        sample_names<- names(input_folders)
        seurat_data<- purrr::map2(seurat_data, sample_names, add_sample_name_to_cell)
        seurat_objects<- purrr::map2(seurat_data, sample_names,
                                     function(x,y) CreateSeuratObject(raw.data = x,
                                                                      project = y,
                                                                      ...))
        #merge to a single seurat object
        merged_seurat<- purrr::reduce(seurat_objects,
                                      function(x,y) {MergeSeurat(x,y,
                                                                 do.normalize = do.normalize)})
}
