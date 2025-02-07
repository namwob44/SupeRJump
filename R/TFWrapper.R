#' Wrapper function to obtain Transcription Factor Activity
#'
#' @param suerat_obj 
#' @param organism 
#' @param model 
#' @param assay_name
#' @return
#' @export
#'
#' @examples
TFWrapper<-function(suerat_obj,assay_name="RNA",organism="mouse",model="TF"){
  mat_normed<-as.matrix(Seurat::GetAssayData(seurat_obj_subset, assay =assay_name, layer = 'data'))

  #rownames(mat_normed)<-rownames(suerat_obj@assays$RNA)
  #colnames(mat_normed)<-colnames(suerat_obj@assays$RNA)
  if(model=="TF"){
    net <- get_collectri(organism=organism, split_complexes=FALSE)
    acts_norm <- run_ulm(mat=mat_normed, net=net, .source='source', .target='target',
                         .mor='mor', minsize = 5)
    suerat_obj[['tfsulm']] <- acts_norm %>%
      tidyr::pivot_wider(id_cols = 'source', names_from = 'condition',
                         values_from = 'score') %>%
      tibble::column_to_rownames('source') %>%
      Seurat::CreateAssayObject(.)
    DefaultAssay(object = suerat_obj) <- "tfsulm"
    suerat_obj <- Seurat::ScaleData(suerat_obj)
    suerat_obj@assays$tfsulm@data <- suerat_obj@assays$tfsulm@scale.data
    return(suerat_obj)
  }
  if(model=="progney"){
    net <- decoupleR::get_progeny(organism = "human", top = 500)
    acts_norm <- run_mlm(mat=mat_normed, net=net, .source='source', .target='target',
                         .mor='weight', minsize = 5)
    suerat_obj[['progeny']] <- acts_norm %>%
      tidyr::pivot_wider(id_cols = 'source', names_from = 'condition',
                         values_from = 'score') %>%
      tibble::column_to_rownames('source') %>%
      Seurat::CreateAssayObject(.)
    DefaultAssay(object = suerat_obj) <- "progeny"
    suerat_obj <- Seurat::ScaleData(suerat_obj)
    suerat_obj@assays$progeny@data <- suerat_obj@assays$progeny@scale.data
    return(suerat_obj)
  }
  if(model=="hallmark"){
    net <- msigdbr::msigdbr(species = "human",category = "H")%>%
      select(source=gs_name,target=gene_symbol)%>%
      mutate(weight=1)%>%distinct()
    acts_norm <- run_mlm(mat=mat_normed, net=net, .source='source', .target='target',
                         minsize = 5)
    suerat_obj[['hallmark']] <- acts_norm %>%
      tidyr::pivot_wider(id_cols = 'source', names_from = 'condition',
                         values_from = 'score') %>%
      tibble::column_to_rownames('source') %>%
      Seurat::CreateAssayObject(.)
    DefaultAssay(object = suerat_obj) <- "hallmark"
    suerat_obj <- Seurat::ScaleData(suerat_obj)
    suerat_obj@assays$hallmark@data <- suerat_obj@assays$hallmark@scale.data
    return(suerat_obj)
  }
}
