#' Wrapper function to obtain Transcription Factor Activity
#'
#' @param seurat_obj the main seurat object to analyze.
#' @param organism this is either "mouse" or "human".
#' @param model this is chosen between "TF", "progeny", or "hallmark"
#' @param assay_name This should be kept at "RNA", however, newer models may allow "SCT"
#' @return the same seurat object with a new assay for TFs
#' @export
#' @import decoupleR
#' @import OmnipathR
#' @import tibble
#' @import Seurat
#' @import tidyr
#' @import msigdbr
#'
#' @examples seurat_obj<-TFWrapper(seurat_obj,assay_name="RNA",organism="mouse",model="TF")
TFWrapper<-function(seurat_obj,assay_name="RNA",organism="mouse",model="TF"){
  mat_normed<-as.matrix(Seurat::GetAssayData(seurat_obj, assay =assay_name))#, layer = 'data'))

  if(model=="TF"){
    net <- decoupleR::get_collectri(organism=organism, split_complexes=FALSE)
    acts_norm <- decoupleR::run_ulm(mat=mat_normed, net=net, .source='source', .target='target',
                         .mor='mor', minsize = 5)
    seurat_obj[['tfsulm']] <- acts_norm %>%
      tidyr::pivot_wider(id_cols = 'source', names_from = 'condition',
                         values_from = 'score') %>%
      tibble::column_to_rownames('source') %>%
      Seurat::CreateAssayObject(.)
    Seurat::DefaultAssay(object = seurat_obj) <- "tfsulm"
    seurat_obj <- Seurat::ScaleData(seurat_obj)
    seurat_obj@assays$tfsulm@data <- seurat_obj@assays$tfsulm@scale.data
    return(seurat_obj)
  }
  if(model=="progney"){
    net <- decoupleR::get_progeny(organism = "human", top = 500)
    acts_norm <- decoupleR::run_mlm(mat=mat_normed, net=net, .source='source', .target='target',
                         .mor='weight', minsize = 5)
    seurat_obj[['progeny']] <- acts_norm %>%
      tidyr::pivot_wider(id_cols = 'source', names_from = 'condition',
                         values_from = 'score') %>%
      tibble::column_to_rownames('source') %>%
      Seurat::CreateAssayObject(.)
    Seurat::DefaultAssay(object = seurat_obj) <- "progeny"
    seurat_obj <- Seurat::ScaleData(seurat_obj)
    seurat_obj@assays$progeny@data <- seurat_obj@assays$progeny@scale.data
    return(seurat_obj)
  }
  if(model=="hallmark"){
    net <- msigdbr::msigdbr(species = "human",category = "H")%>%
      dplyr::select(source=gs_name,target=gene_symbol)%>%
      dplyr::mutate(weight=1)%>%
      dplyr::distinct()
    acts_norm <- decoupleR::run_mlm(mat=mat_normed, net=net, .source='source', .target='target',
                         minsize = 5)
    seurat_obj[['hallmark']] <- acts_norm %>%
      tidyr::pivot_wider(id_cols = 'source', names_from = 'condition',
                         values_from = 'score') %>%
      tibble::column_to_rownames('source') %>%
      Seurat::CreateAssayObject(.)
    Seurat::DefaultAssay(object = seurat_obj) <- "hallmark"
    seurat_obj <- Seurat::ScaleData(seurat_obj)
    seurat_obj@assays$hallmark@data <- seurat_obj@assays$hallmark@scale.data
    return(seurat_obj)
  }
}



#' Get Basic multivariate linear model of fates and and TFs
#' @description
#' This model is used to correlate fates with TFs to nominate candidate hypotheses for mechanisms driving the fates.
#'
#' @param seurat_obj The seurat object which contains lineage fates and the tfsulm assay
#' @param lineages_to_use this is a list of the rownames to use in the lineage_fates assay for the model fit.
#'
#' @return it returns a dataframe containing the statistical results of comparison between the lineage_fates and TFs assay
#' @export
#' @import dplyr
#' @import stats
#' @import tibble
#' @import Seurat
#'
#'
#' @examples model_fit_df<-GetFateTFModel(seurat_obj,c("B_cells","Neu","Mono","Ery"))
GetFateTFModel<-function(seurat_obj,lineages_to_use){

  Seurat::DefaultAssay(seurat_obj)<-"tfsulm"
  a <- rownames(seurat_obj)
  temp_df<-t(Seurat::GetAssayData(seurat_obj,assay = "tfsulm"))%>%
    as.data.frame%>%
    tibble::rownames_to_column(var="row_names")%>%
    dplyr::inner_join(t(Seurat::GetAssayData(seurat_obj,assay = "lineage_fates",layer="scale.data"))%>%
                        as.data.frame%>%
                        tibble::rownames_to_column(var="row_names"),
                      by="row_names")

  # check if more data points than TFs or Motifs.
  # keep TFs
  swapped=FALSE
  if(dim(temp_df)[1]<dim(temp_df[,a])[2]){
    print("swapping")
    temp<-lineages_to_use
    lineages_to_use<-a
    a<-temp
    swapped=TRUE
  }

  res_all<-lm(temp_df[,lineages_to_use]%>%as.matrix~temp_df[,colnames(temp_df)[which(colnames(temp_df)%in%a)]]%>%as.matrix)%>%summary()
  names(res_all) <- colnames(temp_df[,lineages_to_use]%>%as.matrix)

  # summary is a list for each condition. Get the info we need:
  res_new <- res_all %>% lapply(X = ., function(fit){

    scores <- as.vector(fit$coefficients[,3][-1])
    pvals <- as.vector(fit$coefficients[,4][-1])
    sources <- colnames(temp_df[,colnames(temp_df)[which(colnames(temp_df)%in%a)]]%>%as.matrix)
    diff_n <- length(sources) - length(scores)
    if (diff_n > 0) {
      print("PANIC")
    }
    tibble::tibble(score=scores, p_value=pvals, source=sources)
  }) %>% dplyr::bind_rows(.id = "condition") %>%
    dplyr::mutate(statistic = "mlm", .before= 1) %>%
    dplyr::select(statistic, source, condition,
                  score, p_value)%>%
    dplyr::mutate(p_value=ifelse(p_value==0,min(p_value[p_value!=0]),p_value))%>%
    dplyr::mutate(sign_score = (-1)*log10(p_value)*sign(score))

  if(swapped==TRUE){
    res_new<-res_new%>%dplyr::select(statistic,condition,source,score,p_value,sign_score)
    colnames(res_new)<-c("statistic","source","condition","score","p_value","sign_score")
  }
  res_new$p_value<-ifelse(res_new$p_value==0,min(res_new$p_value[res_new$p_value!=0]),res_new$p_value)

  return(res_new)
}
