#' Wrapper function for non-biased and biased cells for ATAC differential expresssion
#'
#' @param suerat_obj 
#' @param poise_df 
#' @param lineages_to_compare 
#' @param foldername 
#'
#' @return
#' @export
#'
#' @examples
RunDifferentialExpressionOnPoiseVsNonPoiseATAC<-function(suerat_obj,poise_df,lineages_to_compare,foldername){
  if(length(as.character(unique(poise_df$Poised)))<2){
    return(NULL)
  }
  if(table(poise_df$Poised)[1]<=2|table(poise_df$Poised)[2]<=2){
    return(NULL)
  }
  if(!dir.exists(foldername)){dir.create(foldername)}
  
  cells_group_1<-rownames(poise_df%>%dplyr::filter(Poised=="poised"))
  cells_group_2<-rownames(poise_df%>%dplyr::filter(Poised=="not"))
  #Idents(suerat_obj,cells=cells_group_1)<-"poised"
  #Idents(suerat_obj,cells=cells_group_2)<-"not"
  
  DefaultAssay(suerat_obj)<-"RNA"
  
  marker_object<-Seurat::FindMarkers(object =suerat_obj@assays[["RNA"]],slot="data",cells.1=cells_group_1,cells.2=cells_group_2)
  
  write.csv(marker_object,file = paste0(foldername,"RNA_",as.character(poise_df$Cluster.Assignment)[1],"_poising_for_",lineages_to_compare,".csv"),sep = ",",append = T)
  DefaultAssay(suerat_obj)<-"chromvar"
  marker_object<-Seurat::FindMarkers(object =suerat_obj@assays[["chromvar"]],cells.1=cells_group_1,cells.2=cells_group_2)
  rownames(marker_object)<-suerat_obj@assays[["ATAC"]]@motifs@motif.names[match(rownames(marker_object),names(suerat_obj@assays[["ATAC"]]@motifs@motif.names))]%>%as.matrix
  
  write.csv(marker_object,file = paste0(foldername,"Motif_",as.character(poise_df$Cluster.Assignment)[1],"_poising_for_",lineages_to_compare,".csv"),sep = ",",append = T)
  return(TRUE)
}
