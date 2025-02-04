
#' A wrapper for finding differential expression on biased cells for ATAC data
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
RunDifferentialExpressionPoiseVsPoiseATAC<-function(suerat_obj,poise_df,lineages_to_compare,foldername){
  
  if(length(as.character(unique(poise_df$Group)))<2){
    print("We up in here")
    return(NULL)
  }
  if(table(poise_df$Group)[1]<=2|table(poise_df$Group)[2]<=2){
    return(NULL)
  }
  if(!dir.exists(foldername)){dir.create(foldername)}
  
  cells_group_1<-rownames(poise_df%>%dplyr::filter(Group==as.character(unique(poise_df$Group)[1])))
  cells_group_2<-rownames(poise_df%>%dplyr::filter(Group==as.character(unique(poise_df$Group)[2])))
  DefaultAssay(suerat_obj)<-"RNA"
  Idents(suerat_obj)<-"Group"
  marker_object<-Seurat::FindMarkers(object =suerat_obj@assays[["RNA"]],slot="data",cells.1=cells_group_1,cells.2=cells_group_2)
  
  write.csv(marker_object%>%as.matrix,file = paste0(foldername,"RNA_",as.character(poise_df$Cluster.Assignment)[1],"_poised_for_",lineages_to_compare,".csv"),sep = ",",append = T)
  DefaultAssay(suerat_obj)<-"chromvar"
  Idents(suerat_obj)<-"Group"
  
  
  marker_object<-Seurat::FindMarkers(object =suerat_obj@assays[["chromvar"]],slot="data",cells.1=cells_group_1,cells.2=cells_group_2)
  rownames(marker_object)<-suerat_obj@assays[["ATAC"]]@motifs@motif.names[match(rownames(marker_object),names(suerat_obj@assays[["ATAC"]]@motifs@motif.names))]%>%as.matrix
  
  write.csv(marker_object%>%as.matrix,file = paste0(foldername,"Motif_",as.character(poise_df$Cluster.Assignment)[1],"_poised_for_",lineages_to_compare,".csv"),sep = ",",append = T)
  return(TRUE)
}
