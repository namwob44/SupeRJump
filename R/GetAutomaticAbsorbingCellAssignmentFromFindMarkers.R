#' Get absorbing states through data-driven means
#'
#' @param sigma_cells 
#' @param sce 
#' @param markers_df 
#'
#' @return
#' @export
#'
#' @examples
GetAutomaticAbsorbingCellAssignmentFromFindMarkers<-function(sigma_cells,sce,markers_df){
  
  unique_cell_types<-levels(sigma_cells$group)
  #as.data.frame%>%
  #group_by(group)%>%
  #pull(group)%>%unique%>%as.character
  temp_matrix<-matrix(data = NA_real_,nrow=ncol(sce),ncol = length(unique_cell_types))
  colnames(temp_matrix)<-unique_cell_types
  rownames(temp_matrix)<-colnames(sce)
  for(iter in 1:length(unique_cell_types)){
    supervised_list_of_list<-markers_df%>%
      dplyr::filter(p_val_adj<0.01)%>%
      dplyr::filter(avg_log2FC>1)%>%
      dplyr::filter(cluster==unique_cell_types[iter])%>%
      pull(gene)%>%
      unique
    print(unique_cell_types[iter])
    temp_scores<-Seurat::AddModuleScore(sce,features=
                                          list("score"=c(gsub("(?<=\\b)([a-z])", "\\U\\1", (supervised_list_of_list), perl=TRUE))),
                                        name="trash",
                                        search=TRUE,)  
    temp_matrix[,iter]<-temp_scores@meta.data$trash1-(min(temp_scores@meta.data$trash1))
    
  }
  
  sink_cells<-sigma_cells%>%
    as.data.frame%>%
    tibble::rownames_to_column(var="row_names")%>%
    dplyr::inner_join(temp_matrix%>%
                        as.data.frame%>%
                        tibble::rownames_to_column(var="row_names")%>%
                        tidyr::pivot_longer(cols=!row_names,names_to="Group_scores",values_to="Marker_scores"),by="row_names")%>%
    dplyr::filter(group==Group_scores)%>%
    dplyr::mutate(Marker_scores=as.numeric(Marker_scores))%>%
    group_by(group)%>%
    dplyr::filter((Marker_scores>=stats::quantile(Marker_scores,0.9)))%>%
    dplyr::select(row_names,group)  
  
  return(sink_cells)
}
