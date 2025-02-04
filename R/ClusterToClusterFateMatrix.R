
#' Turning cell by lineage matrix, to a cluster by lineage matrix
#'
#' @param Bmat 
#' @param cosine_data 
#' @param sink_cells 
#'
#' @return
#' @export
#'
#' @examples
ClusterToClusterFateMatrix<-function(Bmat,cosine_data,sink_cells){
  
  
  # WE. CAN NORMALIZE AND DETERMINE CUSTOM THINGS HERE!!
  sink_groups <-colnames(Bmat)
  Fate_matrix <-matrix(data=NA,nrow=length(unique(cosine_data$group)),ncol=length(sink_groups))
  rownames(Fate_matrix)<-as.character(unique(cosine_data$group))
  colnames(Fate_matrix)<-colnames(Bmat)
  
  for(current_lineage_iter in 1:length(unique(cosine_data$group))){
    current_source_cells<-cosine_data%>%
      dplyr::filter(group==as.character(unique(cosine_data$group)[current_lineage_iter]))%>%
      tibble::rownames_to_column(var="row_names")%>%
      pull(row_names)
    for(sink_lineage_iter in 1:length(sink_groups)){
      Fate_matrix[current_lineage_iter,sink_lineage_iter]<-sum(Bmat[rownames(Bmat)%in%current_source_cells,sink_lineage_iter])
    }
  }
  Fate_matrix_norm<-Fate_matrix#apply(Fate_matrix,2,function(x){x/sum(x)})
  return(Fate_matrix_norm)
}
