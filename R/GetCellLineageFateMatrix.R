#' Aggregating the Fates from absorbing states to lineages
#'
#' @param Bmat 
#' @param sink_cells 
#'
#' @return
#' @export
#'
#' @examples
GetCellLineageFateMatrix<-function(Bmat,sink_cells){
  sink_groups<-c(as.character(unique(sink_cells$group)))
  Bmat_condensed<-matrix(data=0,nrow=nrow(Bmat),ncol=length(sink_groups))
  rownames(Bmat_condensed)<-rownames(Bmat)
  colnames(Bmat_condensed)<-sink_groups
  for(sink_iter in 1:length(sink_groups)){
    lineage_columns<-sink_cells%>%
      dplyr::filter(group==sink_groups[sink_iter])%>%
      pull(row_names)
    if(length(lineage_columns)>1){
      Bmat_condensed[,sink_iter]<-rowSums(Bmat[,colnames(Bmat)%in%lineage_columns],na.rm=T)
    }else{
      Bmat_condensed[,sink_iter]<-Bmat[,colnames(Bmat)%in%lineage_columns]
    }
  }
  return(Bmat_condensed)
}
