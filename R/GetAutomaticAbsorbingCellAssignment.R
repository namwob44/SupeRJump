#' Get a model based approach to determine Absorbing States
#'
#' @param sigma_cells 
#' @param sink_cell_names 
#' @param sink_cell_nodes 
#' @param get_AVG 
#'
#' @return
#' @export
#'
#' @examples
GetAutomaticAbsorbingCellAssignment<-function(sigma_cells,sink_cell_names=NULL,sink_cell_nodes=NULL,get_AVG=F){
  
  if(is.null(sink_cell_names)){
    if(is.null(sink_cell_nodes)){
      # Automatic way to get sink cell names for all groups
      if(get_AVG){
        sink_cells<-sigma_cells%>%
          as.data.frame%>%
          tibble::rownames_to_column(var="row_names")%>%
          mutate(PT_score=as.numeric(PT_score))%>%
          group_by(group)%>%
          dplyr::filter(abs(PT_score-quantile(PT_score,0.5))<=1e-3)%>%
          dplyr::select(row_names,group)
        
      }else{
        sink_cells<-sigma_cells%>%
          as.data.frame%>%
          tibble::rownames_to_column(var="row_names")%>%
          mutate(PT_score=as.numeric(PT_score))%>%
          group_by(group)%>%
          dplyr::filter((PT_score>=stats::quantile(PT_score,0.9)))%>%
          dplyr::select(row_names,group)  
      }
      
      
    }else{
      # Automatic way to get sink cell names for SPECIFIC groups
      sink_cells<-sigma_cells%>%
        as.data.frame%>%
        tibble::rownames_to_column(var="row_names")%>%
        mutate(PT_score=as.numeric(PT_score))%>%
        dplyr::filter(group%in%sink_cell_nodes)%>%
        group_by(group)%>%
        #dplyr::filter((PT_score>=stats::quantile(PT_score,0.9)))%>%
        dplyr::select(row_names,group)
    }
  }else{
    # They gave us cell names to use.
    sink_cells<-sigma_cells%>%
      as.data.frame%>%
      tibble::rownames_to_column(var="row_names")%>%
      mutate(PT_score=as.numeric(PT_score))%>%
      dplyr::filter(row_names%in%sink_cell_names)%>%
      dplyr::select(row_names,group)
  }
  return(sink_cells)
}
