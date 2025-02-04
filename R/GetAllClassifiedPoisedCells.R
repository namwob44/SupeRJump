#' Get Prefentially Biased Cells
#'
#' @param data_input 
#' @param lineages_to_compare 
#'
#' @return
#' @export
#'
#' @examples
GetAllClassifiedPoisedCells<-function(data_input,lineages_to_compare){
  
  namevector <- paste0("Poised_",lineages_to_compare)
  data_input[ , namevector] <- NA
  
  for(cell_iter in 1:length(as.character(unique(data_input$Cluster.Assignment)))){
    for(iter2 in 1:length(lineages_to_compare)){
      #for(iter in 1:length(lineages_to_compare)){
      
      if(dim(data_input%>%dplyr::filter(Cluster.Assignment==as.character(unique(data_input$Cluster.Assignment))[cell_iter]))[1]>1){
        
        poised_vector<-ClassifyPoisedCellsSingleLineage(data_input%>%dplyr::filter(Cluster.Assignment==as.character(unique(data_input$Cluster.Assignment))[cell_iter]),
                                                        c(lineages_to_compare[iter2]))
        
        print(colnames(poised_vector)[(dim(poised_vector)[2]-2):dim(poised_vector)[2]])
        
        colnames(poised_vector)[dim(poised_vector)[2]]<-paste0("Poised_",as.character(lineages_to_compare[iter2]),"_",as.character(lineages_to_compare[iter2]))
        data_input[rownames(poised_vector),paste0("Poised_",as.character(lineages_to_compare[iter2]))]<-poised_vector[,dim(poised_vector)[2]]
      }
      #}
    }
  }
  return(data_input)
}