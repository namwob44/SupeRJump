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
  
  # let's create placeholder columns of interest, along with the outlier scores
  namevector <- paste0("Biased_",lineages_to_compare)
  data_input[ , namevector] <- NA
  namevector <- paste0("Outlier_Score_",lineages_to_compare)
  data_input[ , namevector] <- NA
  namevector <- paste0("Outlier_Score_loss_",lineages_to_compare)
  data_input[ , namevector] <- NA

  for(cluster_iter in 1:length(as.character(unique(data_input$Cluster.Assignment)))){
    for(lineage_iter in 1:length(lineages_to_compare)){
      
        current_group=as.character(unique(data_input$Cluster.Assignment))[cluster_iter]
        poised_vector<-ClassifyPoisedCellsSingleLineage(data_subset = data_input%>%
                                                      filter(Cluster.Assignment==as.character(unique(data_input$Cluster.Assignment))[cluster_iter]),
                                                    c(lineages_to_compare[lineage_iter],
                                                      paste0(current_group,"_membership"),
                                                      paste0(lineages_to_compare[lineage_iter],"_membership")))
      colnames(poised_vector)
      data_input[rownames(poised_vector),paste0("Biased_",as.character(lineages_to_compare[lineage_iter]))]<-poised_vector[,"Biased"]
      data_input[rownames(poised_vector),paste0("Outlier_Score_loss_",as.character(lineages_to_compare[lineage_iter]))]<-poised_vector[,"Outlier_Score_losses"]
      data_input[rownames(poised_vector),paste0("Outlier_Score_",as.character(lineages_to_compare[lineage_iter]))]<-poised_vector[,"Outlier_Score"]
          
    }
  }
  return(data_input)
}
