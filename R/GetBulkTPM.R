#' Crude way to get bulk transition probablity matrix (deprecated)
#'
#' @param TPM_matrix 
#' @param cosine_data 
#'
#' @return
#' @export
#'
#' @examples
GetBulkTPM<-function(TPM_matrix,cosine_data){
  reduced_data<-cosine_data[rownames(TPM_matrix),]
  bulk_TPM_blocks<-matrix(NA_real_,nrow=length(unique(reduced_data$group)),ncol=length(unique(reduced_data$group)))
  rownames(bulk_TPM_blocks)<-colnames(bulk_TPM_blocks)<-as.character(unique(reduced_data$group))
  for(sink_iter in 1:length(unique(reduced_data$group))){
    for(sink_other_iter in 1:length(unique(reduced_data$group))){
      row_cells<-reduced_data%>%
        filter(group==as.character(unique(reduced_data$group)[sink_iter]))%>%
        tibble::rownames_to_column(var="row_names")%>%
        pull(row_names)
      col_cells<-reduced_data%>%
        filter(group==as.character(unique(reduced_data$group)[sink_other_iter]))%>%
        tibble::rownames_to_column(var="row_names")%>%
        pull(row_names)
      
      temp_block<-TPM_matrix[rownames(TPM_matrix)%in%row_cells,
                             colnames(TPM_matrix)%in%col_cells]
      prob_check <-ifelse(temp_block==0,NA,temp_block)
      bulk_TPM_blocks[sink_iter,sink_other_iter]<-mean(rowMeans(temp_block,na.rm = T),na.rm=T)*(table(is.na(prob_check))[1]/sum(table(is.na(prob_check))))
    }
  }
  return(bulk_TPM_blocks)
}
