
#' An alternative to splitting your TPMs to conditions (deprecated)
#'
#' @param TPM_matrix 
#' @param cosine_data 
#'
#' @return
#' @export
#'
#' @examples
SplitTPMByCondition<-function(TPM_matrix,cosine_data){
  split_TPM_list<-list()
  for(condition_iter in 1:length(as.character(unique(cosine_data$type)))){
    cell_ids<-cosine_data%>%
      filter(type==as.character(unique(cosine_data$type)[condition_iter]))%>%
      tibble::rownames_to_column(var="row_names")%>%
      pull(row_names)
    Reduced_matrix<-TPM_matrix[rownames(TPM_matrix)%in%cell_ids,colnames(TPM_matrix)%in%cell_ids]
    split_TPM_list[[condition_iter]]<-Reduced_matrix
  }
  names(split_TPM_list)<-paste0(as.character(unique(cosine_data$type)),"_TPM")
  return(split_TPM_list)
}