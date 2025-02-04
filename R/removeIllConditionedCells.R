
#' Check for removing non-connected cells, necessary for Fate and Visitation
#'
#' @param TPM_matrix 
#'
#' @return
#' @export
#'
#' @examples
removeIllConditionedCells<-function(TPM_matrix){
  ill_condition_check<-ifelse(TPM_matrix!=0,1,0)
  print(which(colSums(ill_condition_check)<=1))
  if(length(which(colSums(ill_condition_check)<=1))>0){
    TPM_matrix<-TPM_matrix[-which(rownames((TPM_matrix))%in%colnames(ill_condition_check)[which(colSums(ill_condition_check)<=1)]),
                           -which(colnames((TPM_matrix))%in%colnames(ill_condition_check)[which(colSums(ill_condition_check)<=1)])]
  }
  return(TPM_matrix)
}