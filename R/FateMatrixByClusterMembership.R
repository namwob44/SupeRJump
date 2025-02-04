
#' Fate Matrix scaled by cell clusters
#'
#' @param TPM_matrix 
#' @param memberships_by_cluster 
#'
#' @return
#' @export
#'
#' @examples
FateMatrixByClusterMembership<-function(TPM_matrix,memberships_by_cluster){
  
  TPM_matrix<-ifelse(is.nan(TPM_matrix),0,TPM_matrix)
  TPM_matrix<-ifelse(is.na(TPM_matrix),0,TPM_matrix)
  
  temp_TPM<-TPM_matrix
  temp_TPM<-temp_TPM/rowSums(temp_TPM)
  
  ill_condition_check<-ifelse(TPM_matrix!=0,1,0)
  print(which(colSums(ill_condition_check)<=1))
  if(length(which(colSums(ill_condition_check)<=1))>0){
    TPM_matrix<-TPM_matrix[-which(rownames((TPM_matrix))%in%colnames(ill_condition_check)[which(colSums(ill_condition_check)<=1)]),-which(colnames((TPM_matrix))%in%colnames(ill_condition_check)[which(colSums(ill_condition_check)<=1)])]
  }
  Qmat<-TPM_matrix
  
  Q_solved<-eigenMatInverse(diag(nrow=dim(Qmat)[1])-Qmat)#this shit goes vroom if it average number of visits to states
  
  memberships_by_cluster<-memberships_by_cluster[rownames(Qmat),]
  
  Cell_To_clust_commitment<-eigenMapMatMult(Q_solved,memberships_by_cluster)
  rownames(Cell_To_clust_commitment)<-rownames(Qmat)
  colnames(Cell_To_clust_commitment)<-colnames(memberships_by_cluster)
  return(Cell_To_clust_commitment)
  
}

