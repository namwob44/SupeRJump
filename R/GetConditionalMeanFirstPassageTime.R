
#' Get Conditional Mean First Passage Time in R parallel strategy
#'
#' @param TPM_matrix 
#' @useDynLib SuperJump
#' @return
#' @export
#'
#' @examples
GetConditionalMeanFirstPassageTime<-function(TPM_matrix){
  TPM_matrix<-ifelse(is.nan(TPM_matrix),0,TPM_matrix)
  TPM_matrix<-ifelse(is.na(TPM_matrix),0,TPM_matrix)
  
  temp_TPM<-TPM_matrix
  temp_TPM<-temp_TPM/rowSums(temp_TPM)
  
  doParallel::registerDoParallel(cores = 8)
  TimeMatrix<-foreach::foreach(current_cell_iter=1:(dim(TPM_matrix)[1]),.combine='rbind')%dopar%{
    temp_vec <-matrix(0,nrow=1,ncol=dim(TPM_matrix)[1])
    Qmat<-temp_TPM[-current_cell_iter,-current_cell_iter]
    Rmat<-temp_TPM[-current_cell_iter,current_cell_iter]%>%as.matrix
    # 
    # Q_solved<-eigenMatInverse(diag(nrow=dim(Qmat)[1])-Qmat)#this shit goes vroom.
    # 
    # Bmat<-eigenMapMatMult(Q_solved,Rmat)
    # diag_mat<-diag(c(Bmat))
    # temp_multiply_matrix<-eigenMapMatMult(eigenMapMatMult(eigenMatInverse(diag_mat),Q_solved),diag_mat)
    
    temp_vec[,-current_cell_iter]<-eigenOutputFlightTime(Qmat,Rmat)
    return(temp_vec)
  }
  stopImplicitCluster()
  
  return(TimeMatrix)
  
}
