#' Calculate Fate and Fundamental Matrices
#'
#' @param TPM_matrix 
#' @param sink_cells 
#' @useDynLib SuperJump
#' @return
#' @export
#'
#' @examples
GetFateMatrix<-function(TPM_matrix,sink_cells){
  TPM_matrix<-ifelse(is.na(TPM_matrix),0,TPM_matrix)
  TPM_matrix<-ifelse(is.nan(TPM_matrix),0,TPM_matrix)
  
  TPM_matrix<-TPM_matrix/rowSums(TPM_matrix)
  # check columns have more than 1 entry
  TPM_matrix<-removeIllConditionedCells(TPM_matrix)
  
  Rmat<-(TPM_matrix)[-which(rownames((TPM_matrix))%in%sink_cells$row_names),which(colnames((TPM_matrix))%in%sink_cells$row_names)]%>%as.matrix # need this matrix conversion to handle single absorbing state
  Qmat<-(TPM_matrix)[-which(rownames((TPM_matrix))%in%sink_cells$row_names),-which(colnames((TPM_matrix))%in%sink_cells$row_names)]
  
  Q_solved<-eigenMatInverse(diag(nrow=dim(Qmat)[1])-Qmat)#this shit goes vroom.
  colnames(Rmat)<-colnames(TPM_matrix)[which(colnames((TPM_matrix))%in%sink_cells$row_names)] # need this in case they only have a single absorbing state
  Bmat<-matrix(data=0,nrow=nrow((TPM_matrix)),ncol=dim(Rmat)[2])
  rownames(Bmat)<-rownames((TPM_matrix))
  colnames(Bmat)<-colnames(Rmat)
  Bmat[-which(rownames((TPM_matrix))%in%sink_cells$row_names),]<-eigenMapMatMult(Q_solved,Rmat)
  if(ncol(Bmat)>1){
    diag(Bmat[which(rownames((TPM_matrix))%in%sink_cells$row_names),])<-1
  }
  rownames(Q_solved)<-colnames(Q_solved)<-rownames(Qmat)#rownames(TPM_matrix)[-which(rownames((TPM_matrix))%in%sink_cells$row_names)]
  return(list("Bmat"=Bmat,
              "Nmat"=Q_solved))
}
