#' Individual lineage and cluster for preferentially biasing
#'
#' @param data_subset 
#' @param lineages_to_compare 
#'
#' @return
#' @export
#'
#' @examples
ClassifyPoisedCellsSingleLineage<-function(data_subset,lineages_to_compare){
  
  left_ecdf_func <- function(data) { 
    Length <- length(data) 
    sorted <- sort(data) 
    ecdf <- rep(0, Length) 
    for (i in 1:Length) { 
      ecdf[i] <- sum(sorted <= data[i]) / Length 
    } 
    return(ecdf) 
  } 
  right_ecdf_func <- function(data) { 
    Length <- length(data) 
    sorted <- sort(data) 
    ecdf <- rep(0, Length) 
    for (i in 1:Length) { 
      ecdf[i] <- sum(sorted >= data[i]) / Length 
    } 
    return(ecdf) 
  } 
  
  measures=data_subset[,lineages_to_compare]%>%as.matrix
  outlier_score<-matrix(NA,nrow=dim(measures)[1],ncol = dim(measures)[2])
  skew_meas <- apply(measures,2,e1071::skewness)
  skew_meas[is.nan(skew_meas)]<-0
  for(pro_iter in 1:dim(measures)[2]){
    right_ecdf<-right_ecdf_func(measures[,pro_iter])
    left_ecdf<-left_ecdf_func(measures[,pro_iter])
    for(cell_iter in 1:dim(measures)[1]){
      
      #if(skew_meas[pro_iter]>0){
      outlier_score[cell_iter,pro_iter]<-right_ecdf[cell_iter]
      #}
      #else{
      #  outlier_score[cell_iter,pro_iter]<-left_ecdf[cell_iter]
      #}
    }
  }
  outlier_score_neglog<-(-1)*log(outlier_score)%*%(abs(skew_meas)/norm(skew_meas,type="2"))
  outlier_cutoff_val<-quantile(outlier_score_neglog)[4]+1.5*IQR(outlier_score_neglog)
  
  data_subset$Outlier_Score<-outlier_score_neglog
  data_subset$Poised<-ifelse(outlier_score_neglog>=outlier_cutoff_val,"poised","not")
  return(data_subset)
}
