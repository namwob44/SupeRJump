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
       
      if(pro_iter==1){ # this get's lineage while other 2 do membership
        outlier_score[cell_iter,pro_iter]<-right_ecdf[cell_iter]
      }
      else if(skew_meas[pro_iter]>0){
        outlier_score[cell_iter,pro_iter]<-right_ecdf[cell_iter]
      }
      else{
        outlier_score[cell_iter,pro_iter]<-left_ecdf[cell_iter]
      }
    }
  }
  # this makes sure the lineage we care about is not over dominated by membership.
  skew_meas<-ifelse(abs(skew_meas[1])<=abs(skew_meas),sign(skew_meas)*skew_meas[1],skew_meas) 
  # this heuristic does the skewing and combination of measures.
  outlier_score_neglog<-(-1)*log(outlier_score)%*%(abs(skew_meas)/norm(skew_meas,type="2"))

  # ok let's now get the distribution of outlier scores and figure out how they are robust 
  mad_value <- mad(outlier_score_neglog)

  # Set delta as a multiple of MAD (you can adjust the multiplier as needed)
  delta <- 1.5 * mad_value

  # Huber loss function
  huber_loss <- function(x, delta) {
    abs_x <- abs(x)
    ifelse(abs_x <= delta, 0.5 * (x^2), delta * (abs_x - 0.5 * delta))
  }

  # Calculate Huber loss for each data point
  losses <- huber_loss(outlier_score_neglog - median(outlier_score_neglog), delta)

  # Identify outliers based on a threshold (e.g., 95th percentile)
  #outlier_cutoff_val <- quantile(losses, 0.95)
  outlier_cutoff_val<-quantile(losses)[4]+1.5*IQR(losses)

  #outliers <- outlier_score_neglog[losses > outlier_cutoff_val]
  
  #outlier_score_neglog<-(-1)*log(outlier_score)%*%(abs(skew_meas)/norm(skew_meas,type="2"))
  #outlier_cutoff_val<-quantile(outlier_score_neglog)[4]+1.5*IQR(outlier_score_neglog)
  
  
  data_subset$Outlier_Score<-outlier_score_neglog
  data_subset$Outlier_Score_losses<-losses
  data_subset$Biased<-ifelse(losses>=outlier_cutoff_val,"biased","not")
  return(data_subset)
}
}
