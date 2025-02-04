#' Get Sink data points
#' this is a summarized group for points and is used for model fits. It will pul all groups
#' @param cosine_data The input matrix for analysis.
#'
#' @return A set of data points for sink points in the jump-drift-diffusion model
#' @export
#'
#' @examples
GetSinkData<-function(cosine_data){
  Y<-cosine_data%>%
    group_by(group)%>%
    summarise(across(everything(), mean))
  # Add random noise to everything so things don't hit singularity
  idx = which(colnames(Y)=="group")
  Y[,-which(colnames(Y)=="group")]<-Y[,-which(colnames(Y)=="group")]+0.001*pracma::randn(n=1)

  return(Y)
}
