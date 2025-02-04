
#' Downsamples the input matrix
#' in case there is too many cells. This uses sigma points to estimate the cells in each group.
#' @param cosine_data The input matrix to downsample
#' @param start_cell if you have a known starting or initiating cell to keep.
#'
#' @return
#' @export
#'
#' @examples
DownSample<-function(cosine_data,start_cell=NULL){

  GetSigmaPoints<-function(clust_data){
    sigma_n <- (dim(clust_data)[2])-1 # number of dimensions but we have group so subtract 1
    sigma_lambda <- 3-(1*sigma_n) # 0.5 seems to get better spread for this data in my opinion. usually 0.5 should be 1.
    sigma_mu <- colMeans(clust_data%>%dplyr::select(-group))
    sigma_cov_mat <- cov(clust_data%>%dplyr::select(-group))
    #if(!is.positive.semi.definite(sigma_cov_mat)){#&!is.positive.definite(sigma_cov_mat)){
    #  sigma_cov_mat <- GetPDFromSPD(sigma_cov_mat)
    #}
    sigma_points <- cbind((sigma_mu),
                          (sigma_mu)+Re(expm::sqrtm((sigma_n+sigma_lambda)*sigma_cov_mat)),
                          (sigma_mu)-Re(expm::sqrtm((sigma_n+sigma_lambda)*sigma_cov_mat)))
    sigma_points[nrow(sigma_points),]<- ifelse(sigma_points[nrow(sigma_points),] < 0, 0, sigma_points[nrow(sigma_points),]) # this ensures time is positive
    sigma_points[nrow(sigma_points),]<- ifelse(sigma_points[nrow(sigma_points),] > 1, 1, sigma_points[nrow(sigma_points),]) # this ensures time doesn't go past our biggest value of 1.
    sigma_points<-as.data.frame(sigma_points)
    rownames(sigma_points)<-colnames(clust_data%>%dplyr::select(-group))
    colnames(sigma_points)<-gsub(pattern = "V",replacement = paste0(as.character((clust_data%>%dplyr::select(group))[1,1]),"_test_point"),x = colnames(sigma_points),ignore.case = T)
    sigma_points[nrow(sigma_points)+1,]<-rep(x = as.character((clust_data%>%dplyr::select(group))[1,1]),times=dim(sigma_points)[2])
    rownames(sigma_points)[nrow(sigma_points)]<-"group"
    sigma_points[nrow(sigma_points)+1,]<-rep(x = "Synthetic",times=dim(sigma_points)[2])
    rownames(sigma_points)[nrow(sigma_points)]<-"type"
    return(sigma_points)
  }
  node_names<-as.character(unique(cosine_data$group))
  sigma_points_full=NULL
  for(clust_iter in 1:length(node_names)){
    temp_data<-cosine_data%>%filter(group==node_names[clust_iter])
    temp_data$PT_score<-cosine_data%>%filter(group==node_names[clust_iter])%>%pull(PT_score)
    sigma_points_temp = GetSigmaPoints(temp_data)
    sigma_points_full<-dplyr::bind_cols(sigma_points_full,sigma_points_temp)
  }
  for_testing_sigma_points<-t(sigma_points_full)%>%as.data.frame
  for_testing_sigma_points[,1:(ncol(for_testing_sigma_points)-2)]<-apply(for_testing_sigma_points[,1:(ncol(for_testing_sigma_points)-2)],2,function(x){as.numeric(x)})

  if(is.null(start_cell)){
    for_testing_sigma_points<-dplyr::bind_rows((cosine_data%>%arrange(PT_score))[1,],for_testing_sigma_points)
  }else{
    for_testing_sigma_points<-dplyr::bind_rows(start_cell,for_testing_sigma_points)

  }


  for_testing_sigma_points[1,(ncol(for_testing_sigma_points)-1)]=0
  for_testing_sigma_points[1,ncol(for_testing_sigma_points)]="Real"
  for_testing_sigma_points<-for_testing_sigma_points%>%arrange(PT_score)
  for_testing_sigma_points$group<-factor(for_testing_sigma_points$group,levels=unique(for_testing_sigma_points$group))
  for_testing_sigma_points$type<-factor(for_testing_sigma_points$type,levels=unique(for_testing_sigma_points$type))

  return(for_testing_sigma_points)
}
