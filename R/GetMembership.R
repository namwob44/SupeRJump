#' Get Membership of cells through gaussian estimation
#'
#' @param raw_data_to_test this is the input matrix to incorporate it needs the column group to function
#'
#' @return It returns a matrix of membership assigment normalized per group
#' @export
#'
#' @examples
GetMembership<-function(raw_data_to_test){
  # raw_data_to_test<-as.data.frame(sce@reductions[["pca"]]@cell.embeddings)
  # raw_data_to_test$group =sce@meta.data$group
  # raw_data_to_test$group<-factor(raw_data_to_test$group,levels=(unique(raw_data_to_test$group)))

  Mean_mus_clust<-raw_data_to_test%>%
    group_by(group)%>%
    summarise(across(everything(), mean))

  var_mus_clust<-raw_data_to_test%>%
    group_by(group)%>%
    summarise(across(everything(), var))
  membership_raw = matrix(data=NA_real_,nrow=dim(raw_data_to_test)[1],ncol =length(unique(raw_data_to_test$group)))
  for(cluster_it in 1:length(unique(raw_data_to_test$group))){
    membership_raw[,cluster_it]<-mvtnorm::dmvnorm(raw_data_to_test%>%
                                                    dplyr::select(-group)%>%
                                                    as.matrix,
                                                  mean = as.matrix(Mean_mus_clust[cluster_it,2:length(Mean_mus_clust)]),
                                                  sigma = diag(as.vector(var_mus_clust[cluster_it,2:length(Mean_mus_clust)])))
  }
  colnames(membership_raw)<-Mean_mus_clust$group
  rownames(membership_raw)<-rownames(raw_data_to_test)
  membership_norm <-membership_raw/rowSums(membership_raw)

  return(membership_norm)
}
