
#' Get Membership from PseudoTime only
#'
#' @param raw_data_to_test 
#'
#' @return
#' @export
#'
#' @examples
GetMembershipWithPTOnly<-function(raw_data_to_test){
  # raw_data_to_test<-as.data.frame(sce@reductions[["pca"]]@cell.embeddings)
  # raw_data_to_test$group =sce@meta.data$group
  # raw_data_to_test$group<-factor(raw_data_to_test$group,levels=(unique(raw_data_to_test$group)))
  
  Mean_mus_clust<-raw_data_to_test%>%
    group_by(group)%>%
    summarise(across(PT_score, mean))
  
  var_mus_clust<-raw_data_to_test%>%
    group_by(group)%>%
    summarise(across(PT_score, var))
  membership_raw = matrix(data=NA_real_,nrow=dim(raw_data_to_test)[1],ncol =length(unique(raw_data_to_test$group)))
  for(cluster_it in 1:length(unique(raw_data_to_test$group))){
    membership_raw[,cluster_it]<-stats::dnorm(raw_data_to_test%>%
                                                dplyr::select(PT_score)%>%
                                                as.matrix,
                                              mean = as.numeric(Mean_mus_clust[cluster_it,2]),
                                              sd = as.numeric(var_mus_clust[cluster_it,2]))
  }
  colnames(membership_raw)<-Mean_mus_clust$group
  rownames(membership_raw)<-rownames(raw_data_to_test)
  membership_norm <-membership_raw/rowSums(membership_raw)
  
  return(membership_norm)
}
