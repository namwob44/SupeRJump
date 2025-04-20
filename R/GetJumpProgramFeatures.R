
#' Get a table of which genes are associated with Jumps
#'
#' @param sce_obj 
#' @param down_Z 
#' @param Y 
#' @param model_fits_list 
#' @param foldername 
#' @param top_n_PCs 
#'
#' @return
#' @export
#'
#' @examples
GetJumpProgramFeatures<-function(sce_obj,down_Z,Y,model_fits_list,foldername,top_n_PCs=10){
  
  list_of_lambdas<-lapply(model_fits_list[["model_fits"]],function(x){
    return(x[["solution"]][1])})%>%as.matrix
  idx_for_jump_models<-(which(list_of_lambdas<max(down_Z$PT_score)&(list_of_lambdas>0)))
  
  long_loading_pc_df_norm<-down_Z%>%
    tibble::rownames_to_column(var="Cell")%>%
    as.data.frame%>%
    tidyr::pivot_longer(cols=!c(Cell,type,group,PT_score),
                        names_to = "PC_name",
                        values_to="score")
  
  sink_levels<-as.character(Y$group)
  
  jump_table<-matrix(FALSE,nrow=length(unique(long_loading_pc_df_norm$PC_name)),ncol=length(unique(long_loading_pc_df_norm$group)))
  rownames(jump_table)<-unique(long_loading_pc_df_norm$PC_name)
  colnames(jump_table)<-levels(long_loading_pc_df_norm$group)
  for(iter in 1:length(idx_for_jump_models)){
    jump_table[idx_for_jump_models[iter]%%nrow(jump_table),sink_levels[ceiling(idx_for_jump_models[iter]/nrow(jump_table))]]<-TRUE
  }
  
  jump_df<-jump_table%>%
    as.data.frame%>%
    tibble::rownames_to_column(var="PC_names")%>%
    tidyr::pivot_longer(cols=!PC_names,names_to = "Sink_Type",values_to="Jump")%>%
    group_by(PC_names)%>%filter(Jump==TRUE)
  write.csv(jump_df,paste0(foldername,"jump_table.csv"))
  PCs_to_look_into<-(jump_df%>%pull(PC_names)%>%unique)
  print(PCs_to_look_into)
  feature_loadings_mat<-sce_obj@reductions[["pca"]]@feature.loadings # *sce_obj@reductions[["pca"]]@stdev
  
  scale_loadings_mat<-abs(apply(feature_loadings_mat,2,function(x){scale(x)}))*(sce_obj@reductions[["pca"]]@stdev/100)
  rownames(scale_loadings_mat)<-rownames(feature_loadings_mat)
  
  candidate_Genes_To_explore<-apply(scale_loadings_mat[,PCs_to_look_into],2,function(x){
    names(x)[which(x>quantile(x)[4]+1.5*IQR(x))]
    sort(x[which(x>quantile(x)[4]+1.5*IQR(x))])
  })
  
  
  #pheatmap::pheatmap(scale_loadings_mat,cluster_rows = T,cluster_cols = F,show_colnames = T,show_rownames = F,)
  
  
  return(candidate_Genes_To_explore)
  
  
  #for(iter in 1:length(PCs_to_look_into)){
  #  PC_importance<-abs(feature_loadings_mat[,PCs_to_look_into[iter]])
  #  
  #}
  
  
}
