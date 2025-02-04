#' Normalize the probability to cut thresholds and preserve some states
#'
#' @param raw_TPM 
#' @param Z 
#'
#' @return
#' @export
#'
#' @examples
NormalizeTPM<-function(raw_TPM,Z){
  Prob_Xt_X0_mat_reduced <- raw_TPM
  
  for(sink_iter in 1:length(unique(Z$group))){
    current_cell_names<-Z%>%
      tibble::rownames_to_column(var="row_names")%>%
      dplyr::filter(group==as.character(unique(Z$group)[sink_iter]))%>%
      pull(row_names)
    for(other_sink_iter in 1:length(unique(Z$group))){
      other_cell_names<-Z%>%
        tibble::rownames_to_column(var="row_names")%>%
        dplyr::filter(group==as.character(unique(Z$group)[other_sink_iter]))%>%
        pull(row_names)
      if(sink_iter!=other_sink_iter){
        temp_block<-Prob_Xt_X0_mat_reduced[rownames(Prob_Xt_X0_mat_reduced)%in%current_cell_names,colnames(Prob_Xt_X0_mat_reduced)%in%other_cell_names]
        temp_block<-ifelse(temp_block<(1/(dim(Prob_Xt_X0_mat_reduced)[1])),NA,temp_block)
        Prob_Xt_X0_mat_reduced[rownames(Prob_Xt_X0_mat_reduced)%in%rownames(temp_block),colnames(Prob_Xt_X0_mat_reduced)%in%colnames(temp_block)]<-temp_block
        
      }
    }
  }
  #Prob_Xt_X0_mat_reduced <- ifelse(Prob_Xt_X0_mat_reduced<(1/dim(Prob_Xt_X0_mat_reduced)[2]),NA,Prob_Xt_X0_mat_reduced)
  #Prob_Xt_X0_mat_reduced<- Prob_Xt_X0_mat_reduced/rowSums(Prob_Xt_X0_mat_reduced,na.rm=T)
  Prob_Xt_X0_mat_reduced<-ifelse(is.na(Prob_Xt_X0_mat_reduced),0,Prob_Xt_X0_mat_reduced)
  return(Prob_Xt_X0_mat_reduced)  
}
