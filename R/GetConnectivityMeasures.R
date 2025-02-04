#' Some connectivity measures from a TPM
#'
#' @param raw_TPM 
#' @param Z 
#'
#' @return
#' @export
#'
#' @examples
GetConnectivityMeasures<-function(raw_TPM,Z){
  Prob_Xt_X0_mat_reduced <- raw_TPM
  connectivitity_alpha<-matrix(NA_real_,ncol=length(unique(Z$group)),nrow=length(unique(Z$group)))
  connectivitity_gamma<-matrix(NA_real_,ncol=length(unique(Z$group)),nrow=length(unique(Z$group)))
  rownames(connectivitity_alpha)<-colnames(connectivitity_alpha)<-as.character(unique(Z$group))
  rownames(connectivitity_gamma)<-colnames(connectivitity_gamma)<-as.character(unique(Z$group))
  
  for(sink_iter in 1:length(unique(Z$group))){
    current_cell_names<-Z%>%
      filter(group==as.character(unique(Z$group)[sink_iter]))%>%
      tibble::rownames_to_column(var="row_names")%>%
      pull(row_names)
    for(other_sink_iter in 1:length(unique(Z$group))){
      other_cell_names<-Z%>%
        filter(group==as.character(unique(Z$group)[other_sink_iter]))%>%
        tibble::rownames_to_column(var="row_names")%>%
        pull(row_names)
      if(sink_iter!=other_sink_iter){
        temp_block<-Prob_Xt_X0_mat_reduced[rownames(Prob_Xt_X0_mat_reduced)%in%current_cell_names,colnames(Prob_Xt_X0_mat_reduced)%in%other_cell_names]
        temp_block<-ifelse(temp_block==0,NA,temp_block)
        check_edges<-ifelse(is.na(temp_block),0,1)
        check_edges_sum <-sum(sum(check_edges))
        check_vertices <-sum(dim(temp_block)[1]+dim(temp_block)[2])
        
        connectivitity_alpha[sink_iter,other_sink_iter] = (check_edges_sum-check_vertices)/( 0.5*(check_vertices*(check_vertices-1))-(check_vertices-1))
        connectivitity_gamma[sink_iter,other_sink_iter] = check_edges_sum/(0.5*check_vertices*(check_vertices-1))
        print(check_edges_sum)
        print(check_vertices)
        print(connectivitity_alpha[sink_iter,other_sink_iter])
        print(connectivitity_gamma[sink_iter,other_sink_iter])
        print(as.character(unique(Z$group)[sink_iter]))
        print(as.character(unique(Z$group)[other_sink_iter]))
      }
    }
  }
  #Prob_Xt_X0_mat_reduced <- ifelse(Prob_Xt_X0_mat_reduced<(1/dim(Prob_Xt_X0_mat_reduced)[2]),NA,Prob_Xt_X0_mat_reduced)
  return(list("alpha_connect"=connectivitity_alpha,
              "gamma_connect"=connectivitity_gamma))  
}