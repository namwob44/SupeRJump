#' Get Basic multivariate linear model of fates and and TFs
#'
#' @param poised_cells 
#' @param model_list_check 
#' @param a 
#'
#' @return
#' @export
#'
#' @examples
GetFateTFModel<-function(poised_cells,model_list_check,a){
  
  # check if more data points than TFs or Motifs.
  # keep TFs
  swapped=FALSE
  if(dim(poised_cells)[1]<dim(poised_cells[,a])[2]){
    print("swapping")
    temp<-model_list_check
    model_list_check<-a
    a<-temp
    swapped=TRUE
  }
  
  
  res_all<-lm(poised_cells[,model_list_check]%>%as.matrix~poised_cells[,colnames(poised_cells)[which(colnames(poised_cells)%in%a)]]%>%as.matrix)%>%summary()
  names(res_all) <- colnames(poised_cells[,model_list_check]%>%as.matrix)
  
  # summary is a list for each condition. Get the info we need: 
  res_new <- res_all %>% lapply(X = ., function(fit){
    
    scores <- as.vector(fit$coefficients[,3][-1])
    pvals <- as.vector(fit$coefficients[,4][-1])
    sources <- colnames(poised_cells[,colnames(poised_cells)[which(colnames(poised_cells)%in%a)]]%>%as.matrix)
    diff_n <- length(sources) - length(scores)
    if (diff_n > 0) {
      print("PANIC")
      #stop(stringr::str_glue('After intersecting mat and network, at least {diff_n} sources in the network are colinear with other sources.
      #Cannot fit a linear model with colinear covariables, please remove them.
      #Please run decoupleR::check_corr to see what regulators are correlated.'))
    }
    tibble(score=scores, p_value=pvals, source=sources)
  }) %>% bind_rows(.id = "condition") %>%
    dplyr::mutate(statistic = "mlm", .before= 1) %>%
    dplyr::select(statistic, source, condition,
                  score, p_value)%>%
    dplyr::mutate(p_value=ifelse(p_value==0,min(p_value[p_value!=0]),p_value))%>%
    dplyr::mutate(sign_score = (-1)*log10(p_value)*sign(score))
  
  if(swapped==TRUE){
    res_new<-res_new%>%dplyr::select(statistic,condition,source,score,p_value,sign_score)
    colnames(res_new)<-c("statistic","source","condition","score","p_value","sign_score")
  }
  res_new$p_value<-ifelse(res_new$p_value==0,min(res_new$p_value[res_new$p_value!=0]),res_new$p_value)
  
  return(res_new)
}
