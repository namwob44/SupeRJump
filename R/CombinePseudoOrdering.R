
#' Combining Pseudo Ordering with Entropy Efficiency
#'
#' @param sce This is the seurat object with the scores for the lineages of interest
#' @param terminal_state_list this is the list of terminal lineages set earlier. Names are irrelevant. We just take the size of the terminal lineages to grab the last N score columns in the metadata
#'
#' @return This returns another column in the metadata of your suerat object called Combined_Ordering
#' @export
#'
#' @examples
CombinePseudoOrdering<-function(sce,terminal_state_list){
  score_mat<-sce@meta.data[,(ncol(sce@meta.data)-length(terminal_state_list)+1):ncol(sce@meta.data)]
  uniform_dist<-rep(x=1,times=length(terminal_state_list))
  test_prob <-uniform_dist/sum(uniform_dist)

  score_mat_dist <- score_mat/rowSums(score_mat)
  Entropy_score <- (-1) * rowSums(score_mat_dist *
                                    log2(score_mat_dist), na.rm = TRUE)
  Entropy_max = -sum(test_prob * log2(test_prob), na.rm = TRUE)
  Efficiency = Entropy_score/Entropy_max

  sce@meta.data[,ncol(sce@meta.data)+1]<-1-Efficiency
  colnames(sce@meta.data)[ncol(sce@meta.data)]<- "Combined_Ordering"
  return(sce)

}
