#' Cluster Assignment for cells
#'# get group assignments from the sce object post seurat work
#' @param cosine_data The transformed PCA matrix we want to add
#' @param sce the seurat object we want to pull information from (i.e. cell type or group)
#' @param column_name The column name we want to get out of the seurat object.
#'
#' @return The input matrix with an addition column added for the cell type group.
#' @export
#'
#' @examples
ClusterAssignColumn <-function(cosine_data,sce,column_name){
  cosine_data$group =sce@meta.data[,column_name]
  cosine_data$group = factor(cosine_data$group,levels=unique(cosine_data$group))
  return(cosine_data)
}

#' Assign Pseudotime Scores
#'
#' @param cosine_data Input matrix to add a new column for.
#' @param sce The seurat object to pull the metadata column about pseudotime from
#' @param column_name column name you want to use for psuedotime
#'
#' @return The input matrix with an addition column added for the pseudotime (PT)
#' @export
#'
#' @examples
AssignPTScores<-function(cosine_data,sce,column_name){
  cosine_data$PT_score<-sce@meta.data[,column_name]
  return(cosine_data)
}
