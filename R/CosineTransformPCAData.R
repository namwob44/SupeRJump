#' Cosine Transform for PCA Data
#' cosine transforming data gives better representation than euclidean distance. We sorta cheat creating cosine dist with the norm below.

#' @param sce Put in the Seurat object. Be sure PCA in suerat has been ran first
#'
#' @return It returns a matrix that transforms PCA into cosine space
#' @export
#'
#' @examples
CosineTransformPCAData<-function(sce){
  raw_data_to_test<-as.data.frame(sce@reductions[["pca"]]@cell.embeddings)
  Z<-raw_data_to_test/wordspace::rowNorms(as.matrix(raw_data_to_test),method = "euclidean",p = 2)
  return(Z)
}
