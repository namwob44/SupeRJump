#' Get Pseudo Ordering
#'
#' @param sce is the seurat object you wish to get pseudo time ordering of
#' @param supervised_list_of_list This is either a list of gene names, or a list of lost of gene names you want scores for
#' @param name_of_score Individual scores can have names
#'
#' @return it returns the suerat object
#' @export
#'
#' @examples
GetPseudoOrdering<-function(sce,supervised_list_of_list,name_of_score){
  # the crazy gsub just makes it proper capitalization instead of all caps.
  # need to add check if already lower case or not.
  if(class(supervised_list_of_list)=="list"){
    for(iter in 1:(length(supervised_list_of_list))){
      sce<-Seurat::AddModuleScore(sce,features=
                                    #list("score"=c(gsub("(?<=\\b)([a-z])", "\\U\\1", (supervised_list_of_list[[iter]]), perl=TRUE))),
                                  list("score"=c(gsub("((?<=\\b.)[[:upper:]])", "\\L\\1", (supervised_list_of_list[[iter]]), perl=TRUE))),
                                  name=paste0("score",iter),
                                  search=TRUE,)
      sce@meta.data[,dim(sce@meta.data)[2]]<-sce@meta.data[,dim(sce@meta.data)[2]]-min(sce@meta.data[,dim(sce@meta.data)[2]]) # this will make lowest value 0.

    }
  }else{
    sce<-Seurat::AddModuleScore(sce,features=
                                  #list("score"=c(gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(supervised_list_of_list), perl=TRUE))),
                                   list("score"=c(gsub("((?<=\\b.)[[:upper:]])", "\\L\\1", (supervised_list_of_list[[iter]]), perl=TRUE))),
                                name=as.character(name_of_score),search = TRUE)
    sce@meta.data[,dim(sce@meta.data)[2]]<-sce@meta.data[,dim(sce@meta.data)[2]]-min(sce@meta.data[,dim(sce@meta.data)[2]]) # this will make lowest value 0.
  }
  return(sce)
}
