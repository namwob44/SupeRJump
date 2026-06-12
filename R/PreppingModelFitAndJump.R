
#' Cosine Transform for PCA Data
#' @description
#' cosine transforming data gives better representation than euclidean distance. We sorta cheat creating cosine dist with the norm below.
#'
#' @param seurat_obj Put in the Seurat object. Be sure PCA in seurat has been ran first since we pull the cell embeddings from the pca layer.
#'
#' @return It returns a matrix that transforms PCA into cosine space
#' @importFrom wordspace rowNorms
#' @import Seurat
#' @examples \dontrun{seurat_obj<-CosineTransformPCAData(seurat_obj)
#' state_grouping_column_name="custom_cell_classes_fine"
#' pseudotime_column_name = "Combined_Ordering"
#' batch_correction_column_name="sample"
#' Z <-seurat_obj[["transformed_pca"]]@cell.embeddings%>%
#' as.data.frame%>%
#' tibble::rownames_to_column(var="row_names")%>%
#' dplyr::inner_join(Seurat::FetchData(seurat_obj,vars=c(state_grouping_column_name,pseudotime_column_name,batch_correction_column_name))%>%
#' tibble::rownames_to_column(var="row_names"),by="row_names")%>%
#' tibble::column_to_rownames(var="row_names")%>%
#' dplyr::ungroup()%>%
#' dplyr::arrange(.data[[pseudotime_column_name]])}
#' @export
CosineTransformPCAData<-function(seurat_obj){
  raw_data_to_test<-as.data.frame(seurat_obj@reductions[["pca"]]@cell.embeddings)
  Z<-raw_data_to_test/wordspace::rowNorms(as.matrix(raw_data_to_test),method = "euclidean",p = 2)

  seurat_obj[["transformed_pca"]]<-Seurat::CreateDimReducObject(embeddings = as.matrix(Z),
                                                                key = "PCT_",
                                                                stdev = seurat_obj[["pca"]]@stdev)

  return(seurat_obj)
}


#' Get Pseudo Ordering for each geneset list.
#'
#' Note that this function will handle matching the names for the genesets to the same format as the seurat rows. This will work regardless human or mouse.
#'
#' @param seurat_obj is the seurat object you wish to get pseudo time ordering of
#' @param supervised_list_of_list This is either a list of gene names, or a list of lists of gene names you want scores for
#' @param name_of_score Individual scores can have names
#'
#' @return it returns the seurat object
#' @export
#' @importFrom Seurat AddModuleScore
#'
#' @examples \dontrun{
#' Idents(seurat_obj)<-"custom_cell_classes_fine"
#' DefaultAssay(seurat_obj)<-"RNA"
#' markers.differences<-FindAllMarkers(seurat_obj)
#' terminal_state_list<-c("early_Erp","MkP","Lymph","Mast", "Neu", "cDC", "Baso")
#' reduced_markers<-markers.differences%>%filter(cluster%in%terminal_state_list)
#' gene_list_from_findallmarkers<-list()
#' for(iter in 1:length(terminal_state_list)){
#'   gene_list_from_findallmarkers[[iter]]<-reduced_markers%>%
#'     filter(p_val_adj<0.01)%>%
#'        filter(avg_log2FC>1)%>%
#'            filter(cluster==terminal_state_list[iter])%>%pull(gene)%>%unique
#'     }
#' seurat_obj<-GetPseudoOrdering(seurat_obj,gene_list_from_findallmarkers,terminal_state_list)
#' seurat_obj<-seurat_obj<-CombinePseudoOrdering(seurat_obj,terminal_state_list)
#' }
#'
GetPseudoOrdering<-function(seurat_obj,supervised_list_of_list,name_of_score){

  # These inner helper functions make sure that dashes, dots, or anyway someone separates some genes to a standard format.
  normalize_names <- function(x) {
    # insert space before camelCase transitions
    x <- gsub("([a-z])([A-Z])", "\\1 \\2", x)
    # replace underscores/dashes/dots with space
    x <- gsub("[_.-]+", " ", x)
    # trim + lowercase + remove spaces
    x <- tolower(trimws(x))
    x <- gsub("\\s+", "", x)
    return(x)
  }
  match_to_reference <- function(reference_names, query_names) {
    ref_norm <- normalize_names(reference_names)
    qry_norm <- normalize_names(query_names)
    # match normalized query to normalized reference
    match_idx <- match(qry_norm, ref_norm)
    # return matched names in original reference format
    matched <- reference_names[match_idx]
    return(matched)
  }

  if(class(supervised_list_of_list)=="list"){

    for(iter in 1:(length(supervised_list_of_list))){
      seurat_obj<- Seurat::AddModuleScore(seurat_obj,
                                          features =list("score"=match_to_reference(reference_names = rownames(seurat_obj),query_names =supervised_list_of_list[[iter]])),
                                          name=paste0("score",iter),
                                          search =TRUE)
      seurat_obj@meta.data[,dim(seurat_obj@meta.data)[2]]<-seurat_obj@meta.data[,dim(seurat_obj@meta.data)[2]]-min(seurat_obj@meta.data[,dim(seurat_obj@meta.data)[2]]) # this will make lowest value 0.
    }

  }else{
    seurat_obj<- Seurat::AddModuleScore(seurat_obj,
                                        features =list("score"=match_to_reference(reference_names = rownames(seurat_obj),query_names =supervised_list_of_list)),
                                        name=as.character(name_of_score),
                                        search =TRUE)

    seurat_obj@meta.data[,dim(seurat_obj@meta.data)[2]]<-seurat_obj@meta.data[,dim(seurat_obj@meta.data)[2]]-min(seurat_obj@meta.data[,dim(seurat_obj@meta.data)[2]]) # this will make lowest value 0.
  }

  score_cols <- grep("^score", colnames(seurat_obj@meta.data), value = TRUE)
  scores_mat <- as.matrix(seurat_obj@meta.data[, score_cols])
  colnames(scores_mat)<-name_of_score
  module_assay <- Seurat::CreateAssayObject(data = t(scores_mat))
  seurat_obj[["ModuleScores"]] <- module_assay

  return(seurat_obj)
}



#' Combining Pseudo Ordering with Entropy Efficiency
#' This strategy deploys a
#' @param seurat_obj This is the seurat object with the scores for the lineages of interest
#' @param terminal_state_list this is the list of terminal lineages set earlier. Names are irrelevant. We just take the size of the terminal lineages to grab the last N score columns in the metadata
#'
#' @return This returns another column in the metadata of your seurat object called Combined_Ordering
#' @export
#' @importFrom Matrix rowSums
#' @import Seurat
#'
#' @examples \dontrun{
#' Idents(seurat_obj)<-"custom_cell_classes_fine"
#' DefaultAssay(seurat_obj)<-"RNA"
#' markers.differences<-FindAllMarkers(seurat_obj)
#' terminal_state_list<-c("early_Erp","MkP","Lymph","Mast", "Neu", "cDC", "Baso")
#' reduced_markers<-markers.differences%>%filter(cluster%in%terminal_state_list)
#' gene_list_from_findallmarkers<-list()
#' for(iter in 1:length(terminal_state_list)){
#'   gene_list_from_findallmarkers[[iter]]<-reduced_markers%>%
#'     filter(p_val_adj<0.01)%>%
#'        filter(avg_log2FC>1)%>%
#'            filter(cluster==terminal_state_list[iter])%>%pull(gene)%>%unique
#'     }
#' seurat_obj<-GetPseudoOrdering(seurat_obj,gene_list_from_findallmarkers,terminal_state_list)
#' seurat_obj<-seurat_obj<-CombinePseudoOrdering(seurat_obj,terminal_state_list)
#' }
#'
CombinePseudoOrdering<-function(seurat_obj,terminal_state_list){
  # obtains all scores generated from GetPseudoOrdering
  Seurat::DefaultAssay(seurat_obj)<-"ModuleScores"
  score_mat<-t(Seurat::GetAssayData(seurat_obj)%>%as.matrix)
  #score_mat<-seurat_obj@meta.data[,grep("score",colnames(seurat_obj@meta.data))]

  # maximum uniform distribution
  uniform_dist<-rep(x=1,times=length(terminal_state_list))
  test_prob <-uniform_dist/sum(uniform_dist)
  # the entropy score for each cell
  score_mat_dist <- score_mat/Matrix::rowSums(score_mat)
  Entropy_score <- (-1) * Matrix::rowSums(score_mat_dist *
                                            log2(score_mat_dist), na.rm = TRUE)
  Entropy_max = -sum(test_prob * log2(test_prob), na.rm = TRUE)
  # the efficiency lets us rank cells that most represent a uniform distribution which represents non-committed lineage specific cells.
  # The subtraction of 1
  Efficiency = Entropy_score/Entropy_max
  if(any(grepl("Combined_Ordering",colnames(seurat_obj@meta.data)))){
    seurat_obj@meta.data[, "Combined_Ordering"] <- 1 - Efficiency
  }else{
    seurat_obj@meta.data[, ncol(seurat_obj@meta.data) + 1] <- 1 - Efficiency
    colnames(seurat_obj@meta.data)[ncol(seurat_obj@meta.data)] <- "Combined_Ordering"
  }

  seurat_obj[["ModuleScores"]] <- Seurat::CreateAssayObject(data = rbind(t(score_mat), Combined_Ordering =  seurat_obj@meta.data$Combined_Ordering))
  return(seurat_obj)

}


#' Get Sink data points
#'
#' this is a summarized group for points and is used for model fits. It will pull all groups that listed as states.
#' @param seurat_obj The seurat object with a transformed PC space.
#' @param state_grouping_column_name the column name to be the states you determine your model by, often it is cell type metadata in seurat.
#'
#' @return A set of data points for eigen cells used for the jump-drift-diffusion model
#' @export
#' @import dplyr
#' @importFrom pracma randn
#' @import Seurat
#' @import tibble
#'
#' @examples \dontrun{
#' Y<-GetSinkData(seurat_obj,state_grouping_column_name="custom_cell_classes_fine",pseudotime_column_name="Combined_Ordering")
#' }
GetSinkData<-function(seurat_obj,state_grouping_column_name="custom_cell_classes_fine",pseudotime_column_name="Combined_Ordering"){


  Y<-seurat_obj[["transformed_pca"]]@cell.embeddings%>%as.data.frame%>%
    tibble::rownames_to_column(var="row_names")%>%
    dplyr::inner_join(Seurat::FetchData(seurat_obj,vars=c(pseudotime_column_name,state_grouping_column_name))%>%
                        tibble::rownames_to_column(var="row_names"),
                      by="row_names")%>%
    tibble::column_to_rownames(var="row_names")%>%
    dplyr::group_by(.data[[state_grouping_column_name]])%>%
    dplyr::summarise(dplyr::across(dplyr::everything(), mean))
  # Add random noise to everything so things don't hit singularity
  Y[,2:ncol(Y)]<-Y[,2:ncol(Y)]+0.001*pracma::randn(n=1)

  return(Y)
}



#' Cluster Assignment for Cells (deprecated)
#'
#'# get group assignments from the seurat object post seurat work
#' @param cosine_data The transformed PCA matrix we want to add another column to
#' @param seurat_obj the seurat object we want to pull information from (i.e. cell type or group)
#' @param column_name The column name we want to get out of the seurat object.
#'
#' @return The input matrix with an addition column added for the cell type group.
#' @export
#'
#' @examples \dontrun{Z<-ClusterAssignColumn(Z,seurat_obj,column_name="custom_cell_classes_fine")}
ClusterAssignColumn <-function(cosine_data,seurat_obj,column_name){
  cosine_data$group =seurat_obj@meta.data[,column_name]
  cosine_data$group = factor(cosine_data$group,levels=unique(cosine_data$group))
  return(cosine_data)
}

#' Assign Pseudotime Scores  (deprecated)
#'
#' This takes the pseudotime column and assigns it to the pseudotime (PT_score) column.
#' @param cosine_data the transformed PCA matrix/dataframe we want to add pseudotime column too.
#' @param seurat_obj The seurat object to pull the metadata column about pseudotime from
#' @param column_name column name you want to use for pseudotime
#'
#' @return The input matrix with an addition column added for the pseudotime (PT)
#' @export
#'
#' @examples \dontrun{Z<-AssignPTScores(Z,seurat_obj,column_name="Combined_Ordering")}
AssignPTScores<-function(cosine_data,seurat_obj,column_name){
  cosine_data$PT_score<-seurat_obj@meta.data[,column_name]
  return(cosine_data)
}

