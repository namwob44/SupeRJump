#' Get a model based approach to determine Absorbing States
#' @description
#' This is a function that identifies the a relatively balanced absorbing state populations for each sample.
#'
#' @param seurat_obj the seurat object with a TPM assay and potentially the CMFPT
#' @param markers_df place the  is used if the mode="markers"
#' @param sink_cell_names a short cut to manually select cell names to be absorbing states default is NULL
#' @param mode if sink cell names is NULL, this selects which approach to use for absorbing states. currently choose 1 of 2 methods: "markers", "CMFPT"
#' @param state_grouping_column_name the column name to be the states you determine your model by, often it is cell type metadata in seurat.
#' @param pseudotime_column_name the column name to be the pseudotime, often it is Combined_Ordering.
#' @param batch_correction_column_name the column that contains sample/mouse/patient differences.
#' @param terminal_state_list this is a character vector specifying which cell states within state_grouping_column_name to keep (terminal states)
#'
#' @return seurat object with an updated metadata column called absorbing_states. NA means not selected, while 1 is selected.
#' @export
#' @import dplyr
#' @import tibble
#' @import Seurat
#'
#' @examples \dontrun{
#' randomly_selected_cell_names<-sample(colnames(subset(seurat_obj,subset= cluster_annotations%in%terminal_state_list)),200,replace = F)
#' seurat_obj<-GetAutomaticAbsorbingCellAssignment(seurat_obj,sink_cell_names = randomly_selected_cell_names)
#' }
#'
GetAutomaticAbsorbingCellAssignment<-function(seurat_obj,
                                              markers_df=NULL,
                                              sink_cell_names=NULL,
                                              mode=c("markers","CMFPT"),
                                              state_grouping_column_name,
                                              pseudotime_column_name,
                                              batch_correction_column_name,
                                              terminal_state_list){

  #choices for sink cell selection:
  # self selection
  # PT
  # topK
  # CMFPT

  if(!is.null(sink_cell_names)){
    temp_data_frame<-data.frame(row_names=colnames(seurat_obj))%>%
      dplyr::mutate(absorbing_states=ifelse(row_names%in%sink_cell_names,1,NA))%>%
      tibble::column_to_rownames(var="row_names")
    seurat_obj<-Seurat::AddMetaData(seurat_obj,metadata = temp_data_frame)
    return(seurat_obj)
  }


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



  # general rule, FOR EACH BATCH given lowest populations(of more than 4 cells), we want absorbing cells per lineage: n_absorbing_cells = min(sqrt(n_lineage),3*sqrt(n_lowest_nonzero))

  #ok so let's first find our numbers!
  lineage_group_numbers_to_pull<-Seurat::FetchData(seurat_obj,vars=c(batch_correction_column_name,state_grouping_column_name))%>%
    data.frame%>%
    dplyr::filter(.data[[state_grouping_column_name]]%in%terminal_state_list)%>%
    dplyr::group_by(.data[[batch_correction_column_name]],.data[[state_grouping_column_name]])%>%
    dplyr::summarise(count=floor(sqrt(n())))%>%
    dplyr::ungroup()%>%
    dplyr::filter(count>=2)%>% # this makes sure we have at least 4 cells in a cluster to consider it terminal...
    dplyr::group_by(.data[[batch_correction_column_name]])%>%
    dplyr::mutate(min_count=min(count))%>%
    dplyr::mutate(absorbing_count = ifelse(count>3*min_count,3*min_count,count))%>%
    dplyr::ungroup()%>%
    dplyr::select(.data[[batch_correction_column_name]],.data[[state_grouping_column_name]],absorbing_count)


  if(mode=="CMFPT"&(c("CMFPT") %in% names(seurat_obj@graphs))){
    message("running CMFPT mode for automatic absorbing cell selection")

      s<-summary(seurat_obj@graphs[["CMFPT"]])
      sink_cell_names<-data.frame(current_state = rownames(seurat_obj@graphs[["CMFPT"]])[s$i],next_state= colnames(seurat_obj@graphs[["CMFPT"]])[s$j],CMFPT = s$x)%>%
        dplyr::filter(!is.na(CMFPT))%>%
      dplyr::inner_join(Seurat::FetchData(seurat_obj,vars=c(state_grouping_column_name))%>%
                          tibble::rownames_to_column(var="next_state")%>%
                          dplyr::filter(.data[[state_grouping_column_name]]%in%terminal_state_list),
                          by="next_state")%>%
      dplyr::rename("cluster_to"=state_grouping_column_name)%>%
      dplyr::inner_join(Seurat::FetchData(seurat_obj,vars=c(state_grouping_column_name,batch_correction_column_name))%>%
                            tibble::rownames_to_column(var="current_state")%>%
                          dplyr::filter(!.data[[state_grouping_column_name]]%in%terminal_state_list)%>%
                          dplyr::select(current_state,.data[[batch_correction_column_name]]),
                        by="current_state")%>%
        dplyr::group_by(cluster_to,.data[[batch_correction_column_name]])%>%
        #dplyr::filter(CMFPT>=stats::quantile(CMFPT,0.8))%>%
        dplyr::mutate(CVaR = mean(CMFPT))%>% # here we would take this many cells instead of display? Still need to scale I think?
        dplyr::ungroup()%>%
        dplyr::group_by(next_state)%>%
        dplyr::mutate(num_of_init_cells_that_are_far=n())%>%
        dplyr::distinct(next_state,cluster_to,.data[[batch_correction_column_name]],CVaR,num_of_init_cells_that_are_far)%>%
        dplyr::ungroup()%>%
        dplyr::mutate(Prob_best_representation=num_of_init_cells_that_are_far/max(num_of_init_cells_that_are_far))%>%
        dplyr::group_by(cluster_to)%>%
        dplyr::mutate(Prob_Far=CVaR/max(CVaR))%>%
        dplyr::ungroup()%>%
        dplyr::mutate(weighted_prob = Prob_best_representation*Prob_Far)%>%
        dplyr::distinct(next_state,cluster_to,.data[[batch_correction_column_name]],weighted_prob)%>%
        dplyr::group_by(next_state)%>%
        dplyr::filter(weighted_prob==max(weighted_prob))%>%
        dplyr::ungroup()%>%
        dplyr::inner_join(lineage_group_numbers_to_pull,by=c(batch_correction_column_name,cluster_to=state_grouping_column_name))%>%
        dplyr::select(next_state,.data[[batch_correction_column_name]],cluster_to,weighted_prob,absorbing_count)%>%
        dplyr::group_by(.data[[batch_correction_column_name]],cluster_to)%>%
        dplyr::group_modify(~ dplyr::slice_sample(.x,n = min(unique(.x$absorbing_count), nrow(.x)) ,weight_by = .x$weighted_prob,replace = FALSE)) %>%
        dplyr::ungroup()%>%
        dplyr::pull(next_state)

      temp_data_frame<-data.frame(row_names=colnames(seurat_obj))%>%
        dplyr::mutate(absorbing_states=ifelse(row_names%in%sink_cell_names,1,NA))%>%
        tibble::column_to_rownames(var="row_names")
      seurat_obj<-Seurat::AddMetaData(seurat_obj,metadata = temp_data_frame)
      return(seurat_obj)
  }
  else{
    if(mode=="CMFPT"&!(c("CMFPT") %in% names(seurat_obj@graphs))){
      warning("CMFPT graph is not found, please run GetBatchCorrectedCMFPT first")
    }
    message("defaulting to pseudotime and markers mode for automatic absorbing cell selection")

    set.seed(303)
    unique_state_types<-Seurat::FetchData(seurat_obj,vars=c(state_grouping_column_name))%>%dplyr::pull(.data[[state_grouping_column_name]])%>%unique
    temp_matrix<-matrix(data = NA_real_,nrow=ncol(seurat_obj),ncol = length(unique_state_types))
    colnames(temp_matrix)<-unique_state_types
    rownames(temp_matrix)<-colnames(seurat_obj)
    Seurat::DefaultAssay(seurat_obj)<-"RNA"
    for(iter in 1:length(unique_state_types)){
      supervised_list_of_list<-markers_df%>%
        dplyr::filter(p_val_adj<0.01)%>%
        dplyr::filter(avg_log2FC>1)%>%
        dplyr::filter(cluster==unique_state_types[iter])%>%
        dplyr::pull(gene)%>%
        unique
      print(unique_state_types[iter])
      temp_scores<-Seurat::AddModuleScore(seurat_obj,
                                          features =list("score"=match_to_reference(reference_names = rownames(seurat_obj),query_names =(supervised_list_of_list))),
                                          name="trash",
                                          search =TRUE)
      temp_matrix[,iter]<-temp_scores@meta.data$trash1-(min(temp_scores@meta.data$trash1))
    }

    # Maybe we just do probability of PT and marker score and call it since I allow custom cells anyway so I can just use my genes only one for paper scripts... Then exchange CMFPT for both PT and marker score?
    # so what columns are needed? rownames of cells, batches, terminal_state_grouping/state_grouping_column_name, num_of_cells per batch-group pair. weighted_probability for a given cell.

    # Ok next step that is downstream and hard to do. We need to get the probabilities setup for stochastically choosing which cells given a "score".
    # I view 2 necessary ways to weight these. We want terminal cells that are outliery far for the maximum number of cells in an initial cluster (for PT) AND
    # we want cells that are furthest away so we can also break ties.
    # ok we need to sample on each group

    sink_cell_names<-Seurat::FetchData(seurat_obj,vars=c(batch_correction_column_name,pseudotime_column_name,state_grouping_column_name))%>%
      as.data.frame%>%
      tibble::rownames_to_column(var="row_names")%>%
      dplyr::inner_join(temp_matrix%>%
                          as.data.frame%>%
                          tibble::rownames_to_column(var="row_names")%>%
                          tidyr::pivot_longer(cols=!c(row_names),names_to="Group_scores",values_to="Marker_scores"),by="row_names")%>%
      dplyr::filter(.data[[state_grouping_column_name]]==Group_scores)%>%
      dplyr::inner_join(lineage_group_numbers_to_pull,by=c(batch_correction_column_name,state_grouping_column_name))%>%
      dplyr::mutate(Marker_scores=as.numeric(Marker_scores))%>%
      #dplyr::group_by(.data[[state_grouping_column_name]])%>%
      #dplyr::filter((Marker_scores>=stats::quantile(Marker_scores,0.9)))%>%
      dplyr::filter(.data[[state_grouping_column_name]]%in%terminal_state_list)%>%
      dplyr::group_by(.data[[batch_correction_column_name]],.data[[state_grouping_column_name]])%>%
      dplyr::mutate(Prob_Self_ID=Marker_scores/max(Marker_scores))%>% # majority voting on who views this as an outlier?
      dplyr::mutate(Prob_Far=.data[[pseudotime_column_name]]/max(.data[[pseudotime_column_name]]))%>%
      dplyr::ungroup()%>%
      dplyr::mutate(weighted_prob = Prob_Self_ID*Prob_Far)%>%
      dplyr::select(row_names,.data[[batch_correction_column_name]],.data[[state_grouping_column_name]],weighted_prob,absorbing_count)%>%
      dplyr::group_by(.data[[batch_correction_column_name]],.data[[state_grouping_column_name]])%>%
      dplyr::group_modify(~ dplyr::slice_sample(.x,n = min(unique(.x$absorbing_count), nrow(.x)) ,weight_by = .x$weighted_prob,replace = FALSE)) %>%
      dplyr::ungroup()%>%
      dplyr::pull(row_names)

      temp_data_frame<-data.frame(row_names=colnames(seurat_obj))%>%
        dplyr::mutate(absorbing_states=ifelse(row_names%in%sink_cell_names,1,NA))%>%
        tibble::column_to_rownames(var="row_names")
      seurat_obj<-Seurat::AddMetaData(seurat_obj,metadata = temp_data_frame)
      return(seurat_obj)


  }

}

