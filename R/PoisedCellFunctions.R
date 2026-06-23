#' Get Preferentially Biased Cells
#' @description
#' This function determines within a cluster group which cells are selectively biased toward the lineages selected. Important to note, we develop an assay we provide both an outlier score in data slot, and a bias score(recommended) in the scale.data slot.
#'
#' @param seurat_obj this must contain dataframe must contain membership and
#' @param lineages_to_compare this is the row names of the lineages_fates assay you wish to include in the study
#'
#' @return it returns a new assay called bias_scores. These are continuous values that dictate a score to use.
#' @export
#' @import dplyr
#' @import tibble
#' @import Seurat
#'
#'
#' @examples \dontrun{seurat_obj<-GetAllClassifiedPoisedCells(seurat_obj,state_grouping_column_name="custom_cell_class_fine",lineages_to_compare=c("Mono","Neu","Ery"))}
GetAllClassifiedPoisedCells<-function(seurat_obj,state_grouping_column_name,lineages_to_compare){



  # we need to make a data_input which contains state_grouping_column, Membership assay and lineage fate assay too.

  data_input <- t(Seurat::GetAssayData(seurat_obj,assay = "Membership"))%>%
    as.data.frame%>%
    dplyr::rename_with(~paste0(.,"_membership"))%>%
    tibble::rownames_to_column(var="row_names")%>%
    dplyr::inner_join(
      t(Seurat::GetAssayData(seurat_obj,assay="lineage_fates",slot="scale.data"))%>%as.data.frame%>%
        tibble::rownames_to_column(var="row_names"),
      by="row_names",suffix=c("_membership",""))%>%
    dplyr::inner_join(Seurat::FetchData(seurat_obj,vars=c(state_grouping_column_name))%>%
                        tibble::rownames_to_column(var="row_names"),
    by="row_names")%>%
    tibble::column_to_rownames(var="row_names")

  colnames(data_input)<-gsub("-","_",colnames(data_input))

  # let's create placeholder columns of interest, along with the outlier scores
  namevector <- paste0("Biased_",lineages_to_compare)
  data_input[ , namevector] <- NA
  namevector <- paste0("Outlier_Score_",lineages_to_compare)
  data_input[ , namevector] <- NA
  namevector <- paste0("Outlier_Score_loss_",lineages_to_compare)
  data_input[ , namevector] <- NA



  for(cluster_iter in 1:length(as.character(unique(data_input[[state_grouping_column_name]])))){
    for(lineage_iter in 1:length(lineages_to_compare)){
      print("in loop")

      current_group=as.character(unique(data_input[[state_grouping_column_name]]))[cluster_iter]
      poised_vector<-ClassifyPoisedCellsSingleLineage(data_subset = data_input%>%
                                                        dplyr::filter(.data[[state_grouping_column_name]]==as.character(unique(data_input[[state_grouping_column_name]]))[cluster_iter]),
                                                      c(lineages_to_compare[lineage_iter],
                                                        paste0(current_group,"_membership"),
                                                        paste0(lineages_to_compare[lineage_iter],"_membership")))
      print("after inner function")

      data_input[rownames(poised_vector),paste0("Biased_",as.character(lineages_to_compare[lineage_iter]))]<-poised_vector[,"Biased"]
      data_input[rownames(poised_vector),paste0("Outlier_Score_loss_",as.character(lineages_to_compare[lineage_iter]))]<-poised_vector[,"Outlier_Score_losses"]
      data_input[rownames(poised_vector),paste0("Outlier_Score_",as.character(lineages_to_compare[lineage_iter]))]<-poised_vector[,"Outlier_Score"]

    }
  }
  #return(data_input)

  outlier_score_data <- data_input[,  paste0("Outlier_Score_", lineages_to_compare), drop = FALSE]
  outlier_loss_data <- data_input[,  paste0("Outlier_Score_loss_", lineages_to_compare), drop = FALSE]
  biased_data <- data_input[, paste0("Biased_", lineages_to_compare), drop = FALSE]
  # Transpose to lineages x cells
  outlier_score_mat <- t(as.matrix(outlier_score_data))

  outlier_loss_mat <- t(as.matrix(outlier_loss_data))

    biased_mat_raw <- t(as.matrix(biased_data))

    rownames(outlier_score_mat) <- rownames(outlier_loss_mat) <- rownames(biased_mat_raw) <- lineages_to_compare

  outlier_loss_mat <- matrix(as.numeric(outlier_loss_mat),
                             nrow = nrow(outlier_loss_mat),
                             ncol = ncol(outlier_loss_mat),
                             dimnames = dimnames(outlier_loss_mat))

  outlier_score_mat <- matrix(as.numeric(outlier_score_mat),
                              nrow = nrow(outlier_score_mat),
                              ncol = ncol(outlier_score_mat),
                              dimnames = dimnames(outlier_score_mat))

  outlier_scores <- Seurat::CreateAssayObject(data = Matrix::Matrix(outlier_score_mat,sparse = TRUE)  )
  outlier_scores@misc$classification <- biased_mat_raw

  seurat_obj[["bias_scores"]] <- outlier_scores
  seurat_obj@assays$bias_scores@scale.data <- outlier_loss_mat

  return(seurat_obj)
}


#' Individual lineage and cluster for preferentially biasing
#' @description
#' This function finds the skewness for a cluster to handle which is an outlier vs inlier
#'
#' @param data_subset This is the raw data subsetted on the initial cluster we want to investigate.
#' @param lineages_to_compare this is the data columns to look into, we use 3 of them: lineage fate, and then membership scores for current population and lineage population
#'
#' @return this returns a dataframe for each input cell with it's nominal outlier score, bias score, and hard cutoff value for biased.
#' @export
#' @import e1071
#' @import stats
#'
#'
#' @examples \dontrun{poised_vector<-ClassifyPoisedCellsSingleLineage(data_subset = data_input%>%
#' filter(.data[[state_grouping_column_name]]==as.character(unique(data_input[[state_grouping_column_name]]))[cluster_iter]),
#' c(lineages_to_compare[lineage_iter],
#'  paste0(current_group,"_membership"),
#'  paste0(lineages_to_compare[lineage_iter],"_membership")))
#' }
ClassifyPoisedCellsSingleLineage<-function(data_subset,lineages_to_compare){

  left_ecdf_func <- function(data) {
    Length <- length(data)
    sorted <- sort(data)
    ecdf <- rep(0, Length)
    for (i in 1:Length) {
      ecdf[i] <- sum(sorted <= data[i]) / Length
    }
    return(ecdf)
  }
  right_ecdf_func <- function(data) {
    Length <- length(data)
    sorted <- sort(data)
    ecdf <- rep(0, Length)
    for (i in 1:Length) {
      ecdf[i] <- sum(sorted >= data[i]) / Length
    }
    return(ecdf)
  }
  measures<-data_subset[,lineages_to_compare]%>%as.matrix
  measures[is.na(measures)]<-0 #This handles missing data calls...
  outlier_score<-matrix(NA,nrow=dim(measures)[1],ncol = dim(measures)[2])
  skew_meas <- apply(measures,2,e1071::skewness)
  print(skew_meas)

  skew_meas[is.nan(skew_meas)]<-0
  skew_meas[is.na(skew_meas)]<-0

  for(pro_iter in 1:dim(measures)[2]){
    right_ecdf<-right_ecdf_func(measures[,pro_iter])
    left_ecdf<-left_ecdf_func(measures[,pro_iter])
    for(cell_iter in 1:dim(measures)[1]){

      if(pro_iter==1){ # this get's lineage while other 2 do membership
        outlier_score[cell_iter,pro_iter]<-right_ecdf[cell_iter]
      }
      else if(skew_meas[pro_iter]>0){
        outlier_score[cell_iter,pro_iter]<-right_ecdf[cell_iter]
      }
      else{
        outlier_score[cell_iter,pro_iter]<-left_ecdf[cell_iter]
      }
    }
  }
  # this makes sure the lineage we care about is not over dominated by membership.
  skew_meas<-ifelse(abs(skew_meas[1])<=abs(skew_meas),sign(skew_meas)*skew_meas[1],skew_meas)
  # this heuristic does the skewing and combination of measures.
  outlier_score_neglog<-(-1)*log(outlier_score)%*%(abs(skew_meas)/norm(skew_meas,type="2"))

  # ok let's now get the distribution of outlier scores and figure out how they are robust
  mad_value <- mad(outlier_score_neglog)

  # Set delta as a multiple of MAD (you can adjust the multiplier as needed)
  delta <- 1.5 * mad_value

  # Huber loss function
  huber_loss <- function(x, delta) {
    abs_x <- abs(x)
    ifelse(abs_x <= delta, 0.5 * (x^2), delta * (abs_x - 0.5 * delta))
  }

  # Calculate Huber loss for each data point
  losses <- huber_loss(outlier_score_neglog - median(outlier_score_neglog), delta)

  # Identify outliers based on a threshold (e.g., 95th percentile)
  #outlier_cutoff_val <- quantile(losses, 0.95)
  outlier_cutoff_val<-quantile(losses)[4]+1.5*IQR(losses)

  #outliers <- outlier_score_neglog[losses > outlier_cutoff_val]

  #outlier_score_neglog<-(-1)*log(outlier_score)%*%(abs(skew_meas)/norm(skew_meas,type="2"))
  #outlier_cutoff_val<-quantile(outlier_score_neglog)[4]+1.5*IQR(outlier_score_neglog)


  data_subset$Outlier_Score<-outlier_score_neglog
  data_subset$Outlier_Score_losses<-losses
  data_subset$Biased<-ifelse(losses>=outlier_cutoff_val,"biased","not")
  return(data_subset)
}
