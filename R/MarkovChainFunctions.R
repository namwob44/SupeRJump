
#' Getting the Fate Matrix as an Assay
#' @description
#' The main wrapper function to place Fates into their own assay, with scaled counts as lineages. The lineage Fates will contain batch corrected values.
#' The Visitation probability will be batch corrected and designated as a graph in the seurat object
#' The CMFPT takes a long time to run, and most analysis does not require this, however, to run the CMFPT measure please run the GetConditionalMeanFirstPassageTime,
#' Once visitation and CMFPT graphs are identified one can run GetWeightedDestinationTime to obtain the WDT which is also stored as a graph.
#' @param seurat_obj the seurat object with a TPM assay and metadata containing a column for absorbing states
#' @param absorbing_state_column_name T This is a metadata column that contains TRUE if it is absorbing state, NA is set for transient cells, used to ensure we get correct names to aggregate
#' @param state_grouping_column_name the column name to be the states you determine your model by, often it is cell type metadata in seurat.
#' @param pseudotime_column_name the column name to be the pseudotime, often it is Combined_Ordering.
#' @param batch_correction_column_name the column that contains sample/mouse/patient differences.
#' @param verbose informs which batch we are currently working.

#'
#' @return It will return a seurat object with updated assays, and graphs.
#' @export
#'
#' @examples seurat_obj<-GetFateMatrixAndMetrics(seurat_obj,state_grouping_column_name="custom_cell_classes_fine",pseudotime_column_name="Combined_Ordering",batch_correction_column_name="sample",absorbing_state_column_name="absorbing_states")
GetFateMatrixAndMetrics <-function(seurat_obj,state_grouping_column_name,pseudotime_column_name,batch_correction_column_name,absorbing_state_column_name,verbose=T){


  all_cells <- colnames(seurat_obj)
  n_global  <- length(all_cells)

  meta <- Seurat::FetchData(seurat_obj,vars = c(state_grouping_column_name,pseudotime_column_name,batch_correction_column_name,absorbing_state_column_name))
  batches <- unique(meta[[batch_correction_column_name]])

  I_all_fate <- list()
  J_all_fate <- list()
  F_all <- list()
  lineage_mats <- list()
  I_all_vp <- list()
  J_all_vp <- list()
  V_all <- list()
  batch_cell_count<-list()
  counter <- 1

  for (b in batches) {
    if (verbose) {
      message("Processing ", batch_correction_column_name, ": ", b)
    }
   # we need to subset seurat objects here and then we need to pump into each Fate matrix and stuff into a list. Once it is placed in the seurat object, let's batch correct

    batch_cells <- rownames(meta)[meta[[batch_correction_column_name]] == b]

    sparse_fate_list<-GetFateMatrix(subset(seurat_obj,cells = batch_cells),absorbing_state_column_name)
    lineage_fate_matrix<-GetCellLineageFateMatrix(sparse_fate_list[["Bmat"]],seurat_obj,absorbing_state_column_name,state_grouping_column_name)
    vistation_metrics<-GetTransientStateMeasures(sparse_fate_list[["Nmat"]])

    # Let's start with Fates graph
    sp_fate <- summary(Matrix::Matrix(sparse_fate_list[["Bmat"]], sparse = TRUE))   # gives i, j, x (1-based)
    fate_source_cells <- rownames(Matrix::Matrix(sparse_fate_list[["Bmat"]], sparse = TRUE))[sp_fate$i]
    fate_dest_cells <- colnames(Matrix::Matrix(sparse_fate_list[["Bmat"]], sparse = TRUE))[sp_fate$j]
    I_all_fate[[counter]] <- match(fate_source_cells,all_cells)
    J_all_fate[[counter]] <- match(fate_dest_cells,all_cells)
    F_all[[counter]] <- sp_fate$x

    # Let's now do Visitation graph

    sparse_vp <- summary(Matrix::Matrix(vistation_metrics[["probability_of_visiting_any_state_before_absorbed"]], sparse = TRUE))

    VP_source_cells <- rownames(Matrix::Matrix(vistation_metrics[["probability_of_visiting_any_state_before_absorbed"]], sparse = TRUE))[sparse_vp$i]
    VP_dest_cells <- colnames(Matrix::Matrix(vistation_metrics[["probability_of_visiting_any_state_before_absorbed"]], sparse = TRUE))[sparse_vp$j]

    I_all_vp[[counter]] <- match(VP_source_cells, all_cells)
    J_all_vp[[counter]] <- match(VP_dest_cells, all_cells)
    V_all[[counter]] <- sparse_vp$x

    # Lineage specific matrices
    # Seurat assay is going to want lineages by cells so transpose this matrix
    lineage_mats[[as.character(b)]]<-t(lineage_fate_matrix)


    batch_cell_count[[counter]]<-batch_cells

    counter <- counter + 1

  }
  I_all_fate <- unlist(I_all_fate, use.names = FALSE)
  J_all_fate <- unlist(J_all_fate, use.names = FALSE)
  F_all <- unlist(F_all, use.names = FALSE)
  batch_cell_count<-unlist(batch_cell_count,use.names = F)
  Fate_sparse <- Matrix::sparseMatrix(i = I_all_fate,j = J_all_fate,x = F_all,dims = c(n_global, n_global))
  rownames(Fate_sparse)<-all_cells
  colnames(Fate_sparse)<- batch_cell_count #this feels incorrect to me.
  seurat_obj@misc[["fates_individual"]] <- Fate_sparse

  # I want "Fates" to be assay name.
  # I want data layer to be sp_fate
  # I want scale.data layer to be sp_lineage
  all_lineages <- unique(unlist(lapply(lineage_mats, rownames)))
  lineage_global <- Matrix::Matrix(0,nrow = length(all_lineages),ncol = n_global,sparse = TRUE)
  rownames(lineage_global) <- all_lineages
  colnames(lineage_global) <- all_cells

  for (m in lineage_mats) {
    lineage_global[rownames(m), colnames(m)] <- m
  }
  lineage_corrected<-BatchCorrectLineages(lineage_global,meta,state_grouping_column_name,batch_correction_column_name)


  fate_assay<-Seurat::CreateAssayObject(data = lineage_global,)
  seurat_obj[["lineage_fates"]] <- fate_assay
  seurat_obj@assays$lineage_fates@scale.data <- lineage_corrected

  # I want visitation probability to be a graph
  I_all_vp <- unlist(I_all_vp, use.names = FALSE)
  J_all_vp <- unlist(J_all_vp, use.names = FALSE)
  V_all <- unlist(V_all, use.names = FALSE)

  VP_sparse <- Matrix::sparseMatrix(i = I_all_vp,j = J_all_vp,x = V_all,dims = c(n_global, n_global))

  rownames(VP_sparse) <- all_cells
  colnames(VP_sparse) <- all_cells
  VP_sparse<-BatchCorrectGraph(VP_sparse,meta,state_grouping_column_name,batch_correction_column_name)

  seurat_obj@graphs[["visitation"]] <- VP_sparse
  # I want to apply mean matching on Visitation graph though. This is done by mean matching the factors of batch_correction_column_name and state_grouping_column_name which we can get from variable meta


  return(seurat_obj)
}




#' Calculate Fate and Fundamental Matrices
#' @description
#' This function is the main workhorse we use absorbing markov chains and partition the TPM  into a Q (transient matrix) and R (absorbing matrix)
#' Afterwards, we solve for the fundamental matrix which is N = (I-Q)^{-1} and find the Fate matrix through B=NR
#'
#' @param seurat_obj The seurat object, please check you have TPM as a graph name.
#' @param absorbing_state_column_name This is a metadata column that contains TRUE if it is absorbing state, NA is set for transient cells
#' @useDynLib SupeRJump, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom Seurat FetchData
#' @importFrom Matrix rowSums
#' @return This will return a list of 2 matrices, a fate matrix B, and a fundamental matrix N.
#' @export
#'
#'
#' @examples sparse_fate_matrices<-GetFateMatrix(seurat_obj,absorbing_state_column_name="absorbing_state_bool")
#'
GetFateMatrix<-function(seurat_obj,absorbing_state_column_name){

  TPM_matrix<-as.matrix(seurat_obj@graphs[["TPM"]])
  sink_cells <- Seurat::FetchData(seurat_obj,vars=c(absorbing_state_column_name))
  sink_cell_names<-rownames(sink_cells)[which(!is.na(sink_cells[[absorbing_state_column_name]]))]

  TPM_matrix<-ifelse(is.na(TPM_matrix),0,TPM_matrix)
  TPM_matrix<-ifelse(is.nan(TPM_matrix),0,TPM_matrix)

  TPM_matrix<-TPM_matrix/Matrix::rowSums(TPM_matrix)
  # check columns have more than 1 entry
  TPM_matrix<-removeIllConditionedCells(TPM_matrix)

  Rmat<-(TPM_matrix)[-which(rownames((TPM_matrix))%in%sink_cell_names),which(colnames((TPM_matrix))%in%sink_cell_names)]%>%as.matrix # need this matrix conversion to handle single absorbing state
  Qmat<-(TPM_matrix)[-which(rownames((TPM_matrix))%in%sink_cell_names),-which(colnames((TPM_matrix))%in%sink_cell_names)]

  Q_solved<-eigenMatInverse(diag(nrow=dim(Qmat)[1])-Qmat) #this goes vroom.
  colnames(Rmat)<-colnames(TPM_matrix)[which(colnames((TPM_matrix))%in%sink_cell_names)] # need this in case they only have a single absorbing state?
  Bmat<-matrix(data=0,nrow=nrow((TPM_matrix)),ncol=dim(Rmat)[2])
  rownames(Bmat)<-rownames((TPM_matrix))
  colnames(Bmat)<-colnames(Rmat)
  Bmat[-which(rownames((TPM_matrix))%in%sink_cell_names),]<-eigenMapMatMult(Q_solved,Rmat)
  if(ncol(Bmat)>1){
    diag(Bmat[which(rownames((TPM_matrix))%in%sink_cell_names),])<-1
  }
  rownames(Q_solved)<-colnames(Q_solved)<-rownames(Qmat)#rownames(TPM_matrix)[-which(rownames((TPM_matrix))%in%sink_cells$row_names)]
  return(list("Bmat"=Bmat,
              "Nmat"=Q_solved))
}


#' Aggregating the Fates from absorbing states to lineages
#' @description
#' The goal of this function is to aggregate the absorbing states by their state columns.
#'
#' @param Bmat This is the fate matrix that was found using GetFateMatrix
#' @param seurat_obj Used so we can extract the two metadata columns below cleanly.
#' @param absorbing_state_column_name T This is a metadata column that contains TRUE if it is absorbing state, NA is set for transient cells, used to ensure we get correct names to aggregate
#' @param state_grouping_column_name This is the grouping function to aggregate the absorbing probabilities by
#'
#' @return It returns a dense matrix of n_cells by l_lineages
#' @export
#' @import dplyr
#' @import Seurat
#'
#' @examples lineage_fate_matrix<-GetCellLineageFateMatrix(Bmat,seurat_obj,absorbing_state_column_name="absorbing_state_bool",state_grouping_column_name="custom_cell_classes_fine")

GetCellLineageFateMatrix<-function(Bmat,seurat_obj,absorbing_state_column_name,state_grouping_column_name){

    sink_cells <- Seurat::FetchData(seurat_obj,vars=c(absorbing_state_column_name,state_grouping_column_name))
    sink_cells <- sink_cells%>%
      dplyr::filter(!is.na(.data[[absorbing_state_column_name]]))%>%
      tibble::rownames_to_column(var="row_names")
    sink_groups<-c(as.character(unique(sink_cells[[state_grouping_column_name]])))

  Bmat_condensed<-matrix(data=0,nrow=nrow(Bmat),ncol=length(sink_groups))
  rownames(Bmat_condensed)<-rownames(Bmat)
  colnames(Bmat_condensed)<-sink_groups
  for(sink_iter in 1:length(sink_groups)){
    lineage_columns<-sink_cells%>%
      dplyr::filter(.data[[state_grouping_column_name]]==sink_groups[sink_iter])%>%
      dplyr::pull(row_names)
    if(length(lineage_columns)>1){
      Bmat_condensed[,sink_iter]<-Matrix::rowSums(Bmat[,colnames(Bmat)%in%lineage_columns],na.rm=T)
    }else{
      Bmat_condensed[,sink_iter]<-Bmat[,colnames(Bmat)%in%lineage_columns]
    }
  }
  return(Bmat_condensed)
}



#' Get Membership of cells through gaussian estimation
#'
#' @param seurat_obj this is the input seurat object to incorporate. It needs the column group to function
#' @param state_grouping_column_name the column name to be the states you determine your model by, often it is cell type metadata in seurat.
#'
#' @return It returns a seurat object assay of membership assignments that are normalized per group
#' @export
#' @import dplyr
#' @import mvtnorm
#' @import Matrix
#' @import Seurat
#' @import tibble
#'
#' @examples \dontrun{seurat_obj<-GetMembership(seurat_obj,state_grouping_column_name="custom_cell_classes_fine")}
GetMembership<-function(seurat_obj,state_grouping_column_name){
  raw_data_to_test<-Seurat::Embeddings(seurat_obj, "pca") %>% as.data.frame # cells x PCs
  added_column <- Seurat::FetchData(seurat_obj,vars = c(state_grouping_column_name))
  raw_data_to_test<- raw_data_to_test%>%
    tibble::rownames_to_column(var = "row_names")%>%
    dplyr::inner_join(added_column%>%
                        tibble::rownames_to_column(var="row_names"),
                      by="row_names")%>%
    tibble::column_to_rownames(var="row_names")


  Mean_mus_clust<-raw_data_to_test%>%
    dplyr::group_by(.data[[state_grouping_column_name]])%>%
    dplyr::summarise( dplyr::across( dplyr::everything(), mean))

  var_mus_clust<-raw_data_to_test%>%
    dplyr::group_by(.data[[state_grouping_column_name]])%>%
    dplyr::summarise( dplyr::across( dplyr::everything(), var))
  membership_raw = matrix(data=NA_real_,nrow=dim(raw_data_to_test)[1],ncol =length(unique(raw_data_to_test[[state_grouping_column_name]])))
  for(cluster_it in 1:length(unique(raw_data_to_test[[state_grouping_column_name]]))){
    membership_raw[,cluster_it]<-mvtnorm::dmvnorm(raw_data_to_test%>%
                                                    dplyr::select(-.data[[state_grouping_column_name]])%>%
                                                    as.matrix,
                                                  mean = as.matrix(Mean_mus_clust[cluster_it,2:length(Mean_mus_clust)]),
                                                  sigma = diag(as.vector(var_mus_clust[cluster_it,2:length(Mean_mus_clust)])))
  }
  colnames(membership_raw)<-Mean_mus_clust[[state_grouping_column_name]]
  rownames(membership_raw)<-rownames(raw_data_to_test)
  membership_norm <-membership_raw/Matrix::rowSums(membership_raw)
  module_assay <- Seurat::CreateAssayObject(data = t(membership_norm))
  seurat_obj[["Membership"]] <- module_assay
  return(seurat_obj)
}




#' Fate Matrix scaled by cell clusters
#' @description
#' This is a utility function for future development. It allows for a metric called committment
#' @param seurat_obj the seurat object with a TPM assay and metadata containing a column for absorbing states
#'
#' @return this returns a matrix that is cells by cluster.
#' @export
#' @import Matrix
#' @useDynLib SupeRJump, .registration=TRUE
#' @importFrom Rcpp evalCpp
#'
#' @examples  \dontrun{cell_to_clust_commitment<-FateMatrixByClusterMembership(seurat_obj)}
FateMatrixByClusterMembership<-function(seurat_obj){

  TPM_matrix<-as.matrix(seurat_obj@graphs[["TPM"]])
  memberships_by_cluster<- seurat_obj@assays[["Membership"]]
  TPM_matrix<-ifelse(is.nan(TPM_matrix),0,TPM_matrix)
  TPM_matrix<-ifelse(is.na(TPM_matrix),0,TPM_matrix)

  temp_TPM<-TPM_matrix
  temp_TPM<-temp_TPM/rowSums(temp_TPM)

  ill_condition_check<-ifelse(TPM_matrix!=0,1,0)
  print(which(Matrix::colSums(ill_condition_check)<=1))
  if(length(which(Matrix::colSums(ill_condition_check)<=1))>0){
    TPM_matrix<-TPM_matrix[-which(rownames((TPM_matrix))%in%colnames(ill_condition_check)[which(Matrix::colSums(ill_condition_check)<=1)]),-which(colnames((TPM_matrix))%in%colnames(ill_condition_check)[which(Matrix::colSums(ill_condition_check)<=1)])]
  }
  Qmat<-TPM_matrix

  Q_solved<-eigenMatInverse(diag(nrow=dim(Qmat)[1])-Qmat)#this shit goes vroom if it average number of visits to states

  memberships_by_cluster<-memberships_by_cluster[rownames(Qmat),]

  Cell_To_clust_commitment<-eigenMapMatMult(Q_solved,memberships_by_cluster)
  rownames(Cell_To_clust_commitment)<-rownames(Qmat)
  colnames(Cell_To_clust_commitment)<-colnames(memberships_by_cluster)
  return(Cell_To_clust_commitment)

}




#' Turning cell by lineage matrix, to a cluster by lineage matrix
#' @description
#' This is a utility function for future use.
#'
#' @param seurat_obj the seurat object with the Fates assay.
#' @param state_grouping_column_name the clusters to summarize and aggregate by
#'
#' @return It produces a cluster by cluster matrix by simply aggregating all the cell states in rows and columns
#' @export
#' @import Seurat
#' @import dplyr
#' @import tibble
#'
#'
#' @examples \dontrun{small_heatmap<-ClusterToClusterFateMatrix(seurat_obj,state_grouping_column_name="custom_cell_classes_fine")
#' pheatmap::pheatmap(small_heatmap,scale="row")
#' }
ClusterToClusterFateMatrix<-function(seurat_obj,state_grouping_column_name){

  Bmat <- Seurat::GetAssayData(seurat_obj,assay = "Fates",layer = "Lineages")
  cosine_data<-Seurat::FetchData(seurat_obj,vars=c(state_grouping_column_name))

  # WE. CAN NORMALIZE AND DETERMINE CUSTOM THINGS HERE!!
  sink_groups <-colnames(Bmat)
  Fate_matrix <-matrix(data=NA,nrow=length(unique(cosine_data[[state_grouping_column_name]])),ncol=length(sink_groups))
  rownames(Fate_matrix)<-as.character(unique(cosine_data[[state_grouping_column_name]]))
  colnames(Fate_matrix)<-colnames(Bmat)

  for(current_lineage_iter in 1:length(unique(cosine_data[[state_grouping_column_name]]))){
    current_source_cells<-cosine_data%>%
      dplyr::filter(.data[[state_grouping_column_name]]==as.character(unique(cosine_data[[state_grouping_column_name]])[current_lineage_iter]))%>%
      tibble::rownames_to_column(var="row_names")%>%
      dplyr::pull(row_names)
    for(sink_lineage_iter in 1:length(sink_groups)){
      Fate_matrix[current_lineage_iter,sink_lineage_iter]<-sum(Bmat[rownames(Bmat)%in%current_source_cells,sink_lineage_iter])
    }
  }
  Fate_matrix_norm<-Fate_matrix#apply(Fate_matrix,2,function(x){x/sum(x)})
  return(Fate_matrix_norm)
}

#' Get Visitation Probability and Expected timesteps for TPM
#' @description
#' This function obtains metrics for transient to transient state networks. It contains 3 measures for now. The first is the expected time to
#'
#' @param Nmat this is the fundamental matrix found before
#' @useDynLib SupeRJump, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @return this returns a list of measures including timesteps to any absorbing state and visitation probability.
#' @export
#'
#' @examples  \dontrun{vistation_metrics<-GetTransientStateMeasures(sparse_fate_list[["Nmat"]])}
#'
GetTransientStateMeasures<-function(Nmat){
  print("Working on Expected Time")
  expected_time_cheat<-eigenMapMatMult(Nmat,matrix(data=1,nrow=nrow(Nmat),ncol=1))
  colnames(expected_time_cheat)<-"expected_time"
  #variance_on_timesteps_until_absorbed<- (((2*Nmat)-(diag(nrow=dim(Nmat)[1])))%*%expected_time_cheat) -(expected_time_cheat*expected_time_cheat)
  #colnames(variance_on_timesteps_until_absorbed)<-"variance_timestep"
  print("Done with Expected time")
  print("Working on Visitation Probability")
  probability_transient_state_Transfer<-eigenMapMatMult((Nmat-diag(nrow=dim(Nmat)[1])),eigenMatInverse(diag(diag(Nmat))))
  rownames(probability_transient_state_Transfer)<-rownames(Nmat)
  colnames(probability_transient_state_Transfer)<-colnames(Nmat)
  print("Finished Visitation Probability")

  # var_time_steps_cell2cell<-eigenMapMatMult(Nmat,2*diag(diag(Nmat))-diag(nrow=dim(Nmat)[1]))-(Nmat*Nmat)

  list_of_measures<-list("timesteps_until_absorbed_by_any_sink"=expected_time_cheat,
                         #  "variance_on_timesteps_until_absorbed"=variance_on_timesteps_until_absorbed,
                         # "variance_on_visits_before_being_absorbed"=var_time_steps_cell2cell,
                         "probability_of_visiting_any_state_before_absorbed"=probability_transient_state_Transfer)
  return(list_of_measures)

}


#' Batch Correct Lineages
#' @description
#' This function performs mean matching batch correction on lineages in a similar fashion as ComBAT.
#'
#' @param matrix_to_correct This can either be the Lineage fates or Visitation Matrix
#' @param metadata This is a dataframe with metadata for the state grouping and batches
#' @param state_grouping_column_name  the clusters to summarize and aggregate by
#' @param batch_correction_column_name the column name to determine how the matrices are subsetted.
#' @param verbose updating where in batch correction process the fit is.
#'
#' @return This returns a lineages by cells matrix with batch corrected values renormalized between 0 and 1
#' @export
#' @import Matrix
#'
#' @examples  \dontrun{lineage_corrected<-BatchCorrectLineages(lineages_matrix,meta,state_grouping_column_name="custom_cell_classes_fine",batch_correction_column_name="sample")}
BatchCorrectLineages<-function(matrix_to_correct,metadata,state_grouping_column_name,batch_correction_column_name,verbose=T){

  prob <- as.matrix(matrix_to_correct)
  adjusted <- prob

  #metadata <- metadata[common_cells, , drop = FALSE]
  lineages <- rownames(prob)

  states <- unique(metadata[[state_grouping_column_name]])
  states <- states[!is.na(states)]

  batches <- unique(metadata[[batch_correction_column_name]])
  batches <- batches[!is.na(batches)]
  for (s in states) {
    state_cells <- rownames(metadata)[metadata[[state_grouping_column_name]] == s]
    state_cells <- intersect(state_cells, colnames(prob))
    if (length(state_cells) == 0) {
      next
    }
    if (verbose) {
      message("Mean matching state: ", s)
    }
    # mean(probability | cluster_init, cluster_to)
    # vector of length n_lineages
    state_lineage_mean <- Matrix::rowMeans(prob[, state_cells, drop = FALSE],na.rm = TRUE)
    for (b in batches) {
      batch_state_cells <- rownames(metadata)[metadata[[state_grouping_column_name]] == s &metadata[[batch_correction_column_name]] == b]
      batch_state_cells <- intersect(batch_state_cells, colnames(prob))
      if (length(batch_state_cells) == 0){
        next
      }
      # mean(probability | cluster_init, cluster_to, batch)
      # vector of length n_lineages
      batch_state_lineage_mean <- Matrix::rowMeans(prob[, batch_state_cells, drop = FALSE],na.rm = TRUE)
      # offset to remove batch effect
      adjustment <- state_lineage_mean - batch_state_lineage_mean
      adjusted[, batch_state_cells] <- sweep(adjusted[, batch_state_cells, drop = FALSE],MARGIN = 1,STATS = adjustment,FUN = "+")
    }
  }
  min_val <- min(adjusted, na.rm = TRUE)
  max_val <- max(adjusted, na.rm = TRUE)
  if (is.finite(min_val) && is.finite(max_val) && max_val > min_val) {
    adjusted <- (adjusted - min_val) / (max_val - min_val)
  } else {
      warning("Global min-max scaling skipped because max <= min.")
    }

  return(adjusted)
}



#' Wrapper for Batch Corrected Conditional Mean First Passage Time (CMFPT)
#' @description
#' This function is the wrapper to the log10 of Conditional Mean First Passage Time (CMFPT). It will perform CMFPT for each batch sequentially. This function will take a long time, so it is advised to
#' For reference, on an M2 Macbook air with 8 threads took 21 hours for 7500 cells. More bench marking will be conducted to improve the speed of this calculation.
#' The Conditional Mean First Passage Time (CMFPT) is a measure that answers how quickly cell i(row) will visit cell j(column) on average before ending up in an absrobing state.
#' The measure is the number of time steps. It is important to batch correct as this is a network size dependent measure, it is ill-advised to compare head to head.

#' @param seurat_obj The seurat object for each batch to determine the CMFPT
#' @param state_grouping_column_name the column name for the states, this is used for batch correction.
#' @param batch_correction_column_name the column name to determine how to subset cell-cell networks with TPMs
#' @param old_flag temporary flag for comparing new methods with one that is slow but mathematically exact.
#' @param targets_in_rows This should be left as TRUE, but optional control for futture downstream changes,
#' @param n_threads left as 1 for safety, consider setting higher based on computer hardware. 8 on an M2 Macbook air took 21hours on about 7500 cells. Further parallel strategies to improve speed will eventually be coming soon.
#' @param verbose Informs whether to update which sample/batch we are currently working on.
#'
#' @return should return a graphs object called CMFPT in the seurat object
#' @export
#' @import Matrix
#' @import Seurat
#'
#' @examples  \dontrun{seurat_obj<-GetBatchCorrectedCMFPT(seurat_obj,batch_correction_column_name="sample",old_flag = F,targets_in_rows = T,n_threads = 8,verbose=T)}

GetBatchCorrectedCMFPT<-function(seurat_obj,state_grouping_column_name,batch_correction_column_name,old_flag = F,targets_in_rows = T,n_threads = 1,verbose=T){

  all_cells <- colnames(seurat_obj)
  n_global  <- length(all_cells)

  meta <- Seurat::FetchData(seurat_obj,vars = c(state_grouping_column_name,batch_correction_column_name))

  batches <- unique(meta[[batch_correction_column_name]])

  # global sparse buffers
  I_all <- list()
  J_all <- list()
  V_all <- list()
  counter <- 1

  for (b in batches) {
    if (verbose) {
      message("Processing ", batch_correction_column_name, ": ", b)
    }
    # ---- batch cells ----
    batch_cells <- rownames(meta)[meta[[batch_correction_column_name]] == b]
    # ---- compute batch TPM ----
    Sparse_CMFPT <- GetConditionalMeanFirstPassageTime(subset(seurat_obj, cells = batch_cells),old_flag=old_flag,targets_in_rows=targets_in_rows,n_threads=n_threads)
    # ---- map local → global indices ----

    # ---- extract sparse entries ----
    sp <- summary(Matrix::Matrix(log10(Sparse_CMFPT+1),sparse = TRUE))   # gives i, j, x (1-based)
    CMFPT_source_cells <- rownames(Matrix::Matrix(Sparse_CMFPT, sparse = TRUE))[sp$i]
    CMFPT_dest_cells <- colnames(Matrix::Matrix(Sparse_CMFPT, sparse = TRUE))[sp$j]


    # map to global indices
    I_all[[counter]] <- match(CMFPT_source_cells, all_cells)
    J_all[[counter]] <-  match(CMFPT_dest_cells, all_cells)
    V_all[[counter]] <- sp$x

    counter <- counter + 1
  }

  # ---- combine all batches ----
  I_all <- unlist(I_all, use.names = FALSE)
  J_all <- unlist(J_all, use.names = FALSE)
  V_all <- unlist(V_all, use.names = FALSE)

  CMFPT_sparse <- Matrix::sparseMatrix(i = I_all,j = J_all,x = V_all,dims = c(n_global, n_global))

  rownames(CMFPT_sparse) <- all_cells
  colnames(CMFPT_sparse) <- all_cells
  print(CMFPT_sparse)

  #return(CMFPT_sparse)

  CMFPT_sparse<-BatchCorrectGraph(CMFPT_sparse,meta,state_grouping_column_name,batch_correction_column_name)

  seurat_obj@graphs[["CMFPT"]] <- CMFPT_sparse

  return(seurat_obj)

}


#' Get Conditional Mean First Passage Time in R parallel strategy
#' @description
#' The Conditional Mean First Passage Time (CMFPT) is a measure that answers how quickly cell i(row) will visit cell j(column) on average before ending up in an absrobing state.
#' The measure is the number of time steps. It is important to batch correct as this is a network size dependent measure, it is ill-advised to compare head to head.
#'
#'
#' @param seurat_obj The seurat object for each batch to determine the CMFPT
#' @param old_flag temporary flag for comparing new methods with one that is slow but mathematically exact.
#' @param targets_in_rows This should be left as TRUE, but optional control for futture downstream changes,
#' @param n_threads left as 1 for safety, consider setting higher based on computer hardware. 8 on an M2 Macbook air took 21hours on about 7500 cells. Further parallel strategies to improve speed will eventually be coming soon.
#' @useDynLib SupeRJump, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @import doParallel
#' @import foreach
#' @import Matrix
#' @return This will return a single square matrix containing how long on average it takes for cell i(row) to visit cell j (column)
#' @export
#'
#' @examples  \dontrun{Sparse_CMFPT <- GetConditionalMeanFirstPassageTime(subset(seurat_obj, cells = batch_cells),old_flag=old_flag,targets_in_rows=targets_in_rows,n_threads=n_threads)}

GetConditionalMeanFirstPassageTime<-function(seurat_obj,old_flag=F,targets_in_rows=T,n_threads=1){
  TPM_matrix <-as.matrix(seurat_obj@graphs[["TPM"]])
  # I don't know if subset cells truncates both rows and columns. So, let's check, but I think rows should get axed for sure.
  # Looks like we are ok, but leaving here just in case for some reason insanity hits.
  if(nrow(TPM_matrix)!=ncol(TPM_matrix)){
    TPM_matrix<- TPM_matrix[rownames(TPM_matrix),rownames(TPM_matrix)]
  }


  if (old_flag==T){
    TPM_matrix<-ifelse(is.nan(TPM_matrix),0,TPM_matrix)
    TPM_matrix<-ifelse(is.na(TPM_matrix),0,TPM_matrix)

    temp_TPM<-TPM_matrix
    temp_TPM<-temp_TPM/rowSums(temp_TPM)

    doParallel::registerDoParallel(cores = 8)
    TimeMatrix<-foreach::foreach(current_cell_iter=1:(dim(TPM_matrix)[1]),.combine='rbind')%dopar%{
      temp_vec <-matrix(0,nrow=1,ncol=dim(TPM_matrix)[1])
      Qmat<-temp_TPM[-current_cell_iter,-current_cell_iter]
      Rmat<-temp_TPM[-current_cell_iter,current_cell_iter]%>%as.matrix
      #
      #   # Q_solved<-eigenMatInverse(diag(nrow=dim(Qmat)[1])-Qmat)#this shit goes vroom.
      #   #
      #   # Bmat<-eigenMapMatMult(Q_solved,Rmat)
      #   # diag_mat<-diag(c(Bmat))
      #   # temp_multiply_matrix<-eigenMapMatMult(eigenMapMatMult(eigenMatInverse(diag_mat),Q_solved),diag_mat)
      #
      temp_vec[,-current_cell_iter]<-eigenOutputFlightTime(Qmat,Rmat)
      return(temp_vec)
    }
    doParallel::stopImplicitCluster()
    rownames(TimeMatrix)<-rownames(TPM_matrix)
    colnames(TimeMatrix)<-colnames(TPM_matrix)
    TimeMatrix<-t(TimeMatrix)
    return(TimeMatrix)
  }
  else{
    TPM_matrix <- as(TPM_matrix, "dgCMatrix")

    bad <- !is.finite(TPM_matrix@x)
    if (any(bad)) {
      TPM_matrix@x[bad] <- 0
    }

    rs <- Matrix::rowSums(TPM_matrix)
    inv_rs <- numeric(length(rs))
    nz <- rs > 0
    inv_rs[nz] <- 1 / rs[nz]

    # Row-stochastic sparse transition matrix
    P <- Matrix::Diagonal(x = inv_rs) %*% TPM_matrix

    out <- flighttime_sparse_targets_cpp(
      P,
      targets = seq_len(nrow(P)),
      targets_in_rows = targets_in_rows,
      hit_tol = 1e-14,
      n_threads = n_threads
    )

    rn <- rownames(TPM_matrix)
    cn <- colnames(TPM_matrix)

    if (!is.null(rn)) {
      if (targets_in_rows) {
        rownames(out) <- rn
        colnames(out) <- cn
      } else {
        rownames(out) <- rn
        colnames(out) <- cn
      }
    }

    out
  }
}



#' Batch Correct Graph
#' @description
#' This function is used to batch correct the input graph either Visitation or CMFPT for each batch. Then it is expanded back out to the graph form.
#'
#' @param graph The graph to be used for the correction, either will be visitation or CMFPT
#' @param meta the meta data columns of the seurat object to access the states and batch corrected columns
#' @param state_grouping_column_name the column name to be the states you determine your model by, often it is cell type metadata in seurat.
#' @param batch_correction_column_name the column that contains sample/mouse/patient differences.
#'
#' @return it returns the graph which should overfit the seurat object.
#' @export
#' @import Matrix
#'
#' @examples   \dontrun{CMFPT_sparse<-BatchCorrectGraph(CMFPT_sparse,meta,state_grouping_column_name,batch_correction_column_name)}
BatchCorrectGraph<-function(graph,meta,state_grouping_column_name,batch_correction_column_name){


    common_source_cells <- intersect(rownames(graph), rownames(meta))
    common_dest_cells <- intersect(colnames(graph), rownames(meta))
    graph <- graph[common_source_cells, common_dest_cells, drop = FALSE]
    source_meta <- meta[rownames(graph), , drop = FALSE]
    dest_meta <- meta[colnames(graph), , drop = FALSE]

    sp <- summary(graph)
    source_cells <- rownames(graph)[sp$i]
    dest_cells <- colnames(graph)[sp$j]

    print(length(source_cells))
    print(length(dest_cells))
    print(length(sp$x))
    print(length(sp$i))
    source_state <- source_meta[source_cells, state_grouping_column_name]
    dest_state <- dest_meta[dest_cells, state_grouping_column_name]
    source_batch <- source_meta[source_cells, batch_correction_column_name]

    triplets <- data.frame(i = sp$i,
      j = sp$j,
      values_to_norm = sp$x,
      source_cell = source_cells,
      dest_cell = dest_cells,
      cluster_init = source_state,
      cluster_to = dest_state,
      init_batch = source_batch,
      stringsAsFactors = FALSE)

    triplets <- triplets[
      !is.na(triplets$cluster_init) &
        !is.na(triplets$cluster_to) &
        !is.na(triplets$init_batch),
      ,
      drop = FALSE]

    reduced_key <- paste(triplets$cluster_init,triplets$cluster_to,sep = "___")

    full_key <- paste(triplets$cluster_init,triplets$cluster_to,triplets$init_batch,sep = "___")

    reduced_mean_lookup <- tapply(triplets$values_to_norm,reduced_key,mean,na.rm = TRUE)

    full_mean_lookup <- tapply(triplets$values_to_norm,full_key,mean,na.rm = TRUE)

    triplets$reduced_mean <- as.numeric(reduced_mean_lookup[reduced_key])
    triplets$full_mean <- as.numeric(full_mean_lookup[full_key])
    triplets$adjusted_values <- triplets$values_to_norm - (triplets$full_mean - triplets$reduced_mean)

    min_val <- min(triplets$adjusted_values, na.rm = TRUE)
    max_val <- max(triplets$adjusted_values, na.rm = TRUE)
    if (is.finite(min_val) && is.finite(max_val) && max_val > min_val) {
      triplets$scaled_values <- (triplets$adjusted_values - min_val) / (max_val - min_val)
    } else {
      warning("Global min-max scaling skipped because max <= min.")
      triplets$scaled_values <- triplets$adjusted_values
    }
    corrected_x <- triplets$scaled_values
    corrected_graph <- Matrix::sparseMatrix(i = triplets$i,j = triplets$j,x = corrected_x,dims = dim(graph))

    rownames(corrected_graph) <- rownames(graph)
    colnames(corrected_graph) <- colnames(graph)

    return(corrected_graph)

}



#' Weighted Destination Time (WDT) function
#' @description
#' This function requires that GetMembership, GetBatchCorrectedCMFPT, and GetFateMatrixAndMetrics were run first.
#' This takes the
#' @param seurat_obj The seurat object for each batch to determine the CMFPT
#' @param preference_lineage_column Which meta data column to appropriately aggregate for a reduced representation, leave NULL if it not interested.
#' @param enrichment_flag Set TRUE if you want the enrichment of the preference lineage column to a null distribution.
#'
#' @return the seurat object with an updated Assay called WDT and 2 new matrices in seurat_obj misc layer
#' @export
#' @import Matrix
#' @import Seurat
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @import tibble
#'
#' @examples  \dontrun{seurat_obj<-GetWeightedDestinationTime(seurat_obj,preference_lineage_column="sgRNA",enrichment_flag=T)}
GetWeightedDestinationTime<-function(seurat_obj,preference_lineage_column=NULL, enrichment_flag=T){


  all_cell_rownames <- intersect(rownames(seurat_obj@graphs[["visitation"]]),rownames(seurat_obj@graphs[["CMFPT"]]))
  all_cell_colnames <- intersect(colnames(seurat_obj@graphs[["visitation"]]),colnames(seurat_obj@graphs[["CMFPT"]]))

  A<-as.matrix(seurat_obj@graphs[["visitation"]][all_cell_rownames,all_cell_colnames])/(as.matrix(seurat_obj@graphs[["CMFPT"]][all_cell_rownames,all_cell_colnames]) +1e-6)
  A_bar <- A/rowSums(A)
  A_bar[!is.finite(A_bar)] <- 0

  memcheck_var <- t(Seurat::GetAssayData(seurat_obj,assay = "Membership"))[all_cell_rownames,]

  S_M <-Matrix::colSums(ifelse(memcheck_var>0,1,0))
  DM_inv <- diag(1/pmax(S_M,1))
  rownames(DM_inv)<-colnames(DM_inv)<-names(S_M)
  WDT_nk <- A_bar %*% memcheck_var[all_cell_rownames,] %*% DM_inv

  WDT_assay<-Seurat::CreateAssayObject(data = t(WDT_nk))
  seurat_obj@assays[["WDT"]] <- WDT_assay

  if(!is.null(preference_lineage_column)){

    meta <- Seurat::FetchData(seurat_obj,vars = c(preference_lineage_column))
    G_mat<-meta%>%
      tibble::rownames_to_column(var="current_state")%>%
      dplyr::select(current_state,.data[[preference_lineage_column]])%>%
      table%>%
      as.data.frame%>%
      tidyr::pivot_wider(id_cols = c(current_state),names_from = preference_lineage_column,values_from = "Freq")%>%
      tibble::column_to_rownames(var="current_state")

    S_L <- Matrix::colSums(G_mat)
    DL_inv <- diag(1/pmax(S_L,1))
    rownames(DL_inv)<-colnames(DL_inv)<-names(S_L)

    WDT_pref_lk <- DL_inv %*% t(G_mat[all_cell_rownames,]) %*% WDT_nk

    WDT_pref_lk_norm<-WDT_pref_lk / rowSums(WDT_pref_lk)
    WDT_pref_lk_norm[!is.finite(WDT_pref_lk_norm)] <- 0

    seurat_obj@misc[["WDT_pref_lineage_by_state"]]<- WDT_pref_lk_norm

    if (enrichment_flag){

      one_vec <- matrix(1, nrow = nrow(A_bar), ncol = 1)
      u_global <- (t(one_vec) %*% A_bar %*% memcheck_var[all_cell_rownames,] %*% DM_inv) / nrow(A_bar)
      p_global <- as.numeric(u_global / sum(u_global))

      E <- log2((WDT_pref_lk_norm + 1e-8) / matrix(p_global, nrow = nrow(WDT_pref_lk_norm), ncol = ncol(WDT_pref_lk_norm), byrow = TRUE))

      seurat_obj@misc[["WDT_pref_enrichment"]]<- E

    }
  }

 return(seurat_obj)
}






#' Check for removing non-connected cells, necessary for Fate and Visitation
#' @description
#' We use this to ensure our graph eliminates singleton cells
#'
#' @param TPM_matrix The transition probability matrix.
#'
#' @return It returns a truncated TPM matrix with the problematic cells eliminated in both rows and columns.
#' @export
#' @import Matrix
#'
#' @examples   \dontrun{TPM_matrix<-removeIllConditionedCells(TPM_matrix)}
removeIllConditionedCells<-function(TPM_matrix){
  ill_condition_check<-ifelse(TPM_matrix!=0,1,0)
  print(which(Matrix::colSums(ill_condition_check)<=1))
  if(length(which(Matrix::colSums(ill_condition_check)<=1))>0){
    TPM_matrix<-TPM_matrix[-which(rownames((TPM_matrix))%in%colnames(ill_condition_check)[which(Matrix::colSums(ill_condition_check)<=1)]),
                           -which(colnames((TPM_matrix))%in%colnames(ill_condition_check)[which(Matrix::colSums(ill_condition_check)<=1)])]
  }
  return(TPM_matrix)
}
