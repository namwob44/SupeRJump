
#' Fitting Jump Model for each Eigen Cell and PC.
#'
#' @param seurat_obj the main seurat object to analyze.
#' @param Total_Eigen_Data  A set of data points for sink points from GetEigenData
#' @param start_point index for which cell is "origin" or earliest cell
#' @param state_grouping_column_name the column name to be the states you determine your model by, often it is cell type metadata in seurat.
#' @param pseudotime_column_name the column name to be the pseudotime, often it is Combined_Ordering.
#' @param batch_correction_column_name the column that contains sample/mouse/patient differences.
#' @param order this is a truncated value to consider the poisson process in the infinite sum.
#' @param maxeval the total number of evaluations for the Negative log-likelihood fit (default is 1000)
#' @param xtol_rel the relative tolerance for the nloptr default is 1e-5
#'
#' @return It returns a list of models for each EigenCell and PC
#' @export
#' @import Seurat
#' @import stats
#'
#' @examples \dontrun{
#' seurat_obj<-CosineTransformPCAData(seurat_obj)
#' Y <- GetSinkData(seurat_obj,state_grouping_column_name="custom_cell_classes_fine",pseudotime_column_name = "Combined_Ordering")
#' model_fits <- Reduced_Jump_Prob_fast(seurat_obj,Y,start_point=1,state_grouping_column_name="custom_cell_classes_fine",pseudotime_column_name="Combined_Ordering",batch_correction_column_name="sample")
#'
#'
#' }
Reduced_Jump_Prob_fast <- function(seurat_obj, Total_Eigen_Data, start_point,
                                   state_grouping_column_name="custom_cell_classes_fine",
                                   pseudotime_column_name="Combined_Ordering",
                                   batch_correction_column_name="sample",
                                   order = 5,
                                   maxeval = 1000,
                                   xtol_rel = 1e-5) {

  # Pull embeddings + metadata (no joins)
  emb <- Seurat::Embeddings(seurat_obj, "transformed_pca")  # cells x PCs
  meta <- Seurat::FetchData(seurat_obj, vars = c(state_grouping_column_name,
                                                 pseudotime_column_name,
                                                 batch_correction_column_name))

  # Align and order by pseudotime once
  meta <- meta[match(rownames(emb), rownames(meta)), , drop = FALSE]
  ord <- order(meta[[pseudotime_column_name]])
  emb <- emb[ord, , drop = FALSE]
  meta <- meta[ord, , drop = FALSE]
  PT  <- as.numeric(meta[[pseudotime_column_name]])

  # Precompute deltat once
  pt0 <- PT[start_point]
  deltat <- abs(PT[-start_point] - pt0)          # length nCells-1
  absdt  <- abs(deltat)

  # Precompute factorials for jump sum
  fac <- factorial(0:order)

  nsink <- nrow(Total_Eigen_Data)
  npc   <- ncol(emb)

  # If Total_Eigen_Data includes PT_score column separately, pull it once
  Total_Eigen_PT <- as.numeric(Total_Eigen_Data[[pseudotime_column_name]])

  fit_list2 <- vector("list", nsink * npc)
  k <- 0

  # Define model + NLL once
  dens_model <- function(pars, Xt, Y, X0, param_t, deltat, absdt, fac, order) {
    lam     <- pars[1]
    sigma_x <- pars[2]
    mu_z    <- pars[3]
    sigma_z <- pars[4]

    # base mean vector (length Xt)
    base_mean <- X0 + param_t * (Y - X0)

    # i = 0 term
    w0 <- exp(-lam * deltat)
    sd0 <- sqrt(sigma_x^2 * absdt)
    dens <- w0 * stats::dnorm(Xt, mean = base_mean, sd = sd0)

    if (order > 0) {
      # i = 1..order terms
      i <- 1:order
      lam_absdt <- lam * absdt
      pow_mat <- outer(lam_absdt, i, `^`)                 # N x order
      wi_mat  <- (w0 * pow_mat) / matrix(fac[i + 1],nrow = length(absdt), ncol = order,byrow = TRUE)
      for (j in seq_along(i)){
        ij <- i[j]
        sdj <- sqrt(sigma_x^2 * absdt + ij * sigma_z^2)
        sdj <- pmax(sdj, 1e-12)
        dens <- dens + wi_mat[, j] * stats::dnorm(Xt, mean = base_mean + ij * mu_z, sd = sdj)
      }
    }
    dens
  }

  nll <- function(pars, Xt, Y, X0, param_t, deltat, absdt, fac, order) {
    d <- dens_model(pars, Xt, Y, X0, param_t, deltat, absdt, fac, order)
    # avoid log(0)
    eps <- 1e-300
    -sum(log(pmax(d, eps)), na.rm = TRUE)
  }

  # Main loops
  for (sink_iter in seq_len(nsink)) {

    # param_t depends only on sink pseudotime scaling
    den <- Total_Eigen_PT[sink_iter] - pt0
    if (!is.finite(den) || den <= 0) {
      param_t <- rep(1, length(deltat))
    } else {
      param_t <- pmin(1, deltat / den)
    }

    # optional warm start per sink
    x0 <- c(lam = 0.5, sigma.x = 1, mu.z = 1, sigma.z = 1)

    for (PC_dim_iter in seq_len(npc)) {

      Xt <- emb[-start_point, PC_dim_iter]
      X0 <- emb[start_point, PC_dim_iter]

      # sink value for this PC: you used PC_dim_iter+1; keep your convention
      Y  <- as.numeric(Total_Eigen_Data[sink_iter, PC_dim_iter + 1])

      fit <- nloptr::nloptr(
        x0     = x0,
        eval_f = function(pars) nll(pars, Xt, Y, X0, param_t, deltat, absdt, fac, order),
        lb     = c(0, 0, -Inf, 0),
        ub     = c(Inf, Inf, Inf, Inf),
        opts   = list(algorithm = "NLOPT_LN_SBPLX",
                      maxeval   = maxeval,
                      xtol_rel  = xtol_rel)
      )

      k <- k + 1
      fit_list2[[k]] <- fit
      x0 <- fit$solution  # warm start for next PC
    }

    message("Finished sink: ", sink_iter)
  }

  list(Probability_matrix = FALSE, model_fits = fit_list2)
}






#' Get a table and genes that associated with Jumps
#'
#' Important to note, we identify genes from the principle components not the transformed PCs
#' @param seurat_obj the seurat object to identify the jumps
#' @param Y  A set of data points for sink points from GetEigenData
#' @param state_grouping_column_name the column name to be the states you determine your model by, often it is cell type metadata in seurat.
#' @param pseudotime_column_name the column name to be the pseudotime, often it is Combined_Ordering.
#' @param batch_correction_column_name the column that contains sample/mouse/patient differences.
#' @param model_fits The model fit equations from Reduced_Jump_Prob_fast
#'
#' @return it returns a list of two entries, 1) is the list of PCs and state clusters that have jumps. 2) is the list of candidate genes to look into that significant contributors to the PC.
#' @export
#' @import Seurat
#' @import tibble
#' @import dplyr
#' @import stats
#' @importFrom tidyr pivot_longer
#'
#' @examples PC_genes_on_jump_processes<-GetJumpProgramFeatures(seurat_obj,Y,model_fits_list, state_grouping_column_name="custom_cell_classes_fine",pseudotime_column_name="Combined_Ordering",batch_correction_column_name="sample")
#'
GetJumpProgramFeatures<-function(seurat_obj,Y,model_fits_list,
                                 state_grouping_column_name="custom_cell_classes_fine",
                                 pseudotime_column_name="Combined_Ordering",
                                 batch_correction_column_name="sample"){

  emb <- Seurat::Embeddings(seurat_obj, "transformed_pca")  # cells x PCs
  meta <- Seurat::FetchData(seurat_obj, vars = c(state_grouping_column_name,
                                                 pseudotime_column_name,
                                                 batch_correction_column_name))

  # Align and order by pseudotime once
  meta <- meta[match(rownames(emb), rownames(meta)), , drop = FALSE]

  down_Z <-emb%>%
    as.data.frame%>%
    tibble::rownames_to_column(var="row_names")%>%
    dplyr::inner_join(meta%>%
                        tibble::rownames_to_column(var="row_names"),
                      by="row_names")%>%
    tibble::column_to_rownames(var="row_names")%>%
    dplyr::ungroup()%>%
    dplyr::arrange(.data[[pseudotime_column_name]])


  list_of_lambdas<-lapply(model_fits_list[["model_fits"]],function(x){
    return(x[["solution"]][1])})%>%as.matrix
  idx_for_jump_models<-(which(list_of_lambdas<max(down_Z[[pseudotime_column_name]])&(list_of_lambdas>0)))

  long_loading_pc_df_norm<-down_Z%>%
    tibble::rownames_to_column(var="Cell")%>%
    as.data.frame%>%
    tidyr::pivot_longer(cols=!c(Cell,state_grouping_column_name,batch_correction_column_name,state_grouping_column_name,pseudotime_column_name),
                        names_to = "PC_name",
                        values_to="score")

  sink_levels<-as.character(Y[[state_grouping_column_name]])

  jump_table<-matrix(FALSE,nrow=length(unique(long_loading_pc_df_norm$PC_name)),ncol=length(unique(long_loading_pc_df_norm[[state_grouping_column_name]])))
  rownames(jump_table)<-unique(long_loading_pc_df_norm$PC_name)
  colnames(jump_table)<-sink_levels
  for(iter in 1:length(idx_for_jump_models)){

    jump_table[idx_for_jump_models[iter]%%nrow(jump_table),sink_levels[ceiling(idx_for_jump_models[iter]/nrow(jump_table))]]<-TRUE
  }

  jump_df<-jump_table%>%
    as.data.frame%>%
    tibble::rownames_to_column(var="PC_names")%>%
    tidyr::pivot_longer(cols=!c(PC_names),names_to = "Sink_Type",values_to="Jump")%>%
    dplyr::group_by(PC_names)%>%
    dplyr::filter(Jump==TRUE)
  PCs_to_look_into<-jump_df%>%dplyr::pull(PC_names)%>%unique
  PCs_to_look_into<-gsub("T","",PCs_to_look_into)
  print(PCs_to_look_into)
  feature_loadings_mat<-seurat_obj@reductions[["pca"]]@feature.loadings # *sce_obj@reductions[["pca"]]@stdev

  scale_loadings_mat<-abs(apply(feature_loadings_mat,2,function(x){scale(x)}))*(seurat_obj@reductions[["pca"]]@stdev/100)
  rownames(scale_loadings_mat)<-rownames(feature_loadings_mat)

  candidate_genes_to_explore<-apply(scale_loadings_mat[,PCs_to_look_into],2,function(x){
    names(x)[which(x>stats::quantile(x)[4]+1.5*stats::IQR(x))]
    sort(x[which(x>stats::quantile(x)[4]+1.5*stats::IQR(x))])
  })

  return(list("jumps_df" =jump_df,
              "candidate_genes"= candidate_genes_to_explore))

}



#' Perform Gene Ontology to find the Relevant Ontology Program Extraction (JumpROPE)
#'
#' Important to note, we identify genes from the principle components not the transformed PCs
#' @param seurat_obj the seurat object to identify the jumps
#' @param jump_gene_list the output list from GetJumpProgramFeatures that contains genes from PCs
#'
#' @return A list of 2 gene ontology sets for each PC where at least a single jump is located.
#' @export
#' @import topGO
#' @import GO.db
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#'
#' @examples GO_terms_in_jump_PCs<-JumpROPE(seurat_obj,jump_gene_list)
#'
JumpROPE<-function(seurat_obj,jump_gene_list){


  Full_gene_list<-jump_gene_list[["candidate_genes"]]

  feature_loadings_mat <- seurat_obj@reductions[["pca"]]@feature.loadings
  all.genes<-rownames(feature_loadings_mat)
  postive_feature_loading_list<-list()
  negative_feature_loading_list<-list()
  for(iter in 1:length(Full_gene_list)){
    full_feature_loading_list<-feature_loadings_mat[names(Full_gene_list[[names(Full_gene_list)[iter]]]),names(Full_gene_list)[iter]]
    postive_feature_loading_list_temp<-full_feature_loading_list[which(full_feature_loading_list>0)]
    negative_feature_loading_list_temp<-full_feature_loading_list[which(full_feature_loading_list<0)]
    postive_feature_loading_list[[iter]]<-postive_feature_loading_list_temp
    negative_feature_loading_list[[iter]]<-negative_feature_loading_list_temp
  }

  pos_genetables_full<-lapply(postive_feature_loading_list,function(x){
    expressed.genes <-names(x)
    geneList <- ifelse(all.genes%in%expressed.genes, 1, 0)
    names(geneList) <- all.genes
    GOdata <- new("topGOdata",
                  ontology = "BP", # use biological process ontology
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  annot = topGO::annFUN.org, mapping = "org.Mm.eg", ID = "symbol")
    # Test for enrichment using Fisher's Exact Test
    resultFisher <- topGO::runTest(GOdata, algorithm = "elim", statistic = "fisher")
    #runTest(GOdata, algorithm = "elim", statistic = "fisher")
    topGO::GenTable(GOdata, Fisher = resultFisher, topNodes = 500, numChar = 1000)
    #GenTable(GOdata,resultFisher,rm.one=TRUE,topNodes=200,pvalCutOff=0.05)

  }	)

  names(pos_genetables_full)<-paste0(names(Full_gene_list),"_positive")

  neg_genetables_full<-lapply(negative_feature_loading_list,function(x){
    expressed.genes <-names(x)
    geneList <- ifelse(all.genes%in%expressed.genes, 1, 0)
    names(geneList) <- all.genes
    GOdata <- new("topGOdata",
                  ontology = "BP", # use biological process ontology
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  annot = topGO::annFUN.org, mapping = "org.Mm.eg", ID = "symbol")
    # Test for enrichment using Fisher's Exact Test
    resultFisher <- topGO::runTest(GOdata, algorithm = "elim", statistic = "fisher")
    #runTest(GOdata, algorithm = "elim", statistic = "fisher")
    topGO::GenTable(GOdata, Fisher = resultFisher, topNodes = 500, numChar = 1000)
    #GenTable(GOdata,resultFisher,rm.one=TRUE,topNodes=200,pvalCutOff=0.05)

  }	)
  names(neg_genetables_full)<-paste0(names(Full_gene_list),"_negative")

  return(list("postive_genes"= pos_genetables_full,
              "negative_genes" = neg_genetables_full))

}



