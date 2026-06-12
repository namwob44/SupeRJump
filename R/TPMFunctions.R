#' This is a summation over the Jump drift diffusion model.
#'
#' @param pars_mat fitted parameters matrix from the Reduced_Jump_Prob
#' @param Xt is all cells we are transitioning towards.
#' @param Yp is the sink/eigen cells we use for particular PC
#' @param deltat this is the absolute value for the pseudo time from the starting cell to all other cells
#' @param y_dts this is the pseudo time of the sink cells Yp
#' @param X0 initial starting cell we are finding transition probabilities from.
#' @param t0 initial starting pseudo time.
#' @param order this is a truncated value to consider the poisson process in the infinite sum.
#' @param eps_sd just a stabilizing control
#'
#' @return probability density of each state Xt
#' @export
#' @import stats
#'
#' @examples \dontrun{
#' for (idx in blocks) {
#' out <- out + JDD_sum_over_sinks_vec(
#'   pars_mat = pars_mat[idx, , drop = FALSE],
#'   Xt       = Xt,
#'   Yp       = Yp[idx],
#'   deltat   = deltat,
#'   y_dts    = y_dts[idx],
#'   X0       = X0,
#'   t0       = t0,
#'   order    = order,
#'   eps_sd   = eps_sd
#'   )
#'   }}
JDD_sum_over_sinks_vec <- function(pars_mat,Xt,Yp,deltat,y_dts,X0,t0,order = 5,eps_sd = 1e-12) {
  # pars_mat: n_sinks x 4
  # columns: lambda, sigma.x, mu.z, sigma.z
  #
  # Xt:      length n_targets
  # Yp:      length n_sinks
  # deltat:  length n_targets
  # y_dts:   length n_sinks
  #
  # returns: length n_targets, summed over sinks

  nsink <- length(Yp)
  ntar  <- length(Xt)

  lam     <- pars_mat[, 1]
  sigma_x <- pars_mat[, 2]
  mu_z    <- pars_mat[, 3]
  sigma_z <- pars_mat[, 4]

  absdt <- abs(deltat)

  # Matrix dimensions:
  # rows = target cells
  # cols = sinks

  dt_mat    <- matrix(deltat, nrow = ntar, ncol = nsink)
  absdt_mat <- matrix(absdt,  nrow = ntar, ncol = nsink)

  lam_mat      <- matrix(lam,     nrow = ntar, ncol = nsink, byrow = TRUE)
  sigx_mat     <- matrix(sigma_x, nrow = ntar, ncol = nsink, byrow = TRUE)
  muz_mat      <- matrix(mu_z,    nrow = ntar, ncol = nsink, byrow = TRUE)
  sigz_mat     <- matrix(sigma_z, nrow = ntar, ncol = nsink, byrow = TRUE)
  Y_mat        <- matrix(Yp,      nrow = ntar, ncol = nsink, byrow = TRUE)
  Xt_mat       <- matrix(Xt,      nrow = ntar, ncol = nsink)
  ydt_mat      <- matrix(y_dts,   nrow = ntar, ncol = nsink, byrow = TRUE)

  den <- ydt_mat - t0

  # Match your original behavior as closely as possible:
  # param_t = min(1, deltat / (y_dts - t0)).
  # But protect against zero / negative / non-finite denominators.
  param_t <- dt_mat / den
  param_t[!is.finite(param_t)] <- 1
  param_t <- pmin(1, param_t)

  base_mean <- X0 + param_t * (Y_mat - X0)

  exp_term <- exp(-lam_mat * dt_mat)

  # i = 0 term
  sd0 <- sqrt(sigx_mat^2 * absdt_mat)
  sd0 <- pmax(sd0, eps_sd)

  dens <- exp_term * stats::dnorm(
    Xt_mat,
    mean = base_mean,
    sd   = sd0
  )

  if (order > 0) {
    for (i in seq_len(order)) {
      weight_i <- exp_term * ((lam_mat * absdt_mat)^i) / factorial(i)

      sdi <- sqrt(sigx_mat^2 * absdt_mat + i * sigz_mat^2)
      sdi <- pmax(sdi, eps_sd)

      dens <- dens + weight_i * stats::dnorm(
        Xt_mat,
        mean = base_mean + i * muz_mat,
        sd   = sdi
      )
    }
  }
  rowSums(dens)
}

#' Extract the model fit parameters into a matrix given sinks and PCs
#'
#' @param model_fits the large list of model filts.
#' @param nsink total number of sink cells we have
#' @param npc The total number of PCs we have (commonly 30-50)
#'
#' @return A matrix of all parameters from the fitted values.
#' @export
#'
#' @examples fit_array <- extract_fit_array(model_fits, nsink = nsink, npc = npc)
#'
extract_fit_array <- function(model_fits, nsink, npc) {
  fit_array <- array(
    NA_real_,
    dim = c(nsink, npc, 4),
    dimnames = list(
      sink = NULL,
      PC = NULL,
      par = c("lambda", "sigma.x", "mu.z", "sigma.z")
    )
  )

  fits <- model_fits[["model_fits"]]

  for (sink_iter in seq_len(nsink)) {
    for (PC_iter in seq_len(npc)) {
      idx <- (sink_iter - 1L) * npc + PC_iter
      fit_array[sink_iter, PC_iter, ] <- fits[[idx]][["solution"]]
    }
  }

  fit_array
}


#' Blocked Structured to solve Jump drift diffusion, helpful for massive amount of cells in Seurat object
#'
#' @param pars_mat fitted parameters matrix from the Reduced_Jump_Prob
#' @param Xt is all cells we are transitioning towards.
#' @param Yp is the sink/eigen cells we use for particular PC
#' @param deltat this is the absolute value for the pseudo time from the starting cell to all other cells
#' @param y_dts this is the pseudo time of the sink cells Yp
#' @param X0 initial starting cell we are finding transition probabilities from.
#' @param t0 initial starting pseudo time.
#' @param order this is a truncated value to consider the poisson process in the infinite sum.
#' @param eps_sd just a stabilizing control
#
#' @return it returns a probability density of the transitions from a starting cell.
#' @export
#'
#' @examples \dontrun{
#' for (PC_iter in seq_len(npc)) {
#' pars_mat <- fit_array[, PC_iter, , drop = FALSE]
#' pars_mat <- matrix(pars_mat, nrow = nsink, ncol = 4)
#' prob_per_pc[, PC_iter] <- JDD_sum_over_sinks_blocked(pars_mat   = pars_mat,
#' Xt         = Xmat[target_idx, PC_iter],
#' Yp         = Ymat[, PC_iter],
#'  deltat     = deltat,
#'  y_dts      = ypt,
#'  X0         = Xmat[starting_point, PC_iter],
#'   t0         = t0,
#'   order      = order,
#'   block_size = 100)
#'   }}
#'
JDD_sum_over_sinks_blocked <- function(pars_mat,Xt,Yp,deltat,y_dts,X0,t0,order = 5,block_size = 100,eps_sd = 1e-12) {
  nsink <- length(Yp)
  out <- numeric(length(Xt))

  blocks <- split(seq_len(nsink), ceiling(seq_len(nsink) / block_size))

  for (idx in blocks) {
    out <- out + JDD_sum_over_sinks_vec(
      pars_mat = pars_mat[idx, , drop = FALSE],
      Xt       = Xt,
      Yp       = Yp[idx],
      deltat   = deltat,
      y_dts    = y_dts[idx],
      X0       = X0,
      t0       = t0,
      order    = order,
      eps_sd   = eps_sd
    )
  }

  out
}




#' Wrapper Function for the cell-to-cell Transition Probability Matrix (TPM) with batch correction and parallelization
#'
#' @param seurat_obj the main seurat object to analyze.
#' @param Y The sink cells identified from GetSinkData
#' @param model_fits The model fit equations from Reduced_Jump_Prob_fast
#' @param state_grouping_column_name the column name to be the states you determine your model by, often it is cell type metadata in seurat.
#' @param pseudotime_column_name the column name to be the pseudotime, often it is Combined_Ordering.
#' @param batch_correction_column_name the column that contains sample/mouse/patient differences.
#' @param n_cores this is how many parallel threads to run at once. If using 8,
#' @param cell_indices_to_use set to NULL, index if only interested in a particular cell.
#' @param order this is a truncated value to consider the poisson process in the infinite sum.
#' @param eps stability value for logprob transform.
#' @param normalize_rows set this to true to turn into a markov process, default is true.
#' @param verbose if want status updates on samples.

#' @return A seurat object with an updated graphs called TPM.
#' @export
#' @import tibble
#' @import dplyr
#' @import tibble
#' @import parallel
#' @import Matrix
#' @useDynLib SupeRJump, .registration = TRUE

#' @examples  \dontrun{seurat_obj<-GetCellToCellTPM_fast(seurat_obj,Y,Prob_reduced,state_grouping_column_name="custom_cell_classes_fine",pseudotime_column_name = "Combined_Ordering",batch_correction_column_name="sample",n_cores = 8,cell_indices_to_use = NULL,order = 5,eps = 1e-300,normalize_rows = TRUE,verbose = TRUE)}
#'
GetCellToCellTPM_fast <- function(seurat_obj,Y,model_fits, state_grouping_column_name,pseudotime_column_name,batch_correction_column_name,n_cores = NULL,cell_indices_to_use = NULL,order = 5,eps = 1e-300,normalize_rows = FALSE,verbose = TRUE) {

  all_cells <- colnames(seurat_obj)
  n_global  <- length(all_cells)

  meta <- Seurat::FetchData(seurat_obj,vars = c(state_grouping_column_name,pseudotime_column_name,batch_correction_column_name))

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
    Sparse_TPM <- Individual_Batch_TPM(subset(seurat_obj, cells = batch_cells),Y,model_fits,state_grouping_column_name,pseudotime_column_name,batch_correction_column_name,n_cores = n_cores,cell_indices_to_use = NULL,order = order,eps = eps,normalize_rows = normalize_rows,verbose = verbose)
    # ---- map local → global indices ----
    local_cells <- rownames(Sparse_TPM)
    local_to_global <- match(local_cells, all_cells)

    # ---- extract sparse entries ----
    sp <- summary(Sparse_TPM)   # gives i, j, x (1-based)

    # map to global indices
    I_all[[counter]] <- local_to_global[sp$i]
    J_all[[counter]] <- local_to_global[sp$j]
    V_all[[counter]] <- sp$x

    counter <- counter + 1
  }

  # ---- combine all batches ----
  I_all <- unlist(I_all, use.names = FALSE)
  J_all <- unlist(J_all, use.names = FALSE)
  V_all <- unlist(V_all, use.names = FALSE)

  TPM_sparse <- Matrix::sparseMatrix(i = I_all,j = J_all,x = V_all,dims = c(n_global, n_global))

  rownames(TPM_sparse) <- all_cells
  colnames(TPM_sparse) <- all_cells

  seurat_obj@graphs[["TPM"]] <- TPM_sparse

  return(seurat_obj)
}


#' The main function to get Transition Probability Matrices.
#'
#' This function gets the TPM for individual batches.
#' @param seurat_obj the main seurat object to analyze.
#' @param Y The sink cells identified from GetSinkData
#' @param model_fits The model fit equations from Reduced_Jump_Prob_fast
#' @param state_grouping_column_name the column name to be the states you determine your model by, often it is cell type metadata in seurat.
#' @param pseudotime_column_name the column name to be the pseudotime, often it is Combined_Ordering.
#' @param batch_correction_column_name the column that contains sample/mouse/patient differences.
#' @param n_cores this is how many parallel threads to run at once. If using 8,
#' @param cell_indices_to_use set to NULL, index if only interested in a particular cell.
#' @param order this is a truncated value to consider the poisson process in the infinite sum.
#' @param eps stability value for logprob transform.
#' @param normalize_rows set this to true to turn into a markov process, default is true.
#' @param verbose if want status updates on samples.
#'
#' @return a transition probability matrix before being added to a seurat graphs section
#' @export
#' @import tibble
#' @import dplyr
#' @import tibble
#' @import parallel
#' @import Matrix
#' @useDynLib SupeRJump, .registration = TRUE
#'
#' @examples   \dontrun{Sparse_TPM <- Individual_Batch_TPM(subset(seurat_obj, cells = batch_cells),Y,model_fits,state_grouping_column_name,pseudotime_column_name,batch_correction_column_name,n_cores = n_cores,cell_indices_to_use = NULL,order = order,eps = eps,normalize_rows = normalize_rows,verbose = verbose)}

Individual_Batch_TPM <- function(seurat_obj,Y,model_fits,state_grouping_column_name,pseudotime_column_name,batch_correction_column_name,n_cores = NULL,cell_indices_to_use = NULL,order = 5,eps = 1e-300,normalize_rows = FALSE,verbose = TRUE) {

  sigma_points <-seurat_obj[["transformed_pca"]]@cell.embeddings%>%
    as.data.frame%>%
    tibble::rownames_to_column(var="row_names")%>%
    dplyr::inner_join(Seurat::FetchData(seurat_obj,vars=c(state_grouping_column_name,pseudotime_column_name,batch_correction_column_name))%>%
                        tibble::rownames_to_column(var="row_names"),
                      by="row_names")%>%
    tibble::column_to_rownames(var="row_names")%>%
    dplyr::ungroup()%>%
    dplyr::arrange(.data[[pseudotime_column_name]])



  ncell <- nrow(sigma_points)
  npc   <- ncol(sigma_points) - 3
  nsink <- nrow(Y)

  if (is.null(cell_indices_to_use)) {
    cell_indices_to_use <- seq_len(ncell)
  }

  # Eigenvalue weights
  eigValues <- (seurat_obj@reductions[["transformed_pca"]]@stdev)^2
  varExplained <- eigValues / sum(eigValues)

  # Use only PCs present in sigma_points
  varExplained <- varExplained[seq_len(npc)]
  varExplained <- varExplained / sum(varExplained)

  # Numeric matrices/vectors once
  Xmat <- as.matrix(sigma_points[, seq_len(npc), drop = FALSE])
  pt   <- as.numeric(sigma_points[, ncol(sigma_points) - 1])

  Ymat <- as.matrix(Y[, 2:(npc + 1), drop = FALSE])
  ypt  <- as.numeric(Y[[pseudotime_column_name]])

  # Extract fitted parameters once:
  # dim = nsink x npc x 4
  fit_array <- extract_fit_array(model_fits, nsink = nsink, npc = npc)

  # Worker for one starting cell
  compute_one_start <- function(starting_point) {
    target_idx <- setdiff(seq_len(ncell), starting_point)

    t0 <- pt[starting_point]
    deltat <- abs(pt[target_idx] - t0)

    # Matrix: target cells x PCs
    prob_per_pc <- matrix(
      NA_real_,
      nrow = length(target_idx),
      ncol = npc
    )

    for (PC_iter in seq_len(npc)) {
      pars_mat <- fit_array[, PC_iter, , drop = FALSE]
      pars_mat <- matrix(pars_mat, nrow = nsink, ncol = 4)

      prob_per_pc[, PC_iter] <- JDD_sum_over_sinks_blocked(pars_mat   = pars_mat,
                                                           Xt         = Xmat[target_idx, PC_iter],
                                                           Yp         = Ymat[, PC_iter],
                                                           deltat     = deltat,
                                                           y_dts      = ypt,
                                                           X0         = Xmat[starting_point, PC_iter],
                                                           t0         = t0,
                                                           order      = order,
                                                           block_size = 100)
    }

    # Normalize densities for each PC across target cells
    col_sums <- colSums(prob_per_pc, na.rm = TRUE)
    col_sums[!is.finite(col_sums) | col_sums <= 0] <- NA_real_

    pc_norm <- sweep(prob_per_pc, 2, col_sums, "/")
    pc_norm[!is.finite(pc_norm)] <- eps
    pc_norm <- pmax(pc_norm, eps)

    # Weighted geometric mean:
    # old: apply(pc_norm, 1, function(x) prod(x ^ varExplained))
    # new: exp(log(pc_norm) %*% varExplained)
    probs <- as.numeric(exp(log(pc_norm) %*% varExplained))

    out <- rep(NA_real_, ncell)
    out[target_idx] <- probs
    out[starting_point] <- NA_real_

    if (verbose) {
      message("Finished start cell: ", starting_point)
    }

    out
  }

  # Parallel over starting cells
  if (!is.null(n_cores) && n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      library(SuperJump)
      NULL
    })
    parallel::clusterExport(
      cl,
      varlist = c(
        "ncell", "npc", "nsink",
        "Xmat", "pt", "Ymat", "ypt",
        "fit_array", "varExplained",
        "order", "eps",
        "JDD_sum_over_sinks_vec",
        "JDD_sum_over_sinks_blocked"
      ),
      envir = environment()
    )

    res_list <- parallel::parLapply(
      cl,
      X = cell_indices_to_use,
      fun = compute_one_start
    )
  } else {
    res_list <- lapply(cell_indices_to_use, compute_one_start)
  }
  message("past the TPM fitting")


  message("starting to fill and sparsify TPM")

  n <- ncell
  groups <- as.character(sigma_points[[state_grouping_column_name]])
  threshold <- 1 / n
  I <- vector("list", length(cell_indices_to_use))
  J <- vector("list", length(cell_indices_to_use))
  V <- vector("list", length(cell_indices_to_use))

  for (idx in seq_along(cell_indices_to_use)){
    i <- cell_indices_to_use[idx]
    row_vals <-res_list[[idx]]
    gi <- groups[i]
    same_group <- (groups==gi)
    keep <- same_group | (row_vals >= threshold)
    if (!any(keep)) next
    cols_to_use <- which(keep)
    I[[idx]] <- rep(i, length(cols_to_use))
    J[[idx]] <- cols_to_use
    V[[idx]] <- row_vals[cols_to_use]

  }
  # Flatten lists
  I_all <- unlist(I, use.names = FALSE)
  J_all <- unlist(J, use.names = FALSE)
  V_all <- unlist(V, use.names = FALSE)
  TPM_sparse <- Matrix::sparseMatrix(i = I_all,j = J_all,x = V_all,dims = c(n, n)  )
  rownames(TPM_sparse) <-  rownames(sigma_points)
  colnames(TPM_sparse) <-  rownames(sigma_points)


  # Optional row normalization if you want each start-cell transition vector
  # to sum to 1 after PC integration.
  if (normalize_rows) {
    message("Converting to Markov")
    rs <- Matrix::rowSums(TPM_sparse, na.rm = TRUE)
    valid <- is.finite(rs) & rs > 0
    TPM_sparse[valid, ] <- TPM_sparse[valid, , drop = FALSE] / rs[valid]
  }

  return(TPM_sparse)
}
