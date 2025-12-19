
#' Parallel R implementation to get Cell-Cell Transition Probability Matrix
#'
#' @param sigma_points 
#' @param Y 
#' @param model_fits 
#' @param sce 
#' @param n_cores 
#' @param cell_indices_to_use 
#'
#' @return
#' @export
#'
#' @examples
GetCellToCellTPM_parallel<-function(sigma_points,Y,model_fits,sce,n_cores=NULL,cell_indices_to_use=NULL){
  
  eigValues = (sce@reductions[["pca"]]@stdev)^2  ## EigenValues
  varExplained = eigValues / sum(eigValues)#total_variance
  Prob_Xt_X0_mat=matrix(data=NA_real_,nrow=dim(sigma_points)[1] ,ncol=dim(sigma_points)[1])
  if(is.null(cell_indices_to_use)){
    cell_indices_to_use = 1:nrow(sigma_points)
  }
  if(!is.null(n_cores)){
    doParallel::registerDoParallel(cores = n_cores)
  }else{
    doParallel::registerDoParallel()
  }
  
  for(starting_point in cell_indices_to_use){
    Prob_X0_given_Xt_per_PC =NULL
    
    Prob_X0_given_Xt_per_PC<-foreach::foreach(PC_iter=1:(ncol(sigma_points)-3),.combine='cbind')%dopar%{
      out<-foreach::foreach(sink_iter=1:nrow(Y),.combine='cbind')%do%{
        
        JDD_multivariate_model2_with_param(pars=model_fits[["model_fits"]][[(sink_iter-1)*(ncol(sigma_points)-3)+PC_iter]][["solution"]],
                                           data=as.numeric(sigma_points[-starting_point,PC_iter]),
                                           sinks=as.numeric(Y[sink_iter,PC_iter+1]),
                                           deltat =abs(as.numeric(sigma_points[-starting_point,(ncol(sigma_points)-1)])-as.numeric(sigma_points[starting_point,(ncol(sigma_points)-1)])),
                                           y_dts=Y$PT_score[sink_iter],
                                           X0=as.numeric(sigma_points[starting_point,PC_iter]),
                                           t0=as.numeric(sigma_points[starting_point,(ncol(sigma_points)-1)]))
        
      }
      rowSums(out)
      #Prob_X0_given_Xt_per_PC = cbind(Prob_X0_given_Xt_per_PC,rowSums(Prob_X0_given_Xt_Yk_for_PC))
    }
    PC_normalized_Prob<-apply(Prob_X0_given_Xt_per_PC,2,function(x){x/sum(x)}) # Normalize the densities for each PC
    Prob_Xt_given_X0 <-apply(PC_normalized_Prob,1,function(x){prod(x^(varExplained))})
    
    Prob_Xt_X0_mat[-starting_point,starting_point]<-Prob_Xt_given_X0
    print(starting_point)
  }
  doParallel::stopImplicitCluster()
  
  colnames(Prob_Xt_X0_mat)<-rownames(sigma_points)
  rownames(Prob_Xt_X0_mat)<-rownames(sigma_points)
  
  return(t(Prob_Xt_X0_mat))
  
  
}

