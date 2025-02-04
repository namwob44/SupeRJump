
#' Fitting Jump Model for each Eigen Cell and PC.
#'
#' @param Z The input matrix for analysis that includes PCs and Pseudotime
#' @param Total_Eigen_Data  A set of data points for sink points from GetEigenData
#' @param start_point index for which cell is "origin" or earliest cell
#'
#' @return It returns a list of models for each EigenCell and PC
#' @export
#'
#' @examples
Reduced_Jump_Prob<-function(Z,Total_Eigen_Data,start_point){
  fit_list2=NULL
  k=0
  
  for(sink_iter in 1:(dim(Total_Eigen_Data)[1])){
    for(PC_dim_iter in 1:(ncol(Z)-3)){
      Xt_deltat= abs(as.numeric(Z[-start_point,(ncol(Z)-1)])-as.numeric(Z[start_point,(ncol(Z)-1)]))
      data_to_use_ind = as.vector(Z[-start_point,PC_dim_iter])
      sink_data = as.numeric(Total_Eigen_Data[sink_iter,PC_dim_iter+1])
      
      X0 = as.numeric(Z[start_point,PC_dim_iter])
      
      JDD_multivariate_model2<-function(pars,data,sinks,deltat){
        Xt<-data
        Y<-sinks
        lam <- pars[1]#exp(pars[1])+10^-10
        #mu.x <-pars[2]
        sigma.x <- pars[2]#exp(pars[3])+10^-10
        mu.z <- pars[3]
        sigma.z <- pars[4]#exp(pars[5])+10^-10
        order<- 5#length(Xt)
        param_t <-ifelse(((deltat)/(as.numeric(Total_Eigen_Data$PT_score[sink_iter])-as.numeric(Z$PT_score[start_point])))>1,
                         1,
                         ((deltat)/(as.numeric(Total_Eigen_Data$PT_score[sink_iter])-as.numeric(Z$PT_score[start_point]))))
        dens <- exp(-lam*deltat)*stats::dnorm(Xt,(X0+((param_t)*as.numeric(Y-X0))),sqrt(sigma.x^2*abs(deltat))) 
        for(i in 1:order){
          dens <-  dens + (exp(-lam*deltat)*((lam*abs(deltat))^i)/factorial(i))*stats::dnorm(Xt,(X0+((param_t)*as.numeric(Y-X0)))+i*mu.z,sqrt(sigma.x^2*abs(deltat)+i*sigma.z^2))
        }
        return(dens)
      }
      NLL = function(pars, data,sinks,deltat) {
        # Negative log-likelihood 
        return( (-1)*sum(log(JDD_multivariate_model2(pars,data,sinks,deltat)),na.rm = T))
        
      }
      par0_ind<-c(lam = (0.5),sigma.x=(1),mu.z=(1),sigma.z=(1.0))
      
      fit<-nloptr::nloptr(x0=par0_ind,
                          eval_f = NLL,
                          data=data_to_use_ind,
                          sinks=sink_data,
                          deltat=Xt_deltat,
                          lb=c(0,0,-Inf,0),
                          ub=c(Inf,Inf,Inf,Inf),
                          opts=list("algorithm"="NLOPT_LN_SBPLX",
                                    "maxeval" = 1000,
                                    "xtol_rel"=1e-6))
      k=k+1
      fit_list2[[k]]<-fit

    }
    print(paste0("We finished sink: ",sink_iter))
  }
  Prob_Xt_fixed_Xs_per_Yk=FALSE
  return(list("Probability_matrix"=Prob_Xt_fixed_Xs_per_Yk,
              "model_fits"=fit_list2))
}
