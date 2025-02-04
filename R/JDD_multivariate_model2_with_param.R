
#' Parameterized function to get cell-cell probabilty
#'
#' @param pars 
#' @param data 
#' @param sinks 
#' @param deltat 
#' @param y_dts 
#' @param X0 
#' @param t0 
#'
#' @return
#' @export
#'
#' @examples
JDD_multivariate_model2_with_param<-function(pars,data,sinks,deltat,y_dts,X0,t0){
  Xt<-data
  Y<-sinks
  lam <- pars[1]#exp(pars[1])+10^-10
  #mu.x <-pars[2]
  sigma.x <- pars[2]#exp(pars[3])+10^-10
  mu.z <- pars[3]
  sigma.z <- pars[4]#exp(pars[5])+10^-10
  order<- 5#length(Xt)
  param_t <-ifelse(((deltat)/(y_dts-as.numeric(t0)))>1,1,((deltat)/(y_dts-as.numeric(t0))))
  
  dens <- exp(-lam*deltat)*stats::dnorm(Xt,(X0+((param_t)*as.numeric(Y-X0))),sqrt(sigma.x^2*abs(deltat))) 
  for(i in 1:order){
    dens <-  dens + (exp(-lam*deltat)*((lam*abs(deltat))^i)/factorial(i))*stats::dnorm(Xt,(X0+((param_t)*as.numeric(Y-X0)))+i*mu.z,sqrt(sigma.x^2*abs(deltat)+i*sigma.z^2))
  }
  return(dens)
}
