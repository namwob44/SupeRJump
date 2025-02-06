#' Get Visitation Probability and Expected timesteps for TPM
#'
#' @param Nmat 
#' @useDynLib SuperJump
#' @return
#' @export
#'
#' @examples
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
