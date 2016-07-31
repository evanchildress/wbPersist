#'Placeholder until model for length at recruitment is available
#'@param river Vector of river locations for recruits
#'@param flow some measure of pre-recruitment flow
#'@param temp some measure of pre-recruitment temp
#'@param biomass density of competitors/dominators
#'@export

recruitLength<-function(river,flow,temp,biomass){
  return(rnorm(length(river),60,10))
}
