#' Returns a vector of stages making transitions for YOY to adults after summer
#' @param stage Vector of starting stages
#' @param season Vector of starting season over which the transition happens
#' @export

stageTransition<-function(stage,season){
  out<-rep(2,length(stage))
  out[which(stage==1 & season %in% c(3,4,1))]<-1
  return(out)
}
