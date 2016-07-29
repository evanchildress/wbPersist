stageTransition<-function(stage,season){
  out<-rep(2,length(stage))
  out[which(stage==1 & season %in% c(3,4,1))]<-1
  return(out)
}
