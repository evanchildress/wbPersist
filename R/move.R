move<-function(moveProb){
  to<-matrix(ncol=ncol(moveProb),nrow=nrow(moveProb))
  for(i in 1:nrow(to)){
    to[i,]<-rmultinom(1,1,moveProb[i,])
  }
  out<-apply(to,1,function(x){which(x==1)})
  return(out)
}