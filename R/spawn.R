#return egg number from fork length based on Letcher et al. 2007
spawn<-function(forkLength){
  totalEggs<-0.00187*forkLength^2.19
  femaleEggs<-round(totalEggs/2)
  return(femaleEggs)
}