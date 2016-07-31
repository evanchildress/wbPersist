#' Returns number of eggs produced from fork lengths based on Letcher et al. 2007
#' @param forkLength Vector of fork lengths to estimate fecundity for
#' @export

spawn<-function(forkLength){
  totalEggs<-0.00187*forkLength^2.19
  femaleEggs<-round(totalEggs/2)
  return(femaleEggs)
}
