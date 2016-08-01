#' Returns number of eggs produced from fork lengths based on Letcher et al. 2007
#' @param forkLength Vector of fork lengths to estimate fecundity for
#' @export

spawn<-function(forkLength,river){
  streamLength<-c(47,17,14)
  streamLength<-streamLength/sum(streamLength)

  totalEggs<-0.00187*forkLength^2.19
  femaleEggs<-totalEggs/2

  eggs<-data.table(river=river,femaleEggs=femaleEggs) %>%
                bind_rows(data.table(river=1:4,femaleEggs=rep(0,4))) %>%
                .[,river2:=river] %>%
                .[river2 %in% 1:3,river2:=1]
  eggs<-eggs[,.(femaleEggs=sum(femaleEggs)),by=river2]
  eggsOut<-c(eggs[river2==1,femaleEggs]*streamLength,eggs[river2==4,femaleEggs])

  return(eggsOut)
}
