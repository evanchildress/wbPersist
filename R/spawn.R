#' Returns number of eggs produced from fork lengths based on Letcher et al. 2007
#' @param forkLength Vector of fork lengths to estimate fecundity for
#' @import data.table
#' @export

spawn<-function(forkLength,river){
  streamLength<-c(5421,715,390)
  streamLength<-streamLength/sum(streamLength)

  totalEggs<-0.00187*forkLength^2.19
  femaleEggs<-totalEggs/2

  eggs<-data.table::data.table(river=river,femaleEggs=femaleEggs) %>%
                list(data.table::data.table(river=1:4,femaleEggs=rep(0,4))) %>%
                data.table::rbindlist() %>%
                .[river %in% 1:3,river:=1]
  eggs<-eggs[,.(femaleEggs=sum(femaleEggs)),by=river]
  eggsOut<-c(eggs[river==1,femaleEggs]*streamLength,eggs[river==4,femaleEggs])

  return(eggsOut)
}
