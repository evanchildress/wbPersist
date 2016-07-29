growthPerformance<-array(dim=c(nSamples-1,nRivers))
for(r in 1:nRivers){
  for(s in 1:(nSamples-1)){
    growthPerformance[s,r]<-predictPerformance(hourlyTemp[sampleHours[s]:sampleHours[s+1]],
                       growthPars$tOpt[r],
                       growthPars$ctMax[r],
                       growthPars$sigma[r]) %>%
      sum()
  }
}

grow<-function(forkLength){
  return(forkLength+20)
}