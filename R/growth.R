#' Estimates growth performance, i.e., the temperature part of the growth model
#'
#'@param sampleHours Breaks between performance periods
#'@param hourlyTemp Matrix of hourly temperatures by river
#'@param growthPars Named list of growth parameters from the model fit
#'@param nRivers how many rivers?
#'
#'@export

getGrowthPerformance<-function(sampleHours,hourlyTemp,growthPars,nRivers=4){
  nSamples<-length(sampleHours)
  growthPerformance<-array(dim=c(nSamples-1,nRivers))

  for(r in 1:nRivers){
    for(s in 1:(nSamples-1)){
      growthPerformance[s,r]<-perform::predictPerformance(hourlyTemp[sampleHours[s]:sampleHours[s+1],r],
                         growthPars$tOpt[r],
                         growthPars$ctMax[r],
                         growthPars$sigma[r]) %>%
        sum()
    }
  }
  return(growthPerformance)
}
