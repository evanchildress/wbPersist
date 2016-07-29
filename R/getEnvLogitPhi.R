
getEnvLogitPhi<-function(nSamples,nDays,sampleDays,phiBeta,nRivers=4,nStages=2){
    #create array of predicted daily logit survival based only on flow and temp
    dailyLogitPhi<-array(dim=c(nDays,nRivers,nStages)) #time,river,stage
    for(r in 1:nRivers){
      for(s in 1:nStages){
        dailyLogitPhi[,r,s]<-phiBeta[1,r,s] +
          phiBeta[2,r,s]*dailyFlow[,r] +
          phiBeta[3,r,s]*dailyTemp[,r] +
          phiBeta[4,r,s]*dailyTemp[,r]*dailyFlow[,r]
      }
    }

    envLogitPhi<-array(dim=c(nSamples,nRivers,nStages))
    for(r in 1:nRivers){
      for(g in 1:nStages){
        for(s in 1:nSamples){
          envLogitPhi[s,r,g]<-prod(plogis(dailyLogitPhi[sampleDays[s]:(sampleDays[s]-1)]))
        }
      }
    }
    return(envLogitPhi)
}
