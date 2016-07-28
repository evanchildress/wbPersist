library(dplyr)
library(data.table)

nRivers<-4
nStages<-2
nYears<-20
nTimes<-365*nYears
nSamples<-nYears*4

sampleTimes<-round(seq(1,nTimes,length.out=nSamples))

#river specific daily environmental variables
hourlyTemp<-rnorm(nTimes*24,0,1) %>%
  cbind(.,.,.,.) %>%
  data.table() %>%
  .[,day:=rep(1:nTimes,each=24)] %>%
  setnames(c("wb","jimmy","mitchell","obear","day"))

dailyTemp<-hourlyTemp[,.(max(wb),
                        max(jimmy),
                        max(mitchell),
                        max(obear)),
                      by=day] %>%
  .[,.(V1,V2,V3,V4)] %>%
  as.matrix()

dailyFlow<-rnorm(nTimes,0,1) %>%
  cbind(.,.,.,.) %>%
  data.table() %>%
  setnames(c("wb","jimmy","mitchell","obear")) %>%
  .[,day:=1:nrow(.)] %>%
  .[,sample:=sum(sampleTimes<=day),by=day]

seasonalFlow<-dailyFlow[,.(mean(wb),
                           mean(jimmy),
                           mean(mitchell),
                           mean(obear)),
                        by=sample] %>%
  .[,.(V1,V2,V3,V4)] %>%
  as.matrix()

hourlyTemp<-as.matrix(hourlyTemp[,.(wb,jimmy,mitchell,obear)])
dailyFlow<-as.matrix(dailyFlow[,.(wb,jimmy,mitchell,obear)])

#return egg number from fork length based on Letcher et al. 2007
spawn<-function(forkLength){
  totalEggs<-0.00187*forkLength^2.19
  femaleEggs<-round(totalEggs/2)
  return(femaleEggs)
}

#betas from survival model
phiBeta<-array( c(6,-1,-1,1,-1,
                  6,-1,-1,1,-1,
                  6,-1,-1,1,-1,
                  6,-1,-1,1,-1,
                  6,-1,-1,1,-1,
                  6,-1,-1,1,-1,
                  6,-1,-1,1,-1,
                  6,-1,-1,1,-1,
                  6,-1,-1,1,-1),
                dim=c(5,4,2))


#create array of predicted daily logit survival based only on flow and temp
dailyLogitPhi<-array(dim=c(nTimes,nRivers,nStages)) #time,river,stage
for(r in 1:nRivers){
  for(s in 1:nStages){
    dailyLogitPhi[,r,s]<-phiBeta[1,r,s] +
                         phiBeta[2,r,s]*flow[,r] +
                         phiBeta[3,r,s]*temp[,r] +
                         phiBeta[4,r,s]*temp[,r]*flow[,r]
  }
}
sourceCpp("src/survive.cpp")
envLogitPhi<-array(dim=c(nSamples,nRivers,nStages))
for(r in 1:nRivers){
  for(g in 1:nStages){
    for(s in 1:nSamples){
      envLogitPhi[s,r,g]<-prod(plogis(dailyLogitPhi[sampleTimes[s]:(sampleTimes[s]-1)]))
    }
  }
}

#get surival by combining daily predictions with size effect
survive<-function(sample,river,stage,forkLength){
    logitSeasonalPhi<-envLogitPhi[sample,river,stage]+phiBeta[5,river,stage]*forkLength
  return(plogis(logitSeasonalPhi))
}




growthPerformance<-array(dim=c(nSamples,nRivers))
for(r in 1:nRivers){
  for(s in 1:nSamples){
    
  }
}

#parameters from growth model
list(tOpt=c(13,13,13,13),
     ctMax=c(20,20,20,20),
     sigma=c(4,4,4,4),
     beta1=c(0.015,0.015,0.015,0.015),
     beta2=c(-0.0006,-0.0006,-0.0006,-0.0006),
     beta3=c(0.01,0.01,0.01,0.01),
     beta4=c(-0.05,-0.05,-0.05,-0.05),
     eps=c(0.00049,0.00049,0.00049,0.00049))

grow<-function(forkLength)

#
core<-createCoreData() %>%
      addTagProperties() %>%
      filter(species=="bkt") %>%
      addSampleProperties() %>%
      data.table()
S<-L<-R<-G<-array(dim=c(100000,nSamples))
initN<-rpois(1,1000)
stageAge<-sample(1:nrow(core[river=="west brook"&season==3&!is.na(observedLength)]),
                 initN)

S[1:initN,1]<-1
L[1:initN,1]<-core[river=="west brook"&season==3&!is.na(observedLength)] %>%
                .[stageAge,observedLength]
G[1:initN,1]<-core[river=="west brook"&season==3&!is.na(observedLength)] %>%
  .[stageAge,ifelse(cohort==year,1,2)]
R[1:initN,1]<-1

alive<-which(S[,1]==1)

phi<-survive(R[alive,1],G[alive,1],L[alive,1],phiBeta[5,,],envLogitPhi[1,,])
