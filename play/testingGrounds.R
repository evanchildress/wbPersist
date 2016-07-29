library(dplyr)
library(data.table)
library(getWBData)
library(Rcpp)

nRivers<-4
nStages<-2
nYears<-20
nTimes<-365*nYears
nSamples<-nYears*4
season<-rep(c(3,4,1,2),nYears)

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
                         phiBeta[2,r,s]*dailyFlow[,r] +
                         phiBeta[3,r,s]*dailyTemp[,r] +
                         phiBeta[4,r,s]*dailyTemp[,r]*dailyFlow[,r]
  }
}

envLogitPhi<-array(dim=c(nSamples,nRivers,nStages))
for(r in 1:nRivers){
  for(g in 1:nStages){
    for(s in 1:nSamples){
      envLogitPhi[s,r,g]<-prod(plogis(dailyLogitPhi[sampleTimes[s]:(sampleTimes[s]-1)]))
    }
  }
}

growthPerformance<-array(dim=c(nSamples,nRivers))
for(r in 1:nRivers){
  for(s in 1:nSamples){
    
  }
}

stageTransition<-function(stage,season){
  out<-rep(2,length(stage))
  out[which(stage==1 & season %in% c(3,4,1))]<-1
  return(out)
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

grow<-function(forkLength){
  return(forkLength+20)
}

moveProb<-array(c(1,0,0,0,
                  0,1,0,0,
                  0,0,1,0,
                  0,0,0,1,
                  1,0,0,0,
                  0,1,0,0,
                  0,0,1,0,
                  0,0,0,1,
                  1,0,0,0,
                  0,1,0,0,
                  0,0,1,0,
                  0,0,0,1,
                  1,0,0,0,
                  0,1,0,0,
                  0,0,1,0,
                  0,0,0,1),
                dim=c(4,nRivers,nRivers))
move<-function(moveProb){
  to<-matrix(ncol=ncol(moveProb),nrow=nrow(moveProb))
  for(i in 1:nrow(to)){
    to[i,]<-rmultinom(1,1,moveProb[i,])
  }
  out<-apply(to,1,function(x){which(x==1)})
  return(out)
}

#
core<-createCoreData() %>%
      addTagProperties() %>%
      filter(species=="bkt") %>%
      addSampleProperties() %>%
      data.table()
A<-L<-R<-G<-array(dim=c(100000,nSamples))
initN<-rpois(1,1000)
maxIndexAlive<-initN

stageAge<-sample(1:nrow(core[river=="west brook"&season==3&!is.na(observedLength)]),
                 initN)
#who is alive?
A[1:initN,1]<-1

#size
L[1:initN,1]<-core[river=="west brook"&season==3&!is.na(observedLength)] %>%
                .[stageAge,observedLength]

#stage
G[1:initN,1]<-core[river=="west brook"&season==3&!is.na(observedLength)] %>%
  .[stageAge,ifelse(cohort==year,1,2)]

#location
R[1:initN,1]<-1

alive<-which(A[,1]==1)

phi<-survive(R[alive,1],G[alive,1],L[alive,1],phiBeta[5,,],envLogitPhi[1,,])

i<-2
#who survives?
S[alive,i]<-rbinom(length(S[alive,i-1]),1,phi)

#grow, change stages, and move, if you survived
alive<-which(S[,i]==1)

L[alive,i]<-grow(L[alive,i-1])
G[alive,i]<-stageTransition(G[alive,i-1],season[i-1])

R[alive,i]<-move(getMoveProb(R[alive,i-1],moveProb[,,season[i-1]]))

#reproduce if it's the right season
eggPhi<-c(0.05,0.05,0.05,0.05) #need to get the YOY model going for this

if(season[i-1]==2){
  eggs<-rep(as.integer(NA),4)
  for(r in 1:nRivers){
    eggs[r]<-sum(spawn(L[alive,i-1][R[alive,i-1]==r]))
  }
  nYoy<-rbinom(nRivers,eggs,eggPhi)
}

#need to add the right number of YOY to the bottom of A and loop
