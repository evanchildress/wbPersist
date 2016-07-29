library(dplyr)
library(data.table)
library(getWBData)
library(Rcpp)

nRivers<-4
nStages<-2
nYears<-20
nDays<-365*nYears
nSamples<-nYears*4
season<-rep(c(3,4,1,2),nYears)


sampleDays<-round(seq(1,nDays,length.out=nSamples))
sampleHours<-round(seq(1,nDays*24,length.out=nSamples))
#river specific daily environmental variables
hourlyTemp<-rnorm(nDays*24,0,1) %>%
  cbind(.,.,.,.) %>%
  data.table() %>%
  .[,day:=rep(1:nDays,each=24)] %>%
  setnames(c("wb","jimmy","mitchell","obear","day"))

dailyTemp<-hourlyTemp[,.(max(wb),
                        max(jimmy),
                        max(mitchell),
                        max(obear)),
                      by=day] %>%
  .[,.(V1,V2,V3,V4)] %>%
  as.matrix()

dailyFlow<-rnorm(nDays,0,1) %>%
  cbind(.,.,.,.) %>%
  data.table() %>%
  setnames(c("wb","jimmy","mitchell","obear")) %>%
  .[,day:=1:nrow(.)] %>%
  .[,sample:=sum(sampleDays<=day),by=day]

seasonalFlow<-dailyFlow[,.(mean(wb),
                           mean(jimmy),
                           mean(mitchell),
                           mean(obear)),
                        by=sample] %>%
  .[,.(V1,V2,V3,V4)] %>%
  as.matrix()

hourlyTemp<-as.matrix(hourlyTemp[,.(wb,jimmy,mitchell,obear)])
dailyFlow<-as.matrix(dailyFlow[,.(wb,jimmy,mitchell,obear)])

#betas from survival model
phiBeta<-array( c(6,-1,-1,1,-1/60,
                  6,-1,-1,1,-1/60,
                  6,-1,-1,1,-1/60,
                  6,-1,-1,1,-1/60,
                  6,-1,-1,1,-1/60,
                  6,-1,-1,1,-1/60,
                  6,-1,-1,1,-1/60,
                  6,-1,-1,1,-1/60),
                dim=c(5,4,2))


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



stageTransition<-function(stage,season){
  out<-rep(2,length(stage))
  out[which(stage==1 & season %in% c(3,4,1))]<-1
  return(out)
}





#parameters from growth model
growthPars<-list(tOpt=c(13,13,13,13),
     ctMax=c(20,20,20,20),
     sigma=c(4,4,4,4),
     betas=matrix(c(0.015,0.015,0.015,0.015,
                    -0.0006,-0.0006,-0.0006,-0.0006,
                    0.01,0.01,0.01,0.01,
                  -0.05,-0.05,-0.05,-0.05),
                  nrow=4,ncol=4,byrow=T),
     eps=c(0.00049,0.00049,0.00049,0.00049))



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


#
core<-createCoreData() %>%
      addTagProperties() %>%
      filter(species=="bkt") %>%
      addSampleProperties() %>%
      data.table()

#set up data structures
A<-L<-R<-G<-array(dim=c(100000,nSamples))
initN<-rpois(1,1000)
maxIndexAlive<-initN

stageAge<-sample(1:nrow(core[river=="west brook"&season==3&!is.na(observedLength)]),
                 initN)
#Intialize the first year
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

for(i in 2:nSamples){
phi<-survive(R[alive,i-1],G[alive,i-1],L[alive,i-1],phiBeta[5,,],envLogitPhi[1,,])

#who survives?
A[alive,i]<-rbinom(length(alive),1,phi)

#grow, change stages, and move, if you survived
alive<-which(A[,i]==1)
if(length(alive)==0) break
L[alive,i]<-grow(L[alive,i-1])
G[alive,i]<-stageTransition(G[alive,i-1],season[i-1])

R[alive,i]<-move(getMoveProb(R[alive,i-1],moveProb[,,season[i-1]]))

#reproduce if it's the right season


if(season[i]==3){
  eggs<-rep(as.integer(NA),4)
  eggPhi<-eggSurvive(dailyFlow,dailyTemp,1,1)
  for(r in 1:nRivers){
    eggs[r]<-sum(spawn(L[alive,i-4][R[alive,i-4]==r]))
  }
  nYoy<-rbinom(nRivers,eggs,eggPhi)
  A[(maxIndexAlive+1):(maxIndexAlive+sum(nYoy)),i]<-1
  R[(maxIndexAlive+1):(maxIndexAlive+sum(nYoy)),i]<-rep(1:4,nYoy)
  G[(maxIndexAlive+1):(maxIndexAlive+sum(nYoy)),i]<-1
  L[(maxIndexAlive+1):(maxIndexAlive+sum(nYoy)),i]<-recruitLength(R[(maxIndexAlive+1):(maxIndexAlive+sum(nYoy)),i],
                                                                 seasonalFlow,
                                                                 dailyTemp,
                                                                 1)
  maxIndexAlive<-maxIndexAlive+sum(nYoy)
  alive<-which(A[,i]==1)
}
}
#need to add the right number of YOY to the bottom of A and loop

N<-colSums(!is.na(L))
