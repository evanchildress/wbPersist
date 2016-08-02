library(wbPersist)
library(dplyr)
library(data.table)
library(getWBData)
library(reshape2)
reconnect()


nRivers<-4
nStages<-2
nYears<-17
nDays<-365*nYears
nSamples<-nYears*4
season<-rep(c(2,3,4,1),nYears)

########################################################################
#load in (or make up data)
########################################################################
sampleDays<-round(seq(1,nDays,length.out=nSamples))
sampleHours<-round(seq(1,nDays*24,length.out=nSamples))
#river specific daily environmental variables
hourlyTemp<-runif(nDays*24,0,10) %>%
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

hourlyTemp<-tbl(conDplyr,"data_hourly_temperature") %>%
            collect(n=Inf) %>%
            data.table() %>%
            .[datetime>as.POSIXct("1999-01-01 00:00:00")] %>%
            .[,.(datetime,temperature,river)] %>%
            .[!duplicated(.[,.(datetime,river)])]

dailyTemp<-hourlyTemp[,.(temperature=max(temperature)),by=.(date=as.Date(datetime),river)] %>%
            melt(id.vars=c("river","date")) %>%
            acast(date~river) %>%
            .[1:nDays,]

hourlyTemp<-melt(hourlyTemp,id.vars=c("river","datetime")) %>%
            acast(datetime~river)

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

#hourlyTemp<-as.matrix(hourlyTemp[,.(wb,jimmy,mitchell,obear)])
dailyFlow<-as.matrix(dailyFlow[,.(wb,jimmy,mitchell,obear)])

#################################################################################
#load in (or make up) parameters for surival, movement, growth, and reproduction
#################################################################################

#betas from survival model
phiBeta<-array( c(5,-1,-1,1,-1/300,
                  5,-1,-1,1,-1/300,
                  5,-1,-1,1,-1/300,
                  5,-1,-1,1,-1/300,
                  5,-1,-1,1,-1/300,
                  5,-1,-1,1,-1/300,
                  5,-1,-1,1,-1/300,
                  5,-1,-1,1,-1/300),
                dim=c(5,4,2))


#parameters from growth model
growthPars<-list(tOpt=c(13,13,13,13),
     ctMax=c(20,20,20,20),
     sigma=c(4,4,4,4),
     betas=matrix(c(0.018,0.018,0.018,0.018,
                    -0.00005,-0.00005,-0.00005,-0.00005,
                    0.004,0.004,0.004,0.004,
                  -0.002,-0.002,-0.002,-0.002),
                  nrow=4,ncol=4,byrow=T),
     eps=c(0.00049,0.00049,0.00049,0.00049))


#movement probabilities from Letcher et al. 2014
moveProb<-movementProbs


#
core<-createCoreData() %>%
      addTagProperties() %>%
      filter(species=="bkt") %>%
      addSampleProperties() %>%
      data.table()
#
# core<-data.table(river="west brook",season=3,
#                  observedLength=runif(1000,60,200)) %>%
#       .[,cohort:=ifelse(observedLength<95,1,2)] %>%
#       .[,year:=1]

##################################################################
#create structures that contain predictions that do not depend on individual traits
##################################################################

envLogitPhi<-getEnvLogitPhi(nSamples,nDays,sampleDays,phiBeta)
growthPerformance<-getGrowthPerformance(sampleHours,
                                        hourlyTemp,
                                        growthPars)

################################################################
#run simulations
################################################################


#set up data structures
A<-L<-R<-G<-array(dim=c(100000,nSamples))
initN<-rpois(1,1000)
maxIndexAlive<-initN

stageAge<-sample(1:nrow(core[river=="west brook"&season==3&!is.na(observedLength)]),
                 initN,replace=T)
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

for(i in 2:5){
phi<-survive(R[alive,i-1],G[alive,i-1],L[alive,i-1],phiBeta[5,,],envLogitPhi[1,,])

#who survives?
A[alive,i]<-rbinom(length(alive),1,phi)

#grow, change stages, and move, if you survived
alive<-which(A[,i]==1)
if(length(alive)==0) break
L[alive,i]<-grow(growthPars$betas,R[alive,i-1],L[alive,i-1],
                 seasonalFlow[i-1,],rep(0,length(alive)),
                 growthPerformance[i-1,])
G[alive,i]<-stageTransition(G[alive,i-1],season[i-1])

R[alive,i]<-move(getMoveProb(R[alive,i-1],moveProb[,,season[i-1]]))

#reproduce if it's the right season


if(season[i]==3){
  eggs<-rep(as.integer(NA),4)
  eggPhi<-eggSurvive(dailyFlow,dailyTemp,1,1)
  #spawning needs to happen in fall! otherwise alive is wrong, so eggs that survive should do that during fall by looking ahead then survive with probability 1 until next fall to save space used by eggs
  eggs<-spawn(L[alive,i][G[alive,i]==2],
                   R[alive,i][G[alive,i]==2]) %>%
        round()


  nYoy<-rbinom(nRivers,eggs,eggPhi)
  if(sum(nYoy)>0){
    A[(maxIndexAlive+1):(maxIndexAlive+sum(nYoy)),i+4]<-1
    R[(maxIndexAlive+1):(maxIndexAlive+sum(nYoy)),i+4]<-rep(1:4,nYoy)
    G[(maxIndexAlive+1):(maxIndexAlive+sum(nYoy)),i+4]<-1
    L[(maxIndexAlive+1):(maxIndexAlive+sum(nYoy)),i+4]<-recruitLength(R[(maxIndexAlive+1):(maxIndexAlive+sum(nYoy)),i],
                                                                   seasonalFlow,
                                                                   dailyTemp,
                                                                   1)
    maxIndexAlive<-maxIndexAlive+sum(nYoy)
    alive<-which(A[,i]==1)
    }
  }
}
#need to add the right number of YOY to the bottom of A and loop

N<-colSums(!is.na(L))
