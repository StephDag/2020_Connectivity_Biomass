
##Maina
##model selection using glmmtmb
##procdure includes testing for VIF
##creating all possible combinations 
##fitting linear models of gamma family

library(usdm)
library(glmmTMB)
library("glmmTMB")
library("bbmle") ## for AICtab 
library("ggplot2")

source('~/Documents/Mygitprojects/localfunctions/functions_analyses_interaction_glmmTMB.R')

all.data<-read.csv(here("_data","fulldatabaseupdated2403.csv"),h=T)

colnames(all.data)

PredictVar<-all.data[,c("temp","grav_total","Age_of_protection","Indegree","btwdegree","Inflow","Outdegree","InflowLR","SelfR", "Class")]

##standrdize x variables
data.std<-data.frame(apply(X = PredictVar[,1:8], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}))

options(stringsAsFactors = FALSE)
data.std$Class<-as.character(all.data$Class)
data.std$biomassarea1<-all.data$biomassarea1
data.std$Larval_behaviour<-as.character(all.data$Larval_behaviour)
data.std$FE<-as.character(all.data$FE)
data.std$ModelMode<-as.character(all.data$ModelMode)

#1. Create list with all possible combinations between predictors 
vifPredCombinations  <-  list()
varnames<-colnames(PredictVar)#

maxCombs  <-  getMaximumNOfCombs(varnames)
for(j in 1:maxCombs) {
  vifPredCombinations  <-  append(runPredCombinations(j, varnames), vifPredCombinations)
}

##2. Create list with all possible combinations between predictors 
vifPredCombinations  <-  list()
varnames<-colnames(PredictVar)#

maxCombs  <-  getMaximumNOfCombs(varnames)
for(j in 1:maxCombs) {
  vifPredCombinations  <-  append(runPredCombinations(j, varnames), vifPredCombinations)
}

##3. filter the combinations above with VIF<1.5
vifPredCombinations_new<- c()
for(con in vifPredCombinations){
  r <- subset(PredictVar, select = con)
  conClasses   <-  unique(sapply(r, class))
  numOfClasses  <-  length(conClasses)
  twoNCols <- ncol(r)==2
  numDf <- r[sapply(r,is.numeric)]
  zeroNumDf<-ncol(numDf)==0
  numeriCols<- ncol(numDf)>1
  
  if(length(con)<2){
    vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
  }else{
    if (zeroNumDf) { 
      vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
    }
    if (numOfClasses==2 && twoNCols) { 
      vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
    }
    if(numeriCols && max(vif(numDf)["VIF"])<=1.5){##vif cutoff
      vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
    }
    next
  }
}

##configure models

modelText.biomass<-mclapply(vifPredCombinations_new, prepareModelText,data.std )


