
##Maina
##model selection using glmmtmb
##procdure includes testing for VIF
##creating all possible combinations 
##fitting linear models of gamma family

library(usdm)
library(glmmTMB)

all.data<-read.csv(here("_data","fulldatabaseupdated2403.csv"),h=T)

colnames(all.data)




PredictVar<-dataBIC[,c("grav_total","LocalRet","Indegree","Outdegree","Btw")]