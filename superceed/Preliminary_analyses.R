#MODELS BIOMASS ~ CONNECTIVITY
library(lme4)
library(mgcv)
library(nlme)
library(lmerTest)

dataBIC<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_on_going/Connectivity data/Fullmatrixdatabiomassformodels.csv",h=T)
summary(dataBIC)
#scaling data
dataBIC$LogIndegree<-log1p(dataBIC$Indegree)
dataBIC$LogLocalRet<-log1p(dataBIC$LocalRet)
dataBIC$LogOutdegree<-log1p(dataBIC$Outdegree)
dataBIC$LogBtw<-log1p(dataBIC$Btw)
dataBIC$logbiomassarea1 <- log(dataBIC$biomassarea1)
dataBIC$logbiomassarea2 <- log(dataBIC$biomassarea2)
dataBIC$logbiomassarea <- log(dataBIC$biomassarea)
dataBIC$loggravity <- log1p(dataBIC$grav_total)

#Connectivity - Differences among groups
fullmodel<-lmer(logbiomassarea2 ~ LogIndegree * LogBtw * loggravity * 
                 Larval_behaviour *Class + (1|locality), data=dataBIC)
summary(fullmodel)


##Fished areas
fishedmodel<-lmer(logbiomassarea2 ~ LogIndegree * LogBtw * loggravity * 
                  Larval_behaviour + (1|locality), data=dataBIC[dataBIC$Class == "Fished",])
summary(fishedmodel)
visreg::visreg(fishedmodel, "LogIndegree")
