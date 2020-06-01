# Script to check general trends in bleaching drivers
# Authors: Steph and Emily
# Date: December 2017

# Clean environment
rm(list=ls())

# Load packages
#source(file="R code/brt.functions.R")
library(here)
library(PerformanceAnalytics)
library(MASS)
library(tree)
library(gbm)
library(rpart)
library(dismo)
require(dplyr)
require(here)
require(forcats)
require(brms)

# Load data
### BRT

# load data
rm(all.data)
all.data<-read.csv(here("_data","fulldatabaseupdated2403.csv"),h=T)

colnames(all.data)

rm(PredictVar)
#PredictVar<-all.data[,c("region","locality","sites","Richness","temp","grav_total","Age_of_protection","Indegree","btwdegree","Inflow","Outdegree","InflowLR","SelfR", "Class")]
PredictVar<-all.data[,c("Richness","temp","grav_total","Age_of_protection","Indegree","btwdegree","Inflow","Outdegree","InflowLR","SelfR", "Class")]

##standrdize x variables
rm(data.std)
data.std<-data.frame(apply(X = PredictVar[,4:13], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}))

data.std$region<-as.factor(all.data$region)
data.std$locality<-as.factor(all.data$locality)
data.std$sites<-as.factor(all.data$sites)
data.std$ModelMode<-as.factor(all.data$ModelMode)
data.std$ModelMode<-as.factor(all.data$ModelMode)
data.std$Class<-as.factor(all.data$Class)
data.std$log_biomassarea1<-log(all.data$biomassarea1)
data.std$Larval_behaviour<-as.factor(all.data$Larval_behaviour)
data.std$FE<-as.factor(all.data$FE)
data.std$ModelMode<-as.factor(all.data$ModelMode)

# Fished as the reference
data.std$Class <- relevel(data.std$Class, ref="Fished")
data.std$Larval_behaviour <- relevel(data.std$Larval_behaviour, ref="passive")
data.std$FE <- relevel(data.std$FE, ref="crypto")

head(data.std)
names(data.std)
summary(data.std)

# filter data with active
data.std <- data.std %>% filter(Larval_behaviour == "active" & ModelMode == "transi15")
dim(data.std)
summary(data.std)

# split
require(caret)
require(klaR)
set.seed(0123)
split=0.90
rm(trainIndex); trainIndex <- createDataPartition(data.std$Richness, p=split, list=FALSE)
rm(data_train); data_train <- data.std[ trainIndex,]
rm(data_test); data_test <- data.std[-trainIndex,]

  # Check for normality of biomass data 
  hist(data.std$log_biomassarea1)

  
### Simple BRT model
  # The response variable we are interested in 
  rm(Resp)
  Resp <- which(colnames(data.std) == "Richness"); Resp

  # Predictors
  rm(Pred)
  Pred=c("temp","grav_total",
         "Age_of_protection","Indegree","btwdegree","Inflow","Outdegree",
         "InflowLR","SelfR", "Class")
  

  # Ncol = the vector containing the column number of each variable use to predict S
  rm(Ncol); Ncol<- which(colnames(data.std) %in% Pred); Ncol 
  # chck
  length(Pred) == length(Ncol)
  
  # 1st Step: Built the full brt model and check for goodness of fit
  rm(BIOM_brt)
  BIOM_brt <- gbm.step(gbm.y=Resp,gbm.x=Ncol,data=data.std,tree.complexity=10,learning.rate=0.0005,bag.fraction=0.7,n.trees=50,family="gaussian",n.folds=10,max.trees=10000)
  
  # output of the model:
  ls(BIOM_brt)
  
  # Contributions of each variable
  BIOM_brt$contributions
  
  # summary brt
  summary(BIOM_brt,las=2)
  abline(v=5,col="red")
  
  # marginal distributions
  windows()
  plot.gbm(BIOM_brt,11)
  
  # CV AUC
  BIOM_brt$cv.statistics$discrimination.mean
  
  # plot observed vs residuals
  plot(data.std$log_biomassarea1,BIOM_brt$fit,
       xlim=c(3,8),ylim=c(3,8),
       xlab="Observed Bleaching Int.",ylab="Predicted Bleaching Int.",pch=16)
  abline(a=0,b=1)
    # linear model between predicted vs fitted
  lm_BlEACH <- lm(BLEACH_brt$fit~BVAR$bleach_intensity)
  abline(lm_BlEACH$coefficients[1],lm_BlEACH$coefficients[2],lty=2,col="red")
  # plot residuals vs observed data
  plot(BVAR$bleach_intensity,BLEACH_brt$residuals,
       xlim=c(0,100),ylim=c(-10,10),
       xlab="Observed Bleaching Int",ylab="Residuals",pch=16)
  abline(h=0)
  
  hist(BLEACH_brt$residuals)
  # plot residuals vs fitted data
  plot(BLEACH_brt$fit,BLEACH_brt$residuals,
       xlim=c(0,100),
       xlab="Fitted Biomass",ylab="Residuals",pch=16)
  abline(h=0)
  
  # 2. Look for interactions
  rm(BIOM_int)
  BIOM_int <- gbm.interactions(BIOM_brt)  
  BIOM_int$rank.list
  # plot the 3D plot between the most interacting variables
  gbm.perspec(BIOM_int, 2, 1,z.range=c(0,60))
  gbm.perspec(BIOM_int, 9, 10,z.range=c(0,60)) # sd spell and avg spell
  gbm.perspec(BIOM_int, 9, 5,z.range=c(0,60)) # bimod. coef and avg spell
  gbm.perspec(BIOM_int, 9, 3,z.range=c(0,60)) # bimod. coef and avg spell
  gbm.perspec(BIOM_int, 8, 3,z.range=c(0,60)) # bimod. coef and avg spell
  #etc.
  
  # 3. Simplify the model
  BIOM_brt_simpl <- gbm.simplify(BIOM_brt,n.drops=10,n.folds=10) # 7 last variables to remove
  
    # New model with only significant explanatory variables
  rm(BLEACH_brt_simpl)
  BLEACH_brt_SIMPL <- gbm.step(gbm.y=Resp,gbm.x=BLEACH_brt_simpl$pred.list[[7]],data=BVAR,tree.complexity=10,learning.rate=0.005,bag.fraction=0.7,n.trees=50,family="gaussian",n.folds=10,max.trees=10000)
  
  # Contributions
  BLEACH_brt_SIMPL$contributions
  
  # Interactions
  rm(BLEACH_int_SIMPL)
  BLEACH_int_SIMPL <- gbm.interactions(BLEACH_brt_SIMPL)  
  BLEACH_int_SIMPL$rank.list
  gbm.perspec(BLEACH_brt_SIMPL, 5, 4,z.range=c(0,60)) # sd spell and avg spell
  gbm.perspec(BLEACH_brt_SIMPL, 6, 3,z.range=c(0,60)) # bimod. coef and avg spell

    
  
