
# load packages
library(here)
require(dplyr)
require(tidyr)
require(forcats)
require(ggplot2)
require(mgcv)
require(visreg)
require(ggpubr)
require(GGally)
require(corrplot)
require(lme4)
source("stdCoef.lmer.R")
require(kableExtra) # nicetable - html
require(scales)
require(itsadug)
## load data
# SEM analysis - Global connectivity 
# author: Steph D'agata
# date: April 2020
# outputs: SEM coefficients

rm(list=ls())

# load data
rm(all.data)
all.data<-read.csv(here("_data","FullDataMay2020Coordinates.csv"),h=T)
# clean first column
all.data$X.1 <- NULL
all.data$X <- NULL

colnames(all.data)

# list of predictors
rm(PredictVar)
PredictVar<-all.data[,c("region","locality","sites","Richness","temp","grav_total",
                        "Age_of_protection","Class","Indegree","ModelMode",         
                        "btwdegree","Inflow","Outdegree","InflowLR",          
                        "SelfR","Larval_behaviour","FE","InflowBR",          
                        "IndegreeBR","CorridorIndegreeBR","grav_neiBR","IndegreeMPABR",     
                        "InflowMPABR","IndegreeNeiBR","InflowNeiBR")]
PredictVar$log_grav_total <- log(PredictVar$grav_total+1)
PredictVar$log_grav_neiBR <- log(PredictVar$grav_neiBR+1)
PredictVar$log_biomassarea1<-log(PredictVar$biomassarea1+1)

##standrdize x variables
rm(data.std)
# standardize
data.std<-data.frame(apply(X = PredictVar[,c(4,5,7,9,11:15,18:20,22:27)], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}))
# add log biomass
data.std$log_biomassarea1<-log(all.data$biomassarea1+1)

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
data.std$Larval_behaviour <- droplevels(data.std$Larval_behaviour)
data.std$FE <- relevel(data.std$FE, ref="crypto")
data.std$FE <- droplevels(data.std$FE)

head(data.std)
names(data.std)
summary(data.std)

# filter data with active
data.std <- data.std %>% filter(Larval_behaviour == "active" & ModelMode == "transi15")
dim(data.std)
summary(data.std)

### Explore SEM to disentangle the effect of coral cover and species richness on total biomass
# load packages
library(ape) #Version 3.3
library(caper) # Vresion 0.5.2
library(nlme) # Version 3.1.122
library(lavaan) # Version 0.5.19
# Load piecewiseSEM from CRAN
library(piecewiseSEM) # Version 1.0.0
library(lmerTest)
library(devtools)
##### create k-folds for cross validation
require(caret)
require(klaR)

# define an 80%/20% train/test split of the dataset
set.seed(123)
split=0.90
rm(trainIndex); trainIndex <- createDataPartition(DATA$n_spp, p=split, list=FALSE)
rm(data_train); data_train <- DATA[ trainIndex,]
rm(data_test); data_test <- DATA[-trainIndex,]

summary(data_train)
str(data_train)
##### Create full piecewise model explaining total biomass 
#set.seed(1234)
rm(Biomass_TOTAL_randomList_train)
Biomass_TOTAL_randomList_train = psem(
  
  # Predicting Number of species
#  lmer(Richness ~ temp + Class +log_grav_total + log_grav_neiBR + IndegreeMPABR + 
#       IndegreeBR + CorridorIndegreeBR + SelfR + (1 |sites/locality/region),data=data.std),
  lmer(Richness ~ temp + Class +log_grav_total + log_grav_neiBR + IndegreeMPABR + (1 |sites/locality/region),data=data.std),
         
  # Predicting total biomass
 #    log_grav_neiBR + IndegreeMPABR + 
  #     IndegreeBR + CorridorIndegreeBR + SelfR+ (1 |sites/locality/region),data=data.std)
    lmer(log_biomassarea1 ~ Richness+temp + Class +log_grav_total +  (1 |sites/locality/region),data=data.std)
           
     
)


plot(Biomass_TOTAL_randomList_train)
# sem fits
summary(Biomass_TOTAL_randomList_train)

new.summary <- summary(Biomass_TOTAL_randomList_train, .progressBar = F)
class(new.summary)
# Old function
sem.missing.paths(Biomass_TOTAL_randomList_train,data=data_train, .progressBar = F)


# prediction with test data
Biomass_TOTAL_randomList_test <- sem.predict(Biomass_TOTAL_randomList_train, 
                                             newdata=data_test, sefit = TRUE)
head(Biomass_TOTAL_randomList_test)

# goodness of fit for each
  # species richness
plot(data_test$n_spp,Biomass_TOTAL_randomList_test$n_spp.fit,xlim=c(10,80),ylim=c(10,80),pch=16)
abline(a=0,b=1)
cor(data_test$n_spp,Biomass_TOTAL_randomList_test$n_spp.fit)
RMSE(data_test$n_spp,Biomass_TOTAL_randomList_test$n_spp.fit)
# biomass
plot(data_test$logTotBiom,Biomass_TOTAL_randomList_test$logTotBiom.fit,xlim=c(1.5,3.5),ylim=c(1.5,3.5),pch=16)
abline(a=0,b=1)
cor(data_test$logTotBiom,Biomass_TOTAL_randomList_test$logTotBiom.fit)
RMSE(data_test$logTotBiom,Biomass_TOTAL_randomList_test$logTotBiom.fit)
SSE = sum((Biomass_TOTAL_randomList_test$Biomass.hard_coral.fit - data_test$logTotBiom)^2)
SST = sum((mean(data_train$logTotBiom) - data_test$logTotBiom) ^ 2)
R2 = 1- SSE/SST;R2
# coral cover
plot(data_test$hard_coral,Biomass_TOTAL_randomList_test$hard_coral.fit,xlim=c(0,100),ylim=c(0,100),pch=16)
abline(a=0,b=1)
cor(data_test$hard_coral,Biomass_TOTAL_randomList_test$hard_coral.fit)
RMSE(data_test$hard_coral,Biomass_TOTAL_randomList_test$Coral.hard_coral.fit)
SSE = sum((Biomass_TOTAL_randomList_test$hard_coral.fit - data_test$hard_coral)^2)
SST = sum((mean(data_train$hard_coral) - data_test$hard_coral) ^ 2)
R2 = 1- SSE/SST;R2

### models separatly
rm(CORAL)
CORAL <- lm(hard_coral ~ logGrav_NC + #reef30km
     mean_depth  + npp+ management_rules+
     we +reef3km+reef30km+#macroalgae +
     reef300km,data=data_train)
# Predicting Number of species
rm(SPECIES)
SPECIES <- lm(n_spp ~ hard_coral+logGrav_NC+management_rules + #NB_SP = 
      mean_depth + reef_zone + we +latitude + 
      reef30km + reef300km,data = data_train) # ,random = ~ 1|country, data = data_train,na.action = na.exclude)

# Predicting total biomass
rm(BIOM)
BIOM <- lm(logTotBiom ~ n_spp + hard_coral+logGrav_NC+management_rules + #NB_SP = 
      mean_depth + reef_zone + reef30km+
      reef300km, data = data_train,na.action = na.exclude) #random = ~ 1|country

# color
library(viridis)
col <- scale_color_viridis(discrete=TRUE)
col[1]
library(scales)
show_col(viridis_pal()(20))
test <- viridis_pal()(20)
### test creating a predictor matrix:
  # 1. high connectivity, low gravity, fore reef, mean latitude, mean depth
rm(HIGH_OA)
# low gravity 
lg <- quantile(DATA$logGrav_NC,0.05)
logGrav_NC <- rep(0.175,200)

# high connectivity 3km
hr3 <- quantile(DATA$reef3km,0.95)
reef3km<- as.numeric(rep(14684774,200))

# high connectivity 30km
hr30 <- quantile(DATA$reef30km,0.95)
reef30km <- rep(266781793,200)

# high connectivity 300km
hr300 <- quantile(DATA$reef300km,0.95)
reef300km <- rep(5912125415,200)

# depth 10 meters
mean_depth <- rep(10,200)

# wave energy
we <- rep(mean(DATA$we),200)

# management rules - Open Acess
management_rules <- as.factor(rep("Open Access",200))
levels(management_rules) <- c(levels(management_rules),"Gear Restriction","No-Take")

# reef zone  - fore reef
reef_zone <- as.factor(rep("fore reef",200))
levels(reef_zone) <- c(levels(reef_zone),"back reef")
# latitude - sample uniformely
latitude <- rep(mean(DATA$latitude),200)

# hard coral - sample uniformely
hard_coral <- round(sort(runif(200,min=0,max=86)),2)

# npp - sample uniformely
npp <- rep(mean(DATA$npp),200)

# n_spp
n_spp <- rep(NA,200)

# hard_coral
hard_coral  <- seq(0,84,0.421)


HIGH.OA <- data.frame(logGrav_NC,reef3km,reef30km,reef300km,reef_zone,management_rules,
                      latitude,mean_depth,npp,n_spp,hard_coral,we)


rm(CORAL.pred)
CORAL.pred <- predict(CORAL,newdata=HIGH.OA)
rm(SPECIES.pred)
SPECIES.pred <- predict(SPECIES,newdata=HIGH.OA)
HIGH.OA$n_spp <- SPECIES.pred

rm(BIOMASS.pred)
BIOMASS.pred <- predict(BIOM,newdata=HIGH.OA)

str(HIGH.OA)
summary(HIGH.OA)

#plot(HIGH.OA$hard_coral,BIOMASS.pred,type="l",ylim=c(1.2,3.8))
pdf(here("_prelim_figures","Pred_biomass.pdf"),width=12,height=12)
par(mfrow=c(2,2),mar=c(4.5,4.5,1.1,1))
plot(HIGH.OA$hard_coral,10^BIOMASS.pred,type="l",ylim=c(0,4000),xlab="Hard Coral cover",
     ylab="Predicted biomass (kg ha-1)",col="#B8DE29FF",lwd=5,main="Sc1: High connectivity - Low gravity",cex.axis=1.4,cex.lab=1.6)
text(x=75,y=540,"500kg ha-1",col="red",cex=1.2,font=2)
# management rules - Gear management
management_rules <- as.factor(rep("Gear Restriction",200))
levels(management_rules) <- c(levels(management_rules),"Open Access","No-Take")
HIGH.OA$management_rules <- management_rules
rm(SPECIES.pred)
SPECIES.pred <- predict(SPECIES,newdata=HIGH.OA)
rm(BIOMASS.pred)
BIOMASS.pred <- predict(BIOM,newdata=HIGH.OA)
lines(HIGH.OA$hard_coral,10^BIOMASS.pred,type="l",col="#39558CFF",lwd=5)

# management rules - No-Take
management_rules <- as.factor(rep("No-Take",200))
levels(management_rules) <- c(levels(management_rules),"Open Access","Gear Restriction")
HIGH.OA$management_rules <- management_rules
rm(SPECIES.pred)
SPECIES.pred <- predict(SPECIES,newdata=HIGH.OA)
rm(BIOMASS.pred)
BIOMASS.pred <- predict(BIOM,newdata=HIGH.OA)
lines(HIGH.OA$hard_coral,10^BIOMASS.pred,type="l",col="#440154FF",lwd=5)
abline(h=500,col="red",lty=2)
       
### test creating a predictor matrix:
# 1. low connectivity, high gravity, fore reef, mean latitude, mean depth
rm(LOW_OA)
# low gravity 
lg <- quantile(DATA$logGrav_NC,0.95)
logGrav_NC <- rep(1.957273,200)

# high connectivity 3km
hr3 <- quantile(DATA$reef3km,0.05)
reef3km<- as.numeric(rep(1354355,200))

# high connectivity 30km
hr30 <- quantile(DATA$reef30km,0.05)
reef30km <- rep(35183124,200)

# high connectivity 300km
hr300 <- quantile(DATA$reef300km,0.05)
reef300km <- rep(373761064,200)

# depth 10 meters
mean_depth <- rep(10,200)

# wave energy
we <- rep(mean(DATA$we),200)

# management rules - Open Acess
management_rules <- as.factor(rep("Open Access",200))
levels(management_rules) <- c(levels(management_rules),"Gear Restriction","No-Take")

# reef zone  - fore reef
reef_zone <- as.factor(rep("fore reef",200))
levels(reef_zone) <- c(levels(reef_zone),"back reef")
# latitude - sample uniformely
latitude <- rep(mean(DATA$latitude),200)

# npp - sample uniformely
npp <- rep(mean(DATA$npp),200)

# n_spp
n_spp <- rep(NA,200)

# hard_coral
hard_coral  <- seq(0,84,0.421)


LOW.OA <- data.frame(logGrav_NC,reef3km,reef30km,reef300km,reef_zone,management_rules,
                      latitude,mean_depth,npp,n_spp,hard_coral,we)


rm(CORAL.pred)
CORAL.pred <- predict(CORAL,newdata=LOW.OA)
rm(SPECIES.pred)
SPECIES.pred <- predict(SPECIES,newdata=LOW.OA)
LOW.OA$n_spp <- SPECIES.pred

rm(BIOMASS.pred)
BIOMASS.pred <- predict(BIOM,newdata=LOW.OA)

str(LOW.OA)
summary(LOW.OA)

#plot(HIGH.OA$hard_coral,BIOMASS.pred,type="l",ylim=c(1.2,3.8))
plot(HIGH.OA$hard_coral,10^BIOMASS.pred,type="l",ylim=c(0,520),xlab="Hard Coral cover",
     ylab="Predicted biomass (kg ha-1)",col="#B8DE29FF",lwd=5,main="Sc2: Low connectivity - High gravity",cex.axis=1.4,cex.lab=1.6)
LOW.OA[which(BIOMASS.pred > 2.69 & BIOMASS.pred <2.7),"hard_coral"]
# management rules - Gear management
management_rules <- as.factor(rep("Gear Restriction",200))
levels(management_rules) <- c(levels(management_rules),"Open Access","No-Take")
LOW.OA$management_rules <- management_rules
rm(SPECIES.pred)
SPECIES.pred <- predict(SPECIES,newdata=LOW.OA)
rm(BIOMASS.pred)
BIOMASS.pred <- predict(BIOM,newdata=LOW.OA)
lines(LOW.OA$hard_coral,10^BIOMASS.pred,type="l",col="#39558CFF",lwd=5)

# management rules - No-Take
management_rules <- as.factor(rep("No-Take",200))
levels(management_rules) <- c(levels(management_rules),"Open Access","Gear Restriction")
LOW.OA$management_rules <- management_rules
rm(SPECIES.pred)
SPECIES.pred <- predict(SPECIES,newdata=LOW.OA)
rm(BIOMASS.pred)
BIOMASS.pred <- predict(BIOM,newdata=LOW.OA)
lines(LOW.OA$hard_coral,10^BIOMASS.pred,type="l",col="#440154FF",lwd=5)
abline(h=500,col="red",lty=2)
text(x=75,y=510,"500kg ha-1",col="red",cex=1.2,font=2)
### test creating a predictor matrix:
# 1. med connectivity, med gravity, fore reef, mean latitude, mean depth
rm(MED_OA)
# low gravity 
lg <- quantile(DATA$logGrav_NC,0.5)
logGrav_NC <- rep(0.88375,200)

# high connectivity 3km
hr3 <- quantile(DATA$reef3km,0.5)
reef3km<- as.numeric(rep(7162614,200))

# high connectivity 30km
hr30 <- quantile(DATA$reef30km,0.5)
reef30km <- rep(117346693,200)

# high connectivity 300km
hr300 <- quantile(DATA$reef300km,0.5)
reef300km <- rep(2489036294,200)

# depth 10 meters
mean_depth <- rep(10,200)

# wave energy
we <- rep(mean(DATA$we),200)

# management rules - Open Acess
management_rules <- as.factor(rep("Open Access",200))
levels(management_rules) <- c(levels(management_rules),"Gear Restriction","No-Take")

# reef zone  - fore reef
reef_zone <- as.factor(rep("fore reef",200))
levels(reef_zone) <- c(levels(reef_zone),"back reef")
# latitude - sample uniformely
latitude <- rep(mean(DATA$latitude),200)

# npp - sample uniformely
npp <- rep(mean(DATA$npp),200)

# n_spp
n_spp <- rep(NA,200)

# hard_coral
hard_coral  <- seq(0,84,0.421)


MED.OA <- data.frame(logGrav_NC,reef3km,reef30km,reef300km,reef_zone,management_rules,
                     latitude,mean_depth,npp,n_spp,hard_coral,we)


rm(CORAL.pred)
CORAL.pred <- predict(CORAL,newdata=MED.OA)
rm(SPECIES.pred)
SPECIES.pred <- predict(SPECIES,newdata=MED.OA)
MED.OA$n_spp <- SPECIES.pred

rm(BIOMASS.pred)
BIOMASS.pred <- predict(BIOM,newdata=MED.OA)

str(MED.OA)
summary(MED.OA)

#plot(HIGH.OA$hard_coral,BIOMASS.pred,type="l",ylim=c(1.2,3.8))
plot(HIGH.OA$hard_coral,10^BIOMASS.pred,type="l",ylim=c(0,1000),xlab="Hard Coral cover",
     ylab="Predicted biomass (kg ha-1)",col="#B8DE29FF",lwd=5,main="Sc3: Medium connectivity - Medium gravity",cex.axis=1.4,cex.lab=1.6)
MED.OA[which(BIOMASS.pred > 2.697 & BIOMASS.pred <2.699),"hard_coral"]
text(x=75,y=520,"500kg ha-1",col="red",cex=1.2,font=2)
# management rules - Gear management
management_rules <- as.factor(rep("Gear Restriction",200))
levels(management_rules) <- c(levels(management_rules),"Open Access","No-Take")
MED.OA$management_rules <- management_rules
rm(SPECIES.pred)
SPECIES.pred <- predict(SPECIES,newdata=MED.OA)
rm(BIOMASS.pred)
BIOMASS.pred <- predict(BIOM,newdata=MED.OA)
lines(LOW.OA$hard_coral,10^BIOMASS.pred,type="l",col="#39558CFF",lwd=5)
MED.OA[which(BIOMASS.pred > 2.697 & BIOMASS.pred <2.699),"hard_coral"]

# management rules - No-Take
management_rules <- as.factor(rep("No-Take",200))
levels(management_rules) <- c(levels(management_rules),"Open Access","Gear Restriction")
MED.OA$management_rules <- management_rules
rm(SPECIES.pred)
SPECIES.pred <- predict(SPECIES,newdata=MED.OA)
rm(BIOMASS.pred)
BIOMASS.pred <- predict(BIOM,newdata=MED.OA)
lines(LOW.OA$hard_coral,10^BIOMASS.pred,type="l",col="#440154FF",lwd=5)
abline(h=500,col="red",lty=2)
segments(42.5,0,42.5, y1 = 500,col="#B8DE29FF",lty=3,lwd=2)
text(40,20,"42.5",col="#B8DE29FF",font=2,cex=1.5,srt=90)
segments(9.3,0, 9.3, y1 = 500,col="#39558CFF",lty=3,lwd=2)
text(7,20,"9.3",col="#39558CFF",font=2,cex=1.5, srt=90)
segments(3,0, 3, y1 = 500,col="#440154FF",lty=3,lwd=2)
text(1,20,"3",col="#440154FF",font=2,cex=1.5, srt=90)

MED.OA[which(BIOMASS.pred > 2.697 & BIOMASS.pred <2.699),"hard_coral"]

# common legend

dev.off()
