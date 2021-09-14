library(gridExtra)
library(ggpubr)
library(ggcorrplot)
library(RColorBrewer)
library(lme4)
library(ggplot2)
library(car)
library(glmmADMB)
library(MASS)
library(MuMIn)
library(lmerTest)
library(mgcv)
library(here)
library(hutilscpp)


#1st part 
#Connectivity Rank----
bigmama %>% rm()
bigmama<-read.csv(here("_data","global_metricsLU.csv"),h=T)
summary(bigmama)

#PASSIVE (ask m=Majambo to have a look in the parental and resident matrices - cor=1 (impossible))
#passive<-bigmama[,c("cryptoIF","pareIF","residIF","transiIF",
#                   "cryptoout","pareout","residout","transiout",
#                   "cryptobtw","parebtw","residbtw","transibtw",
#                   "cryptoIF","pareIF","residIF","transiIF",
#                   "cryptoLR","pareLR","residLR","transiLR")]

active<-bigmama[,c("ID","crypto5in","pare5in","resid15in","transi15in",
                    "crypto5out","pare5out","resid15out","transi15out",
                    "crypto5btw","pare5btw","resid15btw","transi15btw",
                    "crypto5IF","pare5IF","resid15IF","transi15IF",
                    "crypto5LR","pare5LR","resid15LR","transi15LR",
                   "crypto5IFS","pare5IFS","resid15IFS","transi15IFS",
                   "crypto5OF","pare5OF","resid15OF","transi15OF")]



####Quantiles ACTIVE
#95%-100% - HIGH (1)

#IFlow degree
active$crypto5IFCAT <- ifelse(active$crypto5IF < quantile(active$crypto5IF, prob=c(0.95)),  0,1)

active$pare5IFCAT <- ifelse(active$pare5IF < quantile(active$pare5IF, prob=c(0.95)),  0,1)

active$resid15IFCAT <- ifelse(active$resid15IF < quantile(active$resid15IF, prob=c(0.95)), 0,1)

active$transi15IFCAT <- ifelse(active$transi15IF < quantile(active$transi15IF, prob=c(0.95)), 0,1)

active$IFrank<-rowSums(active[,c("pare5IFCAT","crypto5IFCAT","resid15IFCAT","transi15IFCAT")])

#Indegree
active$crypto5inCAT <- ifelse(active$crypto5in < quantile(active$crypto5in, prob=c(0.95)),  0,1)

active$pare5inCAT <- ifelse(active$pare5in < quantile(active$pare5in, prob=c(0.95)),  0,1)

active$resid15inCAT <- ifelse(active$resid15in < quantile(active$resid15in, prob=c(0.95)), 0,1)

active$transi15inCAT <- ifelse(active$transi15in < quantile(active$transi15in, prob=c(0.95)), 0,1)

active$inrank<-rowSums(active[,c("pare5inCAT","crypto5inCAT","resid15inCAT","transi15inCAT")])

#betweenness
active$crypto5btwCAT <- ifelse(active$crypto5btw < quantile(active$crypto5btw, prob=c(0.95)),  0,1)

active$pare5btwCAT <- ifelse(active$pare5btw < quantile(active$pare5btw, prob=c(0.95)),  0,1)

active$resid15btwCAT <- ifelse(active$resid15btw < quantile(active$resid15btw, prob=c(0.95)), 0,1)

active$transi15btwCAT <- ifelse(active$transi15btw < quantile(active$transi15btw, prob=c(0.95)), 0,1)

active$btwrank<-rowSums(active[,c("pare5btwCAT","crypto5btwCAT","resid15btwCAT","transi15btwCAT")])


#Local Retention
active$crypto5LRCAT <- ifelse(active$crypto5LR < quantile(active$crypto5LR, prob=c(0.95)),  0,1)

active$pare5LRCAT <- ifelse(active$pare5LR < quantile(active$pare5LR, prob=c(0.95)),  0,1)

active$resid15LRCAT <- ifelse(active$resid15LR < quantile(active$resid15LR, prob=c(0.95)), 0,1)

active$transi15LRCAT <- ifelse(active$transi15LR < quantile(active$transi15LR, prob=c(0.95)), 0,1)

active$LRrank<-rowSums(active[,c("pare5LRCAT","crypto5LRCAT","resid15LRCAT","transi15LRCAT")])

#IF + SR
active$crypto5IFSCAT <- ifelse(active$crypto5IFS < quantile(active$crypto5IFS, prob=c(0.95)),  0,1)

active$pare5IFSCAT <- ifelse(active$pare5IFS < quantile(active$pare5IFS, prob=c(0.95)),  0,1)

active$resid15IFSCAT <- ifelse(active$resid15IFS < quantile(active$resid15IFS, prob=c(0.95)), 0,1)

active$transi15IFSCAT <- ifelse(active$transi15IFS < quantile(active$transi15IFS, prob=c(0.95)), 0,1)

active$IFSrank<-rowSums(active[,c("pare5IFSCAT","crypto5IFSCAT","resid15IFSCAT","transi15IFSCAT")])

#OutFLow
active$crypto5OFCAT <- ifelse(active$crypto5OF < quantile(active$crypto5OF, prob=c(0.95)),  0,1)

active$pare5OFCAT <- ifelse(active$pare5OF < quantile(active$pare5OF, prob=c(0.95)),  0,1)

active$resid15OFCAT <- ifelse(active$resid15OF < quantile(active$resid15OF, prob=c(0.95)), 0,1)

active$transi15OFCAT <- ifelse(active$transi15OF < quantile(active$transi15OF, prob=c(0.95)), 0,1)

active$OFrank<-rowSums(active[,c("pare5OFCAT","crypto5OFCAT","resid15OFCAT","transi15OFCAT")])

#Outdegree
active$crypto5outCAT <- ifelse(active$crypto5out < quantile(active$crypto5out, prob=c(0.95)),  0,1)

active$pare5outCAT <- ifelse(active$pare5out < quantile(active$pare5out, prob=c(0.95)),  0,1)

active$resid15outCAT <- ifelse(active$resid15out < quantile(active$resid15out, prob=c(0.95)), 0,1)

active$transi15outCAT <- ifelse(active$transi15out < quantile(active$transi15out, prob=c(0.95)), 0,1)

active$outrank<-rowSums(active[,c("pare5outCAT","crypto5outCAT","resid15outCAT","transi15outCAT")])

activeranks<-active[,c("pare5outCAT","crypto5outCAT","resid15outCAT","transi15outCAT",
                       "pare5inCAT","crypto5inCAT","resid15inCAT","transi15inCAT",
                       "pare5IFCAT","crypto5IFCAT","resid15IFCAT","transi15IFCAT",
                       "pare5LRCAT","crypto5LRCAT","resid15LRCAT","transi15LRCAT",
                       "pare5btwCAT","crypto5btwCAT","resid15btwCAT","transi15btwCAT",
                       "pare5OFCAT","crypto5OFCAT","resid15OFCAT","transi15OFCAT")]


head(activeranks)
summary(activeranks)
#CORRELATION AMONG HOTSPOTS 
correlation_matrix <- cor(activeranks)
colors <- brewer.pal(n = 3, name = "RdYlBu")
p <- ggcorrplot(correlation_matrix , type = "upper", hc.order = TRUE, colors = brewer.pal(n = 3, name = "RdYlBu"))
p <- p + scale_fill_gradient2(limit = c(0.5,1), low = "blue", high = "red", mid="orange", midpoint = 0.95)
p


#MAP hotspots----
#Plotting MAP
library("ggmap")
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())
library(maptools)
library(maps)
mp <- NULL
coord<-read.csv(here("_data","IDs.csv"))
head(coord)
#View(coord)
##plotting biomass #
globalbiot<-merge(active,coord[,c("ID","Lon","Lat")], by ="ID",all=F)
summary(globalbiot)
#PLOT- HOTSPOTS
mapWorld <- borders("world", colour="gray", fill="gray") # create a layer of borders

mp <- ggplot() +   mapWorld + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank()) +
  coord_map(ylim=c(-50, 50), xlim = c(-200,200))
####
#GROUPS total congruence
rm(perc)
perc <- round(100*length(which(globalbiot$IFrank == 3))/length(globalbiot$IFrank),2)
perc
mpIF <- mp + geom_point(aes(x=globalbiot$Lon, y=globalbiot$Lat, fill=as.factor(globalbiot$IFrank)), 
                       alpha=0.5, size=1, shape=21, colour="white") +
  scale_fill_manual(values=c("blue","blue","blue","blue","red")) + labs(x="",y="") + theme(legend.position = "none") +
  geom_point(aes(x=globalbiot[globalbiot$IFrank == 3,]$Lon, y=globalbiot[globalbiot$IFrank == 3,]$Lat, 
                 fill=as.factor(globalbiot[globalbiot$IFrank == 3,]$IFrank)), 
            alpha=0.5, size=1, shape=21, colour="red")


rm(perc)
perc <- round(100*length(which(globalbiot$inrank == 3))/length(globalbiot$inrank),2)
perc
mpin <- mp + geom_point(aes(x=globalbiot$Lon, y=globalbiot$Lat, fill=as.factor(globalbiot$inrank)), 
                        alpha=0.5, size=1, shape=21, colour="white") +
  scale_fill_manual(values=c("blue","blue","blue","blue","red")) + labs(x="",y="") + theme(legend.position = "none") +
  geom_point(aes(x=globalbiot[globalbiot$inrank == 3,]$Lon, y=globalbiot[globalbiot$inrank == 3,]$Lat, 
                 fill=as.factor(globalbiot[globalbiot$inrank == 3,]$inrank)), 
             alpha=0.5, size=1, shape=21, colour="red")


rm(perc)
perc <- round(100*length(which(globalbiot$btwrank == 3))/length(globalbiot$btwrank),2)
perc
mpbtw <- mp + geom_point(aes(x=globalbiot$Lon, y=globalbiot$Lat, fill=as.factor(globalbiot$btwrank)), 
                        alpha=0.5, size=1, shape=21, colour="white") +
  scale_fill_manual(values=c("blue","blue","blue","blue","red")) + labs(x="",y="") + theme(legend.position = "none") +
  geom_point(aes(x=globalbiot[globalbiot$btwrank == 3,]$Lon, y=globalbiot[globalbiot$btwrank == 3,]$Lat, 
                 fill=as.factor(globalbiot[globalbiot$btwrank == 3,]$btwrank)), 
             alpha=0.5, size=1, shape=21, colour="red")

rm(perc)
perc <- round(100*length(which(globalbiot$LRrank == 3))/length(globalbiot$LRrank),2)
perc
mpLR <- mp + geom_point(aes(x=globalbiot$Lon, y=globalbiot$Lat, fill=as.factor(globalbiot$LRrank)), 
                         alpha=0.5, size=1, shape=21, colour="white") +
  scale_fill_manual(values=c("blue","blue","blue","blue","red")) + labs(x="",y="") + theme(legend.position = "none") +
  geom_point(aes(x=globalbiot[globalbiot$LRrank == 3,]$Lon, y=globalbiot[globalbiot$LRrank == 3,]$Lat, 
                 fill=as.factor(globalbiot[globalbiot$LRrank == 3,]$LRrank)), 
             alpha=0.5, size=1, shape=21, colour="red")

rm(perc)
perc <- round(100*length(which(globalbiot$OFrank == 3))/length(globalbiot$OFrank),2)
perc
mpOF <- mp + geom_point(aes(x=globalbiot$Lon, y=globalbiot$Lat, fill=as.factor(globalbiot$OFrank)), 
                        alpha=0.5, size=1, shape=21, colour="white") +
  scale_fill_manual(values=c("blue","blue","blue","blue","red")) + labs(x="",y="") + theme(legend.position = "none") +
  geom_point(aes(x=globalbiot[globalbiot$OFrank == 3,]$Lon, y=globalbiot[globalbiot$OFrank == 3,]$Lat, 
                 fill=as.factor(globalbiot[globalbiot$OFrank == 3,]$OFrank)), 
             alpha=0.5, size=1, shape=21, colour="red")

rm(perc)
perc <- round(100*length(which(globalbiot$outrank == 3))/length(globalbiot$outrank),2)
perc
mpout <- mp + geom_point(aes(x=globalbiot$Lon, y=globalbiot$Lat, fill=as.factor(globalbiot$outrank)), 
                        alpha=0.5, size=1, shape=21, colour="white") +
  scale_fill_manual(values=c("blue","blue","blue","blue","red")) + labs(x="",y="") + theme(legend.position = "none") +
  geom_point(aes(x=globalbiot[globalbiot$outrank == 3,]$Lon, y=globalbiot[globalbiot$outrank == 3,]$Lat, 
                 fill=as.factor(globalbiot[globalbiot$outrank == 3,]$outrank)), 
             alpha=0.5, size=1, shape=21, colour="red")

ggarrange(mpIF,mpin,mpOF,mpout,mpLR,mpbtw,labels=c("Inflow","Indegree","Outflow",
                                                   "Outdegree","Self-recruitment",
                                                   "Betweenness"),nrow=3,ncol=2)



#BIOMASS, RICHNESS ~ CONNECTIVITY---- 
#Overlay points/ Preparing data ----
nodesID<-read.csv(here("_data","IDs.csv"),h=T)
head(nodesID) #14804 id's
colnames(nodesID) <- c("ID","ID2","lon","lat","territory","other")

dataBIC<-read.csv(here("_data","databiomassFull.csv"))

dt<-match_nrst_haversine(dataBIC$lat, dataBIC$lon, nodesID$lat, nodesID$lon, nodesID$ID,
                         close_enough = 0.1)
dataBIC$ID<-dt$pos
dataBIC$distprox<-dt$dist
head(dataBIC)
firstat<-merge(dataBIC,active,by=c("ID"),all=F)
dim(firstat) #452 sites
filtertre<-firstat[firstat$distprox<10,] #273 sites -> apply this one (only biomass sites inside reef cells) 

#Sum metrics
filtertre$sumIN<-rowSums(filtertre[,c("crypto5in","pare5in","resid15in","transi15in")])
filtertre$sumIF<-rowSums(filtertre[,c("crypto5IF","pare5IF","resid15IF","transi15IF")])
filtertre$sumLR<-rowSums(filtertre[,c("crypto5LR","pare5LR","resid15LR","transi15LR")])
filtertre$sumbtw<-rowSums(filtertre[,c("crypto5btw","pare5btw","resid15btw","transi15btw")])


filtertre$LogB <- log1p(filtertre$biomassarea)
filtertre$LogG <- log1p(filtertre$grav_total)
filtertre$LogIn <- log1p(filtertre$sumIN)
filtertre$inrank <-as.factor(filtertre$inrank)
filtertre$LogIF <- log1p(filtertre$sumIF) 
filtertre$IFrank <-as.factor(filtertre$IFrank)
filtertre$Logbtw <- log1p(filtertre$sumbtw) 
filtertre$btwrank <-as.factor(filtertre$btwrank)
filtertre$LogLR <- log1p(filtertre$sumLR) 
filtertre$LRrank <-as.factor(filtertre$LRrank)
filtertre$logGrav <- log1p(filtertre$grav_total)

#MODELS ----

bioric<-gam(LogB  ~ richness, family="poisson",data=filtertre)
visreg::visreg(bioric)
summary(bioric) #r=0.06

#Species richness
test1R<-gam(richness ~ Class*LogG + temp, family="poisson",data=filtertre)
test2R<-gam(richness ~ Class*LogG + temp + LogIF + LogIn +
            LogLR + Logbtw   + LogLR, family="poisson",data=filtertre)
test0R<-gam(richness ~ Class*LogG + temp  + LogIF  + LogIn +
              LogLR + Logbtw , family="poisson",data=filtertre)
vif(test1R) #huge VIF

test3R<-gam(richness ~ Class*LogG + temp + LogIn , family="poisson",data=filtertre)
test4R<-gam(richness ~ Class*LogG + temp + LogLR, family="poisson",data=filtertre)
test5R<-gam(richness ~ Class*LogG + temp + Logbtw, family="poisson",data=filtertre) #btw 1st
test6R<-gam(richness ~ Class*LogG + temp+ LogIF, family="poisson",data=filtertre)

model.sel(test0R,test1R,test3R,test4R,test5R,test6R) # always compare nested models ,removing variable with smaller AIC:
#5, 3 and 6 (all the models with conectivity improve AIC) - btw better explain diversity
summary(test0R) # 26.9% explained deviance
visreg::visreg(test0R) 
summary(test1R) # 16.2
summary(test2R)
summary(test5R) #10% improvement in Dev.Expl  25.3
summary(test3R) #9%  # 25.4
anova(test0R)

visreg::visreg(test5RN,"btwrank")
visreg::visreg(test3RN,"inrank")

test5RN<-gam(richness ~ LogG + temp + Class + Logbtw + btwrank,data=filtertre) #btw 1st
test3RN<-gam(richness ~ LogG + temp + Class + LogIn + inrank,data=filtertre)


summary(test1R) #14%
summary(test5R) #24%
summary(test3R) #24%


plot(log(fitted(test5R)), residuals(test5R), xlab = "Fitted Values (log)", ylab = "Residuals", las=1)
abline(h=0, lty=2 ,col='red')
qqPlot(resid(test5R), grid=F, col.lines="grey", ylab="", las=1) #ok

par(mfrow=c(1,2))
visreg::visreg(test5RN,"Logbtw")
visreg::visreg(test3RN,"LogIn")


#BIOMASS
testB1<-gam(LogB ~ LogG + temp + Class + richness,data=filtertre)
testB2<-gam(LogB ~ LogG + temp + Class + LogIF + IFrank + LogIn*inrank +
              LogLR*LRrank + Logbtw + btwrank + richness,data=filtertre)
vif(testB2) #huge VIF

testB3<-gam(LogB ~ LogG + temp + Class + LogIn * inrank+ richness,data=filtertre) #less worse
testB4<-gam(LogB ~ LogG + temp + Class + LogLR * LRrank+ richness,data=filtertre)
testB5<-gam(LogB ~ LogG + temp + Class + Logbtw * btwrank+ richness,data=filtertre)
testB6<-gam(LogB ~ LogG + temp + Class + LogIF * IFrank+ richness,data=filtertre)

model.sel(testB1,testB3,testB4,testB5,testB6) # removing variable with smaller AIC:
#connectivity does not explain biomass
summary(testB3) #14%

par(mfrow=c(2,2))
visreg::visreg(test5RN,"Logbtw")
visreg::visreg(test3RN,"LogIn")

visreg::visreg(testB5,"Logbtw")
visreg::visreg(testB3,"LogIn")





