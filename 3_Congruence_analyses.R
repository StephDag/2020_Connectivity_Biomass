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
library(tidyr)
library(plyr)
library(tidyverse)

#1st part 
#Connectivity Rank----
bigmama %>% rm()
bigmama<-read.csv("~/Documents/GIt/2020_Connectivity_Biomass/_data/global_metrics.csv",h=T)
summary(bigmama)


#PASSIVE

passive<-bigmama[,c("cryptoIF","pareIF","residIF","transiIF",
                   "cryptoout","pareout","residout","transiout",
                   "cryptobtw","parebtw","residbtw","transibtw",
                   "cryptoOF","pareOF","residOF","transiOF",
                   "cryptoLR","pareLR","residLR","transiLR")]

dim(passive)

active<-bigmama[,c("ID","crypto5in","pare5in","resid15in","transi15in",
                    "crypto5out","pare5out","resid15out","transi15out",
                    "crypto5btw","pare5btw","resid15btw","transi15btw",
                    "crypto5IF","pare5IF","resid15IF","transi15IF",
                    "crypto5LR","pare5LR","resid15LR","transi15LR",
                   "crypto5IFS","pare5IFS","resid15IFS","transi15IFS",
                   "crypto5OF","pare5OF","resid15OF","transi15OF")]


#passive corraltion
correlation_matrix <- cor(active[,2:29])
colors <- brewer.pal(n = 3, name = "RdYlBu")
p <- ggcorrplot(correlation_matrix , type = "upper", hc.order = TRUE, colors = brewer.pal(n = 3, name = "RdYlBu"))
p <- p + scale_fill_gradient2(limit = c(0.5,1), low = "blue", high = "red", mid="orange", midpoint = 0.75)
p

#active correlation
correlation_matrix <- cor(passive[,2:20])
colors <- brewer.pal(n = 3, name = "RdYlBu")
p <- ggcorrplot(correlation_matrix , type = "upper", hc.order = TRUE, colors = brewer.pal(n = 3, name = "RdYlBu"))
p <- p + scale_fill_gradient2(limit = c(0.5,1), low = "blue", high = "red", mid="orange", midpoint = 0.75)
p


####Quantiles ACTIVE
#75%-100% - HIGH (1)

#IFlow degree
active$crypto5IFCAT <- ifelse(active$crypto5IF < quantile(active$crypto5IF, prob=c(0.75)),  0,1)

active$pare5IFCAT <- ifelse(active$pare5IF < quantile(active$pare5IF, prob=c(0.75)),  0,1)

active$resid15IFCAT <- ifelse(active$resid15IF < quantile(active$resid15IF, prob=c(0.75)), 0,1)

active$transi15IFCAT <- ifelse(active$transi15IF < quantile(active$transi15IF, prob=c(0.75)), 0,1)

active$IFrank<-rowSums(active[,c("pare5IFCAT","crypto5IFCAT","resid15IFCAT","transi15IFCAT")])

#Indegree
active$crypto5inCAT <- ifelse(active$crypto5in < quantile(active$crypto5in, prob=c(0.75)),  0,1)

active$pare5inCAT <- ifelse(active$pare5in < quantile(active$pare5in, prob=c(0.75)),  0,1)

active$resid15inCAT <- ifelse(active$resid15in < quantile(active$resid15in, prob=c(0.75)), 0,1)

active$transi15inCAT <- ifelse(active$transi15in < quantile(active$transi15in, prob=c(0.75)), 0,1)

active$inrank<-rowSums(active[,c("pare5inCAT","crypto5inCAT","resid15inCAT","transi15inCAT")])

#betweenness
active$crypto5btwCAT <- ifelse(active$crypto5btw < quantile(active$crypto5btw, prob=c(0.75)),  0,1)

active$pare5btwCAT <- ifelse(active$pare5btw < quantile(active$pare5btw, prob=c(0.75)),  0,1)

active$resid15btwCAT <- ifelse(active$resid15btw < quantile(active$resid15btw, prob=c(0.75)), 0,1)

active$transi15btwCAT <- ifelse(active$transi15btw < quantile(active$transi15btw, prob=c(0.75)), 0,1)

active$btwrank<-rowSums(active[,c("pare5btwCAT","crypto5btwCAT","resid15btwCAT","transi15btwCAT")])


#Local Retention
active$crypto5LRCAT <- ifelse(active$crypto5LR < quantile(active$crypto5LR, prob=c(0.75)),  0,1)

active$pare5LRCAT <- ifelse(active$pare5LR < quantile(active$pare5LR, prob=c(0.75)),  0,1)

active$resid15LRCAT <- ifelse(active$resid15LR < quantile(active$resid15LR, prob=c(0.75)), 0,1)

active$transi15LRCAT <- ifelse(active$transi15LR < quantile(active$transi15LR, prob=c(0.75)), 0,1)

active$LRrank<-rowSums(active[,c("pare5LRCAT","crypto5LRCAT","resid15LRCAT","transi15LRCAT")])

#IF + SR
active$crypto5IFSCAT <- ifelse(active$crypto5IFS < quantile(active$crypto5IFS, prob=c(0.75)),  0,1)

active$pare5IFSCAT <- ifelse(active$pare5IFS < quantile(active$pare5IFS, prob=c(0.75)),  0,1)

active$resid15IFSCAT <- ifelse(active$resid15IFS < quantile(active$resid15IFS, prob=c(0.75)), 0,1)

active$transi15IFSCAT <- ifelse(active$transi15IFS < quantile(active$transi15IFS, prob=c(0.75)), 0,1)

active$IFSrank<-rowSums(active[,c("pare5IFSCAT","crypto5IFSCAT","resid15IFSCAT","transi15IFSCAT")])

#OutFLow
active$crypto5OFCAT <- ifelse(active$crypto5OF < quantile(active$crypto5OF, prob=c(0.75)),  0,1)

active$pare5OFCAT <- ifelse(active$pare5OF < quantile(active$pare5OF, prob=c(0.75)),  0,1)

active$resid15OFCAT <- ifelse(active$resid15OF < quantile(active$resid15OF, prob=c(0.75)), 0,1)

active$transi15OFCAT <- ifelse(active$transi15OF < quantile(active$transi15OF, prob=c(0.75)), 0,1)

active$OFrank<-rowSums(active[,c("pare5OFCAT","crypto5OFCAT","resid15OFCAT","transi15OFCAT")])

#Outdegree
active$crypto5outCAT <- ifelse(active$crypto5out < quantile(active$crypto5out, prob=c(0.75)),  0,1)

active$pare5outCAT <- ifelse(active$pare5out < quantile(active$pare5out, prob=c(0.75)),  0,1)

active$resid15outCAT <- ifelse(active$resid15out < quantile(active$resid15out, prob=c(0.75)), 0,1)

active$transi15outCAT <- ifelse(active$transi15out < quantile(active$transi15out, prob=c(0.75)), 0,1)

active$outrank<-rowSums(active[,c("pare5outCAT","crypto5outCAT","resid15outCAT","transi15outCAT")])

activeranks<-active[,c("pare5outCAT","crypto5outCAT","resid15outCAT","transi15outCAT",
                       "pare5inCAT","crypto5inCAT","resid15inCAT","transi15inCAT",
                       "pare5IFCAT","crypto5IFCAT","resid15IFCAT","transi15IFCAT",
                       "pare5LRCAT","crypto5LRCAT","resid15LRCAT","transi15LRCAT",
                       "pare5btwCAT","crypto5btwCAT","resid15btwCAT","transi15btwCAT",
                       "pare5OFCAT","crypto5OFCAT","resid15OFCAT","transi15OFCAT")]


#CORRELATION AMONG HOTSPOTS 
correlation_matrix <- cor(activeranks)
colors <- brewer.pal(n = 3, name = "RdYlBu")
p <- ggcorrplot(correlation_matrix , type = "upper", hc.order = TRUE, colors = brewer.pal(n = 3, name = "RdYlBu"))
p <- p + scale_fill_gradient2(limit = c(0.5,1), low = "blue", high = "red", mid="orange", midpoint = 0.75)
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
coord<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_on_going/Connectivity analyses/IDs.csv")
head(coord)
#View(coord)
##plotting biomass #
globalbiot<-merge(active,coord[,c("ID","Lon","Lat")], by ="ID",all=F)

#PLOT- HOTSPOTS
mapWorld <- borders("world", colour="gray", fill="gray") # create a layer of borders

mp <- ggplot() +   mapWorld + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank()) +
  coord_map(ylim=c(-50, 50), xlim = c(-200,200))
####
#GROUPS

mpIF <- mp + geom_point(aes(x=globalbiot$Lon, y=globalbiot$Lat, fill=as.factor(globalbiot$IFrank)), 
                       alpha=0.5, size=1, shape=21, colour="white") +
  scale_fill_manual(values=c("blue","blue","blue","blue","red")) + labs(x="",y="") + theme(legend.position = "none") +
  geom_point(aes(x=globalbiot[globalbiot$IFrank == 4,]$Lon, y=globalbiot[globalbiot$IFrank == 4,]$Lat, 
                 fill=as.factor(globalbiot[globalbiot$IFrank == 4,]$IFrank)), 
            alpha=0.5, size=1, shape=21, colour="red")


mpin <- mp + geom_point(aes(x=globalbiot$Lon, y=globalbiot$Lat, fill=as.factor(globalbiot$inrank)), 
                        alpha=0.5, size=1, shape=21, colour="white") +
  scale_fill_manual(values=c("blue","blue","blue","blue","red")) + labs(x="",y="") + theme(legend.position = "none") +
  geom_point(aes(x=globalbiot[globalbiot$inrank == 4,]$Lon, y=globalbiot[globalbiot$inrank == 4,]$Lat, 
                 fill=as.factor(globalbiot[globalbiot$inrank == 4,]$inrank)), 
             alpha=0.5, size=1, shape=21, colour="red")


mpbtw <- mp + geom_point(aes(x=globalbiot$Lon, y=globalbiot$Lat, fill=as.factor(globalbiot$btwrank)), 
                        alpha=0.5, size=1, shape=21, colour="white") +
  scale_fill_manual(values=c("blue","blue","blue","blue","red")) + labs(x="",y="") + theme(legend.position = "none") +
  geom_point(aes(x=globalbiot[globalbiot$btwrank == 4,]$Lon, y=globalbiot[globalbiot$btwrank == 4,]$Lat, 
                 fill=as.factor(globalbiot[globalbiot$btwrank == 4,]$btwrank)), 
             alpha=0.5, size=1, shape=21, colour="red")


mpLR <- mp + geom_point(aes(x=globalbiot$Lon, y=globalbiot$Lat, fill=as.factor(globalbiot$LRrank)), 
                         alpha=0.5, size=1, shape=21, colour="white") +
  scale_fill_manual(values=c("blue","blue","blue","blue","red")) + labs(x="",y="") + theme(legend.position = "none") +
  geom_point(aes(x=globalbiot[globalbiot$LRrank == 4,]$Lon, y=globalbiot[globalbiot$LRrank == 4,]$Lat, 
                 fill=as.factor(globalbiot[globalbiot$LRrank == 4,]$LRrank)), 
             alpha=0.5, size=1, shape=21, colour="red")

mpOF <- mp + geom_point(aes(x=globalbiot$Lon, y=globalbiot$Lat, fill=as.factor(globalbiot$OFrank)), 
                        alpha=0.5, size=1, shape=21, colour="white") +
  scale_fill_manual(values=c("blue","blue","blue","blue","red")) + labs(x="",y="") + theme(legend.position = "none") +
  geom_point(aes(x=globalbiot[globalbiot$OFrank == 4,]$Lon, y=globalbiot[globalbiot$OFrank == 4,]$Lat, 
                 fill=as.factor(globalbiot[globalbiot$OFrank == 4,]$OFrank)), 
             alpha=0.5, size=1, shape=21, colour="red")

mpout <- mp + geom_point(aes(x=globalbiot$Lon, y=globalbiot$Lat, fill=as.factor(globalbiot$outrank)), 
                        alpha=0.5, size=1, shape=21, colour="white") +
  scale_fill_manual(values=c("blue","blue","blue","blue","red")) + labs(x="",y="") + theme(legend.position = "none") +
  geom_point(aes(x=globalbiot[globalbiot$outrank == 4,]$Lon, y=globalbiot[globalbiot$outrank == 4,]$Lat, 
                 fill=as.factor(globalbiot[globalbiot$outrank == 4,]$outrank)), 
             alpha=0.5, size=1, shape=21, colour="red")

ggarrange(mpIF,mpin,mpOF,mpout,mpLR,mpbtw,labels=c("Inflow","Indegree","Outflow",
                                                   "Outdegree","Self-recruitment",
                                                   "Betweenness"),nrow=3,ncol=2)





#BIOMASS, RICHNESS ~ CONNECTIVITY---- 

#Preparing data ----
active<-bigmama[,c("ID","crypto5IF","pare5IF","resid15IF","transi15IF",
                   "crypto5in","pare5in","resid15in","transi15in",
                   "crypto5OF","pare5OF","resid15OF","transi15OF")]

active$totalIFA<-rowSums(active[,c("crypto5IF","pare5IF","resid15IF","transi15IF")])
active$totalinA<-rowSums(active[,c("crypto5in","pare5in","resid15in","transi15in")])
active$totalOFA<-rowSums(active[,c("crypto5OF","pare5OF","resid15OF","transi15OF")])


passive<-bigmama[,c("ID","cryptoIF","pareIF","residIF","transiIF",
                    "cryptoOF","pareOF","residOF","transiOF")]

passive$totalIF<-rowSums(passive[,c("cryptoIF","pareIF","residIF","transiIF")])
passive$totalOF<-rowSums(passive[,c("cryptoOF","pareOF","residOF","transiOF")])


nodesID<-read.csv("~/Documents/GIt/2020_Connectivity_Biomass/_data/IDs.csv",h=T)
head(nodesID) #14804 id's
colnames(nodesID) <- c("ID","ID2","lon","lat","territory","other")

dataBIC<-read.csv("~/Documents/GIt/2020_Connectivity_Biomass/_data/databiomassFull.csv",h=T)
dt<-match_nrst_haversine(dataBIC$lat, dataBIC$lon, nodesID$lat, nodesID$lon, nodesID$ID,
                         close_enough = 0.1)
dataBIC$ID<-dt$pos
dataBIC$distprox<-dt$dist
head(dataBIC)
firstat<-merge(dataBIC,active,by=c("ID"),all=F)
dim(firstat) #452 sites
filtertre<-firstat[firstat$distprox<10,] #272 (3 on same ID) sites -> apply this one (only biomass sites inside reef cells) 
filtertre<-left_join(filtertre,passive,by="ID")
summary(filtertre)
#Sum metrics
filtertre$LogB <- log1p(filtertre$biomassarea1)
filtertre$LogG <- log1p(filtertre$grav_total)

#Delta outflow
head(filtertre)
#Transform g/m to Kg/hectare (*10)
#Delta - real/potential OF weighted by biomass (sum OF)
filtertre$deltaOF<-((10*filtertre$biomassarea1)*filtertre$totalOF) - (1000*filtertre$totalOF)
filtertre$deltaOFA<-((10*filtertre$biomassarea1)*filtertre$totalOFA) - (1000*filtertre$totalOFA)
hist(filtertre$deltaOFA)


#MODELS ----
#1st.	How human impact affect larval connectivity?
filtertre$LogOF<-log1p(filtertre$totalOF)
filtertre$LogOFA<-log1p(filtertre$totalOFA)

#Relationship between Outflow, biomass and gravity
test1<-gam(LogOF ~ LogG, data=filtertre)
test1A<-gam(LogOFA ~ LogG, data=filtertre)
summary(test1A) 
par(mfrow=c(1,2))
visreg::visreg(test1) #no relatinship between Outflow and Gravity
visreg::visreg(test1A) #no relatinship between Outflow active model and Gravity
#modelBio<-gam(LogB ~ LogG + Class + LogOFA, data=filtertre)
#modelOF<-gam(LogOFA ~ LogG + Class + LogB, data=filtertre)
#visreg::visreg(modelBio,"LogG",by="Class")
#visreg::visreg(modelOF,"LogG")
#visreg::visreg(modelOF,"LogB","Class")


#Delta outflow (little difference between active and passive)
#passive model Delta
outflowH1<-gam(deltaOF ~ LogG + Class, data=filtertre)
summary(outflowH1) #Deviance explained 11%
#active model Delta
outflowAH1<-gam(deltaOFA ~ LogG * Class, data=filtertre) #14.4
outflowAH1B<-gam(deltaOFA ~ LogG + Class, data=filtertre)

model.sel(outflowAH1,outflowAH1B) #with an interaction - deltaAIC(2.11)
visreg::visreg(outflowAH1,"Class")
visreg::visreg(outflowAH1,"LogG",by="Class")
visreg::visreg(outflowAH1,"LogG",by="Class",overlay=T)
abline(h=0,lty=2)

#Does larval connectivity mitigate human impacts on reefs? (Full model?)
#Quantile regression (.75 and .25) = Logbiomass ~ Inflow(0.75high;0.25low) + LogG + Class (test for Active and Passive)



