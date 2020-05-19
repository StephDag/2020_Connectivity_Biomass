#---Metadata
#ID - node ID from global connectivity matrix 
#Region - Biogeographic region of each site/locality/site (site|locality|region)
#temp - temprature (Methods from Barneche et al 2019 - GEB)
#Model Mode  - 12 models based on different reproductive strategies and larvae life-history traits (PLD, sensory system, release)
#Larval_behavior - active (larvae with sensory system) or passive (larvae can only settle at the end of its PLD) models
#FE - Functional entity

#Response variables
#Richness - species richness for each site (estimation based on Barneche et al 2019 - GEB)
#biomassarea1 - reef fish biomass weiigh/unit of area 

#Explanatory variables
#(BR after variable name means that connections between reefs were filtered by common biogeographical region)
#This means, using Kulbicki et al 2013, we have assumed 0 probability of connection between reefs
#from distinct biogeographical regions once they do not share reef fish species. 

#INTRINSIC FATORS 
#grav-total - distance of each site from nearest market (Cinner et al 2014 - PNAS)
#Class - fished, restricted gear of fully protected reef area
#Age_of_protection - age of protection (if any) of each MPA/localy managed reef area 
#Indegree - number of connections
#btwdegree - the number of particules passing trought the node - Ecological corridor of larvae
#Inflow - number of particules arriving from other reef cells
#Outdegree - the numbr of outward connections from a reef cell
#InflowLR - self-recruitment + number of particles arriving from other reef cells
#SelfR - self-recruitment
#IndegreeBR - the number of inward connections from cells
#IndegreeMPABR - the number of inward connections from cells inside locally managed areas

#Outdegree - the number of outward connections (seeding reefs)
#Outflow - the number of particles reef cells are seeding other reefs

#EXTRINSIC FACTORS (NEIGHBOURS - POTENTIAL LARVAL SOURCE)
#INdegreeNEiBR - the average inward number of connections of neighbours 
#InflowNEiBR - the average inward number of particules arriving to neighbours 
#CorridorIndegreeBR - the number of connections with corridors (high corridor (quantile(>0.75)))
#grav_neiBR <- average gravity of neighbours


###Setting factors for each reef (matrix node)###
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(maptools)
library(maps)
library(hutilscpp) 

#Prparing database
nodesprop<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Grav-Ecoregions.csv",h=T)
kulbickiregions<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/ecoregions_kulbicki.csv",h=T)
colnames(kulbickiregions) <- c("N","PROVINCE", "Kulbicki")
nodesprop<-left_join(nodesprop,kulbickiregions[,2:3],by="PROVINCE")

#Level Protection
nodesIDMPA<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/IDs_with_MPAs.csv",h=T)
nodesIDMPA<-nodesIDMPA[,c("ID","Lon","Lat","Class","General")] #genral is Locally managed (MPA or any gear restriction vs fished/no restrictions)
nodesFinal<-left_join(nodesprop,nodesIDMPA[,c("ID","Class","General")],by="ID")
head(nodesFinal) #node Id, coordinates, gravity, ecoregion, province, bioregion, class, and general

bigmama<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/global_metrics.csv",h=T)
global_factors<-left_join(bigmama,nodesFinal,by="ID")
#MPA/locally managed (1) or lack of management (0)
global_factors$MPA<-ifelse(global_factors$General == "MPA",1,0)
head(global_factors)

#To set whether a point is a high corridor (quantile(>0.75))
#Crypto ("cryptobtw",crypto5btw,"crypto15btw)

global_factors$corridorcrypto<-ifelse(global_factors$cryptobtw > quantile(global_factors$cryptobtw,0.75), 1,0)
global_factors$corridorcrypto5<-ifelse(global_factors$crypto5btw > quantile(global_factors$crypto5btw,0.75), 1,0)
global_factors$corridorcrypto15<-ifelse(global_factors$crypto15btw > quantile(global_factors$crypto15btw,0.75), 1,0)


#pare ("parebtw",pare5btw,"pare15btw,"pare35btw")

global_factors$corridorpare<-ifelse(global_factors$parebtw > quantile(global_factors$parebtw,0.75), 1,0)
global_factors$corridorpare5<-ifelse(global_factors$pare5btw > quantile(global_factors$pare5btw,0.75), 1,0)
global_factors$corridorpare15<-ifelse(global_factors$pare15btw > quantile(global_factors$pare15btw,0.75), 1,0)

#resid ("residbtw",resid5btw,"resid15btw,"resid35btw")

global_factors$corridorresid<-ifelse(global_factors$residbtw > quantile(global_factors$residbtw,0.75), 1,0)
global_factors$corridorresid15<-ifelse(global_factors$resid15btw > quantile(global_factors$resid15btw,0.75), 1,0)
global_factors$corridorresid35<-ifelse(global_factors$resid35btw > quantile(global_factors$resid35btw,0.75), 1,0)


#transi ("transibtw","transi15btw,"transi35btw")

global_factors$corridortransi<-ifelse(global_factors$transibtw > quantile(global_factors$transibtw,0.75), 1,0)
global_factors$corridortransi15<-ifelse(global_factors$transi15btw > quantile(global_factors$transi15btw,0.75), 1,0)
global_factors$corridortransi35<-ifelse(global_factors$transi35btw > quantile(global_factors$transi35btw,0.75), 1,0)

#Connections crypto
myconnections<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Connectivity/crypto.csv")
head(myconnections)
dataconec<-myconnections[,2:4]
colnames(dataconec) <- c("From","To","Inflow")
head(dataconec)

#Kulbicki
crypto_BioRegion <- dataconec %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('From' = 'ID'), all=F) %>%
  dplyr::rename(Source_Region = Kulbicki) %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('To' = 'ID'), all=F) %>%
  dplyr::rename(Sink_Region = Kulbicki) %>% 
  inner_join(nodesFinal %>% dplyr::select(ID, General), by = c('From' = 'ID'), all=F) %>%
  inner_join(global_factors %>% dplyr::select(ID, cryptoin), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Indegree = cryptoin) %>%
  inner_join(global_factors %>% dplyr::select(ID, cryptoIF), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Inflow = cryptoIF) %>% 
  inner_join(global_factors %>% dplyr::select(ID, corridorcrypto), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(CorridorIN = corridorcrypto) %>% 
  inner_join(global_factors %>% dplyr::select(ID, Grav_tot), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(GravNei = Grav_tot)
  
head(crypto_BioRegion)

#filter connections within Biogeo regions
crypto_BioRegion<-crypto_BioRegion[crypto_BioRegion$Source_Region == crypto_BioRegion$Sink_Region,]
crypto_BioRegion$Indegree <- ifelse(crypto_BioRegion$Inflow > 0,1,0)
#
cryptoInflowBR <- aggregate(Inflow ~ To,crypto_BioRegion,sum)
cryptoIDegreeBR <- aggregate(Indegree ~ To, crypto_BioRegion,sum)
cryptoCorridorIndegreeBR <- aggregate(CorridorIN ~ To, crypto_BioRegion,sum)
cryptograv_neiBR <- aggregate(GravNei ~ To, crypto_BioRegion,sum)

cryptoBR<-left_join(cryptoInflowBR,cryptoIDegreeBR, by="To")
cryptoBR<-left_join(cryptoBR,cryptoCorridorIndegreeBR, by="To")
cryptoBR<-left_join(cryptoBR,cryptograv_neiBR, by="To")
head(cryptoBR)

#MPA Indegree/Inflow
crypto_BioRegion$InflowMPA<-ifelse(crypto_BioRegion$General == "MPA",crypto_BioRegion$Inflow,0)
crypto_BioRegion$IndegreeMPA<-ifelse(crypto_BioRegion$General == "MPA",crypto_BioRegion$Indegree,0)
cryptoInflowMPABR <- aggregate(InflowMPA ~ To,crypto_BioRegion,sum)
cryptoIDegreeMPABR <- aggregate(IndegreeMPA ~ To, crypto_BioRegion,sum)
cryptoBR<-left_join(cryptoBR,cryptoIDegreeMPABR, by="To")
cryptoBR<-left_join(cryptoBR,cryptoInflowMPABR, by="To")
head(cryptoBR)

#Average indegree of neighbours/
cryptoInflowNBBR <- aggregate(Source_Indegree ~ To,crypto_BioRegion,mean)
cryptoIDegreeNBBR <- aggregate(Source_Inflow ~ To, crypto_BioRegion,mean)
cryptoBR<-left_join(cryptoBR,cryptoInflowNBBR, by="To")
cryptoBR<-left_join(cryptoBR,cryptoIDegreeNBBR, by="To")
dim(cryptoBR) #merge this data later with other ones
cryptoBR$ModelMode <- ifelse(cryptoBR$Inflow > 0, "crypto","")
head(cryptoBR)

#--------- Crypto 5
#Connections crypto5
myconnections<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Connectivity/crypto5.csv")
head(myconnections)
dataconec<-myconnections[,2:4]
colnames(dataconec) <- c("From","To","Inflow")
head(dataconec)

#Kulbicki
crypto5_BioRegion <- dataconec %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('From' = 'ID'), all=F) %>%
  dplyr::rename(Source_Region = Kulbicki) %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('To' = 'ID'), all=F) %>%
  dplyr::rename(Sink_Region = Kulbicki) %>% 
  inner_join(nodesFinal %>% dplyr::select(ID, General), by = c('From' = 'ID'), all=F) %>%
  inner_join(global_factors %>% dplyr::select(ID, crypto5in), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Indegree = crypto5in) %>%
  inner_join(global_factors %>% dplyr::select(ID, crypto5IF), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Inflow = crypto5IF) %>% 
  inner_join(global_factors %>% dplyr::select(ID, corridorcrypto5), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(CorridorIN = corridorcrypto5) %>% 
  inner_join(global_factors %>% dplyr::select(ID, Grav_tot), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(GravNei = Grav_tot)

head(crypto5_BioRegion)

#filter connections within Biogeo regions
crypto5_BioRegion<-crypto5_BioRegion[crypto5_BioRegion$Source_Region == crypto5_BioRegion$Sink_Region,]
crypto5_BioRegion$Indegree <- ifelse(crypto5_BioRegion$Inflow > 0,1,0)
#
crypto5InflowBR <- aggregate(Inflow ~ To,crypto5_BioRegion,sum)
crypto5IDegreeBR <- aggregate(Indegree ~ To, crypto5_BioRegion,sum)
crypto5CorridorIndegreeBR <- aggregate(CorridorIN ~ To, crypto5_BioRegion,sum)
crypto5grav_neiBR <- aggregate(GravNei ~ To, crypto5_BioRegion,sum)

crypto5BR<-left_join(crypto5InflowBR,crypto5IDegreeBR, by="To")
crypto5BR<-left_join(crypto5BR,crypto5CorridorIndegreeBR, by="To")
crypto5BR<-left_join(crypto5BR,crypto5grav_neiBR, by="To")
head(crypto5BR)

cor(crypto5BR$Inflow,crypto5BR$Indegree) #0.78 correlation

#MPA Indegree/Inflow
crypto5_BioRegion$InflowMPA<-ifelse(crypto5_BioRegion$General == "MPA",crypto5_BioRegion$Inflow,0)
crypto5_BioRegion$IndegreeMPA<-ifelse(crypto5_BioRegion$General == "MPA",crypto5_BioRegion$Indegree,0)
crypto5InflowMPABR <- aggregate(InflowMPA ~ To,crypto5_BioRegion,sum)
crypto5IDegreeMPABR <- aggregate(IndegreeMPA ~ To, crypto5_BioRegion,sum)
crypto5BR<-left_join(crypto5BR,crypto5IDegreeMPABR, by="To")
crypto5BR<-left_join(crypto5BR,crypto5InflowMPABR, by="To")
head(crypto5BR)

#Average indegree of neighbours/
crypto5InflowNBBR <- aggregate(Source_Indegree ~ To,crypto5_BioRegion,mean)
crypto5IDegreeNBBR <- aggregate(Source_Inflow ~ To, crypto5_BioRegion,mean)
crypto5BR<-left_join(crypto5BR,crypto5InflowNBBR, by="To")
crypto5BR<-left_join(crypto5BR,crypto5IDegreeNBBR, by="To")
dim(crypto5BR) #merge this data later with other ones
crypto5BR$ModelMode <- ifelse(crypto5BR$Inflow > 0, "crypto5","")
head(crypto5BR)

#--------- Crypto 15
#Connections crypto
myconnections<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Connectivity/crypto15.csv")
head(myconnections)
dataconec<-myconnections[,2:4]
colnames(dataconec) <- c("From","To","Inflow")
head(dataconec)

#Kulbicki
crypto15_BioRegion <- dataconec %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('From' = 'ID'), all=F) %>%
  dplyr::rename(Source_Region = Kulbicki) %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('To' = 'ID'), all=F) %>%
  dplyr::rename(Sink_Region = Kulbicki) %>% 
  inner_join(nodesFinal %>% dplyr::select(ID, General), by = c('From' = 'ID'), all=F) %>%
  inner_join(global_factors %>% dplyr::select(ID, crypto15in), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Indegree = crypto15in) %>%
  inner_join(global_factors %>% dplyr::select(ID, crypto15IF), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Inflow = crypto15IF) %>% 
  inner_join(global_factors %>% dplyr::select(ID, corridorcrypto15), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(CorridorIN = corridorcrypto15) %>% 
  inner_join(global_factors %>% dplyr::select(ID, Grav_tot), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(GravNei = Grav_tot)

head(crypto15_BioRegion)

#filter connections within Biogeo regions
crypto15_BioRegion<-crypto15_BioRegion[crypto15_BioRegion$Source_Region == crypto15_BioRegion$Sink_Region,]
crypto15_BioRegion$Indegree <- ifelse(crypto15_BioRegion$Inflow > 0,1,0)
#
crypto15InflowBR <- aggregate(Inflow ~ To,crypto15_BioRegion,sum)
crypto15IDegreeBR <- aggregate(Indegree ~ To, crypto15_BioRegion,sum)
crypto15CorridorIndegreeBR <- aggregate(CorridorIN ~ To, crypto15_BioRegion,sum)
crypto15grav_neiBR <- aggregate(GravNei ~ To, crypto15_BioRegion,sum)

crypto15BR<-left_join(crypto15InflowBR,crypto15IDegreeBR, by="To")
crypto15BR<-left_join(crypto15BR,crypto15CorridorIndegreeBR, by="To")
crypto15BR<-left_join(crypto15BR,crypto15grav_neiBR, by="To")
head(crypto15BR)

cor(crypto15BR$Inflow,crypto15BR$Indegree) #0.78 correlation

#MPA Indegree/Inflow
crypto15_BioRegion$InflowMPA<-ifelse(crypto15_BioRegion$General == "MPA",crypto15_BioRegion$Inflow,0)
crypto15_BioRegion$IndegreeMPA<-ifelse(crypto15_BioRegion$General == "MPA",crypto15_BioRegion$Indegree,0)
crypto15InflowMPABR <- aggregate(InflowMPA ~ To,crypto15_BioRegion,sum)
crypto15IDegreeMPABR <- aggregate(IndegreeMPA ~ To, crypto15_BioRegion,sum)
crypto15BR<-left_join(crypto15BR,crypto15IDegreeMPABR, by="To")
crypto15BR<-left_join(crypto15BR,crypto15InflowMPABR, by="To")
head(crypto15BR)

#Average indegree of neighbours/
crypto15InflowNBBR <- aggregate(Source_Indegree ~ To,crypto15_BioRegion,mean)
crypto15IDegreeNBBR <- aggregate(Source_Inflow ~ To, crypto15_BioRegion,mean)
crypto15BR<-left_join(crypto15BR,crypto15InflowNBBR, by="To")
crypto15BR<-left_join(crypto15BR,crypto15IDegreeNBBR, by="To")
dim(crypto15BR) #merge this data later with other ones
crypto15BR$ModelMode <- ifelse(crypto15BR$Inflow > 0, "crypto15","")
head(crypto15BR)

#--------- Pare
#Connections pare
myconnections<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Connectivity/pare.csv")
head(myconnections)
dataconec<-myconnections[,2:4]
colnames(dataconec) <- c("From","To","Inflow")
head(dataconec)

#Kulbicki
pare_BioRegion <- dataconec %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('From' = 'ID'), all=F) %>%
  dplyr::rename(Source_Region = Kulbicki) %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('To' = 'ID'), all=F) %>%
  dplyr::rename(Sink_Region = Kulbicki) %>% 
  inner_join(nodesFinal %>% dplyr::select(ID, General), by = c('From' = 'ID'), all=F) %>%
  inner_join(global_factors %>% dplyr::select(ID, parein), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Indegree = parein) %>%
  inner_join(global_factors %>% dplyr::select(ID, pareIF), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Inflow = pareIF) %>% 
  inner_join(global_factors %>% dplyr::select(ID, corridorpare), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(CorridorIN = corridorpare) %>% 
  inner_join(global_factors %>% dplyr::select(ID, Grav_tot), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(GravNei = Grav_tot)

head(pare_BioRegion)

#filter connections within Biogeo regions
pare_BioRegion<-pare_BioRegion[pare_BioRegion$Source_Region == pare_BioRegion$Sink_Region,]
pare_BioRegion$Indegree <- ifelse(pare_BioRegion$Inflow > 0,1,0)
#
pareInflowBR <- aggregate(Inflow ~ To,pare_BioRegion,sum)
pareIDegreeBR <- aggregate(Indegree ~ To, pare_BioRegion,sum)
pareCorridorIndegreeBR <- aggregate(CorridorIN ~ To, pare_BioRegion,sum)
paregrav_neiBR <- aggregate(GravNei ~ To, pare_BioRegion,sum)

pareBR<-left_join(pareInflowBR,pareIDegreeBR, by="To")
pareBR<-left_join(pareBR,pareCorridorIndegreeBR, by="To")
pareBR<-left_join(pareBR,paregrav_neiBR, by="To")
head(pareBR)

cor(pareBR$Inflow,pareBR$Indegree) #0.78 correlation

#MPA Indegree/Inflow
pare_BioRegion$InflowMPA<-ifelse(pare_BioRegion$General == "MPA",pare_BioRegion$Inflow,0)
pare_BioRegion$IndegreeMPA<-ifelse(pare_BioRegion$General == "MPA",pare_BioRegion$Indegree,0)
pareInflowMPABR <- aggregate(InflowMPA ~ To,pare_BioRegion,sum)
pareIDegreeMPABR <- aggregate(IndegreeMPA ~ To, pare_BioRegion,sum)
pareBR<-left_join(pareBR,pareIDegreeMPABR, by="To")
pareBR<-left_join(pareBR,pareInflowMPABR, by="To")
head(pareBR)

#Average indegree of neighbours/
pareInflowNBBR <- aggregate(Source_Indegree ~ To,pare_BioRegion,mean)
pareIDegreeNBBR <- aggregate(Source_Inflow ~ To, pare_BioRegion,mean)
pareBR<-left_join(pareBR,pareInflowNBBR, by="To")
pareBR<-left_join(pareBR,pareIDegreeNBBR, by="To")
dim(pareBR) #merge this data later with other ones
pareBR$ModelMode <- ifelse(pareBR$Inflow > 0, "pare","")
head(pareBR)

#--------- Pare 5
myconnections<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Connectivity/pare5.csv")
head(myconnections)
dataconec<-myconnections[,2:4]
colnames(dataconec) <- c("From","To","Inflow")
head(dataconec)

#Kulbicki
pare5_BioRegion <- dataconec %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('From' = 'ID'), all=F) %>%
  dplyr::rename(Source_Region = Kulbicki) %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('To' = 'ID'), all=F) %>%
  dplyr::rename(Sink_Region = Kulbicki) %>% 
  inner_join(nodesFinal %>% dplyr::select(ID, General), by = c('From' = 'ID'), all=F) %>%
  inner_join(global_factors %>% dplyr::select(ID, pare5in), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Indegree = pare5in) %>%
  inner_join(global_factors %>% dplyr::select(ID, pare5IF), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Inflow = pare5IF) %>% 
  inner_join(global_factors %>% dplyr::select(ID, corridorpare5), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(CorridorIN = corridorpare5) %>% 
  inner_join(global_factors %>% dplyr::select(ID, Grav_tot), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(GravNei = Grav_tot)

head(pare5_BioRegion)

#filter connections within Biogeo regions
pare5_BioRegion<-pare5_BioRegion[pare5_BioRegion$Source_Region == pare5_BioRegion$Sink_Region,]
pare5_BioRegion$Indegree <- ifelse(pare5_BioRegion$Inflow > 0,1,0)
#
pare5InflowBR <- aggregate(Inflow ~ To,pare5_BioRegion,sum)
pare5IDegreeBR <- aggregate(Indegree ~ To, pare5_BioRegion,sum)
pare5CorridorIndegreeBR <- aggregate(CorridorIN ~ To, pare5_BioRegion,sum)
pare5grav_neiBR <- aggregate(GravNei ~ To, pare5_BioRegion,sum)

pare5BR<-left_join(pare5InflowBR,pare5IDegreeBR, by="To")
pare5BR<-left_join(pare5BR,pare5CorridorIndegreeBR, by="To")
pare5BR<-left_join(pare5BR,pare5grav_neiBR, by="To")
head(pare5BR)

cor(pare5BR$Inflow,pare5BR$Indegree) #0.78 correlation

#MPA Indegree/Inflow
pare5_BioRegion$InflowMPA<-ifelse(pare5_BioRegion$General == "MPA",pare5_BioRegion$Inflow,0)
pare5_BioRegion$IndegreeMPA<-ifelse(pare5_BioRegion$General == "MPA",pare5_BioRegion$Indegree,0)
pare5InflowMPABR <- aggregate(InflowMPA ~ To,pare5_BioRegion,sum)
pare5IDegreeMPABR <- aggregate(IndegreeMPA ~ To, pare5_BioRegion,sum)
pare5BR<-left_join(pare5BR,pare5IDegreeMPABR, by="To")
pare5BR<-left_join(pare5BR,pare5InflowMPABR, by="To")
head(pare5BR)

#Average indegree of neighbours/
pare5InflowNBBR <- aggregate(Source_Indegree ~ To,pare5_BioRegion,mean)
pare5IDegreeNBBR <- aggregate(Source_Inflow ~ To, pare5_BioRegion,mean)
pare5BR<-left_join(pare5BR,pare5InflowNBBR, by="To")
pare5BR<-left_join(pare5BR,pare5IDegreeNBBR, by="To")
dim(pare5BR) #merge this data later with other ones
pare5BR$ModelMode <- ifelse(pare5BR$Inflow > 0, "pare5","")
head(pare5BR)

#--------- Pare 15
myconnections<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Connectivity/pare15.csv")
head(myconnections)
dataconec<-myconnections[,2:4]
colnames(dataconec) <- c("From","To","Inflow")
head(dataconec)

#Kulbicki
pare15_BioRegion <- dataconec %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('From' = 'ID'), all=F) %>%
  dplyr::rename(Source_Region = Kulbicki) %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('To' = 'ID'), all=F) %>%
  dplyr::rename(Sink_Region = Kulbicki) %>% 
  inner_join(nodesFinal %>% dplyr::select(ID, General), by = c('From' = 'ID'), all=F) %>%
  inner_join(global_factors %>% dplyr::select(ID, pare15in), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Indegree = pare15in) %>%
  inner_join(global_factors %>% dplyr::select(ID, pare15IF), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Inflow = pare15IF) %>% 
  inner_join(global_factors %>% dplyr::select(ID, corridorpare15), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(CorridorIN = corridorpare15) %>% 
  inner_join(global_factors %>% dplyr::select(ID, Grav_tot), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(GravNei = Grav_tot)

head(pare15_BioRegion)

#filter connections within Biogeo regions
pare15_BioRegion<-pare15_BioRegion[pare15_BioRegion$Source_Region == pare15_BioRegion$Sink_Region,]
pare15_BioRegion$Indegree <- ifelse(pare15_BioRegion$Inflow > 0,1,0)
#
pare15InflowBR <- aggregate(Inflow ~ To,pare15_BioRegion,sum)
pare15IDegreeBR <- aggregate(Indegree ~ To, pare15_BioRegion,sum)
pare15CorridorIndegreeBR <- aggregate(CorridorIN ~ To, pare15_BioRegion,sum)
pare15grav_neiBR <- aggregate(GravNei ~ To, pare15_BioRegion,sum)

pare15BR<-left_join(pare15InflowBR,pare15IDegreeBR, by="To")
pare15BR<-left_join(pare15BR,pare15CorridorIndegreeBR, by="To")
pare15BR<-left_join(pare15BR,pare15grav_neiBR, by="To")
head(pare15BR)

cor(pare15BR$Inflow,pare15BR$Indegree) #0.78 correlation

#MPA Indegree/Inflow
pare15_BioRegion$InflowMPA<-ifelse(pare15_BioRegion$General == "MPA",pare15_BioRegion$Inflow,0)
pare15_BioRegion$IndegreeMPA<-ifelse(pare15_BioRegion$General == "MPA",pare15_BioRegion$Indegree,0)
pare15InflowMPABR <- aggregate(InflowMPA ~ To,pare15_BioRegion,sum)
pare15IDegreeMPABR <- aggregate(IndegreeMPA ~ To, pare15_BioRegion,sum)
pare15BR<-left_join(pare15BR,pare15IDegreeMPABR, by="To")
pare15BR<-left_join(pare15BR,pare15InflowMPABR, by="To")
head(pare15BR)

#Average indegree of neighbours/
pare15InflowNBBR <- aggregate(Source_Indegree ~ To,pare15_BioRegion,mean)
pare15IDegreeNBBR <- aggregate(Source_Inflow ~ To, pare15_BioRegion,mean)
pare15BR<-left_join(pare15BR,pare15InflowNBBR, by="To")
pare15BR<-left_join(pare15BR,pare15IDegreeNBBR, by="To")
dim(pare15BR) #merge this data later with other ones
pare15BR$ModelMode <- ifelse(pare15BR$Inflow > 0, "pare15","")
head(pare15BR)
#--------- Resid
myconnections<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Connectivity/resid.csv")
head(myconnections)
dataconec<-myconnections[,2:4]
colnames(dataconec) <- c("From","To","Inflow")
head(dataconec)

#Kulbicki
resid_BioRegion <- dataconec %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('From' = 'ID'), all=F) %>%
  dplyr::rename(Source_Region = Kulbicki) %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('To' = 'ID'), all=F) %>%
  dplyr::rename(Sink_Region = Kulbicki) %>% 
  inner_join(nodesFinal %>% dplyr::select(ID, General), by = c('From' = 'ID'), all=F) %>%
  inner_join(global_factors %>% dplyr::select(ID, residin), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Indegree = residin) %>%
  inner_join(global_factors %>% dplyr::select(ID, residIF), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Inflow = residIF) %>% 
  inner_join(global_factors %>% dplyr::select(ID, corridorresid), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(CorridorIN = corridorresid) %>% 
  inner_join(global_factors %>% dplyr::select(ID, Grav_tot), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(GravNei = Grav_tot)

head(resid_BioRegion)

#filter connections within Biogeo regions
resid_BioRegion<-resid_BioRegion[resid_BioRegion$Source_Region == resid_BioRegion$Sink_Region,]
resid_BioRegion$Indegree <- ifelse(resid_BioRegion$Inflow > 0,1,0)
#
residInflowBR <- aggregate(Inflow ~ To,resid_BioRegion,sum)
residIDegreeBR <- aggregate(Indegree ~ To, resid_BioRegion,sum)
residCorridorIndegreeBR <- aggregate(CorridorIN ~ To, resid_BioRegion,sum)
residgrav_neiBR <- aggregate(GravNei ~ To, resid_BioRegion,sum)

residBR<-left_join(residInflowBR,residIDegreeBR, by="To")
residBR<-left_join(residBR,residCorridorIndegreeBR, by="To")
residBR<-left_join(residBR,residgrav_neiBR, by="To")
head(residBR)

cor(residBR$Inflow,residBR$Indegree) #0.78 correlation

#MPA Indegree/Inflow
resid_BioRegion$InflowMPA<-ifelse(resid_BioRegion$General == "MPA",resid_BioRegion$Inflow,0)
resid_BioRegion$IndegreeMPA<-ifelse(resid_BioRegion$General == "MPA",resid_BioRegion$Indegree,0)
residInflowMPABR <- aggregate(InflowMPA ~ To,resid_BioRegion,sum)
residIDegreeMPABR <- aggregate(IndegreeMPA ~ To, resid_BioRegion,sum)
residBR<-left_join(residBR,residIDegreeMPABR, by="To")
residBR<-left_join(residBR,residInflowMPABR, by="To")
head(residBR)

#Average indegree of neighbours/
residInflowNBBR <- aggregate(Source_Indegree ~ To,resid_BioRegion,mean)
residIDegreeNBBR <- aggregate(Source_Inflow ~ To, resid_BioRegion,mean)
residBR<-left_join(residBR,residInflowNBBR, by="To")
residBR<-left_join(residBR,residIDegreeNBBR, by="To")
dim(residBR) #merge this data later with other ones
residBR$ModelMode <- ifelse(residBR$Inflow > 0, "resid","")
head(residBR)
#--------- Resid 15
myconnections<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Connectivity/resid15.csv")
head(myconnections)
dataconec<-myconnections[,2:4]
colnames(dataconec) <- c("From","To","Inflow")
head(dataconec)

#Kulbicki
resid15_BioRegion <- dataconec %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('From' = 'ID'), all=F) %>%
  dplyr::rename(Source_Region = Kulbicki) %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('To' = 'ID'), all=F) %>%
  dplyr::rename(Sink_Region = Kulbicki) %>% 
  inner_join(nodesFinal %>% dplyr::select(ID, General), by = c('From' = 'ID'), all=F) %>%
  inner_join(global_factors %>% dplyr::select(ID, resid15in), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Indegree = resid15in) %>%
  inner_join(global_factors %>% dplyr::select(ID, resid15IF), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Inflow = resid15IF) %>% 
  inner_join(global_factors %>% dplyr::select(ID, corridorresid15), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(CorridorIN = corridorresid15) %>% 
  inner_join(global_factors %>% dplyr::select(ID, Grav_tot), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(GravNei = Grav_tot)

head(resid15_BioRegion)

#filter connections within Biogeo regions
resid15_BioRegion<-resid15_BioRegion[resid15_BioRegion$Source_Region == resid15_BioRegion$Sink_Region,]
resid15_BioRegion$Indegree <- ifelse(resid15_BioRegion$Inflow > 0,1,0)
#
resid15InflowBR <- aggregate(Inflow ~ To,resid15_BioRegion,sum)
resid15IDegreeBR <- aggregate(Indegree ~ To, resid15_BioRegion,sum)
resid15CorridorIndegreeBR <- aggregate(CorridorIN ~ To, resid15_BioRegion,sum)
resid15grav_neiBR <- aggregate(GravNei ~ To, resid15_BioRegion,sum)

resid15BR<-left_join(resid15InflowBR,resid15IDegreeBR, by="To")
resid15BR<-left_join(resid15BR,resid15CorridorIndegreeBR, by="To")
resid15BR<-left_join(resid15BR,resid15grav_neiBR, by="To")
head(resid15BR)

cor(resid15BR$Inflow,resid15BR$Indegree) #0.78 correlation

#MPA Indegree/Inflow
resid15_BioRegion$InflowMPA<-ifelse(resid15_BioRegion$General == "MPA",resid15_BioRegion$Inflow,0)
resid15_BioRegion$IndegreeMPA<-ifelse(resid15_BioRegion$General == "MPA",resid15_BioRegion$Indegree,0)
resid15InflowMPABR <- aggregate(InflowMPA ~ To,resid15_BioRegion,sum)
resid15IDegreeMPABR <- aggregate(IndegreeMPA ~ To, resid15_BioRegion,sum)
resid15BR<-left_join(resid15BR,resid15IDegreeMPABR, by="To")
resid15BR<-left_join(resid15BR,resid15InflowMPABR, by="To")
head(resid15BR)

#Average indegree of neighbours/
resid15InflowNBBR <- aggregate(Source_Indegree ~ To,resid15_BioRegion,mean)
resid15IDegreeNBBR <- aggregate(Source_Inflow ~ To, resid15_BioRegion,mean)
resid15BR<-left_join(resid15BR,resid15InflowNBBR, by="To")
resid15BR<-left_join(resid15BR,resid15IDegreeNBBR, by="To")
dim(resid15BR) #merge this data later with other ones
resid15BR$ModelMode <- ifelse(resid15BR$Inflow > 0, "resid15","")
head(resid15BR)

#--------- Resid 35
myconnections<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Connectivity/resid35.csv")
head(myconnections)
dataconec<-myconnections[,2:4]
colnames(dataconec) <- c("From","To","Inflow")
head(dataconec)

#Kulbicki
resid35_BioRegion <- dataconec %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('From' = 'ID'), all=F) %>%
  dplyr::rename(Source_Region = Kulbicki) %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('To' = 'ID'), all=F) %>%
  dplyr::rename(Sink_Region = Kulbicki) %>% 
  inner_join(nodesFinal %>% dplyr::select(ID, General), by = c('From' = 'ID'), all=F) %>%
  inner_join(global_factors %>% dplyr::select(ID, resid35in), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Indegree = resid35in) %>%
  inner_join(global_factors %>% dplyr::select(ID, resid35IF), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Inflow = resid35IF) %>% 
  inner_join(global_factors %>% dplyr::select(ID, corridorresid35), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(CorridorIN = corridorresid35) %>% 
  inner_join(global_factors %>% dplyr::select(ID, Grav_tot), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(GravNei = Grav_tot)

head(resid35_BioRegion)

#filter connections within Biogeo regions
resid35_BioRegion<-resid35_BioRegion[resid35_BioRegion$Source_Region == resid35_BioRegion$Sink_Region,]
resid35_BioRegion$Indegree <- ifelse(resid35_BioRegion$Inflow > 0,1,0)
#
resid35InflowBR <- aggregate(Inflow ~ To,resid35_BioRegion,sum)
resid35IDegreeBR <- aggregate(Indegree ~ To, resid35_BioRegion,sum)
resid35CorridorIndegreeBR <- aggregate(CorridorIN ~ To, resid35_BioRegion,sum)
resid35grav_neiBR <- aggregate(GravNei ~ To, resid35_BioRegion,sum)

resid35BR<-left_join(resid35InflowBR,resid35IDegreeBR, by="To")
resid35BR<-left_join(resid35BR,resid35CorridorIndegreeBR, by="To")
resid35BR<-left_join(resid35BR,resid35grav_neiBR, by="To")
head(resid35BR)

cor(resid35BR$Inflow,resid35BR$Indegree) #0.78 correlation

#MPA Indegree/Inflow
resid35_BioRegion$InflowMPA<-ifelse(resid35_BioRegion$General == "MPA",resid35_BioRegion$Inflow,0)
resid35_BioRegion$IndegreeMPA<-ifelse(resid35_BioRegion$General == "MPA",resid35_BioRegion$Indegree,0)
resid35InflowMPABR <- aggregate(InflowMPA ~ To,resid35_BioRegion,sum)
resid35IDegreeMPABR <- aggregate(IndegreeMPA ~ To, resid35_BioRegion,sum)
resid35BR<-left_join(resid35BR,resid35IDegreeMPABR, by="To")
resid35BR<-left_join(resid35BR,resid35InflowMPABR, by="To")
head(resid35BR)

#Average indegree of neighbours/
resid35InflowNBBR <- aggregate(Source_Indegree ~ To,resid35_BioRegion,mean)
resid35IDegreeNBBR <- aggregate(Source_Inflow ~ To, resid35_BioRegion,mean)
resid35BR<-left_join(resid35BR,resid35InflowNBBR, by="To")
resid35BR<-left_join(resid35BR,resid35IDegreeNBBR, by="To")
dim(resid35BR) #merge this data later with other ones
resid35BR$ModelMode <- ifelse(resid35BR$Inflow > 0, "resid35","")
head(resid35BR)

#--------- Transi
myconnections<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Connectivity/transi.csv")
head(myconnections)
dataconec<-myconnections[,2:4]
colnames(dataconec) <- c("From","To","Inflow")
head(dataconec)

#Kulbicki
transi_BioRegion <- dataconec %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('From' = 'ID'), all=F) %>%
  dplyr::rename(Source_Region = Kulbicki) %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('To' = 'ID'), all=F) %>%
  dplyr::rename(Sink_Region = Kulbicki) %>% 
  inner_join(nodesFinal %>% dplyr::select(ID, General), by = c('From' = 'ID'), all=F) %>%
  inner_join(global_factors %>% dplyr::select(ID, transiin), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Indegree = transiin) %>%
  inner_join(global_factors %>% dplyr::select(ID, transiIF), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Inflow = transiIF) %>% 
  inner_join(global_factors %>% dplyr::select(ID, corridortransi), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(CorridorIN = corridortransi) %>% 
  inner_join(global_factors %>% dplyr::select(ID, Grav_tot), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(GravNei = Grav_tot)

head(transi_BioRegion)

#filter connections within Biogeo regions
transi_BioRegion<-transi_BioRegion[transi_BioRegion$Source_Region == transi_BioRegion$Sink_Region,]
transi_BioRegion$Indegree <- ifelse(transi_BioRegion$Inflow > 0,1,0)
#
transiInflowBR <- aggregate(Inflow ~ To,transi_BioRegion,sum)
transiIDegreeBR <- aggregate(Indegree ~ To, transi_BioRegion,sum)
transiCorridorIndegreeBR <- aggregate(CorridorIN ~ To, transi_BioRegion,sum)
transigrav_neiBR <- aggregate(GravNei ~ To, transi_BioRegion,sum)

transiBR<-left_join(transiInflowBR,transiIDegreeBR, by="To")
transiBR<-left_join(transiBR,transiCorridorIndegreeBR, by="To")
transiBR<-left_join(transiBR,transigrav_neiBR, by="To")
head(transiBR)

cor(transiBR$Inflow,transiBR$Indegree) #0.78 correlation

#MPA Indegree/Inflow
transi_BioRegion$InflowMPA<-ifelse(transi_BioRegion$General == "MPA",transi_BioRegion$Inflow,0)
transi_BioRegion$IndegreeMPA<-ifelse(transi_BioRegion$General == "MPA",transi_BioRegion$Indegree,0)
transiInflowMPABR <- aggregate(InflowMPA ~ To,transi_BioRegion,sum)
transiIDegreeMPABR <- aggregate(IndegreeMPA ~ To, transi_BioRegion,sum)
transiBR<-left_join(transiBR,transiIDegreeMPABR, by="To")
transiBR<-left_join(transiBR,transiInflowMPABR, by="To")
head(transiBR)

#Average indegree of neighbours/
transiInflowNBBR <- aggregate(Source_Indegree ~ To,transi_BioRegion,mean)
transiIDegreeNBBR <- aggregate(Source_Inflow ~ To, transi_BioRegion,mean)
transiBR<-left_join(transiBR,transiInflowNBBR, by="To")
transiBR<-left_join(transiBR,transiIDegreeNBBR, by="To")
dim(transiBR) #merge this data later with other ones
transiBR$ModelMode <- ifelse(transiBR$Inflow > 0, "transi","")
head(transiBR)
#--------- transi 15
myconnections<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Connectivity/transi15.csv")
head(myconnections)
dataconec<-myconnections[,2:4]
colnames(dataconec) <- c("From","To","Inflow")
head(dataconec)

#Kulbicki
transi15_BioRegion <- dataconec %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('From' = 'ID'), all=F) %>%
  dplyr::rename(Source_Region = Kulbicki) %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('To' = 'ID'), all=F) %>%
  dplyr::rename(Sink_Region = Kulbicki) %>% 
  inner_join(nodesFinal %>% dplyr::select(ID, General), by = c('From' = 'ID'), all=F) %>%
  inner_join(global_factors %>% dplyr::select(ID, transi15in), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Indegree = transi15in) %>%
  inner_join(global_factors %>% dplyr::select(ID, transi15IF), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Inflow = transi15IF) %>% 
  inner_join(global_factors %>% dplyr::select(ID, corridortransi15), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(CorridorIN = corridortransi15) %>% 
  inner_join(global_factors %>% dplyr::select(ID, Grav_tot), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(GravNei = Grav_tot)

head(transi15_BioRegion)

#filter connections within Biogeo regions
transi15_BioRegion<-transi15_BioRegion[transi15_BioRegion$Source_Region == transi15_BioRegion$Sink_Region,]
transi15_BioRegion$Indegree <- ifelse(transi15_BioRegion$Inflow > 0,1,0)
#
transi15InflowBR <- aggregate(Inflow ~ To,transi15_BioRegion,sum)
transi15IDegreeBR <- aggregate(Indegree ~ To, transi15_BioRegion,sum)
transi15CorridorIndegreeBR <- aggregate(CorridorIN ~ To, transi15_BioRegion,sum)
transi15grav_neiBR <- aggregate(GravNei ~ To, transi15_BioRegion,sum)

transi15BR<-left_join(transi15InflowBR,transi15IDegreeBR, by="To")
transi15BR<-left_join(transi15BR,transi15CorridorIndegreeBR, by="To")
transi15BR<-left_join(transi15BR,transi15grav_neiBR, by="To")
head(transi15BR)

cor(transi15BR$Inflow,transi15BR$Indegree) #0.78 correlation

#MPA Indegree/Inflow
transi15_BioRegion$InflowMPA<-ifelse(transi15_BioRegion$General == "MPA",transi15_BioRegion$Inflow,0)
transi15_BioRegion$IndegreeMPA<-ifelse(transi15_BioRegion$General == "MPA",transi15_BioRegion$Indegree,0)
transi15InflowMPABR <- aggregate(InflowMPA ~ To,transi15_BioRegion,sum)
transi15IDegreeMPABR <- aggregate(IndegreeMPA ~ To, transi15_BioRegion,sum)
transi15BR<-left_join(transi15BR,transi15IDegreeMPABR, by="To")
transi15BR<-left_join(transi15BR,transi15InflowMPABR, by="To")
head(transi15BR)

#Average indegree of neighbours/
transi15InflowNBBR <- aggregate(Source_Indegree ~ To,transi15_BioRegion,mean)
transi15IDegreeNBBR <- aggregate(Source_Inflow ~ To, transi15_BioRegion,mean)
transi15BR<-left_join(transi15BR,transi15InflowNBBR, by="To")
transi15BR<-left_join(transi15BR,transi15IDegreeNBBR, by="To")
dim(transi15BR) #merge this data later with other ones
transi15BR$ModelMode <- ifelse(transi15BR$Inflow > 0, "transi15","")
head(transi15BR)

#--------- transi 35
myconnections<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Connectivity/transi35.csv")
head(myconnections)
dataconec<-myconnections[,2:4]
colnames(dataconec) <- c("From","To","Inflow")
head(dataconec)

#Kulbicki
transi35_BioRegion <- dataconec %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('From' = 'ID'), all=F) %>%
  dplyr::rename(Source_Region = Kulbicki) %>%
  inner_join(nodesFinal %>% dplyr::select(ID, Kulbicki), by = c('To' = 'ID'), all=F) %>%
  dplyr::rename(Sink_Region = Kulbicki) %>% 
  inner_join(nodesFinal %>% dplyr::select(ID, General), by = c('From' = 'ID'), all=F) %>%
  inner_join(global_factors %>% dplyr::select(ID, transi35in), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Indegree = transi35in) %>%
  inner_join(global_factors %>% dplyr::select(ID, transi35IF), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(Source_Inflow = transi35IF) %>% 
  inner_join(global_factors %>% dplyr::select(ID, corridortransi35), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(CorridorIN = corridortransi35) %>% 
  inner_join(global_factors %>% dplyr::select(ID, Grav_tot), by=c('From' = 'ID'),all=F) %>%
  dplyr::rename(GravNei = Grav_tot)

head(transi35_BioRegion)

#filter connections within Biogeo regions
transi35_BioRegion<-transi35_BioRegion[transi35_BioRegion$Source_Region == transi35_BioRegion$Sink_Region,]
transi35_BioRegion$Indegree <- ifelse(transi35_BioRegion$Inflow > 0,1,0)
#
transi35InflowBR <- aggregate(Inflow ~ To,transi35_BioRegion,sum)
transi35IDegreeBR <- aggregate(Indegree ~ To, transi35_BioRegion,sum)
transi35CorridorIndegreeBR <- aggregate(CorridorIN ~ To, transi35_BioRegion,sum)
transi35grav_neiBR <- aggregate(GravNei ~ To, transi35_BioRegion,sum)

transi35BR<-left_join(transi35InflowBR,transi35IDegreeBR, by="To")
transi35BR<-left_join(transi35BR,transi35CorridorIndegreeBR, by="To")
transi35BR<-left_join(transi35BR,transi35grav_neiBR, by="To")
head(transi35BR)

cor(transi35BR$Inflow,transi35BR$Indegree) #0.78 correlation

#MPA Indegree/Inflow
transi35_BioRegion$InflowMPA<-ifelse(transi35_BioRegion$General == "MPA",transi35_BioRegion$Inflow,0)
transi35_BioRegion$IndegreeMPA<-ifelse(transi35_BioRegion$General == "MPA",transi35_BioRegion$Indegree,0)
transi35InflowMPABR <- aggregate(InflowMPA ~ To,transi35_BioRegion,sum)
transi35IDegreeMPABR <- aggregate(IndegreeMPA ~ To, transi35_BioRegion,sum)
transi35BR<-left_join(transi35BR,transi35IDegreeMPABR, by="To")
transi35BR<-left_join(transi35BR,transi35InflowMPABR, by="To")
head(transi35BR)

#Average indegree of neighbours/
transi35InflowNBBR <- aggregate(Source_Indegree ~ To,transi35_BioRegion,mean)
transi35IDegreeNBBR <- aggregate(Source_Inflow ~ To, transi35_BioRegion,mean)
transi35BR<-left_join(transi35BR,transi35InflowNBBR, by="To")
transi35BR<-left_join(transi35BR,transi35IDegreeNBBR, by="To")
dim(transi35BR) #merge this data later with other ones
transi35BR$ModelMode <- ifelse(transi35BR$Inflow > 0, "transi35","")
head(transi35BR)
##ALl databases:
library(plyr)
bioregionsconnectivity = Reduce(function(...) merge(..., by=c("ModelMode","To"), all=T), list(cryptoBR,crypto5BR,crypto15BR,pareBR,pare5BR,
                                                                                              pare15BR, residBR,resid15BR,resid35BR,transiBR,
                                                                                              transi15BR, transi35BR))

bioregionsconnectivity$To <- as.factor(bioregionsconnectivity$To)
unique(bioregionsconnectivity$To) #250 sites
#sum(length(crypto5BR$To),length(cryptoBR$To),
#    length(pare5BR$To),length(pareBR$To),
#    length(resid15BR$To),length(residBR$To),
#    length(transi15BR$To),length(transiBR$To)) #1837

bioregionsconnectivity[is.na(bioregionsconnectivity)] <- 0
bioregionfactors<- bioregionsconnectivity[,c(1:10)]
head(bioregionfactors)

colnames(bioregionfactors) <- c("ModelMode","ID","InflowBR","IndegreeBR","CorridorIndegreeBR",
                                "grav_neiBR","IndegreeMPABR","InflowMPABR",
                                "IndegreeNeiBR","InflowNeiBR")
head(bioregionfactors)

#including all nodes 
fulldata<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/fulldatabaseUpdated2403.csv")
dim(fulldata) #previous dataset

biogeofull<-left_join(fulldata,bioregionfactors,by=c("ID","ModelMode"),all=T)
biogeofull[is.na(biogeofull)] <- 0

#checking some attributes
dim(biogeofull) #3276 rows OK
head(biogeofull)
unique(biogeofull$ID) #261 ID's Ok
unique (biogeofull$sites) #272 sites OK

write.csv(biogeofull, "FullDataMay2020_LF.csv")


