##Implementation GAPS

library(ggplot2)
library(ggridges)
theme_set(theme_minimal())


#DATA
#Preparing database
nodesprop<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Data_building_factors/Grav-Ecoregions.csv",h=T)
nodes_grav_correction <- read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Data_building_factors/gravGlobal.csv",h=T)
kulbickiregions<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Data_building_factors/ecoregions_kulbicki.csv",h=T)
colnames(kulbickiregions) <- c("N","PROVINCE", "Kulbicki")
nodesprop<-left_join(nodesprop,kulbickiregions[,2:3],by="PROVINCE")
nodesprop<-left_join(nodesprop[,-6],nodes_grav_correction[,c("ID","Grav_tot")],by="ID")
length(nodesprop$ID) #ok 14804 ID nodes 

#Level Protection
nodesIDMPA<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Data_building_factors/IDs_with_MPAs.csv",h=T)
nodesIDMPA<-nodesIDMPA[,c("ID","Lon","Lat","Class","General")] #general is Locally managed (MPA vs fished/no restrictions)
nodesFinal<-left_join(nodesprop,nodesIDMPA[,c("ID","Class","General")],by="ID")
head(nodesFinal) #node Id, coordinates, gravity, ecoregion, province, bioregion, class, and general
unique(as.factor(nodesFinal$ID)) #14804 levels (ok)

#bigmama<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Data_building_factors/global_metrics.csv",h=T)
#BRbigmama<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Data_building_factors/global_metrics_BRmiss.csv",h=T)
bigmama<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/Data_building_factors/global_metrics_BR_and_total.csv",h=T)
bigmama[is.na(bigmama)] <-0
unique(bigmama$ID) #ok 14804 levels

global_factors<-left_join(bigmama,nodesFinal,by="ID")
unique(as.factor(global_factors$ID)) #14804 levels ok
#MPA/locally managed (1) or lack of management (0)
global_factors$MPA<-ifelse(global_factors$General == "MPA",1,0)
global_factors$ID <- as.factor(global_factors$ID)
head(global_factors)

netflowData<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_Ongoing/Connectivity analyses/MJDatawithNetflow.csv",h=T)
head(netflowData)
netflowData<-left_join(netflowData,global_factors[,c("ID","Kulbicki")],by="ID")
dim(netflowData[netflowData$Kulbicki == "Central Pacific",])


#Plot 4 regions #Central indo_pacific, central pacific, Western indian, Western Atlantic (TEP and Easten Atlantic low number of points)
dataHist<-netflowData[netflowData$Kulbicki %in% c("Central Indo_Pacific","Central Pacific","Western Atlantic", "Western Indian"),]
head(dataHist)
dataProp<-dataHist[,c("ID","Kulbicki","General","transi15Netflow","crypto5BROF","crypto5BRIF")]
head(dataProp)
dataProp$All<- as.factor(ifelse(dataProp$General == "MPA", "yes", "yes"))

#general plot
ggplot(dataProp, aes(x = transi15Netflow, y = All)) +
  geom_density_ridges(aes(fill="navy"), jittered_points = TRUE, quantile_lines = TRUE, scale = 0.5, alpha = 0.5,
                      vline_color = "gray25", point_color="#958D8D", point_size = 0.3, point_alpha = 0.001,quantiles = c(0.10, 0.90),
                      position = position_raincloud(adjust_vlines = TRUE)) +
  stat_density_ridges(quantile_lines = TRUE, vline_color = "gray25", quantiles = c(0.10, 0.90), alpha = 0.2, scale=0.5) + 
  geom_density_ridges((aes(point_color=General)), jittered_points = TRUE, quantile_lines = FALSE, scale = 0.5, alpha = 0.5,
                      vline_color = "gray25", point_size = 0.3, point_alpha = 0.25,quantiles = c(0.10, 0.90),
                      position = position_raincloud(adjust_vlines = TRUE)) +
  stat_density_ridges(quantile_lines = FALSE, vline_color = "gray25", quantiles = c(0.10, 0.90), alpha = 0.2, scale=0.5) +
  scale_discrete_manual(aesthetics = "point_colour", values = c("#958D8D", "#18AF3B")) + theme(legend.position = "none")

#by region plot
ggplot(dataProp, aes(x = transi15Netflow, y = Kulbicki)) +
  geom_density_ridges(aes(fill="navy"), jittered_points = TRUE, quantile_lines = TRUE, scale = 0.5, alpha = 0.5,
                      vline_color = "gray25", point_color="#958D8D", point_size = 0.3, point_alpha = 0.001,quantiles = c(0.10, 0.90),
                      position = position_raincloud(adjust_vlines = TRUE)) +
  stat_density_ridges(quantile_lines = TRUE, vline_color = "gray25", quantiles = c(0.10, 0.90), alpha = 0.2, scale=0.5) + 
  geom_density_ridges((aes(point_color=General)), jittered_points = TRUE, quantile_lines = FALSE, scale = 0.5, alpha = 0.5,
                      vline_color = "gray25", point_size = 0.3, point_alpha = 0.25,quantiles = c(0.10, 0.90),
                      position = position_raincloud(adjust_vlines = TRUE)) +
  stat_density_ridges(quantile_lines = FALSE, vline_color = "gray25", quantiles = c(0.10, 0.90), alpha = 0.2, scale=0.5) +
  scale_discrete_manual(aesthetics = "point_colour", values = c("#958D8D", "#18AF3B")) + theme(legend.position = "none")



dataProp<-na.omit(dataProp)
#14341 cells that either receive or seed larvae (NA's mean 0 Inflow and Outflow - highly isolated reef grids)

dim(dataProp[dataProp$Kulbicki == "Western Atlantic",]) #1385 (total number of grids) 
dim(dataProp[dataProp$Kulbicki == "Central Pacific",]) #4216 (total number of grids)
dim(dataProp[dataProp$Kulbicki == "Western Indian",]) #1830 (total number of grids)
dim(dataProp[dataProp$Kulbicki == "Central Indo_Pacific",]) #6910 (total number of grids)


#Figure 4 Top sinks and sources per region

#Top SINKS (general)

#total number of grids within MPAs
#dim(dataProp[dataProp$General == "MPA",]) #3483 reef cells 
#3483/14805 #23% grids within protected areas IUCN

dim(dataProp[dataProp$transi15Netflow < quantile(dataProp$transi15Netflow, 0.10) & dataProp$General == "MPA",])/
  dim(dataProp[dataProp$transi15Netflow < quantile(dataProp$transi15Netflow, 0.10),]) #399 MPAs / 1434 total top sinks
#27.8% top sinks protected

#by region
#Central Indo-Pacific
dim(dataProp[dataProp$transi15Netflow < quantile(dataProp[dataProp$Kulbicki == "Central Indo_Pacific",]$transi15Netflow, 0.10)
             & dataProp$General == "MPA" & 
               dataProp$Kulbicki == "Central Indo_Pacific",])/
  dim(dataProp[dataProp$transi15Netflow < quantile(dataProp[dataProp$Kulbicki == "Central Indo_Pacific",]$transi15Netflow, 0.10)
               & dataProp$Kulbicki == "Central Indo_Pacific",])
#7.9% of top sinks protected

##Central Pacific
dim(dataProp[dataProp$transi15Netflow < quantile(dataProp[dataProp$Kulbicki == "Central Pacific",]$transi15Netflow, 0.10) & dataProp$General == "MPA" & 
               dataProp$Kulbicki == "Central Pacific",])/
  dim(dataProp[dataProp$transi15Netflow < quantile(dataProp[dataProp$Kulbicki == "Central Pacific",]$transi15Netflow, 0.10) & dataProp$Kulbicki == "Central Pacific",]) 
#44% of top sinks protected

##Western Indian
dim(dataProp[dataProp$transi15Netflow < quantile(dataProp[dataProp$Kulbicki == "Western Indian",]$transi15Netflow, 0.10) & dataProp$General == "MPA" & 
               dataProp$Kulbicki == "Western Indian",])/
  dim(dataProp[dataProp$transi15Netflow < quantile(dataProp[dataProp$Kulbicki == "Western Indian",]$transi15Netflow, 0.10) & dataProp$Kulbicki == "Western Indian",])
#25.6% of top sinks protected

##Western Atlantic
dim(dataProp[dataProp$transi15Netflow < quantile(dataProp[dataProp$Kulbicki == "Western Atlantic",]$transi15Netflow, 0.10) & dataProp$General == "MPA" & 
               dataProp$Kulbicki == "Western Atlantic",])/
  dim(dataProp[dataProp$transi15Netflow < quantile(dataProp[dataProp$Kulbicki == "Western Atlantic",]$transi15Netflow, 0.10) & dataProp$Kulbicki == "Western Atlantic",])
#25.1% of top sources protected


#Top SOURCES
dim(dataProp[dataProp$transi15Netflow > quantile(dataProp$transi15Netflow, 0.90) & dataProp$General == "MPA",])/
  dim(dataProp[dataProp$transi15Netflow > quantile(dataProp$transi15Netflow, 0.90),])
#32% top sources protected
471/3483 #13.5 sources
399/3483 #11.4% sinks

#by region
#Central Indo-Pacific
dim(dataProp[dataProp$transi15Netflow > quantile(dataProp[dataProp$Kulbicki == "Central Indo_Pacific",]$transi15Netflow, 0.90)
             & dataProp$General == "MPA" & 
               dataProp$Kulbicki == "Central Indo_Pacific",])/
  dim(dataProp[dataProp$transi15Netflow > quantile(dataProp[dataProp$Kulbicki == "Central Indo_Pacific",]$transi15Netflow, 0.90)
               & dataProp$Kulbicki == "Central Indo_Pacific",])
#5.4% of top sources protected

##Central Pacific
dim(dataProp[dataProp$transi15Netflow > quantile(dataProp[dataProp$Kulbicki == "Central Pacific",]$transi15Netflow, 0.90) & dataProp$General == "MPA" & 
               dataProp$Kulbicki == "Central Pacific",])/
  dim(dataProp[dataProp$transi15Netflow > quantile(dataProp[dataProp$Kulbicki == "Central Pacific",]$transi15Netflow, 0.90) & dataProp$Kulbicki == "Central Pacific",]) 
#53% of top sources protected

##Western Indian
dim(dataProp[dataProp$transi15Netflow > quantile(dataProp[dataProp$Kulbicki == "Western Indian",]$transi15Netflow, 0.90) & dataProp$General == "MPA" & 
               dataProp$Kulbicki == "Western Indian",])/
  dim(dataProp[dataProp$transi15Netflow > quantile(dataProp[dataProp$Kulbicki == "Western Indian",]$transi15Netflow, 0.90) & dataProp$Kulbicki == "Western Indian",])
#12.5% of top sinks protected

##Western Atlantic
dim(dataProp[dataProp$transi15Netflow > quantile(dataProp[dataProp$Kulbicki == "Western Atlantic",]$transi15Netflow, 0.90) & dataProp$General == "MPA" & 
               dataProp$Kulbicki == "Western Atlantic",])/
  dim(dataProp[dataProp$transi15Netflow > quantile(dataProp[dataProp$Kulbicki == "Western Atlantic",]$transi15Netflow, 0.90) & dataProp$Kulbicki == "Western Atlantic",])
#17.9% of top sources protected



#Fig 3
#Ecological corridors (top 10%)

#general
ggplot(dataProp[dataProp$crypto5BROF>0 & dataProp$crypto5BRIF>0,], aes(x = log1p(crypto5BRIF), y =All)) +
  geom_density_ridges(aes(fill="navy"), jittered_points = TRUE, quantile_lines = TRUE, scale = 0.5, alpha = 0.5,
                      vline_color = "gray25", point_color="#958D8D", point_size = 0.3, point_alpha = 0.001,quantiles =0.90,
                      position = position_raincloud(adjust_vlines = TRUE)) +
  stat_density_ridges(quantile_lines = TRUE, vline_color = "gray25", quantiles =  0.90, alpha = 0.2, scale=0.5) + 
  geom_density_ridges((aes(point_color=General)), jittered_points = TRUE, quantile_lines = FALSE, scale = 0.5, alpha = 0.5,
                      vline_color = "gray25", point_size = 0.3, point_alpha = 0.25,quantiles = 0.90,
                      position = position_raincloud(adjust_vlines = TRUE)) +
  stat_density_ridges(quantile_lines = FALSE, vline_color = "gray25", quantiles = 0.90, alpha = 0.2, scale=0.5) +
  scale_discrete_manual(aesthetics = "point_colour", values = c("#958D8D", "#18AF3B")) + theme(legend.position = "none")


#by region
ggplot(dataProp[dataProp$crypto5BROF>0 & dataProp$crypto5BRIF>0,], aes(x = log1p(crypto5BRIF), y = Kulbicki)) +
  geom_density_ridges(aes(fill="navy"), jittered_points = TRUE, quantile_lines = TRUE, scale = 0.5, alpha = 0.5,
                      vline_color = "gray25", point_color="#958D8D", point_size = 0.3, point_alpha = 0.001,quantiles =0.90,
                      position = position_raincloud(adjust_vlines = TRUE)) +
  stat_density_ridges(quantile_lines = TRUE, vline_color = "gray25", quantiles =  0.90, alpha = 0.2, scale=0.5) + 
  geom_density_ridges((aes(point_color=General)), jittered_points = TRUE, quantile_lines = FALSE, scale = 0.5, alpha = 0.5,
                      vline_color = "gray25", point_size = 0.3, point_alpha = 0.25,quantiles = 0.90,
                      position = position_raincloud(adjust_vlines = TRUE)) +
  stat_density_ridges(quantile_lines = FALSE, vline_color = "gray25", quantiles = 0.90, alpha = 0.2, scale=0.5) +
  scale_discrete_manual(aesthetics = "point_colour", values = c("#958D8D", "#18AF3B")) + theme(legend.position = "none")



#Ecological corridors % protection
EcoCordata<-dataProp[dataProp$crypto5BROF>0 & dataProp$crypto5BRIF>0,] #IF and OF > 0 - means that a reefs is either receiving and seeding larvae

#general
dim(EcoCordata[EcoCordata$crypto5BRIF > quantile(EcoCordata$crypto5BRIF, 0.90) & EcoCordata$General == "MPA",])/
  dim(EcoCordata[EcoCordata$crypto5BRIF > quantile(EcoCordata$crypto5BRIF, 0.90),])
#29.7 corridors within MPAs

#by region
#Central IP
dim(EcoCordata[EcoCordata$crypto5BRIF > quantile(EcoCordata$crypto5BRIF, 0.90) & EcoCordata$General == "MPA" &
                 EcoCordata$Kulbicki == "Central Indo_Pacific",])/
  dim(EcoCordata[EcoCordata$crypto5BRIF > quantile(EcoCordata$crypto5BRIF, 0.90) & EcoCordata$Kulbicki == "Central Indo_Pacific",])
#8.1% MPAs

#Central Pacific
dim(EcoCordata[EcoCordata$crypto5BRIF > quantile(EcoCordata$crypto5BRIF, 0.90) & EcoCordata$General == "MPA" &
                 EcoCordata$Kulbicki == "Central Pacific",])/
  dim(EcoCordata[EcoCordata$crypto5BRIF > quantile(EcoCordata$crypto5BRIF, 0.90) & EcoCordata$Kulbicki == "Central Pacific",])
#77.9% MPAs

#Western Indian
dim(EcoCordata[EcoCordata$crypto5BRIF > quantile(EcoCordata$crypto5BRIF, 0.90) & EcoCordata$General == "MPA" &
                 EcoCordata$Kulbicki == "Western Indian",])/
  dim(EcoCordata[EcoCordata$crypto5BRIF > quantile(EcoCordata$crypto5BRIF, 0.90) & EcoCordata$Kulbicki == "Western Indian",])
#11.9%

#Western Atlantic
dim(EcoCordata[EcoCordata$crypto5BRIF > quantile(EcoCordata$crypto5BRIF, 0.90) & EcoCordata$General == "MPA" &
                 EcoCordata$Kulbicki == "Western Atlantic",])/
  dim(EcoCordata[EcoCordata$crypto5BRIF > quantile(EcoCordata$crypto5BRIF, 0.90) & EcoCordata$Kulbicki == "Western Atlantic",])
#22.1 



##HISTOGRAMS

#corridors
ggplot(dataHist[dataHist$crypto5BROF>0 & dataHist$crypto5BRIF>10,], aes(x=log1p(crypto5BRIF), color=General)) +
  geom_histogram(fill="white", alpha=0.5, position="identity") + theme_classic()

ggplot(dataHist[dataHist$crypto5BROF>0 & dataHist$crypto5BRIF>10,], aes(x=log1p(crypto5BRIF), color=General)) +
  geom_histogram(fill="white", alpha=0.5, position="identity") + facet_wrap(~Kulbicki)+ theme_classic()

median(log1p(dataHist[dataHist$crypto5BROF>0,]$crypto5BRIF))
max(log1p(dataHist[dataHist$crypto5BROF>0,]$crypto5BRIF))
min(log1p(dataHist[dataHist$crypto5BROF>0,]$crypto5BRIF))

median(dataHist[dataHist$crypto5BROF>0,]$crypto5BRIF)
min(dataHist[dataHist$crypto5BROF>0,]$crypto5BRIF)
max(dataHist[dataHist$crypto5BROF>0,]$crypto5BRIF)


#sinks and sources
ggplot(dataHist, aes(x=resid15Netflow, color=General)) +
  geom_histogram(fill="white", alpha=0.5, position="identity") + theme_classic()

ggplot(dataHist, aes(x=transi15Netflow, color=General)) +
  geom_histogram(fill="white", alpha=0.5, position="identity") + facet_wrap(~Kulbicki) + theme_classic()

ggplot(dataHist, aes(x=transi15Netflow, color=General)) +
  geom_histogram(fill="white", alpha=0.5, position="identity") + theme_classic()


