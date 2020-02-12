#Trying to overlay points
nodesID<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_on_going/Connectivity analyses/IDs.csv",h=T)
head(nodesID) #14804 id's
colnames(nodesID) <- c("ID","ID2","lon","lat","territory","other")
sites<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_on_going/Connectivity analyses/databiomassFull.csv",h=T)
###
#library(ggmap)
library(tidyverse)
library(RColorBrewer)
#theme_set(theme_bw())
library(maptools)
library(maps)
library(hutilscpp) 
View(filtertre)
#crossing nodes ID with sites 
dt<-match_nrst_haversine(sites$lat, sites$lon, nodesID$lat, nodesID$lon, nodesID$ID,
                         close_enough = 0.1)
sites$ID<-dt$pos
sites$distprox<-dt$dist
head(sites) #ok 
summary(bigmama)
library(igraph)

#checking<-merge(sites,nodesID,by=c("ID"),all=F)
dim(sites) #452 sites
#View(checking[,c("sites","locality","territory")])

##reading big mama - Global metrics for 14805 reef cells 
bigmama<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_on_going/Connectivity analyses/Connectivity/global_metrics.csv",h=T)

##
firstat<-merge(sites,bigmama,by=c("ID"),all=F)
dim(firstat) #ok
#distance from point (matching ID's)
filterone<-firstat[firstat$distprox<30,] #403 sites
filtertwo<-firstat[firstat$distprox<20,] #374 sites
filtertre<-firstat[firstat$distprox<10,] #258 sites

unique(filterone$sites) #404
unique(filtertwo$sites) #378
unique(filtertre$sites) #272

#creating matrix full - functional groups as factors
testandodata <- gather(filtertre, key = "Modeltype", "INdegree", c("parein","pare5in","pare15in","pare25in",
                                                                    "cryptoin","crypto5in","crypto15in","crypto35in",
                                                                    "residin","resid5in","resid15in","resid35in",
                                                                    "transiin","transi5in","transi15in","transi35in"))
testandodata$ModelMode<-gsub("in", " ", testandodata$Modeltype)
dim(testandodata)

testandodata2 <- gather(filtertre, key = "Modeltype", "Outdegree", c("pareout","pare5out","pare15out","pare25out",
                                                                    "cryptoout","crypto5out","crypto15out","crypto35out",
                                                                    "residout","resid5out","resid15out","resid35out",
                                                                   "transiout","transi5out","transi15out","transi35out"))
testandodata2$ModelMode<-gsub("out", " ", testandodata2$Modeltype)
head(testandodata2)

testandodata3 <- gather(filtertre, key = "Modeltype", "Btwdegree", c("parebtw","pare5btw","pare15btw","pare25btw",
                                                                      "cryptobtw","crypto5btw","crypto15btw","crypto35btw",
                                                                      "residbtw","resid5btw","resid15btw","resid35btw",
                                                                      "transibtw","transi5btw","transi15btw","transi35btw"))

testandodata3$ModelMode<-gsub("btw", " ", testandodata3$Modeltype)

testandodata4 <- gather(filtertre, key = "Modeltype", "LocalRet", c("pareLR","pare5LR","pare15LR","pare25LR",
                                                                      "cryptoLR","crypto5LR","crypto15LR","crypto35LR",
                                                                      "residLR","resid5LR","resid15LR","resid35LR",
                                                                      "transiLR","transi5LR","transi15LR","transi35LR"))
testandodata4$ModelMode<-gsub("LR", " ", testandodata4$Modeltype)
head(testandodata4)

unique(testandodata$ID)
##LR
dd<-(testandodata4[,c("ID","region","locality","sites","lat","lon","biomassarea","biomassarea1","biomassarea2","grav_total","ModelMode","LocalRet")])
#dd$ModelMode<-gsub("LR", " ", dd$Modeltype)
dim(dd) #ok
#in
ww<-(testandodata[,c("ID","region","locality","sites","lat","lon","biomassarea","biomassarea1","biomassarea2","grav_total","ModelMode","INdegree")])
#ww$ModelMode<-gsub("in", " ", ww$Modeltype)
dim(ww)
#out
cc<-(testandodata2[,c("ID","region","locality","sites","lat","lon","biomassarea","biomassarea1","biomassarea2","grav_total","ModelMode","Outdegree")])
#cc$ModelMode<-gsub("out", " ", cc$Modeltype)
dim(cc)
#btw
ee<-(testandodata3[,c("ID","region","locality","sites","lat","lon","biomassarea","biomassarea1","biomassarea2","grav_total","ModelMode","Btwdegree")])
dim(ee)
#ee$ModelMode<-gsub("btw", " ", ee$Modeltype)

#merge daaframes dd (LR), ww (in), cc (out), ee (btw)
identical(ee[,1:11],dd[,1:11]) #ok
identical(ww[,1:11],dd[,1:11]) #ok
identical(cc[,1:11],dd[,1:11]) #ok
#obs - whe I try to merge by columns, for some reason, it creates more rows. COnsidering that all 4 database are identical, its ok to paste additional columns.

dd$Indegree <- ww$INdegree
dd$Outdegree<-cc$Outdegree
dd$Btw<-ee$Btwdegree

mainmatrix<-dd
summary(mainmatrix)

#write.csv(mainmatrix, "mainmatrix.csv")

unique(mainmatrix$sites) #272 sites
unique(mainmatrix$ID) #261 ID (cells)

#plot some maps----

#Active transien IN degree
library(mapproj)
mapcor<-merge(filtertre,nodesID, by="ID")
dim(mapcor) #ok

mapWorld <- borders("world", colour="gray", fill="gray") # create a layer of borders

mp <- ggplot() +   mapWorld + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank()) +
  coord_map(ylim=c(-50, 50))

#Global patterns
Biomassmap <- mp + geom_point(aes(x=filtertre$lon, y=filtertre$lat, color=filtertre$biomassarea2), 
                       alpha=0.5, size=3) +
  scale_color_gradient(low = "blue",high = "red",name="Biomass (kg)/area (m2)") +
  labs(x="",y="")
Biomassmap
