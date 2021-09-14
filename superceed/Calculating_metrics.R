#Trying to overlay points
nodesID<-read.csv("~/OneDrive - Macquarie University/Chapters/Chapter_04_on_going/Connectivity data/IDs.csv",h=T)
dim(nodesID) #14804 id's
colnames(nodesID) <- c("ID2","ID","lon","lat","territory","sites")
sites<-read.csv("~/OneDrive - Macquarie University/Chapters/Chapter_04_on_going/Connectivity data/biomass_data.csv",h=T)
View(sites)
#####-------- plotting MAP
library(ggmap)
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())
library(maptools)
library(maps)
mp <- NULL

#mapWorld <- borders("world", colour="gray", fill="gray") # create a layer of borders
#mp <- ggplot() +   mapWorld + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                    panel.background = element_blank()) +  coord_map(ylim=c(-50, 50))

#mp3 <- mp + geom_point(aes(x=nodesID[nodesID$ID2 == "I7566",]$lon, y=nodesID[nodesID$ID2 == "I7566",]$lat), 
#                       , size=2, color="blue") + labs(x="",y="") +
#  geom_point(aes(x=sites[sites$sites == "cocos_1",]$lon, y=sites[sites$sites == "cocos_1",]$lat), 
#              size=2, color="red")

library(hutilscpp) #crossing nodes ID with sites 
dt<-match_nrst_haversine(sites$lat, sites$lon, nodesID$lat, nodesID$lon, nodesID$ID2,
                 close_enough = 0.1)
sites$nodeID<-dt$pos
sites$distprox<-dt$dist
head(sites) #ok 
View(nodesID)



#1st group = Cryptobenthic group (passive, act 5, act 15 and act 30)
cryptopas<-read.csv("/Volumes/LuisaDrive/Model_outputs/Cryptobenthic/Passive/con_file_pre30.csv",h=FALSE)
cryptoact5<-read.csv("/Volumes/LuisaDrive/Model_outputs/Cryptobenthic/Active/con_file_pre5.csv",h=FALSE)
cryptoact15<-read.csv("/Volumes/LuisaDrive/Model_outputs/Cryptobenthic/Active/con_file_pre15.csv",h=FALSE)
cryptoact30<-read.csv("/Volumes/LuisaDrive/Model_outputs/Cryptobenthic/Active/con_file_pre30.csv",h=FALSE)

#need to replane NaN -> 0
crypto <- replace(cryptopas, cryptopas == "NaN", 0)
rownames(crypto) <-colnames(crypto)
dim(crypto) #dim ok
crypto5 <- replace(cryptoact5, cryptoact5 == "NaN", 0)
rownames(crypto5) <-colnames(crypto5)
crypto15 <- replace(cryptoact15, cryptoact15 == "NaN", 0)
rownames(crypto15) <-colnames(crypto15)
crypto35 <- replace(cryptoact30, cryptoact30 == "NaN", 0)
rownames(crypto35) <-colnames(crypto35)

#Igraph - metrics
library(igraph)

cryptoi<-graph.adjacency(as.matrix(crypto),diag=T,mode="directed",weighted = TRUE)
crypto5i<-graph.adjacency(as.matrix(crypto5),diag=T,mode="directed",weighted = TRUE)
crypto15i<-graph.adjacency(as.matrix(crypto15),diag=T,mode="directed",weighted = TRUE)
crypto35i<-graph.adjacency(as.matrix(crypto35),diag=T,mode="directed",weighted = TRUE)

#In degree

cryptoin <-degree(cryptoi, mode="in")
cryptoinN <-degree(cryptoi, mode="in", normalized = T)

crypto5in <-degree(crypto5i, mode="in")
crypto5inN <-degree(crypto5i, mode="in", normalized = T)

crypto15in <-degree(crypto15i, mode="in")
crypto15inN <-degree(crypto15i, mode="in", normalized = T)

crypto35in <-degree(crypto35i, mode="in")
crypto35inN <-degree(crypto35i, mode="in", normalized = T)

#Out degree

cryptoout <-degree(cryptoi, mode="out")
cryptooutN <-degree(cryptoi, mode="out", normalized = T)

crypto5out <-degree(crypto5i, mode="out")
crypto5outN <-degree(crypto5i, mode="out", normalized = T)

crypto15out <-degree(crypto15i, mode="out")
crypto15outN <-degree(crypto15i, mode="out", normalized = T)

crypto35out <-degree(crypto35i, mode="out")
crypto35outN <-degree(crypto35i, mode="out", normalized = T)

#Between
cryptobtw <-betweenness(cryptoi)
#cryptobtw <-betweenness(cryptoi, normalized = T)

crypto5btw <-betweenness(crypto5i)
#crypto5btwN <-betweenness(crypto5i, normalized = T)

crypto15btw <-betweenness(crypto15i)
#crypto15btwN <-betweenness(crypto15i, normalized = T)

crypto35btw <-betweenness(crypto35i)
#crypto35btwN <-betweenness(crypto35i, normalized = T)


#Diagonal - Local retention
cryptoLR <-diag(as.matrix(crypto))
crypto5LR <-diag(as.matrix(crypto5))
crypto15LR <-diag(as.matrix(crypto15))
crypto35LR <-diag(as.matrix(crypto35))

#ok; for active 5,15 and 30 = out, in, btw and lr ready
#1.3 One object with nodes as column 
cryptometrics<-as.data.frame(colnames(crypto5))
colnames(cryptometrics) <- c("ID")
dim(cryptometrics)

#IN
cryptometrics$cryptoin <- cryptoin
cryptometrics$crypto5in <- crypto5in
cryptometrics$crypto15in <- crypto15in
cryptometrics$crypto35in <- crypto35in

cryptometrics$cryptoinN <- cryptoinN
cryptometrics$crypto5inN <- crypto5inN
cryptometrics$crypto15inN<- crypto15inN
cryptometrics$crypto35inN <- crypto35inN

#OUT
cryptometrics$cryptoout <- cryptoout
cryptometrics$crypto5out <- crypto5out
cryptometrics$crypto15out <- crypto15out
cryptometrics$crypto35out <- crypto35out

cryptometrics$cryptooutN <- cryptooutN
cryptometrics$crypto5outN <- crypto5outN
cryptometrics$crypto15outN<- crypto15outN
cryptometrics$crypto35outN <- crypto35outN

#BTW
cryptometrics$cryptobtw <- cryptobtw
cryptometrics$crypto5btw <- crypto5btw
cryptometrics$crypto15btw <- crypto15btw
cryptometrics$crypto35btw <- crypto35btw

#LR
cryptometrics$cryptoLR <- cryptoLR
cryptometrics$crypto5LR <- crypto5LR
cryptometrics$crypto15LR <- crypto15LR
cryptometrics$crypto35LR <- crypto35LR

dim(cryptometrics) #ok, only missing for Passive model yohoo

#2nd group = Parental care group  - (passive csv is missing - empty folder)
parepas<-read.csv("/Volumes/LuisaDrive/Model_outputs/Parental/Passive/con_file_pre25.csv",h=FALSE)
parenact5<-read.csv("/Volumes/LuisaDrive/Model_outputs/Parental/Active/con_file_pre5.csv",h=FALSE)
pareact15<-read.csv("/Volumes/LuisaDrive/Model_outputs/Parental/Active/con_file_pre15.csv",h=FALSE)
pareact25<-read.csv("/Volumes/LuisaDrive/Model_outputs/Parental/Active/con_file_pre25.csv",h=FALSE)


#need to replane NaN -> 0
pare <- replace(parepas, parepas == "NaN", 0)
rownames(pare) <-colnames(pare)
dim(pare) #dim ok
pare5 <- replace(parenact5, parenact5 == "NaN", 0)
rownames(pare5) <-colnames(pare5)
dim(pare5) #ok
pare15 <- replace(pareact15, pareact15 == "NaN", 0)
rownames(pare15) <-colnames(pare15)
dim(pare15) #ok
pare15 <- replace(pareact15, pareact15 == "NaN", 0)
rownames(pare15) <-colnames(pare15)
dim(pare15) #ok
pare25 <- replace(pareact25, pareact25 == "NaN", 0)
rownames(pare25) <-colnames(pare25)
dim(pare25) #ok


#Igraph - metrics

#missing the passive
parei<-graph.adjacency(as.matrix(pare),diag=T,mode="directed",weighted = TRUE)
pare5i<-graph.adjacency(as.matrix(pare5),diag=T,mode="directed",weighted = TRUE)
pare15i<-graph.adjacency(as.matrix(pare15),diag=T,mode="directed",weighted = TRUE)
pare25i<-graph.adjacency(as.matrix(pare25),diag=T,mode="directed",weighted = TRUE)

#In degree
parein <-degree(parei, mode="in")
pareinN <-degree(parei, mode="in", normalized = T)

pare5in <-degree(pare5i, mode="in")
pare5inN <-degree(pare5i, mode="in", normalized = T)

pare15in <-degree(pare15i, mode="in")
pare15inN <-degree(pare15i, mode="in", normalized = T)

pare25in <-degree(pare25i, mode="in")
pare25inN <-degree(pare25i, mode="in", normalized = T)

#Out degree
pareout <-degree(parei, mode="out")
pareoutN <-degree(parei, mode="out", normalized = T)

pare5out <-degree(pare5i, mode="out")
pare5outN <-degree(pare5i, mode="out", normalized = T)

pare15out <-degree(pare15i, mode="out")
pare15outN <-degree(pare15i, mode="out", normalized = T)

pare25out <-degree(pare25i, mode="out")
pare25outN <-degree(pare25i, mode="out", normalized = T)

#Between
parebtw <-betweenness(parei)

pare5btw <-betweenness(pare5i)
#pare5btwN <-betweenness(pare5i, normalized = T)

pare15btw <-betweenness(pare15i)
#pare15btwN <-betweenness(pare15i, normalized = T)

pare25btw <-betweenness(pare25i)
#pare25btwN <-betweenness(pare25i, normalized = T)


#Diagonal - Local retention
pareLR <-diag(as.matrix(pare))
pare5LR <-diag(as.matrix(pare5))
pare15LR <-diag(as.matrix(pare15))
pare25LR <-diag(as.matrix(pare25))

#ok; for active 5,15 and 30 = out, in, btw and lr ready
#1.3 One object with nodes as column 

paremetrics<-as.data.frame(colnames(pare5))
colnames(paremetrics) <- c("ID")

dim(paremetrics) #ok 25 columns
dim(cryptometrics) #ok 25 columns

#IN
paremetrics$parein <- parein
paremetrics$pare5in <- pare5in
paremetrics$pare15in <- pare15in
paremetrics$pare25in <- pare25in

paremetrics$pareinN <- pareinN
paremetrics$pare5inN <- pare5inN
paremetrics$pare15inN<- pare15inN
paremetrics$pare25inN<- pare25inN


#OUT
paremetrics$pareout <- pareout
paremetrics$pare5out <- pare5out
paremetrics$pare15out <- pare15out
paremetrics$pare25out <- pare25out

paremetrics$pareoutN <- pareoutN
paremetrics$pare5outN <- pare5outN
paremetrics$pare15outN<- pare15outN
paremetrics$pare25outN <- pare25outN

#BTW
paremetrics$parebtw <- parebtw
paremetrics$pare5btw <- pare5btw
paremetrics$pare15btw <- pare15btw
paremetrics$pare25btw <- pare25btw

#LR
paremetrics$pare5LR <- pare5LR
paremetrics$pare15LR <- pare15LR
paremetrics$pare25LR <- pare25LR
paremetrics$pareLR <- pareLR

dim(paremetrics) #ok, only missing for PASsive model yohoo

#3rd Resident group - (missing csv files for active/missing passive)
residpas<-read.csv("/Volumes/LuisaDrive/Model_outputs/Resident/Passive/con_file_pre35.csv",h=FALSE)
residact5<-read.csv("/Volumes/LuisaDrive/Model_outputs/Resident/Active/con_file_pre5.csv",h=FALSE)
residact15<-read.csv("/Volumes/LuisaDrive/Model_outputs/Resident/Active/con_file_pre15.csv",h=FALSE)
residact30<-read.csv("/Volumes/LuisaDrive/Model_outputs/Resident/Active/con_file_pre30.csv",h=FALSE)

#need to replane NaN -> 0
resid <- replace(residpas, residpas == "NaN", 0)
rownames(resid) <-colnames(resid)
dim(resid) #ok
resid5 <- replace(residact5, residact5 == "NaN", 0)
rownames(resid5) <-colnames(resid5)
dim(resid5) #ok
resid15 <- replace(residact15, residact15 == "NaN", 0)
rownames(resid15) <-colnames(resid15)
dim(resid15) #ok
resid35 <- replace(residact30, residact30 == "NaN", 0)
rownames(resid35) <-colnames(resid35)
dim(resid35) #ok

#Igraph - metrics

residi<-graph.adjacency(as.matrix(resid),diag=T,mode="directed",weighted = TRUE)
resid5i<-graph.adjacency(as.matrix(resid5),diag=T,mode="directed",weighted = TRUE)
resid15i<-graph.adjacency(as.matrix(resid15),diag=T,mode="directed",weighted = TRUE)
resid35i<-graph.adjacency(as.matrix(resid35),diag=T,mode="directed",weighted = TRUE)

#In degree
residin <-degree(residi, mode="in")
residinN <-degree(residi, mode="in", normalized = T)

resid5in <-degree(resid5i, mode="in")
resid5inN <-degree(resid5i, mode="in", normalized = T)

resid15in <-degree(resid15i, mode="in")
resid15inN <-degree(resid15i, mode="in", normalized = T)

resid35in <-degree(resid35i, mode="in")
resid35inN <-degree(resid35i, mode="in", normalized = T)

#Out degree
residout <-degree(residi, mode="out")
residoutN <-degree(residi, mode="out", normalized = T)

resid5out <-degree(resid5i, mode="out")
resid5outN <-degree(resid5i, mode="out", normalized = T)

resid15out <-degree(resid15i, mode="out")
resid15outN <-degree(resid15i, mode="out", normalized = T)

resid35out <-degree(resid35i, mode="out")
resid35outN <-degree(resid35i, mode="out", normalized = T)

#Between

residbtw <-betweenness(residi)
#residbtwN <-betweenness(residi, normalized = T)

resid5btw <-betweenness(resid5i)
#resid5btwN <-betweenness(resid5i, normalized = T)

resid15btw <-betweenness(resid15i)
#resid15btwN <-betweenness(resid15i, normalized = T)

resid35btw <-betweenness(resid35i)
#resid35btwN <-betweenness(resid35i, normalized = T)


#Diagonal - Local retention
residLR <-diag(as.matrix(resid))
resid5LR <-diag(as.matrix(resid5))
resid15LR <-diag(as.matrix(resid15))
resid35LR <-diag(as.matrix(resid35))

#ok; for active 5,15 and 30 = out, in, btw and lr ready
#1.3 One object with nodes as column 

residmetrics<-as.data.frame(colnames(resid5))
colnames(residmetrics) <- c("ID")
dim(residmetrics)

#IN
residmetrics$residin <- residin
residmetrics$resid5in <- resid5in
residmetrics$resid15in <- resid15in
residmetrics$resid35in <- resid35in

residmetrics$residinN <- residinN
residmetrics$resid5inN <- resid5inN
residmetrics$resid15inN<- resid15inN
residmetrics$resid35inN <- resid35inN

#OUT
residmetrics$residout <- residout
residmetrics$resid5out <- resid5out
residmetrics$resid15out <- resid15out
residmetrics$resid35out <- resid35out

residmetrics$residoutN <- residoutN
residmetrics$resid5outN <- resid5outN
residmetrics$resid15outN<- resid15outN
residmetrics$resid35outN <- resid35outN

#BTW
residmetrics$residbtw <- residbtw
residmetrics$resid5btw <- resid5btw
residmetrics$resid15btw <- resid15btw
residmetrics$resid35btw <- resid35btw

#LR
residmetrics$residLR <- residLR
residmetrics$resid5LR <- resid5LR
residmetrics$resid15LR <- resid15LR
residmetrics$resid35LR <- resid35LR

dim(residmetrics) #ok 25 columns

####


#4th Transient group - (passive csv missing)
transipas<-read.csv("/Volumes/LuisaDrive/Model_outputs/Transient/Passive/con_file_pre30.csv",h=FALSE)
transiact5<-read.csv("/Volumes/LuisaDrive/Model_outputs/Transient/Active/con_file_pre5.csv",h=FALSE)
transiact15<-read.csv("/Volumes/LuisaDrive/Model_outputs/Transient/Active/con_file_pre15.csv",h=FALSE)
transiact30<-read.csv("/Volumes/LuisaDrive/Model_outputs/Transient/Active/con_file_pre30.csv",h=FALSE)

#need to replane NaN -> 0
transi <- replace(transipas, transipas == "NaN", 0)
rownames(transi) <-colnames(transi)
dim(transi) #ok
transi5 <- replace(transiact5, transiact5 == "NaN", 0)
rownames(transi5) <-colnames(transi5)
dim(transi5) #ok
transi15 <- replace(transiact15, transiact15 == "NaN", 0)
rownames(transi15) <-colnames(transi15)
dim(transi15) #ok
transi35 <- replace(transiact30, transiact30 == "NaN", 0)
rownames(transi35) <-colnames(transi35)
dim(transi35) #ok

#Igraph - metrics

transii<-graph.adjacency(as.matrix(transi),diag=T,mode="directed",weighted = TRUE)
transi5i<-graph.adjacency(as.matrix(transi5),diag=T,mode="directed",weighted = TRUE)
transi15i<-graph.adjacency(as.matrix(transi15),diag=T,mode="directed",weighted = TRUE)
transi35i<-graph.adjacency(as.matrix(transi35),diag=T,mode="directed",weighted = TRUE)

#In degree
transiin <-degree(transii, mode="in")
transiinN <-degree(transii, mode="in", normalized = T)

transi5in <-degree(transi5i, mode="in")
transi5inN <-degree(transi5i, mode="in", normalized = T)

transi15in <-degree(transi15i, mode="in")
transi15inN <-degree(transi15i, mode="in", normalized = T)

transi35in <-degree(transi35i, mode="in")
transi35inN <-degree(transi35i, mode="in", normalized = T)

#Out degree
transiout <-degree(transii, mode="out")
transioutN <-degree(transii, mode="out", normalized = T)

transi5out <-degree(transi5i, mode="out")
transi5outN <-degree(transi5i, mode="out", normalized = T)

transi15out <-degree(transi15i, mode="out")
transi15outN <-degree(transi15i, mode="out", normalized = T)

transi35out <-degree(transi35i, mode="out")
transi35outN <-degree(transi35i, mode="out", normalized = T)

#Between

transibtw <-betweenness(transii)
#transibtwN <-betweenness(transii, normalized = T)

transi5btw <-betweenness(transi5i)
#transi5btwN <-betweenness(transi5i, normalized = T)

transi15btw <-betweenness(transi15i)
#transi15btwN <-betweenness(transi15i, normalized = T)

transi35btw <-betweenness(transi35i)
#transi35btwN <-betweenness(transi35i, normalized = T)


#Diagonal - Local retention
transiLR <-diag(as.matrix(transi))
transi5LR <-diag(as.matrix(transi5))
transi15LR <-diag(as.matrix(transi15))
transi35LR <-diag(as.matrix(transi35))

#ok; for active 5,15 and 30 = out, in, btw and lr ready
#1.3 One object with nodes as column 

transimetrics<-as.data.frame(colnames(transi5))
colnames(transimetrics) <- c("ID")
dim(transimetrics)

#IN
transimetrics$transiin <- transiin
transimetrics$transi5in <- transi5in
transimetrics$transi15in <- transi15in
transimetrics$transi35in <- transi35in

transimetrics$transiinN <- transiinN
transimetrics$transi5inN <- transi5inN
transimetrics$transi15inN<- transi15inN
transimetrics$transi35inN <- transi35inN

#OUT
transimetrics$transiout <- transiout
transimetrics$transi5out <- transi5out
transimetrics$transi15out <- transi15out
transimetrics$transi35out <- transi35out

transimetrics$transioutN <- transioutN
transimetrics$transi5outN <- transi5outN
transimetrics$transi15outN<- transi15outN
transimetrics$transi35outN <- transi35outN

#BTW
transimetrics$transibtw <- transibtw
transimetrics$transi5btw <- transi5btw
transimetrics$transi15btw <- transi15btw
transimetrics$transi35btw <- transi35btw

#LR
transimetrics$transiLR <- transiLR
transimetrics$transi5LR <- transi5LR
transimetrics$transi15LR <- transi15LR
transimetrics$transi35LR <- transi35LR

dim(transimetrics) #ok 25 columns

####

#Crypto x Parental (max 5cm TL and demersal eggs - variation in spawning season)
#IN
cor(paremetrics$pare5in,cryptometrics$crypto5in)
plot(paremetrics$pare5in,cryptometrics$crypto5in)

#OUT
cor(paremetrics$pare5out,cryptometrics$crypto5out)
plot(paremetrics$pare5out,cryptometrics$crypto5out)

#BTW
cor(paremetrics$pare5btw,cryptometrics$crypto5btw)
plot(paremetrics$pare5btw,cryptometrics$crypto5btw)

#LR
cor(paremetrics$pare15LR,cryptometrics$crypto15LR)
plot(paremetrics$pare15LR,cryptometrics$crypto15LR)

#Transient x Parental (max 5cm TL and demersal eggs - variation in spawning season)
#IN
cor(paremetrics$pare5in,transimetrics$transi5in)
plot(paremetrics$pare5in,transimetrics$transi5in)

#OUT
cor(paremetrics$pare5out,transimetrics$transi5out)
plot(paremetrics$pare5out,transimetrics$transi5out)


#Build massive database
globalmetricT <- do.call("rbind", list(cryptometrics,paremetrics,residmetrics,transimetrics))


global_metrics<-merge(paremetrics,cryptometrics,by=c("ID")) 
global_metrics<-merge(global_metrics,transimetrics,by=c("ID"))
global_metrics<-merge(global_metrics,residmetrics,by=c("ID"))

dim(global_metrics)
write.csv(global_metrics,"global_metrics.csv")


