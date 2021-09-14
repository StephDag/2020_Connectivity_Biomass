#Chapter 4 - Connectivity project
library(data.table)

##PLD distributions----
#1st PLD database
traits<-read.csv("~/OneDrive - Macquarie University/Chapters/Chapter_04_on_going/R analyses/Chapter_4/PLD_LATEST.csv", h=T)
head(traits)
pld<-traits
library(ggplot2)

#------Cryptobenthic - species smaller than 5cm and demersal eggs (Brand definition)
pld_crypto<-traits[which(traits$Body_size <5 & traits$Spawn == "DEM"),]
unique(pld_crypto$Family) #definition of cryptobenthic by Simon Branld

#removing pomacentridae
pld_crypto<-pld_crypto[pld_crypto$Family %in% c("APOGONIDAE","CHAENOPSIDAE",
                                                "DACTYLOSCOPIDAE","GOBIIDAE",
                                                "TRIPTERYGIIDAE"),] 


mini<-qplot(PLD_mean, data=pld_crypto, fill="PLD_mean", geom="density", alpha=I(.5)) +
  theme(legend.position="none") + theme_classic() + labs(y="Density", x="PLD") + ggtitle("Cryptobenthic") +
  annotate("text", x = 50, y = 0.08, label = "Mean PLD 32+-11 days ") +
  scale_fill_manual(values=c("#999999")) + xlim(0,100)


mean(pld_crypto$PLD_mean) #32 days
sd(pld_crypto$PLD_mean) # sd +- 11 days
unique(pld_crypto$Genus_species) #19 species

#--------------------------------------Small body size, demersal eggs, parental care
pld_parental<-pld[pld$Group =="Parental care",]
pld_parental<-pld_parental[pld_parental$Body_size>6,]
pld_parental<-pld_parental[-which(pld_parental$Family == "DACTYLOSCOPIDAE"),]
unique(pld_parental$Genus_species) #166 species  
max(pld_parental$Body_size) #6-31cm body size variation


parental<-qplot(PLD_mean, data=pld_parental, fill="PLD_mean", geom="density", alpha=I(.5)) +
  theme(legend.position="none") + theme_classic() + labs(y="Density", x="PLD") + ggtitle("Parental care > 5cm TL") +
  annotate("text", x = 50, y = 0.08, label = "Mean PLD 23+-7 days ") +
  scale_fill_manual(values=c("#E69F00")) + xlim(0,100)



mean(pld_parental$PLD_mean) # 23
sd(pld_parental$PLD_mean) # sd +- 7


#Large transient species - spawnning agreggation limited
head(pld)
pld_transient<-pld[pld$Group == "Transient",]
pld_transient <-pld_transient[pld_transient$Body_size > 39,]
mean(pld_transient$PLD)#34 dias
sd(pld_transient$PLD) #+-12
unique(pld_transient$Genus_species) #23 species


trans<-qplot(PLD_mean, data=pld_transient, fill="PLD_mean", geom="density", alpha=I(.5)) +
  theme(legend.position="none") + theme_classic() + labs(y="Density", x="PLD") + ggtitle("Transient spawnning") +
  annotate("text", x = 50, y = 0.08, label = "Mean PLD 36+-12 days ") +
  scale_fill_manual(values=c("#56B4E9")) + xlim(0,100)


#Resident and pelagic

pld_resident<-pld[pld$Group == "Resident",]
mean(pld_resident$PLD_mean) #37 
sd(pld_resident$PLD_mean) # +-15.8
View(pld_resident)

resi<-qplot(PLD_mean, data=pld_resident, fill="PLD_mean", geom="density", alpha=I(.5)) +
  theme(legend.position="none") + theme_classic() + labs(y="Density", x="PLD") + ggtitle("Resident spawnning") +
  annotate("text", x = 50, y = 0.08, label = "Mean PLD 37+-15 days ") +
  scale_fill_manual(values=c("red")) + xlim(0,100)

ggarrange(mini,parental,resi,trans,ncol = 4, nrow = 1,
          common.legend = TRUE, legend = "none")




####Biomass----
databiomass<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_on_going/Prior biomass and PLD analyses/Chapter_4/data_geb_paper_for_luisa.csv", h=T)
head(databiomass)
databiomass$individual_biomass_g  <-  databiomass$coefa * (databiomass$size_cm * databiomass$ratio)^databiomass$coefb
#View(new.data[new.data$region == "westernindian",])
#Filtering data
#Fique atenta pq nem todos os censos contém abundância 
#e classe de tamanho. Caso vc queira filtrar só pra esses, 
#vc precisa usar os que são 'quant_class' na coluna type
databiomassF<-databiomass[databiomass$type == "quant_class",]
droplevels(databiomassF$type)
#calculating biomass
databiomassF$totalbiom<-databiomassF$abun*databiomassF$individual_biomass_g

#filtering genus according McNeil
genus_filter<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_on_going/Prior biomass and PLD analyses/Chapter_4/genus_list_mcneil.csv", h=T)

unique(genus_filter$genus) #270 genus (try to fit the model using the same genus used for McNeil paper)
levels(genus_filter$genus) <- tolower(levels(genus_filter$genus))
#filtering by genus (followin Aaron paper Nature)
library(dplyr)
library(stringr)
genuslist <- unique(genus_filter$genus)
genuslist<-as.vector(genuslist) #268 genus
#creating a column genus
databiomassF$species<-gsub("_", " ", databiomassF$species)
unique(databiomassF$species)
databiomassF$genus<- word(databiomassF$species, 1)
unique(databiomassF$genus) #540 total

new.data<-subset(databiomassF, genus %in% genuslist)
new.data2<-subset(databiomassF, genus %in% genuslist)

droplevels(new.data)
droplevels(new.data2)

#optional
#remove potential aggregations of large species and mobbing behavior ( bias - following mcneil)
new.data <- new.data[ -which( new.data$size_cm >30 &  new.data$abun >  100) , ]

richness<-as.data.frame(setDT(new.data)[, .(count = uniqueN(species)), by = sites])
richness2<-as.data.frame(setDT(new.data2)[, .(count = uniqueN(species)), by = sites])



#1st  - Biomass filtered (genus and schooling)
##Sum of biomass per site/total sampling area (sum of biomass site/sum of area per site) 
totalbiosites<-aggregate(totalbiom ~ sites + locality, new.data, sum)
areatotal<-aggregate(area~transect_id + sites + locality, new.data,max)
areabysite<-aggregate(area ~ sites + locality,areatotal,sum)
totalbiomass<-merge(totalbiosites,areabysite,by=c("locality","sites"))
totalbiomass$biomassarea1<-totalbiomass$totalbiom/totalbiomass$area
head(totalbiomass) #460 sites
unique(totalbiomass$sites)

#2st  - Biomassfiltered by genus but NOT filtered by schooling)
##Sum of biomass per site/total sampling area (sum of biomass site/sum of area per site) 
totalbiosites2<-aggregate(totalbiom ~ sites + locality, new.data2, sum)
areatotal2<-aggregate(area~transect_id + sites + locality, new.data2,max)
areabysite2<-aggregate(area ~ sites + locality,areatotal2,sum)
totalbiomass2<-merge(totalbiosites2,areabysite2,by=c("locality","sites"))
totalbiomass2$biomassarea2<-totalbiomass2$totalbiom/totalbiomass2$area
unique(totalbiomass2$sites) #460 sites


#Create a biomass dataset with 2 biomass estimatives (merge)

BiomassData<-merge(totalbiomass,totalbiomass2, by=c("locality","sites","area"))
BiomassDataFull<-BiomassData[,c("locality","sites","area","biomassarea1","biomassarea2")]
head(BiomassDataFull) #460 sites

library(dplyr)
#Include coordinates and age of protection 
##Age of protection
protection_mpa <- read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_on_going/Prior biomass and PLD analyses/Chapter_4/age_of_protection_Sites.csv", h=T)
gravitydata <- read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_on_going/Prior biomass and PLD analyses/Chapter_4/coordinates_gravity_chpt4.csv", h=T)

dataG<-merge(BiomassDataFull,gravitydata, by=c("locality","sites"))
dataGlobal<-merge(dataG,protection_mpa, by=c("locality","sites","region"))
unique(dataGlobal$sites) #456 sites 
dataGlobal<-left_join(dataGlobal,richness,by="sites")
dataGlobal<-left_join(dataGlobal,richness2,by="sites")
unique(dataGlobal$sites)
colnames(dataGlobal)[12] <- c("richness1")
colnames(dataGlobal)[13] <- c("richness2")
head(dataGlobal)


#Outline figures----
testando <- transectbio %>% 
  group_by(sites) %>% 
  summarise(count = n())
hist(testando[testando$count < 100,]$count, 60)


#testing the effect bewteen number of transects and biomass/area
transecteffect<-merge(transectbio, testando, by="sites",all=T)
head(transecteffect)

#mean biomass/area - site
globalbio<-aggregate(biomassarea ~ sites + locality + region + count , transecteffect,mean)

globalbio$logbiomassarea<-log(globalbio$biomassarea)
globalbio$logcount<-log(globalbio$count)

test<-lm(logbiomassarea ~ logcount, data=globalbio)
summary(test)
visreg::visreg(test)

#unique(globalbio[globalbio$count<50,]$sites) #401 if less than 50 transects
#no effect then (p >0.05 - not-significant)

#distribution of total numberof transects
par(mfrow=c(1,2))
visreg::visreg(test, "logcount",yaxt ="n",
               xlab="Log Total number of transects/site", ylab="Log Biomass/area",main="",cex.lab=1.3,
               points=list(col="palevioletred4",cex=0.6),
               line=list(col="palevioletred4"),
               fill=list(col="pink"))
hist(globalbio$count, 100, col="palevioletred4", main="", xlab="Total number of transects/site")
max(globalbio$count) #618
min(globalbio$count) #3

#
##
bioandmpa<-merge(globalbio,maisum[,c("region","locality","sites","Age_of_protection","Class")],
by=c("region","sites","locality"), all=T)
unique(bioandmpa$Class) #98
head(bioandmpa)
bioandmpa <- bioandmpa[which(bioandmpa$Class %in% c("Closed","Fished","Restricted")) , ]


##boxplot

boxlogbio<-ggplot(data=bioandmpa, aes(y=logbiomassarea,x=as.factor(Class), fill=as.factor(Class))) + geom_boxplot(alpha=I(1)) +
  scale_fill_manual(values=c("forestgreen","red","navy"),name  ="Status of protection",
                    breaks=c("Closed", "Fished", "Restricted"),
                    labels=c("Closed", "Fished", "Restricted")) + labs(title="", y="Log Total Biomass/area (g/m2)", 
                                                             x="") +
  theme( panel.grid.major.y = element_line(colour = "white", size = NULL, linetype = NULL,  # horizontale Linien
                                           lineend = NULL)
         ,panel.grid.minor.y = element_line(colour = "white", size = NULL, linetype = NULL,
                                            lineend = NULL)
         ,panel.grid.major.x = element_blank()           # vertikale Linien
         ,panel.grid.minor.x = element_blank()
         ,legend.background = element_rect(fill = "white", colour = "white") # Legende 
         ,legend.key = element_rect(fill = "white", colour = "white")
         ,panel.background = element_rect(fill = "white", colour = "white", size = NULL, # Panel Hintergrund
                                          linetype = NULL)
         ,axis.line = element_line(colour = "black", size=.5)
         ,axis.title.x = element_text(family = NULL, face = "bold", size = 11,vjust=0.1)
         ,axis.title.y = element_text(family = NULL, face = "bold", size = 11,vjust=0.1)
         ,axis.text=element_text(colour="black")
         ,legend.title = element_text(family = NULL, face = "plain", size = 11)
         ,legend.text = element_text(family = NULL, face = "plain", size = 9))



library(ggplot2)
library(mgcv)
head(bioandmpa)




##
#meanbio<-aggregate(biomassarea ~ locality + region + Age_of_protection + Class, bioandmpa, mean)
#unique(meanbio$locality)
#head(meanbio)
library(mgcv)
test2<-gam(logbiomassarea ~ log1p(Age_of_protection),REML=F, data=bioandmpa)
test3<-gam(logbiomassarea ~ log1p(Age_of_protection),REML=F, data=bioandmpa[bioandmpa$Class == c("Closed"),])
test4<-gam(logbiomassarea ~ log1p(Age_of_protection),REML=F, data=bioandmpa[bioandmpa$Class == c("Restricted"),])


visreg::visreg(test2, "Age_of_protection")
visreg::visreg(test3, "Age_of_protection")
visreg::visreg(test4, "Age_of_protection")

par(mfrow=c(1,3))
#All
visreg::visreg(test2, "Age_of_protection",
               xlab="Age of protection", ylab="Log Total Fish biomass (g/m2)",main="",cex.lab=1.3,
               points=list(col=rgb(0,0,0.2,0.2),cex=0.6),
               line=list(col=rgb(0.2,0,1,1)),
               fill=list(col=rgb(0.1,0.5,1,0.5)))

#Closed MAP's areas
visreg::visreg(test3, "Age_of_protection",
               xlab="Age of protection of No-take MPAs", ylab="Log Total Fish biomass (g/m2)",main="",cex.lab=1.3,
               points=list(col=rgb(0,0,0.2,0.2),cex=0.6),
               line=list(col=rgb(0.2,0,1,1)),
               fill=list(col=rgb(0.1,0.5,1,0.5)))

#Restricted MPA's
visreg::visreg(test4, "Age_of_protection",
               xlab="Age of protection of MPAs locally managed (Restricted gear)", ylab="Log Total Fish biomass (g/m2)",main="",cex.lab=1.3,
               points=list(col=rgb(0,0,0.2,0.2),cex=0.6),
               line=list(col=rgb(0.2,0,1,1)),
               fill=list(col=rgb(0.1,0.5,1,0.5)))




#write.csv(mpaclosed,"mpabiosites.csv")
##





#####-------- plotting MAP
library("ggmap")
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())
library(maptools)
library(maps)
mp <- NULL
coord<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_on_going/R analyses/Data_chapter_4/geogrVar.csv", h=T)
#View(coord)
##plotting biomass #
globalbiot<-merge(maisum,coord[,c("sites","lon","lat")], by ="sites",all=F)
database<-merge(bioandmpa,coord[,c("sites","lon","lat")], by ="sites",all=F )
head(database)
write.csv(bioan, "sites_coordinates.csv") #send it to Majambo

#####Plotting Map
mapWorld <- borders("world", colour="gray", fill="gray") # create a layer of borders

mp <- ggplot() +   mapWorld + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank()) +
  coord_map(ylim=c(-50, 50), xlim = c(-200,200))
####
mp3 <- mp + geom_point(aes(x=globalbiot$lon, y=globalbiot$lat, fill=factor(globalbiot$Class)), 
                       alpha=0.35, size=3, shape=21, colour="black") +
  scale_fill_manual(values=c("forestgreen","red","navy","gold")) + labs(x="",y="")

####
ggplot(globalbiot, aes(x=Age_of_protection, fill=Class, alpha=I(.9))) + geom_histogram(binwidth = 10) +
  facet_wrap(~region, ncol=3) + theme_classic() + scale_fill_manual(values=c("forestgreen","red","navy","gold")) + labs(y="Density", x="Age of Protection")

ggplot(globalbiot, aes(x=Age_of_protection, fill=Class, alpha=I(.9))) + geom_histogram(binwidth = 5) +
  theme_classic() + scale_fill_manual(values=c("forestgreen","red","navy")) + labs(y="Density", x="Age of Protection")






#Functional diversity----
datatrait<-read.csv("~/Desktop/OneDrive - Macquarie University/Chapters/Chapter_04_on_going/Prior biomass and PLD analyses/Chapter_4/checktraits.csv")
head(datatrait) #dim 1900 rows
library(FD)
library(cluster)
library(ape)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)

#Schooling F IM L M P S
schooling <- sapply(datatrait$Schooling, function(x) {
  if (x=="S") {"1"} else
    if (x=="F") {"2"} else
      if (x=="P") {"3"} else
        if (x=="IM") {"4"} else
          if (x=="M") {"4"} else
            if (x=="L") {"5"}
})

datatrait$schooling <- schooling
unique(datatrait$schooling)

trait<-datatrait[,c("Diet_2012","Size_class","schooling")]
rownames(trait) <- datatrait$species
colnames(trait) <- c("diet","size","school")
trait$size <- as.ordered(trait$size)
trait$school <- as.ordered(trait$school)

#1.4 Reordering rownames alphabetically
trait2 <- trait[ order(row.names(trait)), ]
head(trait2)

#1.5 Abundance matrix
#reading data
head(new.data)
fishbio<-aggregate(abun ~ species + sites, new.data, sum)
head(fishbio)
#1.6 biomass matrix per year and site 
fish3<-spread(fishbio,species,abun)
dim(fish3)
fish3[is.na(fish3)] <- 0
rownames(fish3)<-fish3[,1]
fish3<-fish3[,-1]
#1.7 Binary matrix for betadiversity computation
fishdatatrait<-as.matrix((fish3 > 0) + 0)
#Setting matrices (traits, pcoa and biomass)#setting datasets (traits and biomass matrix)
#####
# 2.1 PCOA with Gower distance with 6 traits
rm(d.traits)
trait2$school<-as.numeric(trait2$school)
trait2$size <- as.numeric(trait2$size)
head(trait2)
d.traits <- daisy(trait2,metric="gower") # compute the distance matrix with GOWER distance with the 3 traits
#pcoa.traits <-pcoa(d.traits, correction = "cailliez") #not working cailliez
pcoa.traitsRAW <- ape::pcoa(d.traits)

pcoa.trait.cor<-pcoa.traitsRAW$vectors[,1:4]

rm(LIZ_FD)
Global_FD <- multidimFD(pcoa.trait.cor,fishdatatrait %>% as.matrix)
fd<-as.data.frame(Global_FD)
fd$sites<-rownames(fd)
head(fd)
fd<-as.data.frame(fd[,c("FRic","sites")])

dataGlobal2<-merge(dataGlobal,fd,by="sites")
unique(dataGlobal2$sites) #456 sites (OK)

cor(dataGlobal2$richness2,dataGlobal2$FRic)# Species richness and Functional Richness - highly correlated 0.85 


#Biomass of target species----
datafulltraits<-merge(new.data,datatrait,by="species",h=T)
head(datafulltraits)

##Sum of biomass per site/total sampling area (sum of biomass site/sum of area per site) 
totalbiositesTR<-aggregate(totalbiom ~ sites + locality, datafulltraits[datafulltraits$Diet_2012 %in% c("FC","HD","HD","IM") & datafulltraits$Size_class>4,], sum)
unique(areabysite$sites) #460
unique(totalbiositesTR$sites)#460

totalbiomassTR<-merge(totalbiositesTR,areabysite,by=c("locality","sites"))
totalbiomassTR$biomassareaTR<-totalbiomassTR$totalbiom/totalbiomassTR$area

lastfulldata<-merge(dataGlobal2,totalbiomassTR,by=c("locality","sites","area"),all=T)
lastfulldata[is.na(lastfulldata)] <- 0
lastfulldataT<-na.omit(lastfulldata)
unique(lastfulldataT$sites) #457 sites
dim(lastfulldataT)
lastfulldataF<-lastfulldataT[!duplicated(lastfulldataT$sites),]
dim(lastfulldataF) #456 ok

head(lastfulldataF)
#write.csv(lastfulldataF, "databiomassFull.csv")


