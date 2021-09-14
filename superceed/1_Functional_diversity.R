# Functional diversity - Global connectivity 
# author: Luisa Fontura/Steph D'agata
# date: March 2020
# outputs: PCOA and functional diversity for global study

rm(list=ls())

# load packages
library(here)
library(FD)
library(cluster)
library(ape)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)

#Functional diversity----
datatrait<-read.csv(here("_data","checktraits.csv"))
head(datatrait) #dim 1900 rows
datatrait %>% summary()

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

# check traits matric ordinal/categorical variables
summary(trait2)

# 2.1 PCOA with Gower distance with 6 traits
rm(d.traits)
#trait2$school<-as.numeric(trait2$school)
#trait2$size <- as.numeric(trait2$size)
head(trait2)
d.traits %>% rm()
d.traits <- daisy(x=trait2,metric="gower") # compute the distance matrix with GOWER distance with the 3 traits
#pcoa.traits <-pcoa(d.traits, correction = "cailliez") #not working cailliez
pcoa.traitsRAW %>% rm()
pcoa.traitsRAW <- ape::pcoa(d.traits)

pcoa.trait.cor<-pcoa.traitsRAW$vectors[,1:4]


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

rm(LIZ_FD)
Global_FD <- multidimFD(pcoa.trait.cor,fishdatatrait %>% as.matrix)
fd<-as.data.frame(Global_FD)
fd$sites<-rownames(fd)
head(fd)
fd<-as.data.frame(fd[,c("FRic","sites")])

dataGlobal2<-merge(dataGlobal,fd,by="sites")
unique(dataGlobal2$sites) #456 sites (OK)

cor(dataGlobal2$richness2,dataGlobal2$FRic)# Species richness and Functional Richness - highly correlated 0.85 

