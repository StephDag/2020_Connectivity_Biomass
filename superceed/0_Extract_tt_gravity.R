# Aim : Extract travel time and gravity for each cell
# Author: Stephanie D'agata
# Date: March 2020
# ** Output: adding gravity and travel time

rm(list=ls())

# packages
library(here)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(raster)
library(rgeos)
library(rgdal)

# load data
# coral data
CORAL <- read.csv(here("data","final","Coral_Data_09032020.csv"))
CORAL %>% head()
CORAL %>% dim()

# fish data
FISH <- read.csv(here("data","final","Fish_Data_09032020.csv"))
FISH %>% head()
FISH %>% dim()

# load travel time 7km
path.tt <- "/Users/stephdagata/OneDrive - Macquarie University/One Drive Steph/OneDrive - Macquarie University/TravelTime"
tt.7 <- raster(paste(path.tt,"TravelTime_7k.tif",sep="/"))
plot(tt.7)

# load travel time 15km
path.tt <- "/Users/stephdagata/OneDrive - Macquarie University/One Drive Steph/OneDrive - Macquarie University/TravelTime"
tt.15 <- raster(paste(path.tt,"TravelTime_15k.tif",sep="/"))
plot(tt.15)

# load gravity time 7km
path.grav <- "/Users/stephdagata/OneDrive - Macquarie University/One Drive Steph/OneDrive - Macquarie University/TravelTime"
grav.7 <- raster(paste(path.grav,"logGRAV_WIO_7k.tif",sep="/"))
plot(grav.7)

# load gravity time 15km
path.grav <- "/Users/stephdagata/OneDrive - Macquarie University/One Drive Steph/OneDrive - Macquarie University/TravelTime"
grav.15 <- raster(paste(path.grav,"logGRAV_WIO_15k.tif",sep="/"))
plot(grav.15)

# load gravity Cinner
path.grav.cinner <- "/Users/stephdagata/Dropbox/MISC DATA/Cinner&al2018_Gravity/PNASGlobalGravity"
grav.15.cinner <- readOGR(dsn=path.grav.cinner,layer="Total Gravity of Coral Reefs 1.0")
plot(grav.15.cinner)
crs(grav.15.cinner)

# Extract value for the CORAL data
  # 1. create a spatial data point for the coordinates
  coords.CORAL %>% rm()
  coords.CORAL = cbind(CORAL$Lon, CORAL$Lat)
  sp.CORAL = SpatialPoints(coords.CORAL)
  points(sp.CORAL)
  crs(sp.CORAL) <- crs(grav.15.cinner)
  # 2. Extract values tt 7km
  tt7k.CORAL.extr <- raster::extract(tt.7,sp.CORAL)
  tt.7k <- tt7k.CORAL.extr
  tt.7k %>% head()
  
  # 3. Extract values tt 15km
  tt15k.CORAL.extr <- raster::extract(tt.15,sp.CORAL)
  tt.15k <- tt15k.CORAL.extr
  tt.15k %>% head()
  
  # 4. Extract values grav 7km
  grav7k.CORAL.extr <- raster::extract(grav.7,sp.CORAL)
  loggrav.7k <- grav7k.CORAL.extr
  grav.7k %>% head()
  
  # 5. Extract values grav 15km
  grav15k.CORAL.extr <- raster::extract(grav.15,sp.CORAL)
  loggrav.15k <- grav15k.CORAL.extr
  grav.15k %>% head()
  
  # 5. Extract values grav 15km cinner
  grav15k.CORAL.cinner <- over(sp.CORAL,grav.15.cinner)
  grav15k.CORAL.cinner %>% head()
  grav15k.CORAL.cinner %>% summary()
  
# bind all for Coral
CORAL <- cbind(CORAL,tt.7k,tt.15k,grav.7k,grav.15k,log1p(grav15k.CORAL.cinner))
CORAL %>% head()  
CORAL %>% summary() 

# save CORAL data
saveRDS(CORAL,here("data","final","Coral_Data.rds"))

# Extract value for the FISH data
# 1. create a spatial data point for the coordinates
coords.FISH %>% rm()
coords.FISH = cbind(FISH$Lon, FISH$Lat)
sp.FISH = SpatialPoints(coords.FISH)
points(sp.FISH)

# 2. Extract values tt 7km
tt7k.FISH.extr <- raster::extract(tt.7,sp.FISH)
tt.7k <- tt7k.FISH.extr
tt.7k %>% head()

# 3. Extract values tt 15km
tt15k.FISH.extr <- raster::extract(tt.15,sp.FISH)
tt.15k <- tt15k.FISH.extr
tt.15k %>% head()

# 4. Extract values grav 7km
grav7k.FISH.extr <- raster::extract(grav.7,sp.FISH)
loggrav.7k <- grav7k.FISH.extr
loggrav.7k %>% head()

# 5. Extract values grav 15km
grav15k.FISH.extr <- raster::extract(grav.15,sp.FISH)
loggrav.15k <- grav15k.FISH.extr
loggrav.15k %>% head()

# 5. Extract values grav 15km cinner
crs(sp.FISH) <- crs(grav.15.cinner)
grav15k.FISH.cinner <- over(sp.FISH,grav.15.cinner)
grav15k.FISH.cinner %>% head()
grav15k.FISH.cinner %>% summary()

# bind all for FISH
FISH %>% rm()
FISH <- cbind(FISH,tt.7k,tt.15k,loggrav.7k,loggrav.15k,log1p(grav15k.FISH.cinner))
FISH %>% head()  
FISH %>% summary()

# save FISH data
saveRDS(FISH,here("data","final","Fish_Data.rds"))


