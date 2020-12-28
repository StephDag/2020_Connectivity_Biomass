library(rgdal)
library(raster)
mypath<-"/Users/josephmaina/Documents/Mygitprojects/2020_Connectivity_Biomass/_data/reefGrids/"

reef.grids<-readOGR(dsn=mypath, "settlement")

wiompagdb<-"~/Dropbox/WIO_MPA_Outlook_Maps_Database/GeoDatabase/WIO_MPA_Outlook_Feb2020.gdb"
subset(ogrDrivers(), grepl("GDB", name))
fc_list <- ogrListLayers(wiompagdb)
print(fc_list)

#wiompa<-readOGR(dsn=wiompagdb, "WIO_MPAs_Combined",require_geomType = 'wkbPolygon')
require(sf)
wiompa <- sf::st_read(wiompagdb, layer = "WIO_MPAs_Combined")
str(wiompa)
nc_sp <- sf:::as_Spatial(nc$geom)
head(wiompa)
levels(wiompa$Designation)


