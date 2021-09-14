library(raster)
library(rgdal)
library(sp)
library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")
#devtools::install_github("ropenscilabs/rnaturalearthdata")
#install.packages("rnaturalearthhires",
#                 repos = "http://packages.ropensci.org",
#                 type = "source")
library(hrbrthemes)
library(ggplot2)
library(maps)
library(colorBlindness)
library(ggthemes)
library(RColorBrewer)
library(cowplot)
library(extrafont)
font_import()
loadfonts(device = "pdf")
library(here)

# Script to run Figure 3 and base maps for Figure 1A and 4A-B
# Majambo Gamoyo and Luisa Fontoura 
# March 2021

#load data required for plotting figures
ID<-read.csv(here("_data","Data_Global_Connectivity","Data_globalFigures.csv"),h=T, stringsAsFactors = F,dec=".")
ID$LonTeste<-ifelse(ID$Lon <0, ID$Lon+360,ID$Lon)
#load all.data.Rds
all.data <- readRDS(here("_data","Connectivity_Biomass_SEMGLMMDATA_March2021.Rds"))
all.data$LonTeste<-ifelse(all.data$Lon <0, all.data$Lon+360,all.data$Lon)
#Calculate netflow for all points
dataNet <- ID
dataNet$crypto5Netflow <- (dataNet$crypto5BROF - dataNet$crypto5BRIF)/(dataNet$crypto5BROF + dataNet$crypto5BRIF)


# rm NAs for netflow for plotting purposes
dataNet <- dataNet %>% filter(!is.na(crypto5Netflow))
#create color for top 10% sinks and sources
dataNet$colorSinkSource<- ifelse(dataNet$crypto5Netflow > quantile(dataNet$crypto5Netflow,0.9),"darkorange",
                                 ifelse(dataNet$crypto5Netflow < quantile(dataNet$crypto5Netflow,0.1),"navyblue","gray25"))
SinkSourceData<-dataNet %>% filter(colorSinkSource %in% c("darkorange","navyblue"))
#head(SinkSourceData)
colourSinks<-as.character(SinkSourceData$colorSinkSource)

#ggplot(dataNet) + 
#  geom_point(aes(x=(crypto5BROF),y=(crypto5BRIF)),color=colourSinks,alpha=0.5)

#Global_data <- read.csv("Global_updated_data.csv",header = T)
#Set World Map
world <- ne_coastline(scale = "large", returnclass = "sf")
world2 = maps::map(wrap=c(0,360), proj="mollweide",plot=FALSE, fill=TRUE,interior = FALSE)
class(world)

##Fig 1A Map####----
theme_set(theme_classic(base_size = 22))
t <- ggplot() +
  geom_sf() +
  borders("world2",fill="gray90", col="gray90", bg="lightblue") +
  #geom_polygon(data = fortify(maps::map("world2",plot=FALSE,fill=TRUE)), aes(x=long, y = lat, group=group)) +
  geom_point(data = ID,
             aes(x = LonTeste, y = Lat), color="gray15",fill="gray55",
             alpha = 0.5, size=0.75) +
  theme_bw() +
  #scale_colour_gradientn(colours = rev(c("DarkRed","orange","white","cyan3","Blue2"))) +
  #scale_fill_discrete(values=Blue2DarkRed12Steps) +
  theme(legend.title = element_blank(),legend.key.height= unit(.5, 'cm'),
        legend.key.width= unit(.5, 'cm')) +
  labs(title = "") +
  labs(x = "Longitude") +
  labs(y = "Latitude") +
  theme(text = element_text(size = 14, family="Arial")) +
  coord_sf(xlim = c(20, 340), ylim = c(-40, 40), expand = FALSE) +
  scale_y_continuous(breaks = seq(-40, 40, by=20))
#Fig 1
fin<-t +
  geom_point(data = all.data, 
             aes(x = LonTeste, y = Lat), color="grey15",fill="steelblue2",
             alpha = 0.5, size=2, shape=21)

fin


##Create transparent
fin+theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA)) + 
  scale_y_continuous(breaks = seq(-40, 40, by = 20))


##################################
################Fig 4$############
##################################
#Netflow # Top 10% sinks and 10% sources dataset (SinkSourceData)
NetTop<-t +
  geom_point(data = SinkSourceData, 
             aes(x = LonTeste, y = Lat, group = colorSinkSource),
             color=SinkSourceData$colorSinkSource,
             alpha = 0.8, size=1, shape=19)

##Create transparent
NetTop+theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA)) + 
  scale_y_continuous(breaks = seq(-40, 40, by = 20))

##Indegree of upstream reefs (outdegree >0) - top dispersal corridors for each biogeographical regions
quantile(dataNet[dataNet$Kulbicki == "Western Indian",]$pare5BRin,0.9) #62
quantile(dataNet[dataNet$Kulbicki == "Western Atlantic",]$pare5BRin,0.9) #60
quantile(dataNet[dataNet$Kulbicki == "Central Indo_Pacific",]$pare5BRin,0.9) #118
quantile(dataNet[dataNet$Kulbicki == "Central Pacific",]$pare5BRin,0.9) #96

InDTop<-t +
  geom_point(data = dataNet[dataNet$Kulbicki == "Western Indian" &
                              dataNet$pare5BRin>61,], 
             aes(x = LonTeste, y = Lat),
             color="darkred",
            alpha=0.45,size=1, shape=19) + 
  geom_point(data = dataNet[dataNet$Kulbicki == "Western Atlantic" &
                              dataNet$pare5BRin>59,], 
             aes(x = LonTeste, y = Lat),
             color="darkred",
             alpha=0.45,size=1, shape=19) +
  geom_point(data = dataNet[dataNet$Kulbicki == "Central Indo_Pacific" &
                              dataNet$pare5BRin>117,], 
             aes(x = LonTeste, y = Lat),
             color="darkred",
             alpha=0.45,size=1, shape=19) +
  geom_point(data = dataNet[dataNet$Kulbicki == "Central Pacific" &
                             dataNet$pare5BRin>95,], 
            aes(x = LonTeste, y = Lat),
            color="darkred",
            alpha=0.45,size=1, shape=19)
  

##Create transparent
InDTop+theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA)) + 
  scale_y_continuous(breaks = seq(-40, 40, by = 20))

#ggsave("Fig1A_basemap.png",type = "cairo-png",bg = "transparent",dpi = 1000,width=10,height=5)
#ggsave("Netflow_crypto5_basemap.png",type = "cairo-png",bg = "transparent",dpi = 1000,width=10,height=5)
#ggsave("IndegreeSources_pare5_basemap.png",type = "cairo-png",bg = "transparent",dpi = 1000,width=10,height=5)
#ggsave("netflowmap.pdf",dpi = 600,width=12,height=7)
#ggsave("nflowmap.png",type = "cairo-png",dpi = 1000,width=10,height=5)


##############
######Fig 3###
##############
head(dataNet)
#filter 

ggplot() +geom_point(data = dataNet,
                     aes(y=crypto5Netflow, x=log(pare5BRin), color=General), alpha=0.5) +
  scale_colour_manual(values = c("seashell","#1F968BFF")) + theme_minimal()


##The stat for adding vetical lines at pecentile locations 
StatPercentileX <- ggproto("StatPercentileX", Stat,
                           compute_group = function(data, scales, probs) {
                             percentiles <- quantile(data$x, probs=probs)
                             data.frame(xintercept=percentiles)
                           },
                           required_aes = c("x")
)

stat_percentile_x <- function(mapping = NULL, data = NULL, geom = "vline",
                              position = "identity", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    stat = StatPercentileX, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

StatPercentileXLabels <- ggproto("StatPercentileXLabels", Stat,
                                 compute_group = function(data, scales, probs) {
                                   percentiles <- quantile(data$x, probs=probs)
                                   data.frame(x=percentiles, y=Inf,
                                              label=paste0("p", probs*100, ": ",
                                                           round(percentiles, digits=3)))
                                 },
                                 required_aes = c("x")
)

stat_percentile_xlab <- function(mapping = NULL, data = NULL, geom = "text",
                                 position = "identity", na.rm = FALSE,
                                 show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    stat = StatPercentileXLabels, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
##The stat for adding vetical lines at pecentile locations 
StatPercentileY <- ggproto("StatPercentileY", Stat,
                           compute_group = function(data, scales, probs) {
                             percentiles <- quantile(data$y, probs=probs)
                             data.frame(yintercept=percentiles)
                           },
                           required_aes = c("y")
)

stat_percentile_y <- function(mapping = NULL, data = NULL, geom = "hline",
                              position = "identity", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    stat = StatPercentileY, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

StatPercentileYLabels <- ggproto("StatPercentileXLabels", Stat,
                                 compute_group = function(data, scales, probs) {
                                   percentiles <- quantile(data$y, probs=probs)
                                   data.frame(y=percentiles, y=Inf,
                                              label=paste0("p", probs*100, ": ",
                                                           round(percentiles, digits=3)))
                                 },
                                 required_aes = c("y")
)

stat_percentile_ylab <- function(mapping = NULL, data = NULL, geom = "text",
                                 position = "identity", na.rm = FALSE,
                                 show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    stat = StatPercentileXLabels, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


#range01 <- function(x){(x-min(x))/(max(x)-min(x))}
#dataNet$crypto5NetflowNOR<-range01(dataNet$crypto5Netflow)
#dataNet$pare5BRinNOR<-range01(dataNet$pare5BRin)
#quantile(dataNet$pare5BRin,probs=seq(0,1,0.10))
library(stats)
percentile <- ecdf(dataNet$crypto5Netflow)
dataNet$crypto5Netflow_percentile<-percentile(dataNet$crypto5Netflow)
percentile <- ecdf(dataNet$pare5BRin)
dataNet$pare5BRin_percentile<-percentile(dataNet$pare5BRin)
head(dataNet)

allRe<- ggplot(dataNet, aes(y=crypto5Netflow_percentile, x=pare5BRin_percentile)) +
  geom_point(alpha=0) + 
  geom_point(data=dataNet %>% filter(General == "MPA"), 
             aes(y=crypto5Netflow_percentile, x=pare5BRin_percentile),alpha=0.25, color="#1F968BFF") +theme_bw() +
  stat_percentile_x(probs=c(0.9), linetype=2) +
  stat_percentile_xlab(probs=c(0.9), hjust=1, vjust=1.5, angle=180) +
  stat_percentile_y(probs=c(0.1,0.90), linetype=2) +
  stat_percentile_ylab(probs=c(0.1,0.90), hjust=1, vjust=1.5, angle=180)

#marginal densities along y axis
xdens <- axis_canvas(allRe, axis = "x")+
  geom_density(data = dataNet %>% filter(General == "MPA"), aes(x = pare5BRin, fill = as.factor(General)),
               alpha = 0.7, size = 0.2)+
  scale_fill_manual(values = c("#1F968BFF")) 
#Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(allRe, axis = "y", coord_flip = TRUE)+
  geom_density(data = dataNet %>% filter(General == "MPA"), aes(x = crypto5Netflow, fill = as.factor(General)),
               alpha = 0.7, size = 0.2)+ coord_flip()+
  scale_fill_manual(values = c("#1F968BFF")) 

p1 <- insert_xaxis_grob(allRe, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
MPa_plot<-ggdraw(p2)

#Fished Areas
allReFi<- ggplot(dataNet, aes(y=crypto5Netflow_percentile, x=pare5BRin_percentile)) +
  geom_point(alpha=0) + 
  geom_point(data=dataNet %>% filter(General == "Fished"), 
             aes(y=crypto5Netflow_percentile, x=pare5BRin_percentile),alpha=0.25, color="gray50") +theme_bw() +
  stat_percentile_x(probs=c(0.9), linetype=2) +
  stat_percentile_xlab(probs=c(0.9), hjust=1, vjust=1.5, angle=180) +
  stat_percentile_y(probs=c(0.1,0.90), linetype=2) +
  stat_percentile_ylab(probs=c(0.1,0.90), hjust=1, vjust=1.5, angle=180)

#marginal densities along y axis
xdens <- axis_canvas(allReFi, axis = "x")+
  geom_density(data = dataNet %>% filter(General == "Fished"), aes(x = pare5BRin, fill = as.factor(General)),
               alpha = 0.7, size = 0.2)+
  scale_fill_manual(values = c("gray50")) 
#Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(allReFi, axis = "y", coord_flip = TRUE)+
  geom_density(data = dataNet %>% filter(General == "Fished"), aes(x = crypto5Netflow, fill = as.factor(General)),
               alpha = 0.7, size = 0.2)+ coord_flip()+
  scale_fill_manual(values = c("gray50")) 

p3 <- insert_xaxis_grob(allReFi, xdens, grid::unit(.2, "null"), position = "top")
p4<- insert_yaxis_grob(p3, ydens, grid::unit(.2, "null"), position = "right")
fished_plot<-ggdraw(p4)

singleplot<- ggplot(dataNet, aes(y=crypto5Netflow_percentile, x=pare5BRin_percentile)) +
  geom_point(alpha=0) + theme_bw() +
  stat_percentile_x(probs=c(0.9), linetype=2) +
  stat_percentile_xlab(probs=c(0.9), hjust=1, vjust=1.5, angle=180) +
  stat_percentile_y(probs=c(0.1,0.90), linetype=2) +
  stat_percentile_ylab(probs=c(0.1,0.90), hjust=1, vjust=1.5, angle=180)
p5 <- insert_xaxis_grob(singleplot, xdens, grid::unit(.2, "null"), position = "top")
p6<- insert_yaxis_grob(p5, ydens, grid::unit(.2, "null"), position = "right")
singleplotF<-ggdraw(p6)
  
library(ggpubr)
relatinnetfin <- ggarrange(singleplotF,fished_plot,
                           MPa_plot, nrow=1,common.legend = T,legend="right",labels = c("A","B","C"))

relatinnetfin_RV <- ggarrange(singleplot,allReFi,
                          allRe, nrow=1,common.legend = T,legend="right",labels = c("A","B","C"))

