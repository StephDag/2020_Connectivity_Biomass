# Preliminary analysis - Global connectivity 
# author: Steph D'agata
# date: February 2020
# outpupts: preliminary plots

rm(list=ls())

# load packages
library(lme4)
library(mgcv)
library(nlme)
library(lmerTest)
library(here)
library(dplyr)
library(ggplot2)
library(GGally)
library(tidyr)
library(stringr)

# load data
dataBIC %>% rm()
dataBIC<-read.csv(here("_data","Fullmatrixdatabiomassformodels.csv"),h=T)
dataBIC %>% head()
dataBIC %>% summary()
dataBIC %>% dim()
dataBIC %>% class()

# create unique ID for the full matrix

# correlation plot

  
crypto %>% rm()
crypto <- dataBIC %>% 
         filter(str_detect(ModelMode, c("crypto15","crypto5","crypto35"))) %>%
        dplyr::select(ID,region,locality,sites,ModelMode,Indegree) %>%
         spread(key=ModelMode,value=Indegree)
crypto <- crypto[-which(is.na(crypto$`crypto15 `)),]

          ggpairs(columns=18:20)
        
crypto %>% summary()
crypto %>% head()
crypto %>% dim()
crypto$ID %>% unique() %>% length()

crypto[which(crypto$ID == "V1272"),]
       
   

  
  require(datasets)
data("swiss")
              require(GGally)
              require(ggplot2)
              
              my_fn <- function(data, mapping, ...){
              p <- ggplot(data = data, mapping = mapping) + 
              geom_point() + 
              geom_smooth(method=loess, fill="red", color="red", ...) +
              geom_smooth(method=lm, fill="blue", color="blue", ...)
              p
              }
              
              g = ggpairs(swiss,columns = 1:4, lower = list(continuous = my_fn))
              g
  
        ggpairs(Indegree, lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.1)))



#%>%
#         spread(key=ModelMode, value=Indegree) %>% 
crypto %>% head()
crypto %>% summary()
crypto %>% dim()
crypto %>% names()
%>%
 %>%
 head()

dataBIC$ID %>% length()
dataBIC$ID %>% unique() %>% length()

  ggplot(color=Larval_behaviour)+
          geom_point() +
          geom_smooth()

ggcorr(nba[, -1],
       label = TRUE,
       label_alpha = TRUE,
       name = "") +
  ggplot2::theme(legend.position = "bottom")  

plot(log1p(dataBIC$grav_total),dataBIC$Indegree)
plot(log1p(dataBIC$grav_total),dataBIC$LocalRet)
plot(log1p(dataBIC$grav_total),dataBIC$Outdegree)


