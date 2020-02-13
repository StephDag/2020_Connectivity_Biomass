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
library(tidyverse)

# load data
dataBIC %>% rm()
dataBIC<-read.csv(here("_data","Fullmatrixdatabiomassformodels.csv"),h=T)
dataBIC %>% head()
dataBIC %>% summary()
dataBIC %>% dim()
dataBIC %>% class()
dataBIC %>% str()
dataBIC$ModelMode <- gsub(" ","",dataBIC$ModelMode)
dataBIC$logGrav <- log1p(dataBIC$grav_total)


#########################
#       IN DEGREE       #
#########################

# all
all %>% rm()
all <- dataBIC %>% 
  mutate(ModelMode = as.character(ModelMode)) %>%
  dplyr::filter(ModelMode %in% c("crypto5","crypto15","crypto35",
                                 "pare5","pare15","pare25",
                                 "transi5","transi15","transi35",
                                 "resid5","resid15","resid35")) %>%
  dplyr::filter(Larval_behaviour %in% c("active")) %>%
  dplyr::group_by(ModelMode) %>% 
  dplyr::select(sites,
                ModelMode,
                Indegree) %>% 
  mutate(grouped_id = row_number()) %>% 
  spread(key=ModelMode,value=Indegree) %>%
  as.data.frame()
dim(all)
head(all)
summary(all)

# fonction to plot correlation with lm and loess function
              my_fn <- function(data, mapping, ...){
              p <- ggplot(data = data, mapping = mapping) + 
              geom_point() + 
              geom_smooth(method=loess, fill="red", color="red", ...) +
              geom_smooth(method=lm, fill="blue", color="blue", ...)
              p
              }
              
# correlation between all functional groups and  scnearios with biplots        
g = ggpairs(all,columns = 3:14, lower = list(continuous = my_fn))
g
ggsave(here("_prelim.figures","correl_scenario_INDEGREE.V1.pdf"),plot=g,width=20,height=20)
# correlation between all functional groups and  scnearios only R2 values 
library(ggcorrplot)
library(RColorBrewer)
correlation_matrix <- cor(all[,3:14])
colors <- brewer.pal(n = 3, name = "RdYlBu")
p <- ggcorrplot(correlation_matrix , type = "upper", hc.order = TRUE, colors = brewer.pal(n = 3, name = "RdYlBu"))
p <- p + scale_fill_gradient2(limit = c(0.5,1), low = "blue", high = "red", mid="orange", midpoint = 0.75)
p
ggsave(here("_prelim.figures","correl_scenario__INDEGREE.V2.pdf"),plot=p)

# explanatory variables and connectivity (in)
all %>% rm()
all <- dataBIC %>% 
  mutate(ModelMode = as.character(ModelMode)) %>%
  dplyr::filter(ModelMode %in% c("crypto15","pare15","transi15","resid15")) %>%
  dplyr::filter(Larval_behaviour %in% c("active")) %>%
  dplyr::group_by(ModelMode) %>% 
  dplyr::select(sites,
                ModelMode,
                Indegree,biomassarea,biomassarea1,biomassarea2,logGrav) %>% 
  mutate(grouped_id = row_number()) %>% 
  spread(key=ModelMode,value=Indegree) %>%
  as.data.frame()
dim(all)
head(all)
summary(all)

g %>% rm()
g = ggpairs(all, columns=c(2:5,7:10), lower = list(continuous = my_fn))
g
ggsave(here("_prelim.figures","correl_scenario_INDEGREE.expla.pdf"),plot=g,width=20,height=20)

# explanatory variables and connectivity (LR)
all %>% rm()
all <- dataBIC %>% 
  mutate(ModelMode = as.character(ModelMode)) %>%
  dplyr::filter(ModelMode %in% c("crypto15","pare15","transi15","resid15")) %>%
  dplyr::filter(Larval_behaviour %in% c("active")) %>%
  dplyr::group_by(ModelMode) %>% 
  dplyr::select(sites,
                ModelMode,
                LocalRet,biomassarea,biomassarea1,biomassarea2,logGrav) %>% 
  mutate(grouped_id = row_number()) %>% 
  spread(key=ModelMode,value=LocalRet) %>%
  as.data.frame()
dim(all)
head(all)
summary(all)

g %>% rm()
g = ggpairs(all, columns=c(2:5,7:10), lower = list(continuous = my_fn))
g
ggsave(here("_prelim.figures","correl_scenario_LR.expla.pdf"),plot=g,width=20,height=20)

