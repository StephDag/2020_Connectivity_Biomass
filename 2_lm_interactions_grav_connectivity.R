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
dataBIC$logB <- log1p(dataBIC$biomassarea)
dataBIC$logB1 <- log1p(dataBIC$biomassarea1)
dataBIC$logB2 <- log1p(dataBIC$biomassarea2)

# explanatory variables and connectivity (LR)
all %>% rm()
all <- dataBIC %>% 
  mutate(ModelMode = as.character(ModelMode)) %>%
  dplyr::filter(ModelMode %in% c("crypto15","pare15","transi15","resid15")) %>%
  dplyr::filter(Larval_behaviour %in% c("active")) %>%
 # dplyr::filter(Age_MPA > 0) %>%
  dplyr::group_by(ModelMode) %>% 
  dplyr::select(sites,
                ModelMode,
                Indegree,logB,logB1,logB2,logGrav,Age_MPA,Class) %>% 
  mutate(grouped_id = row_number()) %>% 
  spread(key=ModelMode,value=Indegree) %>%
  as.data.frame()
dim(all)
head(all)
summary(all)

#all$quantresid <- ifelse(all$resid15 < 17, "LOW", ifelse(
#                        all$resid15 >= 17 & all$resid15<80, "MEDIUM.low",ifelse(
#                        all$resid15>= 80 & all$resid15 < 200, 'MEDIUM.high',"HIGH")))

#all$quantresid <- ifelse(all$resid15 < 17, "LOW", ifelse(
#  all$resid15 >= 17 & all$resid15<80, "MEDIUM","HIGH"))

all$quantresid <- ifelse(all$resid15 < 20, "LOW","HIGH")
all$quantresid <- as.factor(all$quantresid)
summary(as.factor(all$quantresid))
quantile(all$resid15, prob = seq(0,1,0.01))
hist(all$resid15,breaks=100)

# linear model
test %>% rm()
test <- lm(logB ~  quantresid*logGrav+ quantresid*Age_MPA,data=all)
summary(test)
anova(test)

library(visreg)
visreg(test)
visreg2d(test,xvar="logGrav",yvar="Age_MPA")
visreg::visreg(test,xvar="logGrav",by="quantresid")
visreg::visreg(test,xvar="logGrav",by="Class")
visreg::visreg(test,xvar="Age_MPA",by="quantresid")


test.2 %>% rm()
test.2 <- lm(logB ~  logGrav+Age_MPA,data=all)
summary(test.2)
anova(test.2)

AIC(test, test.2)

# gam
library(mgcv)
test.gam <- mgcv::gam(logB ~ quantresid+  s(logGrav,by = quantresid) + 
                        s(logGrav,by = Class),data=all,method="REML")
summary(test.gam)
anova(test.gam)

visreg(test.gam)
visreg::visreg(test.gam,xvar="logGrav",by="quantresid")
visreg::visreg(test.gam,xvar="logGrav",by="Class")
