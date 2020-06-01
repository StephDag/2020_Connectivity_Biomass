# SEM analysis - Global connectivity 
# author: Steph D'agata
# date: April 2020
# outputs: SEM coefficients

rm(list=ls())

# brms
install.packages("brms")
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("paul-buerkner/brms")
require(brms)
# stan
remove.packages("rstan")
if (file.exists(".RData")) file.remove(".RData")
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)

# packages
require(dplyr)
require(here)
require(forcats)
require(brms)
require(stan)
# load data
rm(all.data)
all.data<-read.csv(here("_data","FullDataMay2020Coordinates.csv"),h=T, stringsAsFactors = F,dec=".")
# clean first column
all.data$X.1 <- NULL
all.data$X <- NULL

# check 
head(all.data)
summary(all.data)

summary(all.data)
str(all.data)
apply(all.data,2,class)
# filter data with active

dim(data.std)
summary(data.std)
colnames(all.data)
head(all.data)
str(all.data)

# list of predictors
rm(PredictVar)
PredictVar<-all.data[,c("region","locality","sites","Richness","biomassarea1","temp","grav_total",
                        "Age_of_protection","Class","Indegree","ModelMode",         
                        "btwdegree","Inflow","Outdegree","InflowLR",          
                        "SelfR","Larval_behaviour","FE","InflowBR",          
                        "IndegreeBR","CorridorIndegreeBR","grav_neiBR","IndegreeMPABR",     
                        "InflowMPABR","IndegreeNeiBR","InflowNeiBR","Lon","Lat")]
PredictVar$log_grav_total <- log(PredictVar$grav_total+1)
PredictVar$log_grav_neiBR <- log(PredictVar$grav_neiBR+1)
PredictVar$log_biomassarea1<-log(PredictVar$biomassarea1+1)

# chage to factor
PredictVar$region <- as.factor(PredictVar$region)
PredictVar$locality <- as.factor(PredictVar$locality)
PredictVar$sites <- as.factor(PredictVar$sites)
PredictVar$Class <- as.factor(PredictVar$Class)
PredictVar$ModelMode <- as.factor(PredictVar$ModelMode)
PredictVar$Larval_behaviour <- as.factor(PredictVar$Larval_behaviour)
PredictVar$FE <- as.factor(PredictVar$FE)

# Transient
TRANSIENT %>% rm()
TRANSIENT <- PredictVar %>% filter(Larval_behaviour == "active" & ModelMode == "transi15") %>% droplevels()

# Fished as the reference
TRANSIENT$Class <- relevel(TRANSIENT$Class, ref="Fished")
summary(TRANSIENT)

library(corrgram)
corrgram(TRANSIENT,order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt)

## standrdize x variables
rm(TRANSIENT.std)
TRANSIENT.std<-data.frame(apply(X = TRANSIENT[,c(4,5,6,7,8,10,12:16,19:31)], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}))
TRANSIENT.std <- cbind(TRANSIENT$region,TRANSIENT.std)
TRANSIENT.std <- cbind(TRANSIENT$Class,TRANSIENT.std)
colnames(TRANSIENT.std)[1:2] <- c("Class","region")
# add log biomass
TRANSIENT.std$log_biomassarea1<-log(TRANSIENT.std$biomassarea1+1)
TRANSIENT.std$log_InflowLR<-log(TRANSIENT.std$InflowLR+1)
TRANSIENT.std$log_IndegreeBR<-log(TRANSIENT.std$IndegreeBR+1)
TRANSIENT.std$log_SelfR<-log(TRANSIENT.std$SelfR+1)
TRANSIENT.std$log_InflowMPABR<-log(TRANSIENT.std$InflowMPABR+1)
TRANSIENT.std$log_IndegreeMPABR<-log(TRANSIENT.std$IndegreeMPABR+1)
TRANSIENT.std$log_InflowNeiBR<-log(TRANSIENT.std$InflowNeiBR+1)
TRANSIENT.std$log_IndegreeNeiBR<-log(TRANSIENT.std$IndegreeNeiBR+1)
TRANSIENT.std$log_CorridorIndegreeBR <-log(TRANSIENT.std$CorridorIndegreeBR+1)

head(TRANSIENT.std)
# built SEM model
# scaled 
  # gravity of neightbour is highly positively correlated with connectivity variables
rm(species_mod_inflow)
species_mod_inflow <- bf(Richness ~ temp + log_grav_total + Class  + 
                    log_CorridorIndegreeBR + log_InflowLR +log_SelfR + 
                      log_InflowMPABR + log_InflowNeiBR + (1+ log_grav_total +log_SelfR |region))

rm(biom_mod_inflow)
biom_mod_inflow <- bf(log_biomassarea1 ~ Richness+temp + log_grav_total + Class  + 
                        log_CorridorIndegreeBR + log_InflowLR +log_SelfR + 
                        log_InflowMPABR + log_InflowNeiBR + 
                        (1+ log_grav_total + log_InflowLR + log_InflowMPABR + log_InflowNeiBR + log_SelfR |region))

# with lat/long
rm(species_mod_inflow_lat)
species_mod_inflow_lat <- bf(Richness ~ temp + log_grav_total+ Class  + 
                           log_CorridorIndegreeBR + log_InflowLR +log_SelfR + 
                           log_InflowMPABR + log_InflowNeiBR + Lat + Lon + (1+ log_grav_total +log_SelfR + Lat + Lon |region))

rm(biom_mod_inflow_lat)
biom_mod_inflow_lat <- bf(log_biomassarea1 ~ Richness+temp + log_grav_total+ Class  + 
                        log_CorridorIndegreeBR + log_InflowLR +log_SelfR +  
                        log_InflowMPABR + log_InflowNeiBR + Lat + Lon +
                        (1+ log_grav_total + log_InflowLR + log_InflowMPABR + log_InflowNeiBR + log_SelfR + Lat + Lon |region))


#rm(species_mod_simp)
#species_mod_simp <- bf(Richness ~ temp + log_grav_total+Class + (1+ log_grav_total |region)) 

#rm(biom_mod_simp)
#biom_mod_simp <- bf(log_biomassarea1 ~ Richness+temp + log_grav_total+Class+ (1+ log_grav_total + Richness |region))

# SEM
## scaled
    # all_fit_brms = connectivity explains both biomass and species richness
all_fit_brms %>% rm()
all_fit_brms <-brm(species_mod_inflow + biom_mod_inflow + set_rescor(FALSE), data=TRANSIENT.std,cores=4, chains = 4,
                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                   prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

all_fit_brms_lat %>% rm()
all_fit_brms_lat <-brm(species_mod_inflow_lat + biom_mod_inflow_lat + set_rescor(FALSE), data=TRANSIENT.std,cores=4, chains = 4,
                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                   prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

pairs(all_fit_brms)
# with no connectivity
all_fit_brms_simp %>% rm()
all_fit_brms_simp <-brm(species_mod_simp + biom_mod_simp + set_rescor(FALSE), data=TRANSIENT.std,cores=4, chains = 4,
                      iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 1),
                      prior = c(prior(normal(0, 10),class = "Intercept"), prior(normal(0, 10), class = "b")))

# with no connectivity
all_fit_brms_NR %>% rm()
all_fit_brms_NR <-brm(biom_mod_NR + set_rescor(FALSE), data=data.std,cores=4, chains = 4,
                        iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.990),
                        prior = c(prior(normal(0, 10),class = "Intercept"), prior(normal(0, 10), class = "b")))


# check the difference between the simple and full model 
LOO(all_fit_brms,all_fit_brms_simp)

# check divergene
pairs(all_fit_brms)
    # the full model is the best model

# summary of the full model
summary(all_fit_brms)
plot(all_fit_brms)
plot(all_fit_brms_FE)
pp_check(all_fit_brms, resp="logbiomassarea1")
pp_check(all_fit_brms, resp="Richness")
bayes_R2(all_fit_brms)
bayes_R2(all_fit_brms_simp)

library(tidyverse)
ggsave(plot=mcmc_plot(all_fit_brms),here("_prelim.figures","SEM_bays_extract.pdf"),width=20,height=30)

#
library("bayesplot")
library("ggplot2")
library("rstanarm")   

posterior <- as.array(all_fit_brms)
dimnames(posterior)

posterior$parameters
dim(posterior)
head(posterior)

Richness_inflow <- mcmc_intervals(posterior, pars = c("b_Richness_Intercept","b_Richness_temp","b_Richness_log_grav_total",                                                 
                                   "b_Richness_ClassClosed","b_Richness_ClassRestricted",                                                  
                                   "b_Richness_log_CorridorIndegreeBR","b_Richness_log_InflowLR",                                                    
                                   "b_Richness_log_SelfR","b_Richness_log_InflowMPABR","b_Richness_log_InflowNeiBR"))
Biomass_inflow <- mcmc_intervals(posterior,c("b_logbiomassarea1_Intercept","b_logbiomassarea1_Richness","b_logbiomassarea1_temp",                                                      
                          "b_logbiomassarea1_log_grav_total","b_logbiomassarea1_ClassClosed",                                               
                          "b_logbiomassarea1_ClassRestricted","b_logbiomassarea1_log_CorridorIndegreeBR",                                    
                          "b_logbiomassarea1_log_InflowLR","b_logbiomassarea1_log_SelfR",                                                 
                          "b_logbiomassarea1_log_InflowMPABR","b_logbiomassarea1_log_InflowNeiBR"))
Biomass
sjstats::rope(all_fit_brms, rope = c(-1, 1))
equi_test(all_fit_brms, out = "plot")
hist(sqrt(data.std$Inflow))

T_cub = sign(data.std$Inflow) * abs(data.std$Inflow)^(1/3)
library(rcompanion)
plotNormalHistogram(T_cub)
# split
require(caret)
require(klaR)
set.seed(0123)
split=0.90
rm(trainIndex); trainIndex <- createDataPartition(data.std$Richness, p=split, list=FALSE)
rm(data_train); data_train <- data.std[ trainIndex,]
rm(data_test); data_test <- data.std[-trainIndex,]


