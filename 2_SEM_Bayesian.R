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

require(dplyr)
require(here)
require(forcats)
require(brms)

# load data
rm(all.data)
all.data<-read.csv(here("_data","fulldatabaseupdated2403.csv"),h=T)

colnames(all.data)

rm(PredictVar)
PredictVar<-all.data[,c("region","locality","sites","Richness","temp","grav_total","Age_of_protection","Indegree","btwdegree","Inflow","Outdegree","InflowLR","SelfR", "Class")]

##standrdize x variables
rm(data.std)
data.std<-data.frame(apply(X = PredictVar[,4:13], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}))

data.std$region<-as.factor(all.data$region)
data.std$locality<-as.factor(all.data$locality)
data.std$sites<-as.factor(all.data$sites)
data.std$ModelMode<-as.factor(all.data$ModelMode)
data.std$ModelMode<-as.factor(all.data$ModelMode)
data.std$Class<-as.factor(all.data$Class)
data.std$log_biomassarea1<-log(all.data$biomassarea1)
data.std$Larval_behaviour<-as.factor(all.data$Larval_behaviour)
data.std$FE<-as.factor(all.data$FE)
data.std$ModelMode<-as.factor(all.data$ModelMode)

# Fished as the reference
data.std$Class <- relevel(data.std$Class, ref="Fished")
data.std$Larval_behaviour <- relevel(data.std$Larval_behaviour, ref="passive")
data.std$FE <- relevel(data.std$FE, ref="crypto")

head(data.std)
names(data.std)
summary(data.std)

# filter data with active
data.std <- data.std %>% filter(Larval_behaviour == "active")
dim(data.std)
summary(data.std)

plot(data.std$Richness,data.std$log_biomassarea1)
test <- lm(log_biomassarea1 ~ poly(Richness,3),data=data.std)
# check simple relationships
test <- lm(Richness ~ poly(temp,2),data=data.std)
visreg::visreg(test)
test <- lm(log_biomassarea1 ~ poly(temp,2),data=data.std)
visreg::visreg(test)
test <- lm(Richness ~ poly(Age_of_protection,2),data=data.std)
visreg::visreg(test)
test <- lm(log_biomassarea1 ~ poly(Age_of_protection,2),data=data.std)
visreg::visreg(test)
test <- lm(Richness ~ sqrt(btwdegree),data=data.std)
plot(sqrt(data.std$btwdegree),data.std$Richness)
visreg::visreg(test)
test <- lm(log_biomassarea1 ~ Inflow,data=data.std)
visreg::visreg(test)

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


# built SEM model
# scaled

rm(species_mod)
species_mod <- bf(Richness ~ poly(temp,3) + grav_total*Class + poly(Age_of_protection,3) + 
                    Indegree+btwdegree+Inflow +SelfR+ModelMode/FE+
                    (1 |sites/locality/region))
#                   (1 |ModelMode/Larval_behaviour/FE))
rm(biom_mod)
biom_mod <- bf(log_biomassarea1 ~ poly(Richness,3)+poly(temp,3) + grav_total*Class+poly(Age_of_protection,3) +
                 Indegree+btwdegree+Inflow+ SelfR+ModelMode/FE+
                 (1 |sites/locality/region))
#                 (1 |ModelMode/Larval_behaviour/FE))
rm(biom_mod_simp)
biom_mod_simp <- bf(log_biomassarea1 ~ poly(Richness,3)+poly(temp,3) + grav_total*Class+poly(Age_of_protection,3)+
                      ModelMode/FE+
                      (1 |sites/locality/region))

random <- bf(log_biomassarea1 ~ (1 |sites/locality/region))
# SEM
## scaled
    # all_fit_brms = connectivity explains both biomass and species richness
all_fit_brms %>% rm()
all_fit_brms <-brm(species_mod + biom_mod + set_rescor(FALSE), data=data.std,cores=4, chains = 3,
                   iter = 2000, warmup = 200,thin = 2, refresh = 0, control = list(adapt_delta = 0.95),
                   prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

    # all_fit_brms_simp = connectivity is mediated through species richness
all_fit_brms_simp <-brm(species_mod + biom_mod_simp + set_rescor(FALSE), data=data.std,cores=4, chains = 3,
                        iter = 2000, warmup = 200,thin = 2, refresh = 0, control = list(adapt_delta = 0.95),
                        prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

# check the difference between the simple and full model 
LOO(all_fit_brms,all_fit_brms_simp)
    # the full model is the best model

# summary of the full model
summary(all_fit_brms)
plot(all_fit_brms_r)
pp_check(all_fit_brms, resp="logbiomassarea1")
pp_check(all_fit_brms, resp="Richness")
bayes_R2(all_fit_brms)
bayes_R2(all_fit_brms_simp)

library(tidyverse)
ggsave(plot=mcmc_plot(all_fit_brms),here("_prelim.figures","SEM_bays_extract.pdf"),width=20,height=30)

## nice summary
library(sjPlot)
library(insight)
library(httr)

tab_model(all_fit_brms)
