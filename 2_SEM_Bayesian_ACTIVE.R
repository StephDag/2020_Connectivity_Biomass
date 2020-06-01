# SEM analysis - Global connectivity 
# author: Steph D'agata
# date: April 2020
# updates: May 2020
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
require(rstan)
# load data
rm(all.data)
all.data<-read.csv(here("_data","FullDataMay2020Coordinates.csv"),h=T)
# clean first column
all.data$X.1 <- NULL
all.data$X <- NULL

# filter data with active
summary(all.data %>% filter(Larval_behaviour == "active" & ModelMode == "transi15"))


colnames(all.data)
head(all.data)

# list of predictors
rm(PredictVar)
PredictVar<-all.data[,c("region","locality","sites","Richness","biomassarea1","temp","grav_total",
                        "Age_of_protection","Class","Indegree","ModelMode",         
                        "btwdegree","Inflow","Outdegree","InflowLR",          
                        "SelfR","Larval_behaviour","FE","InflowBR",          
                        "IndegreeBR","CorridorIndegreeBR","grav_neiBR","IndegreeMPABR",     
                        "InflowMPABR","IndegreeNeiBR","InflowNeiBR")]
PredictVar$log_grav_total <- log(PredictVar$grav_total+1)
PredictVar$log_grav_neiBR <- log(PredictVar$grav_neiBR+1)
PredictVar$log_biomassarea1<-log(PredictVar$biomassarea1+1)


hist(sqrt(PredictVar$Richness))

transformed <- abs(PredictVar$Richness - mean(PredictVar$Richness))
shapiro.test(PredictVar$log_biomassarea1)

shapiro.test(PredictVar$biomassarea1)
hist(PredictVar$log_biomassarea1,breaks=100)
hist(PredictVar$Richness,breaks=100)

demo2 <- rnorm(40, mean = 50, sd = 10)
shapiro.test(demo2)

summary(PredictVar)
names(PredictVar)
##standrdize x variables
rm(data.std)
# standardize
data.std<-data.frame(apply(X = PredictVar[,c(4,5,7,9,11:15,18:20,22:27)], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}))
# add log biomass
data.std$log_biomassarea1<-log(all.data$biomassarea1+1)

library(corrgram)
corrgram(data.std,order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt)

#IndegreeBR - the number of inward connections from cells
#IndegreeMPABR
plot(data.std$IndegreeBR,data.std$IndegreeMPABR)

plot(data.std$IndegreeBR,data.std$log_biomassarea1)
plot(data.std$IndegreeMPABR,data.std$log_biomassarea1)
  # add factor variables
data.std$region<-as.factor(all.data$region)
data.std$locality<-as.factor(all.data$locality)
data.std$sites<-as.factor(all.data$sites)
data.std$ModelMode<-as.factor(all.data$ModelMode)
data.std$Class<-as.factor(all.data$Class)
data.std$Larval_behaviour<-as.factor(all.data$Larval_behaviour)
data.std$FE<-as.factor(all.data$FE)

# Fished as the reference
data.std$Class <- relevel(data.std$Class, ref="Fished")
data.std$Larval_behaviour <- relevel(data.std$Larval_behaviour, ref="passive")
data.std$FE <- relevel(data.std$FE, ref="crypto")

head(data.std)
names(data.std)
summary(data.std)

# filter data with active
data.std <- data.std %>% filter(Larval_behaviour == "active" & ModelMode == "transi15") %>% droplevels()
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
species_mod <- bf(Richness ~ temp + log_grav_total+Class  + 
                    CorridorIndegreeBR+InflowBR +SelfR+ InflowMPABR+log_grav_neiBR+(1 |region))

rm(biom_mod)
biom_mod <- bf(log_biomassarea1 ~ Richness+ temp + log_grav_total+Class  + 
                 CorridorIndegreeBR+InflowBR +SelfR+ InflowMPABR+log_grav_neiBR+(1 |region))


rm(species_mod_simp)
species_mod_simp <- bf(Richness ~ temp + log_grav_total+Class + (1 |region)) 

rm(biom_mod_simp)
biom_mod_simp <- bf(log_biomassarea1 ~ Richness+temp + log_grav_total+Class+ (1 |region))

# SEM
## scaled
    # all_fit_brms = connectivity explains both biomass and species richness
all_fit_brms %>% rm()
all_fit_brms <-brm(species_mod + biom_mod + set_rescor(FALSE), data=data.std,cores=4, chains = 4,
                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                   prior = c(prior(normal(0, 10),class = "Intercept"), prior(normal(0, 10), class = "b")))
saveRDS(all_fit_brms,here("_prelim.figures","Species_Biomass_connectivity.RDS"))
pairs(all_fit_brms)
# with no connectivity
all_fit_brms_simp %>% rm()
all_fit_brms_simp <-brm(species_mod_simp + biom_mod_simp + set_rescor(FALSE), data=data.std,cores=4, chains = 4,
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
brms::pp_check(all_fit_brms, resp="logbiomassarea1",nsamples=100)
brms::pp_check(all_fit_brms, resp="Richness",nsamples=100)
bayes_R2(all_fit_brms)
bayes_R2(all_fit_brms_simp)

library(tidyverse)
ggsave(plot=mcmc_plot(all_fit_brms),here("_prelim.figures","SEM_bays_extract.pdf"),width=20,height=30)

## nice summary
library(sjPlot)
library(insight)
library(httr)

tab_model(all_fit_brms)
