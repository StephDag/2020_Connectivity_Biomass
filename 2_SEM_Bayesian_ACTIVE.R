# SEM analysis - Global connectivity 
# author: Steph D'agata
# date: April 2020
# outputs: SEM coefficients

#rm(list=ls())

# brms

#if (!requireNamespace("remotes")) {
#  install.packages("remotes")
#}
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
#remotes::install_github("stan-dev/rstan", ref = "develop", subdir = "rstan/rstan", build_opts = "")
 
# stan
#remove.packages("rstan")
#if (file.exists(".RData")) file.remove(".RData")
#pkgbuild::has_build_tools(debug = TRUE)
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)

# packages
require(dplyr)
require(here)
require(forcats)
require(brms)
require(rstan)

# load data
rm(all.data)
all.data<-read.csv(here("_data","Connectivity_Biomass_SEMGLMMDATA.csv"),h=T, stringsAsFactors = F,dec=".")
# clean first column
all.data$X <- NULL

# check 
head(all.data)
summary(all.data)
dim(all.data)
summary(all.data)
str(all.data)
names(all.data)
apply(all.data,2,class)

dim(all.data)

# log all the data
all.data$log_grav_total <- log(all.data$grav_total+1)
all.data$log_grav_neiBR <- log(all.data$grav_nei+1)
all.data$log_biomassarea <-log(all.data$biomassarea1+1)

# chage to factor
all.data$region <- as.factor(all.data$region)
all.data$locality <- as.factor(all.data$locality)
all.data$sites <- as.factor(all.data$sites)
all.data$Class <- as.factor(all.data$Class)
all.data$ModelMode <- as.factor(all.data$ModelMode)
all.data$Larval_behaviour <- as.factor(all.data$Larval_behaviour)
all.data$FE <- as.factor(all.data$FE)

# Transform connection with MPAs
all.data$IndegreeMPABR.bin <- ifelse(all.data$IndegreeMPABR == 0,0,1) %>% as.factor()
all.data$IndegreeMPABR.bin %>% summary()
# Transient
TRANSIENT %>% rm()
TRANSIENT <- all.data %>% filter(Larval_behaviour == "active" & ModelMode == "transi15") %>% droplevels()

# Fished as the reference
TRANSIENT$Class <- relevel(TRANSIENT$Class, ref="Fished")
summary(TRANSIENT)

library(corrgram)
corrgram(TRANSIENT,order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt)

connectivity <- all.data[,c("SelfR","InflowBR","IndegreeBR",        
"CorridorIndegreeBR","IndegreeMPABR","InflowMPABR",       
"IndegreeNeiBR","InflowNeiBR","OutFlow","Outdegree","btwdegree")] 
corr.connectivity <- corrgram(connectivity,order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.cor, text.panel=panel.txt)

ggsave(here("_prelim.figures","Corr_connectivity.pdf"),corr.connectivity,width=20,height=12)

## standrdize x variables
rm(TRANSIENT.std)
TRANSIENT.std<-data.frame(apply(X = TRANSIENT[,c(5,6,12:18,19:23)], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))
#TRANSIENT.std<-data.frame(apply(X = TRANSIENT[,c(5,6,14:17,19:28)], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))
TRANSIENT.std <- cbind(TRANSIENT$region,TRANSIENT.std)
TRANSIENT.std <- cbind(TRANSIENT$Class,TRANSIENT.std)

colnames(TRANSIENT.std)[c(1,2)] <- c("Class","region")
names(TRANSIENT.std)# add log biomass

TRANSIENT.std$log_InflowBR<-log(TRANSIENT.std$Inflow+1)
TRANSIENT.std$log_IndegreeBR<-log(TRANSIENT.std$Indegree+1)
TRANSIENT.std$log_SelfR<-log(TRANSIENT.std$SelfR+1)
TRANSIENT.std$log_InflowMPABR<-log(TRANSIENT.std$InflowMPA+1)
TRANSIENT.std$log_IndegreeMPABR<-log(TRANSIENT.std$IndegreeMPA+1)
TRANSIENT.std$log_InflowNeiBR<-log(TRANSIENT.std$InflowNei+1)
TRANSIENT.std$log_CorridorIndegreeBR <-log(TRANSIENT.std$CorridorIndegree+1)
TRANSIENT.std$log_btwdegree <-log(TRANSIENT.std$btwdegree+1)
TRANSIENT.std$log_outdegree <-log(TRANSIENT.std$Outdegree+1)
TRANSIENT.std$Class <- relevel(TRANSIENT.std$Class, ref="Fished")
TRANSIENT.std$Netflow<- TRANSIENT$Netflow
TRANSIENT.std$log_biomassarea <- TRANSIENT$log_biomassarea
TRANSIENT.std$log_grav_total <- TRANSIENT$log_grav_total
TRANSIENT.std$log_grav_neiBR  <- TRANSIENT$log_grav_neiBR
head(TRANSIENT.std)
dim(TRANSIENT.std)
summary(TRANSIENT.std)

# Parental
PARENTAL %>% rm()
PARENTAL <- all.data %>% filter(Larval_behaviour == "active" & ModelMode == "pare5") %>% droplevels()

# Fished as the reference
PARENTAL$Class <- relevel(PARENTAL$Class, ref="Fished")
summary(PARENTAL.std)

## standrdize x variables
rm(PARENTAL.std)
PARENTAL.std<-data.frame(apply(X = PARENTAL[,c(5,6,12:18,19:23)], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))
#PARENTAL.std<-data.frame(apply(X = PARENTAL[,c(5,6,14:17,19:28)], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))
PARENTAL.std <- cbind(PARENTAL$region,PARENTAL.std)
PARENTAL.std <- cbind(PARENTAL$Class,PARENTAL.std)

colnames(PARENTAL.std)[c(1,2)] <- c("Class","region")
names(PARENTAL.std)# add log biomass

PARENTAL.std$log_InflowBR<-log(PARENTAL.std$Inflow+1)
PARENTAL.std$log_IndegreeBR<-log(PARENTAL.std$Indegree+1)
PARENTAL.std$log_SelfR<-log(PARENTAL.std$SelfR+1)
PARENTAL.std$log_InflowMPABR<-log(PARENTAL.std$InflowMPA+1)
PARENTAL.std$log_IndegreeMPABR<-log(PARENTAL.std$IndegreeMPA+1)
PARENTAL.std$log_InflowNeiBR<-log(PARENTAL.std$InflowNei+1)
PARENTAL.std$log_CorridorIndegreeBR <-log(PARENTAL.std$CorridorIndegree+1)
PARENTAL.std$log_btwdegree <-log(PARENTAL.std$btwdegree+1)
PARENTAL.std$log_outdegree <-log(PARENTAL.std$Outdegree+1)
PARENTAL.std$Class <- relevel(PARENTAL.std$Class, ref="Fished")
PARENTAL.std$Netflow<- PARENTAL$Netflow
PARENTAL.std$log_biomassarea <- PARENTAL$log_biomassarea
PARENTAL.std$log_grav_total <- PARENTAL$log_grav_total
PARENTAL.std$log_grav_neiBR  <- PARENTAL$log_grav_neiBR
head(PARENTAL.std)
dim(PARENTAL.std)
summary(PARENTAL.std)

# Cryptic
CRYPTIC %>% rm()
CRYPTIC <- all.data %>% filter(Larval_behaviour == "active" & ModelMode == "crypto5") %>% droplevels()

# Fished as the reference
CRYPTIC$Class <- relevel(CRYPTIC$Class, ref="Fished")
summary(CRYPTIC)


## standrdize x variables
rm(CRYPTIC.std)
CRYPTIC.std<-data.frame(apply(X = CRYPTIC[,c(5,6,12:18,19:23)], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))
#CRYPTIC.std<-data.frame(apply(X = CRYPTIC[,c(5,6,14:17,19:28)], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))
CRYPTIC.std <- cbind(CRYPTIC$region,CRYPTIC.std)
CRYPTIC.std <- cbind(CRYPTIC$Class,CRYPTIC.std)

colnames(CRYPTIC.std)[c(1,2)] <- c("Class","region")
names(CRYPTIC.std)# add log biomass

CRYPTIC.std$log_InflowBR<-log(CRYPTIC.std$Inflow+1)
CRYPTIC.std$log_IndegreeBR<-log(CRYPTIC.std$Indegree+1)
CRYPTIC.std$log_SelfR<-log(CRYPTIC.std$SelfR+1)
CRYPTIC.std$log_InflowMPABR<-log(CRYPTIC.std$InflowMPA+1)
CRYPTIC.std$log_IndegreeMPABR<-log(CRYPTIC.std$IndegreeMPA+1)
CRYPTIC.std$log_InflowNeiBR<-log(CRYPTIC.std$InflowNei+1)
CRYPTIC.std$log_CorridorIndegreeBR <-log(CRYPTIC.std$CorridorIndegree+1)
CRYPTIC.std$log_btwdegree <-log(CRYPTIC.std$btwdegree+1)
CRYPTIC.std$log_outdegree <-log(CRYPTIC.std$Outdegree+1)
CRYPTIC.std$Class <- relevel(CRYPTIC.std$Class, ref="Fished")
CRYPTIC.std$Netflow<- CRYPTIC$Netflow
CRYPTIC.std$log_biomassarea <- CRYPTIC$log_biomassarea
CRYPTIC.std$log_grav_total <- CRYPTIC$log_grav_total
CRYPTIC.std$log_grav_neiBR  <- CRYPTIC$log_grav_neiBR
head(CRYPTIC.std)
dim(CRYPTIC.std)
summary(CRYPTIC.std)

# Resident
RESID %>% rm()
RESID <- all.data %>% filter(Larval_behaviour == "active" & ModelMode == "resid15") %>% droplevels()

# Fished as the reference
RESID$Class <- relevel(RESID$Class, ref="Fished")
summary(RESID)

## standrdize x variables
rm(RESID.std)
RESID.std<-data.frame(apply(X = RESID[,c(5,6,12:18,19:23)], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))
#RESID.std<-data.frame(apply(X = RESID[,c(5,6,14:17,19:28)], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))
RESID.std <- cbind(RESID$region,RESID.std)
RESID.std <- cbind(RESID$Class,RESID.std)

colnames(RESID.std)[c(1,2)] <- c("Class","region")
names(RESID.std)# add log biomass

RESID.std$log_InflowBR<-log(RESID.std$Inflow+1)
RESID.std$log_IndegreeBR<-log(RESID.std$Indegree+1)
RESID.std$log_SelfR<-log(RESID.std$SelfR+1)
RESID.std$log_InflowMPABR<-log(RESID.std$InflowMPA+1)
RESID.std$log_IndegreeMPABR<-log(RESID.std$IndegreeMPA+1)
RESID.std$log_InflowNeiBR<-log(RESID.std$InflowNei+1)
RESID.std$log_CorridorIndegreeBR <-log(RESID.std$CorridorIndegree+1)
RESID.std$log_btwdegree <-log(RESID.std$btwdegree+1)
RESID.std$log_outdegree <-log(RESID.std$Outdegree+1)
RESID.std$Class <- relevel(RESID.std$Class, ref="Fished")
RESID.std$Netflow<- RESID$Netflow
RESID.std$log_biomassarea <- RESID$log_biomassarea
RESID.std$log_grav_total <- RESID$log_grav_total
RESID.std$log_grav_neiBR  <- RESID$log_grav_neiBR
head(RESID.std)
dim(RESID.std)
summary(RESID.std)

# built SEM model
# scaled 
  # gravity of neightbour is highly positively correlated with connectivity variables
rm(species_mod_inflow.intr) # ENV + CON
species_mod_inflow.intr <- bf(Richness ~ temp + log_grav_total*Class*Netflow  + 
                    log_InflowBR +log_SelfR + log_btwdegree +
                      (1+ log_grav_total + log_InflowBR +log_SelfR + log_btwdegree + Netflow |region))

rm(species_mod_inflow.intr.extr) # ENV + CON
species_mod_inflow.intr.extr <- bf(Richness ~ temp + log_grav_total*Class*Netflow+ 
                                log_InflowBR +log_SelfR + log_btwdegree + 
                                log_grav_neiBR + log_InflowMPABR + log_InflowNeiBR +
                                log_CorridorIndegreeBR + (1+ log_grav_total + log_InflowBR +log_SelfR + log_btwdegree + log_grav_neiBR + log_InflowMPABR + log_InflowNeiBR +
                                   log_CorridorIndegreeBR + Netflow |region))
rm(species_mod_inflow.intr.extr.simp) # ENV + CON
species_mod_inflow.intr.extr.simp <- bf(Richness ~ temp + log_grav_total*Class*Netflow+ 
                                     log_InflowBR +log_SelfR +  
                                     log_grav_neiBR +  log_InflowNeiBR  + 
                                     (1+ log_grav_total + log_InflowBR +log_SelfR +  log_grav_neiBR + log_InflowNeiBR  + Netflow |region))

rm(biom_mod_inflow.intr) # ENV + CON intr +S 
biom_mod_inflow.intr <- bf(log_biomassarea ~ Richness + temp + log_grav_total*Class*Netflow  + 
                             log_InflowBR +log_SelfR + log_btwdegree +  
                             (1+ log_grav_total + log_InflowBR +log_SelfR + log_btwdegree+ Netflow |region))

rm(biom_mod_inflow.intr.extr) # ENV + CON
biom_mod_inflow.intr.extr <- bf(log_biomassarea ~ Richness +temp + log_grav_total*Class*Netflow+ 
                                     log_InflowBR +log_SelfR + log_btwdegree + 
                                     log_grav_neiBR + log_InflowMPABR + log_InflowNeiBR +
                                     log_CorridorIndegreeBR + (1+ log_grav_total + log_InflowBR +log_SelfR + log_btwdegree + log_grav_neiBR + log_InflowMPABR + log_InflowNeiBR +
                                                                 log_CorridorIndegreeBR + Netflow |region))

rm(biom_mod_inflow.intr.extr.simp) # ENV + CON
biom_mod_inflow.intr.extr.simp <- bf(log_biomassarea ~ Richness + temp + log_grav_total*Class*Netflow+ 
                                       log_InflowBR +log_SelfR +  
                                       log_grav_neiBR +  log_InflowNeiBR  + 
                                       (1+ log_grav_total + log_InflowBR +log_SelfR +  log_grav_neiBR + log_InflowNeiBR  + Netflow |region))
rm(biom_mod_nocon.S) # ENV + S
biom_mod_nocon.S <- bf(log_biomassarea ~ Richness+temp + log_grav_total*Class  + (1+ log_grav_total|region))

rm(S_mod_nocon) # ENV
S_mod_nocon <- bf(Richness ~ temp + log_grav_total*Class  + (1+ log_grav_total|region))

# SEM
## scaled
    # all_fit_brms = connectivity explains both biomass and species richness
#TRANSIENT total
all_fit_brms.tot.TRANSIENT.intr %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.TRANSIENT.intr <-brm(species_mod_inflow.intr + biom_mod_inflow.intr + set_rescor(FALSE), data=TRANSIENT.std,cores=4,chains = 4,
                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                   prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.tot.TRANSIENT.intr,"all_fit_brms.tot.TRANSIENT.intr.Rds")
all_fit_brms.tot.TRANSIENT.intr.extr %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.TRANSIENT.intr.extr <-brm(species_mod_inflow.intr.extr + biom_mod_inflow.intr.extr + set_rescor(FALSE), data=TRANSIENT.std,cores=4,chains = 4,
                                      iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                      prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.tot.TRANSIENT.intr.extr,"all_fit_brms.tot.TRANSIENT.intr.extr.Rds")

all_fit_brms.tot.TRANSIENT.intr.extr.simp %>% rm()  # Simple - connectivity through both S and B
all_fit_brms.tot.TRANSIENT.intr.extr.simp <-brm(species_mod_inflow.intr.extr.simp + biom_mod_inflow.intr.extr.simp + set_rescor(FALSE), data=TRANSIENT.std,cores=4,chains = 4,
                                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.tot.TRANSIENT.intr.extr.simp,"all_fit_brms.tot.TRANSIENT.intr.extr.simp.Rds")


#TRANSIENT species mediated
all_fit_brms.nocon.S.TRANSIENT.intr %>% rm() # connectivity only through S + ENV
all_fit_brms.nocon.S.TRANSIENT.intr <-brm(species_mod_inflow.intr + biom_mod_nocon.S + set_rescor(FALSE), data=TRANSIENT.std,cores=4,chains = 4,
                                     iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                     prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.nocon.S.TRANSIENT.intr,"all_fit_brms.nocon.S.TRANSIENT.intr.Rds")

all_fit_brms.nocon.S.TRANSIENT.intr.extr %>% rm() # connectivity only through S + ENV
all_fit_brms.nocon.S.TRANSIENT.intr.extr <-brm(species_mod_inflow.intr.extr + biom_mod_nocon.S + set_rescor(FALSE), data=TRANSIENT.std,cores=4,chains = 4,
                                          iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                          prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.nocon.S.TRANSIENT.intr.extr,"all_fit_brms.nocon.S.TRANSIENT.intr.extr.Rds")

#
all_fit_brms.nocon.TRANSIENT %>% rm() # no connectivity for S and B, only environment
all_fit_brms.nocon.TRANSIENT <-brm(S_mod_nocon + biom_mod_nocon.S + set_rescor(FALSE), data=TRANSIENT.std,cores=4,chains = 4,
                                  iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                  prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.nocon.TRANSIENT,"all_fit_brms.nocon.TRANSIENT.Rds")


TRANSIENT.weight <- model_weights(all_fit_brms.tot.TRANSIENT.intr,
                                  all_fit_brms.tot.TRANSIENT.intr.extr,
                                  all_fit_brms.nocon.S.TRANSIENT.intr,
                                  all_fit_brms.nocon.S.TRANSIENT.intr.extr,
                                  all_fit_brms.tot.TRANSIENT.intr.extr.simp,
                                  all_fit_brms.nocon.TRANSIENT,weights = "loo")

TRANSIENT.LOO <- LOO(all_fit_brms.tot.TRANSIENT.intr,
                     all_fit_brms.tot.TRANSIENT.intr.extr,
                     all_fit_brms.nocon.S.TRANSIENT.intr,
                     all_fit_brms.nocon.S.TRANSIENT.intr.extr,
                     all_fit_brms.tot.TRANSIENT.intr.extr.simp,
                     all_fit_brms.nocon.TRANSIENT)
WAIC(all_fit_brms.tot.TRANSIENT.intr,
     all_fit_brms.tot.TRANSIENT.intr.extr,
     all_fit_brms.nocon.S.TRANSIENT.intr,
     all_fit_brms.nocon.S.TRANSIENT.intr.extr,
     all_fit_brms.tot.TRANSIENT.intr.extr.simp,
     all_fit_brms.nocon.TRANSIENT)

TRANSIENT.R2 <- rbind(bayes_R2(all_fit_brms.tot.TRANSIENT.intr),
bayes_R2(all_fit_brms.tot.TRANSIENT.intr.extr),
bayes_R2(all_fit_brms.nocon.S.TRANSIENT.intr),
bayes_R2(all_fit_brms.nocon.S.TRANSIENT.intr.extr),
bayes_R2(all_fit_brms.tot.TRANSIENT.intr.extr.simp),
bayes_R2(all_fit_brms.nocon.TRANSIENT))

mcmc_plot(all_fit_brms.tot.TRANSIENT.intr)
mcmc_plot(all_fit_brms.tot.TRANSIENT.intr.extr)
mcmc_plot(all_fit_brms.nocon.S.TRANSIENT.intr)
mcmc_plot(all_fit_brms.nocon.S.TRANSIENT.intr.extr)
mcmc_plot(all_fit_brms.nocon.TRANSIENT)

describe_posterior(all_fit_brms.nocon.TRANSIENT)



library("bayesplot")
library("ggplot2")
library("rstanarm")   
require(ggpubr)

posterior.TRANSIENT.tot <- as.array(all_fit_brms.tot.TRANSIENT.intr.extr)
posterior.TRANSIENT.simp <- as.array(all_fit_brms.simp.TRANSIENT)

# RESID
all_fit_brms.tot.RESID.intr %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.RESID.intr <-brm(species_mod_inflow.intr + biom_mod_inflow.intr + set_rescor(FALSE), data=RESID.std,cores=4,chains = 4,
                                      iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                      prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.tot.RESID.intr,"all_fit_brms.tot.RESID.intr.Rds")


all_fit_brms.tot.RESID.intr.extr %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.RESID.intr.extr <-brm(species_mod_inflow.intr.extr + biom_mod_inflow.intr.extr + set_rescor(FALSE), data=RESID.std,cores=4,chains = 4,
                                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.tot.RESID.intr.extr,"all_fit_brms.tot.RESID.intr.extr.Rds")

all_fit_brms.tot.RESID.intr.extr.simp %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.RESID.intr.extr.simp <-brm(species_mod_inflow.intr.extr.simp + biom_mod_inflow.intr.extr.simp + set_rescor(FALSE), data=RESID.std,cores=4,chains = 4,
                                                iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                                prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.tot.RESID.intr.extr.simp,"all_fit_brms.tot.RESID.intr.extr.simp.Rds")


#RESID species mediated
all_fit_brms.nocon.S.RESID.intr %>% rm() # connectivity only through S + ENV
all_fit_brms.nocon.S.RESID.intr <-brm(species_mod_inflow.intr + biom_mod_nocon.S + set_rescor(FALSE), data=RESID.std,cores=4,chains = 4,
                                          iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                          prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.nocon.S.RESID.intr,"all_fit_brms.nocon.S.RESID.intr.Rds")

all_fit_brms.nocon.S.RESID.intr.extr %>% rm() # connectivity only through S + ENV
all_fit_brms.nocon.S.RESID.intr.extr <-brm(species_mod_inflow.intr.extr + biom_mod_nocon.S + set_rescor(FALSE), data=RESID.std,cores=4,chains = 4,
                                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.nocon.S.RESID.intr.extr,"all_fit_brms.nocon.S.RESID.intr.extr.Rds")

#
all_fit_brms.nocon.RESID %>% rm() # no connectivity for S and B, only environment
all_fit_brms.nocon.RESID <-brm(S_mod_nocon + biom_mod_nocon.S + set_rescor(FALSE), data=RESID.std,cores=4,chains = 4,
                                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                   prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.nocon.RESID,"all_fit_brms.nocon.RESID.Rds")


RESID.weight <- model_weights(all_fit_brms.tot.RESID.intr,
                                  all_fit_brms.tot.RESID.intr.extr,
                                  all_fit_brms.nocon.S.RESID.intr,
                                  all_fit_brms.nocon.S.RESID.intr.extr,
                                  all_fit_brms.tot.RESID.intr.extr.simp,
                                  all_fit_brms.nocon.RESID,weights = "loo")

RESID.LOO <- LOO(all_fit_brms.tot.RESID.intr,
                     all_fit_brms.tot.RESID.intr.extr,
                     all_fit_brms.nocon.S.RESID.intr,
                     all_fit_brms.nocon.S.RESID.intr.extr,
                     all_fit_brms.tot.RESID.intr.extr.simp,
                     all_fit_brms.nocon.RESID)
WAIC(all_fit_brms.tot.RESID.intr,
     all_fit_brms.tot.RESID.intr.extr,
     all_fit_brms.nocon.S.RESID.intr,
     all_fit_brms.nocon.S.RESID.intr.extr,
     all_fit_brms.tot.RESID.intr.extr.simp,
     all_fit_brms.nocon.RESID)

RESID.R2 <- rbind(bayes_R2(all_fit_brms.tot.RESID.intr),
                      bayes_R2(all_fit_brms.tot.RESID.intr.extr),
                      bayes_R2(all_fit_brms.nocon.S.RESID.intr),
                      bayes_R2(all_fit_brms.nocon.S.RESID.intr.extr),
                      bayes_R2(all_fit_brms.tot.RESID.intr.extr.simp),
                      bayes_R2(all_fit_brms.nocon.RESID))

mcmc_plot(all_fit_brms.tot.RESID.intr)
mcmc_plot(all_fit_brms.tot.RESID.intr.extr)
mcmc_plot(all_fit_brms.nocon.S.RESID.intr)
mcmc_plot(all_fit_brms.nocon.S.RESID.intr.extr)
mcmc_plot(all_fit_brms.nocon.RESID)

# PARENTAL
all_fit_brms.tot.PARENTAL.intr %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.PARENTAL.intr <-brm(species_mod_inflow.intr + biom_mod_inflow.intr + set_rescor(FALSE), data=PARENTAL.std,cores=4,chains = 4,
                                      iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                      prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.tot.PARENTAL.intr,"all_fit_brms.tot.PARENTAL.intr.Rds")

all_fit_brms.tot.PARENTAL.intr.extr %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.PARENTAL.intr.extr <-brm(species_mod_inflow.intr.extr + biom_mod_inflow.intr.extr + set_rescor(FALSE), data=PARENTAL.std,cores=4,chains = 4,
                                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.tot.PARENTAL.intr.extr,"all_fit_brms.tot.PARENTAL.intr.extr.Rds")

all_fit_brms.tot.PARENTAL.intr.extr.simp %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.PARENTAL.intr.extr.simp <-brm(species_mod_inflow.intr.extr.simp + biom_mod_inflow.intr.extr.simp + set_rescor(FALSE), data=PARENTAL.std,cores=4,chains = 4,
                                                iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                                prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.tot.PARENTAL.intr.extr.simp,"all_fit_brms.tot.PARENTAL.intr.extr.simp.Rds")


#PARENTAL species mediated
all_fit_brms.nocon.S.PARENTAL.intr %>% rm() # connectivity only through S + ENV
all_fit_brms.nocon.S.PARENTAL.intr <-brm(species_mod_inflow.intr + biom_mod_nocon.S + set_rescor(FALSE), data=PARENTAL.std,cores=4,chains = 4,
                                          iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                          prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.nocon.S.PARENTAL.intr,"all_fit_brms.nocon.S.PARENTAL.intr.Rds")

all_fit_brms.nocon.S.PARENTAL.intr.extr %>% rm() # connectivity only through S + ENV
all_fit_brms.nocon.S.PARENTAL.intr.extr <-brm(species_mod_inflow.intr.extr + biom_mod_nocon.S + set_rescor(FALSE), data=PARENTAL.std,cores=4,chains = 4,
                                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.nocon.S.PARENTAL.intr.extr,"all_fit_brms.nocon.S.PARENTAL.intr.extr.Rds")

#
all_fit_brms.nocon.PARENTAL %>% rm() # no connectivity for S and B, only environment
all_fit_brms.nocon.PARENTAL <-brm(S_mod_nocon + biom_mod_nocon.S + set_rescor(FALSE), data=PARENTAL.std,cores=4,chains = 4,
                                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                   prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.nocon.PARENTAL,"all_fit_brms.nocon.PARENTAL.Rds")


PARENTAL.weight <- model_weights(all_fit_brms.tot.PARENTAL.intr,
                                  all_fit_brms.tot.PARENTAL.intr.extr,
                                  all_fit_brms.nocon.S.PARENTAL.intr,
                                  all_fit_brms.nocon.S.PARENTAL.intr.extr,
                                  all_fit_brms.tot.PARENTAL.intr.extr.simp,
                                  all_fit_brms.nocon.PARENTAL,weights = "loo")

PARENTAL.LOO <- LOO(all_fit_brms.tot.PARENTAL.intr,
                     all_fit_brms.tot.PARENTAL.intr.extr,
                     all_fit_brms.nocon.S.PARENTAL.intr,
                     all_fit_brms.nocon.S.PARENTAL.intr.extr,
                     all_fit_brms.tot.PARENTAL.intr.extr.simp,
                     all_fit_brms.nocon.PARENTAL)
WAIC(all_fit_brms.tot.PARENTAL.intr,
     all_fit_brms.tot.PARENTAL.intr.extr,
     all_fit_brms.nocon.S.PARENTAL.intr,
     all_fit_brms.nocon.S.PARENTAL.intr.extr,
     all_fit_brms.tot.PARENTAL.intr.extr.simp,
     all_fit_brms.nocon.PARENTAL)

PARENTAL.R2 <- rbind(bayes_R2(all_fit_brms.tot.PARENTAL.intr),
                      bayes_R2(all_fit_brms.tot.PARENTAL.intr.extr),
                      bayes_R2(all_fit_brms.nocon.S.PARENTAL.intr),
                      bayes_R2(all_fit_brms.nocon.S.PARENTAL.intr.extr),
                      bayes_R2(all_fit_brms.tot.PARENTAL.intr.extr.simp),
                      bayes_R2(all_fit_brms.nocon.PARENTAL))

mcmc_plot(all_fit_brms.tot.PARENTAL.intr)
mcmc_plot(all_fit_brms.tot.PARENTAL.intr.extr)
mcmc_plot(all_fit_brms.nocon.S.PARENTAL.intr)
mcmc_plot(all_fit_brms.nocon.S.PARENTAL.intr.extr)
mcmc_plot(all_fit_brms.nocon.PARENTAL)

# CRYPTIC
all_fit_brms.tot.CRYPTIC.intr %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.CRYPTIC.intr <-brm(species_mod_inflow.intr + biom_mod_inflow.intr + set_rescor(FALSE), data=CRYPTIC.std,cores=4,chains = 4,
                                      iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                      prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

all_fit_brms.tot.CRYPTIC.intr.extr %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.CRYPTIC.intr.extr <-brm(species_mod_inflow.intr.extr + biom_mod_inflow.intr.extr + set_rescor(FALSE), data=CRYPTIC.std,cores=4,chains = 4,
                                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(all_fit_brms.tot.CRYPTIC.intr.extr,"all_fit_brms.tot.CRYPTIC.intr.extr.Rds")

all_fit_brms.tot.CRYPTIC.intr.extr.simp %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.CRYPTIC.intr.extr.simp <-brm(species_mod_inflow.intr.extr.simp + biom_mod_inflow.intr.extr.simp + set_rescor(FALSE), data=CRYPTIC.std,cores=4,chains = 4,
                                                iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                                prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))


#CRYPTIC species mediated
all_fit_brms.nocon.S.CRYPTIC.intr %>% rm() # connectivity only through S + ENV
all_fit_brms.nocon.S.CRYPTIC.intr <-brm(species_mod_inflow.intr + biom_mod_nocon.S + set_rescor(FALSE), data=CRYPTIC.std,cores=4,chains = 4,
                                          iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                          prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
all_fit_brms.nocon.S.CRYPTIC.intr.extr %>% rm() # connectivity only through S + ENV
all_fit_brms.nocon.S.CRYPTIC.intr.extr <-brm(species_mod_inflow.intr.extr + biom_mod_nocon.S + set_rescor(FALSE), data=CRYPTIC.std,cores=4,chains = 4,
                                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
#
all_fit_brms.nocon.CRYPTIC %>% rm() # no connectivity for S and B, only environment
all_fit_brms.nocon.CRYPTIC <-brm(S_mod_nocon + biom_mod_nocon.S + set_rescor(FALSE), data=CRYPTIC.std,cores=4,chains = 4,
                                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                   prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))


CRYPTIC.weight <- model_weights(all_fit_brms.tot.CRYPTIC.intr,
                                  all_fit_brms.tot.CRYPTIC.intr.extr,
                                  all_fit_brms.nocon.S.CRYPTIC.intr,
                                  all_fit_brms.nocon.S.CRYPTIC.intr.extr,
                                  all_fit_brms.tot.CRYPTIC.intr.extr.simp,
                                  all_fit_brms.nocon.CRYPTIC,weights = "loo")

CRYPTIC.LOO <- LOO(all_fit_brms.tot.CRYPTIC.intr,
                     all_fit_brms.tot.CRYPTIC.intr.extr,
                     all_fit_brms.nocon.S.CRYPTIC.intr,
                     all_fit_brms.nocon.S.CRYPTIC.intr.extr,
                     all_fit_brms.tot.CRYPTIC.intr.extr.simp,
                     all_fit_brms.nocon.CRYPTIC)
WAIC(all_fit_brms.tot.CRYPTIC.intr,
     all_fit_brms.tot.CRYPTIC.intr.extr,
     all_fit_brms.nocon.S.CRYPTIC.intr,
     all_fit_brms.nocon.S.CRYPTIC.intr.extr,
     all_fit_brms.tot.CRYPTIC.intr.extr.simp,
     all_fit_brms.nocon.CRYPTIC)

CRYPTIC.R2 <- rbind(bayes_R2(all_fit_brms.tot.CRYPTIC.intr),
                      bayes_R2(all_fit_brms.tot.CRYPTIC.intr.extr),
                      bayes_R2(all_fit_brms.nocon.S.CRYPTIC.intr),
                      bayes_R2(all_fit_brms.nocon.S.CRYPTIC.intr.extr),
                      bayes_R2(all_fit_brms.tot.CRYPTIC.intr.extr.simp),
                      bayes_R2(all_fit_brms.nocon.CRYPTIC))

mcmc_plot(all_fit_brms.tot.CRYPTIC.intr)
mcmc_plot(all_fit_brms.tot.CRYPTIC.intr.extr)
mcmc_plot(all_fit_brms.nocon.S.CRYPTIC.intr)
mcmc_plot(all_fit_brms.nocon.S.CRYPTIC.intr.extr)
mcmc_plot(all_fit_brms.nocon.CRYPTIC)

#
library("bayesplot")
library("ggplot2")
library("rstanarm")   
require(ggpubr)

# classify explanatory variables based on categories


# TRANSIENT TOT
a <- mcmc_intervals(all_fit_brms.tot.TRANSIENT.intr.extr)

# no classification
rich.var <- as.data.frame(a$data)[grep("b_Richness",as.data.frame(a$data)[,"parameter"]),"parameter"]
biom.var <- as.data.frame(a$data)[grep("b_logbiomassarea",as.data.frame(a$data)[,"parameter"]),"parameter"]

Richness_inflow_TRANSIENT_tot <- mcmc_intervals(all_fit_brms.tot.TRANSIENT.intr.extr, pars = as.character(rich.var))
Biomass_inflow_TRANSIENT_tot <- mcmc_intervals(all_fit_brms.tot.TRANSIENT.intr.extr,pars = as.character(biom.var))

rm(a)

# Classification
INT <- "Intercept"
INTR <- c("Netflow","log_InflowBR","log_SelfR","log_btwdegree")
EXTR <- c("log_InflowMPABR","log_InflowNeiBR","log_CorridorIndegreeBR","log_grav_neiBR")
HUMAN <- c("_ClassClosed","_ClassRestricted","log_grav_total","log_grav_total:ClassClosed","log_grav_total:ClassRestricted")
ENV <- "temp"
HUMAN_INTR <- c("log_grav_total:Netflow","ClassClosed:Netflow","ClassRestricted:Netflow","log_grav_total:ClassClosed:Netflow","log_grav_total:ClassRestricted:Netflow")

Richness_inflow_TRANSIENT_tot_data <- mcmc_intervals_data(all_fit_brms.tot.TRANSIENT.intr.extr, pars = as.character(rich.var))
Biomass_inflow_TRANSIENT_tot_data <- mcmc_intervals_data(all_fit_brms.tot.TRANSIENT.intr.extr,pars = as.character(biom.var))


library(dotwhisker)
library(broom)
library(dplyr)
#m1_df <- tidy(b) %>% filter(term != "Intercept") %>% mutate(model = "Model 1")

PARAM <- c("Intercept","Richness","Temperature","Tot.Gravity","No-Take","Restricted gears",
"Netflow","Inflow","Self Recruit.","Betweeness Centr.",
"Gravity Neighb.","Inflow MPA","Inflow Neighb.","Corridor Indegree",                      
"Tot.Gravity x No-Take","Tot.Gravity x Restricted gears",
"Tot.Gravity x Netflow",
"No-Take x Netflow","Restricted gears x Netflow",
"Tot.Gravity x No-Take x Netflow","Tot.Gravity x Restricted gears x Netflow")

CAT <- c("Intercept","Richness","Human/Env.","Human/Env.","Human/Env.","Human/Env.",
           "Intrinsic Connect.","Intrinsic Connect.","Intrinsic Connect.","Intrinsic Connect.",
           "Extrinsic Connect.","Extrinsic Connect.","Extrinsic Connect.","Extrinsic Connect.",                      
           "Human/Env.","Human/Env.",
           "Intrinsic Connect.",
           "Intrinsic Connect.","Intrinsic Connect.",
           "Intrinsic Connect.","Intrinsic Connect.")


DesiredOrder <- c("Intercept",
                  "Richness",
                  "Temperature",
                  "Tot.Gravity",
                  "No-Take",
                  "Restricted gears",
                  "Tot.Gravity x No-Take","Tot.Gravity x Restricted gears",
                  "Netflow","Inflow","Self Recruit.","Betweeness Centr.",
                  "Gravity Neighb.","Inflow MPA","Inflow Neighb.","Corridor Indegree",                      
                  "Tot.Gravity x Netflow",
                  "No-Take x Netflow","Restricted gears x Netflow",
                  "Tot.Gravity x No-Take x Netflow","Tot.Gravity x Restricted gears x Netflow")
# Transient
rm(m1_df_transient)
m1_df_transient <- Richness_inflow_TRANSIENT_tot_data %>%  mutate(FE = "TRANSIENT") %>%  
  mutate(model = "Richness") %>% 
  add_row(point_est="median",ll = 0, l = 0,m=0,h=0,hh=0, FE="TRANSIENT",model="Richness",parameter="b_Richness_Richness",.before = 2)
m1_df_transient <- cbind(PARAM,m1_df_transient,CAT)
m1_df_transient$PARAM <- factor(m1_df_transient$PARAM, levels = rev(DesiredOrder))
m1_df_transient$CAT <- factor(m1_df_transient$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)

rm(m2_df_transient)
m2_df_transient <- Biomass_inflow_TRANSIENT_tot_data %>%  mutate(FE = "TRANSIENT") %>%  mutate(model = "Biomass")
m2_df_transient <- cbind(PARAM,m2_df_transient,CAT)
m2_df_transient$PARAM <- factor(m2_df_transient$PARAM, levels = rev(DesiredOrder))
m2_df_transient$CAT <- factor(m2_df_transient$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)

two_models <- rbind(m1_df_transient, m2_df_transient)
two_models$FE <- as.factor(two_models$FE)
two_models$model <- as.factor(two_models$model)
two_models$CAT <- as.factor(two_models$CAT)

rm(Transient.SEM)
Transient.SEM <- ggplot(two_models %>% filter(PARAM != "Intercept"), aes(colour = model)) +
geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
geom_linerange(aes(x = PARAM, ymin = ll,
                                ymax = hh,colour = model),lwd = 1/2, position = position_dodge(width = 1/2)) +
geom_point(aes(x = PARAM, y = m,colour = model), lwd = 1.5, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE") +
geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = model),lwd = 1, position = position_dodge(width = 1/2)) +
facet_grid(CAT ~ ., scales = "free", space = "free") +
#scale_x_discrete(labels= PARAM)  +
coord_flip() + theme_bw() +
labs(x="",y="Standardized coefficients") +
ggtitle("Transient") 
Transient.SEM  # The trick to these is position_dodge()
rm(two_models)
# RESIDENT
a <- mcmc_intervals(all_fit_brms.tot.RESID.intr.extr)
rich.var <- as.data.frame(a$data)[grep("b_Richness",as.data.frame(a$data)[,"parameter"]),"parameter"]
biom.var <- as.data.frame(a$data)[grep("b_logbiomassarea",as.data.frame(a$data)[,"parameter"]),"parameter"]

Richness_inflow_RESIDENT_tot <- mcmc_intervals(all_fit_brms.tot.RESID.intr.extr, pars = as.character(rich.var))
Biomass_inflow_RESIDENT_tot <- mcmc_intervals(all_fit_brms.tot.RESID.intr.extr,pars = as.character(biom.var))

Richness_inflow_RESIDENT_tot_data <- mcmc_intervals_data(all_fit_brms.tot.RESID.intr.extr, pars = as.character(rich.var))
Biomass_inflow_RESIDENT_tot_data <- mcmc_intervals_data(all_fit_brms.tot.RESID.intr.extr,pars = as.character(biom.var))

rm(m1_df_resident)
m1_df_resident <- Richness_inflow_RESIDENT_tot_data %>%  mutate(FE = "RESIDENT") %>%  
  mutate(model = "Richness") %>% 
  add_row(point_est="median",ll = 0, l = 0,m=0,h=0,hh=0, FE="RESIDENT",model="Richness",parameter="b_Richness_Richness",.before = 2)
m1_df_resident <- cbind(PARAM,m1_df_resident,CAT)
m1_df_resident$PARAM <- factor(m1_df_resident$PARAM, levels = rev(DesiredOrder))
m1_df_resident$CAT <- factor(m1_df_resident$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)

rm(m2_df_resident)
m2_df_resident <- Biomass_inflow_RESIDENT_tot_data %>%  mutate(FE = "RESIDENT") %>%  mutate(model = "Biomass")
m2_df_resident <- cbind(PARAM,m2_df_resident,CAT)
m2_df_resident$PARAM <- factor(m2_df_resident$PARAM, levels = rev(DesiredOrder))
m2_df_resident$CAT <- factor(m2_df_resident$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)


two_models <- rbind(m1_df_resident, m2_df_resident)
two_models$FE <- as.factor(two_models$FE)
two_models$model <- as.factor(two_models$model)
two_models$CAT <- as.factor(two_models$CAT)

rm(Resident.SEM)
Resident.SEM <- ggplot(two_models %>% filter(PARAM != "Intercept"), aes(colour = model)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = model),lwd = 1/2, position = position_dodge(width = 1/2)) +
  geom_point(aes(x = PARAM, y = m,colour = model), lwd = 1.5, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = model),lwd = 1, position = position_dodge(width = 1/2)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") +
  ggtitle("Resident") 
Resident.SEM  # The trick to these is position_dodge()
rm(two_models)
# PARENTAL
a <- mcmc_intervals(all_fit_brms.tot.PARENTAL.intr.extr)
rich.var <- as.data.frame(a$data)[grep("b_Richness",as.data.frame(a$data)[,"parameter"]),"parameter"]
biom.var <- as.data.frame(a$data)[grep("b_logbiomassarea",as.data.frame(a$data)[,"parameter"]),"parameter"]

Richness_inflow_PARENTAL_tot <- mcmc_intervals(all_fit_brms.tot.PARENTAL.intr.extr, pars = as.character(rich.var))
Biomass_inflow_PARENTAL_tot <- mcmc_intervals(all_fit_brms.tot.PARENTAL.intr.extr,pars = as.character(biom.var))

Richness_inflow_PARENTAL_tot_data <- mcmc_intervals_data(all_fit_brms.tot.PARENTAL.intr.extr, pars = as.character(rich.var))
Biomass_inflow_PARENTAL_tot_data <- mcmc_intervals_data(all_fit_brms.tot.PARENTAL.intr.extr,pars = as.character(biom.var))

rm(m1_df_parental)
m1_df_parental <- Richness_inflow_PARENTAL_tot_data %>%  mutate(FE = "PARENTAL") %>%  
  mutate(model = "Richness") %>% 
  add_row(point_est="median",ll = 0, l = 0,m=0,h=0,hh=0, FE="PARENTAL",model="Richness",parameter="b_Richness_Richness",.before = 2)
m1_df_parental <- cbind(PARAM,m1_df_parental,CAT)
m1_df_parental$PARAM <- factor(m1_df_parental$PARAM, levels = rev(DesiredOrder))
m1_df_parental$CAT <- factor(m1_df_parental$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)

rm(m2_df_parental)
m2_df_parental <- Biomass_inflow_PARENTAL_tot_data %>%  mutate(FE = "PARENTAL") %>%  mutate(model = "Biomass")
m2_df_parental <- cbind(PARAM,m2_df_parental,CAT)
m2_df_parental$PARAM <- factor(m2_df_parental$PARAM, levels = rev(DesiredOrder))
m2_df_parental$CAT <- factor(m2_df_parental$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)


two_models <- rbind(m1_df_parental, m2_df_parental)
two_models$FE <- as.factor(two_models$FE)
two_models$model <- as.factor(two_models$model)
two_models$CAT <- as.factor(two_models$CAT)

rm(Parental.SEM)
Parental.SEM <- ggplot(two_models %>% filter(PARAM != "Intercept"), aes(colour = model)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = model),lwd = 1/2, position = position_dodge(width = 1/2)) +
  geom_point(aes(x = PARAM, y = m,colour = model), lwd = 1.5, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = model),lwd = 1, position = position_dodge(width = 1/2)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") +
  ggtitle("Parental") 
Parental.SEM  # The trick to these is position_dodge()

rm(two_models)
# CRYPTIC
a <- mcmc_intervals(all_fit_brms.tot.CRYPTIC.intr.extr)
rich.var <- as.data.frame(a$data)[grep("b_Richness",as.data.frame(a$data)[,"parameter"]),"parameter"]
biom.var <- as.data.frame(a$data)[grep("b_logbiomassarea",as.data.frame(a$data)[,"parameter"]),"parameter"]

Richness_inflow_CRYPTIC_tot <- mcmc_intervals(all_fit_brms.tot.CRYPTIC.intr.extr, pars = as.character(rich.var))
Biomass_inflow_CRYPTIC_tot <- mcmc_intervals(all_fit_brms.tot.CRYPTIC.intr.extr,pars = as.character(biom.var))

Richness_inflow_CRYPTIC_tot_data <- mcmc_intervals_data(all_fit_brms.tot.CRYPTIC.intr.extr, pars = as.character(rich.var))
Biomass_inflow_CRYPTIC_tot_data <- mcmc_intervals_data(all_fit_brms.tot.CRYPTIC.intr.extr,pars = as.character(biom.var))

rm(m1_df_cryptic)
m1_df_cryptic <- Richness_inflow_CRYPTIC_tot_data %>%  mutate(FE = "CRYPTIC") %>%  
  mutate(model = "Richness") %>% 
  add_row(point_est="median",ll = 0, l = 0,m=0,h=0,hh=0, FE="CRYPTIC",model="Richness",parameter="b_Richness_Richness",.before = 2)
m1_df_cryptic <- cbind(PARAM,m1_df_cryptic,CAT)
m1_df_cryptic$PARAM <- factor(m1_df_cryptic$PARAM, levels = rev(DesiredOrder))
m1_df_cryptic$CAT <- factor(m1_df_cryptic$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)

rm(m2_df_cryptic)
m2_df_cryptic <- Biomass_inflow_CRYPTIC_tot_data %>%  mutate(FE = "CRYPTIC") %>%  mutate(model = "Biomass")
m2_df_cryptic <- cbind(PARAM,m2_df_cryptic,CAT)
m2_df_cryptic$PARAM <- factor(m2_df_cryptic$PARAM, levels = rev(DesiredOrder))
m2_df_cryptic$CAT <- factor(m2_df_cryptic$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)


two_models <- rbind(m1_df_cryptic, m2_df_cryptic)
two_models$FE <- as.factor(two_models$FE)
two_models$model <- as.factor(two_models$model)
two_models$CAT <- as.factor(two_models$CAT)

rm(Cryptic.SEM)
Cryptic.SEM <- ggplot(two_models %>% filter(PARAM != "Intercept"), aes(colour = model)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = model),lwd = 1/2, position = position_dodge(width = 1/2)) +
  geom_point(aes(x = PARAM, y = m,colour = model), lwd = 1.5, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = model),lwd = 1, position = position_dodge(width = 1/2)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") +
  ggtitle("Cryptic") 
Cryptic.SEM  # The trick to these is position_dodge()

# total

SEM_inflow_tot <- ggarrange(Richness_inflow_TRANSIENT_tot,Richness_inflow_RESIDENT_tot,
                            Richness_inflow_PARENTAL_tot,Richness_inflow_CRYPTIC_tot, 
                            Biomass_inflow_TRANSIENT_tot,
                            Biomass_inflow_PARENTAL_tot,
                            Biomass_inflow_CRYPTIC_tot,
                            Biomass_inflow_RESIDENT_tot,
                            ncol=4,nrow=2,labels = c("A_Richness_Transient","B_Richness_Resident","C_Richness_Parental","D_Richness_Cryptic",
                                                     "E_Biomass_Transient","F_Biomass_Resident","G_Biomass_Parental","H_Biomass_Cryptic"),align="hv")
ggsave(here("_prelim.figures","SEM_bayes_inflow_coef_total_FE.pdf"),SEM_inflow_tot,width=30,height=15)

SEM_inflow_tot_weight <- round(rbind(TRANSIENT.weight[c(1,2,5,3,4,6)],RESID.weight[c(1,2,5,3,4,6)],
                                     PARENTAL.weight[c(1,2,5,3,4,6)],CRYPTIC.weight[c(1,2,5,3,4,6)]),2)
colnames(SEM_inflow_tot_weight) <- c("Model1_Full_intr","Model1_Full_intr_extr","Model1_Full_intr_extr_simp","Model2_Smed_intr","Model2_Smed_intr_extr","Model3_noCon")
rownames(SEM_inflow_tot_weight) <- c("TRANSIENT","RESIDENT","PARENTAL","CRYPTIC")

require(kableExtra)
SEM_inflow_tot_weight.table <- SEM_inflow_tot_weight  %>%
kable(align="c") %>%
  kable_styling() %>%
  save_kable("_prelim.figures/SEM_bayes_inflow_weight_FE.png")

SEM_inflow_tot_R2 <- round(rbind(TRANSIENT.R2[c(1,2,3,4,9,10,5,6,7,8,11,12),1],
                                 RESID.R2[c(1,2,3,4,9,10,5,6,7,8,11,12),1],
                                 PARENTAL.R2[c(1,2,3,4,9,10,5,6,7,8,11,12),1],
                                 CRYPTIC.R2[c(1,2,3,4,9,10,5,6,7,8,11,12),1]),2)
colnames(SEM_inflow_tot_R2) <- c("Model1_Full_S_intr","Model1_Full_B_intr",
                                 "Model1_Full_S_intr_extr","Model1_Full_B_intr_extr",
                                 "Model1_Full_S_intr_extr_simp","Model1_Full_B_intr_extr_simp",
                                 "Model2_Smed_S_intr","Model2_Smed_B_intr",
                                 "Model2_Smed_S_intr_extr","Model2_Smed_B_intr_extr",
                                 "Model3_noCon_S","Model3_noCon_B")
rownames(SEM_inflow_tot_R2) <- c("TRANSIENT",
                                     "RESIDENT",
                                     "PARENTAL",
                                     "CRYPTIC")
SEM_inflow_tot_R2.table <- SEM_inflow_tot_R2  %>%
  kable(align="c") %>%
  kable_styling() %>%
  save_kable("_prelim.figures/SEM_bayes_inflow_R2_FE.png")

#####
SEM_coef_tot <- ggarrange(Transient.SEM,Resident.SEM,Parental.SEM,Cryptic.SEM,
                            ncol=2,nrow=2,labels = c("A","B","C","D"),align="hv",common.legend = T,legend="bottom")
ggsave(here("_prelim.figures","SEM_bayes_coef_total_FE.pdf"),SEM_coef_tot,width=15,height=15)

#### two plots 1) Richness 2) Biomass

  # richness data.frame
four_FE_Richness <- rbind(m1_df_transient,m1_df_resident,m1_df_parental,m1_df_cryptic)
four_FE_Richness$FE <- as.factor(four_FE_Richness$FE)
four_FE_Richness$CAT <- as.factor(four_FE_Richness$CAT)
four_FE_Richness$model <- as.factor(four_FE_Richness$model)
summary(four_FE_Richness)
  # biomass data.frame
four_FE_Biomass <- rbind(m2_df_transient,m2_df_resident,m2_df_parental,m2_df_cryptic)
four_FE_Biomass$FE <- as.factor(four_FE_Biomass$FE)
four_FE_Biomass$CAT <- as.factor(four_FE_Biomass$CAT)
four_FE_Biomass$model <- as.factor(four_FE_Biomass$model)
summary(four_FE_Biomass)

# plot richness
rm(Richness.FE)
Richness.FE <- ggplot(four_FE_Richness %>% filter(PARAM != "Intercept" & PARAM != "Richness"), aes(colour = FE)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = FE),lwd = 1/2, position = position_dodge(width = 0.8)) +
  geom_point(aes(x = PARAM, y = m,colour = FE), lwd = 1.5, position = position_dodge(width = 0.8),shape = 21, fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = FE),lwd = 1, position = position_dodge(width = 0.8)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") +
  ggtitle("Richness") 
Richness.FE  # The trick to these is position_dodge()

# plot biomass
rm(Biomass.FE)
Biomass.FE <- ggplot(four_FE_Biomass %>% filter(PARAM != "Intercept"), aes(colour = FE)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = FE),lwd = 1/2, position = position_dodge(width = 0.8)) +
  geom_point(aes(x = PARAM, y = m,colour = FE), lwd = 1.5, position = position_dodge(width = 0.8),shape = 21, fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = FE),lwd = 1, position = position_dodge(width = 0.8)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") +
  ggtitle("Biomass") 
Biomass.FE  # The trick to these is position_dodge()

SEM_coef_tot_rich_biom <- ggarrange(Richness.FE,Biomass.FE,
                          ncol=2,nrow=1,labels = c("A","B"),align="hv",common.legend = T,legend="bottom")
ggsave(here("_prelim.figures","SEM_bayes_coef_total_FE_rich_biom.pdf"),SEM_coef_tot_rich_biom,width=15,height=15)

# save parameters for each model
write.csv(four_FE_Biomass,here("_prelim.figures","SEM_bayes_coef_total_FE_Biomass.csv"))
write.csv(four_FE_Richness,here("_prelim.figures","SEM_bayes_coef_total_FE_Richness.csv"))


###### SIMPLE MODELS ######

# TRANSIENT TOT
a <- mcmc_intervals(all_fit_brms.tot.TRANSIENT.intr.extr.simp)

# no classification
rich.var <- as.data.frame(a$data)[grep("b_Richness",as.data.frame(a$data)[,"parameter"]),"parameter"]
biom.var <- as.data.frame(a$data)[grep("b_logbiomassarea",as.data.frame(a$data)[,"parameter"]),"parameter"]

Richness_inflow_TRANSIENT_simp <- mcmc_intervals(all_fit_brms.tot.TRANSIENT.intr.extr.simp, pars = as.character(rich.var))
Biomass_inflow_TRANSIENT_simp <- mcmc_intervals(all_fit_brms.tot.TRANSIENT.intr.extr.simp,pars = as.character(biom.var))

rm(a)

# Classification
INT <- "Intercept"
INTR <- c("Netflow","log_InflowBR","log_SelfR")
EXTR <- c("log_InflowNeiBR","log_grav_neiBR")
HUMAN <- c("_ClassClosed","_ClassRestricted","log_grav_total","log_grav_total:ClassClosed","log_grav_total:ClassRestricted")
ENV <- "temp"
HUMAN_INTR <- c("log_grav_total:Netflow","ClassClosed:Netflow","ClassRestricted:Netflow","log_grav_total:ClassClosed:Netflow","log_grav_total:ClassRestricted:Netflow")

Richness_inflow_TRANSIENT_simp_data <- mcmc_intervals_data(all_fit_brms.tot.TRANSIENT.intr.extr.simp, pars = as.character(rich.var))
Biomass_inflow_TRANSIENT_simp_data <- mcmc_intervals_data(all_fit_brms.tot.TRANSIENT.intr.extr.simp,pars = as.character(biom.var))


library(dotwhisker)
library(broom)
library(dplyr)
#m1_df <- tidy(b) %>% filter(term != "Intercept") %>% mutate(model = "Model 1")

PARAM <- c("Intercept","Richness","Temperature","Tot.Gravity","No-Take","Restricted gears",
           "Netflow","Inflow","Self Recruit.",
           "Gravity Neighb.","Inflow Neighb.",                      
           "Tot.Gravity x No-Take","Tot.Gravity x Restricted gears",
           "Tot.Gravity x Netflow",
           "No-Take x Netflow","Restricted gears x Netflow",
           "Tot.Gravity x No-Take x Netflow","Tot.Gravity x Restricted gears x Netflow")

CAT <- c("Intercept","Richness","Human/Env.","Human/Env.","Human/Env.","Human/Env.",
         "Intrinsic Connect.","Intrinsic Connect.","Intrinsic Connect.",
         "Extrinsic Connect.","Extrinsic Connect.",                      
         "Human/Env.","Human/Env.",
         "Intrinsic Connect.",
         "Intrinsic Connect.","Intrinsic Connect.",
         "Intrinsic Connect.","Intrinsic Connect.")


DesiredOrder <- c("Intercept",
                  "Richness",
                  "Temperature",
                  "Tot.Gravity",
                  "No-Take",
                  "Restricted gears",
                  "Tot.Gravity x No-Take","Tot.Gravity x Restricted gears",
                  "Netflow","Inflow","Self Recruit.",
                  "Gravity Neighb.","Inflow Neighb.",                      
                  "Tot.Gravity x Netflow",
                  "No-Take x Netflow","Restricted gears x Netflow",
                  "Tot.Gravity x No-Take x Netflow","Tot.Gravity x Restricted gears x Netflow")
# Transient
rm(m1_df_transient)
m1_df_transient <- Richness_inflow_TRANSIENT_simp_data %>%  mutate(FE = "TRANSIENT") %>%  
  mutate(model = "Richness") %>% 
  add_row(point_est="median",ll = 0, l = 0,m=0,h=0,hh=0, FE="TRANSIENT",model="Richness",parameter="b_Richness_Richness",.before = 2)
m1_df_transient <- cbind(PARAM,m1_df_transient,CAT)
m1_df_transient$PARAM <- factor(m1_df_transient$PARAM, levels = rev(DesiredOrder))
m1_df_transient$CAT <- factor(m1_df_transient$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)

rm(m2_df_transient)
m2_df_transient <- Biomass_inflow_TRANSIENT_simp_data %>%  mutate(FE = "TRANSIENT") %>%  mutate(model = "Biomass")
m2_df_transient <- cbind(PARAM,m2_df_transient,CAT)
m2_df_transient$PARAM <- factor(m2_df_transient$PARAM, levels = rev(DesiredOrder))
m2_df_transient$CAT <- factor(m2_df_transient$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)

two_models <- rbind(m1_df_transient, m2_df_transient)
two_models$FE <- as.factor(two_models$FE)
two_models$model <- as.factor(two_models$model)
two_models$CAT <- as.factor(two_models$CAT)

rm(Transient.SEM)
Transient.SEM <- ggplot(two_models %>% filter(PARAM != "Intercept"), aes(colour = model)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = model),lwd = 1/2, position = position_dodge(width = 1/2)) +
  geom_point(aes(x = PARAM, y = m,colour = model), lwd = 1.5, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = model),lwd = 1, position = position_dodge(width = 1/2)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") +
  ggtitle("Transient") 
Transient.SEM  # The trick to these is position_dodge()
rm(two_models)
# RESIDENT
rm(a)
a <- mcmc_intervals(all_fit_brms.tot.RESID.intr.extr.simp)
rich.var <- as.data.frame(a$data)[grep("b_Richness",as.data.frame(a$data)[,"parameter"]),"parameter"]
biom.var <- as.data.frame(a$data)[grep("b_logbiomassarea",as.data.frame(a$data)[,"parameter"]),"parameter"]

Richness_inflow_RESIDENT_simp <- mcmc_intervals(all_fit_brms.tot.RESID.intr.extr.simp, pars = as.character(rich.var))
Biomass_inflow_RESIDENT_simp <- mcmc_intervals(all_fit_brms.tot.RESID.intr.extr.simp,pars = as.character(biom.var))

Richness_inflow_RESIDENT_simp_data <- mcmc_intervals_data(all_fit_brms.tot.RESID.intr.extr.simp, pars = as.character(rich.var))
Biomass_inflow_RESIDENT_simp_data <- mcmc_intervals_data(all_fit_brms.tot.RESID.intr.extr.simp,pars = as.character(biom.var))

rm(m1_df_resident)
m1_df_resident <- Richness_inflow_RESIDENT_simp_data %>%  mutate(FE = "RESIDENT") %>%  
  mutate(model = "Richness") %>% 
  add_row(point_est="median",ll = 0, l = 0,m=0,h=0,hh=0, FE="RESIDENT",model="Richness",parameter="b_Richness_Richness",.before = 2)
m1_df_resident <- cbind(PARAM,m1_df_resident,CAT)
m1_df_resident$PARAM <- factor(m1_df_resident$PARAM, levels = rev(DesiredOrder))
m1_df_resident$CAT <- factor(m1_df_resident$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)

rm(m2_df_resident)
m2_df_resident <- Biomass_inflow_RESIDENT_simp_data %>%  mutate(FE = "RESIDENT") %>%  mutate(model = "Biomass")
m2_df_resident <- cbind(PARAM,m2_df_resident,CAT)
m2_df_resident$PARAM <- factor(m2_df_resident$PARAM, levels = rev(DesiredOrder))
m2_df_resident$CAT <- factor(m2_df_resident$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)


two_models <- rbind(m1_df_resident, m2_df_resident)
two_models$FE <- as.factor(two_models$FE)
two_models$model <- as.factor(two_models$model)
two_models$CAT <- as.factor(two_models$CAT)

rm(Resident.SEM)
Resident.SEM <- ggplot(two_models %>% filter(PARAM != "Intercept"), aes(colour = model)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = model),lwd = 1/2, position = position_dodge(width = 1/2)) +
  geom_point(aes(x = PARAM, y = m,colour = model), lwd = 1.5, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = model),lwd = 1, position = position_dodge(width = 1/2)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") +
  ggtitle("Resident") 
Resident.SEM  # The trick to these is position_dodge()
rm(two_models)
# PARENTAL
rm(a)
a <- mcmc_intervals(all_fit_brms.tot.PARENTAL.intr.extr.simp)
rich.var <- as.data.frame(a$data)[grep("b_Richness",as.data.frame(a$data)[,"parameter"]),"parameter"]
biom.var <- as.data.frame(a$data)[grep("b_logbiomassarea",as.data.frame(a$data)[,"parameter"]),"parameter"]

Richness_inflow_PARENTAL_simp <- mcmc_intervals(all_fit_brms.tot.PARENTAL.intr.extr.simp, pars = as.character(rich.var))
Biomass_inflow_PARENTAL_simp <- mcmc_intervals(all_fit_brms.tot.PARENTAL.intr.extr.simp,pars = as.character(biom.var))

Richness_inflow_PARENTAL_simp_data <- mcmc_intervals_data(all_fit_brms.tot.PARENTAL.intr.extr.simp, pars = as.character(rich.var))
Biomass_inflow_PARENTAL_simp_data <- mcmc_intervals_data(all_fit_brms.tot.PARENTAL.intr.extr.simp,pars = as.character(biom.var))

rm(m1_df_parental)
m1_df_parental <- Richness_inflow_PARENTAL_simp_data %>%  mutate(FE = "PARENTAL") %>%  
  mutate(model = "Richness") %>% 
  add_row(point_est="median",ll = 0, l = 0,m=0,h=0,hh=0, FE="PARENTAL",model="Richness",parameter="b_Richness_Richness",.before = 2)
m1_df_parental <- cbind(PARAM,m1_df_parental,CAT)
m1_df_parental$PARAM <- factor(m1_df_parental$PARAM, levels = rev(DesiredOrder))
m1_df_parental$CAT <- factor(m1_df_parental$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)

rm(m2_df_parental)
m2_df_parental <- Biomass_inflow_PARENTAL_simp_data %>%  mutate(FE = "PARENTAL") %>%  mutate(model = "Biomass")
m2_df_parental <- cbind(PARAM,m2_df_parental,CAT)
m2_df_parental$PARAM <- factor(m2_df_parental$PARAM, levels = rev(DesiredOrder))
m2_df_parental$CAT <- factor(m2_df_parental$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)

two_models <- rbind(m1_df_parental, m2_df_parental)
two_models$FE <- as.factor(two_models$FE)
two_models$model <- as.factor(two_models$model)
two_models$CAT <- as.factor(two_models$CAT)

rm(Parental.SEM)
Parental.SEM <- ggplot(two_models %>% filter(PARAM != "Intercept"), aes(colour = model)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = model),lwd = 1/2, position = position_dodge(width = 1/2)) +
  geom_point(aes(x = PARAM, y = m,colour = model), lwd = 1.5, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = model),lwd = 1, position = position_dodge(width = 1/2)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") +
  ggtitle("Parental") 
Parental.SEM  # The trick to these is position_dodge()

rm(two_models)
# CRYPTIC
rm(a)
a <- mcmc_intervals(all_fit_brms.tot.CRYPTIC.intr.extr.simp)
rich.var <- as.data.frame(a$data)[grep("b_Richness",as.data.frame(a$data)[,"parameter"]),"parameter"]
biom.var <- as.data.frame(a$data)[grep("b_logbiomassarea",as.data.frame(a$data)[,"parameter"]),"parameter"]

Richness_inflow_CRYPTIC_simp <- mcmc_intervals(all_fit_brms.tot.CRYPTIC.intr.extr.simp, pars = as.character(rich.var))
Biomass_inflow_CRYPTIC_simp <- mcmc_intervals(all_fit_brms.tot.CRYPTIC.intr.extr.simp,pars = as.character(biom.var))

Richness_inflow_CRYPTIC_simp_data <- mcmc_intervals_data(all_fit_brms.tot.CRYPTIC.intr.extr.simp, pars = as.character(rich.var))
Biomass_inflow_CRYPTIC_simp_data <- mcmc_intervals_data(all_fit_brms.tot.CRYPTIC.intr.extr.simp,pars = as.character(biom.var))

rm(m1_df_cryptic)
m1_df_cryptic <- Richness_inflow_CRYPTIC_simp_data %>%  mutate(FE = "CRYPTIC") %>%  
  mutate(model = "Richness") %>% 
  add_row(point_est="median",ll = 0, l = 0,m=0,h=0,hh=0, FE="CRYPTIC",model="Richness",parameter="b_Richness_Richness",.before = 2)
m1_df_cryptic <- cbind(PARAM,m1_df_cryptic,CAT)
m1_df_cryptic$PARAM <- factor(m1_df_cryptic$PARAM, levels = rev(DesiredOrder))
m1_df_cryptic$CAT <- factor(m1_df_cryptic$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)

rm(m2_df_cryptic)
m2_df_cryptic <- Biomass_inflow_CRYPTIC_simp_data %>%  mutate(FE = "CRYPTIC") %>%  mutate(model = "Biomass")
m2_df_cryptic <- cbind(PARAM,m2_df_cryptic,CAT)
m2_df_cryptic$PARAM <- factor(m2_df_cryptic$PARAM, levels = rev(DesiredOrder))
m2_df_cryptic$CAT <- factor(m2_df_cryptic$CAT , levels = c("Intercept","Richness","Human/Env.","Intrinsic Connect.","Extrinsic Connect."),ordered = TRUE)


two_models <- rbind(m1_df_cryptic, m2_df_cryptic)
two_models$FE <- as.factor(two_models$FE)
two_models$model <- as.factor(two_models$model)
two_models$CAT <- as.factor(two_models$CAT)

rm(Cryptic.SEM)
Cryptic.SEM <- ggplot(two_models %>% filter(PARAM != "Intercept"), aes(colour = model)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = model),lwd = 1/2, position = position_dodge(width = 1/2)) +
  geom_point(aes(x = PARAM, y = m,colour = model), lwd = 1.5, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = model),lwd = 1, position = position_dodge(width = 1/2)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") +
  ggtitle("Cryptic") 
Cryptic.SEM  # The trick to these is position_dodge()

# total

SEM_inflow_simp <- ggarrange(Richness_inflow_TRANSIENT_simp,Richness_inflow_RESIDENT_simp,
                            Richness_inflow_PARENTAL_simp,Richness_inflow_CRYPTIC_simp, 
                            Biomass_inflow_TRANSIENT_simp,
                            Biomass_inflow_PARENTAL_simp,
                            Biomass_inflow_CRYPTIC_simp,
                            Biomass_inflow_RESIDENT_simp,
                            ncol=4,nrow=2,labels = c("A_Richness_Transient","B_Richness_Resident","C_Richness_Parental","D_Richness_Cryptic",
                                                     "E_Biomass_Transient","F_Biomass_Resident","G_Biomass_Parental","H_Biomass_Cryptic"),align="hv")
ggsave(here("_prelim.figures","SEM_bayes_inflow_coef_simp_FE.pdf"),SEM_inflow_simp,width=30,height=15)

#### two plots 1) Richness 2) Biomass

# richness data.frame
four_FE_Richness <- rbind(m1_df_transient,m1_df_resident,m1_df_parental,m1_df_cryptic)
four_FE_Richness$FE <- as.factor(four_FE_Richness$FE)
four_FE_Richness$CAT <- as.factor(four_FE_Richness$CAT)
four_FE_Richness$model <- as.factor(four_FE_Richness$model)
summary(four_FE_Richness)
# biomass data.frame
four_FE_Biomass <- rbind(m2_df_transient,m2_df_resident,m2_df_parental,m2_df_cryptic)
four_FE_Biomass$FE <- as.factor(four_FE_Biomass$FE)
four_FE_Biomass$CAT <- as.factor(four_FE_Biomass$CAT)
four_FE_Biomass$model <- as.factor(four_FE_Biomass$model)
summary(four_FE_Biomass)

# plot richness
rm(Richness.FE)
Richness.FE <- ggplot(four_FE_Richness %>% filter(PARAM != "Intercept" & PARAM != "Richness"), aes(colour = FE)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = FE),lwd = 1/2, position = position_dodge(width = 0.8)) +
  geom_point(aes(x = PARAM, y = m,colour = FE), lwd = 1.5, position = position_dodge(width = 0.8),shape = 21, fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = FE),lwd = 1, position = position_dodge(width = 0.8)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") +
  ggtitle("Richness") 
Richness.FE  # The trick to these is position_dodge()

# plot biomass
rm(Biomass.FE)
Biomass.FE <- ggplot(four_FE_Biomass %>% filter(PARAM != "Intercept"), aes(colour = FE)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = FE),lwd = 1/2, position = position_dodge(width = 0.8)) +
  geom_point(aes(x = PARAM, y = m,colour = FE), lwd = 1.5, position = position_dodge(width = 0.8),shape = 21, fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = FE),lwd = 1, position = position_dodge(width = 0.8)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") +
  ggtitle("Biomass") 
Biomass.FE  # The trick to these is position_dodge()

SEM_coef_simp_rich_biom <- ggarrange(Richness.FE,Biomass.FE,
                                    ncol=2,nrow=1,labels = c("A","B"),align="hv",common.legend = T,legend="bottom")
ggsave(here("_prelim.figures","SEM_bayes_coef_simp_FE_rich_biom.pdf"),SEM_coef_simp_rich_biom,width=15,height=15)

# save parameters for each model
write.csv(four_FE_Biomass,here("_prelim.figures","SEM_bayes_coef_simp_FE_Biomass.csv"))
write.csv(four_FE_Richness,here("_prelim.figures","SEM_bayes_coef_simp_FE_Richness.csv"))


# END





# with latitude
posterior_lat <- as.array(all_fit_brms_lat)
dimnames(posterior_lat)
posterior_lat$parameters
dim(posterior_lat)
head(posterior_lat)

Richness_inflow_lat <- mcmc_intervals(posterior_lat, pars = c("b_Richness_Intercept","b_Richness_temp","b_Richness_log_grav_total",                                                 
                                                      "b_Richness_ClassClosed","b_Richness_ClassRestricted",                                                  
                                                      "b_Richness_log_CorridorIndegreeBR","b_Richness_log_InflowLR",                                                    
                                                      "b_Richness_log_SelfR","b_Richness_log_InflowMPABR","b_Richness_log_InflowNeiBR",
                                                      "b_Richness_Lat","b_Richness_Lon"))
Biomass_inflow_lat <- mcmc_intervals(posterior_lat,c("b_logbiomassarea1_Intercept","b_logbiomassarea1_Richness","b_logbiomassarea1_temp",                                                      
                                             "b_logbiomassarea1_log_grav_total","b_logbiomassarea1_ClassClosed",                                               
                                             "b_logbiomassarea1_ClassRestricted","b_logbiomassarea1_log_CorridorIndegreeBR",                                    
                                             "b_logbiomassarea1_log_InflowLR","b_logbiomassarea1_log_SelfR",                                                 
                                             "b_logbiomassarea1_log_InflowMPABR","b_logbiomassarea1_log_InflowNeiBR","b_logbiomassarea1_Lat","b_logbiomassarea1_Lon"))

require(ggpubr)
SEM_inflow_lat <- ggarrange(Richness_inflow_lat,Biomass_inflow_lat,ncol=1,nrow=2,labels = c("A","B"),align="hv")
ggsave(here("_prelim.figures","SEM_bayes_inflow_coef_latlong.pdf"),SEM_inflow_lat,width=10,height=8)

###### indegree
# with lat/long

rm(species_mod_indegree_lat)
species_mod_indegree_lat <- bf(Richness ~ temp + log_grav_total+ Class  + 
                               log_CorridorIndegreeBR + log_IndegreeBR +log_SelfR + 
                                 log_IndegreeMPABR + log_IndegreeNeiBR + Lat + Lon + (1+ log_grav_total +log_SelfR + Lat + Lon |region))

rm(biom_mod_indegree_lat)
biom_mod_indegree_lat <- bf(log_biomassarea1 ~ Richness+temp + log_grav_total+ Class  + 
                            log_CorridorIndegreeBR + log_IndegreeBR +log_SelfR + 
                              log_IndegreeMPABR + log_IndegreeNeiBR + Lat + Lon +
                            (1+ log_grav_total + log_IndegreeBR +log_SelfR + 
                               log_IndegreeMPABR + log_IndegreeNeiBR + Lat + Lon |region))


all_fit_brms_lat_indegree %>% rm()
all_fit_brms_lat_indegree <-brm(species_mod_indegree_lat + biom_mod_indegree_lat + set_rescor(FALSE), data=TRANSIENT.std,cores=4, chains = 4,
                       iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                       prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

# posterior indegree
posterior_lat_ind_indegree <- as.array(all_fit_brms_lat_indegree)
dimnames(posterior_lat_ind_indegree)

dim(posterior_lat_ind_indegree)
head(posterior_lat_ind_indegree)

Richness_inflow_lat_ind_indegree <- mcmc_intervals(posterior_lat_ind_indegree, pars = c("b_Richness_Intercept","b_Richness_temp","b_Richness_log_grav_total",                                                 
                                                              "b_Richness_ClassClosed","b_Richness_ClassRestricted",                                                  
                                                              "b_Richness_log_CorridorIndegreeBR",                                               
                                                              "b_Richness_log_IndegreeBR",                                                       
                                                              "b_Richness_log_SelfR",                                                            
                                                              "b_Richness_log_IndegreeMPABR",                                                    
                                                              "b_Richness_log_IndegreeNeiBR",
                                                              "b_Richness_Lat","b_Richness_Lon")) +
                                                              ggplot2::labs(title = "Species Richness Posterior distributions",
                                                              subtitle = "with medians and 90% intervals")
Biomass_inflow_lat_ind_indegree <- mcmc_intervals(posterior_lat_ind_indegree,c("b_logbiomassarea1_Intercept","b_logbiomassarea1_Richness","b_logbiomassarea1_temp",                                                      
                                                     "b_logbiomassarea1_log_grav_total","b_logbiomassarea1_ClassClosed",                                               
                                                     "b_logbiomassarea1_ClassRestricted","b_logbiomassarea1_log_CorridorIndegreeBR",                                        
                                                     "b_logbiomassarea1_log_IndegreeBR",                                                
                                                     "b_logbiomassarea1_log_SelfR",                                                     
                                                     "b_logbiomassarea1_log_IndegreeMPABR",                                             
                                                     "b_logbiomassarea1_log_IndegreeNeiBR","b_logbiomassarea1_Lat","b_logbiomassarea1_Lon"))

SEM_inflow_lat_indegree <- ggarrange(Richness_inflow_lat_ind_indegree,Biomass_inflow_lat_ind,ncol=1,nrow=2,labels = c("A","B"),align="hv")
ggsave(here("_prelim.figures","SEM_bayes_inflow_coef_latlong_indegree.pdf"),SEM_inflow_lat_indegree,width=10,height=8)


### marginal plots
library(ggeffects)
plot(conditional_effects(all_fit_brms.tot.TRANSIENT.intr.extr))
conditions <- data.frame(Netflow = c(-0.53, 0.12, 0.77))
plot(conditional_effects(all_fit_brms.tot.TRANSIENT.intr.extr, effects = "log_grav_total:Class:Netflow",conditions=conditions))

conditions <- data.frame(Netflow = c(-0.53, 0.12, 0.77))
plot(conditional_effects(all_fit_brms.tot.TRANSIENT.intr.extr, effects = "log_grav_total:Class:Netflow",
                         conditions = conditions))


# 3-way interactions
Netflow <- c(-0.53, 0.12, 0.77)
Class <- c("Closed","Fished","Restricted")

#model prediction
pred <- expand.grid(Netflow=c(-0.53, 0.12, 0.77),log_grav_total=c(min(TRANSIENT.std$log_grav_total),0,max(TRANSIENT.std$log_grav_total)),Class=factor(Class),
                    log_InflowBR = mean(TRANSIENT.std$log_InflowBR),log_SelfR= mean(TRANSIENT.std$log_SelfR),
                    log_btwdegree = mean(TRANSIENT.std$btwdegree),log_grav_neiBR = mean(TRANSIENT.std$log_grav_neiBR),
                    log_InflowMPABR = mean(TRANSIENT.std$log_InflowMPABR),log_InflowNeiBR = mean(TRANSIENT.std$log_InflowNeiBR),
                    log_CorridorIndegreeBR = mean(TRANSIENT.std$log_CorridorIndegreeBR),region="sw_pacific",
                    Richness = mean(TRANSIENT.std$Richness),temp = mean(TRANSIENT.std$temp))

pred$log_biomassarea <- predict(all_fit_brms.tot.TRANSIENT.intr.extr,pred)
pred$Netflow <- as.factor(pred$Netflow)

matrix(pred$log_biomassarea, ncol = 1)

print(pred$log_biomassarea[,,2][,1])

# carefull, biomass is an array
# Species richness
three.way.Richness <- ggplot(TRANSIENT.std,aes(x=log_grav_total,y=Richness,color=Netflow))+
  geom_point() +
  facet_grid(~Class) +
  geom_line(data=pred,aes(y=log_biomassarea[,,1][,1],x=log_grav_total,group=Netflow))
three.way.Richness  

  # biomass
three.way.Biomass <- ggplot(TRANSIENT.std,aes(x=log_grav_total,y=log_biomassarea,color=Netflow))+
  geom_point() +
  facet_grid(~Class) +
  geom_line(data=pred,aes(y=log_biomassarea[,,2][,1],x=log_grav_total,group=Netflow))
three.way.Biomass 

# carefull, biomass is an array with MPA
# Species richness MPA
three.way.Richness.MPA <- ggplot(TRANSIENT.std,aes(x=log_grav_total,y=Richness,color=Class))+
  geom_point() +
  facet_grid(~Netflow) +
  geom_line(data=pred,aes(y=log_biomassarea[,,1][,1],x=log_grav_total,group=Netflow))
three.way.Richness  

# biomass
three.way.Biomass <- ggplot(TRANSIENT.std,aes(x=log_grav_total,y=log_biomassarea,color=Netflow))+
  geom_point() +
  facet_grid(~Class) +
  geom_line(data=pred,aes(y=log_biomassarea[,,2][,1],x=log_grav_total,group=Netflow))
three.way.Biomass  


  geom_text(label, size = 6, colour = "white")
?geom_line

for(w in 1:3){
  # define the subset of the original data
  dt <- TRANSIENT[TRANSIENT.std$Class == Class[w], ]
  # defining our new data
  nd <- tibble(Class = Class[w], Netflow = Netflow)
  # use our sampling skills, like before
  f <- 
    fitted(all_fit_brms.tot.TRANSIENT.intr.extr, newdata = nd) %>%
    as_tibble() %>%
    bind_cols(nd)
  
  # specify our custom plot
  fig <- 
    ggplot() +
    geom_smooth(data = f,
                aes(x = shade_c, y = Estimate, ymin = Q2.5, ymax = Q97.5),
                stat = "identity", 
                fill = "#CC79A7", color = "#CC79A7", alpha = 1/5, size = 1/2) +
    geom_point(data = dt, 
               aes(x = shade_c, y = blooms),
               shape = 1, color = "#CC79A7") +
    coord_cartesian(xlim = range(d$shade_c), 
                    ylim = range(d$blooms)) +
    scale_x_continuous("Shade (centered)", breaks = c(-1, 0, 1)) +
    labs("Blooms", 
         title = paste("Water (centered) =", w)) +
    theme_pander() + 
    theme(text = element_text(family = "Times"))
  
  # plot that joint
  plot(fig)
}



library(emmeans)
# 1. estimate conditional means
em_ <- emmeans(all_fit_brms.tot.TRANSIENT.intr.extr, ~ log_grav_total + Class + Netflow,
               cov.red = unique)

# 2. estimate the diff of diffs between `visit * Trt` conditional on `zAge`
c_ <- contrast(em_, interaction = c("pairwise", "pairwise"), by = "zAge")

# 3. Plot
emmip(c_, ~ zAge, CIs = TRUE)

plot(conditional_effects(all_fit_brms.tot.RESID.intr.extr))
plot(conditional_effects(all_fit_brms.tot.PARENTAL.intr.extr))
plot(conditional_effects(all_fit_brms.tot.CRYPTIC.intr.extr))






all_fit_brms_lat %>% rm()
all_fit_brms_lat <-brm(species_mod_inflow_lat + biom_mod_inflow_lat + set_rescor(FALSE), data=TRANSIENT.std,cores=4, chains = 4,
                       iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                       prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))


all_fit_brms_lat_un %>% rm()
all_fit_brms_lat_un <-brm(species_mod_inflow_lat + biom_mod_simp + set_rescor(FALSE), data=TRANSIENT.std,cores=4, chains = 4,
                          iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                          prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))


pairs(all_fit_brms)
# with no connectivity
all_fit_brms_simp %>% rm()
all_fit_brms_simp <-brm(species_mod_simp + biom_mod_simp + set_rescor(FALSE), data=TRANSIENT.std,cores=4, chains = 4,
                        iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                        prior = c(prior(normal(0, 10),class = "Intercept"), prior(normal(0, 10), class = "b")))

# with no connectivity
all_fit_brms_NR %>% rm()
all_fit_brms_NR <-brm(biom_mod_NR + set_rescor(FALSE), data=data.std,cores=4, chains = 4,
                      iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.990),
                      prior = c(prior(normal(0, 10),class = "Intercept"), prior(normal(0, 10), class = "b")))


# check the difference between the simple and full model 
LOO(all_fit_brms_lat,all_fit_brms_lat_un,all_fit_brms_simp)
model_weights(all_fit_brms_lat,all_fit_brms_lat_un,all_fit_brms_simp, weights = "loo")

model_weights(all_fit_brms_lat,all_fit_brms_lat_un,all_fit_brms_simp, weights = "waic")

# check divergene
pairs(all_fit_brms)
# the full model is the best model

# summary of the full model
summary(all_fit_brms)
plot(all_fit_brms)
plot(all_fit_brms_FE)
pp_check(all_fit_brms, resp="logbiomassarea1")
pp_check(all_fit_brms, resp="Richness")
bayes_R2(all_fit_brms_lat)
bayes_R2(all_fit_brms_lat_un)
bayes_R2(all_fit_brms_simp)

# with latitude
bayes_R2(all_fit_brms_lat)
pp_check(all_fit_brms_lat, resp="logbiomassarea1")
pp_check(all_fit_brms_lat, resp="Richness")

library(tidyverse)
ggsave(plot=mcmc_plot(all_fit_brms),here("_prelim.figures","SEM_bays_extract.pdf"),width=20,height=30)






#### GLOBAL RESULTS
ALL.FE <- ggarrange(Richness_inflow_lat,Biomass_inflow_lat,ncol=4,nrow=2,labels = c("A","B","C","D","E",'F'),align="hv")
ggsave(here("_Figures","SEM_bayes_coef_latlong.pdf"),SEM_inflow_lat,width=10,height=8)


#####################


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


