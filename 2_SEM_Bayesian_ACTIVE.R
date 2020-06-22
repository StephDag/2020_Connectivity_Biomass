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
all.data<-read.csv(here("_data","FullDataBiomassJune2020.csv"),h=T, stringsAsFactors = F,dec=".")
# clean first column
all.data$X.1 <- NULL
all.data$X <- NULL

# check 
head(all.data)
summary(all.data)
dim(all.data)
summary(all.data)
str(all.data)
apply(all.data,2,class)
# filter data with active

dim(data.std)
summary(data.std)
colnames(all.data)
head(all.data)
str(all.data)
names(all.data)
# list of predictors
rm(PredictVar)
PredictVar<-all.data[,c("ID","region","locality","sites",             
                        "temp","Richness","grav_total","Age_of_protection", 
                        "Class","ModelMode","InflowLR","SelfR",             
                        "Larval_behaviour","FE","InflowBR","IndegreeBR",        
                        "CorridorIndegreeBR","grav_neiBR","IndegreeMPABR","InflowMPABR",       
                        "IndegreeNeiBR","InflowNeiBR","LonID","LatID",             
                        "Lon","Lat","InflowLRBR","biomassarea","crypto","pare","resid","transi")]
PredictVar$log_grav_total <- log(PredictVar$grav_total+1)
PredictVar$log_grav_neiBR <- log(PredictVar$grav_neiBR+1)
PredictVar$log_biomassarea <-log(PredictVar$biomassarea+1)

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

connectivity <- all.data[,c("InflowLR","SelfR","InflowBR","IndegreeBR",        
                         "CorridorIndegreeBR","IndegreeMPABR","InflowMPABR",       
                         "IndegreeNeiBR","InflowNeiBR","InflowLRBR","Outdegree")] 
corr.connectivity <- corrgram(connectivity,order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.cor, text.panel=panel.txt)
ggsave(here("_prelim.figures","Corr_connectivity.pdf"),corr.connectivity,width=20,height=12)

## standrdize x variables
rm(TRANSIENT.std)
TRANSIENT.std<-data.frame(apply(X = TRANSIENT[,c(5,6,11,12,15:22,27:35)], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))
TRANSIENT.std <- cbind(TRANSIENT$region,TRANSIENT.std)
TRANSIENT.std <- cbind(TRANSIENT$Class,TRANSIENT.std)
colnames(TRANSIENT.std)[1:2] <- c("Class","region")
# add log biomass
#TRANSIENT.std$log_biomassarea1<-log(TRANSIENT.std$biomassarea1+1)
TRANSIENT.std$log_InflowLR<-log(TRANSIENT.std$InflowLR+1)
TRANSIENT.std$log_IndegreeBR<-log(TRANSIENT.std$IndegreeBR+1)
TRANSIENT.std$log_SelfR<-log(TRANSIENT.std$SelfR+1)
TRANSIENT.std$log_InflowMPABR<-log(TRANSIENT.std$InflowMPABR+1)
TRANSIENT.std$log_IndegreeMPABR<-log(TRANSIENT.std$IndegreeMPABR+1)
TRANSIENT.std$log_InflowNeiBR<-log(TRANSIENT.std$InflowNeiBR+1)
TRANSIENT.std$log_IndegreeNeiBR<-log(TRANSIENT.std$IndegreeNeiBR+1)
TRANSIENT.std$log_CorridorIndegreeBR <-log(TRANSIENT.std$CorridorIndegreeBR+1)

head(TRANSIENT.std)


# Parental
PARENTAL %>% rm()
PARENTAL <- PredictVar %>% filter(Larval_behaviour == "active" & ModelMode == "pare5") %>% droplevels()

# Fished as the reference
PARENTAL$Class <- relevel(PARENTAL$Class, ref="Fished")
summary(PARENTAL)

library(corrgram)
corrgram(PARENTAL,order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt)

## standrdize x variables
rm(PARENTAL.std)
PARENTAL.std<-data.frame(apply(X = PARENTAL[,c(5,6,11,12,15:22,27:35)], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}))
PARENTAL.std <- cbind(PARENTAL$region,PARENTAL.std)
PARENTAL.std <- cbind(PARENTAL$Class,PARENTAL.std)
colnames(PARENTAL.std)[1:2] <- c("Class","region")
# add log biomass
PARENTAL.std$log_biomassarea1<-log(PARENTAL.std$biomassarea1+1)
PARENTAL.std$log_InflowLR<-log(PARENTAL.std$InflowLR+1)
PARENTAL.std$log_IndegreeBR<-log(PARENTAL.std$IndegreeBR+1)
PARENTAL.std$log_SelfR<-log(PARENTAL.std$SelfR+1)
PARENTAL.std$log_InflowMPABR<-log(PARENTAL.std$InflowMPABR+1)
PARENTAL.std$log_IndegreeMPABR<-log(PARENTAL.std$IndegreeMPABR+1)
PARENTAL.std$log_InflowNeiBR<-log(PARENTAL.std$InflowNeiBR+1)
PARENTAL.std$log_IndegreeNeiBR<-log(PARENTAL.std$IndegreeNeiBR+1)
PARENTAL.std$log_CorridorIndegreeBR <-log(PARENTAL.std$CorridorIndegreeBR+1)

head(PARENTAL.std)

# Cryptic
CRYPTIC %>% rm()
CRYPTIC <- PredictVar %>% filter(Larval_behaviour == "active" & ModelMode == "crypto15") %>% droplevels()

# Fished as the reference
CRYPTIC$Class <- relevel(CRYPTIC$Class, ref="Fished")
summary(CRYPTIC)

library(corrgram)
corrgram(CRYPTIC,order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt)

## standrdize x variables
rm(CRYPTIC.std)
CRYPTIC.std<-data.frame(apply(X = CRYPTIC[,c(5,6,11,12,15:22,27:35)], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}))
CRYPTIC.std <- cbind(CRYPTIC$region,CRYPTIC.std)
CRYPTIC.std <- cbind(CRYPTIC$Class,CRYPTIC.std)
colnames(CRYPTIC.std)[1:2] <- c("Class","region")
# add log biomass
#CRYPTIC.std$log_biomassarea1<-log(CRYPTIC.std$biomassarea1+1)
CRYPTIC.std$log_InflowLR<-log(CRYPTIC.std$InflowLR+1)
CRYPTIC.std$log_IndegreeBR<-log(CRYPTIC.std$IndegreeBR+1)
CRYPTIC.std$log_SelfR<-log(CRYPTIC.std$SelfR+1)
CRYPTIC.std$log_InflowMPABR<-log(CRYPTIC.std$InflowMPABR+1)
CRYPTIC.std$log_IndegreeMPABR<-log(CRYPTIC.std$IndegreeMPABR+1)
CRYPTIC.std$log_InflowNeiBR<-log(CRYPTIC.std$InflowNeiBR+1)
CRYPTIC.std$log_IndegreeNeiBR<-log(CRYPTIC.std$IndegreeNeiBR+1)
CRYPTIC.std$log_CorridorIndegreeBR <-log(CRYPTIC.std$CorridorIndegreeBR+1)

head(CRYPTIC.std)

# Resident
RESID %>% rm()
RESID <- PredictVar %>% filter(Larval_behaviour == "active" & ModelMode == "resid15") %>% droplevels()

# Fished as the reference
RESID$Class <- relevel(RESID$Class, ref="Fished")
summary(RESID)

library(corrgram)
corrgram(RESID,order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt)

## standrdize x variables
rm(RESID.std)
RESID.std<-data.frame(apply(X = RESID[,c(5,6,11,12,15:22,27:35)], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}))
RESID.std <- cbind(RESID$region,RESID.std)
RESID.std <- cbind(RESID$Class,RESID.std)
colnames(RESID.std)[1:2] <- c("Class","region")
# add log biomass

RESID.std$log_InflowLR<-log(RESID.std$InflowLR+1)
RESID.std$log_IndegreeBR<-log(RESID.std$IndegreeBR+1)
RESID.std$log_SelfR<-log(RESID.std$SelfR+1)
RESID.std$log_InflowMPABR<-log(RESID.std$InflowMPABR+1)
RESID.std$log_IndegreeMPABR<-log(RESID.std$IndegreeMPABR+1)
RESID.std$log_InflowNeiBR<-log(RESID.std$InflowNeiBR+1)
RESID.std$log_IndegreeNeiBR<-log(RESID.std$IndegreeNeiBR+1)
RESID.std$log_CorridorIndegreeBR <-log(RESID.std$CorridorIndegreeBR+1)

head(RESID.std)

# Adjancy matrix to include spatial autocorrelation
library(spdep)
library(sp)

TRANSIENT.std <- cbind(TRANSIENT.std,TRANSIENT$Lat,TRANSIENT$Lon)
colnames(TRANSIENT.std)[c(35:36)] <- c("Latitude","Longitude")
dim(TRANSIENT.std)
head(TRANSIENT.std)

TRANSIENT_sf <- st_as_sf(TRANSIENT.std, coords = c("Longitude", "Latitude"))
str(TRANSIENT_sf)
class(TRANSIENT_sf)
coord(TRANSIENT_sf)
ls(TRANSIENT_sf)

library(raster)
library(ecodist)
pts <- TRANSIENT.std[,c("Longitude", "Latitude")]
pts.sf <- SpatialPoints(pts)
head(pts)
summary(pts)

pts.dist <- knn2nb(knearneigh(pts.sf, k = 1))

library(sp)
rm(grid)
grid <- expand.grid(Longitude = seq(-179,178, l = 100),
                     Latitude = seq(-32, 27, l = 100))
head(grid)
K <- nrow(grid)

# set up distance and neighbourhood matrices
distance <- as.matrix(dist(grid))
W <- array(0, c(K, K))
W[distance == 1] <- 1 


gdis <- pointDistance(pts, lonlat=TRUE,allpairs=T)
gdis
dim(gdis)
gdis <- as.matrix(gdis)
full(gdis)

W <- 1 / gdis
W.ad <- nb2mat(W,style="B")
W.ad
head(W.ad)
round(W, 4)
head(W)[1:5,1:5]
class(W)

car(W)


# built SEM model
# scaled 
  # gravity of neightbour is highly positively correlated with connectivity variables
rm(species_mod_inflow) # ENV + CON
species_mod_inflow <- bf(Richness ~ temp + log_grav_total + Class  + 
                    log_InflowLR +log_SelfR + log_CorridorIndegreeBR +(1+ log_grav_total +log_SelfR +log_CorridorIndegreeBR |region))

rm(species_mod_inflow.nocor) # ENV + CON
species_mod_inflow.nocor <- bf(Richness ~ temp + log_grav_total + Class  + 
                           log_InflowLR +log_SelfR +(1+ log_grav_total +log_SelfR  |region))

rm(biom_mod_inflow.tot) # ENV + CON +S 
biom_mod_inflow.tot <- bf(log_biomassarea ~ Richness+temp + log_grav_total + Class  + 
                        log_InflowLR +log_SelfR + log_CorridorIndegreeBR +
                        (1+ log_grav_total + log_InflowLR +  log_SelfR +log_CorridorIndegreeBR  |region))

rm(biom_mod_inflow.tot.nocor) # ENV + CON +S 
biom_mod_inflow.tot.nocor <- bf(log_biomassarea ~ Richness+temp + log_grav_total + Class  + 
                            log_InflowLR +log_SelfR +
                            (1+ log_grav_total + log_InflowLR +  log_SelfR |region))

rm(biom_mod_inflow.simp) # ENV + CON
biom_mod_inflow.simp <- bf(log_biomassarea ~ temp + log_grav_total + Class  + 
                             log_InflowLR +log_SelfR + log_CorridorIndegreeBR 
                             (1+ log_grav_total + log_InflowLR + log_SelfR +log_CorridorIndegreeBR  |region))

rm(biom_mod_nocon) # ENV
biom_mod_nocon <- bf(log_biomassarea ~ temp + log_grav_total + Class  + (1+ log_grav_total|region))

rm(biom_mod_nocon.S) # ENV + S
biom_mod_nocon.S <- bf(log_biomassarea ~ Richness+temp + log_grav_total + Class  + (1+ log_grav_total|region))

rm(S_mod_nocon) # ENV
S_mod_nocon <- bf(Richness ~ temp + log_grav_total + Class  + (1+ log_grav_total|region))

# with lat/long
rm(species_mod_inflow_lat)
species_mod_inflow_lat <- bf(Richness ~ temp + log_grav_total+ Class  + 
                           log_CorridorIndegreeBR + log_InflowLR +log_SelfR + 
                           log_InflowMPABR + log_InflowNeiBR + Latitude + Longitude + (1+ log_grav_total +log_SelfR + Lat + Lon |region))

rm(biom_mod_inflow_lat)
biom_mod_inflow_lat <- bf(log_biomassarea ~ Richness+temp + log_grav_total+ Class  + 
                        log_CorridorIndegreeBR + log_InflowLR +log_SelfR +  
                        log_InflowMPABR + log_InflowNeiBR + Latitude + Longitude +
                        (1+ log_grav_total + log_InflowLR + log_InflowMPABR + log_InflowNeiBR + log_SelfR + Lat + Lon |region))


rm(species_mod_simp)
species_mod_simp <- bf(Richness ~ temp + log_grav_total+Class + Latitude + Longitude+ (1+ log_grav_total |region)) 

rm(biom_mod_simp)
biom_mod_simp <- bf(log_biomassarea ~ Richness+temp + log_grav_total+Class+ Latitude + Longitude+ (1+ log_grav_total + Richness |region))


# SEM
## scaled
    # all_fit_brms = connectivity explains both biomass and species richness
#TRANSIENT
all_fit_brms.tot.TRANSIENT %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.TRANSIENT <-brm(species_mod_inflow + biom_mod_inflow.tot + set_rescor(FALSE), data=TRANSIENT.std,cores=4,chains = 4,
                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                   prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

all_fit_brms.tot.TRANSIENT.nocor %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.TRANSIENT.nocor <-brm(species_mod_inflow.nocor + biom_mod_inflow.tot.nocor + set_rescor(FALSE), data=TRANSIENT.std,cores=4,chains = 4,
                                 iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                 prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

all_fit_brms.nocon.S.TRANSIENT %>% rm() # connectivity only through S + ENV
all_fit_brms.nocon.S.TRANSIENT <-brm(species_mod_inflow + biom_mod_nocon.S + set_rescor(FALSE), data=TRANSIENT.std,cores=4,chains = 4,
                                     iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                     prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

all_fit_brms.nocon.TRANSIENT %>% rm() # no connectivity for S and B, only environment
all_fit_brms.nocon.TRANSIENT <-brm(S_mod_nocon + biom_mod_nocon.S + set_rescor(FALSE), data=TRANSIENT.std,cores=4,chains = 4,
                                  iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                  prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

TRANSIENT.weight <- model_weights(all_fit_brms.tot.TRANSIENT,
              all_fit_brms.nocon.S.TRANSIENT,
              all_fit_brms.nocon.TRANSIENT,weights = "loo")

TRANSIENT.LOO <- LOO(all_fit_brms.tot.TRANSIENT,
              all_fit_brms.nocon.S.TRANSIENT,
              all_fit_brms.nocon.TRANSIENT)
WAIC(all_fit_brms.tot.TRANSIENT,
    all_fit_brms.nocon.S.TRANSIENT,
    all_fit_brms.nocon.TRANSIENT)

TRANSIENT.R2 <- rbind(bayes_R2(all_fit_brms.tot.TRANSIENT),
bayes_R2(all_fit_brms.nocon.S.TRANSIENT),
bayes_R2(all_fit_brms.nocon.TRANSIENT))

mcmc_plot(all_fit_brms.nocon.S.TRANSIENT)
mcmc_plot(all_fit_brms.tot.TRANSIENT)

library("bayesplot")
library("ggplot2")
library("rstanarm")   
require(ggpubr)

posterior.TRANSIENT.tot <- as.array(all_fit_brms.tot.TRANSIENT)
posterior.TRANSIENT.simp <- as.array(all_fit_brms.simp.TRANSIENT)

# RESID
all_fit_brms.tot.RESID %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.RESID <-brm(species_mod_inflow + biom_mod_inflow.tot + set_rescor(FALSE), data=RESID.std,cores=4,chains = 4,
                                 iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                 prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

all_fit_brms.tot.RESID.nocor %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.RESID.nocor <-brm(species_mod_inflow.nocor + biom_mod_inflow.tot.nocor + set_rescor(FALSE), data=RESID.std,cores=4,chains = 4,
                                       iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                       prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

all_fit_brms.nocon.S.RESID %>% rm() # connectivity only through S + ENV
all_fit_brms.nocon.S.RESID <-brm(species_mod_inflow + biom_mod_nocon.S + set_rescor(FALSE), data=RESID.std,cores=4,chains = 4,
                                     iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                     prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

all_fit_brms.nocon.RESID %>% rm() # no connectivity for S and B, only environment
all_fit_brms.nocon.RESID <-brm(S_mod_nocon + biom_mod_nocon.S + set_rescor(FALSE), data=RESID.std,cores=4,chains = 4,
                                   iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                   prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

RESID.weight <- model_weights(all_fit_brms.tot.RESID,
              all_fit_brms.nocon.S.RESID,
              all_fit_brms.nocon.RESID,weights = "loo")
model_weights(all_fit_brms.tot.RESID,
              all_fit_brms.nocon.S.RESID,
              all_fit_brms.nocon.RESID,weights = "WAIC")

RESID.LOO <- LOO(all_fit_brms.tot.RESID,
    all_fit_brms.nocon.S.RESID,
    all_fit_brms.nocon.RESID)
WAIC(all_fit_brms.tot.RESID,
     all_fit_brms.nocon.S.RESID,
     all_fit_brms.nocon.RESID)
mcmc_plot(all_fit_brms.tot.RESID)

RESID.R2 <- rbind(bayes_R2(all_fit_brms.tot.RESID),
bayes_R2(all_fit_brms.nocon.S.RESID),
bayes_R2(all_fit_brms.nocon.RESID))

posterior.RESID.tot <- as.array(all_fit_brms.tot.RESID)

# PARENTAL
all_fit_brms.tot.PARENTAL %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.PARENTAL <-brm(species_mod_inflow + biom_mod_inflow.tot + set_rescor(FALSE), data=PARENTAL.std,cores=4,chains = 4,
                             iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                             prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

all_fit_brms.tot.PARENTAL.nocor %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.PARENTAL.nocor <-brm(species_mod_inflow.nocor + biom_mod_inflow.tot.nocor + set_rescor(FALSE), data=PARENTAL.std,cores=4,chains = 4,
                                       iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                       prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))


all_fit_brms.nocon.S.PARENTAL %>% rm() # connectivity only through S + ENV
all_fit_brms.nocon.S.PARENTAL <-brm(species_mod_inflow + biom_mod_nocon.S + set_rescor(FALSE), data=PARENTAL.std,cores=4,chains = 4,
                                 iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                 prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

all_fit_brms.nocon.PARENTAL %>% rm() # no connectivity for S and B, only environment
all_fit_brms.nocon.PARENTAL <-brm(S_mod_nocon + biom_mod_nocon.S + set_rescor(FALSE), data=PARENTAL.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

PARENTAL.weight <- model_weights(all_fit_brms.tot.PARENTAL,
              all_fit_brms.nocon.S.PARENTAL,
              all_fit_brms.nocon.PARENTAL,weights = "loo")

PARENTAL.LOO <- LOO(all_fit_brms.tot.PARENTAL,
    all_fit_brms.nocon.S.PARENTAL,
    all_fit_brms.nocon.PARENTAL)
WAIC(all_fit_brms.tot.PARENTAL,
     all_fit_brms.nocon.S.PARENTAL,
     all_fit_brms.nocon.PARENTAL)
mcmc_plot(all_fit_brms.tot.PARENTAL)

PARENTAL.R2 <- rbind(bayes_R2(all_fit_brms.tot.PARENTAL),
bayes_R2(all_fit_brms.nocon.S.PARENTAL),
bayes_R2(all_fit_brms.nocon.PARENTAL))

posterior.PARENTAL.tot <- as.array(all_fit_brms.tot.PARENTAL)
posterior.PARENTAL.simp <- as.array(all_fit_brms.simp.PARENTAL)

# CRYPTIC
all_fit_brms.tot.CRYPTIC %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.CRYPTIC <-brm(species_mod_inflow + biom_mod_inflow.tot + set_rescor(FALSE), data=CRYPTIC.std,cores=4,chains = 4,
                                iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

all_fit_brms.tot.CRYPTIC.nocor %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.CRYPTIC.nocor <-brm(species_mod_inflow.nocor + biom_mod_inflow.tot.nocor + set_rescor(FALSE), data=CRYPTIC.std,cores=4,chains = 4,
                                       iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                       prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

all_fit_brms.nocon.S.CRYPTIC %>% rm() # connectivity only through S + ENV
all_fit_brms.nocon.S.CRYPTIC <-brm(species_mod_inflow + biom_mod_nocon.S + set_rescor(FALSE), data=CRYPTIC.std,cores=4,chains = 4,
                                    iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                    prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

all_fit_brms.nocon.CRYPTIC %>% rm() # no connectivity for S and B, only environment
all_fit_brms.nocon.CRYPTIC <-brm(S_mod_nocon + biom_mod_nocon.S + set_rescor(FALSE), data=CRYPTIC.std,cores=4,chains = 4,
                                  iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                  prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

CRYPTIC.weight <- model_weights(all_fit_brms.tot.CRYPTIC,
              all_fit_brms.nocon.S.CRYPTIC,
              all_fit_brms.nocon.CRYPTIC,weights = "loo")



bayes_R2(all_fit_brms.tot.CRYPTIC.nocor)
mcmc_plot(all_fit_brms.tot.CRYPTIC.nocor)

CRYPTIC.LOO <- LOO(all_fit_brms.tot.CRYPTIC,
    all_fit_brms.nocon.S.CRYPTIC,
    all_fit_brms.nocon.CRYPTIC,all_fit_brms.tot.CRYPTIC.nocor)

WAIC(all_fit_brms.tot.CRYPTIC,
     all_fit_brms.nocon.S.CRYPTIC,
     all_fit_brms.nocon.CRYPTIC,all_fit_brms.tot.CRYPTIC.nocor)
mcmc_plot(all_fit_brms.tot.CRYPTIC)
mcmc_plot(all_fit_brms.nocon.S.CRYPTIC)

CRYPTIC.R2 <-rbind(bayes_R2(all_fit_brms.tot.CRYPTIC),
bayes_R2(all_fit_brms.nocon.S.CRYPTIC),
bayes_R2(all_fit_brms.nocon.CRYPTIC))

posterior.CRYPTIC.tot <- as.array(all_fit_brms.tot.CRYPTIC)
posterior.CRYPTIC.simp <- as.array(all_fit_brms.simp.CRYPTIC)

#
library("bayesplot")
library("ggplot2")
library("rstanarm")   
require(ggpubr)

# TRANSIENT TOT
Richness_inflow_TRANSIENT_tot <- mcmc_intervals(posterior.TRANSIENT.tot, pars = c("b_Richness_Intercept","b_Richness_temp","b_Richness_log_grav_total",                                                 
                                   "b_Richness_ClassClosed","b_Richness_ClassRestricted",                                                  
                                   "b_Richness_log_CorridorIndegreeBR","b_Richness_log_InflowLR",                                                    
                                   "b_Richness_log_SelfR"))
Biomass_inflow_TRANSIENT_tot <- mcmc_intervals(posterior.TRANSIENT.tot,c("b_logbiomassarea_Intercept","b_logbiomassarea_Richness","b_logbiomassarea_temp",                                                      
                          "b_logbiomassarea_log_grav_total","b_logbiomassarea_ClassClosed",                                               
                          "b_logbiomassarea_ClassRestricted","b_logbiomassarea_log_CorridorIndegreeBR",                                    
                          "b_logbiomassarea_log_InflowLR","b_logbiomassarea_log_SelfR"))

# RESIDENT
Richness_inflow_RESIDENT_tot <- mcmc_intervals(posterior.RESID.tot, pars = c("b_Richness_Intercept","b_Richness_temp","b_Richness_log_grav_total",                                                 
                                                                                  "b_Richness_ClassClosed","b_Richness_ClassRestricted",                                                  
                                                                                  "b_Richness_log_CorridorIndegreeBR","b_Richness_log_InflowLR",                                                    
                                                                                  "b_Richness_log_SelfR"))
Biomass_inflow_RESIDENT_tot <- mcmc_intervals(posterior.RESID.tot,c("b_logbiomassarea_Intercept","b_logbiomassarea_Richness","b_logbiomassarea_temp",                                                      
                                                                         "b_logbiomassarea_log_grav_total","b_logbiomassarea_ClassClosed",                                               
                                                                         "b_logbiomassarea_ClassRestricted","b_logbiomassarea_log_CorridorIndegreeBR",                                    
                                                                         "b_logbiomassarea_log_InflowLR","b_logbiomassarea_log_SelfR"))
# PARENTAL
Richness_inflow_PARENTAL_tot <- mcmc_intervals(posterior.PARENTAL.tot, pars = c("b_Richness_Intercept","b_Richness_temp","b_Richness_log_grav_total",                                                 
                                                                             "b_Richness_ClassClosed","b_Richness_ClassRestricted",                                                  
                                                                             "b_Richness_log_CorridorIndegreeBR","b_Richness_log_InflowLR",                                                    
                                                                             "b_Richness_log_SelfR"))
Biomass_inflow_PARENTAL_tot <- mcmc_intervals(posterior.PARENTAL.tot,c("b_logbiomassarea_Intercept","b_logbiomassarea_Richness","b_logbiomassarea_temp",                                                      
                                                                    "b_logbiomassarea_log_grav_total","b_logbiomassarea_ClassClosed",                                               
                                                                    "b_logbiomassarea_ClassRestricted","b_logbiomassarea_log_CorridorIndegreeBR",                                    
                                                                    "b_logbiomassarea_log_InflowLR","b_logbiomassarea_log_SelfR"))
# CRYPTIC
Richness_inflow_CRYPTIC_tot <- mcmc_intervals(posterior.CRYPTIC.tot, pars = c("b_Richness_Intercept","b_Richness_temp","b_Richness_log_grav_total",                                                 
                                                                                "b_Richness_ClassClosed","b_Richness_ClassRestricted",                                                  
                                                                                "b_Richness_log_CorridorIndegreeBR","b_Richness_log_InflowLR",                                                    
                                                                                "b_Richness_log_SelfR"))
Biomass_inflow_CRYPTIC_tot <- mcmc_intervals(posterior.CRYPTIC.tot,c("b_logbiomassarea_Intercept","b_logbiomassarea_Richness","b_logbiomassarea_temp",                                                      
                                                                       "b_logbiomassarea_log_grav_total","b_logbiomassarea_ClassClosed",                                               
                                                                       "b_logbiomassarea_ClassRestricted","b_logbiomassarea_log_CorridorIndegreeBR",                                    
                                                                       "b_logbiomassarea_log_InflowLR","b_logbiomassarea_log_SelfR"))

SEM_inflow_tot <- ggarrange(Richness_inflow_TRANSIENT_tot,Richness_inflow_RESIDENT_tot,
                            Richness_inflow_PARENTAL_tot,Richness_inflow_CRYPTIC_tot, 
                            Biomass_inflow_TRANSIENT_tot,
                            Biomass_inflow_PARENTAL_tot,
                            Biomass_inflow_CRYPTIC_tot,
                            Biomass_inflow_RESIDENT_tot,
                            ncol=4,nrow=2,labels = c("A_Richness_Transient","B_Richness_Resident","C_Richness_Parental","D_Richness_Cryptic",
                                                     "E_Biomass_Transient","F_Biomass_Resident","G_Biomass_Parental","H_Biomass_Cryptic"),align="hv")
ggsave(here("_prelim.figures","SEM_bayes_inflow_coef_total_FE.pdf"),SEM_inflow_tot,width=20,height=12)

SEM_inflow_tot_weight <- round(rbind(TRANSIENT.weight,RESID.weight,PARENTAL.weight,CRYPTIC.weight),2)
colnames(SEM_inflow_tot_weight) <- c("Model1_Full","Model2_Smed","Model3_noCon")
rownames(SEM_inflow_tot_weight) <- c("TRANSIENT","RESIDENT","PARENTAL","CRYPTIC")

require(kableExtra)
SEM_inflow_tot_weight.table <- SEM_inflow_tot_weight  %>%
kable(align="c") %>%
  kable_styling() %>%
  save_kable("_prelim.figures/SEM_bayes_inflow_weight_FE.png")

SEM_inflow_tot_R2 <- round(rbind(TRANSIENT.R2[,1],RESID.R2[,1],PARENTAL.R2[,1],CRYPTIC.R2[,1]),2)
colnames(SEM_inflow_tot_R2) <- c("Model1_Full_S","Model1_Full_B","Model2_Smed_S","Model2_Smed_B","Model3_noCon_S","Model3_noCon_B")
rownames(SEM_inflow_tot_R2) <- c("TRANSIENT",
                                     "RESIDENT",
                                     "PARENTAL",
                                     "CRYPTIC")
SEM_inflow_tot_R2.table <- SEM_inflow_tot_R2  %>%
  kable(align="c") %>%
  kable_styling() %>%
  save_kable("_prelim.figures/SEM_bayes_inflow_R2_FE.png")












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


