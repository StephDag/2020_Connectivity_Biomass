# SEM analysis - Global connectivity 
# author: Steph D'agata
# date: April 2020
# outputs: Power of the coefficients

rm(list=ls())

# packages
library(dotwhisker)
library(broom)
library(dplyr)
library("bayesplot")
library("ggplot2")
library("rstanarm")   
require(ggpubr)

# load models
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
biom_mod_inflow.intr.extr <- bf(log_biomassarea ~ Richness + temp + log_grav_total*Class*Netflow+ 
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

# load original data
  # use the 0.Load_DATA_biomass.R script to load data

# full data for each functional group
all_fit_brms.tot.TRANSIENT.intr.extr %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.TRANSIENT.intr.extr <-brm(species_mod_inflow.intr.extr + biom_mod_inflow.intr.extr + set_rescor(FALSE), data=TRANSIENT.std,cores=24,chains = 4,
                                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))

all_fit_brms.tot.RESID.intr.extr %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.RESID.intr.extr <-brm(species_mod_inflow.intr.extr + biom_mod_inflow.intr.extr + set_rescor(FALSE), data=RESID.std,cores=24,chains = 4,
                                       iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                       prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
all_fit_brms.tot.PARENTAL.intr.extr %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.PARENTAL.intr.extr <-brm(species_mod_inflow.intr.extr + biom_mod_inflow.intr.extr + set_rescor(FALSE), data=PARENTAL.std,cores=24,chains = 4,
                                          iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                          prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
all_fit_brms.tot.CRYPTIC.intr.extr %>% rm()  # Full - connectivity through both S and B
all_fit_brms.tot.CRYPTIC.intr.extr <-brm(species_mod_inflow.intr.extr + biom_mod_inflow.intr.extr + set_rescor(FALSE), data=CRYPTIC.std,cores=24,chains = 4,
                                         iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                         prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
# newdata 
post.pred.TRANSIENT <- posterior_predict(all_fit_brms.tot.TRANSIENT.intr.extr) # 8000*256*2 (fish richness and biomass)
post.pred.RESID <- posterior_predict(all_fit_brms.tot.RESID.intr.extr) # 8000*256*2 (fish richness and biomass)
post.pred.PARENTAL <- posterior_predict(all_fit_brms.tot.PARENTAL.intr.extr) # 8000*256*2 (fish richness and biomass)
post.pred.CRYPTIC <- posterior_predict(all_fit_brms.tot.CRYPTIC.intr.extr) # 8000*256*2 (fish richness and biomass)

# Number of posterior to draw from each matrix = 20% = 1600
post.draws <- sample(1:8000,1600) # sample randomly 1600 draws from the 8000 from each model

# post.pred: 1st list = Richness
Rich.TRANSIENT <- post.pred.TRANSIENT[,,1] 
Biom.TRANSIENT <- post.pred.TRANSIENT[,,2]
Rich.RESID <- post.pred.RESID[,,1] 
Biom.RESID <- post.pred.RESID[,,2]
Rich.PARENTAL <- post.pred.PARENTAL[,,1] 
Biom.PARENTAL <- post.pred.PARENTAL[,,2]
Rich.CRYPTIC <- post.pred.CRYPTIC[,,1] 
Biom.CRYPTIC <- post.pred.CRYPTIC[,,2]

rm(newdata.TRANSIENT)
newdata.TRANSIENT <- TRANSIENT.std[-which(is.na(TRANSIENT.std$Netflow)),] # rm NA
rm(newdata.RESID)
newdata.RESID <- RESID.std[-which(is.na(RESID.std$Netflow)),] # rm NA
rm(newdata.PARENTAL)
newdata.PARENTAL <- PARENTAL.std[-which(is.na(PARENTAL.std$Netflow)),] # rm NA
rm(newdata.CRYPTIC)
newdata.CRYPTIC <- CRYPTIC.std[-which(is.na(CRYPTIC.std$Netflow)),] # rm NA

# updates the model with each raw of post.pred and extract the coefficient values
  # create list to store the coefficient for each of the 30 model
  rm(coef.rich.TRANSIENT); rm(coef.biom.TRANSIENT)
  coef.rich.TRANSIENT <- list()
  coef.biom.TRANSIENT <- list()

  rm(coef.rich.RESID); rm(coef.biom.RESID)
  coef.rich.RESID <- list()
  coef.biom.RESID <- list()
  
  rm(coef.rich.PARENTAL); rm(coef.biom.PARENTAL)
  coef.rich.PARENTAL <- list()
  coef.biom.PARENTAL <- list()
  
  rm(coef.rich.CRYPTIC); rm(coef.biom.CRYPTIC)
  coef.rich.CRYPTIC <- list()
  coef.biom.CRYPTIC <- list()
  
for (i in 1:length(post.draws)){
  paste("i=",i,sep="_")
  
  Todraw <- post.draws[i]
  
  # TRANSIENT
  newdata.TRANSIENT$Richness <- Rich.TRANSIENT[Todraw,]
  newdata.TRANSIENT$log_biomassarea <- Biom.TRANSIENT[Todraw,]
  temp.fit.TRANSIENT <- update(all_fit_brms.tot.TRANSIENT.intr.extr, newdata=newdata.TRANSIENT)
  # extract coefficients for each TRANSIENT model
  a.TRANSIENT <- mcmc_intervals(temp.fit.TRANSIENT)
  rich.var.TRANSIENT <- as.data.frame(a.TRANSIENT$data)[grep("b_Richness",as.data.frame(a.TRANSIENT$data)[,"parameter"]),"parameter"]
  biom.var.TRANSIENT <- as.data.frame(a.TRANSIENT$data)[grep("b_logbiomassarea",as.data.frame(a.TRANSIENT$data)[,"parameter"]),"parameter"]

  coef.rich.TRANSIENT[[Todraw]] <- mcmc_intervals_data(temp.fit.TRANSIENT, pars = as.character(rich.var.TRANSIENT))
  coef.biom.TRANSIENT[[Todraw]]  <- mcmc_intervals_data(temp.fit.TRANSIENT,pars = as.character(biom.var.TRANSIENT))
  
  # RESIDENT
  newdata.RESID$Richness <- Rich.RESID[Todraw,]
  newdata.RESID$log_biomassarea <- Biom.RESID[Todraw,]
  temp.fit.RESID <- update(all_fit_brms.tot.RESID.intr.extr, newdata=newdata.RESID)
  # extract coefficients for each TRANSIENT model
  a.RESID <- mcmc_intervals(temp.fit.RESID)
  rich.var.RESID <- as.data.frame(a.RESID$data)[grep("b_Richness",as.data.frame(a.RESID$data)[,"parameter"]),"parameter"]
  biom.var.RESID <- as.data.frame(a.RESID$data)[grep("b_logbiomassarea",as.data.frame(a.RESID$data)[,"parameter"]),"parameter"]
  
  coef.rich.RESID[[Todraw]] <- mcmc_intervals_data(temp.fit.RESID, pars = as.character(rich.var.RESID))
  coef.biom.RESID[[Todraw]]  <- mcmc_intervals_data(temp.fit.RESID,pars = as.character(biom.var.RESID))
  
  # PARENTAL
  newdata.PARENTAL$Richness <- Rich.PARENTAL[Todraw,]
  newdata.PARENTAL$log_biomassarea <- Biom.PARENTAL[Todraw,]
  temp.fit.PARENTAL <- update(all_fit_brms.tot.PARENTAL.intr.extr, newdata=newdata.PARENTAL)
  # extract coefficients for each TRANSIENT model
  a.PARENTAL <- mcmc_intervals(temp.fit.PARENTAL)
  rich.var.PARENTAL <- as.data.frame(a.PARENTAL$data)[grep("b_Richness",as.data.frame(a.PARENTAL$data)[,"parameter"]),"parameter"]
  biom.var.PARENTAL <- as.data.frame(a.PARENTAL$data)[grep("b_logbiomassarea",as.data.frame(a.PARENTAL$data)[,"parameter"]),"parameter"]
  
  coef.rich.PARENTAL[[Todraw]] <- mcmc_intervals_data(temp.fit.PARENTAL, pars = as.character(rich.var.PARENTAL))
  coef.biom.PARENTAL[[Todraw]]  <- mcmc_intervals_data(temp.fit.PARENTAL,pars = as.character(biom.var.PARENTAL))
  
  # CRYPTIC
  newdata.CRYPTIC$Richness <- Rich.CRYPTIC[Todraw,]
  newdata.CRYPTIC$log_biomassarea <- Biom.CRYPTIC[Todraw,]
  temp.fit.CRYPTIC <- update(all_fit_brms.tot.CRYPTIC.intr.extr, newdata=newdata.CRYPTIC)
  # extract coefficients for each TRANSIENT model
  a.CRYPTIC <- mcmc_intervals(temp.fit.CRYPTIC)
  rich.var.CRYPTIC <- as.data.frame(a.CRYPTIC$data)[grep("b_Richness",as.data.frame(a.CRYPTIC$data)[,"parameter"]),"parameter"]
  biom.var.CRYPTIC <- as.data.frame(a.CRYPTIC$data)[grep("b_logbiomassarea",as.data.frame(a.CRYPTIC$data)[,"parameter"]),"parameter"]
  
  coef.rich.CRYPTIC[[Todraw]] <- mcmc_intervals_data(temp.fit.CRYPTIC, pars = as.character(rich.var.CRYPTIC))
  coef.biom.CRYPTIC[[Todraw]]  <- mcmc_intervals_data(temp.fit.CRYPTIC,pars = as.character(biom.var.CRYPTIC))
  
  # clean temporary object at the end of the loop
  rm(temp.fit.TRANSIENT); rm(temp.fit.RESID);rm(temp.fit.PARENTAL);rm(temp.fit.CRYPTIC)
  rm(a.TRANSIENT); rm(a.RESID); rm(a.PARENTAL); rm(a.CRYPTIC)
  rm(rich.var.TRANSIENT);  rm(rich.var.RESID);  rm(rich.var.PARENTAL);  rm(rich.var.CRYPTIC)
  rm(biom.var.TRANSIENT);  rm(biom.var.RESID);  rm(biom.var.PARENTAL);  rm(biom.var.CRYPTIC)

}

rm(temp.rich.TRANSIENT); rm(temp.biom.TRANSIENT)
temp.rich.TRANSIENT <-  coef.rich.TRANSIENT[[1]]$parameter
temp.biom.TRANSIENT <-  coef.biom.TRANSIENT[[1]]$parameter
rm(temp.rich.RESID); rm(temp.biom.RESID)
temp.rich.RESID <-  coef.rich.RESID[[1]]$parameter
temp.biom.RESID <-  coef.biom.RESID[[1]]$parameter
rm(temp.rich.PARENTAL); rm(temp.biom.TRANSIENT)
temp.rich.PARENTAL <-  coef.rich.PARENTAL[[1]]$parameter
temp.biom.PARENTAL <-  coef.biom.PARENTAL[[1]]$parameter
rm(temp.rich.CRYPTIC); rm(temp.biom.CRYPTIC)
temp.rich.CRYPTIC <-  coef.rich.CRYPTIC[[1]]$parameter
temp.biom.CRYPTIC <-  coef.biom.CRYPTIC[[1]]$parameter

for (j in 1:length(post.draws)){
  temp.rich.TRANSIENT <- cbind(temp.rich.TRANSIENT,as.vector(as.data.frame(coef.rich.TRANSIENT[[j]][,"m"])))
  temp.biom.TRANSIENT <- cbind(temp.biom.TRANSIENT,as.vector(as.data.frame(coef.biom.TRANSIENT[[j]][,"m"])))
  temp.rich.PARENTAL <- cbind(temp.rich.PARENTAL,as.vector(as.data.frame(coef.rich.PARENTAL[[j]][,"m"])))
  temp.biom.PARENTAL <- cbind(temp.biom.PARENTAL,as.vector(as.data.frame(coef.biom.PARENTAL[[j]][,"m"])))
  temp.rich.RESID <- cbind(temp.rich.RESID,as.vector(as.data.frame(coef.rich.RESID[[j]][,"m"])))
  temp.biom.RESID <- cbind(temp.biom.RESID,as.vector(as.data.frame(coef.biom.RESID[[j]][,"m"])))
  temp.rich.CRYPTIC <- cbind(temp.rich.CRYPTIC,as.vector(as.data.frame(coef.rich.CRYPTIC[[j]][,"m"])))
  temp.biom.CRYPTIC <- cbind(temp.biom.CRYPTIC,as.vector(as.data.frame(coef.biom.CRYPTIC[[j]][,"m"])))
  
}
# TRANSIENT
 # richness
 rm(df.rich.TRANSIENT)
 df.rich.TRANSIENT <- temp.rich.TRANSIENT 
 rownames(df.rich.TRANSIENT) <- temp.rich.TRANSIENT[,1]; df.rich.TRANSIENT$temp.rich <- NULL
 saveRDS(df.rich.TRANSIENT,"Richness.power.pred.TRANSIENT")
 # biomass
 rm(df.biom.TRANSIENT)
 df.biom.TRANSIENT <- temp.biom.TRANSIENT 
 rownames(df.biom.TRANSIENT) <- temp.biom.TRANSIENT[,1]; df.biom.TRANSIENT$temp.biom <- NULL
 saveRDS(df.biom.TRANSIENT,"Biomass.power.pred.TRANSIENT")

 # RESIDENT
 # richness
 rm(df.rich.RESIDENT)
 df.rich.RESIDENT <- temp.rich.RESIDENT 
 rownames(df.rich.RESIDENT) <- temp.rich.RESIDENT[,1]; df.rich.RESIDENT$temp.rich <- NULL
 saveRDS(df.rich.RESIDENT,"Richness.power.pred.RESIDENT")
 # biomass
 rm(df.biom.RESIDENT)
 df.biom.RESIDENT <- temp.biom.RESIDENT 
 rownames(df.biom.RESIDENT) <- temp.biom.RESIDENT[,1]; df.biom.RESIDENT$temp.biom <- NULL
 saveRDS(df.biom.RESIDENT,"Biomass.power.pred.RESIDENT")
 
 # PARENTAL
 # richness
 rm(df.rich.PARENTAL)
 df.rich.PARENTAL <- temp.rich.PARENTAL 
 rownames(df.rich.PARENTAL) <- temp.rich.PARENTAL[,1]; df.rich.PARENTAL$temp.rich <- NULL
 saveRDS(df.rich.PARENTAL,"Richness.power.pred.PARENTAL")
 # biomass
 rm(df.biom.PARENTAL)
 df.biom.PARENTAL <- temp.biom.PARENTAL 
 rownames(df.biom.PARENTAL) <- temp.biom.PARENTAL[,1]; df.biom.PARENTAL$temp.biom <- NULL
 saveRDS(df.biom.PARENTAL,"Biomass.power.pred.PARENTAL")
 
 # CRYPTIC
 # richness
 rm(df.rich.CRYPTIC)
 df.rich.CRYPTIC <- temp.rich.CRYPTIC 
 rownames(df.rich.CRYPTIC) <- temp.rich.CRYPTIC[,1]; df.rich.CRYPTIC$temp.rich <- NULL
 saveRDS(df.rich.CRYPTIC,"Richness.power.pred.CRYPTIC")
 # biomass
 rm(df.biom.CRYPTIC)
 df.biom.CRYPTIC <- temp.biom.CRYPTIC 
 rownames(df.biom.CRYPTIC) <- temp.biom.CRYPTIC[,1]; df.biom.CRYPTIC$temp.biom <- NULL
 saveRDS(df.biom.CRYPTIC,"Biomass.power.pred.CRYPTIC")
 
 
 
 
 

  
 
 
 
 
 

