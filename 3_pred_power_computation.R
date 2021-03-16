# SEM analysis - Global connectivity 
# author: Steph D'agata
# date: September 2020
# outputs: Power of the coefficients for each variables

rm(list=ls())

# load packages
library(dplyr)
library(here)
library("bayesplot")
library("ggplot2")
library("rstanarm")   
require(ggpubr)
require(bayestestR)

# load data of power
Richness.power.pred.CRYPTIC <- readRDS("~/Documents/GitHub/2020_Connectivity_Biomass_V2/Richness.power.pred.CRYPTIC.rds")
Richness.power.pred.PARENTAL <- readRDS("~/Documents/GitHub/2020_Connectivity_Biomass_V2/Richness.power.pred.PARENTAL.rds")
Richness.power.pred.RESIDENT <- readRDS("~/Documents/GitHub/2020_Connectivity_Biomass_V2/Richness.power.pred.RESIDENT.rds")
Richness.power.pred.TRANSIENT <- readRDS("~/Documents/GitHub/2020_Connectivity_Biomass_V2/Richness.power.pred.TRANSIENT.rds")

Biomass.power.pred.CRYPTIC <- readRDS("~/Documents/GitHub/2020_Connectivity_Biomass_V2/Biomass.power.pred.CRYPTIC.rds")
Biomass.power.pred.PARENTAL <- readRDS("~/Documents/GitHub/2020_Connectivity_Biomass_V2/Biomass.power.pred.PARENTAL.rds")
Biomass.power.pred.RESIDENT <- readRDS("~/Documents/GitHub/2020_Connectivity_Biomass_V2/Biomass.power.pred.RESIDENT.rds")
Biomass.power.pred.TRANSIENT <- readRDS("~/Documents/GitHub/2020_Connectivity_Biomass_V2/Biomass.power.pred.TRANSIENT.rds")

# load models
all_fit_brms.tot.TRANSIENT.intr.extr <- readRDS("~/Documents/GitHub/2020_Connectivity_Biomass_V2/all_fit_brms.tot.TRANSIENT.intr.extr.Rds")
all_fit_brms.tot.PARENTAL.intr.extr <- readRDS("~/Documents/GitHub/2020_Connectivity_Biomass_V2/all_fit_brms.tot.PARENTAL.intr.extr.Rds")
all_fit_brms.tot.CRYPTIC.intr.extr <- readRDS("~/Documents/GitHub/2020_Connectivity_Biomass_V2/all_fit_brms.tot.CRYPTIC.intr.extr.Rds")
all_fit_brms.tot.RESIDENT.intr.extr <- readRDS("~/Documents/GitHub/2020_Connectivity_Biomass_V2/all_fit_brms.tot.RESID.intr.extr.Rds")

# TRANSIENT TOT
a.transient <- mcmc_intervals(all_fit_brms.tot.TRANSIENT.intr.extr)
hdi_transient <- hdi(all_fit_brms.tot.TRANSIENT.intr.extr,ci=0.95)

# no classification
rich.var <- as.data.frame(a.transient$data)[grep("b_Richness",as.data.frame(a.transient$data)[,"parameter"]),"parameter"]
biom.var <- as.data.frame(a.transient$data)[grep("b_logbiomassarea",as.data.frame(a.transient$data)[,"parameter"]),"parameter"]

Richness_inflow_TRANSIENT_tot <- mcmc_intervals_data(all_fit_brms.tot.TRANSIENT.intr.extr, pars = as.character(rich.var))
Biomass_inflow_TRANSIENT_tot <- mcmc_intervals_data(all_fit_brms.tot.TRANSIENT.intr.extr,pars = as.character(biom.var))

# PARENTAL
a.parental <- mcmc_intervals(all_fit_brms.tot.PARENTAL.intr.extr)
hdi_parental <- hdi(all_fit_brms.tot.PARENTAL.intr.extr,ci=0.95)

# no classification
rich.var <- as.data.frame(a.parental$data)[grep("b_Richness",as.data.frame(a.parental$data)[,"parameter"]),"parameter"]
biom.var <- as.data.frame(a.parental$data)[grep("b_logbiomassarea",as.data.frame(a.parental$data)[,"parameter"]),"parameter"]

Richness_inflow_PARENTAL_tot <- mcmc_intervals_data(all_fit_brms.tot.PARENTAL.intr.extr, pars = as.character(rich.var))
Biomass_inflow_PARENTAL_tot <- mcmc_intervals_data(all_fit_brms.tot.PARENTAL.intr.extr,pars = as.character(biom.var))

# CRYPTIC
a.cryptic <- mcmc_intervals(all_fit_brms.tot.CRYPTIC.intr.extr)
hdi_cryptic <- hdi(all_fit_brms.tot.CRYPTIC.intr.extr,ci=0.95)

# no classification
rich.var <- as.data.frame(a.cryptic$data)[grep("b_Richness",as.data.frame(a.cryptic$data)[,"parameter"]),"parameter"]
biom.var <- as.data.frame(a.cryptic$data)[grep("b_logbiomassarea",as.data.frame(a.cryptic$data)[,"parameter"]),"parameter"]

Richness_inflow_CRYPTIC_tot <- mcmc_intervals_data(all_fit_brms.tot.CRYPTIC.intr.extr, pars = as.character(rich.var))
Biomass_inflow_CRYPTIC_tot <- mcmc_intervals_data(all_fit_brms.tot.CRYPTIC.intr.extr,pars = as.character(biom.var))


# RESIDENT
a.resident <- mcmc_intervals(all_fit_brms.tot.RESIDENT.intr.extr)
hdi_resident <- hdi(all_fit_brms.tot.RESIDENT.intr.extr,ci=0.95)

# no classification
rich.var <- as.data.frame(a.resident$data)[grep("b_Richness",as.data.frame(a.resident$data)[,"parameter"]),"parameter"]
biom.var <- as.data.frame(a.resident$data)[grep("b_logbiomassarea",as.data.frame(a.resident$data)[,"parameter"]),"parameter"]

Richness_inflow_RESIDENT_tot <- mcmc_intervals_data(all_fit_brms.tot.RESIDENT.intr.extr, pars = as.character(rich.var))
Biomass_inflow_RESIDENT_tot <- mcmc_intervals_data(all_fit_brms.tot.RESIDENT.intr.extr,pars = as.character(biom.var))

##### compute power analysis for richness
rm(power.rich.transient)
power.rich.transient <- matrix()
power.rich.parental <- matrix()
power.rich.cryptic <- matrix()
power.rich.resident <- matrix()

for (i in 1:20){
 power.rich.transient[i] <-  length(which(between(as.numeric(Richness.power.pred.TRANSIENT[i,]),as.data.frame(Richness_inflow_TRANSIENT_tot[i,"ll"]),as.data.frame(Richness_inflow_TRANSIENT_tot[i,"hh"]))))
 power.rich.parental[i] <-  length(which(between(as.numeric(Richness.power.pred.PARENTAL[i,]),as.data.frame(Richness_inflow_PARENTAL_tot[i,"ll"]),as.data.frame(Richness_inflow_PARENTAL_tot[i,"hh"]))))
 power.rich.cryptic[i] <-  length(which(between(as.numeric(Richness.power.pred.CRYPTIC[i,]),as.data.frame(Richness_inflow_CRYPTIC_tot[i,"ll"]),as.data.frame(Richness_inflow_CRYPTIC_tot[i,"hh"]))))
 power.rich.resident[i] <-  length(which(between(as.numeric(Richness.power.pred.RESIDENT[i,]),as.data.frame(Richness_inflow_RESIDENT_tot[i,"ll"]),as.data.frame(Richness_inflow_RESIDENT_tot[i,"hh"]))))
}

power.rich.transient <- round(power.rich.transient/106,2)
power.rich.parental <- round(power.rich.parental/106,2)
power.rich.cryptic <- round(power.rich.cryptic/106,2)
power.rich.resident <- round(power.rich.resident/106,2)

power.rich <- cbind(as.character(Richness_inflow_TRANSIENT_tot$parameter),power.rich.transient,power.rich.parental,power.rich.cryptic,power.rich.resident)
power.rich <- power.rich[-1,]
write.csv(power.rich ,"Pred_powr_richness.csv")
##### compute power analysis for biomass
rm(power.biom.transient)
power.biom.transient <- matrix()
power.biom.parental <- matrix()
power.biom.cryptic <- matrix()
power.biom.resident <- matrix()

for (i in 1:21){
  power.biom.transient[i] <-  length(which(between(as.numeric(Biomass.power.pred.TRANSIENT[i,]),as.data.frame(Biomass_inflow_TRANSIENT_tot[i,"ll"]),as.data.frame(Biomass_inflow_TRANSIENT_tot[i,"hh"]))))
  power.biom.parental[i] <-  length(which(between(as.numeric(Biomass.power.pred.PARENTAL[i,]),as.data.frame(Biomass_inflow_PARENTAL_tot[i,"ll"]),as.data.frame(Biomass_inflow_PARENTAL_tot[i,"hh"]))))
  power.biom.cryptic[i] <-  length(which(between(as.numeric(Biomass.power.pred.CRYPTIC[i,]),as.data.frame(Biomass_inflow_CRYPTIC_tot[i,"ll"]),as.data.frame(Biomass_inflow_CRYPTIC_tot[i,"hh"]))))
  power.biom.resident[i] <-  length(which(between(as.numeric(Biomass.power.pred.RESIDENT[i,]),as.data.frame(Biomass_inflow_RESIDENT_tot[i,"ll"]),as.data.frame(Biomass_inflow_RESIDENT_tot[i,"hh"]))))
}
power.biom.transient <- round(power.biom.transient/106,2)
power.biom.parental <- round(power.biom.parental/106,2)
power.biom.cryptic <- round(power.biom.cryptic/106,2)
power.biom.resident <- round(power.biom.resident/106,2)

power.biom <- cbind(as.character(Biomass_inflow_TRANSIENT_tot$parameter),power.biom.transient,power.biom.parental,power.biom.cryptic,power.biom.resident)
power.biom <- power.biom[-c(1,2),]
write.csv(power.biom,"Pred_powr_biomass.csv")
tail(power.biom)
