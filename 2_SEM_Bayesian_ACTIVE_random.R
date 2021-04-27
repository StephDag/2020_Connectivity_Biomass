# SEM analysis - Global connectivity 
# author: Steph D'agata
# date: April 2020
# updated: March 2021
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

#Sys.setenv(MAKEFLAGS = "-j10") # 10 cores used
#install.packages(c("Rcpp", "RcppEigen", "RcppParallel", "StanHeaders"), type = "source")
#install.packages("rstan", type = "source")

# packages
require(dplyr)
require(here)
require(forcats)
require(brms)
require(rstan)
require(Rserve)
require(rJava)

# load data
#rm(all.data)
#all.data <- readRDS(here::here("_data","DataBiomassConnectivityBR_March2021.rds"))
#all.data<-read.csv(here("_data","Connectivity_Biomass_SEMGLMMDATA_March2021.csv"),h=T, stringsAsFactors = F,dec=".")

# check 
#head(all.data)
#summary(all.data)
#dim(all.data)
#summary(all.data)
#str(all.data)
#names(all.data)
#apply(all.data,2,class)
#dim(all.data)

# built SEM model
# scaled 
  # gravity of neightbour is highly positively correlated with connectivity variables

### 26/03/2021
  # SPECIES RICHNESS
  # rm "-way interactions for species richness in the full and alternative models (non comparable if not)
  # rm temperature (the 50% interval crosses the 0 line for all functional group)

  # BIOMASS
  # rm inflow from neighbour in the simplified model
  # rm productivity (the 50% interval crosses the 0 line for all functional group)

#### SPECIES RICHNESS = Environment + Connectivity
rm(MOD_1_S_null) # ENV + CON
MOD_1_S_null <- bf(Richness ~ 1 + (1 |region))

# models with connectivity
rm(MOD_1_S) # ENV + CON
MOD_1_S <- bf(Richness ~ temp + prod.annual + 
                log_btwdegree + log_CorridorIn + 
                log_InflowMPA + log_InflowNei + Netflow + 
                log_grav_total + Class + (1 + temp + prod.annual + 
                                            log_btwdegree + log_CorridorIn + 
                                            log_InflowMPA + log_InflowNei + Netflow + 
                                            log_grav_total  |region))

rm(MOD_2_S) # ENV + CON
MOD_2_S <- bf(Richness ~ temp + prod.annual + 
                log_btwdegree + log_Outflow + 
                log_Inflow + log_Indegree_Neigh + 
                log_grav_total + Class + (1 + temp + prod.annual + 
                                            log_btwdegree + log_Outflow + 
                                            log_Inflow + log_Indegree_Neigh + 
                                            log_grav_total |region))

rm(MOD_3_S) # ENV + CON
MOD_3_S <- bf(Richness ~ temp + prod.annual + log_Outdegree + 
                log_btwdegree +
                log_InflowNei + log_Indegree_MPA + 
                log_grav_total + Class + (1 + temp + prod.annual + log_Outdegree + 
                                            log_btwdegree +
                                            log_InflowNei + log_Indegree_MPA + 
                                            log_grav_total |region))

rm(MOD_4_S) # ENV + CON
MOD_4_S <- bf(Richness ~ temp + prod.annual + log_SelfR + 
                log_btwdegree + log_Indegree + Netflow + 
                log_grav_total + Class + 
                (1 + temp + prod.annual + log_SelfR + 
                   log_btwdegree + log_Indegree + Netflow + 
                   log_grav_total |region))

#### SPECIES RICHNESS = Environment
S_mod_nocon <- bf(Richness ~ temp +  prod.annual + Class + log_grav_total + (1 +temp +  prod.annual + Class + log_grav_total |region))
                  
#### BIOMASS = Environment + Connectivity
  # null model
rm(MOD_B.null) # ENV + CON
MOD_B.null <- bf(log_biomassarea ~ 1 + (1 |region))

# full model
rm(MOD_1_B) # ENV + CON
MOD_1_B <- bf(log_biomassarea ~ Richness + 
                                  temp + prod.annual + 
                                  log_grav_total*Class*Netflow + 
                                  log_InflowMPA + log_InflowNei + Class + log_SelfR +
                (1 + log_InflowMPA + log_InflowNei + Richness + 
                   temp + prod.annual + 
                   log_grav_total +log_SelfR  |region))

rm(MOD_2_B) # ENV + CON
MOD_2_B <- bf(log_biomassarea ~ Richness + 
                temp + prod.annual + 
                log_grav_total*Class*Netflow + 
                log_SelfR + 
                log_Indegree + 
                (1 + Richness + 
                   temp + prod.annual + 
                   log_grav_total + Netflow + 
                   log_SelfR + 
                   log_Indegree |region))

rm(MOD_3_B) # ENV + CON
MOD_3_B <- bf(log_biomassarea ~ Richness + 
                temp + prod.annual + 
                log_grav_total*Class*Netflow + 
                log_SelfR + 
                log_Indegree_MPA  + 
                (1 + Richness + 
                   temp + prod.annual + 
                   log_grav_total+Netflow + 
                   log_SelfR + 
                   log_Indegree_MPA|region))

rm(MOD_4_B) # ENV + CON
MOD_4_B <- bf(log_biomassarea ~ Richness + 
                temp + prod.annual + 
                log_grav_total*Class*Netflow + 
                log_SelfR +  
                (1 + Richness + 
                   temp + prod.annual + 
                   log_grav_total+Netflow + 
                   log_SelfR |region))

rm(MOD_5_B) # ENV + CON
MOD_5_B <- bf(log_biomassarea ~ Richness + 
                temp + prod.annual + 
                log_grav_total*Class*Netflow + 
                (1 +Richness + 
                   temp + prod.annual + 
                   log_grav_total+Netflow|region))
#### Biomass = Environment
B_mod_nocon <- bf(log_biomassarea ~ temp +  prod.annual + Class*log_grav_total + 
                    (1 + temp +  prod.annual + log_grav_total  |region))

### run species models
  # null model
MOD_S_1_run_TRANSIENT_null %>% rm()  # Full - connectivity through both S and B
MOD_S_1_run_TRANSIENT_null <-brm(MOD_1_S_null , data=TRANSIENT.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2)
saveRDS(MOD_S_1_run_TRANSIENT,"ACTIVE_models/Full/MOD_S_1_run_TRANSIENT_null.Rds")

# full model
# null model
#MOD_S_full_run_TRANSIENT %>% rm()  # Full - connectivity through both S and B
#MOD_S_full_run_TRANSIENT <-brm(MOD_full_B , data=TRANSIENT.std,cores=4,chains = 4,
#                                 iter = 5000, warmup = 1000,thin = 2)
#saveRDS(MOD_S_full_run_TRANSIENT,"ACTIVE_models/Full/MOD_S_full_run_TRANSIENT.Rds")

  # model 1
MOD_S_1_run_TRANSIENT %>% rm()  # Full - connectivity through both S and B
MOD_S_1_run_TRANSIENT <-brm(MOD_1_S , data=TRANSIENT.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_S_1_run_TRANSIENT,"ACTIVE_models/Full/MOD_S_1_run_TRANSIENT.Rds")
  # model 2
MOD_S_2_run_TRANSIENT %>% rm()  # Full - connectivity through both S and B
MOD_S_2_run_TRANSIENT <-brm(MOD_2_S , data=TRANSIENT.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_S_2_run_TRANSIENT,"ACTIVE_models/Full/MOD_S_2_run_TRANSIENT.Rds")

# model 3
MOD_S_3_run_TRANSIENT %>% rm()  # Full - connectivity through both S and B
MOD_S_3_run_TRANSIENT <-brm(MOD_3_S , data=TRANSIENT.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_S_3_run_TRANSIENT,"ACTIVE_models/Full/MOD_S_3_run_TRANSIENT.Rds")

# model 4
MOD_S_4_run_TRANSIENT %>% rm()  # Full - connectivity through both S and B
MOD_S_4_run_TRANSIENT <-brm(MOD_4_S , data=TRANSIENT.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_S_4_run_TRANSIENT,"ACTIVE_models/Full/MOD_S_4_run_TRANSIENT.Rds")

# model S env
MOD_S_env %>% rm()  # Full - connectivity through both S and B
MOD_S_env <-brm(S_mod_nocon , data=TRANSIENT.std,cores=4,chains = 4,
                            iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                            prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_S_env,"ACTIVE_models/Full/MOD_S_env.Rds")

#### test all S models
r2_bayes(MOD_S_1_run_TRANSIENT_null)
r2_bayes(MOD_S_1_run_TRANSIENT)
r2_bayes(MOD_S_2_run_TRANSIENT)
r2_bayes(MOD_S_3_run_TRANSIENT)
r2_bayes(MOD_S_4_run_TRANSIENT)
r2_bayes(MOD_S_env)

sort(model_weights(MOD_S_1_run_TRANSIENT_null,MOD_S_1_run_TRANSIENT,
              MOD_S_2_run_TRANSIENT,MOD_S_3_run_TRANSIENT,
              MOD_S_4_run_TRANSIENT,MOD_S_env,weights="loo"),decreasing = T)

summary(MOD_S_2_run_TRANSIENT)

plot(MOD_S_4_run_TRANSIENT)

loo(MOD_S_1_run_TRANSIENT_null,MOD_S_1_run_TRANSIENT,MOD_S_2_run_TRANSIENT,
    MOD_S_3_run_TRANSIENT,MOD_S_4_run_TRANSIENT, MOD_S_env)

### 

fit2 <- update(MOD_BIOM_1_run_TRANSIENT, formula. = . ~  Richness + 
                 temp + prod.annual + 
                 log_grav_total*Class*Netflow + 
                 log_InflowMPA + (1 |region),newdata=TRANSIENT.std)

loo(MOD_BIOM_1_run_TRANSIENT,fit2)
r2_bayes(MOD_BIOM_1_run_TRANSIENT)
r2_bayes(fit2)
summary(MOD_BIOM_1_run_TRANSIENT)


#fit3 <- update(MOD_S_1_run_TRANSIENT, formula. = . ~ temp + prod.annual + 
#                 log_btwdegree + log_CorridorIn + log_grav_total + (1 |region),newdata=TRANSIENT.std)
#mcmc_areas(as.matrix(fit3),prob_outer = .99)
#loo(MOD_S_1_run_TRANSIENT_null,MOD_S_1_run_TRANSIENT,fit2,fit3,fit4,fit5)
#plot(MOD_S_1_run_TRANSIENT)
#pp_check(MOD_S_1_run_TRANSIENT, resp="log_biomassarea")
#pp_check(fit2, resp="Richness", nsamples = 1000)

library(performance)
r2_bayes(fit5.env)
r2_bayes(MOD_S_1_run_TRANSIENT_null)

### run biomass models

## transient
MOD_BIOM_1_run_TRANSIENT.null %>% rm()  # Full - connectivity through both S and B
MOD_BIOM_1_run_TRANSIENT.null <-brm(MOD_B.null , data=TRANSIENT.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30))
saveRDS(MOD_BIOM_1_run_TRANSIENT,"ACTIVE_models/Full/MOD_BIOM_1_run_TRANSIENT.Rds")



MOD_BIOM_1_run_TRANSIENT %>% rm()  # Full - connectivity through both S and B
MOD_BIOM_1_run_TRANSIENT <-brm(MOD_1_B , data=TRANSIENT.std,cores=4,chains = 4,
                                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_BIOM_1_run_TRANSIENT,"ACTIVE_models/Full/MOD_BIOM_1_run_TRANSIENT.Rds")


MOD_BIOM_2_run_TRANSIENT %>% rm()  # Full - connectivity through both S and B
MOD_BIOM_2_run_TRANSIENT <-brm(MOD_2_B , data=TRANSIENT.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_BIOM_2_run_TRANSIENT,"ACTIVE_models/Full/MOD_BIOM_2_run_TRANSIENT.Rds")


MOD_BIOM_3_run_TRANSIENT %>% rm()  # Full - connectivity through both S and B
MOD_BIOM_3_run_TRANSIENT <-brm(MOD_3_B , data=TRANSIENT.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_BIOM_3_run_TRANSIENT,"ACTIVE_models/Full/MOD_BIOM_3_run_TRANSIENT.Rds")

MOD_BIOM_4_run_TRANSIENT %>% rm()  # Full - connectivity through both S and B
MOD_BIOM_4_run_TRANSIENT <-brm(MOD_4_B , data=TRANSIENT.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_BIOM_4_run_TRANSIENT,"ACTIVE_models/Full/MOD_BIOM_4_run_TRANSIENT.Rds")

MOD_BIOM_5_run_TRANSIENT %>% rm()  # Full - connectivity through both S and B
MOD_BIOM_5_run_TRANSIENT <-brm(MOD_5_B , data=TRANSIENT.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_BIOM_5_run_TRANSIENT,"ACTIVE_models/Full/MOD_BIOM_5_run_TRANSIENT.Rds")

# model B env
MOD_B_env %>% rm()  # Full - connectivity through both S and B
MOD_B_env <-brm(B_mod_nocon , data=TRANSIENT.std,cores=4,chains = 4,
                iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_B_env,"ACTIVE_models/Full/MOD_B_env.Rds")

#### test all S models
r2_bayes(MOD_BIOM_1_run_TRANSIENT.null)
r2_bayes(MOD_BIOM_1_run_TRANSIENT)
r2_bayes(MOD_BIOM_2_run_TRANSIENT)
r2_bayes(MOD_BIOM_3_run_TRANSIENT)
r2_bayes(MOD_BIOM_4_run_TRANSIENT)
r2_bayes(MOD_BIOM_5_run_TRANSIENT)
r2_bayes(MOD_B_env)
summary(MOD_B_env)

bayes_R2(MOD_BIOM_1_run_TRANSIENT)

loo(MOD_BIOM_1_run_TRANSIENT.null,MOD_BIOM_1_run_TRANSIENT,MOD_BIOM_2_run_TRANSIENT,
    MOD_BIOM_3_run_TRANSIENT,MOD_BIOM_4_run_TRANSIENT,MOD_BIOM_5_run_TRANSIENT,MOD_B_env)

sort(model_weights(MOD_BIOM_1_run_TRANSIENT.null,MOD_BIOM_1_run_TRANSIENT,MOD_BIOM_2_run_TRANSIENT,
              MOD_BIOM_3_run_TRANSIENT,MOD_BIOM_4_run_TRANSIENT,MOD_BIOM_5_run_TRANSIENT,MOD_B_env,weights="loo"),decreasing=T)

plot(conditional_effects(MOD_BIOM_2_run_TRANSIENT))

mcmc_plot(MOD_BIOM_1_run_TRANSIENT, pars = "^b_")
mcmc_plot(MOD_B_env, pars = "^b_")

mcmc_plot(MOD_S_2_run_TRANSIENT, pars = "^b_")
mcmc_plot(MOD_S_env, pars = "^b_")

posterior <- as.array(MOD_BIOM_1_run_TRANSIENT)
dim(posterior)

if (requireNamespace("hexbin", quietly = TRUE)) {
  mcmc_hex(posterior, pars = c("b_log_InflowMPA", "b_log_InflowNei"))
}

mcmc_pairs(posterior, pars = c("b_Netflow","b_log_InflowMPA", "b_log_InflowNei",
                               "b_log_grav_total","b_log_SelfR"), diag_fun = "dens", off_diag_fun = "hex")
color_scheme_set("pink")
mcmc_pairs(posterior, pars = c("Netflow", "log_InflowMPA"),
           off_diag_args = list(size = 1.5))
#
color_scheme_set("pink")
(p3 <- mcmc_scatter(posterior, pars = c("b_log_InflowMPA", "b_log_InflowNei"), alpha = 0.25, size = 3))
p3 + geom_smooth(method = "lm", se = FALSE, color = "gray20",
                 size = .75, linetype = 2)

(p3 <- mcmc_scatter(posterior, pars = c("b_Netflow", "b_log_InflowNei"), alpha = 0.25, size = 3))
p3 + geom_smooth(method = "lm", se = FALSE, color = "gray20",
                 size = .75, linetype = 2)

(p3 <- mcmc_scatter(posterior, pars = c("b_Netflow", "b_log_InflowNei"), alpha = 0.25, size = 3))
p3 + geom_smooth(method = "lm", se = FALSE, color = "gray20",
                 size = .75, linetype = 2)

(p3 <- mcmc_scatter(posterior, pars = c("b_log_InflowNei", "b_log_InflowMPA"), alpha = 0.25, size = 3))
p3 + geom_smooth(method = "lm", se = FALSE, color = "gray20",
                 size = .75, linetype = 2)
(p3 <- mcmc_scatter(posterior, pars = c("b_Netflow", "b_log_SelfR"), alpha = 0.25, size = 3))
p3 + geom_smooth(method = "lm", se = FALSE, color = "gray20",
                 size = .75, linetype = 2)



#> `geom_smooth()` using formula 'y ~ x'

summary(MOD_BIOM_1_run_TRANSIENT)
summary(MOD_BIOM_2_run_TRANSIENT)
summary(MOD_BIOM_3_run_TRANSIENT)
summary(MOD_BIOM_4_run_TRANSIENT)
summary(MOD_B_env)
plot(MOD_BIOM_1_run_TRANSIENT)
rm(fit2)
fit2 <- update(MOD_BIOM_1_run_TRANSIENT, formula. = . ~ Richness + 
                 temp + prod.annual + 
                 log_grav_total*Class*Netflow +  (1 |region),newdata=TRANSIENT.std)
rm(fit3)
fit3 <- update(MOD_BIOM_1_run_TRANSIENT, formula. = . ~ Richness + 
                 prod.annual + 
                 log_grav_total*Class*Netflow +  (1 |region),newdata=TRANSIENT.std)

fit.env <- update(MOD_BIOM_1_run_TRANSIENT, formula. = . ~  prod.annual + 
                     temp  +   log_grav_total*Class + (1 |region),newdata=TRANSIENT.std)


mcmc_areas(as.matrix(fit3),prob_outer = .99)
summary(MOD_BIOM_1_run_TRANSIENT)
summary(fit3)
summary(fit4)
summary(fit5)
loo(MOD_BIOM_1_run_TRANSIENT.null,MOD_BIOM_1_run_TRANSIENT,fit2)
#,fit3,fit4,fit5)

plot(MOD_S_1_run_TRANSIENT)

plot(fit4)
plot(fit5)
plot(fit5.env)

brms::pp_check(MOD_BIOM_1_run_TRANSIENT, resp="log_biomassarea", nsamples = 100)
pp_check(fit2, resp="Richness", nsamples = 1000)
bayes_R2(MOD_BIOM_1_run_TRANSIENT)
bayes_R2(fit2)
bayes_R2(fit3)
bayes_R2(fit4)
bayes_R2(fit5)
bayes_R2(fit.env)

library(performance)
r2_bayes(MOD_BIOM_1_run_TRANSIENT.null)
r2_bayes(MOD_BIOM_1_run_TRANSIENT)
r2_bayes(fit2)
r2_bayes(fit3)
r2_bayes(fit.env)

bayes_R2(MOD_S_1_run_TRANSIENT_null)


## cryptic


## transient
MOD_BIOM_1_run_CRYPTIC.null %>% rm()  # Full - connectivity through both S and B
MOD_BIOM_1_run_CRYPTIC.null <-brm(MOD_B.null , data=CRYPTIC.std,cores=4,chains = 4,
                                    iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30))
saveRDS(MOD_BIOM_1_run_CRYPTIC,"ACTIVE_models/Full/MOD_BIOM_1_run_CRYPTIC.Rds")



MOD_BIOM_1_run_CRYPTIC %>% rm()  # Full - connectivity through both S and B
MOD_BIOM_1_run_CRYPTIC <-brm(MOD_1_B , data=CRYPTIC.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_BIOM_1_run_CRYPTIC,"ACTIVE_models/Full/MOD_BIOM_1_run_CRYPTIC.Rds")


MOD_BIOM_2_run_CRYPTIC %>% rm()  # Full - connectivity through both S and B
MOD_BIOM_2_run_CRYPTIC <-brm(MOD_2_B , data=CRYPTIC.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_BIOM_2_run_CRYPTIC,"ACTIVE_models/Full/MOD_BIOM_2_run_CRYPTIC.Rds")


MOD_BIOM_3_run_CRYPTIC %>% rm()  # Full - connectivity through both S and B
MOD_BIOM_3_run_CRYPTIC <-brm(MOD_3_B , data=CRYPTIC.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_BIOM_3_run_TRANSIENT,"ACTIVE_models/Full/MOD_BIOM_3_run_CRYPTIC.Rds")

MOD_BIOM_4_run_TRANSIENT %>% rm()  # Full - connectivity through both S and B
MOD_BIOM_4_run_TRANSIENT <-brm(MOD_4_B , data=TRANSIENT.std,cores=4,chains = 4,
                               iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                               prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(MOD_BIOM_4_run_TRANSIENT,"ACTIVE_models/Full/MOD_BIOM_4_run_CRYPTIC.Rds")

#### test all S models
r2_bayes(MOD_BIOM_1_run_TRANSIENT.null)
r2_bayes(MOD_BIOM_1_run_TRANSIENT)
r2_bayes(MOD_BIOM_2_run_TRANSIENT)
r2_bayes(MOD_BIOM_3_run_TRANSIENT)
r2_bayes(MOD_BIOM_4_run_TRANSIENT)

loo(MOD_BIOM_1_run_TRANSIENT.null,MOD_BIOM_1_run_TRANSIENT,MOD_BIOM_2_run_TRANSIENT,
    MOD_BIOM_3_run_TRANSIENT,MOD_BIOM_4_run_TRANSIENT)




#### END OF MODELS

#
library("bayesplot")
library("ggplot2")
library("rstanarm")   
require(ggpubr)

# classify explanatory variables based on categories


# TRANSIENT TOT
a <- mcmc_intervals(all_fit_brms.tot.TRANSIENT.intr)

# no classification
rich.var <- as.data.frame(a$data)[grep("b_Richness",as.data.frame(a$data)[,"parameter"]),"parameter"]
biom.var <- as.data.frame(a$data)[grep("b_logbiomassarea",as.data.frame(a$data)[,"parameter"]),"parameter"]

Richness_inflow_TRANSIENT_tot <- mcmc_intervals(all_fit_brms.tot.TRANSIENT.intr, pars = as.character(rich.var))
Biomass_inflow_TRANSIENT_tot <- mcmc_intervals(all_fit_brms.tot.TRANSIENT.intr,pars = as.character(biom.var))

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


