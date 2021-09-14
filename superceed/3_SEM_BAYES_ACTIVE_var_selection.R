
library(rstanarm)
options(mc.cores = parallel::detectCores())
library(loo)
library(tidyverse)
library(GGally)
library(bayesplot)
theme_set(bayesplot::theme_default())
library(projpred)


Full <- all_fit_brms.tot.TRANSIENT.intr.extr

summary(Full)


mcmc_areas(as.matrix(Full),prob_outer = .99)

# full
TRANSIENT.S.full %>% rm()  # Full - connectivity through both S and B
TRANSIENT.S.full <-brm(Richness ~ temp + prod.annual + 
                         log_btwdegree + log_CorridorIn + 
                         log_InflowMPA + log_InflowNei + log_Inflow + log_Indegree +
                         log_Indegree_Neigh + log_Indegree_MPA + 
                         Netflow +
                         log_grav_total + log_grav_neiBR + log_Indegree + 
                         Age_of_pro +  Class + (1 |region), data=TRANSIENT.std,cores=4,chains = 4,
                                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30)
)

TRANSIENT.B.full %>% rm()  # Full - connectivity through both S and B
TRANSIENT.B.full <-brm(log_biomassarea ~ Richness + temp + prod.annual + 
                         log_btwdegree + log_CorridorIn + 
                         log_InflowMPA + log_InflowNei + log_Inflow + log_Indegree +
                         log_Indegree_Neigh + log_Indegree_MPA +
                         log_grav_total+Class+Netflow + log_grav_neiBR + log_Indegree + 
                        (1 |region), data=TRANSIENT.std,cores=4,chains = 4,
                       iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30)
)




saveRDS(all_fit_brms.tot.TRANSIENT.intr.extr,"ACTIVE_models/Full/all_fit_brms.tot.TRANSIENT.intr.extr.Rds")

# full CRYPTIC
CRYPTIC.S.full %>% rm()  # Full - connectivity through both S and B
CRYPTIC.S.full <-brm(Richness ~ temp + prod.annual + 
                         log_btwdegree + log_CorridorIn + 
                         log_InflowMPA + log_InflowNei + log_Inflow + log_Indegree +
                         log_Indegree_Neigh + log_Indegree_MPA + 
                         Netflow +
                         log_grav_total + log_grav_neiBR + log_Indegree + 
                         Age_of_pro +  Class + (1 |region), data=CRYPTIC.std,cores=4,chains = 4,
                       iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30)
)

CRYPTIC.B.full %>% rm()  # Full - connectivity through both S and B
CRYPTIC.B.full <-brm(log_biomassarea ~ Richness + temp + prod.annual + 
                         log_btwdegree + log_CorridorIn + 
                         log_InflowMPA + log_InflowNei + log_Inflow + log_Indegree +
                         log_Indegree_Neigh + log_Indegree_MPA + 
                         log_grav_total+Class+Netflow + log_grav_neiBR + log_Indegree + 
                         (1 |region), data=CRYPTIC.std,cores=4,chains = 4,
                       iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30)
)



# fit null model
TRANSIENT.S.null %>% rm()  # Full - connectivity through both S and B
TRANSIENT.S.null <-brm(Richness ~ 1, data=TRANSIENT.std,cores=4,chains = 4,
                       iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30)
)

TRANSIENT.B.null %>% rm()  # Full - connectivity through both S and B
TRANSIENT.B.null <-brm(log_biomassarea ~ 1, data=TRANSIENT.std,cores=4,chains = 4,
                       iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30)
)
#
#TRANSIENT full
TRANSIENT.null %>% rm()  # Full - connectivity through both S and B
TRANSIENT.null <-brm(S.null + B.null + set_rescor(FALSE), data=TRANSIENT.std,cores=4,chains = 4,
                                           iter = 5000, warmup = 1000,thin = 2, refresh = 0, control = list(adapt_delta = 0.99999,max_treedepth = 30),
                                           prior = c(prior(normal(0, 100),class = "Intercept"), prior(normal(0, 100), class = "b")))
saveRDS(TRANSIENT.null,"ACTIVE_models/Full/TRANSIENT.null.Rds")

mcmc_pairs(as.matrix(Full))


# posterior checks



# check loo
(looFull.S <- loo(TRANSIENT.S.full))
(looFull.B <- loo(TRANSIENT.B.full))
(looNull.S <- loo(TRANSIENT.S.null))
(looNull.B <- loo(TRANSIENT.B.null))

loo_compare(looFull.S, looNull.S)
loo_compare(looFull.B, looNull.B)

# cross validation 
fitg_cv.S <- projpred::cv_varsel(MOD_S_1_run_TRANSIENT, method='forward', cv_method='LOO')
plot(fitg_cv.S )
(nv <- fitg_cv.S$suggested_size)
nv <- 10
projg <- project(fitg_cv.B, nv = nv, ns = 4000)
round(colMeans(as.matrix(projg)),1)
round(posterior_interval(as.matrix(projg)),1)
fitg_cv.S$solution_terms

mcmc_areas(as.matrix(projg), prob_outer=0.99)
           pars = c('(Intercept)',fitg_cv.S$solution_terms[1:nv]), )

plot(fitg_cv.S)
projpred::varsel_plot(fitg_cv.S, stats = c('elpd', 'rmse'))

varsel_plot(vs, stats=c('elpd', 'rmse'))




# Biomass transient
rm(fitg_cv.B)
fitg_cv.B <- projpred::cv_varsel(MOD_BIOM_1_run_TRANSIENT, method='forward', cv_method='LOO')
plot(fitg_cv.B)

fitg.B <- projpred::varsel(MOD_BIOM_1_run_TRANSIENT, method='forward', cv_method='LOO')
fitg.B$suggested_size

projg.B.bis <- projpred::project(fitg.B,fitg.B$suggested_size)
round(colMeans(as.matrix(projg.B.bis)),1)
round(posterior_interval(as.matrix(projg.B.bis)),1)

fitg_cv.B.cryp<- projpred::cv_varsel(CRYPTIC.B.full, method='forward', cv_method='LOO')
plot(fitg_cv.B.cryp)
(nv.B <- fitg_cv.B$suggested_size)

projg.B <- projpred::project(fitg_cv.B, nv.B)
round(colMeans(as.matrix(projg.B)),1)
round(posterior_interval(as.matrix(projg.B)),1)
fitg_cv.B$solution_terms

mcmc_areas(as.matrix(projg.B), prob_outer=0.99)
pars = c('(Intercept)',fitg_cv.S$solution_terms[1:nv]), )

# Biomass cryptic
fitg_cv.B <- projpred::cv_varsel(CRYPTIC.B.full, method='forward', cv_method='LOO')

(nv.B <- fitg_cv.B$suggested_size)


projg.B.C <- projpred::project(fitg_cv.B.cryp, nterms = nv.B, ndraws = 4000)



projg.B <- project(fitg_cv.B.cryp, nv = 10, ns = 4000)
round(colMeans(as.matrix(projg.B)),1)
round(posterior_interval(as.matrix(projg.B)),1)
fitg_cv.S$solution_terms

mcmc_areas(as.matrix(projg.B), prob_outer=0.99)
pars = c('(Intercept)',fitg_cv.S$solution_terms[1:nv]), )

plot(TRANSIENT.std$log_grav_total,TRANSIENT.std$log_biomassarea)
plot(TRANSIENT.std$Netflow,TRANSIENT.std$log_biomassarea)
plot(TRANSIENT.std$log_btwdegree,TRANSIENT.std$Netflow)

