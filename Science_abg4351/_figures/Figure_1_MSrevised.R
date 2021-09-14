# packages
require(dplyr)
require(here)
require(forcats)
require(brms)
require(rstan)
library("bayesplot")
library("ggplot2")
library("rstanarm")   
require(ggpubr)
### marginal plots
library(ggeffects)
library(here)

#Figure 1 - Part B and C
#BIOMASS models
#GLobal model
MOD_BIOM_1_run_ACTIVE.1 <- readRDS(here("Models output","GLOBAL","ACTIVE 1","Biomass","MOD_BIOM_1_run_ACTIVE.1_V2_May.Rds"))

#Resident 15
MOD_BIOM_1_run_RESID <- readRDS(here("Models output","ACTIVE_1","Resident","Biomass","MOD_BIOM_1_run_RESID_V2_May.Rds"))

#Cryptic 5
MOD_BIOM_1_run_CRYPTIC <- readRDS(here("Models output","ACTIVE_1","CRYPTIC","Biomass","MOD_BIOM_1_run_CRYPTIC_V2_May.Rds"))


#RICHNESS models
#GLobal model
MOD_S_1_run_ACTIVE.1 <- readRDS(here("Models output","GLOBAL","ACTIVE 1","Richness","MOD_S_1_run_ACTIVE.1_V2_May.Rds"))

#Transient 15
MOD_S_1_run_TRANSIENT <- readRDS(here("Models output","ACTIVE_1","Transient","Richness","MOD_S_1_run_TRANSIENT_V2_May.Rds"))

#Parental 5
MOD_S_1_run_PARENTAL <- readRDS(here("Models output","ACTIVE_1","PARENTAL","Richness","MOD_S_1_run_PARENTAL_V2_May.Rds"))

#
# Figure 1 MS - Global Active 1, Transient and Cryptic##----
# Full model # ENV + CON (Netflow, Inflow of upstream reefs, self-retention)
# Plot MOD_B_1
library(dotwhisker)
library(broom)
library(dplyr)

# TRANSIENT TOT
rm(a)
a <- mcmc_intervals(MOD_BIOM_1_run_ACTIVE.1)
# no classification
rm(biom.var)
biom.var <- a$data[1:10,]$parameter

Biomass_inflow_GLOBAL_simp_data <- mcmc_intervals_data(MOD_BIOM_1_run_ACTIVE.1,pars = as.character(biom.var))
Biomass_inflow_RESID_simp_data <- mcmc_intervals_data(MOD_BIOM_1_run_RESID,pars = as.character(biom.var))
Biomass_inflow_CRYPTO_simp_data <- mcmc_intervals_data(MOD_BIOM_1_run_CRYPTIC,pars = as.character(biom.var))

rm(a)


# Classification
INT <- "Intercept"
CON <- c("log_SelfR","Netflow","log_InflowNei")
HUMAN <- c("log_grav_total","_ClassClosed","_ClassRestricted")
ENV <- c("Richness","temp","prod.annual")

PARAM <- c("Intercept","Richness","Temperature","Productivity",
           "Tot.Gravity","No-Take","Restricted gears",
            "Netflow","Local Recruit.","Exogenous Inflow")

CAT <- c("Intercept","Human/Env.","Human/Env.","Human/Env.","Human/Env.","Human/Env.",
         "Human/Env.",
         "Connectivity",
         "Connectivity",
         "Connectivity")


DesiredOrder <- c("Intercept",
                  "Richness",
                  "Temperature",
                  "Productivity",
                  "Tot.Gravity",
                  "No-Take",
                  "Restricted gears",
                  "Netflow",
                  "Local Recruit.",
                  "Exogenous Inflow")

#GLOBAL
rm(m2_df_global)
m2_df_global <- Biomass_inflow_GLOBAL_simp_data %>%  mutate(FE = "GLOBAL") %>%  mutate(model = "Biomass")
m2_df_global <- cbind(PARAM,m2_df_global,CAT)
m2_df_global$PARAM <- factor(m2_df_global$PARAM, levels = rev(DesiredOrder))
m2_df_global$CAT <- factor(m2_df_global$CAT , levels = c("Intercept","Human/Env.","Connectivity"),ordered = TRUE)

#TRANSIENT
rm(m2_df_resident)
m2_df_resident <- Biomass_inflow_RESID_simp_data %>%  mutate(FE = "Resident") %>%  mutate(model = "Biomass")
m2_df_resident <- cbind(PARAM,m2_df_resident,CAT)
m2_df_resident$PARAM <- factor(m2_df_resident$PARAM, levels = rev(DesiredOrder))
m2_df_resident$CAT <- factor(m2_df_resident$CAT , levels = c("Intercept","Human/Env.","Connectivity"),ordered = TRUE)

#CRYPTIC
rm(m2_df_crypto)
m2_df_crypto <- Biomass_inflow_CRYPTO_simp_data %>%  mutate(FE = "Crypto") %>%  mutate(model = "Biomass")
m2_df_crypto <- cbind(PARAM,m2_df_crypto,CAT)
m2_df_crypto$PARAM <- factor(m2_df_crypto$PARAM, levels = rev(DesiredOrder))
m2_df_crypto$CAT <- factor(m2_df_crypto$CAT , levels = c("Intercept","Human/Env.","Connectivity"),ordered = TRUE)


# biomass data.frame
rm(four_FE_Biomass)
four_FE_Biomass <- rbind(m2_df_global,m2_df_resident,m2_df_crypto)
four_FE_Biomass$FE <- factor(four_FE_Biomass$FE, levels=c("Crypto","Resident","GLOBAL"))
four_FE_Biomass$CAT <- as.factor(four_FE_Biomass$CAT)
four_FE_Biomass$model <- as.factor(four_FE_Biomass$model)
summary(four_FE_Biomass)

library(viridis)
# plot biomass
rm(Biomass.FE)
Biomass.FE <- ggplot(four_FE_Biomass %>% filter(PARAM != "Intercept"), aes(colour = FE)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = FE),lwd = 1/2, position = position_dodge(width = 0.8)) +
  geom_point(aes(x = PARAM, y = m,colour = FE,size=6), lwd = 1.5, position = position_dodge(width = 0.8),shape = 21,fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = FE),lwd = 1, position = position_dodge(width = 0.8)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") +  scale_color_viridis(discrete=TRUE) +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") + theme(text = element_text(size = 10, family="Arial")) +
  ggtitle("Biomass") 
Biomass.FE  # The trick to these is position_dodge()


#Species richness PLOT FIG 1 Parental and Transient (highest model weight)
# 

rm(a)
a <- mcmc_intervals(MOD_S_1_run_ACTIVE.1)
# no classification
rm(rich.var)
rich.var <- a$data[1:10,]$parameter

Richness_inflow_GLOBAL_simp_data <- mcmc_intervals_data(MOD_S_1_run_ACTIVE.1,pars = as.character(rich.var))
Richness_inflow_TRANSIENT_simp_data <- mcmc_intervals_data(MOD_S_1_run_TRANSIENT,pars = as.character(rich.var))
Richness_inflow_CRYPTO_simp_data <- mcmc_intervals_data(MOD_S_1_run_PARENTAL,pars = as.character(rich.var))

rm(a)


# Classification
INT <- "Intercept"
CON <- c("log_btwdegree","Netflow","log_SelfR","log_Indegree_Neigh")
HUMAN <- c("log_grav_total","_ClassClosed","_ClassRestricted")
ENV <- c("temp","prod.annual")

PARAM <- c("Intercept","Temperature","Productivity",
           "Betweenness","Netflow","Local Recruit.","Exogenous Indegree", 
           "Tot.Gravity","No-Take","Restricted gears")

CAT <- c("Intercept","Human/Env.","Human/Env.",
         "Connectivity",
         "Connectivity",
         "Connectivity",
         "Connectivity",
         "Human/Env.","Human/Env.",
         "Human/Env.")


DesiredOrder <- c("Intercept",
                  "Temperature",
                  "Productivity",
                  "Tot.Gravity",
                  "No-Take",
                  "Restricted gears",
                  "Betweenness",
                  "Netflow",
                  "Local Recruit.",
                  "Exogenous Indegree")

#GLOBAL
rm(m2_df_global)
m2_df_global <- Richness_inflow_GLOBAL_simp_data %>%  mutate(FE = "GLOBAL") %>%  mutate(model = "Richness")
m2_df_global <- cbind(PARAM,m2_df_global,CAT)
m2_df_global$PARAM <- factor(m2_df_global$PARAM, levels = rev(DesiredOrder))
m2_df_global$CAT <- factor(m2_df_global$CAT , levels = c("Intercept","Human/Env.","Connectivity"),ordered = TRUE)

#TRANSIENT
rm(m2_df_transient)
m2_df_transient <- Richness_inflow_TRANSIENT_simp_data %>%  mutate(FE = "TRANSIENT") %>%  mutate(model = "Richness")
m2_df_transient <- cbind(PARAM,m2_df_transient,CAT)
m2_df_transient$PARAM <- factor(m2_df_transient$PARAM, levels = rev(DesiredOrder))
m2_df_transient$CAT <- factor(m2_df_transient$CAT , levels = c("Intercept","Human/Env.","Connectivity"),ordered = TRUE)

#CRYPTIC
rm(m2_df_crypto)
m2_df_crypto <- Richness_inflow_CRYPTO_simp_data %>%  mutate(FE = "PARENTAL") %>%  mutate(model = "Richness")
m2_df_crypto <- cbind(PARAM,m2_df_crypto,CAT)
m2_df_crypto$PARAM <- factor(m2_df_crypto$PARAM, levels = rev(DesiredOrder))
m2_df_crypto$CAT <- factor(m2_df_crypto$CAT , levels = c("Intercept","Human/Env.","Connectivity"),ordered = TRUE)


# brichness data.frame
rm(four_FE_Richness)
four_FE_Richness <- rbind(m2_df_global,m2_df_transient,m2_df_crypto)
four_FE_Richness$FE <- factor(four_FE_Richness$FE, levels=c("PARENTAL","TRANSIENT","GLOBAL"))
four_FE_Richness$CAT <- as.factor(four_FE_Richness$CAT)
four_FE_Richness$model <- as.factor(four_FE_Richness$model)
summary(four_FE_Richness)

library(viridis)
# plot Richness
rm(Richness.FE)
Richness.FE <- ggplot(four_FE_Richness %>% filter(PARAM != "Intercept"), aes(colour = FE)) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_linerange(aes(x = PARAM, ymin = ll,
                     ymax = hh,colour = FE),lwd = 1/2, position = position_dodge(width = 0.8)) +
  geom_point(aes(x = PARAM, y = m,colour = FE, size=6), lwd = 1.5, position = position_dodge(width = 0.8),shape = 21,fill = "WHITE") +
  geom_linerange(aes(x = PARAM, ymin = l,
                     ymax = h,colour = FE),lwd = 1, position = position_dodge(width = 0.8)) +
  facet_grid(CAT ~ ., scales = "free", space = "free") + scale_color_viridis(discrete=TRUE) +
  #scale_x_discrete(labels= PARAM)  +
  coord_flip() + theme_bw() +
  labs(x="",y="Standardized coefficients") + theme(text = element_text(size = 10, family="Arial")) +
  ggtitle("Richness") 
Richness.FE  # The trick to these is position_dodge()


BAYE_coef_tot_rich_biom <- ggarrange(Biomass.FE,Richness.FE,
                                     ncol=2,nrow=1,labels = c("A","B"),align="hv",common.legend = F,legend="bottom")


#part C and D 


#Netflow and indegree of neightbours #back transform netflow, species richness and indegree of neighbours
conditional_effects(MOD_S_1_run_ACTIVE.1, "log_Indegree_Neigh")
##back transform richness
FUNSTEPH = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}
# z = (x - mean)/ (1 * sd (x))
meanRich<-mean(ACTIVE.1.sub.V2$Richness,na.rm=T)
sdRich<- 1*sd(ACTIVE.1.sub.V2$Richness,na.rm=T)
datanorm<- (FUNSTEPH(ACTIVE.1.sub.V2$Richness))
#head(CRYPTIC.std$Richness) #duble checking
#now back transform
FUNINVLU = function(x){(x * sdRich) + meanRich}
head(FUNINVLU(datanorm)) #ok
head(ACTIVE.1.sub.V2$Richness) #ok

##back transform InDEGREE neighbors
meanINfN<-mean(ACTIVE.1.sub.V2$IndegreeNe,na.rm=T)
sdINfN<- 1*sd(ACTIVE.1.sub.V2$IndegreeNe,na.rm=T)
datanormINfN<- (FUNSTEPH(ACTIVE.1.sub.V2$IndegreeNe))
head(datanormINfN)
head(FUNSTEPH(ACTIVE.1.sub.V2$IndegreeNe))
#now back transform
FUNINVLUIF = function(x){(x * sdINfN) + meanINfN}
head(FUNINVLUIF(datanormINfN)) #ok
head(ACTIVE.1.sub.V2$IndegreeNe) #ok

#marginal plot (richness vs extrinsic indegree )

margCry<-plot(conditional_effects(MOD_S_1_run_ACTIVE.1, "log_Indegree_Neigh"))
marRich_Cor<-margCry$log_Indegree_Neigh
try2<-margCry$log_Indegree_Neigh$plot_env$plots$log_Indegree_Neigh$plot_env$x$log_Indegree_Neigh

rm(margPlot)
margPlot<-try2[,c("log_Indegree_Neigh","estimate__","lower__","upper__")]
margPlot$rwRich<-FUNINVLU(margPlot$estimate__)
margPlot$rwRichlow<-FUNINVLU(margPlot$lower__)
margPlot$rwRichupp<-FUNINVLU(margPlot$upper__)
margPlot$unscaled_Indegree_Neigh <- FUNINVLUIF(expm1(margPlot$log_Indegree_Neigh))

PartC<-ggplot(margPlot, aes(y = rwRich, x = log_Indegree_Neigh)) +
  geom_point(size=0.01) +
  geom_ribbon( aes(ymin = rwRichlow, ymax = rwRichupp), alpha = .15) +
  geom_line( aes(y = rwRich), size = 1.5, color="blue") + theme_classic() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))
PartC

##Netflow
conditional_effects(MOD_BIOM_1_run_ACTIVE.1, "Netflow")

##back transform Netflow
meanINfN<-mean(ACTIVE.1.sub.V2$Netflow,na.rm=T)
sdINfN<- 1*sd(ACTIVE.1.sub.V2$Netflow,na.rm=T)
datanormINfN<- (FUNSTEPH(ACTIVE.1.sub.V2$Netflow))
head(datanormINfN)
head(FUNSTEPH(ACTIVE.1.sub.V2$Netflow))
#now back transform
FUNINVLUIF = function(x){(x * sdINfN) + meanINfN}
head(FUNINVLUIF(datanormINfN)) #ok
head(ACTIVE.1.sub.V2$Netflow) #ok

#marginal plot (log biomass vs Netfflow - back transform Netflow )
rm(margCry)
margCry<-plot(conditional_effects(MOD_BIOM_1_run_ACTIVE.1, "Netflow"))
try2<-margCry$Netflow$plot_env$plots$Netflow$plot_env$x$Netflow

rm(margPlot)
margPlot<-try2[,c("Netflow","estimate__","lower__","upper__")]
margPlot$unscaled_Netflow <- FUNINVLUIF(margPlot$Netflow)

rm(PartD)
PartD<-ggplot(margPlot, aes(y = estimate__, x = unscaled_Netflow)) +
  geom_point(size=0.01) +
  geom_ribbon( aes(ymin = lower__, ymax =upper__), alpha = .15) +
  geom_line( aes(y = estimate__), size = 1.5, color="purple") + theme_classic() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
PartD

Netflow_Indegree <- ggarrange(PartD,PartC,
                                     ncol=2,nrow=1,labels = c("A","B"),align="hv")


