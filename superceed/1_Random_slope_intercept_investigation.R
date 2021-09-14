


# packages
require(dplyr)
require(here)
require(forcats)
require(brms)
require(rstan)
require(ggplot2)
require(ggpubr)
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

# investigate raw relationships gravity and region for random intercept and slopes
  # intercept
SP_region_gravity <- ggplot(data=PredictVar,aes(x=log_grav_total,y=Richness,col=region)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="Total gravity - Intercept & Slope")
SP_region_gravity
  # intercept and slope
Biom_region_gravity <- ggplot(data=PredictVar,aes(x=log_grav_total,y=log_biomassarea1,col=region)) +
  geom_point() +
  geom_smooth(method="lm")+
  labs(title="Total gravity - Intercept & Slope")
Biom_region_gravity

# investigate raw relationships gravity neighbour and region for random intercept and slopes
# intercept
SP_region_gravity.ne <- ggplot(data=PredictVar,aes(x=log_grav_neiBR,y=Richness,col=region)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="Gravity neighbour - Intercept")
SP_region_gravity.ne

# intercept and slope
Biom_region_gravity.ne <- ggplot(data=PredictVar,aes(x=log_grav_neiBR,y=log_biomassarea1,col=region)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="Gravity neighbour - Intercept & slope")
Biom_region_gravity.ne

# investigate raw relationships InflowBR and region for random intercept and slopes
  # intercept
SP_region_Inflow <- ggplot(data=PredictVar,aes(x=log(InflowBR+1),y=Richness,col=region)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="Inflow - Intercept")
SP_region_Inflow

  # slope and intercept
Biom_region_Inflow <- ggplot(data=PredictVar,aes(x=log(InflowBR+1),y=log_biomassarea1,col=region)) +
  geom_point() +
  geom_smooth(method="lm")  +
  labs(title="Inflow -Slope & Intercept")
Biom_region_Inflow

# investigate raw relationships InflowBR MPA and region for random intercept and slopes
# intercept
SP_region_Inflow_MPA <- ggplot(data=PredictVar,aes(x=log(InflowMPABR+1),y=Richness,col=region)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="Inflow MPA - Intercept")
SP_region_Inflow_MPA

# slope and intercept
Biom_region_Inflow_MPA <- ggplot(data=PredictVar,aes(x=log(InflowMPABR+1),y=log_biomassarea1,col=region)) +
  geom_point() +
  geom_smooth(method="lm")  +
  labs(title="Inflow MPA -Slope & Intercept")
Biom_region_Inflow_MPA

# investigate raw relationships Indegree and region for random intercept and slopes
# intercept only
SP_region_Indegree <- ggplot(data=PredictVar,aes(x=log(IndegreeBR+1),y=Richness,col=region)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="Indegree -Intercept")
SP_region_Indegree

# intercept and slope
Biom_region_Indegree <- ggplot(data=PredictVar,aes(x=log(IndegreeBR+1),y=log_biomassarea1,col=region)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="Indegree -Intercept & Slope")
Biom_region_Indegree

# investigate raw relationships Selfrecruit and region for random intercept and slopes
  # intercept only
SP_region_SR <- ggplot(data=PredictVar,aes(x=log(SelfR+1),y=Richness,col=region)) +
  geom_point() +
  geom_smooth(method="lm")  +
  labs(title="Self Retention -Intercept")
SP_region_SR

  # intercept and slope
Biom_region_SR <- ggplot(data=PredictVar,aes(x=log(SelfR+1),y=log_biomassarea1,col=region)) +
  geom_point() +
  geom_smooth(method="lm")  +
  labs(title="Self Retention -Slope and Intercept")
Biom_region_SR

# investigate raw relationships Indegree MPA and region for random intercept and slopes
# intercept only
SP_region_Indegree.MPA <- ggplot(data=PredictVar,aes(x=log(IndegreeMPABR+1),y=Richness,col=region)) +
  geom_point() +
  geom_smooth(method="lm")  +
  labs(title="Indegree MPA -Intercept")
SP_region_Indegree.MPA

# intercept and slope
Biom_region_Indegree.MPA <- ggplot(data=PredictVar,aes(x=log(IndegreeMPABR+1),y=log_biomassarea1,col=region)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="Indegree MPA -Slope and Intercept")
Biom_region_Indegree.MPA

# investigate raw relationships Indegree NEIG and region for random intercept and slopes
# intercept only
SP_region_Indegree.nei <- ggplot(data=PredictVar,aes(x=log(IndegreeNeiBR+1),y=Richness,col=region)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="Indegree Nei - Intercept")
SP_region_Indegree.nei

# intercept and slope
Biom_region_Indegree.nei <- ggplot(data=PredictVar,aes(x=log(IndegreeNeiBR+1),y=log_biomassarea1,col=region)) +
  geom_point() +
  geom_smooth(method="lm")  +
  labs(title="Indegree Nei - Slope and Intercept")
Biom_region_Indegree.nei

# investigate raw relationships Inflow NEIG and region for random intercept and slopes
# intercept only
SP_region_Inflow.nei <- ggplot(data=PredictVar,aes(x=log(InflowNeiBR+1),y=Richness,col=region)) +
  geom_point() +
  geom_smooth(method="lm")  +
  labs(title="Inflow Nei - Intercept")
SP_region_Inflow.nei

# intercept and slope
Biom_region_Inflow.nei <- ggplot(data=PredictVar,aes(x=log(InflowNeiBR+1),y=log_biomassarea1,col=region)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="Inflow Nei - Slope and Intercept")
Biom_region_Inflow.nei

# investigate raw relationships gravity and region for random intercept and slopes

# intercept and slope
Biom_region_gravity_inflow <- ggplot(data=PredictVar,aes(x=log_grav_total,y=log(InflowBR+1),col=region)) +
  geom_point() +
  geom_smooth(method="lm")+
  labs(title="Total gravity - Intercept & Slope")
Biom_region_gravity_inflow 

#####################
#     SAVE PLOTS    #
#####################

# save gravity plots
Gravity <- ggarrange(SP_region_gravity,Biom_region_gravity,
                     SP_region_gravity.ne,Biom_region_gravity.ne,
                     ncol=2,nrow=2,
                     labels = c("A","B","C","D"),common.legend = T,legend="right")
ggsave(here("_prelim.figures","Gravity_random_slope_intercept.pdf"),Gravity,width=10,height=8)

# save inflow plots
Inflow <- ggarrange(SP_region_Inflow,Biom_region_Inflow,
                     SP_region_Inflow.nei,Biom_region_Inflow.nei,
                    SP_region_Inflow_MPA,Biom_region_Inflow_MPA,
                     ncol=2,nrow=3,
                     labels = c("A","B","C","D","E","F"),common.legend = T,legend="right")
ggsave(here("_prelim.figures","Inflow_random_slope_intercept.pdf"),Inflow,width=10,height=8)

# save indegtree plots
Indegree <- ggarrange(SP_region_Indegree,Biom_region_Indegree,
                    SP_region_Indegree.nei,Biom_region_Indegree.nei,
                    SP_region_Indegree.MPA,Biom_region_Indegree.MPA,
                    ncol=2,nrow=3,
                    labels = c("A","B","C","D","E","F"),common.legend = T,legend="right")
ggsave(here("_prelim.figures","Indegree_random_slope_intercept.pdf"),Indegree,width=10,height=8)

# self retention
SR <- ggarrange(SP_region_SR,Biom_region_SR,
                      ncol=2,nrow=1,
                      labels = c("A","B"),common.legend = T,legend="right")
ggsave(here("_prelim.figures","SelfRetention_random_slope_intercept.pdf"),SR,width=10,height=8)




# plot all variables together
library(GGally)
library(tidyr)
library(stringr)
library(tidyverse)
# fonction to plot correlation with lm and loess function
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}

# correlation between all functional groups and  scnearios with biplots  
nums <- unlist(lapply(PredictVar, is.numeric))  
PredictVar.num <- PredictVar[,nums]
g = ggpairs(PredictVar.num,lower = list(continuous = my_fn))
g
ggsave(here("_prelim.figures","correl_connectivity_gravity.pdf"),plot=g,width=20,height=20)



