### script to load data for the connectivity project when necessary


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
#install.packages(c("dplyr","here","forcats","corrgram","ggpubr"))

require(dplyr)
require(here)
require(forcats)
library(corrgram) # for corrgram
library(ggpubr)


# load data
rm(all.data)
all.data<-read.csv(here("_data","Connectivity_Biomass_SEMGLMMDATA_March2021.csv"),h=T, stringsAsFactors = F,dec=".")

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

# log all the data
all.data$log_grav_total <- log(all.data$grav_total+1)
all.data$log_grav_neiBR <- log(all.data$grav_nei+1)
all.data$log_biomassarea <-log(all.data$biomassare+1)

# chage to factor
all.data$region <- as.factor(all.data$region)
all.data$locality <- as.factor(all.data$locality)
all.data$sites <- as.factor(all.data$sites)
all.data$Class <- as.factor(all.data$Class)
all.data$ModelMode <- as.factor(all.data$ModelMode)
all.data$Larval_behaviour <- as.factor(all.data$Larval_beh)
all.data$FE <- as.factor(all.data$FE)

# Transform connection with MPAs
all.data$IndegreeMPA.bin <- ifelse(all.data$IndegreeMP == 0,0,1) %>% as.factor()
all.data$IndegreeMPA.bin %>% summary()

  # summarize productivity to annual productivity
all.data <- all.data%>%   
  rowwise() %>%
  mutate(prod.annual = mean(c(prod.Jan:prod.Dec),na.rm=T)) %>%
  as.data.frame()

# rm NAs for netflow
all.data <- all.data %>% filter(!is.na(Netflow))

# log transformed skewed data
all.data$log_btwdegree <- log(all.data$btwdegree+1)
all.data$log_SelfR <- log(all.data$SelfR+1)
all.data$log_CorridorIn <- log(all.data$CorridorIn+1)
all.data$log_InflowMPA <- log(all.data$InflowMPA+1)
all.data$log_InflowNei <- log(all.data$InflowNei+1)
all.data$log_Inflow <- log(all.data$Inflow+1)
all.data$log_annual_prod <- log(all.data$prod.annual+1)

# save all.data file
saveRDS(all.data,here::here("_data","Connectivity_Biomass_SEMGLMMDATA_March2021.rds"))

# correlation between connectivity attributes
connectivity <- all.data[,c("SelfR","Inflow","Indegree",        
                            "CorridorIn","IndegreeMP","InflowMPA",       
                            "IndegreeNe","InflowNei","OutFlow","Outdegree","btwdegree","Netflow")] 
corr.connectivity <- corrgram(connectivity,order=TRUE, lower.panel=panel.shade,
                              upper.panel=panel.cor, text.panel=panel.txt)

#ggexport(corr.connectivity,filename=here("_prelim.figures","Correlations","Corr_connectivity.pdf"),width=20,height=12)

# correlation between environmental + human attributes
env_human <- all.data[,c("log_grav_total","log_grav_neiBR","temp","prod.annual","Age_of_pro")] 
corr.env_human <- corrgram(env_human,order=TRUE, lower.panel=panel.shade,
                              upper.panel=panel.cor, text.panel=panel.txt)

#ggexport(corr.env_human,filename=here("_prelim.figures","Correlations","Corr_env_human.pdf"),width=20,height=12)

ggplot(all.data,aes(x=log_grav_total,y=log_grav_neiBR)) +
  geom_point() +
  geom_smooth()

# Transient
TRANSIENT %>% rm()
TRANSIENT <- all.data %>% filter(Larval_behaviour == "active" & ModelMode == "transi15") %>% droplevels()

# Fished as the reference
TRANSIENT$Class <- relevel(TRANSIENT$Class, ref="Fished")
summary(TRANSIENT)

## standrdize x variables
rm(TRANSIENT.std)
TRANSIENT.std<-data.frame(apply(X = TRANSIENT[,c("temp","Age_of_pro","prod.annual",
                                                 "Netflow","log_grav_total","log_grav_neiBR",
                                                 "log_btwdegree","log_SelfR","log_CorridorIn","log_Inflow",  
                                                 "log_InflowMPA","log_InflowNei")], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))
# add management and region
TRANSIENT.std <- cbind(TRANSIENT$region,TRANSIENT.std)
TRANSIENT.std <- cbind(TRANSIENT$Class,TRANSIENT.std)
colnames(TRANSIENT.std)[c(1,2)] <- c("Class","region")

# add biomass and richness
TRANSIENT.std$log_biomassarea <- TRANSIENT$log_biomassarea
TRANSIENT.std$Richness <- TRANSIENT$Richness


ggplot(TRANSIENT.std,aes(x=log_InflowNei,y=log_biomassarea)) +
  geom_point() +
  geom_smooth(method="gam")


head(TRANSIENT.std)
dim(TRANSIENT.std)
summary(TRANSIENT.std)
names(TRANSIENT.std)# add log biomass


# Parental
PARENTAL %>% rm()
PARENTAL <- all.data %>% filter(Larval_behaviour == "active" & ModelMode == "pare5") %>% droplevels()

# Fished as the reference
PARENTAL$Class <- relevel(PARENTAL$Class, ref="Fished")
summary(PARENTAL)

## standardize x variables
rm(PARENTAL.std)
PARENTAL.std<-data.frame(apply(X = PARENTAL[,c("temp","Age_of_pro","prod.annual",
                                               "Netflow","log_grav_total","log_grav_neiBR",
                                               "log_btwdegree","log_SelfR","log_CorridorIn","log_Inflow",  
                                               "log_InflowMPA","log_InflowNei")], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))

# add management and region
PARENTAL.std <- cbind(PARENTAL$region,PARENTAL.std)
PARENTAL.std <- cbind(PARENTAL$Class,PARENTAL.std)
colnames(PARENTAL.std)[c(1,2)] <- c("Class","region")

# add biomass and richness
PARENTAL.std$log_biomassarea <- PARENTAL$log_biomassarea
PARENTAL.std$Richness <- PARENTAL$Richness

# Cryptic
CRYPTIC %>% rm()
CRYPTIC <- all.data %>% filter(Larval_behaviour == "active" & ModelMode == "crypto5") %>% droplevels()

# Fished as the reference
CRYPTIC$Class <- relevel(CRYPTIC$Class, ref="Fished")
summary(CRYPTIC)

## standrdize x variables
rm(CRYPTIC.std)
CRYPTIC.std<-data.frame(apply(X = CRYPTIC[,c("temp","Age_of_pro","prod.annual",
                                             "Netflow","log_grav_total","log_grav_neiBR",
                                             "log_btwdegree","log_SelfR","log_CorridorIn","log_Inflow",  
                                             "log_InflowMPA","log_InflowNei")], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))

# add management and region
CRYPTIC.std <- cbind(CRYPTIC$region,CRYPTIC.std)
CRYPTIC.std <- cbind(CRYPTIC$Class,CRYPTIC.std)
colnames(CRYPTIC.std)[c(1,2)] <- c("Class","region")

# add biomass and richness
CRYPTIC.std$log_biomassarea <- CRYPTIC$log_biomassarea
CRYPTIC.std$Richness <- CRYPTIC$Richness

# Resident
RESID %>% rm()
RESID <- all.data %>% filter(Larval_behaviour == "active" & ModelMode == "resid15") %>% droplevels()

# Fished as the reference
RESID$Class <- relevel(RESID$Class, ref="Fished")
summary(RESID)

## standardize x variables
rm(RESID.std)
RESID.std<-data.frame(apply(X = RESID[,c("temp","Age_of_pro","prod.annual",
                                         "Netflow","log_grav_total","log_grav_neiBR",
                                         "log_btwdegree","log_SelfR","log_CorridorIn","log_Inflow",  
                                         "log_InflowMPA","log_InflowNei")], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (1*sd(x,na.rm=T))}))

# add management and region
RESID.std <- cbind(RESID$region,RESID.std)
RESID.std <- cbind(RESID$Class,RESID.std)
colnames(RESID.std)[c(1,2)] <- c("Class","region")

# add biomass and richness
RESID.std$log_biomassarea <- RESID$log_biomassarea
RESID.std$Richness <- RESID$Richness

