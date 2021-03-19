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

# save all.data file
saveRDS(all.data,here::here("_data","Connectivity_Biomass_SEMGLMMDATA_March2021.rds"))

# correlation between connectivity attributes
connectivity <- all.data[,c("SelfR","Inflow","Indegree",        
                            "CorridorIn","IndegreeMP","InflowMPA",       
                            "IndegreeNe","InflowNei","OutFlow","Outdegree","btwdegree","Netflow")] 
corr.connectivity <- corrgram(connectivity,order=TRUE, lower.panel=panel.shade,
                              upper.panel=panel.cor, text.panel=panel.txt)

ggexport(corr.connectivity,filename=here("_prelim.figures","Correlations","Corr_connectivity.pdf"),width=20,height=12)

# correlation between environmental + human attributes
env_human <- all.data[,c("log_grav_total","log_grav_neiBR","temp","prod.annual","Age_of_pro")] 
corr.env_human <- corrgram(env_human,order=TRUE, lower.panel=panel.shade,
                              upper.panel=panel.cor, text.panel=panel.txt)

ggexport(corr.env_human,filename=here("_prelim.figures","Correlations","Corr_env_human.pdf"),width=20,height=12)

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

