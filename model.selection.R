##Maina
##model selection using glmmtmb
##procdure includes testing for VIF
##creating all possible combinations 
##fitting linear models of gamma family

rm(list=ls())
library(usdm)
library(glmmTMB)
library(ggplot2)
library(data.table)
library(dplyr)
library(MuMIn)
library(parallel)

library(ape) #Version 3.3
library(caper) # Vresion 0.5.2
library(nlme) # Version 3.1.122
library(lavaan) # Version 0.5.19
# Load piecewiseSEM from CRAN
library(piecewiseSEM) # Version 1.0.0


#_____________________________________
#NOTES: model construction follows this convention:
#Random effects are specified as x|g, where x is an effect and g is a grouping factor (which must be a fac- tor variable, or a nesting of/interaction among factor variables). For example, the formula would be 1|block for a random-intercept model or time|block for a model with random variation in slopes through time across groups specified by block. A model of nested random effects (block within site) would be 1|site/block; a model of crossed random effects (block and year) would be (1|block)+(1|year).
#_____________________________________

source('functions_analyses_glmmTMB.R')

all.data<-read.csv("_data/DataBiomassConnectivityBR.csv")

colnames(all.data)
options(stringsAsFactors = FALSE)

PredictVar<-all.data[,c("temp","Richness","grav_total","Age_of_protection","SelfR","InflowBR", "IndegreeBR","CorridorIndegreeBR","grav_neiBR", "IndegreeMPABR","InflowMPABR","IndegreeNeiBR","InflowNeiBR","OutFlow","Outdegree","btwdegree","Class","FE")]

PredictVar_R<-all.data[,c("temp","grav_total","Age_of_protection","SelfR","InflowBR", "IndegreeBR","CorridorIndegreeBR","grav_neiBR", "IndegreeMPABR","InflowMPABR","IndegreeNeiBR","InflowNeiBR","OutFlow","Outdegree","btwdegree","Class","FE")]

PredictVar[,17:18] <- sapply(PredictVar[,17:18], function(x) if (is.factor(x)) as.character(x) else x)
PredictVar_R[,16:17] <- sapply(PredictVar_R[,16:17], function(x) if (is.factor(x)) as.character(x) else x)

##standrdize 16 predicctor variables

data.std<-data.frame(apply(X = PredictVar[,1:16], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}))
#data.std<-data.frame(apply(X = PredictVar_R[,1:15], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}))


#convert facrtor to character
#bob[] <- lapply(bob, as.character)

#add spatial covariance
all.data$pos <- numFactor(all.data$Lon, all.data$Lat)

all.data$group <- factor(rep(1, nrow(all.data)))

all.dataChar<-all.data[,c("region","Class","Larval_behaviour","FE","ModelMode")]

all.dataChar[] <- sapply(all.dataChar, function(x) if (is.factor(x)) as.character(x) else x)

data.std1<-cbind(data.std,all.data[,c("Richness","biomassarea1")],all.dataChar)

colnames(data.std1)[17]<-"Richness_resp"
colnames(data.std1)[18]<-"Biomass_resp"

#1. Create list with all possible combinations between predictors 
vifPredCombinations  <-  list()
varnames<-colnames(PredictVar)
#varnames<-colnames(PredictVar_R)
##
maxCombs  <-  getMaximumNOfCombs(varnames)
for(j in 1:maxCombs) {
vifPredCombinations  <-  append(runPredCombinations(j, varnames), vifPredCombinations)
}

##2. Create list with all possible combinations between predictors 
vifPredCombinations  <-  list()
varnames<-colnames(PredictVar)#
#varnames<-colnames(PredictVar_R)

maxCombs  <-  getMaximumNOfCombs(varnames)
for(j in 1:maxCombs) {
vifPredCombinations  <-  append(runPredCombinations(j, varnames), vifPredCombinations)
}

##3. filter the combinations above with VIF<1.5
vifPredCombinations_new<- c()
for(con in vifPredCombinations){
r <- subset(PredictVar, select = con)
conClasses   <-  unique(sapply(r, class))
numOfClasses  <-  length(conClasses)
twoNCols <- ncol(r)==2
numDf <- r[sapply(r,is.numeric)]
zeroNumDf<-ncol(numDf)==0
numeriCols<- ncol(numDf)>1

if(length(con)<2){
vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
}else{
if (zeroNumDf) { 
vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
}
if (numOfClasses==2 && twoNCols) { 
vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
}
if(numeriCols && max(vif(numDf)["VIF"])<=1.5){##vif cutoff
vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
}
next
}
}

##configure models
vifPredCombinations_biom<-vifPredCombinations_new
#vifPredCombinations_rich<-vifPredCombinations_new
rm(vifPredCombinations_new)

modelText.biom<-lapply(vifPredCombinations_biom, prepareModelText,data.std1 )

modelText.rich<-lapply(vifPredCombinations_rich, prepareModelText,data.std1 )

#set reference level for categorical variable
data.std1$Class<-relevel( as.factor(data.std1$Class), ref="Fished" )

detectCores()

##run using single core
#modList.biom<-lapply(modelText.biom, evalTextModel)
#modList.rich<-lapply(modelText.rich, evalTextModel)

##run using multicore
system.time({modList.biom<-mclapply(modelText.biom,mc.cores=24,evalTextModel)})

#system.time({modList.rich<-mclapply(modelText.rich,mc.cores=24,evalTextModel)})

findNonConverge<-lapply(modList.biom, AIC)
nonconv.index<-which(is.na(findNonConverge))
modList.biom<- modList.biom[-nonconv.index]
#modList2<- modList1[-1]

#modelSel<-model.sel(modList1, rank.args = list(REML = FALSE), extra =c(AIC, BIC))
#modelSel1<-model.sel(modList2, rank.args = list(REML = FALSE),extra = list(AIC, BIC,R2 = function(x) r.squaredGLMM(x, fmnull)["delta", ]))
modelSel.biom<-model.sel(modList.biom, rank.args = list(REML = FALSE),extra = list(AIC, BIC,R2 = function(x) r.squaredGLMM(x, fmnull)["delta", ]))
write.csv(modelSel.biom, 'modSelSel.biomass_june_24.csv')


#top.model<-get.models(modelSel, subset=delta<2)
top.model.biom<-get.models(modelSel.biom, subset=delta<2)

topModelAve.b<-model.avg(top.model.biom) 


mA<-summary(topModelAve.b) #pulling out model averages
df1<-as.data.frame(mA$coefmat.full) #selecting full model coefficient averages

CI <- as.data.frame(confint(topModelAve.b, full=T)) # get confidence intervals for full model
df1$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df1$CI.max <-CI$`97.5 %`# order of coeffients same in both, so no mixups; but should check anyway
setDT(df1, keep.rownames = "coefficient") #put rownames into column
names(df1) <- gsub(" ", "", names(df1)) # remove spaces from column headers
df1$coefficient<-gsub("cond\\(|)","",x)#remove brackets around predictor names

myPath<-"~/Documents/Connectvity_Biomass/2020_Connectivity_Biomass/_prelim.figures/"
#pdf(file = myPath, onefile = F, width = 4, height = 8.5)

df1[2:14,] %>% mutate(Color = ifelse( Estimate> 0, "blue", "red")) %>%
ggplot(aes(x=coefficient, y=Estimate,color = Color))+ #again, excluding intercept because estimates so much larger
geom_hline(yintercept=0, color = "black",linetype="dashed", lwd=1.5)+ #add dashed line at zero
geom_errorbar(aes(ymin=CI.min, ymax=CI.max), colour="black", #CI
width=0, lwd=1.5) + coord_flip()+ # flipping x and y axes
geom_point(size=8)+theme_classic(base_size = 20)+ ylab("Coefficient") +
scale_color_identity()

#geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error), colour ="red", # SE
#             width=.2, lwd=3) +
# geom_errorbar(aes(ymin=Estimate-AdjustedSE, ymax=Estimate+AdjustedSE), colour="blue", #adj SE
#               width=.2, lwd=2) +
# geom_errorbar(aes(ymin=CI.min, ymax=CI.max), colour="pink", # CIs
#              width=.2,lwd=1) 

ggsave("TopModelAvgCoef_biomass_june24.pdf",path = myPath,width = 8, height = 8)
unlink("TopModelAvgCoef_biomass_june24.pdf")

write.csv(df1, 'model.averaged.coefficients.csv')

#save worksopace to Luisa's drive
save.image("/Volumes/LuisaDrive/ModelSel/modelSelection_biomass.RData")


load('~/Documents/Mygitprojects/localfiles/modelSelection_biomass.RData')



