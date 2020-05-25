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
#_____________________________________
#NOTES: model construction follows this convention:
#Random effects are specified as x|g, where x is an effect and g is a grouping factor (which must be a fac- tor variable, or a nesting of/interaction among factor variables). For example, the formula would be 1|block for a random-intercept model or time|block for a model with random variation in slopes through time across groups specified by block. A model of nested random effects (block within site) would be 1|site/block; a model of crossed random effects (block and year) would be (1|block)+(1|year).
#_____________________________________

source('functions_analyses_glmmTMB.R')')

all.data<-read.csv("_data/FullDataMay2020_LF.csv")

colnames(all.data)
options(stringsAsFactors = FALSE)
#PredictVar<-all.data[,c("temp","Richness","grav_total","Age_of_protection","Indegree","btwdegree","Inflow","Outdegr#ee","InflowLR","SelfR","Class","FE")]

PredictVar<-all.data[,c("temp","Richness","grav_neiBR","Age_of_protection","IndegreeBR","btwdegree","InflowBR","Outdegree","InflowLR","SelfR","IndegreeMPABR","CorridorIndegreeBR","grav_neiBR","InflowMPABR","IndegreeNeiBR","InflowNeiBR","Class","FE")]

##standrdize 16 predicctor variables
data.std<-data.frame(apply(X = PredictVar[,1:16], MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}))

#convert facrtor to character
#bob[] <- lapply(bob, as.character)


data.std1<-cbind(data.std,all.data[,"biomassarea1" ],all.data[,c("Class","Larval_behaviour","FE","ModelMode")])

colnames(data.std1)[17]<-"biomassarea1"

#1. Create list with all possible combinations between predictors 
vifPredCombinations  <-  list()
varnames<-colnames(PredictVar)#

maxCombs  <-  getMaximumNOfCombs(varnames)
for(j in 1:maxCombs) {
vifPredCombinations  <-  append(runPredCombinations(j, varnames), vifPredCombinations)
}

##2. Create list with all possible combinations between predictors 
vifPredCombinations  <-  list()
varnames<-colnames(PredictVar)#

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
modelText.biomass<-lapply(vifPredCombinations_new, prepareModelText,data.std1 )

#set reference level for categorical variable
data.std1$Class<-relevel( as.factor(data.std1$Class), ref="Fished" )

modList<-lapply(modelText.biomass, evalTextModel)

findNonConverge<-lapply(modList, AIC)
nonconv.index<-which(is.na(findNonConverge))
modList1<- modList[-nonconv.index]
modList2<- modList1[-1]

#modelSel<-model.sel(modList1, rank.args = list(REML = FALSE), extra =c(AIC, BIC))
modelSel1<-model.sel(modList2, rank.args = list(REML = FALSE),extra = list(AIC, BIC,R2 = function(x) r.squaredGLMM(x, fmnull)["delta", ]))
write.csv(modelSel1, 'modelSel.biom1.csv')


#top.model<-get.models(modelSel, subset=delta<2)
top.model<-get.models(modelSel1, subset=delta<2)

topModelAve<-model.avg(top.model) 


mA<-summary(topModelAve) #pulling out model averages
df1<-as.data.frame(mA$coefmat.full) #selecting full model coefficient averages

CI <- as.data.frame(confint(topModelAve, full=T)) # get confidence intervals for full model
df1$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df1$CI.max <-CI$`97.5 %`# order of coeffients same in both, so no mixups; but should check anyway
setDT(df1, keep.rownames = "coefficient") #put rownames into column
names(df1) <- gsub(" ", "", names(df1)) # remove spaces from column headers
df1$coefficient<-gsub("cond\\(|)","",x)#remove brackets around predictor names

myPath<-"/Users/josephmaina/Documents/Mygitprojects/2020_Connectivity_Biomass/_prelim.figures/"
#pdf(file = myPath, onefile = F, width = 4, height = 8.5)

df1[2:13,] %>% mutate(Color = ifelse( Estimate> 0, "blue", "red")) %>%
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



ggsave("TopModelAvgCoef.pdf",path = myPath,width = 8, height = 8)
unlink("TopModelAvgCoef.pdf")

save.image("modelSelection.RData")



