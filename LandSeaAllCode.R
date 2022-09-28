############
### HMSC ###
############

library(Hmsc)
library(parallel)
library(corrplot)
library(viridis)
library(ape)
library(Matrix)
library(dplyr)
library(fuzzySim)
library(knitr)

# define environment
setwd("C:/Users/lynds/Dropbox/Seabird Shit moves from land to sea/Final data and code/")

# load phylogenetic tree of study species
phyloTree <- ape::read.tree("Merc_algae_final.tre")
plot(phyloTree)

# load data
Data=read.csv("LandSeaData.csv")
head(Data)

# remove rows with zeros for all species (empty quadrats)
Data <- Data[apply(Data[,23:46], 1, function(x) !all(x==0)),]

# remove duplicate rows for fixed variables
Data <- Data[!duplicated(Data[,c('season','treatment','depth','fetch','runoff','soilN15','m2burrow','aspect','compaction','algaeSiteID','slope')]),]

# remove rows with NAs for fixed variables
Data <- Data[complete.cases(Data[,c('season','treatment','depth','fetch','runoff','soilN15','m2burrow','aspect','compaction','algaeSiteID','slope')]),]

# create fixed variable dataframe with original values
XData_orig = subset(Data, select=c(algaeSiteID,transect,samplePair,treatment,depth,fetch,runoff,soilN15,m2burrow,aspect,compaction,slope,season))
colnames(XData_orig) = c("Plot","Transect","SamplePair","Treatment","Depth","Wave","Runoff","Soild15N","BurrowDensity","Aspect","Compaction","Slope","Season")
head(XData_orig)

# create fixed variable dataframe with mean and centered values
XData <- XData_orig
XData$Depth <- scale(XData$Depth)
XData$Wave <- scale(XData$Wave)
XData$Runoff <- scale(XData$Runoff)
XData$Soild15N <- scale(XData$Soild15N)
XData$BurrowDensity <- scale(XData$BurrowDensity)
XData$Aspect <- scale(XData$Aspect)
XData$Compaction <- scale(XData$Compaction)
XData$Slope <- scale(XData$Slope)
head(XData)

# create community dataframe
Y = as.matrix(Data[,23:46])
head(Y)

# load trait data, with row names as the first column (species names)
TrData = read.csv("Merc_traits.csv", header=T, row.names=1)
head(TrData)

# define covariate formula (include squares of numeric variables for possible intermediate optimums
XFormula = ~  
  poly(Soild15N, degree = 2, raw = TRUE) +
  poly(BurrowDensity, degree = 2, raw = TRUE) +
  Treatment + 
  poly(Runoff, degree = 2, raw = TRUE) + 
  poly(Aspect, degree = 2, raw = TRUE) + 
  poly(Compaction, degree = 2, raw = TRUE) +
  poly(Slope, degree = 2, raw = TRUE) +
  poly(Depth, degree = 2, raw = TRUE) +
  poly(Wave, degree = 2, raw = TRUE) +
  Season 

# define trait formula
TrFormula = ~ Group

# define study design (dataframe of random effects - algae and soil plot pairs)
studyDesign <- as.data.frame(XData$SamplePair)
colnames(studyDesign) = c("Plot")
studyDesign$Plot <- as.factor(studyDesign$Plot)
head(studyDesign)

# define Hmsc random level to be sample pair names without spatial info (set min and max)
NrL = HmscRandomLevel(units=levels(studyDesign[,1]))
NrL$nfMin = 5
NrL$nfMax = 10
NrL

# convert all variables to numeric or factors
str(XData)
XData$Plot <- as.factor(XData$Plot)
XData$Transect <- as.factor(XData$Transect)
XData$SamplePair <- as.factor(XData$SamplePair)
XData$Treatment <- as.factor(XData$Treatment)
XData$Season <- as.factor(XData$Season)
XData$Depth <- as.numeric(XData$Depth)
XData$Wave <- as.numeric(XData$Wave)
XData$Runoff <- as.numeric(XData$Runoff)
XData$Soild15N <- as.numeric(XData$Soild15N)
XData$BurrowDensity <- as.numeric(XData$BurrowDensity)
XData$Aspect <- as.numeric(XData$Aspect)
XData$Compaction <- as.numeric(XData$Compaction)
XData$Slope <- as.numeric(XData$Slope)
str(XData)

# construct model
m = Hmsc(Y = Y, 
         XData = XData, XFormula = XFormula,
         TrData = TrData, TrFormula = TrFormula,
         phyloTree = phyloTree,
         studyDesign = studyDesign,
         ranLevels = list("Plot"=NrL))

# run MCMC and save model
thin = 10
samples = 100
nChains = 3
verbose = 0
m = sampleMcmc(m, samples = samples, thin = thin,
               adaptNf = rep(ceiling(0.4*samples*thin),1),
               transient = ceiling(0.5*samples*thin),
               nChains = nChains, nParallel = nChains,
               verbose = verbose)
save(m, "916_model_3_chains_100_samples_10_thin")

# visualize the phylogenetic correlation matrix C to detect any patterns
image(Matrix(m$C))

# evaluate MCMC convergence
# Beta=species niche
# Gamma=influence of traits on species
# Omega=residual species associations
# Rho=strength of phylogenetic signal
mpost = convertToCodaObject(m)
par(mar=c(5,6,4,1)+.1)
par(mfrow=c(3,2))
ess.beta = effectiveSize(mpost$Beta)
psrf.beta = gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
hist(ess.beta)
hist(psrf.beta)
ess.gamma = effectiveSize(mpost$Gamma)
psrf.gamma = gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf
hist(ess.gamma)
hist(psrf.gamma)
ns=ncol(Y)
sppairs = matrix(sample(x = 1:ns^2, size = 100))
tmp = mpost$Omega[[1]]
for (chain in 1:length(tmp)){
  tmp[[chain]] = tmp[[chain]][,sppairs]
}
ess.omega = effectiveSize(tmp)
psrf.omega = gelman.diag(tmp, multivariate=FALSE)$psrf
hist(ess.omega)
hist(psrf.omega)

# convergence diagnostics
# Reasonably good MCMC convergence: most potential scale reduction factors (psrf) are close to one and effective sample sizes (ess) relatively high
print("ess.rho:")
effectiveSize(mpost$Rho)
print("psrf.rho:")
gelman.diag(mpost$Rho)$psrf

# model fit - measure explanatory power of model (R2)
preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=preds)
MF$R2
hist(MF$R2, xlim = c(0,1), main=paste0("Mean = ", round(mean(MF$R2),2)))

# model R2 by species
r2.sp <- data.frame(colnames(m$Y), MF$R2)
colnames(r2.sp) <- c("Species", "modelR2")
r2.sp
r2.sp %>% arrange(desc(modelR2))

# variance partitioning
# examine design matrix X and assign groups of variables
head(m$X)
VP = computeVariancePartitioning(m)
VP
VP = computeVariancePartitioning(m, group = c(1,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10), groupnames = c("Soil d15N","Burrow Density","Island (eradication history)","Runoff","Aspect","Compaction","Slope","Depth","Wave action","Season"))

# extract 'vals' dataframe from VP to reorder based on seabird values and transpose
vals=VP$vals
valst=as.data.frame(t(vals))
head(valst)

# create and sort new dataframe from sums of seabirds and terrestrial variables
sums=as.data.frame(rowSums(valst[,c("Soil d15N","Burrow Density","Island (eradication history)","Runoff","Aspect","Compaction","Slope")]))
colnames(sums)=c("sums")
sums.sort=sums[order(-sums$sums), , drop=FALSE]
head(sums.sort)

# reorder vals dataframe based on order of sums dataframe and transpose back
valst.sort=valst[match(rownames(sums.sort), rownames(valst)), ]
vals.sort=as.matrix(t(valst.sort))
head(vals.sort)
vals.sort

# replace old vals dataframe with new sorted dataframe within VP object
VP2=VP
VP2$vals=vals.sort
VP2

# plot variance partitioning
# abbreviated names, bold study species, decreasing seabird influence
barcolor=c("#482677FF","#482677FF","#482677FF","#3CBB75FF","#3CBB75FF","#3CBB75FF","#3CBB75FF","#238A8DFF","#238A8DFF","#DCE319FF","grey50")
longname=row.names(sums.sort)
shortname=spCodes(longname, nchar.gen=1, nchar.sp=3, sep.species="_", sep.spcode=". ")
par(mar=c(7,4,3,1) +0.1)
plotVariancePartitioning(m, VP = VP2, cols=barcolor, 
                         legend.text=F,
                         names.arg=shortname, 
                         border="white",
                         cex.axis=2,
                         cex.names=2,
                         cex.main=2,
                         col.lab=adjustcolor("white",alpha.f=0),
                         main="All islands and seasons",
                         las=2)

# table of how much the traits explain out of the variation among species in their responses to environmental covariates
kable(VP$R2T$Beta)

# how the influence of traits on species niches propagates into influence of traits on species abundances
VP$R2T$Y

# plot species' environmental responses mapped to phylogeny
par(mar=c(20,4,3,1) +0.1)
postBeta = getPostEstimate(m, parName = "Beta")
plotBeta(m, post = postBeta, param = "Support",
         plotTree = TRUE, supportLevel = 0.95, split=.4, spNamesNumbers = c(F,F))

# plot how species traits influence env responses
par(mar=c(20,4,3,1) +0.1)
postGamma = getPostEstimate(m, parName = "Gamma")
plotGamma(m, post=postGamma, param="Support", supportLevel = 0.95)

# plot estimated residual associations among species
OmegaCor = computeAssociations(m)
supportLevel = 0.95
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
corrplot(toPlot, method = "color",
         col=colorRampPalette(c("blue","white","red"))(200),
         tl.cex=.6, tl.col="black",
         title=paste("random effect level:", m$rLNames[1]), mar=c(0,0,1,0))

# calculate strength of phylogentic signal (if rho contains zero, no strong evidence for phylo signal)
summary(mpost$Rho)

# plot variation over env gradients
# example: depth gradient
Gradient = constructGradient(m,focalVariable = "Depth",
                             non.focalVariables = list("Treatment"=list(3, "Old Erad"), 
                                                       "Soild15N"=list(3,0), 
                                                       "BurrowDensity"=list(3,0),
                                                       "Runoff"=list(3,0),
                                                       "Aspect"=list(3,0),
                                                       "Compaction"=list(3,0),
                                                       "Slope"=list(3,0),
                                                       "Wave"=list(3,0),
                                                       "Season"=list(3, "winter")))
Gradient$XDataNew
predY = predict(m, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,
                ranLevels=Gradient$rLNew, expected=TRUE)
# S=species abundance
plotGradient(m, Gradient, pred=predY, measure="S", showData = TRUE, cicol=adjustcolor("#238A8DFF",alpha.f=0.7))
# Y=individual species, by index (m$Y)
plotGradient(m, Gradient, pred=predY, measure="Y", index=10, showData = TRUE)

# example: island gradient
Gradient = constructGradient(m,focalVariable = "Treatment",
                             non.focalVariables = list("Soild15N"=list(1), 
                                                       "BurrowDensity"=list(1),
                                                       "Runoff"=list(1),
                                                       "Aspect"=list(1),
                                                       "Compaction"=list(1),
                                                       "Slope"=list(1),
                                                       "Wave"=list(1),
                                                       "Depth"=list(1),
                                                       "Season"=list(1)))
Gradient$XDataNew
predY = predict(m, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,
                ranLevels=Gradient$rLNew, expected=TRUE)
# S=species abundance
plotGradient(m, Gradient, pred=predY, measure="S", showData = TRUE, jigger = 0.2)
# Y=individual species, by index (m$Y)
plotGradient(m, Gradient, pred=predY, measure="Y", index=1, showData = TRUE, jigger = 0.2)




############
### LMMs ###
############

library(nlme)
library(MuMIn)
library(AICcmodavg)

# define environment
setwd("C:/Users/lynds/Dropbox/Seabird Shit moves from land to sea/Final data and code/")

# load data
Data=read.csv("LandSeaData.csv")
Data <- as.data.frame(Data[,1:22])
head(Data)

# remove rows with NAs for fixed variables
Data <- Data[complete.cases(Data[,c('season','treatment','depth','fetch','runoff','soilN15','m2burrow','aspect','compaction','algaeSiteID','slope')]),]

# mean and center
Data$fetch <- scale(Data$fetch)
Data$runoff <- scale(Data$runoff)
Data$soilN15 <- scale(Data$soilN15)
Data$m2burrow <- scale(Data$m2burrow)
Data$aspect <- scale(Data$aspect)
Data$compaction <- scale(Data$compaction)
Data$slope <- scale(Data$slope)
head(Data)

# specify season and transect to factor
Data$season <-as.factor(Data$season)
Data$transect <-as.factor(Data$transect)

# subset by species
Data.Erad <- Data[Data$algaeSpecies=="E. radiata",]
Data.Xchon <- Data[Data$algaeSpecies=="X. chondrophylla",]
Data.Cmas <- Data[Data$algaeSpecies=="C. maschalocarpum",]
Data.Cflex <- Data[Data$algaeSpecies=="C. flexuosum",]
Data.Cplum <- Data[Data$algaeSpecies=="C. plumosum",]
Data.Vcol <- Data[Data$algaeSpecies=="V. colensoi",]

ctrl <- lmeControl(opt='optim')

# determine weighted variables
weight.none <- gls(algaeN15 ~ soilN15 + season + runoff + depth + fetch + m2burrow + slope + compaction + aspect + treatment,
                   data=Data.Erad, control=ctrl)
plot(weight.none, select=c(1), pch=16)
weight.treatment <- gls(algaeN15 ~ soilN15 + season + runoff + depth + fetch + m2burrow + slope + compaction + aspect + treatment,
                        data=Data.Erad, control=ctrl, weights=varIdent(form= ~ 1 | treatment))
plot(weight.treatment, select=c(1), pch=16)
weight.season <- gls(algaeN15 ~ soilN15 + season + runoff + depth + fetch + m2burrow + slope + compaction + aspect + treatment,
                     data=Data.Erad, control=ctrl, weights=varIdent(form= ~ 1 | season))
plot(weight.season, select=c(1), pch=16)
AIC(weight.none, weight.treatment, weight.season)

# determine random variables
random.none <- gls(algaeN15 ~ soilN15 + season + runoff + depth + fetch + m2burrow + slope + compaction + aspect + treatment,
                   data=Data.Erad, control=ctrl, weights=varIdent(form= ~ 1 | treatment))
plot(random.none, select=c(1), pch=16)
random.transect <- lme(algaeN15 ~ soilN15 + season + runoff + depth + fetch + m2burrow + slope + compaction + aspect + treatment,
                       data=Data.Erad, control=ctrl, random= ~1 | transect, weights=varIdent(form= ~ 1 | treatment), method="REML")
plot(random.transect, select=c(1), pch=16)
random.season <- lme(algaeN15 ~ soilN15 + season + runoff + depth + fetch + m2burrow + slope + compaction + aspect + treatment,
                     data=Data.Erad, control=ctrl, random= ~1 | season, weights=varIdent(form= ~ 1 | treatment), method="REML")
plot(random.season, select=c(1), pch=16)
random.season.transect <- lme(algaeN15 ~ soilN15 + season + runoff + depth + fetch + m2burrow + slope + compaction + aspect + treatment,
                              data=Data.Erad, control=ctrl, random= ~1 + season | transect, weights=varIdent(form= ~ 1 | treatment), method="REML")
plot(random.season.transect, select=c(1), pch=16)
AIC(random.none, random.transect, random.season, random.season.transect)

# final model
LM.add.Data.Erad <- lme(algaeN15 ~ soilN15 + season + runoff + depth + fetch + m2burrow + slope + compaction + aspect + treatment,
                       data=Data.Erad, control=ctrl, random= ~1 + season | transect, weights=varIdent(form= ~ 1 | treatment), method="ML")
plot(LM.add.Data.Erad, select=c(1), pch=16)

MinsTabList <- dredge(LM.add.Data.Erad)
MinsTabList
CandList <- subset(MinsTabList, weight>=0.01)
CandList
CandList <- subset(CandList, delta<=6)
CandList
CandMod <- get.models(CandList, subset=TRUE)
CandMod
summary(model.avg(CandMod))
confint(model.avg(CandMod))
CandList
CandModAvg <- model.avg(CandMod)
summary(model.avg(CandMod))
confint(model.avg(CandMod))
aictab(CandMod)

Bopt <- lme(algaeN15 ~ soilN15 + season + runoff + depth + fetch + m2burrow + slope + compaction + aspect + treatment,
            data=Data.Erad, control=ctrl, random= ~1 + season | transect, weights=varIdent(form= ~ 1 | treatment), method="ML")

Data.Erad.dry <- Data.Erad[Data.Erad$season=="Dry",]
Data.Erad.rainy <- Data.Erad[Data.Erad$season=="Rainy",]

newdatB.dry <- expand.grid(soilN15=c(min(Data.Erad.dry$soilN15), max(Data.Erad.dry$soilN15)), 
                           season="Dry", 
                           runoff=c(min(Data.Erad.dry$runoff), max(Data.Erad.dry$runoff)), 
                           depth=c(min(Data.Erad.dry$depth), max(Data.Erad.dry$depth)), 
                           fetch=c(min(Data.Erad.dry$fetch), max(Data.Erad.dry$fetch)),
                           m2burrow=c(min(Data.Erad.dry$m2burrow), max(Data.Erad.dry$m2burrow)),
                           slope=c(min(Data.Erad.dry$slope), max(Data.Erad.dry$slope)),
                           compaction=c(min(Data.Erad.dry$compaction), max(Data.Erad.dry$compaction)),
                           aspect=c(min(Data.Erad.dry$aspect), max(Data.Erad.dry$aspect)),
                           treatment=unique(Data.Erad.dry$treatment))
newdatB.rainy <- expand.grid(soilN15=c(min(Data.Erad.rainy$soilN15), max(Data.Erad.rainy$soilN15)), 
                             season="Rainy", 
                             runoff=c(min(Data.Erad.rainy$runoff), max(Data.Erad.rainy$runoff)), 
                             depth=c(min(Data.Erad.rainy$depth), max(Data.Erad.rainy$depth)), 
                             fetch=c(min(Data.Erad.rainy$fetch), max(Data.Erad.rainy$fetch)),
                             m2burrow=c(min(Data.Erad.rainy$m2burrow), max(Data.Erad.rainy$m2burrow)),
                             slope=c(min(Data.Erad.rainy$slope), max(Data.Erad.rainy$slope)),
                             compaction=c(min(Data.Erad.rainy$compaction), max(Data.Erad.rainy$compaction)),
                             aspect=c(min(Data.Erad.rainy$aspect), max(Data.Erad.rainy$aspect)),
                             treatment=unique(Data.Erad.rainy$treatment))
newdatB <- rbind(newdatB.dry, newdatB.rainy)

newdatB$pred <- predict(CandModAvg, newdatB, level=0)
DesignmatB=model.matrix(~newdatB$soilN15 + newdatB$season + newdatB$runoff + newdatB$depth + newdatB$fetch + 
                          newdatB$m2burrow + newdatB$slope + newdatB$compaction + newdatB$aspect + newdatB$treatment)
predvarB=diag(DesignmatB%*%Bopt$varFix%*%t(DesignmatB))
newdatB$SE=sqrt(predvarB)
newdatB$SE2=sqrt(predvarB+Bopt$sigma^2)
newdatB$LCI <-newdatB$pred-(2* newdatB$SE)
newdatB$UCI <-newdatB$pred+(2* newdatB$SE)
newdatB
newdatB[,"Species"] <- NA
newdatB$Species <- rep("E. radiata",nrow(newdatB))
Data.Erad.pred <- newdatB
head(Data.Erad.pred)















