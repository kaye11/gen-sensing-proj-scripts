#read data
trackdata <- read.table ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/final track data.csv", 
                         sep=";", header=T)

library(ggplot2)
library(gtable)
library(ggthemes)
library(mgcv)
library(ggplot2)
library(grid)
s=mgcv:::s


source("AED.R")
source("vif.R")
source("summarySE.R")
source("resizewin.R")

resize.win(15,12)

#treat time as factor
trackdata$timef <- as.factor(trackdata$time)
trackdata$time <- as.numeric(trackdata$time)

#combine A+B
trackdata$newbin <- trackdata$bin
levels(trackdata$newbin) <- c("BinsAB", "BinsAB", "BinC")

#time2
tm=seq(0, 3600, by = 120)
trackdata$time2 <- cut(trackdata$T, tm, labels=paste(tail(tm, -1L)))

trackdata$time2 = factor (trackdata$time2, levels=c(levels(trackdata$time2), 0))

trackdata$time2[is.na(trackdata$time2)] <- 0 #replace NAs with 0. some data time points have the starting point as NA because they started at 0
trackdata$time2n <- as.numeric (trackdata$time2)*120

#summaries
speed.sumV <- summarySE(trackdata, measurevar="V", groupvars=c("treatment", "newbin", "time2"))
speed.sumVlog <- summarySE(trackdata, measurevar="Vlog", groupvars=c("treatment", "newbin"))


#Bin A
BinC= subset (trackdata, newbin=='binC')

BinC <- BinC [! BinC$time=="0",  ]

qplot(timef, Vlog, color = treatment, data = BinC,  geom = "boxplot") + facet_grid(treatment~., scales="free") 

expA=as.data.frame(data.table(cbind(treatment=as.numeric(as.factor(BinC$treatment)), T=as.numeric(BinC$time), ID=as.numeric(as.factor(BinC$ID)))))

cor(expA, method = "spearman")

vif_func(in_frame=expA,thresh=5,trace=T)

pairs(expA, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

SpeedBinC <- summarySE(BinC, measurevar="V", groupvars=c("treatment", "time", "bin"))

#boxplots
op=par(mfrow=c(2,2))
boxplot(Vlog~treatment, data=BinC)
boxplot(Vlog~ID, data=BinC)
boxplot (Vlog~time, data=BinC)

#levene
library(lawstat)
levene.test(BinC$Vlog, group=BinC$ID, location="mean") #unequal
levene.test(BinC$Vlog, group=BinC$time, location="mean") #unequal
levene.test(BinC$Vlog, group=BinC$treatment, location="mean") #unequal

#gamm
BA <- gamm (Vlog~s(time, by=treatment, bs="fs"), method="REML", data = BinC)

#smoothing splines for factors can are 3 types cr (cubic regression), cs (shrinkage version of cr) and cc (cyclic cubic)
BA1 <- gamm (Vlog~s(time, by=treatment, bs="fs", xt="cr"), method="REML", data = BinC) 
BA2 <- gamm (Vlog~s(time, by=treatment, bs="fs", xt="cs"), method="REML", data = BinC) 
BA2.1 <- gamm (Vlog~s(time, by=treatment, bs="fs", xt="cc"), method="REML", data = BinC) #best but use cr because it yields better results with corAR1

anova(BA$lme, BA1$lme, BA2$lme, BA2.1$lme)

#make random factor and correlations
fBinC <- Vlog~s(time, by=treatment, bs="fs", xt="cr")

BA3 <- gamm (fBinC, method="REML",  random=list(ID=~1), data = BinC) 
BA4 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), data = BinC) #BEST
BA5 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (), data = BinC) #same with BA4

anova(BA$lme, BA1$lme, BA2$lme, BA2.1$lme, BA3$lme, BA4$lme, BA5$lme)

#make variance structures
#BA6 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varIdent(form=~1| time), data = BinC) #no convergence

#BA7 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varIdent(form=~1| ID), data = BinC) #no convergence

BA8 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), 
             weights = varIdent(form=~1| treatment), data = BinC) 

BA9 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varExp(form=~fitted(.)), data = BinC) 

BA10 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varExp(0.2, form=~treatment), data = BinC) 


#not possible BA10 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights=varComb(varIdent(form=~1|treatment), varExp(form=~fitted(.))), data = BinC) 

#BA9 worked but kind of weirdly

anova(BA$lme, BA1$lme, BA2$lme, BA2.1$lme, BA3$lme, BA4$lme, BA5$lme, BA8$lme)

AIC(BA$lme, BA1$lme, BA2$lme, BA2.1$lme, BA3$lme, BA4$lme, BA5$lme, BA8$lme, BA9$lme)
