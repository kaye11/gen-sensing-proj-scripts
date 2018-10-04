#libraries
library(ggplot2)
library(grid)
library(gtable)
library(ggthemes)
library(gridExtra)
library(mgcv)
library(data.table)
library(nlme)
library(lawstat)

#source files
source("vif.R")
source ("AED.R")
source("resizewin.R")
source("summarySE.R")
resize.win (9,6)

#read in data
starveday <- read.csv("D:/Karen's/PhD/R program/initial data/Data in csv/starveday2.csv", sep=";")

#make new variables
starveday$reptreat <- as.factor(paste(starveday$replicate, starveday$treatment, sep = "-"))
starveday$reptreatbin <- as.factor(paste(starveday$reptreat, starveday$Bin, sep = "-"))
starveday$reptreatbindays <- as.factor(paste(starveday$reptreatbin, starveday$starvedays, sep = "-"))
starveday$reptreatdays <- as.factor(paste(starveday$reptreat, starveday$starvedays, sep = "-"))
starveday$treatdays <- as.factor(paste(starveday$treatment, starveday$starvedays, sep = "-"))

#correct time variable
starveday$T= starveday$time*60
starveday$T.factor= as.factor(starveday$T)

#initial plots
qplot(Bin,Cells, color = treatment, data = starveday,  geom = "boxplot") + facet_wrap(~starvedays) 
qplot(T.factor,Cells, color = treatment, data = starveday,  geom = "boxplot") + facet_wrap(starvedays~Bin, scales="free") 

#standardization per treatment and starveday

starveday$CellsS=NA
k=split(starveday, starveday$treatdays)
starvedaystd <- lapply(k, function (x) scale(x[,c("Cells")], center=T, scale=T))
starveday$CellsS=unsplit(starvedaystd, starveday$treatdays)

#baselining to 0 at time point 0
NT<-data.table(starveday, key=c("reptreatbindays"))

t1=NT[,list(treatment=treatment, starvedays=starvedays, Bin=Bin, reptreatdays=reptreatdays, treatdays=treatdays, T=T, Cells=Cells, CellsS=CellsS,
            CellsBase=(CellsS-CellsS[1])), by=c("reptreatbindays")]

starvedaybase <- t1 #DATA IS NOW CALLED starvedayBASE

starvedaybase$T.factor <- as.factor (starvedaybase$T)

qplot(T.factor, CellsBase, color = treatment, data = starvedaybase,  geom = "boxplot") + facet_wrap(starvedays~Bin) 

#summaries

starvedaybase.sum <- summarySE(starvedaybase, measurevar="CellsBase", groupvars=c("T", "Bin", "treatment", "starvedays"))

ggplot(data=starvedaybase.sum, aes(x=T, y=CellsBase, shape=treatment, color=treatment)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=CellsBase-se, ymax=CellsBase+se), width=10, size=1) + facet_grid(starvedays~Bin)

##BINS

#Bin A
BinA= subset (starvedaybase, Bin=='A')
expA=as.data.frame(data.table(cbind(treatment=BinA$treatment, T=BinA$T, starvedays=BinA$starvedays, ID=BinA$reptreatdays)))
cor(expA, method = "spearman")


vif_func(in_frame=expA,thresh=5,trace=T)

pairs(expA, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(BinA$CellsBase, group=BinA$reptreatdays, location="mean") #significant
levene.test(BinA$CellsBase, group=BinA$T, location="mean") # significant
levene.test(BinA$CellsBase, group=BinA$treatment, location="mean") #significant
levene.test(BinA$CellsBase, group=BinA$starveday, location="mean") #not significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(CellsBase~treatment, data=BinA)
boxplot(CellsBase~reptreatdays, data=BinA)
boxplot (CellsBase~T, data=BinA)
boxplot (CellsBase~starvedays, data=BinA)


#fit a gls
Form <- formula (CellsBase ~ treatment*T*starvedays)
BinA.gls<- gls(Form, data=BinA)


#####NLME models######


#random structures #winner: BinA1.lme
BinA1.lme <- lme (Form, random = ~1|reptreatdays, method="REML", data=BinA) 

BinA2.lme <- lme (Form, random = ~1|treatment/reptreatdays, method="REML", data=BinA)

anova(BinA.gls, BinA1.lme, BinA2.lme)


#random structures+correlation structures #all produce the same AIC values even BinA5.lme AIC=-80.53

BinA3.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (form=~1|reptreatdays/treatment), method="REML", data=BinA) 

BinA4.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), method="REML", data=BinA) 

BinA5.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (form=~1|reptreatdays/starvedays), method="REML", data=BinA) 

anova(BinA.gls, BinA1.lme, BinA2.lme, BinA3.lme, BinA4.lme, BinA5.lme)


#random structures+correlation structures+weights #BinA10.lme wins AIC=-100.29

#BinA6.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varIdent(form=~1|reptreatdays), method="REML", data=BinA) 

BinA7.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varIdent(form=~1|starvedays), method="REML", data=BinA) 

BinA8.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varIdent(form=~1|treatment), method="REML", data=BinA) 

#BinA9.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varComb(varIdent(form=~1|treatment), varIdent(form=~1|reptreatdays)), method="REML", data=BinA) 

BinA10.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varComb(varIdent(form=~1|treatment), varIdent(form=~1|starvedays)), 
                  method="REML", data=BinA)

anova(BinA.gls, BinA1.lme, BinA2.lme, BinA3.lme, BinA4.lme, BinA5.lme, BinA7.lme, BinA8.lme, BinA10.lme)

#if I only check the share of weight structures #BinA13.lme wins

BinA11.lme <- lme (Form, random = ~1|reptreatdays, weights=varIdent(form=~1|starvedays), method="REML", data=BinA) 

BinA12.lme <- lme (Form, random = ~1|reptreatdays,weights=varIdent(form=~1|treatment), method="REML", data=BinA) 

BinA13.lme <- lme (Form, random = ~1|reptreatdays, weights=varComb(varIdent(form=~1|treatment), varIdent(form=~1|starvedays)), 
                   method="REML", data=BinA)

anova(BinA.gls, BinA1.lme, BinA2.lme, BinA3.lme, BinA4.lme, BinA5.lme, BinA7.lme, BinA8.lme, BinA10.lme, BinA11.lme, BinA12.lme, BinA13.lme)

#OVERALL WINNER: Bin A10.lme

#residuals
BinA.E2<-resid(BinA10.lme,type="normalized")
BinA.F2<-fitted(BinA10.lme)
op<-par(mfrow=c(3,2),mar=c(5,5,4,3))
MyYlab="Residuals"

plot(x=BinA.F2,y=BinA.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(BinA.E2~treatment,data=BinA, main="Treatment",ylab=MyYlab)
boxplot(BinA.E2~starvedays, data=BinA, main="Length starvation",ylab=MyYlab,xlab="Length starvation (days)")
boxplot(BinA.E2~T, data=BinA, main="Time",ylab=MyYlab,xlab="Time(sec)")
hist(BinA.E2, main="Residuals", xlab=MyYlab)
plot(fitted(BinA10.lme), residuals(BinA10.lme),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(BinA10.lme), residuals(BinA10.lme)))
par(op)

summary(BinA10.lme)
anova(BinA10.lme)

#let's plot this!

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 
library (AICcmodavg)



#BinA fit
BinA.sum <- summarySE(BinA, measurevar="CellsBase", groupvars=c("T", "treatment", "starvedays"))

BinA.fit <- as.data.frame(predictSE.lme(BinA10.lme, BinA, se.fit = TRUE, level = 0,
                                        print.matrix = FALSE))

BinA.fit$upr <- BinA.fit$fit + (1.96 * BinA.fit$se)
BinA.fit$lwr <- BinA.fit$fit - (1.96 * BinA.fit$se)

BinA.fit.combdata <- cbind(BinA, BinA.fit)

ggplot(data=BinA.sum, aes(x=T, y=CellsBase, shape=treatment, color=treatment)) + geom_point(size=5)+
  geom_errorbar(aes(ymin=CellsBase-se, ymax=CellsBase+se), width=30, size=1) +
  geom_smooth(data=BinA.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, fill=treatment), method="lm", stat="identity", alpha=0.1)+ 
  scale_colour_manual(values = c(Control="lightcoral", dSi="steelblue2"), name="Treatment") + facet_wrap(~starvedays) +
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") + 
  labs(list(x = "Time (s)", y = "Normalized cell count", title="Bin A"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_blank(), 
        axis.title.x=element_blank(),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 

##Bin B

#Bin B
BinB= subset (starvedaybase, Bin=='B')
expA=as.data.frame(data.table(cbind(treatment=BinB$treatment, T=BinB$T, starvedays=BinB$starvedays, ID=BinB$reptreatdays)))
cor(expA, method = "spearman")


vif_func(in_frame=expA,thresh=5,trace=T)

pairs(expA, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(BinB$CellsBase, group=BinB$reptreatdays, location="mean") #significant
levene.test(BinB$CellsBase, group=BinB$T, location="mean") # significant
levene.test(BinB$CellsBase, group=BinB$treatment, location="mean") #significant
levene.test(BinB$CellsBase, group=BinB$starveday, location="mean") #significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(CellsBase~treatment, data=BinB)
boxplot(CellsBase~reptreatdays, data=BinB)
boxplot (CellsBase~T, data=BinB)
boxplot (CellsBase~starvedays, data=BinB)


#fit a gls
Form <- formula (CellsBase ~ treatment*T*starvedays)
BinB.gls<- gls(Form, data=BinB)


#####NLME models######


#random structures #winner: BinB1.lme
BinB1.lme <- lme (Form, random = ~1|reptreatdays, method="REML", data=BinB) 

BinB2.lme <- lme (Form, random = ~1|treatment/reptreatdays, method="REML", data=BinB)

anova(BinB.gls, BinB1.lme, BinB2.lme)


#random structures+correlation structures #all produce the same AIC values even BinB5.lme AIC=174.13

BinB3.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (form=~1|reptreatdays/treatment), method="REML", data=BinB) 

BinB4.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), method="REML", data=BinB) 

BinB5.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (form=~1|reptreatdays/starvedays), method="REML", data=BinB) 

anova(BinB.gls, BinB1.lme, BinB2.lme, BinB3.lme, BinB4.lme, BinB5.lme)


#random structures+correlation structures+weights #BinB10.lme wins AIC=134.38

BinB6.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varIdent(form=~1|reptreatdays), method="REML", data=BinB) 

BinB7.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varIdent(form=~1|starvedays), method="REML", data=BinB) 

BinB8.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varIdent(form=~1|treatment), method="REML", data=BinB) 

BinB9.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varComb(varIdent(form=~1|treatment), varIdent(form=~1|reptreatdays)), method="REML", data=BinB) 

BinB10.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varComb(varIdent(form=~1|treatment), varIdent(form=~1|starvedays)), 
                   method="REML", data=BinB)

anova(BinB.gls, BinB1.lme, BinB2.lme, BinB3.lme, BinB4.lme, BinB5.lme, BinB6.lme, BinB7.lme, BinB8.lme, BinB9.lme, BinB10.lme)

#if I only check the share of weight structures #BinB14.lme wins

BinB11.lme <- lme (Form, random = ~1|reptreatdays, weights=varIdent(form=~1|starvedays), method="REML", data=BinB) 

BinB12.lme <- lme (Form, random = ~1|reptreatdays, weights=varIdent(form=~1|treatment), method="REML", data=BinB) 

BinB13.lme <- lme (Form, random = ~1|reptreatdays, weights=varComb(varIdent(form=~1|treatment), varIdent(form=~1|starvedays)), 
                   method="REML", data=BinB)

BinB14.lme <- lme (Form, random = ~1|reptreatdays, weights=varIdent(form=~1|reptreatdays), method="REML", data=BinB) 

anova(BinB.gls, BinB1.lme, BinB2.lme, BinB3.lme, BinB4.lme, BinB5.lme, BinB6.lme, BinB7.lme, BinB8.lme, BinB9.lme, BinB10.lme,
      BinB11.lme, BinB12.lme, BinB13.lme, BinB14.lme)

#OVERALL WINNER: BinB6.lme but fit is ridiculous so a better fitting and a more logical model is BinB10.lme

#residuals
BinB.E2<-resid(BinB10.lme,type="normalized")
BinB.F2<-fitted(BinB10.lme)
op<-par(mfrow=c(3,2),mar=c(5,5,4,3))
MyYlab="Residuals"

plot(x=BinB.F2,y=BinB.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(BinB.E2~treatment,data=BinB, main="Treatment",ylab=MyYlab)
boxplot(BinB.E2~starvedays, data=BinB, main="Length starvation",ylab=MyYlab,xlab="Length starvation (days)")
boxplot(BinB.E2~T, data=BinB, main="Time",ylab=MyYlab,xlab="Time(sec)")
hist(BinB.E2, main="Residuals", xlab=MyYlab)
plot(fitted(BinB10.lme), residuals(BinB10.lme),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(BinB10.lme), residuals(BinB10.lme)))
par(op)

summary(BinB10.lme)
anova(BinB10.lme)

#BinB fit
BinB.sum <- summarySE(BinB, measurevar="CellsBase", groupvars=c("T", "treatment", "starvedays"))

BinB.fit <- as.data.frame(predictSE.lme(BinB10.lme, BinB, se.fit = TRUE, level = 0,
                                        print.matrix = FALSE))

BinB.fit$upr <- BinB.fit$fit + (1.96 * BinB.fit$se)
BinB.fit$lwr <- BinB.fit$fit - (1.96 * BinB.fit$se)

BinB.fit.combdata <- cbind(BinB, BinB.fit)

ggplot(data=BinB.sum, aes(x=T, y=CellsBase, shape=treatment, color=treatment)) + geom_point(size=5)+
  geom_errorbar(aes(ymin=CellsBase-se, ymax=CellsBase+se), width=30, size=1) +
  geom_smooth(data=BinB.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, fill=treatment), method="lm", stat="identity", alpha=0.1)+ 
  scale_colour_manual(values = c(Control="lightcoral", dSi="steelblue2"), name="Treatment") + facet_wrap(~starvedays) +
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") + 
  labs(list(x = "Time (s)", y = "Normalized cell count", title="Bin B"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_blank(), 
        axis.title.x=element_blank(),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 


#Bin C

#Bin B
BinC= subset (starvedaybase, Bin=='C')
expA=as.data.frame(data.table(cbind(treatment=BinC$treatment, T=BinC$T, starvedays=BinC$starvedays, ID=BinC$reptreatdays)))
cor(expA, method = "spearman")


vif_func(in_frame=expA,thresh=5,trace=T)

pairs(expA, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(BinC$CellsBase, group=BinC$reptreatdays, location="mean") #significant
levene.test(BinC$CellsBase, group=BinC$T, location="mean") # significant
levene.test(BinC$CellsBase, group=BinC$treatment, location="mean") #significant
levene.test(BinC$CellsBase, group=BinC$starveday, location="mean") #significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(CellsBase~treatment, data=BinC)
boxplot(CellsBase~reptreatdays, data=BinC)
boxplot (CellsBase~T, data=BinC)
boxplot (CellsBase~starvedays, data=BinC)


#fit a gls
Form <- formula (CellsBase ~ treatment*T*starvedays)
BinC.gls<- gls(Form, data=BinC)


#####NLME models######


#random structures #winner: BinC1.lme
BinC1.lme <- lme (Form, random = ~1|reptreatdays, method="REML", data=BinC) 

BinC2.lme <- lme (Form, random = ~1|treatment/reptreatdays, method="REML", data=BinC)

anova(BinC.gls, BinC1.lme, BinC2.lme)


#random structures+correlation structures #all produce the same AIC values even BinC5.lme AIC=160.64

BinC3.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (form=~1|reptreatdays/treatment), method="REML", data=BinC) 

BinC4.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), method="REML", data=BinC) 

BinC5.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (form=~1|reptreatdays/starvedays), method="REML", data=BinC) 

anova(BinC.gls, BinC1.lme, BinC2.lme, BinC3.lme, BinC4.lme, BinC5.lme)


#random structures+correlation structures+weights #BinC10.lme wins AIC=142.58

BinC6.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varIdent(form=~1|reptreatdays), method="REML", data=BinC) 

BinC7.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varIdent(form=~1|starvedays), method="REML", data=BinC) 

BinC8.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varIdent(form=~1|treatment), method="REML", data=BinC) 

BinC9.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varComb(varIdent(form=~1|treatment), varIdent(form=~1|reptreatdays)), method="REML", data=BinC) 

BinC10.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varComb(varIdent(form=~1|treatment), varIdent(form=~1|starvedays)), 
                   method="REML", data=BinC)

anova(BinC.gls, BinC1.lme, BinC2.lme, BinC3.lme, BinC4.lme, BinC5.lme, BinC6.lme, BinC7.lme, BinC8.lme, BinC9.lme, BinC10.lme)

#if I only check the share of weight structures #BinC14.lme wins

BinC11.lme <- lme (Form, random = ~1|reptreatdays, weights=varIdent(form=~1|starvedays), method="REML", data=BinC) 

BinC12.lme <- lme (Form, random = ~1|reptreatdays, weights=varIdent(form=~1|treatment), method="REML", data=BinC) 

BinC13.lme <- lme (Form, random = ~1|reptreatdays, weights=varComb(varIdent(form=~1|treatment), varIdent(form=~1|starvedays)), 
                   method="REML", data=BinC)

BinC14.lme <- lme (Form, random = ~1|reptreatdays, weights=varIdent(form=~1|reptreatdays), method="REML", data=BinC) 

anova(BinC.gls, BinC1.lme, BinC2.lme, BinC3.lme, BinC4.lme, BinC5.lme, BinC6.lme, BinC7.lme, BinC8.lme, BinC9.lme, BinC10.lme,
      BinC11.lme, BinC12.lme, BinC13.lme, BinC14.lme)

#OVERALL WINNER: BinC6.lme but fit is ridiculous so a better fitting and a more logical model is BinC10.lme

#residuals
BinC.E2<-resid(BinC10.lme,type="normalized")
BinC.F2<-fitted(BinC10.lme)
op<-par(mfrow=c(3,2),mar=c(5,5,4,3))
MyYlab="Residuals"

plot(x=BinC.F2,y=BinC.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(BinC.E2~treatment,data=BinC, main="Treatment",ylab=MyYlab)
boxplot(BinC.E2~starvedays, data=BinC, main="Length starvation",ylab=MyYlab,xlab="Length starvation (days)")
boxplot(BinC.E2~T, data=BinC, main="Time",ylab=MyYlab,xlab="Time(sec)")
hist(BinC.E2, main="Residuals", xlab=MyYlab)
plot(fitted(BinC10.lme), residuals(BinC10.lme),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(BinC10.lme), residuals(BinC10.lme)))
par(op)

summary(BinC10.lme)
anova(BinC10.lme)

#BinC fit
BinC.sum <- summarySE(BinC, measurevar="CellsBase", groupvars=c("T", "treatment", "starvedays"))

BinC.fit <- as.data.frame(predictSE.lme(BinC10.lme, BinC, se.fit = TRUE, level = 0,
                                        print.matrix = FALSE))

BinC.fit$upr <- BinC.fit$fit + (1.96 * BinC.fit$se)
BinC.fit$lwr <- BinC.fit$fit - (1.96 * BinC.fit$se)

BinC.fit.combdata <- cbind(BinC, BinC.fit)

ggplot(data=BinC.sum, aes(x=T, y=CellsBase, shape=treatment, color=treatment)) + geom_point(size=5)+
  geom_errorbar(aes(ymin=CellsBase-se, ymax=CellsBase+se), width=30, size=1) +
  geom_smooth(data=BinC.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, fill=treatment), method="lm", stat="identity", alpha=0.1)+ 
  scale_colour_manual(values = c(Control="lightcoral", dSi="steelblue2"), name="Treatment") + facet_wrap(~starvedays) +
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") + 
  labs(list(x = "Time (s)", y = "Normalized cell count", title="Bin C"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_blank(), 
        axis.title.x=element_blank(),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 


###FITTING ALL DATA

allbins.fitdata = rbind (BinA.fit.combdata, BinB.fit.combdata, BinC.fit.combdata)
BinA.sum$Bin="A"
BinB.sum$Bin="B"
BinC.sum$Bin="C"

allbins.sum = rbind (BinA.sum, BinB.sum, BinC.sum)

mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="Bin") { 
    value[value=="A"] <- "Bin A"
    value[value=="B"]   <- "Bin B"
    value[value=="C"] <- "Bin C"
  }
  return(value)
}

ggplot(data=allbins.sum, aes(x=T, y=CellsBase, shape=treatment, color=treatment)) + geom_point(size=5)+ 
  facet_grid(Bin~starvedays, labeller=mf_labeller, scales="free")+
  geom_errorbar(aes(ymin=CellsBase-se, ymax=CellsBase+se), width=20, size=1) +
  geom_smooth(data=allbins.fitdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, fill=treatment), method="lm", stat="identity", alpha=0.1)+ 
  scale_colour_manual(values = c(Control="lightcoral", dSi="steelblue2"), name="Treatment") +
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") + 
  labs(list(x = "Time (s)", y = "Normalized cell count"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (1, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 


