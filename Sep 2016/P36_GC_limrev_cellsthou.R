library(ggplot2)
library(gtable)
library(ggthemes)
library(grLabelExtra)
library(mgcv)
library(ggplot2)
library(grid)
library(data.table)

source("AED.R")
source("vif.R")
source("summarySE.R")
source("resizewin.R")
resize.win(12,9)


#read in data

P36_GC_limrev <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/growth curves/P36_limitation and recovery_GC_longformat.csv", sep=";")

#make new factors

P36_GC_limrev$medtreat <- as.factor (paste(P36_GC_limrev$med, P36_GC_limrev$treatment, sep="-"))
P36_GC_limrev$medtreatreplet <- as.factor (paste(P36_GC_limrev$medtreat, P36_GC_limrev$replet, sep="-"))
P36_GC_limrev$daynum <- (as.numeric(P36_GC_limrev$day))-1
P36_GC_limrev$cellsmlthou <- (P36_GC_limrev$cellsml)/1000

#plotting
qplot(as.factor(daynum), cellsmlthou, data = P36_GC_limrev,  geom = "boxplot") + facet_grid(medtreat~., scales="free") 

source("summarySE.R")
GC.sum <- summarySE(P36_GC_limrev, measurevar="cellsmlthou", groupvars=c("medtreat", "daynum"))

#plot with log normal cells R
ggplot(data=GC.sum, aes(x=daynum, y=cellsmlthou)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsmlthou-se, ymax=cellsmlthou+se), width=0.5, size=1) + facet_wrap(~medtreat, nrow=2)+
  geom_smooth(method="loess", aes(group=1))

ggplot(data=GC.sum, aes(x=daynum, y=cellsmlthou)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsmlthou-se, ymax=cellsmlthou+se), width=0.5, size=1) + facet_wrap(~medtreat, nrow=2)+
  geom_smooth(method="lm", aes(group=1))

#P36_GC_limrev

source ("AED.R")
source ("vif.R")
library (lawstat)

expGC=as.data.frame(data.table(cbind(medtreat=P36_GC_limrev$medtreat, daynum=P36_GC_limrev$daynum, ID=P36_GC_limrev$medtreatreplet)))
cor(expGC, method = "spearman")


vif_func(in_frame=expGC,thresh=5,trace=T)

pairs(expGC, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(P36_GC_limrev$cellsmlthou, group=P36_GC_limrev$medtreatreplet, location="mean") #significant
levene.test(P36_GC_limrev$cellsmlthou, group=P36_GC_limrev$daynum, location="mean") # not significant
levene.test(P36_GC_limrev$cellsmlthou, group=P36_GC_limrev$medtreat, location="mean") # significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(cellsmlthou~medtreat, data=P36_GC_limrev)
boxplot(cellsmlthou~medtreatreplet, data=P36_GC_limrev)
boxplot (cellsmlthou~daynum, data=P36_GC_limrev)


#fit a gls
Form <- formula (cellsmlthou ~ medtreat*daynum)
P36_GC_limrev.gls<- gls(Form, data=P36_GC_limrev)


#nlme model

#random factor
P36_GC_limrev1.lme <- lme (Form, random = ~1|medtreatreplet, method="REML", data=P36_GC_limrev) #BEST

P36_GC_limrev2.lme <- lme (Form, random = ~1|medtreat, method="REML", data=P36_GC_limrev)

P36_GC_limrev3.lme <- lme (Form, random = ~1|medtreatreplet/medtreat, method="REML", data=P36_GC_limrev)

anova(P36_GC_limrev.gls, P36_GC_limrev1.lme, P36_GC_limrev2.lme, P36_GC_limrev3.lme)

#variance structure

#P36_GC_limrev4.lme <- lme (Form, random = ~1|medtreatreplet,  weights=varIdent(form=~1|medtreatreplet), method="REML", data=P36_GC_limrev)

P36_GC_limrev5.lme <- lme (Form, random = ~1|medtreatreplet,  weights=varIdent(form=~1|medtreat), method="REML", data=P36_GC_limrev) #BEST but with no variance structure is also okay

#P36_GC_limrev6.lme <- lme (Form, random = ~1|medtreatreplet,  weights=varIdent(form=~1|medtreatreplet/medtreat), method="REML", data=P36_GC_limrev) 

anova(P36_GC_limrev.gls, P36_GC_limrev1.lme, P36_GC_limrev2.lme, P36_GC_limrev3.lme, P36_GC_limrev5.lme)

#correlation structures

P36_GC_limrev7.lme <- lme (Form, random = ~1|medtreatreplet, correlation=corAR1(), method="REML", data=P36_GC_limrev) #BEST use this

P36_GC_limrev8.lme <- lme (Form, random = ~1|medtreatreplet,  
                           correlation=corAR1 (form=~1|medtreatreplet), method="REML", data=P36_GC_limrev) 

P36_GC_limrev9.lme <- lme (Form, random = ~1|medtreatreplet,  
                           correlation=corAR1 (form=~1|medtreatreplet/medtreat), method="REML", data=P36_GC_limrev) 


anova(P36_GC_limrev.gls, P36_GC_limrev1.lme, P36_GC_limrev2.lme, P36_GC_limrev3.lme, P36_GC_limrev5.lme, P36_GC_limrev7.lme, P36_GC_limrev8.lme, P36_GC_limrev9.lme)

#models 7-9 are the same are the same because it treats the correlation structures the same! HAHA

summary(P36_GC_limrev9.lme)
anova(P36_GC_limrev9.lme)

#multiple comparisons
library (multcomp)

summary(glht(P36_GC_limrev9.lme, linfct=mcp(medtreat="Tukey", covariate_average=TRUE, interaction_average = TRUE)))

library(lsmeans)

pairs(lsmeans(P36_GC_limrev9.lme, ~medtreat|daynum))

library(multcompView)
cld(lsmeans(P36_GC_limrev9.lme, ~medtreat), alpha=0.05)

#residuals
P36_GC_limrev.E2<-resid(P36_GC_limrev1.lme,type="normalized")
P36_GC_limrev.F2<-fitted(P36_GC_limrev1.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=P36_GC_limrev.F2,y=P36_GC_limrev.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(P36_GC_limrev.E2~medtreat,data=P36_GC_limrev, main="medtreat",ylab=MyYlab)
plot(x=P36_GC_limrev$daynum,y=P36_GC_limrev.E2,main="Time",ylab=MyYlab,xlab="daynum")
par(op)

xyplot (P36_GC_limrev.E2 ~ daynum| medtreat, data=P36_GC_limrev, ylab="Residuals", xlab="daynum", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})


#fitting data

library(AICcmodavg)

#DPR fit
P36_GC_limrev.fit <- as.data.frame(predictSE.lme(P36_GC_limrev9.lme, P36_GC_limrev, se.fit = TRUE, level = 0,
                                                 print.matrix = FALSE))

P36_GC_limrev.fit$upr <- P36_GC_limrev.fit$fit + (1.96 * P36_GC_limrev.fit$se)
P36_GC_limrev.fit$lwr <- P36_GC_limrev.fit$fit - (1.96 * P36_GC_limrev.fit$se)

P36_GC_limrev.fit.combdata <- cbind(P36_GC_limrev, P36_GC_limrev.fit)


#plot

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 

ggplot(data=GC.sum, aes(x=daynum, y=cellsmlthou)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsmlthou-se, ymax=cellsmlthou+se), width=0.5, size=1) +  facet_wrap(~medtreat, nrow=2)+
  geom_smooth(data=P36_GC_limrev.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, group=1), color="black", method="lm", stat="identity", alpha=0.2)+ 
  labs(list(x = "daynum", y = "cells thou per cm2"))



