library(ggplot2)
library(gtable)
library(ggthemes)
library(mgcv)
library(ggplot2)
library(grid)
library(data.table)
library(lattice)
library(plyr)
s=mgcv:::s

source("AED.R")
source("vif.R")
source("tsDiagGamm.R")
source("summarySE.R")
source("resizewin.R")
source("lang.R")


##read data

trackdata <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/choice/final track data.csv", sep=";")

#make another copy of the trackdata just in case

trackdata2 <- trackdata

#make time as min

trackdata2$timemin <- trackdata2$time/60

trackdata2$timeminfac <- as.factor(trackdata2$time/60)


#drop bin_both

trackdata2 <- trackdata2 [! trackdata2$bin=="bin_both",  ]

trackdata2$bin <- droplevels(trackdata2$bin)

#calculate distance between cell position and bead
trackdata2$dist_dSi=sqrt(((trackdata2$x_dSi-trackdata2$X)^2)+((trackdata2$y_dSi-trackdata2$Y)^2))
trackdata2$dist_DPR=sqrt(((trackdata2$x_DPR-trackdata2$X)^2)+((trackdata2$y_DPR-trackdata2$Y)^2))


#calculate the sumdist for all. REMEMBER THIS IS YOUR DATA FOR THIS ANALYSIS. NOT SUMMARRY SE RESULTS BUT SUMDIST!!!
#dSi
dSi_dist <- ddply(trackdata2, c("wellvid", "treatment", "ID", "timeminfac", "timemin"), summarise,
                  N    = length(ID),
                  mean = mean(dist_dSi, na.rm=TRUE),
                  sumdist= sum(dist_dSi, na.rm=TRUE), 
                  sd   = sd(dist_dSi, na.rm=TRUE),
                  se   = sd / sqrt(N))

dSi_dist$distcm = dSi_dist$sumdist/10000

dSi_dist$bead ="dSibead"

dSi_dist.sum <- summarySE(dSi_dist, measurevar="mean", groupvars=c("treatment", "bead", "timemin", "timeminfac"))

qplot(timeminfac, mean, color = treatment, data = dSi_dist.sum)+ geom_point(size=5) + facet_grid(~treatment, scales="free") + geom_smooth(aes(group=1))

#DPR
DPR_dist <- ddply(trackdata2, c("wellvid", "treatment",  "ID", "timeminfac", "timemin"), summarise,
                  N    = length(ID),
                  mean = mean(dist_DPR, na.rm=TRUE),
                  sumdist= sum(dist_DPR, na.rm=TRUE), 
                  sd   = sd(dist_DPR, na.rm=TRUE),
                  se   = sd / sqrt(N))

DPR_dist$distcm = DPR_dist$sumdist/10000

DPR_dist$bead ="DPRbead"

DPR_dist.sum <- summarySE(DPR_dist, measurevar="mean", groupvars=c("treatment", "bead", "timemin", "timeminfac"))

qplot(timeminfac, mean, color = treatment, data = DPR_dist.sum)+ geom_point(size=5) + facet_grid(~treatment, scales="free") + geom_smooth(aes(group=1))

#binding

trackdata.dist = rbind(dSi_dist, DPR_dist)

trackdata.dist$meandist <- trackdata.dist$mean

trackdata.dist$bead <- as.factor(trackdata.dist$bead)

trackdata.dist.sum = rbind (dSi_dist.sum, DPR_dist.sum)

qplot(timemin, mean, color = treatment, data = trackdata.dist.sum)+ geom_point(size=5) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5, size=1) +
  facet_grid(treatment~bead, scales="free") + stat_smooth(aes(outfit=fit<<-..y..)) # fitted data was extracted as fit


##subsetting
induced.dist = subset (trackdata.dist, trackdata.dist$treatment=="Si_induced")

notinduced.dist = subset (trackdata.dist, trackdata.dist$treatment=="Si_notinduced")


#LME: induced.dist

qplot(timeminfac, meandist, color = bead, data = induced.dist,  geom = "boxplot") + facet_grid(bead~., scales="free") +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=1))

exp.induced=as.data.frame(data.table(cbind(wellvid=as.numeric(as.factor(induced.dist$wellvid)),bead=as.numeric(as.factor(induced.dist$bead)), 
                                           T=as.numeric(induced.dist$timeminfac), ID=as.numeric(as.factor(induced.dist$ID)))))

cor(exp.induced, method = "spearman")

vif_func(in_frame=exp.induced,thresh=5,trace=T)

pairs(exp.induced, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

#boxplots
op=par(mfrow=c(2,2))
boxplot(meandist~wellvid, data=induced.dist)
boxplot(meandist~bead, data=induced.dist)
boxplot(meandist~ID, data=induced.dist)
boxplot (meandist~timeminfac, data=induced.dist)

#levene
library(lawstat)
levene.test(induced.dist$meandist, group=induced.dist$wellvid, location="mean") #significant
levene.test(induced.dist$meandist, group=induced.dist$ID, location="mean") # significant
levene.test(induced.dist$meandist, group=induced.dist$timeminfac, location="mean")  # not significant
levene.test(induced.dist$meandist, group=induced.dist$bead, location="mean") # significant

#fit a gls
Form.induced <- formula (meandist ~ bead*timemin)
induced.gls<- gls(Form.induced, data=induced.dist)

#nlme model

#random factor
induced1.lme <- lme (Form.induced, random = ~1|ID, method="REML", data=induced.dist) 

induced2.lme <- lme (Form.induced, random = ~1|bead, method="REML", data=induced.dist)

induced2a.lme <- lme (Form.induced, random=list(ID=~1), method="REML", data=induced.dist) #second best

induced3.lme <- lme (Form.induced, random = ~1|ID/bead, method="REML", data=induced.dist) #best

induced3a.lme <- lme (Form.induced, random=list(ID=~1|bead), method="REML", data=induced.dist)

induced4.lme <- lme (Form.induced, random = ~1|wellvid, method="REML", data=induced.dist)

anova(induced.gls, induced1.lme, induced2.lme, induced2a.lme, induced3.lme, induced3a.lme, induced4.lme) #2a and 3a the same


#variance structure

#induced5.lme <- lme (Form.induced, random = ~1|ID/bead,  weights=varIdent(form=~1|ID), method="REML", data=induced.dist)

induced6.lme <- lme (Form.induced, random = ~1|ID/bead,  weights=varIdent(form=~1|bead), method="REML", data=induced.dist) #best

induced6a.lme <- lme (Form.induced,  random=list(ID=~1),  weights=varIdent(form=~1|bead), method="REML", data=induced.dist) 

#induced7.lme <- lme (Form.induced, random = ~1|ID/bead,  weights=varIdent(form=~1|ID/bead), method="REML", data=induced.dist) 

#induced8.lme <- lme (Form.induced, random = ~1|ID/bead,  weights=varIdent(form=~1|wellvid), method="REML", data=induced.dist) 

induced8a.lme <- lme (Form.induced, random=list(ID=~1),  weights=varIdent(form=~1|wellvid), method="REML", data=induced.dist) 

anova(induced.gls, induced1.lme, induced2.lme, induced2a.lme, induced3.lme, induced3a.lme, induced4.lme,
      induced6.lme, induced6a.lme, induced8a.lme)

#6 and 6a, 8 and 8a produce the same result. but 6 is lower

#correlation structures

ctrl <- lmeControl (maxIter = 100, msMaxIter = 100, msMaxEval = 400)

induced9.lme <- lme (Form.induced,random = ~1|ID/bead,  weights=varIdent(form=~1|bead), control = ctrl, 
                     correlation=corAR1(), method="REML", data=induced.dist) #BEST use this

#induced10.lme <- lme (Form.induced, random = ~1|ID/bead, weights=varIdent(form=~1|bead), control = ctrl, correlation=corAR1 (form=~1|bead), method="REML", data=induced.dist) 

#induced11.lme <- lme (Form.induced, random = ~1|ID,  weights=varIdent(form=~1|bead), control = ctrl, correlation=corAR1 (form=~1|bead), method="REML", data=induced.dist) 


anova(induced.gls, induced1.lme, induced2.lme, induced2a.lme, induced3.lme, induced3a.lme, induced4.lme,
      induced6.lme, induced6a.lme, induced8a.lme, induced9.lme)

#best induced9.lme, AIC=188984.8 

summary(induced9.lme)
anova(induced9.lme)

#multiple comparisons
#library(lsmeans)

#pairs(lsmeans(induced11.lme, ~bead|timeminfac))

#library(multcompView)
#cld(lsmeans(induced11.lme, ~bead), alpha=0.05)

library(multcomp)

#glht doesnt work
summary(glht(induced9.lme, linfct=mcp(bead="Tukey", covariate_average=TRUE, interaction_average = TRUE)))

#residuals
induced.E2<-resid(induced9.lme,type="normalized")
induced.F2<-fitted(induced9.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=induced.F2,y=induced.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(induced.E2~bead,data=induced.dist, main="bead",ylab=MyYlab)
plot(x=induced.dist$timemin,y=induced.E2,main="Time",ylab=MyYlab,xlab="Time (min)")
boxplot(induced.E2~wellvid,data=induced.dist, main="wellvid",ylab=MyYlab)
par(op)

xyplot (induced.E2 ~ timemin| bead, data=induced.dist, ylab="Residuals", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})

#LME: notinduced.dist

qplot(timeminfac, meandist, color = bead, data = notinduced.dist,  geom = "boxplot") + facet_grid(bead~., scales="free") +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=1))

exp.notinduced=as.data.frame(data.table(cbind(wellvid=as.numeric(as.factor(notinduced.dist$wellvid)),bead=as.numeric(as.factor(notinduced.dist$bead)), 
                                              T=as.numeric(notinduced.dist$timeminfac), ID=as.numeric(as.factor(notinduced.dist$ID)))))

cor(exp.notinduced, method = "spearman")

vif_func(in_frame=exp.notinduced,thresh=5,trace=T)

pairs(exp.notinduced, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

#boxplots
op=par(mfrow=c(2,2))
boxplot(meandist~wellvid, data=notinduced.dist)
boxplot(meandist~bead, data=notinduced.dist)
boxplot(meandist~ID, data=notinduced.dist)
boxplot (meandist~timeminfac, data=notinduced.dist)

#levene
levene.test(notinduced.dist$meandist, group=notinduced.dist$wellvid, location="mean") #significant
levene.test(notinduced.dist$meandist, group=notinduced.dist$ID, location="mean") # significant
levene.test(notinduced.dist$meandist, group=notinduced.dist$timeminfac, location="mean")  # not significant
levene.test(notinduced.dist$meandist, group=notinduced.dist$bead, location="mean") # significant

#fit a gls
Form.notinduced <- formula (meandist ~ bead*timemin)
notinduced.gls<- gls(Form.notinduced, data=notinduced.dist)

#nlme model

#random factor
notinduced1.lme <- lme (Form.notinduced, random = ~1|ID, method="REML", data=notinduced.dist) 

notinduced2.lme <- lme (Form.notinduced, random = ~1|bead, method="REML", data=notinduced.dist)

notinduced2a.lme <- lme (Form.notinduced, random=list(ID=~1), method="REML", data=notinduced.dist) #second best

notinduced3.lme <- lme (Form.notinduced, random = ~1|ID/bead, method="REML", data=notinduced.dist) #best

notinduced3a.lme <- lme (Form.notinduced, random=list(ID=~1|bead), method="REML", data=notinduced.dist)

notinduced4.lme <- lme (Form.notinduced, random = ~1|wellvid, method="REML", data=notinduced.dist)

anova(notinduced.gls, notinduced1.lme, notinduced2.lme, notinduced2a.lme, notinduced3.lme, notinduced3a.lme, notinduced4.lme) #2a and 3a the same


#variance structure

#notinduced5.lme <- lme (Form.notinduced, random = ~1|ID/bead,  weights=varIdent(form=~1|ID), method="REML", data=notinduced.dist)

notinduced6.lme <- lme (Form.notinduced, random = ~1|ID/bead,  weights=varIdent(form=~1|bead), method="REML", data=notinduced.dist) #best

notinduced6a.lme <- lme (Form.notinduced,  random=list(ID=~1),  weights=varIdent(form=~1|bead), method="REML", data=notinduced.dist) 

#notinduced7.lme <- lme (Form.notinduced, random = ~1|ID/bead,  weights=varIdent(form=~1|ID/bead), method="REML", data=notinduced.dist) 

#notinduced8.lme <- lme (Form.notinduced, random = ~1|ID/bead,  weights=varIdent(form=~1|wellvid), method="REML", data=notinduced.dist) 

notinduced8a.lme <- lme (Form.notinduced, random=list(ID=~1),  weights=varIdent(form=~1|wellvid), method="REML", data=notinduced.dist) 

anova(notinduced.gls, notinduced1.lme, notinduced2.lme, notinduced2a.lme, notinduced3.lme, notinduced3a.lme, notinduced4.lme,
      notinduced6.lme, notinduced6a.lme, notinduced8a.lme)

#6 and 6a, 8 and 8a produce the same result. but 6 is lower

#correlation structures

ctrl <- lmeControl (maxIter = 100, msMaxIter = 100, msMaxEval = 400)

notinduced9.lme <- lme (Form.notinduced,random = ~1|ID/bead,  weights=varIdent(form=~1|bead), control = ctrl, 
                        correlation=corAR1(), method="REML", data=notinduced.dist) #BEST use this

notinduced10.lme <- lme (Form.notinduced, random = ~1|ID/bead, weights=varIdent(form=~1|bead), control = ctrl, 
                         correlation=corAR1 (form=~1|bead), method="REML", data=notinduced.dist) 

notinduced11.lme <- lme (Form.notinduced, random = ~1|ID,  weights=varIdent(form=~1|bead), control = ctrl, 
                         correlation=corAR1 (form=~1|bead), method="REML", data=notinduced.dist) 

anova(notinduced.gls, notinduced1.lme, notinduced2.lme, notinduced2a.lme, notinduced3.lme, notinduced3a.lme, notinduced4.lme,
      notinduced6.lme, notinduced6a.lme, notinduced8a.lme, notinduced9.lme)

#best notinduced9.lme, AIC=194217.8

summary(notinduced9.lme)
anova(notinduced9.lme)

#multiple comparisons
#library(lsmeans)

#pairs(lsmeans(notinduced11.lme, ~bead|timeminfac))

#library(multcompView)
#cld(lsmeans(notinduced11.lme, ~bead), alpha=0.05)

library(multcomp)

#glht doesnt work
summary(glht(notinduced9.lme, linfct=mcp(bead="Tukey", covariate_average=TRUE, interaction_average = TRUE)))

#residuals
notinduced.E2<-resid(notinduced9.lme,type="normalized")
notinduced.F2<-fitted(notinduced9.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=notinduced.F2,y=notinduced.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(notinduced.E2~bead,data=notinduced.dist, main="bead",ylab=MyYlab)
plot(x=notinduced.dist$timemin,y=notinduced.E2,main="Time",ylab=MyYlab,xlab="Time (min)")
boxplot(notinduced.E2~wellvid,data=notinduced.dist, main="wellvid",ylab=MyYlab)
par(op)

xyplot (notinduced.E2 ~ timemin| bead, data=notinduced.dist, ylab="Residuals", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})

#plot

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 

scaleFUN <- function(x) sprintf("%.1f", x)

trackdata.dist.sum$bead2 <- factor(trackdata.dist.sum$bead, levels=c("DPRbead", "dSibead"), labels =c ("Diproline bead", "dSi bead"))
trackdata.dist.sum$treatment2 <-  factor(trackdata.dist.sum$treatment, levels=c("Si_induced", "Si_notinduced"), labels =c ("induced", "not induced"))

resize.win (12,9)
ggplot(data=trackdata.dist.sum, aes(x=timemin, y=mean)) + geom_point(size=5)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5, size=1) +
  facet_grid(treatment2~bead2, scale="free") +
  labs(list(x = "Time (min)", y = "Mean distance relative to the bead (µm/s)"))+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5),
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text = element_text(size=15), legend.title=text, legend.text=text, panel.margin=unit (1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_x_continuous (breaks=c(0, 2, 4, 6, 8, 10))+
  scale_y_continuous(labels=scaleFUN)

