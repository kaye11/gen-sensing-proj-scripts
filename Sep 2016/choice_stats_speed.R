library(ggplot2)
library(gtable)
library(ggthemes)
library(mgcv)
library(ggplot2)
library(grid)
library(data.table)
library(lattice)
s=mgcv:::s

source("AED.R")
source("vif.R")
source("tsDiagGamm.R")
source("summarySE.R")
source("resizewin.R")
source("lang.R")


##read data

trackdata <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/choice/final track data.csv", sep=";")

#induced.trackdata <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/choice/induced track data.csv", sep=";")

#make another copy of the trackdata just in case

trackdata2 <- trackdata

#make time as min

trackdata2$timemin <- trackdata2$time/60

trackdata2$timeminfac <- as.factor(trackdata2$time/60)


#drop bin_both

trackdata2 <- trackdata2 [! trackdata2$bin=="bin_both",  ]

trackdata2$bin <- droplevels(trackdata2$bin)

#make summaries

trackdata.sum <- summarySE(trackdata2, measurevar="Vlog", groupvars=c("wellvid",  "ID", "treatment", "bin","timemin", "timeminfac"), na.rm=TRUE)

trackdata.sumV <- summarySE(trackdata2, measurevar="V", groupvars=c("wellvid",  "ID", "treatment", "bin","timemin", "timeminfac"), na.rm=TRUE)


trackdata.sum2 <- ddply(trackdata2, c("wellvid",  "ID", "treatment", "bin","timemin", "timeminfac"), summarise,
                  N    = sum(!is.na(Vlog)),
                  mean = mean(Vlog, na.rm=TRUE),
                  sd   = sd(Vlog, na.rm=TRUE),
                  se   = sd / sqrt(N))

#same results of sum and sum2

trackdata.sum.noID <- summarySE(trackdata2, measurevar="Vlog", groupvars=c("treatment", "bin", "timemin", "timeminfac"))

trackdata.sum.noID.V <- summarySE(trackdata2, measurevar="V", groupvars=c("treatment", "bin", "timemin", "timeminfac"))

# plot

ggplot(data=trackdata.sum.noID, aes(x=timemin, y=Vlog, shape=bin, color=bin)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=Vlog-se, ymax=Vlog+se), width=0.5, size=1) + facet_grid(bin~treatment)

ggplot(data=trackdata.sum.noID.V, aes(x=timemin, y=V, shape=bin, color=bin)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=V-se, ymax=V+se), width=0.5, size=1) + facet_grid(bin~treatment)

#stats on summary because its impossible to run the whole dataset :( R crashes if you try

induced.sum <- subset(trackdata.sum, trackdata.sum$treatment=="Si_induced")

notinduced.sum <- subset(trackdata.sum, trackdata.sum$treatment=="Si_notinduced")



#LME: induced.sum

qplot(timeminfac, Vlog, color = bin, data = induced.sum,  geom = "boxplot") + facet_grid(bin~., scales="free") +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=1))

exp.induced=as.data.frame(data.table(cbind(wellvid=as.numeric(as.factor(induced.sum$wellvid)),bin=as.numeric(as.factor(induced.sum$bin)), 
                                           T=as.numeric(induced.sum$timeminfac), ID=as.numeric(as.factor(induced.sum$ID)))))

cor(exp.induced, method = "spearman")

vif_func(in_frame=exp.induced,thresh=5,trace=T)

pairs(exp.induced, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

#boxplots
op=par(mfrow=c(2,2))
boxplot(Vlog~wellvid, data=induced.sum)
boxplot(Vlog~bin, data=induced.sum)
boxplot(Vlog~ID, data=induced.sum)
boxplot (Vlog~timemin, data=induced.sum)

#levene
library(lawstat)
levene.test(induced.sum$Vlog, group=induced.sum$wellvid, location="mean") #significant
levene.test(induced.sum$Vlog, group=induced.sum$ID, location="mean") # significant
levene.test(induced.sum$Vlog, group=induced.sum$timeminfac, location="mean")  # significant
levene.test(induced.sum$Vlog, group=induced.sum$bin, location="mean") # significant

#fit a gls
Form.induced <- formula (Vlog ~ bin*timeminfac)
induced.gls<- gls(Form.induced, data=induced.sum)

#nlme model

#random factor
induced1.lme <- lme (Form.induced, random = ~1|ID, method="REML", data=induced.sum) #BEST

induced2.lme <- lme (Form.induced, random = ~1|bin, method="REML", data=induced.sum)

induced2a.lme <- lme (Form.induced, random=list(ID=~1), method="REML", data=induced.sum) #second best

induced3.lme <- lme (Form.induced, random = ~1|ID/bin, method="REML", data=induced.sum) #best

induced3a.lme <- lme (Form.induced, random=list(ID=~1|bin), method="REML", data=induced.sum)

induced4.lme <- lme (Form.induced, random = ~1|wellvid, method="REML", data=induced.sum)

anova(induced.gls, induced1.lme, induced2.lme, induced2a.lme, induced3.lme, induced3a.lme, induced4.lme) #2a and 3a the same


#variance structure

#induced5.lme <- lme (Form.induced, random = ~1|ID,  weights=varIdent(form=~1|ID), method="REML", data=induced.sum)

induced6.lme <- lme (Form.induced, random = ~1|ID,  weights=varIdent(form=~1|bin), method="REML", data=induced.sum) 

induced6a.lme <- lme (Form.induced,  random=list(ID=~1),  weights=varIdent(form=~1|bin), method="REML", data=induced.sum) 

#induced7.lme <- lme (Form.induced, random = ~1|ID,  weights=varIdent(form=~1|ID/bin), method="REML", data=induced.sum) 

induced8.lme <- lme (Form.induced, random = ~1|ID,  weights=varIdent(form=~1|wellvid), method="REML", data=induced.sum) 

induced8a.lme <- lme (Form.induced, random=list(ID=~1),  weights=varIdent(form=~1|wellvid), method="REML", data=induced.sum) 

anova(induced.gls, induced1.lme, induced2.lme, induced2a.lme, induced3.lme, induced3a.lme, induced4.lme,
      induced6.lme, induced6a.lme, induced8.lme, induced8a.lme)

#6 and 6a, 8 and 8a produce the same result. but 6 is lower

#correlation structures

induced9.lme <- lme (Form.induced, random = ~1|ID,  weights=varIdent(form=~1|bin), correlation=corAR1(), method="REML", data=induced.sum) #BEST use this

induced10.lme <- lme (Form.induced, random = ~1|ID,  weights=varIdent(form=~1|bin), 
                      correlation=corAR1 (form=~1|ID), method="REML", data=induced.sum) 

induced11.lme <- lme (Form.induced, random = ~1|ID,  weights=varIdent(form=~1|bin), 
                      correlation=corAR1 (form=~1|ID/bin), method="REML", data=induced.sum) 


anova(induced.gls, induced1.lme, induced2.lme, induced2a.lme, induced3.lme, induced3a.lme, induced4.lme,
      induced6.lme, induced6a.lme, induced8.lme, induced8a.lme, induced9.lme, induced10.lme, induced11.lme)

#best induced11.lme, AIC=4727.753

summary(induced11.lme)
anova(induced11.lme)

#multiple comparisons
#library(lsmeans)

#pairs(lsmeans(induced11.lme, ~bin|timeminfac))

#library(multcompView)
#cld(lsmeans(induced11.lme, ~bin), alpha=0.05)

#library(multcomp)

#glht doesnt work
summary(glht(induced11.lme, linfct=mcp(bin="Tukey", covariate_average=TRUE, interaction_average = TRUE)))

#residuals
induced.E2<-resid(induced11.lme,type="normalized")
induced.F2<-fitted(induced11.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=induced.F2,y=induced.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(induced.E2~bin,data=induced.sum, main="bin",ylab=MyYlab)
plot(x=induced.sum$timemin,y=induced.E2,main="Time",ylab=MyYlab,xlab="Time (min)")
boxplot(induced.E2~wellvid,data=induced.sum, main="wellvid",ylab=MyYlab)
par(op)

xyplot (induced.E2 ~ timemin| bin, data=induced.sum, ylab="Residuals", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})



#LME: notinduced.sum

qplot(timeminfac, Vlog, color = bin, data = notinduced.sum,  geom = "boxplot") + facet_grid(bin~., scales="free") +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=1))

exp.notinduced=as.data.frame(data.table(cbind(wellvid=as.numeric(as.factor(notinduced.sum$wellvid)),bin=as.numeric(as.factor(notinduced.sum$bin)), 
                                              T=as.numeric(notinduced.sum$timeminfac), ID=as.numeric(as.factor(notinduced.sum$ID)))))

cor(exp.notinduced, method = "spearman")

vif_func(in_frame=exp.notinduced,thresh=5,trace=T)

pairs(exp.notinduced, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

#boxplots
op=par(mfrow=c(2,2))
boxplot(Vlog~wellvid, data=notinduced.sum)
boxplot(Vlog~bin, data=notinduced.sum)
boxplot(Vlog~ID, data=notinduced.sum)
boxplot (Vlog~timemin, data=notinduced.sum)

#levene
library(lawstat)
levene.test(notinduced.sum$Vlog, group=notinduced.sum$wellvid, location="mean") #significant
levene.test(notinduced.sum$Vlog, group=notinduced.sum$ID, location="mean") # significant
levene.test(notinduced.sum$Vlog, group=notinduced.sum$timeminfac, location="mean")  # significant
levene.test(notinduced.sum$Vlog, group=notinduced.sum$bin, location="mean") # significant

#fit a gls
Form.notinduced <- formula (Vlog ~ bin*timeminfac)
notinduced.gls<- gls(Form.notinduced, data=notinduced.sum)

#nlme model

#random factor
notinduced1.lme <- lme (Form.notinduced, random = ~1|ID, method="REML", data=notinduced.sum) #BEST

notinduced2.lme <- lme (Form.notinduced, random = ~1|bin, method="REML", data=notinduced.sum)

notinduced2a.lme <- lme (Form.notinduced, random=list(ID=~1), method="REML", data=notinduced.sum) #second best

notinduced3.lme <- lme (Form.notinduced, random = ~1|ID/bin, method="REML", data=notinduced.sum) #best

notinduced3a.lme <- lme (Form.notinduced, random=list(ID=~1|bin), method="REML", data=notinduced.sum)

notinduced4.lme <- lme (Form.notinduced, random = ~1|wellvid, method="REML", data=notinduced.sum)

anova(notinduced.gls, notinduced1.lme, notinduced2.lme, notinduced2a.lme, notinduced3.lme, notinduced3a.lme, notinduced4.lme) #2a and 3a the same


#variance structure

#notinduced5.lme <- lme (Form.notinduced, random = ~1|ID,  weights=varIdent(form=~1|ID), method="REML", data=notinduced.sum)

notinduced6.lme <- lme (Form.notinduced, random = ~1|ID,  weights=varIdent(form=~1|bin), method="REML", data=notinduced.sum) 

notinduced6a.lme <- lme (Form.notinduced,  random=list(ID=~1),  weights=varIdent(form=~1|bin), method="REML", data=notinduced.sum) 

#notinduced7.lme <- lme (Form.notinduced, random = ~1|ID,  weights=varIdent(form=~1|ID/bin), method="REML", data=notinduced.sum) 

notinduced8.lme <- lme (Form.notinduced, random = ~1|ID,  weights=varIdent(form=~1|wellvid), method="REML", data=notinduced.sum) 

notinduced8a.lme <- lme (Form.notinduced, random=list(ID=~1),  weights=varIdent(form=~1|wellvid), method="REML", data=notinduced.sum) 

anova(notinduced.gls, notinduced1.lme, notinduced2.lme, notinduced2a.lme, notinduced3.lme, notinduced3a.lme, notinduced4.lme,
      notinduced6.lme, notinduced6a.lme, notinduced8.lme, notinduced8a.lme)

#6 and 6a, 8 and 8a produce the same result. but 6 is lower

#correlation structures

notinduced9.lme <- lme (Form.notinduced, random = ~1|ID,  weights=varIdent(form=~1|bin), correlation=corAR1(), method="REML", data=notinduced.sum) #BEST use this

notinduced10.lme <- lme (Form.notinduced, random = ~1|ID,  weights=varIdent(form=~1|bin), 
                         correlation=corAR1 (form=~1|ID), method="REML", data=notinduced.sum) 

notinduced11.lme <- lme (Form.notinduced, random = ~1|ID,  weights=varIdent(form=~1|bin), 
                         correlation=corAR1 (form=~1|ID/bin), method="REML", data=notinduced.sum) 


anova(notinduced.gls, notinduced1.lme, notinduced2.lme, notinduced2a.lme, notinduced3.lme, notinduced3a.lme, notinduced4.lme,
      notinduced6.lme, notinduced6a.lme, notinduced8.lme, notinduced8a.lme, notinduced9.lme, notinduced10.lme, notinduced11.lme)

#best notinduced11.lme, AIC=6063.907

summary(notinduced11.lme)
anova(notinduced11.lme)

#multiple comparisons
#library(lsmeans)

#pairs(lsmeans(notinduced11.lme, ~bin|timeminfac))

#library(multcompView)
#cld(lsmeans(notinduced11.lme, ~bin), alpha=0.05)

#library(multcomp)

#glht doesnt work
summary(glht(notinduced11.lme, linfct=mcp(bin="Tukey", covariate_average=TRUE, interaction_average = TRUE)))

#residuals
notinduced.E2<-resid(notinduced11.lme,type="normalized")
notinduced.F2<-fitted(notinduced11.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=notinduced.F2,y=notinduced.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(notinduced.E2~bin,data=notinduced.sum, main="bin",ylab=MyYlab)
plot(x=notinduced.sum$timemin,y=notinduced.E2,main="Time",ylab=MyYlab,xlab="Time (min)")
boxplot(notinduced.E2~wellvid,data=notinduced.sum, main="wellvid",ylab=MyYlab)
par(op)

xyplot (notinduced.E2 ~ timemin| bin, data=notinduced.sum, ylab="Residuals", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})

# data not fitted because V should be shown and not Vlog

#plot

mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="bin") { 
    value[value=="bin_DPR"] <- "DPR bead"
    value[value=="bin_dSi"]   <- "dSi bead"
    value[value=="bin_out"] <- "outside"
  }
  return(value)
}

scaleFUN <- function(x) sprintf("%.1f", x)

trackdata.sumV.noID <- summarySE(trackdata.sumV, measurevar="V", groupvars=c("treatment", "bin","timemin", "timeminfac"), na.rm=TRUE)

trackdata.sumV.noID$bin2 <- factor(trackdata.sumV.noID$bin, levels=c("bin_DPR", "bin_dSi", "bin_out"), labels =c ("Diproline bead", "dSi bead", "outside beads"))
trackdata.sumV.noID$treatment2 <-  factor(trackdata.sumV.noID$treatment, levels=c("Si_induced", "Si_notinduced"), labels =c ("induced", "not induced"))

#bw
resize.win(12,9)

ggplot(data=trackdata.sumV.noID, aes(x=timemin, y=V)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=V-se, ymax=V+se), width=0.5, size=1) +
  facet_grid(treatment2~bin2, scale="free") +
  labs(list(x = "Time (min)", y = "Speed (µm/s)"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text = element_text(size=15), legend.title=text, legend.text=text, panel.margin=unit (1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + 
  scale_x_continuous (breaks=c(0, 2, 4, 6, 8, 10))+
  scale_y_continuous(labels=scaleFUN)



##no fitting


library(AICcmodavg)

trackdata.sumV <- summarySE(trackdata2, measurevar="V", groupvars=c("wellvid", "treatment", "bin", "ID", "timemin", "timeminfac"), na.rm=TRUE)

induced.sumV <- subset(trackdata.sumV, trackdata.sum$treatment=="Si_induced")

notinduced.sumV <- subset(trackdata.sumV, trackdata.sum$treatment=="Si_notinduced")


#induced fit
induced.fit <- as.data.frame(predictSE.lme(induced9.lme, induced.sum, se.fit = TRUE, level = 0,
                                           print.matrix = FALSE))

induced.fit$upr <- induced.fit$fit + (1.96 * induced.fit$se.fit)
induced.fit$lwr <- induced.fit$fit - (1.96 * induced.fit$se.fit)

induced.fit.combdata <- cbind(induced.sum, induced.fit)

#notinduced fit
notinduced.fit <- as.data.frame(predictSE.lme(notinduced9.lme, notinduced.sum, se.fit = TRUE, level = 0,
                                              print.matrix = FALSE))

notinduced.fit$upr <- notinduced.fit$fit + (1.96 * notinduced.fit$se.fit)
notinduced.fit$lwr <- notinduced.fit$fit - (1.96 * notinduced.fit$se.fit)

notinduced.fit.combdata <- cbind(notinduced.sum, notinduced.fit)

#summaries, Vlog is used because its summarized data but only for now

#summaries

induced.sum.noID <- summarySE(induced.sum, measurevar="Vlog", groupvars=c("treatment", "bin", "timeminfac", "timemin"), na.rm=TRUE)

notinduced.sum.noID <- summarySE(notinduced.sum, measurevar="Vlog", groupvars=c("treatment", "bin", "timeminfac", "timemin"), na.rm=TRUE)





