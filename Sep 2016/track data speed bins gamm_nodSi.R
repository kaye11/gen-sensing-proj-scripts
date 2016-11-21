#read data
trackdata3 <- read.table ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/final track data.csv", 
                         sep=";", header=T)

library(ggplot2)
library(gtable)
library(ggthemes)
library(mgcv)
library(ggplot2)
library(grid)
library(data.table)
library(gdata)
s=mgcv:::s

source("AED.R")
source("vif.R")
source("tsDiagGamm.R")
source("summarySE.R")
source("resizewin.R")

resize.win(12,9)

#drop dSi
trackdata2 <- trackdata3 [! trackdata3$treatment=="Si bead", ]
trackdata <- trackdata2
trackdata$treatment <- trackdata$treatment [, drop=TRUE]

#treat time as factor
trackdata$timef <- as.factor(trackdata$time)


#combine A+B
trackdata$newbin <- trackdata$bin
levels(trackdata$newbin) <- c("BinsAB", "BinsAB", "BinC")

#time2
tm=seq(0, 3600, by = 120)
trackdata$time2 <- cut(trackdata$T, tm, labels=paste(tail(tm, -1L)))

trackdata$time2 = factor (trackdata$time2, levels=c(levels(trackdata$time2), 0))

trackdata$time2[is.na(trackdata$time2)] <- 0 #replace NAs with 0. some data time points have the starting point as NA because they started at 0
trackdata$time2n <- as.numeric (trackdata$time2)*120

trackdata <- trackdata [! trackdata$time=="0",  ]
trackdata$timemin <- as.numeric(trackdata$time2n)/120*2


speed.trackdata <- summarySE(trackdata,  measurevar="Vlog", groupvars=c("newbin", "time2", "ID"))
trackdata.slow <- speed.trackdata [speed.trackdata$Vlog< 1, ]


A3control <- trackdata [trackdata$wellvid =="A3-004", ]
speed.A3control <- summarySE(A3control,  measurevar="Vlog", groupvars=c("newbin", "time2", "ID"))
A3control.slow <- speed.A3control [speed.A3control$V< 1, ]

B4control <- trackdata [trackdata$wellvid =="B4-026", ]
speed.B4control <- summarySE(B4control,  measurevar="Vlog", groupvars=c("newbin", "time2", "ID"))
B4control.slow <- speed.B4control [speed.B4control$V< 1, ]

B1control <- trackdata [trackdata$wellvid =="B1-032", ]
speed.B1control <- summarySE(B1control,  measurevar="Vlog", groupvars=c("newbin", "time2", "ID"))
B1control.slow <- speed.B1control [speed.B1control$V< 1, ]


#8 tracks deleted from control, control now has 37 tracks
#7929 spots (control), 14770 (P bead)
trackdata <- trackdata [! trackdata$ID=="A3-004-0", ]
trackdata <- trackdata [! trackdata$ID=="A3-004-469", ]
trackdata <- trackdata [! trackdata$ID=="A3-004-946", ]
trackdata <- trackdata [! trackdata$ID=="A3-004-6", ]

trackdata <- trackdata [! trackdata$ID=="B4-026-739", ]

trackdata <- trackdata [! trackdata$ID=="B1-032-119", ]
trackdata <- trackdata [! trackdata$ID=="B1-032-532", ]
trackdata <- trackdata [! trackdata$ID=="B1-032-34", ]

#trackdata <- trackdata [! trackdata$ID=="B4-026-126", ]

#plotting

qplot(as.factor(timemin), Vlog, color = treatment, data = trackdata,  geom = "boxplot") + facet_grid(treatment~newbin, scales="free") 

#summaries
speed.sumV <- summarySE(trackdata, measurevar="V", groupvars=c("treatment", "newbin", "timemin"))
speed.sumVlog <- summarySE(trackdata, measurevar="Vlog", groupvars=c("treatment", "newbin", "timemin"))

#Bin AB
BinsAB= subset (trackdata, newbin=='BinsAB')

BinsAB <- BinsAB [! BinsAB$time2n=="0",  ]

qplot(as.factor(timemin), Vlog, color = treatment, data = BinsAB,  geom = "boxplot") + facet_grid(treatment~., scales="free") 

expAB=as.data.frame(data.table(cbind(treatment=as.numeric(as.factor(BinsAB$treatment)), T=as.numeric(BinsAB$time2), ID=as.numeric(as.factor(BinsAB$ID)))))

cor(expAB, method = "spearman")

vif_func(in_frame=expAB,thresh=5,trace=T)

pairs(expAB, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

SpeedBinsAB <- summarySE(BinsAB, measurevar="V", groupvars=c("treatment", "timemin", "ID"))

#boxplots
op=par(mfrow=c(2,2))
boxplot(Vlog~treatment, data=BinsAB)
boxplot(Vlog~ID, data=BinsAB)
boxplot (Vlog~timemin, data=BinsAB)

#levene
library(lawstat)
levene.test(BinsAB$Vlog, group=BinsAB$ID, location="mean") #unequal
levene.test(BinsAB$Vlog, group=BinsAB$timemin, location="mean") #unequal
levene.test(BinsAB$Vlog, group=BinsAB$treatment, location="mean") #unequal

#gamm
BA <- gamm (Vlog~s(timemin, by=treatment, bs="fs"), method="REML", data = BinsAB)

#smoothing splines for factors can are 3 types cr (cubic regression), cs (shrinkage version of cr) and cc (cyclic cubic)
BA1 <- gamm (Vlog~s(timemin, by=treatment, bs="fs", xt="cr"), method="REML", data = BinsAB) 
BA2 <- gamm (Vlog~s(timemin, by=treatment, bs="fs", xt="cs"), method="REML", data = BinsAB) 
BA2.1 <- gamm (Vlog~s(timemin, by=treatment, bs="fs", xt="cc"), method="REML", data = BinsAB) #best but use cr because it yields better results with corAR1
anova(BA$lme, BA1$lme, BA2$lme, BA2.1$lme)

#make random factor and correlations
fBinsAB <- Vlog~s(timemin, by=treatment, bs="fs", xt="cr")

BA3 <- gamm (fBinsAB, method="REML",  random=list(ID=~1), data = BinsAB) 
BA4 <- gamm (fBinsAB, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), data = BinsAB) #BEST
BA5 <- gamm (fBinsAB, method="REML", random=list(ID=~1), correlation= corAR1 (), data = BinsAB) #same with BA4

anova(BA$lme, BA1$lme, BA2$lme, BA2.1$lme, BA3$lme, BA4$lme, BA5$lme)

#make variance structures
#BA6 <- gamm (fBinsAB, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varIdent(form=~1| timemin), data = BinsAB) #no convergence

#BA7 <- gamm (fBinsAB, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varIdent(form=~1| ID), data = BinsAB) #no convergence

BA8 <- gamm (fBinsAB, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), 
             weights = varIdent(form=~1| treatment), data = BinsAB) 

BA8a <- gamm (Vlog~s(timemin, by=treatment, bs="fs", xt="cc"), method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), 
              weights = varIdent(form=~1| treatment), data = BinsAB) 

BA9 <- gamm (fBinsAB, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varExp(form=~fitted(.)), data = BinsAB) 

#BA10 <- gamm (fBinsAB, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varExp(0.3, form=~timemin), data = BinsAB) 

#BA9 worked but kind of weirdly

anova(BA$lme, BA1$lme, BA2$lme, BA2.1$lme, BA3$lme, BA4$lme, BA5$lme, BA8$lme, BA8a$lme)

AIC(BA$lme, BA1$lme, BA2$lme, BA2.1$lme, BA3$lme, BA4$lme, BA5$lme, BA8$lme, BA9$lme)

#Take BA8
with(BinsAB, tsDiagGamm(BA8, timevar=timemin, observed=Vlog)) #visual checking of residuals

plot(BA8$lme)

summary(BA8$gam)
summary(BA8$lme)
anova(BA8$gam)
anova(BA8$lme)

#checking the model produced by GAMM
resize.win(9,12)
op=par(mfrow=c(3,1), mar=c(4.5,4.5,1.5,1.5))
#plot(BA8$gam, cex.lab=1.1, cex.axis=1.1, xlab ="Time (s)")
plot(BA8$gam, cex.lab=1.5, cex.axis=1.5, xlab ="Time (min)")


#Bin C
BinC= subset (trackdata, newbin=='BinC')

BinC <- BinC [! BinC$timemin=="0",  ]

qplot(as.factor(timemin), Vlog, color = treatment, data = BinC,  geom = "boxplot") + facet_grid(treatment~., scales="free") 

expC=as.data.frame(data.table(cbind(treatment=as.numeric(as.factor(BinC$treatment)), T=as.numeric(BinC$timemin), ID=as.numeric(as.factor(BinC$ID)))))

cor(expC, method = "spearman")

vif_func(in_frame=expC,thresh=5,trace=T)

pairs(expC, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

SpeedBinC <- summarySE(BinC, measurevar="V", groupvars=c("treatment", "timemin", "ID"))

#boxplots
op=par(mfrow=c(2,2))
boxplot(Vlog~treatment, data=BinC)
boxplot(Vlog~ID, data=BinC)
boxplot (Vlog~timemin, data=BinC)

#levene
library(lawstat)
levene.test(BinC$Vlog, group=BinC$ID, location="mean") #unequal
levene.test(BinC$Vlog, group=BinC$time2, location="mean") #unequal
levene.test(BinC$Vlog, group=BinC$treatment, location="mean") #unequal


#gamm
BC <- gamm (Vlog~s(timemin, by=treatment, bs="fs"), method="REML", data = BinC)

#smoothing splines for factors can are 3 types cr (cubic regression), cs (shrinkage version of cr) and cc (cyclic cubic)
BC1 <- gamm (Vlog~s(timemin, by=treatment, bs="fs", xt="cr"), method="REML", data = BinC) 
BC2 <- gamm (Vlog~s(timemin, by=treatment, bs="fs", xt="cs"), method="REML", data = BinC) 
BC2.1 <- gamm (Vlog~s(timemin, by=treatment, bs="fs", xt="cc"), method="REML", data = BinC) #best but use cr because it yields better results with corAR1
anova(BC$lme, BC1$lme, BC2$lme, BC2.1$lme)

#make random factor and correlations
fBinC <- Vlog~s(timemin, by=treatment, bs="fs", xt="cr")

BC3 <- gamm (fBinC, method="REML",  random=list(ID=~1), data = BinC) 
BC4 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), data = BinC) #BEST
BC5 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (), data = BinC) #same with BC4

anova(BC$lme, BC1$lme, BC2$lme, BC2.1$lme, BC3$lme, BC4$lme, BC5$lme)

#make variance structures
#BC6 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varIdent(form=~1| time2), data = BinC) #no convergence

#BC7 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varIdent(form=~1| ID), data = BinC) #no convergence

BC8 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), 
             weights = varIdent(form=~1| treatment), data = BinC) 

BC8a <- gamm (Vlog~s(timemin, by=treatment, bs="fs", xt="cc"), method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), 
              weights = varIdent(form=~1| treatment), data = BinC) 

BC9 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varExp(form=~fitted(.)), data = BinC) 

#BC10 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varExp(0.3, form=~timemin), data = BinC) 

#BC9 worked but kind of weirdly

anova(BC$lme, BC1$lme, BC2$lme, BC2.1$lme, BC3$lme, BC4$lme, BC5$lme, BC8$lme, BC8a$lme)

AIC(BC$lme, BC1$lme, BC2$lme, BC2.1$lme, BC3$lme, BC4$lme, BC5$lme, BC8$lme, BC9$lme)

#Take BC8
with(BinC, tsDiagGamm(BC8, timevar=timemin, observed=Vlog)) #visual checking of residuals

plot(BC8$lme)

summary(BC8$gam)
summary(BC8$lme)
anova(BC8$gam)
anova(BC8$lme)

#checking the model produced by GAMM
resize.win (9,12)
op=par(mfrow=c(3,1), mar=c(4.5,4.5,1.5,1.5))
#plot(BC8$gam, cex.lab=1.1, cex.axis=1.1, xlab ="Time (s)")
plot(BC8$gam, cex.lab=1.5, cex.axis=1.5, xlab ="Time (min)")

##PLOTTING 

#plotting
grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 

resize.win(12,9)

#bw

speed.sumV <- summarySE(trackdata, measurevar="V", groupvars=c("treatment", "newbin", "timemin"))

speed.sumV$treatlabels <- factor(speed.sumV$treatment, levels= c("control bead", "P bead"), labels = c(" control bead  ", " dP bead  "))
speed.sumV$newbinlabels <- factor(speed.sumV$newbin, levels=c ("BinsAB", "BinC"), labels =c ("Bins A+B", "Bin C"))


ggplot(data=speed.sumV, aes(x=timemin, y=V, shape=treatlabels)) +  geom_errorbar(aes(ymin=V-se, ymax=V+se), width=2, size=1) +
  geom_point(size=5, shape = 21, color='black', aes(fill = treatlabels)) + 
  facet_grid(newbinlabels~treatlabels)+
  labs(y = expression("Mean cell speed"~("µm"~s^-1)), x = "Time (min)",
       title = "Speed of dP-starved cells in response to control and dP-loaded beads" ) +
  scale_fill_manual(values = c('white', 'black')) +
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20, vjust=1.5), 
        axis.title.x=element_text(size=20, vjust=-0.5),
        plot.title = element_text(size =24), axis.text=text,  legend.position="bottom", legend.title=element_blank(),
        strip.text.x = element_blank(), strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))

##save trackdata

write.table (trackdata, "d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/final track data with preprocessing_nodSi.csv", 
             sep=";", col.names=T, row.names=F)

