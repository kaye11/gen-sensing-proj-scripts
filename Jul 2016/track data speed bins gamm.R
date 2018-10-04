#read data
trackdata <- read.table ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/final track data.csv", 
                         sep=";", header=T)

library(ggplot2)
library(gtable)
library(ggthemes)
library(mgcv)
library(ggplot2)
library(grid)
library(data.table)
s=mgcv:::s

source("AED.R")
source("vif.R")
source("tsDiagGamm.R")
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

trackdata <- trackdata [! trackdata$time=="0",  ]

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
#7929 spots (control), 14770 (P bead), 11907 (Si bead)
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

qplot(time2, Vlog, color = treatment, data = trackdata,  geom = "boxplot") + facet_grid(treatment~newbin, scales="free") 

#s to min
trackdata$timemin <- as.numeric(trackdata$time2)*120/60


#summaries
speed.sumV <- summarySE(trackdata, measurevar="V", groupvars=c("treatment", "newbin", "time2"))
speed.sumVlog <- summarySE(trackdata, measurevar="Vlog", groupvars=c("treatment", "newbin", "time2"))

#Bin AB
BinsAB= subset (trackdata, newbin=='BinsAB')

BinsAB <- BinsAB [! BinsAB$time2n=="0",  ]

qplot(time2, Vlog, color = treatment, data = BinsAB,  geom = "boxplot") + facet_grid(treatment~., scales="free") 

expAB=as.data.frame(data.table(cbind(treatment=as.numeric(as.factor(BinsAB$treatment)), T=as.numeric(BinsAB$time2), ID=as.numeric(as.factor(BinsAB$ID)))))

cor(expAB, method = "spearman")

vif_func(in_frame=expAB,thresh=5,trace=T)

pairs(expAB, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

SpeedBinsAB <- summarySE(BinsAB, measurevar="V", groupvars=c("treatment", "time2", "ID"))

#boxplots
op=par(mfrow=c(2,2))
boxplot(Vlog~treatment, data=BinsAB)
boxplot(Vlog~ID, data=BinsAB)
boxplot (Vlog~time2, data=BinsAB)

#levene
library(lawstat)
levene.test(BinsAB$Vlog, group=BinsAB$ID, location="mean") #unequal
levene.test(BinsAB$Vlog, group=BinsAB$time2, location="mean") #unequal
levene.test(BinsAB$Vlog, group=BinsAB$treatment, location="mean") #unequal


#gamm
BA <- gamm (Vlog~s(time2n, by=treatment, bs="fs"), method="REML", data = BinsAB)

#smoothing splines for factors can are 3 types cr (cubic regression), cs (shrinkage version of cr) and cc (cyclic cubic)
BA1 <- gamm (Vlog~s(time2n, by=treatment, bs="fs", xt="cr"), method="REML", data = BinsAB) 
BA2 <- gamm (Vlog~s(time2n, by=treatment, bs="fs", xt="cs"), method="REML", data = BinsAB) 
BA2.1 <- gamm (Vlog~s(time2n, by=treatment, bs="fs", xt="cc"), method="REML", data = BinsAB) #best but use cr because it yields better results with corAR1
anova(BA$lme, BA1$lme, BA2$lme, BA2.1$lme)

#make random factor and correlations
fBinsAB <- Vlog~s(time2n, by=treatment, bs="fs", xt="cr")
fBinsABmin <- Vlog~s(timemin, by=treatment, bs="fs", xt="cr")

BA3 <- gamm (fBinsAB, method="REML",  random=list(ID=~1), data = BinsAB) 
BA4 <- gamm (fBinsAB, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), data = BinsAB) #BEST
BA5 <- gamm (fBinsAB, method="REML", random=list(ID=~1), correlation= corAR1 (), data = BinsAB) #same with BA4

anova(BA$lme, BA1$lme, BA2$lme, BA2.1$lme, BA3$lme, BA4$lme, BA5$lme)

#make variance structures
#BA6 <- gamm (fBinsAB, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varIdent(form=~1| time2), data = BinsAB) #no convergence

#BA7 <- gamm (fBinsAB, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varIdent(form=~1| ID), data = BinsAB) #no convergence

BA8 <- gamm (fBinsAB, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), 
             weights = varIdent(form=~1| treatment), data = BinsAB) 

BA8min <- gamm (fBinsABmin, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), 
             weights = varIdent(form=~1| treatment), data = BinsAB) 

BA9 <- gamm (fBinsAB, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varExp(form=~fitted(.)), data = BinsAB) 

#BA10 <- gamm (fBinsAB, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varExp(0.3, form=~time2n), data = BinsAB) 

#BA9 worked but kind of weirdly

anova(BA$lme, BA1$lme, BA2$lme, BA2.1$lme, BA3$lme, BA4$lme, BA5$lme, BA8$lme)

AIC(BA$lme, BA1$lme, BA2$lme, BA2.1$lme, BA3$lme, BA4$lme, BA5$lme, BA8$lme, BA9$lme)

#Take BA8
with(BinsAB, tsDiagGamm(BA8, timevar=time2n, observed=Vlog)) #visual checking of residuals

plot(BA8$lme)

summary(BA8$gam)
summary(BA8$lme)
anova(BA8$gam)
anova(BA8$lme)

#checking the model produced by GAMM
resize.win(9,12)
op=par(mfrow=c(3,1), mar=c(4.5,4.5,1.5,1.5))
#plot(BA8$gam, cex.lab=1.1, cex.axis=1.1, xlab ="Time (s)")
plot(BA8min$gam, cex.lab=1.5, cex.axis=1.5, xlab ="Time (min)")


#Bin C
BinC= subset (trackdata, newbin=='BinC')

BinC <- BinC [! BinC$time2n=="0",  ]

qplot(time2, Vlog, color = treatment, data = BinC,  geom = "boxplot") + facet_grid(treatment~., scales="free") 

expC=as.data.frame(data.table(cbind(treatment=as.numeric(as.factor(BinC$treatment)), T=as.numeric(BinC$time2), ID=as.numeric(as.factor(BinC$ID)))))

cor(expC, method = "spearman")

vif_func(in_frame=expC,thresh=5,trace=T)

pairs(expC, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

SpeedBinC <- summarySE(BinC, measurevar="V", groupvars=c("treatment", "time2", "ID"))

#boxplots
op=par(mfrow=c(2,2))
boxplot(Vlog~treatment, data=BinC)
boxplot(Vlog~ID, data=BinC)
boxplot (Vlog~time2, data=BinC)

#levene
library(lawstat)
levene.test(BinC$Vlog, group=BinC$ID, location="mean") #unequal
levene.test(BinC$Vlog, group=BinC$time2, location="mean") #unequal
levene.test(BinC$Vlog, group=BinC$treatment, location="mean") #unequal


#gamm
BC <- gamm (Vlog~s(time2n, by=treatment, bs="fs"), method="REML", data = BinC)

#smoothing splines for factors can are 3 types cr (cubic regression), cs (shrinkage version of cr) and cc (cyclic cubic)
BC1 <- gamm (Vlog~s(time2n, by=treatment, bs="fs", xt="cr"), method="REML", data = BinC) 
BC2 <- gamm (Vlog~s(time2n, by=treatment, bs="fs", xt="cs"), method="REML", data = BinC) 
BC2.1 <- gamm (Vlog~s(time2n, by=treatment, bs="fs", xt="cc"), method="REML", data = BinC) #best but use cr because it yields better results with corAR1
anova(BC$lme, BC1$lme, BC2$lme, BC2.1$lme)

#make random factor and correlations
fBinC <- Vlog~s(time2n, by=treatment, bs="fs", xt="cr")
fBinCmin <- Vlog~s(timemin, by=treatment, bs="fs", xt="cr")

BC3 <- gamm (fBinC, method="REML",  random=list(ID=~1), data = BinC) 
BC4 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), data = BinC) #BEST
BC5 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (), data = BinC) #same with BC4

anova(BC$lme, BC1$lme, BC2$lme, BC2.1$lme, BC3$lme, BC4$lme, BC5$lme)

#make variance structures
#BC6 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varIdent(form=~1| time2), data = BinC) #no convergence

#BC7 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varIdent(form=~1| ID), data = BinC) #no convergence

BC8 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), 
             weights = varIdent(form=~1| treatment), data = BinC) 

BC8min <- gamm (fBinCmin, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), 
             weights = varIdent(form=~1| treatment), data = BinC) 


BC9 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varExp(form=~fitted(.)), data = BinC) 

#BC10 <- gamm (fBinC, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varExp(0.3, form=~time2n), data = BinC) 

#BC9 worked but kind of weirdly

anova(BC$lme, BC1$lme, BC2$lme, BC2.1$lme, BC3$lme, BC4$lme, BC5$lme, BC8$lme)

AIC(BC$lme, BC1$lme, BC2$lme, BC2.1$lme, BC3$lme, BC4$lme, BC5$lme, BC8$lme, BC9$lme)

#Take BC8
with(BinC, tsDiagGamm(BC8, timevar=time2n, observed=Vlog)) #visual checking of residuals

plot(BC8$lme)

summary(BC8$gam)
summary(BC8$lme)
anova(BC8$gam)
anova(BC8$lme)

#checking the model produced by GAMM
resize.win (6,9)
op=par(mfrow=c(3,1), mar=c(4.5,4.5,1.5,1.5))
#plot(BC8$gam, cex.lab=1.1, cex.axis=1.1, xlab ="Time (s)")
plot(BC8min$gam, cex.lab=1.5, cex.axis=1.5, xlab ="Time (min)")

##PLOTTING 

#plotting
grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 

resize.win(12,9)


mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="newbin") { 
    value[value=="BinsAB"] <- "Bins A+B"
    value[value=="BinC"]   <- "Bin C"
    }
  return(value)
}


cbPalette <- c("#999999",  "#009E73", "#56B4E9","#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#set time in min
speed.sumV$timemin <- as.numeric(speed.sumV$time2)*120/60

ggplot(data=speed.sumV, aes(x=time2, y=V, shape=treatment, color=treatment, group=treatment)) + geom_point(size=5)+
  geom_errorbar(aes(ymin=V-se, ymax=V+se), width=0.5, size=1) + facet_grid(newbin~., labeller=mf_labeller)+
  scale_color_manual(values = cbPalette, name="Treatment") +
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") +
  labs(list(x = "Time (s)", y = "Mean cell speed (µm/s)"))+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + 
  scale_x_discrete (breaks=c(seq(120, 3600, 720))) 

#+ geom_line(size=1.5)


#facet_grid

ggplot(data=speed.sumV, aes(x=time2, y=V, shape=treatment, color=treatment, group=treatment)) + geom_point(size=5)+
  geom_errorbar(aes(ymin=V-se, ymax=V+se), width=0.5, size=1) + facet_grid(treatment~newbin, labeller=mf_labeller, scale="free")+
  scale_color_manual(values = cbPalette, name="Treatment") +
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") +
  labs(list(x = "Time (s)", y = "Mean cell speed (µm/s)"))+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + 
  scale_x_discrete (breaks=c(seq(120, 3600, 720)))

#set time in min
speed.sumV$timemin <- as.numeric(speed.sumV$time2)*120/60

ggplot(data=speed.sumV, aes(x=timemin, y=V, shape=treatment, color=treatment, group=treatment)) + geom_point(size=5)+
  geom_errorbar(aes(ymin=V-se, ymax=V+se), width=0.5, size=1) + facet_grid(newbin~., labeller=mf_labeller)+
  scale_color_manual(values = cbPalette, name="Treatment") +
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") +
  labs(list(x = "Time (min)", y = "Mean cell speed (µm/s)"))+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))

ggplot(data=speed.sumV, aes(x=timemin, y=V, shape=treatment, color=treatment, group=treatment)) + geom_point(size=5)+
  geom_errorbar(aes(ymin=V-se, ymax=V+se), width=2, size=1) + facet_grid(treatment~newbin, labeller=mf_labeller, scale="free")+
  scale_color_manual(values = cbPalette, name="Treatment") +
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") +
  labs(list(x = "Time (min)", y = "Mean cell speed (µm/s)"))+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))


#bw

speed.sumV$treatlabels <- factor(speed.sumV$treatment, levels= c("control bead", "P bead", "Si bead"), labels = c("control bead", "dP bead", "dSi bead"))


ggplot(data=speed.sumV, aes(x=timemin, y=V)) + geom_point(size=5)+
  geom_errorbar(aes(ymin=V-se, ymax=V+se), width=2, size=1) + facet_grid(treatlabels~newbin, labeller=mf_labeller, scale="free")+
  labs(list(x = "Time (min)", y = "Mean cell speed (µm/s)"))+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))

##save trackdata

write.table (trackdata, "d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/final track data with preprocessing.csv", 
             sep=";", col.names=T, row.names=F)