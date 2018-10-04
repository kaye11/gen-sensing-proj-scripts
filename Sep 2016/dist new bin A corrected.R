
#read data
#trackdata <- read.table ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/final track data with preprocessing.csv", sep=";", header=T)

trackdata <- read.table ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/final track data.csv", 
                         sep=";", header=T)

coordinates <- read.table ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/coordinates.csv", 
                         sep=";", header=T)


#treat time as factor
trackdata$timef <- as.factor(trackdata$time)
trackdata$time <- as.numeric(trackdata$time)

library (CircStats)

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


trackdata$timemin <- as.numeric(trackdata$time2)*120/60

library(plyr)
library(dplyr)
library(ggplot2)

#merging coordinates and data set (trackdata)
trackdata <- merge(trackdata, coordinates, by="wellvid")

#calculate distance between the last cell position and bead
trackdata$dist=sqrt(((trackdata$bead.X-trackdata$X)^2)+((trackdata$bead.Y-trackdata$Y)^2))

##distance correction of bin

trackdata$distc <- trackdata$dist

binA <- subset (trackdata, trackdata$bin=="binA", )
binB <- subset (trackdata, trackdata$bin=="binB", )
binC <- subset (trackdata, trackdata$bin=="binC", )

binA$distc = binA$dist - binA$radius

trackdata.corrected <- rbind(binA, binB, binC)

#compute whole distance

trial.dist <- aggregate(trackdata.corrected$distc, by = list(ID = trackdata.corrected$ID, treatment=trackdata.corrected$treatment, bin=trackdata.corrected$bin), last)

qplot(treatment, x, data=trial.dist, geom="boxplot")+facet_grid(~bin)

#summary
source("summarySE.R")

trial.dist$dist <- trial.dist$x

#dist.sum <- summarySE(trial.dist, measurevar="dist", groupvars=c("treatment", "bin"), na.rm=TRUE)

dist.sum <- ddply(trial.dist, c("treatment", "bin"), summarise,
                  N    = length(ID),
                  mean = mean(dist, na.rm=TRUE),
                  sumdist= sum(dist, na.rm=TRUE), 
                  sd   = sd(dist, na.rm=TRUE),
                  se   = sd / sqrt(N))

ggplot(data=dist.sum, aes(x=treatment, y=mean, shape=treatment, color=treatment)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1) + facet_grid(bin~., scales="free")


source("vif.R")
source ("AED.R")
source("lang.R")
source("summarySE.R")
source("lang.R")
source("tsDiagGamm.R")
source("resizewin.R")
resize.win(12,9)

#no bins will be used since data points are not that much
#data will be compared accordingly as a whole


expdist=as.data.frame(data.table(cbind(treatment=trial.dist$treatment, bin=trial.dist$bin, ID=trial.dist$ID)))
cor(expdist, method = "spearman")

vif_func(in_frame=expdist,thresh=5,trace=T)

pairs(expdist, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

library(lawstat)
levene.test(trial.dist$dist, group=trial.dist$ID, location="mean") #not significant
levene.test(trial.dist$dist, group=trial.dist$bin, location="mean") # significant
levene.test(trial.dist$dist, group=trial.dist$treatment, location="mean") #not significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(dist~treatment, data=trial.dist)
boxplot(dist~ID, data=trial.dist)
boxplot (dist~bin, data=trial.dist)

#fit a gls
Form <- formula (dist ~ treatment*bin)
dist.gls<- gls(Form, data=trial.dist)


#nlme model
dist1.lme <- lme (Form, random = ~1|ID, method="REML", data=trial.dist)

#dist2.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|ID), method="REML", data=trial.dist)

#dist21.lme <- lme (Form, random = ~1|treatment/ID,  weights=varIdent(form=~1|ID), method="REML", data=trial.dist) 

#dist3.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|ID), correlation=corAR1 (form=~1|ID/treatment), method="REML", data=trial.dist) 

dist4.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|bin), correlation=corAR1 (), method="REML", data=trial.dist) #BEST

dist5.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (), method="REML", data=trial.dist) 

dist6.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (form=~1|ID), method="REML", data=trial.dist) #same as 5


anova(dist.gls, dist1.lme, dist4.lme, dist5.lme, dist6.lme)

summary(dist4.lme)
anova(dist4.lme)

#residuals
dist.E2<-resid(dist4.lme,type="normalized")
dist.F2<-fitted(dist4.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=dist.F2,y=dist.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(dist.E2~treatment,data=trial.dist, main="Treatment",ylab=MyYlab)
plot(x=trial.dist$bin,y=dist.E2,main="Bin",ylab=MyYlab,xlab="Bin")
par(op)


#multimple comparisons
#use lstrends - lstrends for estimating and comparing the slopes of fitted lines (or curves).
#https://cran.r-project.org/web/packages/lsmeans/vignettes/using-lsmeans.pdf

library(lsmeans)

#lstrends(dist4.lme, pairwise ~treatment, var="as.numeric(bin)")
#lstrends only work for numeric not for factors, therefore use lsmeans


lsmeans(dist4.lme, list(pairwise ~ treatment|bin))
pairs(lsmeans(dist4.lme, ~treatment))


#plot

grid.newpage()
text <- element_text(size = 15) #change the size of the axes
theme_set(theme_bw()) 

cbPalette <- c("#999999",  "#009E73", "#56B4E9","#F0E442", "#0072B2", "#D55E00", "#CC79A7")

mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="bin") { 
    value[value=="binA"] <- "Bin A"
    value[value=="binB"]   <- "Bin B"
    value[value=="binC"]   <- "Bin C"
  }
  return(value)
}

dist.sum$treatlabels <- factor(dist.sum$treatment, levels= c("control bead", "P bead", "Si bead"), labels = c("control bead", "dP bead", "dSi bead"))

resize.win(6,9)

#bw

dist.sum$sigbars <- rep(c(1:2, 1), each=3)

ggplot(data=dist.sum, aes(x=treatlabels, y=mean, shape=factor(sigbars))) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1) + facet_grid(bin~., scales="free", labeller=mf_labeller) +
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") + 
  labs(list(x="Treatment", y = "Mean distance relative to the bead (µm)"))+ 
  theme(axis.text=text, axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))




ggplot(data=dist.sum, aes(x=treatment, y=mean, shape=treatment, color=treatment)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1) + facet_grid(bin~., scales="free", labeller=mf_labeller) +
  scale_color_manual(values = cbPalette, name="Treatment") + 
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") + 
  labs(list(x="Treatment", y = "Mean distance relative to the bead (µm)"))+ 
  theme(axis.text=text, axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))











