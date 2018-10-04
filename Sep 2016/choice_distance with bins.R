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

dSi.dist <- aggregate(trackdata2$dist_dSi, by = list(ID = trackdata2$ID, treatment=trackdata2$treatment), last)
DPR.dist <- aggregate(trackdata2$dist_DPR, by = list(ID = trackdata2$ID, treatment=trackdata2$treatment), last)

dSi.dist$bead = as.factor ("dSibead")
DPR.dist$bead = as.factor("DPRbead")

trackdata.dist <- rbind (dSi.dist, DPR.dist)

trackdata.dist$dist = trackdata.dist$x

dist.sum <- ddply(trackdata.dist, c("treatment", "bead"), summarise,
                  N    = length(ID),
                  mean = mean(dist, na.rm=TRUE),
                  sumdist= sum(dist, na.rm=TRUE), 
                  sd   = sd(dist, na.rm=TRUE),
                  se   = sd / sqrt(N))

ggplot(data=dist.sum, aes(x=bead, y=mean)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1) + facet_grid(treatment~., scales="free")

