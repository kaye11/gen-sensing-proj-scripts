
#read data
#trackdata <- read.table ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/final track data with preprocessing.csv", sep=";", header=T)

trackdata <- read.table ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/final track data.csv", 
                         sep=";", header=T)

coordinates <- read.table ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/coordinates.csv", 
                         sep=";", header=T)


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

#trackdata <- trackdata [! trackdata$ID=="A3-004-0", ]
#trackdata <- trackdata [! trackdata$ID=="A3-004-469", ]
trackdata <- trackdata [! trackdata$ID=="A3-004-946", ] #up to here bins AB is okay, just delete this track
#trackdata <- trackdata [! trackdata$ID=="A3-004-6", ]

#trackdata <- trackdata [! trackdata$ID=="B4-026-739", ]

#trackdata <- trackdata [! trackdata$ID=="B1-032-119", ]
#trackdata <- trackdata [! trackdata$ID=="B1-032-532", ]
#trackdata <- trackdata [! trackdata$ID=="B1-032-34", ]

trackdata <- trackdata [! trackdata$ID=="B4-026-126", ]



trackdata$timemin <- as.numeric(trackdata$time2)*120/60

library(plyr)
library(dplyr)
library(ggplot2)

#merging coordinates and data set (trackdata)
trackdata <- merge(trackdata, coordinates, by="wellvid")

#calculate distance between cell position and bead
trackdata$dist=sqrt(((trackdata$bead.X-trackdata$X)^2)+((trackdata$bead.Y-trackdata$Y)^2))

#calculate the sumdist for all. REMEMBER THIS IS YOUR DATA FOR THIS ANALYSIS. NOT SUMMARRY SE RESULTS BUT SUMDIST!!!
phosdist <- ddply(trackdata, c("timemin", "treatment", "newbin"), summarise,
                 N    = length(dist),
                 mean = mean(dist, na.rm=TRUE),
                 sumdist= sum(dist, na.rm=TRUE), 
                 sd   = sd(dist, na.rm=TRUE),
                 se   = sd / sqrt(N))

phosdist$distmm <- phosdist$sumdist/1000

qplot(as.factor(timemin), distmm, color = treatment, data = phosdist)+ stat_smooth(aes(group=treatment), method="lm")+facet_grid(newbin~., scales="free")

#smoothing is different if you have collated data beforehand or if you ask ggplot to collate it for you. better use collated data for consistency
#not collated data is good for detecting outliers and tracks that has only few point contributions
#max data pts per ID per time = 60

phosdist.ID <- ddply(trackdata, c("timemin", "treatment", "newbin", "wellvid", "ID"), summarise,
                  N    = length(dist),
                  mean = mean(dist, na.rm=TRUE),
                  sumdist= sum(dist, na.rm=TRUE), 
                  sd   = sd(dist, na.rm=TRUE),
                  se   = sd / sqrt(N))

phosdist.ID$distmm <- phosdist.ID$sumdist/1000

qplot(as.factor(timemin), distmm, color = treatment, data = phosdist.ID,  geom = "boxplot") + facet_grid(treatment~newbin, scales="free") +
  stat_smooth(aes(group=treatment), method="lm")

qplot(as.factor(timemin), distmm, color = treatment, data = phosdist.ID)+ stat_smooth(aes(group=treatment), method="lm")+
  facet_grid(newbin~wellvid, scales="free")


qplot(as.factor(timemin), distmm, data = phosdist.ID [phosdist.ID$wellvid=="A1-031", ], color=ID)+ 
  stat_smooth(aes(group=ID), method="lm") + facet_grid(newbin~ID, scales="free")

#phosdist.ID.sub <- phosdist.ID [phosdist.ID$N > 9, ]

qplot(as.factor(timemin), distmm, color = treatment, data = phosdistsub)+ stat_smooth(aes(group=treatment), method="lm")+facet_grid(newbin~wellvid, scales="free")

qplot(as.factor(timemin), distmm, color = treatment, data = phosdistsub)+ stat_smooth(aes(group=treatment), method="lm")+facet_grid(~newbin, scales="free")

qplot(as.factor(timemin), distmm, color = treatment, data = phosdist)+ stat_smooth(aes(group=treatment), method="lm")+facet_grid(~newbin, scales="free")

phosdistout <- ddply(phosdist.ID.sub, c("timemin", "treatment", "newbin"), summarise,
                  N    = sum(N),
                  sumdist= sum(distmm, na.rm=TRUE), 
                  sd   = sd(distmm, na.rm=TRUE),
                  se   = sd / sqrt(N))

qplot(as.factor(timemin), sumdist, color = treatment, data = phosdistout)+ stat_smooth(aes(group=treatment), method="lm")+facet_grid(~newbin, scales="free")


write.table (phosdist, "d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/distance complete minus one track.csv", 
             sep=";", col.names=T, row.names=F)

