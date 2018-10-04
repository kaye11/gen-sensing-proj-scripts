
#read data
#trackdata <- read.table ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/final track data with preprocessing.csv", sep=";", header=T)

trackdata <- read.table ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate automatic/final track data.csv", 
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

trackdata$timemin <- as.numeric(trackdata$time2)*120/60

source("summarySE.R")
speed.trackdata <- summarySE(trackdata,  measurevar="V", groupvars=c("ID"))

trackdata.slow <- speed.trackdata [speed.trackdata$V< 1, ]

trackdata_out <- trackdata[ ! trackdata$ID %in% trackdata.slow$ID, ] #deleting variables based on another dataframe 

speed.trackdata_out <- summarySE(trackdata_out,  measurevar="V", groupvars=c("ID"))

trackdata.recheck <- speed.trackdata_out [speed.trackdata_out$V<1, ]


library(plyr)
library(dplyr)
library(ggplot2)

#merging coordinates and data set (trackdata)
trackdata <- merge(trackdata, coordinates, by="wellvid")



#calculate distance between cell position and bead
trackdata$dist=sqrt(((trackdata$bead.X-trackdata$X)^2)+((trackdata$bead.Y-trackdata$Y)^2))

#calculate the sumdist for all. REMEMBER THIS IS YOUR DATA FOR THIS ANALYSIS. NOT SUMMARRY SE RESULTS BUT SUMDIST!!!
phosdist <- ddply(trackdata, c("timemin", "treatment", "bin"), summarise,
                 N    = length(dist),
                 mean = mean(dist, na.rm=TRUE),
                 sumdist= sum(dist, na.rm=TRUE), 
                 sd   = sd(dist, na.rm=TRUE),
                 se   = sd / sqrt(N))

phosdist$distmm <- phosdist$sumdist/1000

#qplot(as.factor(timemin), distmm, color = treatment, data = phosdist)+ stat_smooth(aes(group=treatment))

qplot(as.factor(timemin), distmm, color = treatment, data = phosdist)+ stat_smooth(aes(group=treatment), method="lm")+facet_grid(bin~., scales="free")

write.table (phosdist, "d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate automatic/distance automatic tracking complete.csv", 
             sep=";", col.names=T, row.names=F)

