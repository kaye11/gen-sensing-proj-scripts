
#read data
trackdata <- read.table ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/final track data with preprocessing.csv", 
                         sep=";", header=T)

coordinates <- read.table ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/coordinates.csv", 
                         sep=";", header=T)

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

qplot(as.factor(timemin), distmm, color = treatment, data = phosdist)+ stat_smooth(aes(group=treatment))

qplot(as.factor(timemin), distmm, color = treatment, data = phosdist)+ stat_smooth(aes(group=treatment), method="lm")+facet_grid(newbin~., scales="free")

write.table (phosdist, "d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/distance preprocessed.csv", 
             sep=";", col.names=T, row.names=F)

