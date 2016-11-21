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

trackdata <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/choice/final track data with ND.csv", sep=";")

count_tracks <- ddply(trackdata, ~T, summarise, ID=length(unique(ID)))

NT=as.data.table(trackdata, key=c("bin"))

count_tracks <- NT[, list(count= length(unique(ID))), by = c("treatment", "bin", "T")]

count_tracks <- count_tracks [!count_tracks$bin=="bin_out", ]
count_tracks <- count_tracks [!count_tracks$bin=="bin_both", ]

track <- ggplot(count_tracks, aes(x=T, y=count, color=treatment)) + geom_line() + 
  labs(list(x = "Time (s)", y = "Number of Tracks"))  + facet_grid(~bin)

induced_DPR <- subset (trackdata, trackdata$treatment=="Si_induced" & trackdata$bin=="bin_DPR")
notinduced_DPR <- subset (trackdata, trackdata$treatment=="Si_notinduced" & trackdata$bin=="bin_DPR")

induced_dSi <- subset (trackdata, trackdata$treatment=="Si_induced" & trackdata$bin=="bin_dSi")
notinduced_dSi <- subset (trackdata, trackdata$treatment=="Si_notinduced" & trackdata$bin=="bin_dSi")

write.table (induced_DPR, "d:/Karen's/PhD/R program/General sensing proj/csv files/choice/induced_DPR.csv", 
             sep=";", col.names=T, row.names=F)

write.table (notinduced_DPR, "d:/Karen's/PhD/R program/General sensing proj/csv files/choice/notinduced_DPR.csv", 
             sep=";", col.names=T, row.names=F)

write.table (induced_dSi, "d:/Karen's/PhD/R program/General sensing proj/csv files/choice/induced_dSi.csv", 
             sep=";", col.names=T, row.names=F)

write.table (notinduced_dSi, "d:/Karen's/PhD/R program/General sensing proj/csv files/choice/notinduced_dSi.csv", 
             sep=";", col.names=T, row.names=F)
