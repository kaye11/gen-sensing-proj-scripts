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

#for max number of tracks with the start of each track forced to 0
t1=trackdata
NT = data.table(t1)
NT2=NT[, T := seq(from = 1L, by = 1L, length.out = .N), by = ID]

trackdata.forced <- NT2

#save
write.table (trackdata.forced, "d:/Karen's/PhD/R program/General sensing proj/csv files/choice/trackdata_forced.csv", 
             sep=";", col.names=T, row.names=F)

#set number of tracks
NT=as.data.table(trackdata.forced, key=c("bin"))

count_tracks <- NT[, list(count= length(unique(ID))), by = c("treatment", "bin", "T")]

count_tracks <- count_tracks [!count_tracks$bin=="bin_out", ]
count_tracks <- count_tracks [!count_tracks$bin=="bin_both", ]

track <- ggplot(count_tracks, aes(x=T, y=count, color=treatment)) + geom_line() + 
  labs(list(x = "Time (s)", y = "Number of Tracks"))  + facet_grid(~bin)


induced_DPR <- subset (trackdata.forced, trackdata.forced$treatment=="Si_induced" & trackdata.forced$bin=="bin_DPR")
notinduced_DPR <- subset (trackdata.forced, trackdata.forced$treatment=="Si_notinduced" & trackdata.forced$bin=="bin_DPR")

induced_dSi <- subset (trackdata.forced, trackdata.forced$treatment=="Si_induced" & trackdata.forced$bin=="bin_dSi")
notinduced_dSi <- subset (trackdata.forced, trackdata.forced$treatment=="Si_notinduced" & trackdata.forced$bin=="bin_dSi")

induced_track <- subset (trackdata.forced, trackdata.forced$treatment=="Si_induced")
notinduced_track <- subset (trackdata.forced, trackdata.forced$treatment=="Si_notinduced")

write.table (induced_DPR, "d:/Karen's/PhD/R program/General sensing proj/csv files/choice/induced_DPR_forced.csv", 
             sep=";", col.names=T, row.names=F)

write.table (notinduced_DPR, "d:/Karen's/PhD/R program/General sensing proj/csv files/choice/notinduced_DPR_forced.csv", 
             sep=";", col.names=T, row.names=F)

write.table (induced_dSi, "d:/Karen's/PhD/R program/General sensing proj/csv files/choice/induced_dSi_forced.csv", 
             sep=";", col.names=T, row.names=F)

write.table (notinduced_dSi, "d:/Karen's/PhD/R program/General sensing proj/csv files/choice/notinduced_dSi_forced.csv", 
             sep=";", col.names=T, row.names=F)

write.table (induced_track, "d:/Karen's/PhD/R program/General sensing proj/csv files/choice/induced_track_forced.csv", 
             sep=";", col.names=T, row.names=F)

write.table (notinduced_track, "d:/Karen's/PhD/R program/General sensing proj/csv files/choice/notinduced_track_forced.csv", 
             sep=";", col.names=T, row.names=F)

