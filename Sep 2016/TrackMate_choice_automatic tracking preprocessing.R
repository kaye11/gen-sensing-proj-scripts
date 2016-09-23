#data set is split into induced and not induced

#STEP 1
dirdata<-"d:/Karen's/PhD/R program/General sensing proj/csv files/choice/raw track data/"


#get dir/file names
ff<-do.call(rbind, lapply(dirdata, function(x) {
  ff<-list.files(x, "\\.xls$", include.dirs = FALSE, full.names = TRUE)
  data.frame(dir=basename(x), file=basename(ff), 
             fullpath=ff, stringsAsFactors=F)
}))

#read into data.frame

source("readstack.R") #run the code readstack.R
library(plyr)
merged.data= read.stack(ff$fullpath, extra=list(file=ff$file))

library(plyr)

#read into data.frame

pattern <- "(\\w{1})(\\d{1})(_)(\\w+)(_)(\\w+)(-)(\\d+)(.xls)"

well <- gsub(pattern,'\\1\\2',merged.data$file)
wellvid <- gsub(pattern,'\\1\\2\\7\\8',merged.data$file)
treatment <- gsub(pattern,'\\4\\5\\6',merged.data$file)

final.data=cbind(merged.data, well, wellvid, treatment)

#make new ID
final.data$ID <- as.factor(paste(final.data$wellvid, final.data$TRACK_ID, sep="-"))

##data table challenge
library (data.table)

##direct input from trackmate 

NT = data.table(final.data, key="ID")
trackdata = NT[, list(treatment=treatment, wellvid=wellvid, A=TRACK_ID, X=POSITION_X, Y=POSITION_Y, T=FRAME), by=c("ID")]

trackdata [, V:= c(0, sqrt(diff(X)^2+diff(Y)^2)), by=c("ID")] 
trackdata [, Vlog:= log(V+1)] 
trackdata [, ND:=sqrt((X-X[1])^2 + (Y-Y[1])^2), by=c("ID")]


#merging coordinates

bead_coords <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/bead coordinates choice.csv", sep=";")

trackdatacompvar <- merge(trackdata, bead_coords, by="wellvid")

#binning spatially
source ("circle.R")

trackdatacompvar <- within(trackdatacompvar,
                           { 
                             bin_dSi <- check_if_in_circle(points = cbind(X, Y), 
                                                           x_dSi, y_dSi, rad_dSi + 115)
                             bin_DPR <- check_if_in_circle(points = cbind(X, Y), 
                                                           x_DPR, y_DPR, rad_DPR + 115)
                             bin <- ifelse(bin_dSi, ifelse(bin_DPR, 'bin_both', 'bin_dSi'),
                                           ifelse(bin_DPR, 'bin_DPR', 'bin_out'))
                           } )

trackdatacompvar$bin <- as.factor(trackdatacompvar$bin)


#binning temporally

tm=seq(0, 600, by = 60)

trackdatacompvar$time <- cut(trackdatacompvar$T, tm, labels=paste(tail(tm, -1L)))

trackdatacompvar$time = factor(trackdatacompvar$time, levels=c(levels(trackdatacompvar$time), 0))

trackdatacompvar$time[is.na(trackdatacompvar$time)] <- 0 #replace NAs with 0. some data time points have the starting point as NA because they started at 0

trackdatacompvar$time <- relevel (trackdatacompvar$time, "0")

#plotting
library(ggplot2)

qplot(time, Vlog, color = treatment, data = trackdatacompvar,  geom = "boxplot") + facet_grid(treatment~bin, scales="free") 

#drop time=0

final.track.data <- trackdatacompvar[trackdatacompvar$time!= "0",]

qplot(time, Vlog, color = treatment, data = final.track.data,  geom = "boxplot") + facet_grid(treatment~bin, scales="free") 


#saving 

write.table (final.track.data, "d:/Karen's/PhD/R program/General sensing proj/csv files/choice/final track data with ND.csv", 
             sep=";", col.names=T, row.names=F)

#split data points to induced and not induced

induced.track.data <- subset (final.track.data, final.track.data$treatment=="Si_induced")

notinduced.track.data <- subset (final.track.data, final.track.data$treatment=="Si_notinduced")

write.table (induced.track.data, "d:/Karen's/PhD/R program/General sensing proj/csv files/choice/induced track data.csv", 
             sep=";", col.names=T, row.names=F)

write.table (notinduced.track.data, "d:/Karen's/PhD/R program/General sensing proj/csv files/choice/notinduced track data.csv", 
             sep=";", col.names=T, row.names=F)