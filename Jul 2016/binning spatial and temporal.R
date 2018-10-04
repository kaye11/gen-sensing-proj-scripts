##THIS CODE ONLY WORKS WITH DATA FRAMES AND NOT WITH DATA TABLES! binning temporally and spatially
#rename data files to t1, plots are saved manually

t1=as.data.frame(t1)

##BINNING
library(plotrix)
library(MASS)
library(data.table)
library(fBasics)
library(TeachingDemos)
library(data.table)

#binning temporally (CHANGE THIS WHEN NEEDED!)

tm=seq(0, 3600, by = 60)
t1$time <- cut(t1$T, tm, labels=paste(tail(tm, -1L)))

t1$time = factor(t1$time, levels=c(levels(t1$time), 0))


t1$time[is.na(t1$time)] <- 0 #replace NAs with 0. some data time points have the starting point as NA because they started at 0

##binning spatially
##What is the xy position of the bead in this video
angStr <- readline("Enter xy position of bead: ");
a <- as.numeric(unlist(strsplit(angStr, ",")));

#define center and radii of the circles
a
center <- a
radii <- seq(112, 336, 112)


# calculate the distance to the center for each row of object t1
distcenter <- apply(t1[, c("X", "Y")], 1, function(rowi) dist(rbind(rowi, center)))
# assign each row of object t1 to a group based on distcenter
group <- cut(distcenter, c(-1, radii, 2*max(radii)), labels=FALSE)


# draw circles in x,y coordinate space, use the eqscplot() function from package MASS
# do not rerun this part of the script, otherwise the points will be messed up
eqscplot(0, 0, type="n", xlim=range(t1$X, center[1]+c(-1, 1)*max(radii)), ylim=rev(range(t1$Y, center[2]+c(-1, 1)*max(radii))))
draw.circle(center[1], center[2], 112)
draw.circle(center[1], center[2], 224)
draw.circle(center[1], center[2], 336)
points(t1$X, t1$Y, col=group)
#save the plots manually



# to subset into bins
binA=t1[group==1, ]
binB=t1[group==2, ]
binC=t1[group==3, ]
outbin=t1[group==4, ]

#t1
t1binA=binA
t1binB=binB
t1binC=binC
t1out=outbin

t1binA$bin=c("binA")
t1binB$bin=c("binB")
t1binC$bin=c("binC")
t1out$bin=c("outbin")

t1=rbind (t1binA, t1binB, t1binC, t1out)

VN<- readline("What data did you analyse? SPECIFY DATA TYPE:")
Vid<-paste ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate automatic/binned track data/",VN,".csv")
write.table(t1, Vid, sep=";", col.names=T, row.names=F)

