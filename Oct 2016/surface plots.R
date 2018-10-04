trackdata2 <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/choice/induced track data.csv", sep=";")

trackdata <- subset(trackdata2, trackdata2$wellvid=="C4-036", )

##from jan
library(plyr)

tracksum <- ddply(trackdata, .(X, Y), summarize, meanV = mean(V, na.rm = TRUE),medianV = median(V, na.rm = TRUE))

library(ggplot2)

ggplot(data = tracksum, aes(X,Y,z =meanV)) + stat_summary_2d()
ggplot(data = tracksum, aes(X,Y,z =medianV)) + stat_summary_2d()

##karen version
ggplot(data=trackdata, aes(X, Y, z=V))+stat_summary_2d(mapping=aes(x=X, y=Y, z=V), 
                                                data=trackdata, fun = mean)


library(spatial)
library(akima)

##medianV
topo.meter.ls3 <- surf.ls(6, tracksum$X,tracksum$Y,tracksum$medianV)
summary(tracksum)
topo.meter.surface3 <- trmat(topo.meter.ls3, 0, max(tracksum$X), 0, max(tracksum$Y), 100)
image(topo.meter.surface3)
contour(topo.meter.surface3, add = TRUE)

##meanV
topo.meter.ls3 <- surf.ls(6, tracksum$X,tracksum$Y,log(tracksum$meanV+1))
summary(tracksum)
topo.meter.surface3 <- trmat(topo.meter.ls3, 0, max(tracksum$X), 0, max(tracksum$Y), 50)
image(topo.meter.surface3)
contour(topo.meter.surface3, add = TRUE)

filled.contour(topo.meter.surface3, col=terrain.colors(25))


library(fields)
##akima
x = tracksum$X
y = tracksum$Y
z = tracksum$meanV

akima <- interp(x, y, z)
image.plot  (akima)
contour(akima.li, add=TRUE)

tracksum2 <- ddply(silab1, .(X, Y), summarize, meanP = mean(P, na.rm = TRUE), medianP = median(P, na.rm = TRUE))
tracksum2[is.na(tracksum2)]<-0

topo.meter.ls3 <- surf.ls(6, tracksum2$X,tracksum2$Y,tracksum2$medianV)
summary(tracksum2)
topo.meter.surface4 <- trmat(topo.meter.ls3, 0, max(tracksum2$X), 0, max(tracksum2$Y), 50)
image(topo.meter.surface4)
contour(topo.meter.surface4, add = TRUE)


