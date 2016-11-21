library(gridExtra)
library(ggplot2)
library(sp)
library(maptools)
library(gstat)
library(plyr)

trackdata2 <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/choice/notinduced track data.csv", sep=";")

trackdata <- subset(trackdata2, trackdata2$wellvid=="B4-015", )

##from jan
library(plyr)

tracksum <- ddply(trackdata, .(X, Y), summarize, meanV = mean(V, na.rm = TRUE),medianV = median(V, na.rm = TRUE))

# convert this basic data frame into a spatial points data frame
coordinates(tracksum) = ~ X + Y

plot(tracksum)

x.range <- as.integer(range(tracksum@coords[, 1]))
y.range <- as.integer(range(tracksum@coords[, 2]))
plot(tracksum)

grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 20), y = seq(from = y.range[1], to = y.range[2], by = 50))
coordinates(grd) <- ~x + y
gridded(grd) <- TRUE

plot(grd, cex = 1.5)
points(tracksum, pch = 1, col = "red", cex = 1)
title("B4-015 notinduced")


#Kriging

semivariog <- variogram(log(meanV+1) ~ 1, tracksum)
plot(semivariog)
semivariog

autovariogram = autofitVariogram(log(meanV+1) ~ 1, tracksum)
plot(autovariogram)

model.variog <- vgm(psill = 4.7, model = "Gau", nugget = 0, range = 5.5)
fit.variog <- fit.variogram(semivariog, model.variog)
plot(semivariog, fit.variog)
bubble(tracksum, zcol='meanV', fill=FALSE, do.sqrt=FALSE, maxsize=2)


x.range <- as.integer(range(tracksum@coords[, 1]))
y.range <- as.integer(range(tracksum@coords[, 2]))
plot(tracksum)


grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 50), y = seq(from = y.range[1], to = y.range[2], by = 50))
coordinates(grd) <- ~x + y
gridded(grd) <- TRUE

krig <- krige(formula = meanV ~ 1, locations = tracksum, newdata = grd, model = model.variog, nmax = 1000)

krig.output = as.data.frame(krig)
names(krig.output)[1:3] <- c("X", "Y", "var1.pred")

plot <- ggplot(data = krig.output, aes(x = X, y = Y))  #start with the base-plot and add the Kriged data to it
layer1 <- c(geom_tile(data = krig.output, aes(fill = var1.pred)))  #then create a tile layer and fill with predicted
KrigControl = plot + layer1  + scale_fill_gradient(low = "yellow", high = "red") + 
  labs(title = "C3-006-induced") + coord_fixed()+  scale_y_reverse()+
  annotate("path",x = 481.2 + 16*cos(seq(0,2*pi,length.out=100)), y=539.6 + 16*sin(seq(0,2*pi,length.out=100))) +
  annotate("path",x = 831.2 + 55*cos(seq(0,2*pi,length.out=100)), y=585.4 + 55*sin(seq(0,2*pi,length.out=100))) 

KrigControl


library(automap)

# Ordinary kriging
kriging_result = autoKrige(meanV~1, tracksum, grd, nmax=100)
plot(kriging_result)
# Universal kriging
kriging_result = autoKrige(zinc~soil+ffreq+dist, meuse, meuse.grid)
plot(kriging_result)
