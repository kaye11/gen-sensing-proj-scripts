NM2 <- read.csv("d:/Karen's/PhD/R program/General sensing proj/csv files/choice/notinduced_DPR.csv", sep=";")

NM <- subset (NM2, NM2$T>401 & NM2$T<601, )
#NM <- NM2

#This computes and saves track parameters from Taylor's equations, compelted data set, RMS data set.
#Rename data file into NM

library(ggplot2)
library(multcomp)
library(nlme)
library(lme4)
library(gdata)
library(grid)
library(plyr)

"%!in%" <- function(x, y) !(x %in% y)

# This makes a unique track identifier
NM$ID = as.factor(NM$ID)

#subset data into time points to save convergence failures
#SD <- subset (NM, NM$T<410, ) 
#SD <- subset (NM, NM$T>120 & NM$T<241)
#SD <- subset (NM, NM$T>240 & NM$T<361)
#SD <- subset (NM, NM$T>360 & NM$T<481)
#SD <- subset (NM, NM$T>480 & NM$T<601)

#NM=SD

# Total Number of tracks
lengthtracks=length(unique(NM$ID))

# Mean swimming speeds
Vmean=mean(NM$V, na.rm = TRUE)
NM$ND2=NM$ND*NM$ND

# make the mean ND2 for each timestep over all tracks
Rdata3 = aggregate(NM$ND2, by = list(time = NM$T), mean, na.action = na.omit)
Rdata3$RMS = sqrt(Rdata3$x)


# Taylor equation, fitting the curve
m4 = nls(RMS ~ ((v^2) * tau * 2 * (time - tau * (1 -
                                                   exp(-(time/tau)))))^0.5,
         data = Rdata3, start = list(v = 5, tau = 1), lower = 0.01,
         algorithm = "port", nls.control(maxiter=1000000))
coef(m4)


##Computing decorrelation length scale
decordistance = coef(m4)[1] * coef(m4)[2]
v = coef(m4)[1]
tau = coef(m4)[2]

##Computing diffusivity (D) (divisor is equal to 2 if video is 2d, 3 if video is 3d)
D = v^2*tau/2

#compute beta kernel 
rad<-as.numeric(readline("What is the radius of the bead?"))
B = 4*pi*rad*D

Rdata3$MF = (((v^2) * tau * 2) * (Rdata3$time - tau * (1 -exp(-(Rdata3$time/tau)))))^0.5

## Export calculated parameters
TA2<- readline("What data did you analyse? This is the summary for trackanalysis data!")
TA<-cbind(TA2, lengthtracks, Vmean, v, tau, decordistance, D, B)
write.table(TA, "d:/Karen's/PhD/R program/General sensing proj/Processed_data/trackanalysis.csv", 
            sep=";", col.names=F ,row.names = F, append=T)

##Export Rdata3
RM2<- readline("What RMS data is this?")
RM3<- paste("d:/Karen's/PhD/R program/General sensing proj/Processed_data/RMS data_choice/",RM2,".csv")
write.table(Rdata3, RM3, sep=";", row.names = F)

## Summary
print(c(lengthtracks, Vmean, v, tau, decordistance, D, B))

#Simple RMS plotting with the axis titles changed
Rdata3$MF = (((v^2) * tau * 2) * (Rdata3$time - tau * (1 -exp(-(Rdata3$time/tau)))))^0.5
RMS.plot=qplot(time, MF, data = Rdata3, na.rm=T)
RMS.plot + labs(list(x = "Time (s)", y = "RMS (µm)"))


# RMS distance travelled, computation and plotting RMS with t, regression line computed and 
#overlayed with data points
Rdata3$MF = (((v^2) * tau * 2) * (Rdata3$time - tau * (1 -exp(-(Rdata3$time/tau)))))^0.5
RMS.plot=qplot(time, MF, data = Rdata3, na.rm=T) + geom_point (shape=21, size=3, color="black")
newtime<- seq(from=min(Rdata3[,1]), to=max(Rdata3[,1]), length=length(Rdata3[,1])) 
reglineY <- predict(m4,list(time = newtime)) 

plot=RMS.plot + labs(list(x = "Time (s)", y = "RMS (µm)")) + 
  geom_line(aes(x=newtime, y= reglineY), colour = "red", lwd=2)


