#read in data

#Pbead
Pbead.rawtrack.simp <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/Pbead track data raw_simplified.csv", 
                          sep=";")

library(MotilityLab)
library(pracma)
#make data frame as list

Pbead.rawtrack.simp.list <- as.tracks(split (Pbead.rawtrack.simp, Pbead.rawtrack.simp$ID)) #splitting and make S3 obj as class track

length(unique(Pbead.rawtrack.simp$ID)) #45 IDs

Pbead.x.coords<- lapply(Pbead.rawtrack.simp.list,"[", , 3)
Pbead.y.coords<- lapply(Pbead.rawtrack.simp.list,"[", , 4)

Pbead.x.coords.results <- do.call(rbind.data.frame, lapply(Pbead.x.coords, hurstexp, display=FALSE)) #use pracma package and make result as dataframe
Pbead.y.coords.results <- do.call(rbind.data.frame, lapply(Pbead.y.coords, hurstexp, display=FALSE))

Pbead.x.coords.results$treatment <- "Pbead.x"
Pbead.y.coords.results$treatment <- "Pbead.y"

#control bead

control.rawtrack.simp <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/control track data raw_simplified.csv", 
                                  sep=";")

control.rawtrack.simp.list <- as.tracks(split (control.rawtrack.simp, control.rawtrack.simp$ID)) 

length(unique(control.rawtrack.simp$ID)) #45 IDs

control.x.coords<- lapply(control.rawtrack.simp.list,"[", , 3)
control.y.coords<- lapply(control.rawtrack.simp.list,"[", , 4)

control.x.coords.results <- do.call(rbind.data.frame, lapply(control.x.coords, hurstexp, display=FALSE))
control.y.coords.results <- do.call(rbind.data.frame, lapply(control.y.coords, hurstexp, display=FALSE))
control.x.coords.results$treatment <- "control.x"
control.y.coords.results$treatment <- "control.y"

hurstdata <- rbind(Pbead.x.coords.results, Pbead.y.coords.results, control.x.coords.results, control.y.coords.results)

library(tibble)
hurstdata <-rownames_to_column(hurstdata, var="ID")

library(reshape)
hurstdata_long <- melt(data=hurstdata, id.var=c("treatment", "ID"), 
                  measure.vars=c("Hs", "Hrs", "He", "Hal", "Ht"),
                  variable.name="HurstExp")

names(hurstdata_long)[names(hurstdata_long)=="value"] <- "HurstExp"
names(hurstdata_long)[names(hurstdata_long)=="variable"] <- "HurstExpType"

source("summarySE.R")

sum.hurst <- summarySE(hurstdata_long, measurevar="HurstExp", groupvars=c("treatment", "HurstExpType"))
#no difference on any HurstExpType when in control or P

##plotting
plot(Pbead.rawtrack.simp.list, dims = c("X", "Y"), add = F,
     pch.start = 1, pch.end = 2, cex = 0.5, xlab="X", ylab="Y", main="dP bead")

plot(control.rawtrack.simp.list, dims = c("X", "Y"), add = F,
     pch.start = 1, pch.end = 2, cex = 0.5, xlab="X", ylab="Y", main="control bead")

##measuring other random things

lapply(control.rawtrack.simp.list, straightness) #divides displacement by tracklength, straight=1 max

lapply(control.rawtrack.simp.list, asphericity) #like straightness but ignores back and forth movement of cells, aspherical=1 max