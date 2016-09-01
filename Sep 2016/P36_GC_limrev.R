library(ggplot2)
library(grid)
library(gtable)
library(ggthemes)
library(gridExtra)
library(mgcv)
library(data.table)
library(nlme)


dev.new(width=6, height=9)
source("resizewin.R")
resize.win(12,9)


#read in data

P36_GC_limrev <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/growth curves/P36_limitation and recovery_GC_longformat.csv", sep=";")

#compare log normal in R and excel
#lncells is from excel

P36_GC_limrev$lncellsR <- log(P36_GC_limrev$cellsml)

#almost the same, only decimal places are different

#make new factors

P36_GC_limrev$medtreat <- as.factor (paste(P36_GC_limrev$med, P36_GC_limrev$treatment, sep="-"))
P36_GC_limrev$medtreatreplet <- as.factor (paste(P36_GC_limrev$medtreat, P36_GC_limrev$replet, sep="-"))

#plotting
qplot(day, cellsml, data = P36_GC_limrev,  geom = "boxplot") + facet_grid(medtreat~., scales="free") 
qplot(day, lncellsR, data = P36_GC_limrev,  geom = "boxplot") + facet_grid(medtreat~., scales="free") 

source("summarySE.R")
GC.sum <- summarySE(P36_GC_limrev, measurevar="lncellsR", groupvars=c("medtreat", "day"))

#plot with log normal cells R
ggplot(data=GC.sum, aes(x=day, y=lncellsR)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=lncellsR-se, ymax=lncellsR+se), width=0.5, size=1) + facet_wrap(~medtreat, nrow=2)+
  geom_smooth(method="loess", aes(group=1))

ggplot(data=GC.sum, aes(x=day, y=lncellsR)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=lncellsR-se, ymax=lncellsR+se), width=0.5, size=1) + facet_wrap(~medtreat, nrow=2)+
  geom_smooth(method="lm", aes(group=1))

GC.sum.raw <- summarySE(P36_GC_limrev, measurevar="cellsml", groupvars=c("medtreat", "day"))

ggplot(data=GC.sum.raw, aes(x=day, y=cellsml)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsml-se, ymax=cellsml+se), width=0.2, size=1) + facet_wrap(~medtreat, nrow=2)+
  geom_smooth(method="loess", aes(group=medtreat))

ggplot(data=GC.sum.raw, aes(x=day, y=cellsml)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsml-se, ymax=cellsml+se), width=0.2, size=1) + facet_wrap(~medtreat, nrow=2)+
  geom_smooth(method="lm", aes(group=medtreat))


#compute some parameters
#Two points, N1 and N2,  at the extremes of this linear phase (see fig below) are taken and substituted into the equation
# Growth rate ;  K' = Ln (N2 / N1) / (t2 - t1)
# Divisions per day ; Div.day-1 = K' / Ln2
# Generation time ; Gen' t  = 1 / Div.day-1

#check days first

#ASW
ASW = data.table(subset (P36_GC_limrev, P36_GC_limrev$medtreat=="ASW-control" ), key=c("medtreatreplet"))

ASW.NT = ASW [ASW$day %in% c("day0", "day3"),  ] #or ASW.NT = subset(ASW, ASW$day=="day1" | ASW$day=="day3", )

ASW.comp=ASW.NT [,list(K= (lncellsR[match("day3", day)]-lncellsR [match ("day0", day)])/3),   #divide according to number of days
                 by=c("medtreatreplet")]

ASW.comp [,div.day:=K/log(2)] #this is always 2
ASW.comp [, genT:= 1/div.day]
ASW.comp [, medtreat:="ASW-control"]

#NO3.lim

NO3.lim = data.table(subset (P36_GC_limrev, P36_GC_limrev$medtreat=="-NO3-limitation" ), key=c("medtreatreplet"))

NO3.lim.NT = NO3.lim [NO3.lim$day %in% c("day0", "day2"),  ] 

NO3.lim.comp=NO3.lim.NT [,list(K= (lncellsR[match("day2", day)]-lncellsR [match ("day0", day)])/2),   #divide according to number of days
                         by=c("medtreatreplet")]

NO3.lim.comp [,div.day:=K/log(2)] #this is always 2
NO3.lim.comp [, genT:= 1/div.day]
NO3.lim.comp  [, medtreat:="-NO3-limitation"]

#NO3.rec

NO3.rec = data.table(subset (P36_GC_limrev, P36_GC_limrev$medtreat=="-NO3-recovery" ), key=c("medtreatreplet"))

NO3.rec.NT = NO3.rec [NO3.rec$day %in% c("day1", "day5"),  ] 

NO3.rec.comp=NO3.rec.NT [,list(K= (lncellsR[match("day5", day)]-lncellsR [match ("day1", day)])/4),   #divide according to number of days
                         by=c("medtreatreplet")]

NO3.rec.comp [,div.day:=K/log(2)] #this is always 2
NO3.rec.comp [, genT:= 1/div.day]
NO3.rec.comp  [, medtreat:="-NO3-recovery"]


#PO4.lim

PO4.lim = data.table(subset (P36_GC_limrev, P36_GC_limrev$medtreat=="-PO4-limitation" ), key=c("medtreatreplet"))

PO4.lim.NT = PO4.lim [PO4.lim$day %in% c("day1", "day3"),  ] 

PO4.lim.comp=PO4.lim.NT [,list(K= (lncellsR[match("day3", day)]-lncellsR [match ("day1", day)])/2),   #divide according to number of days
                         by=c("medtreatreplet")]

PO4.lim.comp [,div.day:=K/log(2)] #this is always 2
PO4.lim.comp [, genT:= 1/div.day]
PO4.lim.comp  [, medtreat:="-PO4-limitation"]


#PO4.rec

PO4.rec = data.table(subset (P36_GC_limrev, P36_GC_limrev$medtreat=="-PO4-recovery" ), key=c("medtreatreplet"))

PO4.rec.NT = PO4.rec [PO4.rec$day %in% c("day0", "day3"),  ] 

PO4.rec.comp=PO4.rec.NT [,list(K= (lncellsR[match("day3", day)]-lncellsR [match ("day0", day)])/3),   #divide according to number of days
                         by=c("medtreatreplet")]

PO4.rec.comp [,div.day:=K/log(2)] #this is always 2
PO4.rec.comp [, genT:= 1/div.day]
PO4.rec.comp  [, medtreat:="-PO4-recovery"]

#Si.lim

Si.lim = data.table(subset (P36_GC_limrev, P36_GC_limrev$medtreat=="-Si-limitation" ), key=c("medtreatreplet"))

Si.lim.NT = Si.lim [Si.lim$day %in% c("day0", "day2"),  ] 

Si.lim.comp=Si.lim.NT [,list(K= (lncellsR[match("day2", day)]-lncellsR [match ("day0", day)])/2),   #divide according to number of days
                       by=c("medtreatreplet")]

Si.lim.comp [,div.day:=K/log(2)] #this is always 2
Si.lim.comp [, genT:= 1/div.day]
Si.lim.comp  [, medtreat:="-Si-limitation"]

#Si.rec

Si.rec = data.table(subset (P36_GC_limrev, P36_GC_limrev$medtreat=="-Si-recovery" ), key=c("medtreatreplet"))

Si.rec.NT = Si.rec [Si.rec$day %in% c("day0", "day4"),  ] 

Si.rec.comp=Si.rec.NT [,list(K= (lncellsR[match("day4", day)]-lncellsR [match ("day0", day)])/4),   #divide according to number of days
                       by=c("medtreatreplet")]

Si.rec.comp [,div.day:=K/log(2)] #this is always 2
Si.rec.comp [, genT:= 1/div.day]
Si.rec.comp  [, medtreat:="-Si-recovery"]


#bind all

GC.comp <- rbind(ASW.comp, NO3.lim.comp, NO3.rec.comp, PO4.lim.comp, PO4.rec.comp, Si.lim.comp, Si.rec.comp)

#summaries

GC.comp.sumK <- summarySE(data=GC.comp, measurevar="K", groupvars="medtreat", na.rm=TRUE)


