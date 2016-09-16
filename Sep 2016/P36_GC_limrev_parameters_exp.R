#check days first

#ASW
ASW = data.table(subset (P36_GC_limrev, P36_GC_limrev$medtreat=="ASW-control" ), key=c("medtreatreplet"))

ASW.NT = ASW [ASW$day %in% c("day1", "day2"),  ] #or ASW.NT = subset(ASW, ASW$day=="day1" | ASW$day=="day3", )

ASW.comp=ASW.NT [,list(K= (lncellsR[match("day2", day)]/lncellsR [match ("day1", day)])/1),   #divide according to number of days
                 by=c("medtreatreplet")]

ASW.comp [,div.day:=K/log(2)] #this is always 2
ASW.comp [, genT:= 1/div.day]
ASW.comp [, medtreat:="ASW-control"]

#NO3.lim

NO3.lim = data.table(subset (P36_GC_limrev, P36_GC_limrev$medtreat=="-NO3-limitation" ), key=c("medtreatreplet"))

NO3.lim.NT = NO3.lim [NO3.lim$day %in% c("day1", "day2"),  ] 

NO3.lim.comp=NO3.lim.NT [,list(K= (lncellsR[match("day2", day)]/lncellsR [match ("day1", day)])/1),   #divide according to number of days
                         by=c("medtreatreplet")]

NO3.lim.comp [,div.day:=K/log(2)] #this is always 2
NO3.lim.comp [, genT:= 1/div.day]
NO3.lim.comp  [, medtreat:="-NO3-limitation"]

#NO3.rec

NO3.rec = data.table(subset (P36_GC_limrev, P36_GC_limrev$medtreat=="-NO3-recovery" ), key=c("medtreatreplet"))

NO3.rec.NT = NO3.rec [NO3.rec$day %in% c("day1", "day5"),  ] 

NO3.rec.comp=NO3.rec.NT [,list(K= (lncellsR[match("day5", day)]/lncellsR [match ("day1", day)])/5),   #divide according to number of days
                         by=c("medtreatreplet")]

NO3.rec.comp [,div.day:=K/log(2)] #this is always 2
NO3.rec.comp [, genT:= 1/div.day]
NO3.rec.comp  [, medtreat:="-NO3-recovery"]


#PO4.lim

PO4.lim = data.table(subset (P36_GC_limrev, P36_GC_limrev$medtreat=="-PO4-limitation" ), key=c("medtreatreplet"))

PO4.lim.NT = PO4.lim [PO4.lim$day %in% c("day1", "day3"),  ] 

PO4.lim.comp=PO4.lim.NT [,list(K= (lncellsR[match("day3", day)]/lncellsR [match ("day1", day)])/2),   #divide according to number of days
                         by=c("medtreatreplet")]

PO4.lim.comp [,div.day:=K/log(2)] #this is always 2
PO4.lim.comp [, genT:= 1/div.day]
PO4.lim.comp  [, medtreat:="-PO4-limitation"]


#PO4.rec

PO4.rec = data.table(subset (P36_GC_limrev, P36_GC_limrev$medtreat=="-PO4-recovery" ), key=c("medtreatreplet"))

PO4.rec.NT = PO4.rec [PO4.rec$day %in% c("day1", "day2"),  ] 

PO4.rec.comp=PO4.rec.NT [,list(K= (lncellsR[match("day2", day)]/lncellsR [match ("day0", day)])/2),   #divide according to number of days
                         by=c("medtreatreplet")]

PO4.rec.comp [,div.day:=K/log(2)] #this is always 2
PO4.rec.comp [, genT:= 1/div.day]
PO4.rec.comp  [, medtreat:="-PO4-recovery"]

#Si.lim

Si.lim = data.table(subset (P36_GC_limrev, P36_GC_limrev$medtreat=="-Si-limitation" ), key=c("medtreatreplet"))

Si.lim.NT = Si.lim [Si.lim$day %in% c("day1", "day7"),  ] 

Si.lim.comp=Si.lim.NT [,list(K= (lncellsR[match("day7", day)]/lncellsR [match ("day1", day)])/6),   #divide according to number of days
                       by=c("medtreatreplet")]

Si.lim.comp [,div.day:=K/log(2)] #this is always 2
Si.lim.comp [, genT:= 1/div.day]
Si.lim.comp  [, medtreat:="-Si-limitation"]

#Si.rec

Si.rec = data.table(subset (P36_GC_limrev, P36_GC_limrev$medtreat=="-Si-recovery" ), key=c("medtreatreplet"))

Si.rec.NT = Si.rec [Si.rec$day %in% c("day0", "day4"),  ] 

Si.rec.comp=Si.rec.NT [,list(K= (lncellsR[match("day4", day)]/lncellsR [match ("day0", day)])/4),   #divide according to number of days
                       by=c("medtreatreplet")]

Si.rec.comp [,div.day:=K/log(2)] #this is always 2
Si.rec.comp [, genT:= 1/div.day]
Si.rec.comp  [, medtreat:="-Si-recovery"]


#bind all

GC.comp <- rbind(ASW.comp, NO3.lim.comp, NO3.rec.comp, PO4.lim.comp, PO4.rec.comp, Si.lim.comp, Si.rec.comp)

#summaries

GC.comp.sumK <- summarySE(data=GC.comp, measurevar="K", groupvars="medtreat", na.rm=TRUE)
GC.comp.sumdiv.day <- summarySE(data=GC.comp, measurevar="div.day", groupvars="medtreat", na.rm=TRUE)
GC.comp.sumgenT <- summarySE(data=GC.comp, measurevar="genT", groupvars="medtreat", na.rm=TRUE)

GC.comp.sum <- cbind(GC.comp.sumK, GC.comp.sumdiv.day, GC.comp.sumgenT)

GC.comp <- na.omit(GC.comp)
GC.comp$medtreat <- as.factor (GC.comp$medtreat)


Kplot= qplot(medtreat, K, data = GC.comp,  geom = "boxplot", na.rm=TRUE) 
divdayplot= qplot(medtreat, div.day, data = GC.comp,  geom = "boxplot", na.rm=TRUE) 
genTplot= qplot(medtreat, genT, data = GC.comp,  geom = "boxplot", na.rm=TRUE) 

source ("multiplot.R")

multiplot(Kplot, divdayplot, genTplot)


##K

#boxplots
op=par(mfrow=c(1,2))
boxplot(K~medtreat, data=GC.comp)
boxplot(K~medtreatreplet, data=GC.comp)

library(lawstat)
levene.test(GC.comp$K, group=GC.comp$medtreat, location="mean") #equal

shapiro.test(GC.comp$K) #not normal
bartlett.test(K~ medtreat, data=GC.comp) #not homogenous


#check normality with models
K.aov <- aov(K~medtreat, data=GC.comp)

op=par(mfrow=c(2,2))
plot(K.aov) #not normaaaal! revert to kruskal wallis

#kruskal
kruskal.test(K ~ medtreat, data=GC.comp)

library(PMCMR)
posthoc.kruskal.nemenyi.test(K ~ medtreat, data=GC.comp, dist="Tukey")

##div.day

#boxplots
op=par(mfrow=c(1,2))
boxplot(div.day~medtreat, data=GC.comp)
boxplot(div.day~medtreatreplet, data=GC.comp)

levene.test(GC.comp$div.day, group=GC.comp$medtreat, location="mean") #equal

shapiro.test(GC.comp$div.day) #not normal
bartlett.test(div.day~ medtreat, data=GC.comp) #not homogenous


#check div.day normality with models
div.day.aov <- aov(div.day~medtreat, data=GC.comp)

op=par(mfrow=c(2,2))
plot(div.day.aov) #not normaaaal! revert to kruskal wallis

#kruskal wallis
kruskal.test(div.day ~ medtreat, data=GC.comp)

posthoc.kruskal.nemenyi.test(div.day ~ medtreat, data=GC.comp, dist="Tukey")


##genT

#boxplots
op=par(mfrow=c(1,2))
boxplot(genT~medtreat, data=GC.comp)
boxplot(genT~medtreatreplet, data=GC.comp)

levene.test(GC.comp$genT, group=GC.comp$medtreat, location="mean") #equal

shapiro.test(GC.comp$genT) #not normal
bartlett.test(genT~ medtreat, data=GC.comp) #homogenous


#check normality with models
genT.aov <- aov(genT~medtreat, data=GC.comp)

op=par(mfrow=c(2,2))
plot(genT.aov) #not normaaaal! revert to kruskal wallis

#kruskal wallis
kruskal.test(genT ~ medtreat, data=GC.comp)

posthoc.kruskal.nemenyi.test(genT ~ medtreat, data=GC.comp, dist="Tukey")






