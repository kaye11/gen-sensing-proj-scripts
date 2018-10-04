library(ggplot2)
library(gtable)
library(ggthemes)
library(grLabelExtra)
library(mgcv)
library(ggplot2)
library(grid)
library(data.table)

source("AED.R")
source("vif.R")
source("summarySE.R")
source("resizewin.R")


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

ggplot(data=GC.sum.raw, aes(x=day, y=cellsml, color=medtreat, shape=medtreat)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsml-se, ymax=cellsml+se), width=0.2, size=1) + 
  geom_smooth(method="loess", aes(group=medtreat))


#compute some parameters
#Two points, N1 and N2,  at the extremes of this linear phase (see fig below) are taken and substituted into the equation
# Growth rate ;  K' = Ln (N2 / N1) / (t2 - t1)
# Divisions per day ; Div.day-1 = K' / Ln2
# Generation time ; Gen' t  = 1 / Div.day-1

##WHAT IF TRY DAYS 0-7
NT = data.table (P36_GC_limrev,  key=c("medtreatreplet" ))

GC.comp=NT [,list(K= (lncellsR[match("day7", day)]/lncellsR [match ("day0", day)])/7),   #divide according to number of days
            by=c("medtreatreplet")]

GC.comp [,div.day:=K/log(2)] #this is always 2
GC.comp [, genT:= 1/div.day]

GC.comp <- na.omit(GC.comp)

#put medtreat as a factor
medtreat <- rep(c ("dN limitation", "dN recovery", "dP limitation", "dP recovery", "dSi limitation", "dSi recovery", "ASW control" ), each=9)

GC.comp$medtreat <- as.factor(medtreat)


#summaries

GC.comp.sumK <- summarySE(data=GC.comp, measurevar="K", groupvars="medtreat")
GC.comp.sumdiv.day <- summarySE(data=GC.comp, measurevar="div.day", groupvars="medtreat")
GC.comp.sumgenT <- summarySE(data=GC.comp, measurevar="genT", groupvars="medtreat")

GC.comp.sum <- cbind(GC.comp.sumK, GC.comp.sumdiv.day, GC.comp.sumgenT)

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
bartlett.test(K~ medtreat, data=GC.comp) #homogenous


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


##more post

GC.kr <- kruskalmc(K~medtreat, GC.comp, probs=0.05)

# sometimes you want to look for "homogenous groups"
# and give each group a code (like a letter)
library(multcompView) # load library (install if necessary)
test <- GC.kr$dif.com$difference # select logical vector
names(test) <- row.names(GC.kr$dif.com)# add comparison names
# create a list with "homogenous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),
                       reversed = FALSE)

boxplot(K ~ medtreat, GC.comp
        ,xlab = "Treatment", ylab = "K"
        , ylim = c(0,max(GC.comp$K,na.rm=T)*1.5)
        , notch = F, pch = ".")
mtext(side=3,text=let$Letters,at=1:length(let$Letters),cex=1.5) # letters at top


