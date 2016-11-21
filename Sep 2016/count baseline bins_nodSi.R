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

#read data
count <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/count data/count_nodSi.csv", sep=";")

#correct time variable
count$T= count$time
count$T.factor= as.factor(count$T)

#make a new factor
count$wellvidbin <- as.factor(paste(count$wellvid, count$bin, sep = "-"))

#check initial plot
qplot(bin,cells, color = treatment, data = count,  geom = "boxplot") + facet_wrap(~treatment) 
qplot(as.factor(T),cells, color = treatment, data = count,  geom = "boxplot") + facet_grid(treatment~bin, scales="free") 

#standardization

count$cellsS=NA
k=split(count, count$treatment)
countstd <- lapply(k, function (x) scale(x[,c("cells")], center=T, scale=T))
count$cellsS=unsplit(countstd, count$treatment)

#baselining to 0 at time point 0
NT<-data.table(count, key=c("wellvidbin"))

t1=NT[,list(treatment=treatment, wellvid=wellvid, bin=bin, T=T, T.factor=as.factor(T), cells=cells, cellsS=cellsS,
            cellsBase=(cellsS-cellsS[1])), by=c("wellvidbin")]

countbase <- t1 #DATA IS NOW CALLED COUNTBASE

qplot(as.factor(T),cellsBase, color = treatment, data = countbase,  geom = "boxplot") + facet_grid(treatment~bin, scales="free") +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=1))

source("summarySE.R")
countbase.sum <- summarySE(countbase, measurevar="cellsBase", groupvars=c("T", "bin", "treatment"))

ggplot(data=countbase.sum, aes(x=T, y=cellsBase, shape=treatment, color=treatment)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=0.5, size=1) + facet_grid(bin~treatment)

#combine treatment and T in 1
countbase$treatT <- interaction(countbase$treatment, countbase$T)


##BINS

source ("AED.R")
source ("vif.R")
library (lawstat)

#bin A
binA= subset (countbase, bin=='BinA')
expA=as.data.frame(data.table(cbind(treatment=binA$treatment, T=binA$T, ID=binA$wellvidbin)))
cor(expA, method = "spearman")


vif_func(in_frame=expA,thresh=5,trace=T)

pairs(expA, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(binA$cellsBase, group=binA$wellvidbin, location="mean") #significant
levene.test(binA$cellsBase, group=binA$T, location="mean") # significant
levene.test(binA$cellsBase, group=binA$treatment, location="mean") #significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(cellsBase~treatment, data=binA)
boxplot(cellsBase~wellvidbin, data=binA)
boxplot (cellsBase~T, data=binA)

#fit a gls
Form <- formula (cellsBase ~ treatment*T)
binA.gls<- gls(Form, data=binA)


#nlme model
binA1.lme <- lme (Form, random = ~1|wellvidbin, method="REML", data=binA)

#binA2.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|wellvidbin), method="REML", data=binA)

#binA21.lme <- lme (Form, random = ~1|treatment/wellvidbin,  weights=varIdent(form=~1|wellvidbin), method="REML", data=binA) 

#binA3.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|wellvidbin), correlation=corAR1 (form=~1|wellvidbin/treatment), method="REML", data=binA) 

binA4.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|T), correlation=corAR1 (), method="REML", data=binA) 

binA5.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (), method="REML", data=binA) #BEST

binA6.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (form=~1|wellvidbin), method="REML", data=binA) #same as 5

binA7.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (form=~1|wellvidbin/treatment), method="REML", data=binA) #same as 5


anova(binA.gls, binA1.lme,  binA5.lme, binA6.lme, binA7.lme)


summary(binA5.lme)
anova(binA5.lme)

#residuals
binA.E2<-resid(binA5.lme,type="normalized")
binA.F2<-fitted(binA5.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=binA.F2,y=binA.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(binA.E2~treatment,data=binA, main="Treatment",ylab=MyYlab)
plot(x=binA$T,y=binA.E2,main="Time",ylab=MyYlab,xlab="Time (min)")
par(op)

xyplot (binA.E2 ~ T| treatment, data=binA, ylab="Residuals", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})

xyplot (binA.F2 ~ T| treatment, data=binA, ylab="Fitted", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})


#multimple comparisons
#use lstrends - lstrends for estimating and comparing the slopes of fitted lines (or curves).
#https://cran.r-project.org/web/packages/lsmeans/vignettes/using-lsmeans.pdf

#library(lsmeans)

library(multcomp)
#glht only works when T is numeric
summary(glht(binA5.lme, linfct=mcp(treatment="Tukey", covariate_average=TRUE, interaction_average = TRUE)))

#lstrends(binA5.lme, pairwise ~treatment, var="as.numeric(T)")
#lstrends only work for numeric not for factors, therefore use lsmeans


#lsmeans(binA5.lme, list(pairwise ~ treatment|T))
#pairs(lsmeans(binA5.lme, ~treatment))

#bin B
binB= subset (countbase, bin=='BinB')
expB=as.data.frame(data.table(cbind(treatment=binB$treatment, T=binB$T, ID=binB$wellvidbin)))
cor(expB, method = "spearman")


vif_func(in_frame=expB,thresh=5,trace=T)

pairs(expB, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(binB$cellsBase, group=binB$wellvidbin, location="mean") #significant
levene.test(binB$cellsBase, group=binB$T, location="mean") # not significant
levene.test(binB$cellsBase, group=binB$treatment, location="mean") #significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(cellsBase~treatment, data=binB)
boxplot(cellsBase~wellvidbin, data=binB)
boxplot (cellsBase~T, data=binB)

#fit a gls
Form <- formula (cellsBase ~ treatment*T)
binB.gls<- gls(Form, data=binB)


#nlme model
binB1.lme <- lme (Form, random = ~1|wellvidbin, method="REML", data=binB)

binB2.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|wellvidbin), method="REML", data=binB)

binB21.lme <- lme (Form, random = ~1|treatment/wellvidbin,  weights=varIdent(form=~1|wellvidbin), method="REML", data=binB) 

binB3.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|wellvidbin), correlation=corAR1 (), method="REML", data=binB) 

binB4.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|T), correlation=corAR1 (), method="REML", data=binB) 

binB5.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (), method="REML", data=binB) #BEST

binB6.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (form=~1|wellvidbin), method="REML", data=binB) #same as 5

binB7.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (form=~1|wellvidbin/treatment), method="REML", data=binB) #same as 5


anova(binB.gls, binB1.lme, binB2.lme, binB21.lme, binB3.lme, binB4.lme, binB5.lme, binB6.lme, binB7.lme)

summary(binB3.lme)
anova(binB3.lme)

#residuals
binB.E2<-resid(binB4.lme,type="normalized")
binB.F2<-fitted(binB4.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=binB.F2,y=binB.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(binB.E2~treatment,data=binB, main="Treatment",ylab=MyYlab)
plot(x=binB$T,y=binB.E2,main="Time",ylab=MyYlab,xlab="Time (min)")
par(op)

xyplot (binB.E2 ~ T| treatment, data=binB, ylab="Residuals", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})

xyplot (binB.F2 ~ T| treatment, data=binB, ylab="Fitted", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})


#multimple comparisons
#use lstrends - lstrends for estimating and comparing the slopes of fitted lines (or curves).
#https://cran.r-project.org/web/packages/lsmeans/vignettes/using-lsmeans.pdf

library(lsmeans)

summary(glht(binB3.lme, linfct=mcp(treatment="Tukey", covariate_average=TRUE, interaction_average = TRUE)))


#lstrends(binB5.lme, pairwise ~treatment, var="as.numeric(T)")
#lstrends only work for numeric not for factors, therefore use lsmeans


#lsmeans(binB5.lme, list(pairwise ~ treatment|T))
#pairs(lsmeans(binB5.lme, ~treatment))

#bin C
binC= subset (countbase, bin=='BinC')
expC=as.data.frame(data.table(cbind(treatment=binC$treatment, T=binC$T, ID=binC$wellvidbin)))
cor(expC, method = "spearman")


vif_func(in_frame=expC,thresh=5,trace=T)

pairs(expC, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(binC$cellsBase, group=binC$wellvidbin, location="mean") #significant
levene.test(binC$cellsBase, group=binC$T, location="mean") # significant
levene.test(binC$cellsBase, group=binC$treatment, location="mean") #not significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(cellsBase~treatment, data=binC)
boxplot(cellsBase~wellvidbin, data=binC)
boxplot (cellsBase~T, data=binC)

#fit a gls
Form <- formula (cellsBase ~ treatment*T)
binC.gls<- gls(Form, data=binC)


#nlme model
binC1.lme <- lme (Form, random = ~1|wellvidbin, method="REML", data=binC)

binC2.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|wellvidbin), method="REML", data=binC)

binC21.lme <- lme (Form, random = ~1|treatment/wellvidbin,  weights=varIdent(form=~1|wellvidbin), method="REML", data=binC) 

binC3.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|wellvidbin), correlation=corAR1 (), method="REML", data=binC) 

binC4.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|T), correlation=corAR1 (), method="REML", data=binC) 

binC5.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (), method="REML", data=binC) #BEST

binC6.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (form=~1|wellvidbin), method="REML", data=binC) #same as 5

binC7.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (form=~1|wellvidbin/treatment), method="REML", data=binC) #same as 5


anova(binC.gls, binC1.lme, binC2.lme, binC21.lme, binC3.lme, binC4.lme, binC5.lme, binC6.lme, binC7.lme)

#use binC3.lme for consistency

summary(binC3.lme)
anova(binC3.lme)


#residuals
binC.E2<-resid(binC3.lme,type="normalized")
binC.F2<-fitted(binC3.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=binC.F2,y=binC.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(binC.E2~treatment,data=binC, main="Treatment",ylab=MyYlab)
plot(x=binC$T,y=binC.E2,main="Time",ylab=MyYlab,xlab="Time (min)")
par(op)

xyplot (binC.E2 ~ T| treatment, data=binC, ylab="Residuals", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})

xyplot (binC.F2 ~ T| treatment, data=binC, ylab="Fitted", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})


#multimple comparisons
#use lstrends - lstrends for estimating and comparing the slopes of fitted lines (or curves).
#https://cran.r-project.org/web/packages/lsmeans/vignettes/using-lsmeans.pdf

library(lsmeans)

#lstrends(binC5.lme, pairwise ~treatment, var="as.numeric(T)")
#lstrends only work for numeric not for factors, therefore use lsmeans

summary(glht(binC3.lme, linfct=mcp(treatment="Tukey", covariate_average=TRUE, interaction_average = TRUE)))

#lsmeans(binC5.lme, list(pairwise ~ treatment|T))
#pairs(lsmeans(binC5.lme, ~treatment))

##LET'S PLOT THIS

library(AICcmodavg)

#Bin A fit

BA.fit <- as.data.frame(predictSE.lme(binA5.lme, binA, se.fit = TRUE, level = 0, print.matrix = FALSE))

BA.fit$upr <- BA.fit$fit + (1.96 * BA.fit$se)
BA.fit$lwr <- BA.fit$fit - (1.96 * BA.fit$se)

BA.fit.combdata <- cbind(binA, BA.fit)

#Bin B fit

BB.fit <- as.data.frame(predictSE.lme(binB3.lme, binB, se.fit = TRUE, level = 0, print.matrix = FALSE))

BB.fit$upr <- BB.fit$fit + (1.96 * BB.fit$se)
BB.fit$lwr <- BB.fit$fit - (1.96 * BB.fit$se)

BB.fit.combdata <- cbind(binB, BB.fit)

#Bin C fit

BC.fit <- as.data.frame(predictSE.lme(binC3.lme, binC, se.fit = TRUE, level = 0, print.matrix = FALSE))

BC.fit$upr <- BC.fit$fit + (1.96 * BC.fit$se)
BC.fit$lwr <- BC.fit$fit - (1.96 * BC.fit$se)

BC.fit.combdata <- cbind(binC, BC.fit)

#summaries
#you cannot have Y as a factor!!! so use T
BinA.sum <- summarySE(binA, measurevar="cellsBase", groupvars=c("T", "treatment"))

BinB.sum <- summarySE(binB, measurevar="cellsBase", groupvars=c("T", "treatment"))

BinC.sum <- summarySE(binC, measurevar="cellsBase", groupvars=c("T", "treatment"))

BinA.sum$bin="BinA"
BinB.sum$bin="BinB"
BinC.sum$bin="BinC"

allbins.sum = rbind (BinA.sum, BinB.sum, BinC.sum)


#plot

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 

cbPalette <- c("#999999",  "#009E73", "#56B4E9","#F0E442", "#0072B2", "#D55E00", "#CC79A7")

allbins.fitdata = rbind (BA.fit.combdata, BB.fit.combdata, BC.fit.combdata)

allbins.fitdata$treatment <- factor(allbins.fitdata$treatment, levels=c("control bead", "P bead"), 
                                    labels =c (" control bead  "," dP bead  "))
allbins.sum$treatment <- factor(allbins.sum$treatment, levels=c("control bead", "P bead"), 
                                labels =c (" control bead  "," dP bead  "))

allbins.sum$bin <- factor(allbins.sum$bin, levels=c("BinA", "BinB", "BinC"), 
                                labels =c ("Bin A", "Bin B", "Bin C"))

allbins.fitdata$bin <- factor(allbins.fitdata$bin, levels=c("BinA", "BinB", "BinC"), 
                          labels =c ("Bin A", "Bin B", "Bin C"))

resize.win(9,12)

ggplot(data=allbins.sum, aes(x=T, y=cellsBase, shape=treatment)) + 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=2, size=1) +
  geom_point(size=5, shape = 21, color='black',
             aes(fill = treatment)) +
  geom_smooth(data=allbins.fitdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr), 
              color="black", method="lm", stat="identity", alpha=0.2)+ 
  facet_grid(bin~.) +
  scale_fill_manual(values = c('white', 'black')) +
  labs(list(x = "Time (min)", y = "Normalized cell count", title = "Attraction of dP-starved cells \nto dP-loaded beads"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20, vjust=1.5), 
        axis.title.x=element_text(size=20, vjust=-0.5),
        plot.title = element_text(size =24), axis.text=text,  legend.position="bottom", legend.title=element_blank(),
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))

