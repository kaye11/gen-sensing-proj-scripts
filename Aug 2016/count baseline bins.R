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
count <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/count data/count.csv", sep=";")

#correct time variable
count$T= count$time
count$T.factor= as.factor(count$T)

#make a new factor
count$wellvidbin <- as.factor(paste(count$wellvid, count$bin, sep = "-"))

#check initial plot
qplot(bin,cells, color = treatment, data = count,  geom = "boxplot") + facet_wrap(~treatment) 
qplot(T.factor,cells, color = treatment, data = count,  geom = "boxplot") + facet_grid(treatment~bin, scales="free") 

#standardization

count$cellsS=NA
k=split(count, count$treatment)
countstd <- lapply(k, function (x) scale(x[,c("cells")], center=T, scale=T))
count$cellsS=unsplit(countstd, count$treatment)

#baselining to 0 at time point 0
NT<-data.table(count, key=c("wellvidbin"))

t1=NT[,list(treatment=treatment, bin=bin, T=T, T.factor=T.factor, cells=cells, cellsS=cellsS,
            cellsBase=(cellsS-cellsS[1])), by=c("wellvidbin")]

countbase <- t1 #DATA IS NOW CALLED COUNTBASE

qplot(as.factor(T),cellsBase, color = treatment, data = countbase,  geom = "boxplot") + facet_grid(treatment~bin, scales="free") +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=1))

source("summarySE.R")
countbase.sum <- summarySE(countbase, measurevar="cellsBase", groupvars=c("T.factor", "bin", "treatment"))

ggplot(data=countbase.sum, aes(x=T.factor, y=cellsBase, shape=treatment, color=treatment)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=0.5, size=1) + facet_grid(bin~treatment)

#combine treatment and T.factor in 1
countbase$treatT <- interaction(countbase$treatment, countbase$T.factor)

##BINS

source ("AED.R")
source ("vif.R")
library (lawstat)

#bin A
binA= subset (countbase, bin=='BinA')
expA=as.data.frame(data.table(cbind(treatment=binA$treatment, T=binA$T.factor, ID=binA$wellvidbin)))
cor(expA, method = "spearman")


vif_func(in_frame=expA,thresh=5,trace=T)

pairs(expA, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(binA$cellsBase, group=binA$wellvidbin, location="mean") #significant
levene.test(binA$cellsBase, group=binA$T.factor, location="mean") # significant
levene.test(binA$cellsBase, group=binA$treatment, location="mean") #significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(cellsBase~treatment, data=binA)
boxplot(cellsBase~wellvidbin, data=binA)
boxplot (cellsBase~T.factor, data=binA)

#fit a gls
Form <- formula (cellsBase ~ treatment*T.factor)
binA.gls<- gls(Form, data=binA)


#nlme model
binA1.lme <- lme (Form, random = ~1|wellvidbin, method="REML", data=binA)

#binA2.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|wellvidbin), method="REML", data=binA)

binA21.lme <- lme (Form, random = ~1|treatment/wellvidbin,  weights=varIdent(form=~1|wellvidbin), method="REML", data=binA) 

binA3.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|wellvidbin), 
                  correlation=corAR1 (form=~1|wellvidbin/treatment), method="REML", data=binA) 

binA4.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|T.factor), correlation=corAR1 (), method="REML", data=binA) 

binA5.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (), method="REML", data=binA) #BEST

binA6.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (form=~1|wellvidbin), method="REML", data=binA) #same as 5


anova(binA.gls, binA1.lme, binA3.lme, binA4.lme, binA5.lme, binA6.lme)


#gls model
BA.gls1 <- gls (Form, method="REML", correlation=corAR1(), data=binA)
BA.gls2 <- gls (Form, method="REML", correlation= corAR1 (form=~1|treatment), data=binA) 
BA.gls3 <- gls (Form, method="REML", correlation= corAR1 (form=~1|treatment/wellvidbin), data=binA) #BEST
BA.gls4 <- gls (Form, method="REML", correlation= corAR1 (form=~1|wellvidbin), data=binA)

BA.gls5 <- gls (Form, weights=varIdent(form=~1|treatment), 
                correlation= corAR1 (form=~1|treatment/wellvidbin), method="REML", data=binA) 

BA.gls6 <- gls (Form, weights=varIdent(form=~1|treatment), method="REML", data=binA) 

anova(BA.gls1, BA.gls2, BA.gls3, BA.gls4, BA.gls5, BA.gls6)

anova(BA.gls3, binA5.lme) # no difference between 2, choose lme because statistically its more correct with the random formulation

summary(binA5.lme)
anova(binA5.lme)

#combine treatment and T. factor into one variable

binA5.lmecomb <- lme (cellsBase ~ treatT, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                      correlation=corAR1 (), method="REML", data=binA) #BEST

anova(binA5.lmecomb) # same pvalue as intercept and interaction with binA5.lme

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

library(lsmeans)

#lstrends(binA5.lme, pairwise ~treatment, var="as.numeric(T.factor)")
#lstrends only work for numeric not for factors, therefore use lsmeans


lsmeans(binA5.lme, list(pairwise ~ treatment|T.factor))
pairs(lsmeans(binA5.lme, ~treatment))

#bin B
binB= subset (countbase, bin=='BinB')
expB=as.data.frame(data.table(cbind(treatment=binB$treatment, T=binB$T.factor, ID=binB$wellvidbin)))
cor(expB, method = "spearman")


vif_func(in_frame=expB,thresh=5,trace=T)

pairs(expB, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(binB$cellsBase, group=binB$wellvidbin, location="mean") #significant
levene.test(binB$cellsBase, group=binB$T.factor, location="mean") # significant
levene.test(binB$cellsBase, group=binB$treatment, location="mean") #significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(cellsBase~treatment, data=binB)
boxplot(cellsBase~wellvidbin, data=binB)
boxplot (cellsBase~T.factor, data=binB)

#fit a gls
Form <- formula (cellsBase ~ treatment*T.factor)
binB.gls<- gls(Form, data=binB)


#nlme model
binB1.lme <- lme (Form, random = ~1|wellvidbin, method="REML", data=binB)

#binB2.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|wellvidbin), method="REML", data=binB)

#binB21.lme <- lme (Form, random = ~1|treatment/wellvidbin,  weights=varIdent(form=~1|wellvidbin), method="REML", data=binB) 

#binB3.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|wellvidbin), correlation=corAR1 (form=~1|wellvidbin/treatment), method="REML", data=binB) 

binB4.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|T.factor), correlation=corAR1 (), method="REML", data=binB) 

binB5.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (), method="REML", data=binB) #BEST

binB6.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (form=~1|wellvidbin), method="REML", data=binB) #same as 5


anova(binB.gls, binB1.lme, binB4.lme, binB5.lme, binB6.lme)

summary(binB5.lme)
anova(binB5.lme)

#combine treatment and T. factor into one variable

binB5.lmecomb <- lme (cellsBase ~ treatT, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                      correlation=corAR1 (), method="REML", data=binB) #BEST

anova(binB5.lmecomb) # same pvalue as intercept and interaction with binB5.lme

#residuals
binB.E2<-resid(binB5.lme,type="normalized")
binB.F2<-fitted(binB5.lme)
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

#lstrends(binB5.lme, pairwise ~treatment, var="as.numeric(T.factor)")
#lstrends only work for numeric not for factors, therefore use lsmeans


lsmeans(binB5.lme, list(pairwise ~ treatment|T.factor))
pairs(lsmeans(binB5.lme, ~treatment))

#bin C
binC= subset (countbase, bin=='BinC')
expC=as.data.frame(data.table(cbind(treatment=binC$treatment, T=binC$T.factor, ID=binC$wellvidbin)))
cor(expC, method = "spearman")


vif_func(in_frame=expC,thresh=5,trace=T)

pairs(expC, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(binC$cellsBase, group=binC$wellvidbin, location="mean") #significant
levene.test(binC$cellsBase, group=binC$T.factor, location="mean") # significant
levene.test(binC$cellsBase, group=binC$treatment, location="mean") #not significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(cellsBase~treatment, data=binC)
boxplot(cellsBase~wellvidbin, data=binC)
boxplot (cellsBase~T.factor, data=binC)

#fit a gls
Form <- formula (cellsBase ~ treatment*T.factor)
binC.gls<- gls(Form, data=binC)


#nlme model
binC1.lme <- lme (Form, random = ~1|wellvidbin, method="REML", data=binC)

binC2.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|wellvidbin), method="REML", data=binC)

binC21.lme <- lme (Form, random = ~1|treatment/wellvidbin,  weights=varIdent(form=~1|wellvidbin), method="REML", data=binC) 

binC3.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|wellvidbin), correlation=corAR1 (form=~1|wellvidbin/treatment), method="REML", data=binC) 

#binC4.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|T.factor), correlation=corAR1 (), method="REML", data=binC) 

binC5.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (), method="REML", data=binC) #BEST

binC6.lme <- lme (Form, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                  correlation=corAR1 (form=~1|wellvidbin), method="REML", data=binC) #same as 5


anova(binC.gls, binC1.lme, binC2.lme, binC21.lme, binC3.lme, binC5.lme, binC6.lme)

#for consistency use binC5.lme

summary(binC5.lme)
anova(binC5.lme)

#combine treatment and T. factor into one variable

binC5.lmecomb <- lme (cellsBase ~ treatT, random = ~1|wellvidbin,  weights=varIdent(form=~1|treatment), 
                      correlation=corAR1 (), method="REML", data=binC) #BEST

anova(binC5.lmecomb) # same pvalue as intercept and interaction with binC5.lme

#residuals
binC.E2<-resid(binC5.lme,type="normalized")
binC.F2<-fitted(binC5.lme)
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

#lstrends(binC5.lme, pairwise ~treatment, var="as.numeric(T.factor)")
#lstrends only work for numeric not for factors, therefore use lsmeans


lsmeans(binC5.lme, list(pairwise ~ treatment|T.factor))
pairs(lsmeans(binC5.lme, ~treatment))

##LET'S PLOT THIS

library(AICcmodavg)

#Bin A fit

BA.fit <- as.data.frame(predictSE.lme(binA5.lme, binA, se.fit = TRUE, level = 0, print.matrix = FALSE))

BA.fit$upr <- BA.fit$fit + (1.96 * BA.fit$se)
BA.fit$lwr <- BA.fit$fit - (1.96 * BA.fit$se)

BA.fit.combdata <- cbind(binA, BA.fit)

#Bin B fit

BB.fit <- as.data.frame(predictSE.lme(binB5.lme, binB, se.fit = TRUE, level = 0, print.matrix = FALSE))

BB.fit$upr <- BB.fit$fit + (1.96 * BB.fit$se)
BB.fit$lwr <- BB.fit$fit - (1.96 * BB.fit$se)

BB.fit.combdata <- cbind(binB, BB.fit)

#Bin C fit

BC.fit <- as.data.frame(predictSE.lme(binC5.lme, binC, se.fit = TRUE, level = 0, print.matrix = FALSE))

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

mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="bin") { 
    value[value=="BinA"] <- "Bin A"
    value[value=="BinB"]   <- "Bin B"
    value[value=="BinC"]   <- "Bin C"
  }
  return(value)
}

resize.win(15,12)

allbins.fitdata$treatment2 <- factor(allbins.fitdata$treatment, levels=c("control bead", "P bead", "Si bead"), 
                                 labels =c ("control bead", "dP bead", "dSi bead"))

ggplot(data=allbins.sum, aes(x=T, y=cellsBase)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=2, size=1) +
  geom_smooth(data=allbins.fitdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr), color="black", method="lm", stat="identity", alpha=0.2)+ 
  facet_grid(treatment2~bin, labeller=mf_labeller, scale="free") +
  labs(list(x = "Time (min)", y = "Normalized cell count"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))




ggplot(data=allbins.sum, aes(x=T, y=cellsBase, shape=treatment, color=treatment)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=2, size=1) +
  geom_smooth(data=allbins.fitdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, fill=treatment), stat="identity", alpha=0.2)+ 
  scale_color_manual(values = cbPalette, name="Treatment") + facet_grid(treatment~bin, labeller=mf_labeller, scale="free") +
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") + 
  labs(list(x = "Time (min)", y = "Normalized cell count"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))

