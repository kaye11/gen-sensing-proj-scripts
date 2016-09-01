
#reading the data

phosdist <- read.table ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/distance preprocessed.csv", 
                         sep=";", header=T)

library(ggplot2)
library(grid)
library(gtable)
library(ggthemes)
library(gridExtra)
library(data.table)
library(plyr)
library(mgcv)
source("vif.R")
source ("AED.R")
source("lang.R")
source("summarySE.R")
source("lang.R")
source("tsDiagGamm.R")
source("resizewin.R")
resize.win(12,9)

phosdist$newbin <- relevel(phosdist$newbin, "BinsAB")

qplot(timemin, distmm, color = treatment, data = phosdist)+ stat_smooth(aes(group=treatment), method="lm")+facet_grid(newbin~., scales="free")

qplot(x=timemin, y=distmm, data=phosdist, color=treatment)+geom_line()+facet_grid(treatment~newbin, scales="free")

#BINS AB

BinsAB= subset (phosdist, newbin=='BinsAB')

##check for collinearity and correlation, this only applies to the explanatory variables!
expA=as.data.frame(data.table(cbind(treatment=BinsAB$treatment, T=BinsAB$timemin)))
cor(expA, method = "spearman")

vif_func(in_frame=expA,thresh=5,trace=T)

pairs(expA, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

#boxplots
op=par(mfrow=c(1,2))
boxplot(distmm~treatment, data=BinsAB)
boxplot(distmm~timemin, data=BinsAB)

#levene
library(lawstat)
levene.test(BinsAB$distmm, group=BinsAB$treatment, location="mean") #unequal
levene.test(BinsAB$distmm, group=BinsAB$timemin, location="mean") #unequal

#try linear model

BA.gls <- gls (distmm~timemin*treatment, method="REML", data=BinsAB)
BA.lm <- lm (distmm~timemin*treatment,  data=BinsAB) 
#BA.lme <- lme (distmm~timemin*treatment, method="REML", correlation=corAR1(), data=BinsAB ) #lme doesn't work
#residuals of BA.lm shows nonlinearity but normal qq plots so choose lm

BA.gls1 <- gls (distmm~timemin*treatment, method="REML", correlation=corAR1(), data=BinsAB)
BA.gls2 <- gls (distmm~timemin*treatment, method="REML", correlation= corAR1 (form=~1|treatment), data=BinsAB) #BEST
BA.gls3 <- gls (distmm~timemin*treatment, method="REML", correlation= corAR1 (form=~1|treatment/timemin), data=BinsAB)
BA.gls4 <- gls (distmm~timemin*treatment, method="REML", correlation= corAR1 (form=~1|timemin), data=BinsAB)

anova(BA.gls, BA.lm, BA.gls1, BA.gls2, BA.gls3, BA.gls4)
#BA.gls2 produced lowest AIC values

summary(BA.gls2)
anova(BA.gls2)

#residuals
BinsAB.E2<-resid(BA.gls2,type="response")
BinsAB.F2<-fitted(BA.gls2)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=BinsAB.F2,y=BinsAB.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(BinsAB.E2~treatment,data=BinsAB, main="treatment",ylab=MyYlab)
plot(x=BinsAB$timemin,y=BinsAB.E2,main="timemin",ylab=MyYlab,xlab="time (min)")
par(op)

xyplot (BinsAB.E2 ~ timemin| treatment, data=BinsAB, ylab="Residuals", xlab="time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})


#multimple comparisons
library(lsmeans)
BinsAB.lsmeans=lsmeans(BA.gls2,  ~ treatment)
summary(BinsAB.lsmeans, adjust="bonf")
contrast(BinsAB.lsmeans, adjust="bonf")
pairs(BinsAB.lsmeans)

library(multcompView)
cld(BinsAB.lsmeans, alpha=0.05)


##BIN C

BinC= subset (phosdist, newbin=='BinC')

##check for collinearity and correlation, this only applies to the explanatory variables!
expC=as.data.frame(data.table(cbind(treatment=BinC$treatment, T=BinC$timemin)))
cor(expC, method = "spearman")

vif_func(in_frame=expC,thresh=5,trace=T)

pairs(expC, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

#boxplots
op=par(mfrow=c(1,2))
boxplot(distmm~treatment, data=BinC)
boxplot(distmm~timemin, data=BinC)

#levene
library(lawstat)
levene.test(BinC$distmm, group=BinC$treatment, location="mean") #unequal
levene.test(BinC$distmm, group=BinC$timemin, location="mean") #equal

#try linear model

BC.gls <- gls (distmm~timemin*treatment, method="REML", data=BinC)
BC.lm <- lm (distmm~timemin*treatment,  data=BinC) 
#BC.lme <- lme (distmm~timemin*treatment, method="REML", correlation=corAR1(), data=BinC ) #lme doesn't work
#residuals of BC.lm shows nonlinearity but normal qq plots so choose lm

BC.gls1 <- gls (distmm~timemin*treatment, method="REML", correlation=corAR1(), data=BinC)
BC.gls2 <- gls (distmm~timemin*treatment, method="REML", correlation= corAR1 (form=~1|treatment), data=BinC) #BEST
BC.gls3 <- gls (distmm~timemin*treatment, method="REML", correlation= corAR1 (form=~1|treatment/timemin), data=BinC)
BC.gls4 <- gls (distmm~timemin*treatment, method="REML", correlation= corAR1 (form=~1|timemin), data=BinC)

anova(BC.gls, BC.lm, BC.gls1, BC.gls2, BC.gls3, BC.gls4)
#BC.gls2 produced lowest AIC values

summary(BC.gls2)
anova(BC.gls2)

#residuals
BinC.E2<-resid(BC.gls2,type="response")
BinC.F2<-fitted(BC.gls2)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=BinC.F2,y=BinC.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(BinC.E2~treatment,data=BinC, main="treatment",ylab=MyYlab)
plot(x=BinC$timemin,y=BinC.E2,main="timemin",ylab=MyYlab,xlab="time (min)")
par(op)

xyplot (BinC.E2 ~ timemin| treatment, data=BinC, ylab="Residuals", xlab="time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})


library(lsmeans)
BinC.lsmeans=lsmeans(BC.gls2,  ~ treatment)
summary(BinC.lsmeans, adjust="bonf")
contrast(BinC.lsmeans, adjust="bonf")
pairs(BinC.lsmeans)

library(multcompView)
cld(BinC.lsmeans, alpha=0.05)

##LET'S PLOT THIS

library(AICcmodavg)

#BA fit

BA.fit <- as.data.frame(predictSE.gls(BA.gls2, BinsAB, se.fit = TRUE, level = 0, print.matrix = FALSE))

#BC.fit <- as.data.frame(predict(BC.gls2, BinC, se.fit = TRUE, level = 0, print.matrix = FALSE))

BA.fit$upr <- BA.fit$fit + (1.96 * BA.fit$se)
BA.fit$lwr <- BA.fit$fit - (1.96 * BA.fit$se)

BA.fit.combdata <- cbind(BinsAB, BA.fit)

#BC fit

BC.fit <- as.data.frame(predictSE.gls(BC.gls2, BinC, se.fit = TRUE, level = 0, print.matrix = FALSE))

#BC.fit <- as.data.frame(predict(BC.gls2, BinC, se.fit = TRUE, level = 0, print.matrix = FALSE))

BC.fit$upr <- BC.fit$fit + (1.96 * BC.fit$se)
BC.fit$lwr <- BC.fit$fit - (1.96 * BC.fit$se)

BC.fit.combdata <- cbind(BinC, BC.fit)

#plot

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 

cbPalette <- c("#999999",  "#009E73", "#56B4E9","#F0E442", "#0072B2", "#D55E00", "#CC79A7")

allbins.fitdata = rbind (BA.fit.combdata, BC.fit.combdata)

mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="newbin") { 
    value[value=="BinsAB"] <- "Bins A+B"
    value[value=="BinC"]   <- "Bin C"
  }
  return(value)
}

resize.win(15,12)

ggplot(data=phosdist, aes(x=timemin, y=distmm, shape=treatment, color=treatment)) + geom_point(size=5)+ 
  geom_smooth(data=allbins.fitdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, fill=treatment), stat="identity", alpha=0.2)+ 
  scale_color_manual(values = cbPalette, name="Treatment") + facet_grid(treatment~newbin, labeller=mf_labeller, scale="free") +
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") + 
  labs(list(x = "Time (min)", y = "Sum distance from the bead (mm)"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))

