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

panel_85AS <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/85AS_panel_longformat.csv", sep=";")

#correct time variable
panel_85AS$T= panel_85AS$timemin
panel_85AS$T.factor= as.factor(panel_85AS$T)

#make a new factor
panel_85AS$condind <- as.factor(paste(panel_85AS$condition, panel_85AS$induction, sep = "-"))
panel_85AS$condindbead <- as.factor(paste(panel_85AS$condind, panel_85AS$bead, sep = "-"))
panel_85AS$wellvidbead <- as.factor(paste(panel_85AS$wellvid, panel_85AS$bead, sep = "-"))
panel_85AS$wellvidcondind <- as.factor(paste(panel_85AS$wellvid, panel_85AS$condind, sep = "-"))

#check initial plot
qplot(condind,cells, color = bead, data = panel_85AS,  geom = "boxplot")
qplot(T.factor,cells, color = bead, data = panel_85AS,  geom = "boxplot") + facet_grid(condition~induction, scales="free")
qplot(T.factor,cells, data = panel_85AS,  geom = "boxplot") + facet_grid(bead~condind, scales="free") 


#standardization

panel_85AS$cellsS=NA
k=split(panel_85AS, panel_85AS$condind)
panel_85ASstd <- lapply(k, function (x) scale(x[,c("cells")], center=T, scale=T))
panel_85AS$cellsS=unsplit(panel_85ASstd, panel_85AS$condind)

#baselining to 0 at time point 0
NT<-data.table(panel_85AS, key=c("wellvidbead"))

t1=NT[,list(condition=condition, induction=induction, condind=condind, bead=bead, T=T, T.factor=T.factor, cells=cells, cellsS=cellsS,
            cellsBase=(cellsS-cellsS[1])), by=c("wellvidbead")]

countbase <- t1 #DATA IS NOW CALLED COUNTBASE

qplot(as.factor(T),cellsBase, color = condind, data = countbase,  geom = "boxplot") + facet_grid(bead~condind, scales="free") +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=1))

qplot(as.factor(T),cellsBase, color = condind, data = countbase,  geom = "boxplot") + facet_grid(bead~.,  scales="free") +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=condind, fill=condind), alpha=0.2)

source ("summarySE.R")
countbase.sum <- summarySE(countbase, measurevar="cellsBase", groupvars=c("T.factor", "bead", "condind"))

countbase.sum.cells <- summarySE(countbase, measurevar="cells", groupvars=c("T.factor", "bead", "condind"))

ggplot(data=countbase.sum, aes(x=T.factor, y=cellsBase, shape=condind, color=condind)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=0.5, size=1) + facet_grid(bead~condind)

ggplot(data=countbase.sum.cells, aes(x=T.factor, y=cells, shape=condind, color=condind)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cells-se, ymax=cells+se), width=0.5, size=1) + facet_grid(bead~condind)

##subset the bead, so that you compare responses of cells in one bead depending on condind

#DPRbead

source ("AED.R")
source ("vif.R")
library (lawstat)

DPRbead= subset (countbase, bead=='DPRbead')
expDPR=as.data.frame(data.table(cbind(condind=DPRbead$condind, T=DPRbead$T.factor, ID=DPRbead$wellvidbead)))
cor(expDPR, method = "spearman")


vif_func(in_frame=expDPR,thresh=5,trace=T)

pairs(expDPR, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(DPRbead$cellsBase, group=DPRbead$wellvidbead, location="mean") #significant
levene.test(DPRbead$cellsBase, group=DPRbead$T.factor, location="mean") # significant
levene.test(DPRbead$cellsBase, group=DPRbead$condind, location="mean") # significant
#levene.test(DPRbead$cellsBase, group=DPRbead$condition, location="mean") # significant
#levene.test(DPRbead$cellsBase, group=DPRbead$induction, location="mean") # not significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(cellsBase~condind, data=DPRbead)
boxplot(cellsBase~wellvidbead, data=DPRbead)
boxplot (cellsBase~T.factor, data=DPRbead)
#boxplot(cellsBase~condition, data=DPRbead)
#boxplot(cellsBase~induction, data=DPRbead)

#fit a gls
Form <- formula (cellsBase ~ condind*T)
DPRbead.gls<- gls(Form, data=DPRbead)


#nlme model

#random factor
DPRbead1.lme <- lme (Form, random = ~1|wellvidbead, method="REML", data=DPRbead) #BEST

DPRbead2.lme <- lme (Form, random = ~1|condind, method="REML", data=DPRbead)

DPRbead3.lme <- lme (Form, random = ~1|wellvidbead/condind, method="REML", data=DPRbead)

anova(DPRbead.gls, DPRbead1.lme, DPRbead2.lme, DPRbead3.lme)

#variance structure

#DPRbead4.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|wellvidbead), method="REML", data=DPRbead)

DPRbead5.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|condind), method="REML", data=DPRbead) #BEST

#DPRbead6.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|wellvidbead/condind), method="REML", data=DPRbead) 

anova(DPRbead.gls, DPRbead1.lme, DPRbead2.lme, DPRbead3.lme, DPRbead5.lme)

#correlation structures

DPRbead7.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|condind), correlation=corAR1(), method="REML", data=DPRbead) #BEST use this

DPRbead8.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|condind), 
                     correlation=corAR1 (form=~1|wellvidbead), method="REML", data=DPRbead) 

DPRbead9.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|condind), 
                     correlation=corAR1 (form=~1|wellvidbead/condind), method="REML", data=DPRbead) 


anova(DPRbead.gls, DPRbead1.lme, DPRbead2.lme, DPRbead3.lme, DPRbead5.lme, DPRbead7.lme, DPRbead8.lme, DPRbead9.lme)

#models 7-9 are the same are the same because it treats the correlation structures the same! HAHA

summary(DPRbead7.lme)
anova(DPRbead7.lme)


#multiple comparisons
library(lsmeans)

pairs(lsmeans(DPRbead7.lme, ~condind|T))

library(multcompView)
cld(lsmeans(DPRbead7.lme, ~condind), alpha=0.05)

#residuals
DPRbead.E2<-resid(DPRbead7.lme,type="normalized")
DPRbead.F2<-fitted(DPRbead7.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=DPRbead.F2,y=DPRbead.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(DPRbead.E2~condind,data=DPRbead, main="condind",ylab=MyYlab)
plot(x=DPRbead$T,y=DPRbead.E2,main="Time",ylab=MyYlab,xlab="Time (min)")
par(op)

xyplot (DPRbead.E2 ~ T| condind, data=DPRbead, ylab="Residuals", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})


#dSi bead

dSibead= subset (countbase, bead=='dSibead')
expdSi=as.data.frame(data.table(cbind(condind=dSibead$condind, T=dSibead$T.factor, ID=dSibead$wellvidbead)))
cor(expdSi, method = "spearman")


vif_func(in_frame=expdSi,thresh=5,trace=T)

pairs(expdSi, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(dSibead$cellsBase, group=dSibead$wellvidbead, location="mean") #significant
levene.test(dSibead$cellsBase, group=dSibead$T.factor, location="mean") # significant
levene.test(dSibead$cellsBase, group=dSibead$condind, location="mean") # significant
#levene.test(dSibead$cellsBase, group=dSibead$condition, location="mean") # significant
#levene.test(dSibead$cellsBase, group=dSibead$induction, location="mean") # not significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(cellsBase~condind, data=dSibead)
boxplot(cellsBase~wellvidbead, data=dSibead)
boxplot (cellsBase~T.factor, data=dSibead)
#boxplot(cellsBase~condition, data=dSibead)
#boxplot(cellsBase~induction, data=dSibead)

#fit a gls
Form <- formula (cellsBase ~ condind*T)
dSibead.gls<- gls(Form, data=dSibead)


#nlme model

#random factor
dSibead1.lme <- lme (Form, random = ~1|wellvidbead, method="REML", data=dSibead) #BEST

dSibead2.lme <- lme (Form, random = ~1|condind, method="REML", data=dSibead)

dSibead3.lme <- lme (Form, random = ~1|wellvidbead/condind, method="REML", data=dSibead)

anova(dSibead.gls, dSibead1.lme, dSibead2.lme, dSibead3.lme)

#variance structure

#dSibead4.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|wellvidbead), method="REML", data=dSibead)

dSibead5.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|condind), method="REML", data=dSibead) #BEST

#dSibead6.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|wellvidbead/condind), method="REML", data=dSibead) 

anova(dSibead.gls, dSibead1.lme, dSibead2.lme, dSibead3.lme, dSibead5.lme)

#correlation structures

dSibead7.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|condind), correlation=corAR1(), method="REML", data=dSibead) #BEST use this

dSibead8.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|condind), 
                     correlation=corAR1 (form=~1|wellvidbead), method="REML", data=dSibead) 

dSibead9.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|condind), 
                     correlation=corAR1 (form=~1|wellvidbead/condind), method="REML", data=dSibead) 


anova(dSibead.gls, dSibead1.lme, dSibead2.lme, dSibead3.lme, dSibead5.lme, dSibead7.lme, dSibead8.lme, dSibead9.lme)

#models 7-9 are the same are the same because it treats the correlation structures the same! HAHA

summary(dSibead7.lme)
anova(dSibead7.lme)


#multiple comparisons
library(lsmeans)

pairs(lsmeans(dSibead7.lme, ~condind|T))

library(multcompView)
cld(lsmeans(dSibead7.lme, ~condind), alpha=0.05)

#residuals
dSibead.E2<-resid(dSibead7.lme,type="normalized")
dSibead.F2<-fitted(dSibead7.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=dSibead.F2,y=dSibead.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(dSibead.E2~condind,data=dSibead, main="condind",ylab=MyYlab)
plot(x=dSibead$T,y=dSibead.E2,main="Time",ylab=MyYlab,xlab="Time (min)")
par(op)

xyplot (dSibead.E2 ~ T| condind, data=dSibead, ylab="Residuals", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})

#fitting data

library(AICcmodavg)

#DPR fit
DPRbead.fit <- as.data.frame(predictSE.lme(DPRbead7.lme, DPRbead, se.fit = TRUE, level = 0,
                                           print.matrix = FALSE))

DPRbead.fit$upr <- DPRbead.fit$fit + (1.96 * DPRbead.fit$se)
DPRbead.fit$lwr <- DPRbead.fit$fit - (1.96 * DPRbead.fit$se)

DPRbead.fit.combdata <- cbind(DPRbead, DPRbead.fit)

#dSi fit
dSibead.fit <- as.data.frame(predictSE.lme(dSibead7.lme, dSibead, se.fit = TRUE, level = 0,
                                           print.matrix = FALSE))

dSibead.fit$upr <- dSibead.fit$fit + (1.96 * dSibead.fit$se)
dSibead.fit$lwr <- dSibead.fit$fit - (1.96 * dSibead.fit$se)

dSibead.fit.combdata <- cbind(dSibead, dSibead.fit)

#summaries
DPRbead.sum <- summarySE(DPRbead, measurevar="cellsBase", groupvars=c("T", "bead", "condind"))

dSibead.sum <- summarySE(dSibead, measurevar="cellsBase", groupvars=c("T", "bead", "condind"))


#plot

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 

cbPalette <- c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plotshapes <- c(15, 16, 0, 1)

allbead.fitdata = rbind (DPRbead.fit.combdata, dSibead.fit.combdata)
allbead.sum <- rbind (DPRbead.sum, dSibead.sum)

mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="condind") { 
    value[value=="ASW-induced"] <- "ASW, induced"
    value[value=="ASW-notinduced"]   <- "ASW, not induced"
    value[value=="dSi-induced"] <- "-Si, induced"
    value[value=="dSi-notinduced"]   <- "-Si, not induced"
  }
  return(value)
}

scaleFUN <- function(x) sprintf("%.1f", x)

allbead.sum$bead2 <- factor(allbead.sum$bead, levels=c("DPRbead", "dSibead"), labels =c ("Diproline bead", "dSi bead"))
allbead.fitdata$bead2 <-  factor(allbead.fitdata$bead, levels=c("DPRbead", "dSibead"), labels =c ("Diproline bead", "dSi bead"))


#bw
resize.win(12,9)

ggplot(data=allbead.sum, aes(x=T, y=cellsBase)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=0.5, size=1) +
  geom_smooth(data=allbead.fitdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr), color="black", method="lm", stat="identity", alpha=0.2)+ 
  facet_grid(bead2~condind, labeller=mf_labeller, scale="free") +
  scale_shape_manual (values=plotshapes, name="Treatment") +
  labs(list(x = "Time (min)", y = "Normalized cell count"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.50, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + 
  scale_x_continuous (breaks=c(0, 2, 4, 6, 8, 10))+
  scale_y_continuous(labels=scaleFUN)



#others

ggplot(data=allbead.sum, aes(x=T, y=cellsBase, shape=condind)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=0.5, size=1) +
  geom_smooth(data=allbead.fitdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr), color="black", method="lm", stat="identity", alpha=0.2)+ 
  facet_grid(bead~condind, labeller=mf_labeller, scale="free") +
  scale_shape_manual (values=plotshapes, name="Treatment") +
  labs(list(x = "Time (min)", y = "Normalized cell count"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.50, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + 
  scale_x_continuous (breaks=c(0, 2, 4, 6, 8, 10))+
  scale_y_continuous(labels=scaleFUN)


ggplot(data=allbead.sum, aes(x=T, y=cellsBase, shape=condind, color=condind)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=0.2, size=1) +
  geom_smooth(data=allbead.fitdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, fill=condind), method="lm", stat="identity", alpha=0.2)+ 
  scale_color_manual(values = cbPalette, name="Treatment") + facet_grid(bead~condind, labeller=mf_labeller, scale="free") +
  scale_shape_manual (values=plotshapes, name="Treatment") +
  scale_fill_manual(values=cbPalette, name="Treatment") + 
  labs(list(x = "Time (min)", y = "Normalized cell count"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))
