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

countchoice_85AS <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/Counts85As_choice.csv", sep=";")

#correct time variable
countchoice_85AS$T= countchoice_85AS$timemin
countchoice_85AS$T.factor= as.factor(countchoice_85AS$T)

#make a new factor
#important is induction because condition is the same for both
countchoice_85AS$indbead <- as.factor(paste(countchoice_85AS$induction, countchoice_85AS$bead, sep = "-"))
countchoice_85AS$wellvidbead <- as.factor(paste(countchoice_85AS$wellvid, countchoice_85AS$bead, sep = "-"))
countchoice_85AS$wellvidind <- as.factor(paste(countchoice_85AS$wellvid, countchoice_85AS$induction, sep = "-"))

#check initial plot
qplot(induction,cells, color = bead, data = countchoice_85AS,  geom = "boxplot")
qplot(T.factor,cells, color = bead, data = countchoice_85AS,  geom = "boxplot") + facet_grid(induction~., scales="free")

#standardization

countchoice_85AS$cellsS=NA
k=split(countchoice_85AS, countchoice_85AS$induction)
countchoice_85ASstd <- lapply(k, function (x) scale(x[,c("cells")], center=T, scale=T))
countchoice_85AS$cellsS=unsplit(countchoice_85ASstd, countchoice_85AS$induction)

#baselining to 0 at time point 0
NT<-data.table(countchoice_85AS, key=c("wellvidbead"))

t1=NT[,list(condition=condition, induction=induction, indbead=indbead, bead=bead, T=T, T.factor=T.factor, cells=cells, cellsS=cellsS,
            cellsBase=(cellsS-cellsS[1]), cellsnorm=(cells-cells[1])), by=c("wellvidbead")]

countbase <- t1 #DATA IS NOW CALLED COUNTBASE

qplot(as.factor(T),cellsBase, color = bead, data = countbase,  geom = "boxplot") + facet_grid(induction~.,  scales="free") +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=bead, fill=bead), alpha=0.2)

qplot(as.factor(T),cellsnorm, color = bead, data = countbase,  geom = "boxplot") + facet_grid(induction~.,  scales="free") +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=bead, fill=bead), alpha=0.2)

source ("summarySE.R")
countbase.sum <- summarySE(countbase, measurevar="cellsBase", groupvars=c("T.factor", "bead", "induction"))

countbase.sum.cellsnorm <- summarySE(countbase, measurevar="cellsnorm", groupvars=c("T.factor", "bead", "induction"))

ggplot(data=countbase.sum, aes(x=T.factor, y=cellsBase, shape=induction, color=induction)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=0.5, size=1) + facet_grid(bead~induction)

ggplot(data=countbase.sum.cellsnorm, aes(x=T.factor, y=cellsnorm, shape=induction, color=induction)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsnorm-se, ymax=cellsnorm+se), width=0.5, size=1) + facet_grid(bead~induction)


##choice: subset induced

source ("AED.R")
source ("vif.R")
library (lawstat)

induced_choice= subset (countbase, induction=="induced", )
expSi_ind=as.data.frame(data.table(cbind(bead=induced_choice$bead, T=induced_choice$T.factor, ID=induced_choice$wellvidbead)))
cor(expSi_ind, method = "spearman")

vif_func(in_frame=expSi_ind,thresh=5,trace=T)

pairs(expSi_ind, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(induced_choice$cellsBase, group=induced_choice$wellvidbead, location="mean") #significant
levene.test(induced_choice$cellsBase, group=induced_choice$T.factor, location="mean") # significant
levene.test(induced_choice$cellsBase, group=induced_choice$bead, location="mean") # significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(cellsBase~bead, data=induced_choice)
boxplot(cellsBase~wellvidbead, data=induced_choice)
boxplot (cellsBase~T.factor, data=induced_choice)

#fit a gls
Form <- formula (cellsBase ~ bead*T)
induced_choice.gls<- gls(Form, data=induced_choice)


#nlme model

#random factor
induced_choice1.lme <- lme (Form, random = ~1|wellvidbead, method="REML", data=induced_choice) #BEST

induced_choice2.lme <- lme (Form, random = ~1|bead, method="REML", data=induced_choice)

induced_choice3.lme <- lme (Form, random = ~1|wellvidbead/bead, method="REML", data=induced_choice)

anova(induced_choice.gls, induced_choice1.lme, induced_choice2.lme, induced_choice3.lme)

#variance structure

induced_choice4.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|wellvidbead), method="REML", data=induced_choice)

induced_choice5.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|bead), method="REML", data=induced_choice) #BEST

#induced_choice6.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|wellvidbead/bead), method="REML", data=induced_choice) 

anova(induced_choice.gls, induced_choice1.lme, induced_choice2.lme, induced_choice3.lme, induced_choice5.lme)

#correlation structures

induced_choice7.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|bead), correlation=corAR1(), method="REML", data=induced_choice) #BEST use this

induced_choice8.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|bead), 
                            correlation=corAR1 (form=~1|wellvidbead), method="REML", data=induced_choice) 

induced_choice9.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|bead), 
                            correlation=corAR1 (form=~1|wellvidbead/bead), method="REML", data=induced_choice) 


anova(induced_choice.gls, induced_choice1.lme, induced_choice2.lme, induced_choice3.lme, induced_choice5.lme, induced_choice7.lme, induced_choice8.lme, induced_choice9.lme)

#induced_choice7-9 same, for consistency use 9

summary(induced_choice9.lme)
anova(induced_choice9.lme)

library(multcomp)
#multiple comparisons
summary(glht(induced_choice9.lme, linfct=mcp(bead="Tukey", covariate_average=TRUE, interaction_average = TRUE)))

#residuals
induced_choice.E2<-resid(induced_choice9.lme,type="normalized")
induced_choice.F2<-fitted(induced_choice9.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=induced_choice.F2,y=induced_choice.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(induced_choice.E2~bead,data=induced_choice, main="bead",ylab=MyYlab)
plot(x=induced_choice$T,y=induced_choice.E2,main="Time",ylab=MyYlab,xlab="Time (min)")
par(op)

xyplot (induced_choice.E2 ~ T| bead, data=induced_choice, ylab="Residuals", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})


##choice: subset not induced

notinduced_choice= subset (countbase, induction=="notinduced", )
expSi_notind=as.data.frame(data.table(cbind(bead=notinduced_choice$bead, T=notinduced_choice$T.factor, ID=notinduced_choice$wellvidbead)))
cor(expSi_notind, method = "spearman")

vif_func(in_frame=expSi_notind,thresh=5,trace=T)

pairs(expSi_notind, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(notinduced_choice$cellsBase, group=notinduced_choice$wellvidbead, location="mean") #significant
levene.test(notinduced_choice$cellsBase, group=notinduced_choice$T.factor, location="mean") # significant
levene.test(notinduced_choice$cellsBase, group=notinduced_choice$bead, location="mean") # significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(cellsBase~bead, data=notinduced_choice)
boxplot(cellsBase~wellvidbead, data=notinduced_choice)
boxplot (cellsBase~T.factor, data=notinduced_choice)

#fit a gls
Form <- formula (cellsBase ~ bead*T)
notinduced_choice.gls<- gls(Form, data=notinduced_choice)


#nlme model

#random factor
notinduced_choice1.lme <- lme (Form, random = ~1|wellvidbead, method="REML", data=notinduced_choice) #BEST

notinduced_choice2.lme <- lme (Form, random = ~1|bead, method="REML", data=notinduced_choice)

notinduced_choice3.lme <- lme (Form, random = ~1|wellvidbead/bead, method="REML", data=notinduced_choice)

anova(notinduced_choice.gls, notinduced_choice1.lme, notinduced_choice2.lme, notinduced_choice3.lme)

#variance structure

notinduced_choice4.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|wellvidbead), method="REML", data=notinduced_choice)

notinduced_choice5.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|bead), method="REML", data=notinduced_choice) #BEST

#notinduced_choice6.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|wellvidbead/bead), method="REML", data=notinduced_choice) 

anova(notinduced_choice.gls, notinduced_choice1.lme, notinduced_choice2.lme, notinduced_choice3.lme, notinduced_choice5.lme)

#correlation structures

notinduced_choice7.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|bead), correlation=corAR1(), method="REML", data=notinduced_choice) #BEST use this

notinduced_choice8.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|bead), 
                               correlation=corAR1 (form=~1|wellvidbead), method="REML", data=notinduced_choice) 

notinduced_choice9.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|bead), 
                               correlation=corAR1 (form=~1|wellvidbead/bead), method="REML", data=notinduced_choice) 


anova(notinduced_choice.gls, notinduced_choice1.lme, notinduced_choice2.lme, notinduced_choice3.lme, notinduced_choice5.lme, notinduced_choice7.lme, notinduced_choice8.lme, notinduced_choice9.lme)

#notinduced_choice7-9 same, for consistency use 9

summary(notinduced_choice9.lme)
anova(notinduced_choice9.lme)

#multiple comparisons
summary(glht(notinduced_choice9.lme, linfct=mcp(bead="Tukey", covariate_average=TRUE, interaction_average = TRUE)))

#residuals
notinduced_choice.E2<-resid(notinduced_choice9.lme,type="normalized")
notinduced_choice.F2<-fitted(notinduced_choice9.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=notinduced_choice.F2,y=notinduced_choice.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(notinduced_choice.E2~bead,data=notinduced_choice, main="bead",ylab=MyYlab)
plot(x=notinduced_choice$T,y=notinduced_choice.E2,main="Time",ylab=MyYlab,xlab="Time (min)")
par(op)

xyplot (notinduced_choice.E2 ~ T| bead, data=notinduced_choice, ylab="Residuals", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})


#fitting data

library(AICcmodavg)

#DPR fit
induced_choice.fit <- as.data.frame(predictSE.lme(induced_choice9.lme, induced_choice, se.fit = TRUE, level = 0,
                                                  print.matrix = FALSE))

induced_choice.fit$upr <- induced_choice.fit$fit + (1.96 * induced_choice.fit$se)
induced_choice.fit$lwr <- induced_choice.fit$fit - (1.96 * induced_choice.fit$se)

induced_choice.fit.combdata <- cbind(induced_choice, induced_choice.fit)

#dSi fit
notinduced_choice.fit <- as.data.frame(predictSE.lme(notinduced_choice9.lme, notinduced_choice, se.fit = TRUE, level = 0,
                                                     print.matrix = FALSE))

notinduced_choice.fit$upr <- notinduced_choice.fit$fit + (1.96 * notinduced_choice.fit$se)
notinduced_choice.fit$lwr <- notinduced_choice.fit$fit - (1.96 * notinduced_choice.fit$se)

notinduced_choice.fit.combdata <- cbind(notinduced_choice, notinduced_choice.fit)

#summaries
induced_choice.sum <- summarySE(induced_choice, measurevar="cellsBase", groupvars=c("T", "bead", "induction"))

notinduced_choice.sum <- summarySE(notinduced_choice, measurevar="cellsBase", groupvars=c("T", "bead", "induction"))


#plot

grid.newpage()
text <- element_text(size = 20, color="black") #change the size of the axes
theme_set(theme_bw()) 

allind.fitdata = rbind (induced_choice.fit.combdata, notinduced_choice.fit.combdata)
allind.sum <- rbind (induced_choice.sum, notinduced_choice.sum)

scaleFUN <- function(x) sprintf("%.1f", x)

allind.sum$bead <- factor(allind.sum$bead, levels=c("DPRbead", "dSibead"), labels =c (" Diproline bead  ", " dSi bead  "))
allind.fitdata$bead <-  factor(allind.fitdata$bead, levels=c("DPRbead", "dSibead"), labels =c (" Diproline bead  ", " dSi bead  "))

label_metrics <- function(x){
  x[x=="induced"] <- "induced"
  x[x=="notinduced"]   <- "not induced"
  x
}

mf_labeller <- ggplot2::as_labeller(label_metrics)

#bw
resize.win(7,8)

ggplot(data=allind.sum, aes(x=T, y=cellsBase, shape=bead)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=0.5, size=1) +
  geom_point(size=5, color='black') +
  geom_ribbon(data=allind.fitdata, aes(ymin=lwr, ymax=upr, linetype=NA), 
              stat="identity", alpha=0.2)+ 
  geom_line (data=allind.fitdata, size=1, aes(y=fit)) +
  facet_grid(~induction, labeller=mf_labeller) +
  labs(list(x = "Time (min)", y = "Normalized cell count",  title="Attraction of dSi-starved cells to beads"))+ 
  theme(axis.title = text,
        plot.title = element_text(size =24, hjust=0.5), axis.text=text,  legend.position="bottom", legend.title=element_blank(),
        strip.text= text, legend.text=text, 
        panel.grid.major = element_blank(), panel.spacing = unit (0.5, "lines") , 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_x_continuous (breaks=c(0, 2, 4, 6, 8, 10))+
  scale_y_continuous(labels=scaleFUN)

ggplot(data=allind.sum, aes(x=T, y=cellsBase, shape=bead)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=0.5, size=1) +
  geom_point(size=5, color='black') +
  geom_ribbon(data=allind.fitdata, aes(ymin=lwr, ymax=upr, linetype=NA), 
              stat="identity", alpha=0.2)+ 
  geom_line (data=allind.fitdata, size=1, aes(y=fit)) +
  facet_grid(induction~., labeller=mf_labeller) +
  labs(list(x = "Time (min)", y = "Normalized cell count",  title="Attraction of dSi-starved cells to beads"))+ 
  theme(axis.title = text,
        plot.title = element_text(size =24, hjust=0.5), axis.text=text,  legend.position="bottom", legend.title=element_blank(),
        strip.text= text, legend.text=text, 
        panel.grid.major = element_blank(), panel.spacing = unit (0.5, "lines") , 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_x_continuous (breaks=c(0, 2, 4, 6, 8, 10))+
  scale_y_continuous(labels=scaleFUN)

##colored

ggplot(data=allind.sum, aes(x=T, y=cellsBase, color=bead)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=0.5, size=1, color="black") +
  geom_point(size=5, shape=21, color='black', aes(fill=bead)) +
  geom_ribbon(data=allind.fitdata, aes(ymin=lwr, ymax=upr, linetype=NA, fill=bead), 
              stat="identity", alpha=0.2)+ 
  geom_line (data=allind.fitdata, size=1, aes(y=fit)) +
  facet_grid(induction~., labeller=mf_labeller) +
  scale_color_manual(values = c("#E69F00", "steelblue2")) +
  scale_fill_manual(values = c("#E69F00", "steelblue2")) +
  labs(list(x = "Time (min)", y = "Normalized cell count", title="Attraction of dSi-starved cells to beads"))+  
  theme(axis.title = text,
        plot.title = element_text(size =24, hjust=0.5), axis.text=text,  legend.position="bottom", legend.title=element_blank(),
        strip.text= text, legend.text=text, 
        panel.grid.major = element_blank(), panel.spacing = unit (0.5, "lines") , 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_x_continuous (breaks=c(0, 2, 4, 6, 8, 10))+
  scale_y_continuous(labels=scaleFUN)
