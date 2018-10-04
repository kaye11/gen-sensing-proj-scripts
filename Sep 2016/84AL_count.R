library(ggplot2)
library(grid)
library(gtable)
library(ggthemes)
library(gridExtra)
library(mgcv)
library(data.table)
library(nlme)

source("resizewin.R")
resize.win(12,9)

#read data

count_84AL <- read.csv("D:/Karen's/PhD/R program/General sensing proj/csv files/84AL_count.csv", sep=";")

#condition is just one (-dSi) because its just really one. BEAD is the treatment variable (dSi vs. control bead)

#correct time variable
count_84AL$T= count_84AL$timemin
count_84AL$T.factor= as.factor(count_84AL$T)

#make a new factor
count_84AL$wellvidbead <- as.factor(paste(count_84AL$wellvid, count_84AL$bead, sep = "-"))
count_84AL$wellvidbead <- as.factor(paste(count_84AL$wellvid, count_84AL$bead, sep = "-"))
count_84AL$wellvidcond <- as.factor(paste(count_84AL$wellvid, count_84AL$cond, sep = "-"))

#check initial plot
qplot(T.factor,cells, data = count_84AL,  geom = "boxplot") + facet_grid(bead~., scales="free") 

#standardization

count_84AL$cellsS=NA
k=split(count_84AL, count_84AL$bead)
count_84ALstd <- lapply(k, function (x) scale(x[,c("cells")], center=T, scale=T))
count_84AL$cellsS=unsplit(count_84ALstd, count_84AL$bead)

#standardization without treatments
count_84ALstd2 <- count_84AL
count_84ALstd2$cellsS <- scale(x=count_84ALstd2$cells, center=T, scale=T)

#baselining to 0 at time point 0
NT<-data.table(count_84AL, key=c("wellvidbead"))

t1=NT[,list(bead=bead, T=T, T.factor=T.factor, cells=cells, cellsS=cellsS,
            cellsBase=(cellsS-cellsS[1]), cellsnorm=(cells-cells[1])), by=c("wellvidbead")]

NT2 <- data.table(count_84ALstd2, key=c("wellvidbead"))

t2=NT2[,list(bead=bead, T=T, T.factor=T.factor, cells=cells, cellsS=cellsS,
            cellsBase=(cellsS-cellsS[1]), cellsnorm=(cells-cells[1])), by=c("wellvidbead")]

countbase <- t1 #DATA IS NOW CALLED COUNTBASE

countbase2 <- t2

countbaseorig<- countbase

countbase <- as.data.frame(countbase)
countbase2 <- as.data.frame(countbase2)

source("outlierKD.R")

outlierKD(countbase, cellsBase)
outlierKD(countbase2, cellsBase)

boxplot.stats(countbase$cellsBase)$out


qplot(as.factor(T),cellsBase, data = na.omit(countbase),  geom = "boxplot") + facet_grid(bead~.) +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=1))

qplot(as.factor(T),cellsBase, data = na.omit(countbase2),  geom = "boxplot") + facet_grid(bead~.) +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=1)) ##use countbase2

ggplotly(a)

qplot(as.factor(T),cellsnorm, data = countbase,  geom = "boxplot") + facet_grid(bead~.) +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=1))

source ("summarySE.R")

countbase.sum <- summarySE(countbase, measurevar="cellsBase", groupvars=c("T.factor", "bead"))

countbase.sum.cellsnorm <- summarySE(countbase, measurevar="cellsnorm", groupvars=c("T.factor", "bead"))

ggplot(data=countbase.sum, aes(x=T.factor, y=cellsBase)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=0.5, size=1) + facet_grid(bead~.)

ggplot(data=countbase.sum.cellsnorm, aes(x=T.factor, y=cellsnorm)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsnorm-se, ymax=cellsnorm+se), width=0.5, size=1) + facet_grid(bead~.)

source ("AED.R")
source ("vif.R")
library (lawstat)

expdSi=as.data.frame(data.table(cbind(bead=countbase$bead, T=countbase$T.factor, ID=countbase$wellvidbead)))
cor(expdSi, method = "spearman")

vif_func(in_frame=expdSi,thresh=5,trace=T)

pairs(expdSi, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(countbase$cellsnorm, group=countbase$wellvidbead, location="mean") #significant
levene.test(countbase$cellsnorm, group=countbase$T, location="mean") # significant
levene.test(countbase$cellsnorm, group=countbase$bead, location="mean") #  significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(cellsnorm~bead, data=countbase)
boxplot(cellsnorm~wellvidbead, data=countbase)
boxplot (cellsnorm~T, data=countbase)

#fit a gls
Form <- formula (cellsBase~ bead*T)
countbase.gls<- gls(Form, data=countbase)


#nlme model

#random factor
countbase1.lme <- lme (Form, random = ~1|wellvidbead, method="REML", data=countbase) #BEST

countbase2.lme <- lme (Form, random = ~1|bead, method="REML", data=countbase)

countbase3.lme <- lme (Form, random = ~1|wellvidbead/bead, method="REML", data=countbase)

anova(countbase.gls, countbase1.lme, countbase2.lme, countbase3.lme)

#variance structure

countbase4.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|wellvidbead), method="REML", data=countbase) 

countbase5.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|bead), method="REML", data=countbase) 

#countbase6.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|wellvidbead/bead), method="REML", data=countbase) 

anova(countbase.gls, countbase1.lme, countbase2.lme, countbase3.lme, countbase4.lme, countbase5.lme)

#correlation structures

countbase7.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|bead), correlation=corAR1(), method="REML", data=countbase) 

countbase8.lme <- lme (Form, random = ~1|wellvidbead,   weights=varIdent(form=~1|wellvidbead), correlation=corAR1(), method="REML", data=countbase) 


anova(countbase.gls, countbase1.lme, countbase2.lme, countbase3.lme, countbase4.lme, countbase5.lme, countbase7.lme, countbase8.lme)

#countbase5.lme= 144.7539

summary(countbase8.lme)
anova(countbase8.lme)

library(multcomp)
#multiple comparisons
summary(glht(countbase5.lme, linfct=mcp(bead="Tukey", covariate_average=TRUE, interaction_average = TRUE)))

#residuals
countbase.E2<-resid(countbase8.lme,type="normalized")
countbase.F2<-fitted(countbase8.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=countbase.F2,y=countbase.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(countbase.E2~bead,data=countbase, main="bead",ylab=MyYlab)
plot(x=countbase$T,y=countbase.E2,main="Time",ylab=MyYlab,xlab="Time (min)")
par(op)

xyplot (countbase.E2 ~ T| bead, data=countbase, ylab="Residuals", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})

#fitting data

library(AICcmodavg)

countbase.fit <- as.data.frame(predictSE.lme(countbase8.lme, countbase, se.fit = TRUE, level = 0,
                                             print.matrix = FALSE))

countbase.fit$upr <- countbase.fit$fit + (1.96 * countbase.fit$se)
countbase.fit$lwr <- countbase.fit$fit - (1.96 * countbase.fit$se)

countbase.fit.combdata <- cbind(countbase, countbase.fit)

#summaries
countbase.sum <- summarySE(countbase, measurevar="cellsBase", groupvars=c("T", "bead"))


#plot

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 

scaleFUN <- function(x) sprintf("%.1f", x)

#bw
resize.win(6,8)

ggplot(data=countbase.sum, aes(x=T, y=cellsBase, shape=bead)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=0.5, size=1) +
  geom_smooth(data=countbase.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr), color="black", 
              method="lm", stat="identity", alpha=0.2)+ 
  scale_shape_discrete(name="Treatment") +
  labs(list(x = "Time (min)", y = "Normalized cell count"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text = element_text(size=15), legend.title=element_blank(), legend.text=element_text(size=15), panel.margin=unit (1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + 
  scale_x_continuous (breaks=c(0, 2, 4, 6, 8, 10))+
  scale_y_continuous(labels=scaleFUN)


#for microscale

ggplot(data=countbase.sum, aes(x=T, y=cellsBase, shape=bead, color=bead)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=cellsBase-se, ymax=cellsBase+se), width=0.2, size=1) + 
  geom_smooth(data=countbase.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, fill=bead), 
              method="lm", stat="identity", alpha=0.1)+ 
  scale_colour_manual(values = c("control bead"="lightcoral", "dSi bead"="steelblue2"), name="Treatment") +
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") + 
  labs(list(x = "Time (s)", y = "Normalized cell count", title="Small cells"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20), 
        axis.title.x=element_text(size=20),
        plot.title = element_text(size =20, face="bold"),  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=element_blank(), legend.text=text, panel.spacing=unit (0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(seq(0, 10, 2)))

