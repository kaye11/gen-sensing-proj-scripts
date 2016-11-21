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

#standardization without treatments
count_84ALstd2 <- count_84AL
count_84ALstd2$cellsS <- scale(x=count_84ALstd2$cells, center=T, scale=T)

#baselining to 0 at time point 0

NT2 <- data.table(count_84ALstd2, key=c("wellvidbead"))

t2=NT2[,list(bead=bead, T=T, T.factor=T.factor, cells=cells, cellsS=cellsS,
            cellsBase=(cellsS-cellsS[1]), cellsnorm=(cells-cells[1])), by=c("wellvidbead")]

#DATA IS NOW CALLED COUNTBASE

countbase2 <- t2

countbaseorig<- countbase2 #keep original data (no transformation and removal of outliers)

countbase2 <- as.data.frame(countbase2)

source("outlierKD.R")

outlierKD(countbase2, cellsBase)

boxplot.stats(countbase2$cellsBase)$out

qplot(as.factor(T),cellsBase, data = na.omit(countbase2),  geom = "boxplot") + facet_grid(bead~.) +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=1)) ##use countbase2


source ("summarySE.R")

countbase.sum<- ddply(countbase2, c("T.factor", "T", "bead"), summarise,
                  N    = length(cellsBase),
                  mean = mean(cellsBase, na.rm=TRUE),
                  sd   = sd(cellsBase, na.rm=TRUE),
                  se   = sd / sqrt(N))


countbase.sum <- summarySE(countbase2, measurevar="cellsBase", groupvars=c("T.factor", "bead"), na.omit=TRUE)

ggplot(data=na.omit(countbase.sum), aes(x=T.factor, y=mean)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5, size=1) + facet_grid(bead~.)


source ("AED.R")
source ("vif.R")
library (lawstat)

expdSi=as.data.frame(data.table(cbind(bead=countbase2$bead, T=countbase2$T.factor, ID=countbase2$wellvidbead)))
cor(expdSi, method = "spearman")

vif_func(in_frame=expdSi,thresh=5,trace=T)

pairs(expdSi, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(countbase2$cellsnorm, group=countbase2$wellvidbead, location="mean") #significant
levene.test(countbase2$cellsnorm, group=countbase2$T, location="mean") # significant
levene.test(countbase2$cellsnorm, group=countbase2$bead, location="mean") #  significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(cellsnorm~bead, data=countbase2)
boxplot(cellsnorm~wellvidbead, data=countbase2)
boxplot (cellsnorm~T, data=countbase2)

#fit a gls
Form <- formula (cellsBase~ bead*T)
countbase2.gls<- gls(Form, data=countbase2, na.action=na.exclude)


#nlme model

#random factor
countbase21.lme <- lme (Form, random = ~1|wellvidbead, method="REML", data=countbase2, na.action=na.exclude) #BEST

countbase22.lme <- lme (Form, random = ~1|bead, method="REML", data=countbase2, na.action=na.exclude)

countbase23.lme <- lme (Form, random = ~1|wellvidbead/bead, method="REML", data=countbase2, na.action=na.exclude)

anova(countbase2.gls, countbase21.lme, countbase22.lme, countbase23.lme)

#variance structure

countbase24.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|wellvidbead), method="REML", 
                        data=countbase2, na.action=na.exclude) 

countbase25.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|bead), method="REML", 
                        na.action=na.exclude, data=countbase2) #BEST

#countbase26.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|wellvidbead/bead), method="REML", na.action=na.exclude, data=countbase2) 

anova(countbase2.gls, countbase21.lme, countbase22.lme, countbase23.lme, countbase24.lme, countbase25.lme)

#correlation structures

countbase27.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|bead), correlation=corAR1(), method="REML", 
                        na.action=na.exclude, data=countbase2) 

countbase28.lme <- lme (Form, random = ~1|wellvidbead,   weights=varIdent(form=~1|wellvidbead), correlation=corAR1(), method="REML", 
                        na.action=na.exclude, data=countbase2) 


anova(countbase2.gls, countbase21.lme, countbase22.lme, countbase23.lme, countbase24.lme, countbase25.lme, countbase27.lme, countbase28.lme)

#countbase25.lme= 35.81855

summary(countbase25.lme)
anova(countbase25.lme)

library(multcomp)
#multiple comparisons
summary(glht(countbase25.lme, linfct=mcp(bead="Tukey", covariate_average=TRUE, interaction_average = TRUE)))

#residuals
countbase2.E2<-resid(countbase25.lme,type="normalized")
countbase2.F2<-fitted(countbase25.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=countbase2.F2,y=countbase2.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(countbase2.E2~bead,data=countbase2, main="bead",ylab=MyYlab)
plot(x=countbase2$T,y=countbase2.E2,main="Time",ylab=MyYlab,xlab="Time (min)")
par(op)

xyplot (countbase2.E2 ~ T| bead, data=countbase2, ylab="Residuals", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})

#fitting data

library(AICcmodavg)

countbase2.fit <- as.data.frame(predictSE.lme(countbase25.lme, countbase2, se.fit = TRUE, level = 0,
                                              print.matrix = FALSE))

countbase2.fit$upr <- countbase2.fit$fit + (1.96 * countbase2.fit$se)
countbase2.fit$lwr <- countbase2.fit$fit - (1.96 * countbase2.fit$se)

countbase2.fit.combdata <- cbind(countbase2, countbase2.fit)

#summaries
#use countbase.sum

#plot

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 

scaleFUN <- function(x) sprintf("%.1f", x)

#bw
resize.win(6,8)

ggplot(data=countbase.sum, aes(x=T, y=mean, shape=bead)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5, size=1) +
  geom_smooth(data=countbase2.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr), color="black", 
              method="lm", stat="identity", alpha=0.2)+ 
  scale_shape_discrete(name="Treatment") +
  labs(list(x = "Time (min)", y = "Normalized cell count"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text = element_text(size=15), legend.title=element_blank(), legend.text=element_text(size=15), panel.margin=unit (1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) 



#for microscale

resize.win(6,7)

ggplot(data=countbase.sum, aes(x=T, y=mean, shape=bead, color=bead)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1) + 
  geom_smooth(data=countbase2.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, fill=bead), 
              method="lm", stat="identity", alpha=0.1)+ 
  scale_colour_manual(values = c("control bead"="lightcoral", "dSi bead"="steelblue2"), name="Treatment") +
  scale_shape_discrete (name="Treatment") +
  scale_fill_discrete(name="Treatment") + 
  labs(list(x = "Time (min)", y = "Normalized cell count", title="Large cells"))+ 
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face="bold"), 
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=element_blank(), legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) +scale_x_continuous (breaks=c(seq(0, 10, 2)))

