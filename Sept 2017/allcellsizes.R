library(ggplot2)
library(grid)
library(gtable)
library(ggthemes)
library(gridExtra)
library(mgcv)
library(data.table)
library(nlme)


#check if accummulation is dependent on cell size

source("resizewin.R")
resize.win(12,9)

library(readr)
allcellsizes <- read_delim("D:/Karen's/PhD/R program/General sensing proj/csv files/allcellsizes.csv", 
                           ";", escape_double = FALSE, trim_ws = TRUE)
View(allcellsizes)

library(ggplot2)

qplot(as.factor(T),cellsBase, color = size, data = allcellsizes,  geom = "boxplot") +  stat_smooth (method="loess", formula=y~x, size=1, aes(group=size))


source ("summarySE.R")

allcellsizes.sum<- ddply(allcellsizes, c("T.factor", "T", "size"), summarise,
                      N    = length(cellsBase),
                      mean = mean(cellsBase, na.rm=TRUE),
                      sd   = sd(cellsBase, na.rm=TRUE),
                      se   = sd / sqrt(N))

allcellsizes.sum.cells <- summarySE(allcellsizes, measurevar="cells", groupvars=c("T.factor", "size"))

ggplot(data=allcellsizes.sum, aes(x=T.factor, y=mean, shape=size, color=size)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5, size=1) 

ggplot(data=allcellsizes.sum.cells, aes(x=T.factor, y=mean, shape=size, color=size)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5, size=1) 


source ("AED.R")
source ("vif.R")
library (lawstat)

expdSi=as.data.frame(data.table(cbind(bead=allcellsizes$size, T=allcellsizes$T.factor, ID=allcellsizes$wellvidbead)))
cor(expdSi, method = "spearman")

vif_func(in_frame=expdSi,thresh=5,trace=T)

pairs(expdSi, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(allcellsizes$cellsBase, group=allcellsizes$wellvidbead, location="mean") #significant
levene.test(allcellsizes$cellsBase, group=allcellsizes$T, location="mean") # significant
levene.test(allcellsizes$cellsBase, group=allcellsizes$size, location="mean") #  not significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(cellsBase~size, data=allcellsizes)
boxplot(cellsBase~wellvidbead, data=allcellsizes)
boxplot (cellsBase~T, data=allcellsizes)

library(nlme)
#fit a gls
Form <- formula (cellsBase~ size*T)
allcellsizes.gls<- gls(Form, data=allcellsizes, na.action=na.exclude)


#nlme model

#random factor
allcellsizes1.lme <- lme (Form, random = ~1|wellvidbead, method="REML", data=allcellsizes, na.action=na.exclude) #BEST

allcellsizes2.lme <- lme (Form, random = ~1|size, method="REML", data=allcellsizes, na.action=na.exclude)

allcellsizes3.lme <- lme (Form, random = ~1|wellvidbead/size, method="REML", data=allcellsizes, na.action=na.exclude)

anova(allcellsizes.gls, allcellsizes1.lme, allcellsizes2.lme, allcellsizes3.lme)

#variance structure

allcellsizes4.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|wellvidbead), method="REML", 
                          data=allcellsizes, na.action=na.exclude) #BEST

allcellsizes5.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|size), method="REML", 
                          na.action=na.exclude, data=allcellsizes) 

#allcellsizes6.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|wellvidbead/size), method="REML", na.action=na.exclude, data=allcellsizes) 

anova(allcellsizes.gls, allcellsizes1.lme, allcellsizes2.lme, allcellsizes3.lme, allcellsizes4.lme, allcellsizes5.lme)

#correlation structures

allcellsizes7.lme <- lme (Form, random = ~1|wellvidbead,  weights=varIdent(form=~1|size), correlation=corAR1(), method="REML", 
                          na.action=na.exclude, data=allcellsizes) 

allcellsizes8.lme <- lme (Form, random = ~1|wellvidbead,   weights=varIdent(form=~1|wellvidbead), correlation=corAR1(), method="REML", 
                          na.action=na.exclude, data=allcellsizes) 


anova(allcellsizes.gls, allcellsizes1.lme, allcellsizes2.lme, allcellsizes3.lme, allcellsizes4.lme, allcellsizes5.lme, allcellsizes7.lme, allcellsizes8.lme)

#allcellsizes8.lme= 81.62901

summary(allcellsizes8.lme)
anova(allcellsizes8.lme)

#no significant differences between sizes

#residuals
allcellsizes.E2<-resid(allcellsizes8.lme,type="normalized")
allcellsizes.F2<-fitted(allcellsizes8.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=allcellsizes.F2,y=allcellsizes.E2,xlab="Fitted values",ylab=MyYlab)
boxplot(allcellsizes.E2~size,data=allcellsizes, main="size",ylab=MyYlab)
plot(x=allcellsizes$T,y=allcellsizes.E2,main="Time",ylab=MyYlab,xlab="Time (min)")
par(op)

xyplot (allcellsizes.E2 ~ T| size, data=allcellsizes, ylab="Residuals", xlab="Time (min)", 
        panel=function(x,y){
          panel.grid(h=-1, v= 2)
          panel.points(x,y,col=1)
          panel.loess(x,y,span=0.5,col=1,lwd=2)})

#fitting data

library(AICcmodavg)

allcellsizes.fit <- as.data.frame(predictSE.lme(allcellsizes8.lme, allcellsizes, se.fit = TRUE, level = 0,
                                                print.matrix = FALSE))

allcellsizes.fit$upr <- allcellsizes.fit$fit + (1.96 * allcellsizes.fit$se)
allcellsizes.fit$lwr <- allcellsizes.fit$fit - (1.96 * allcellsizes.fit$se)

allcellsizes.fit.combdata <- cbind(allcellsizes, allcellsizes.fit)

#summaries
#use allcellsizes.sum

#plot

grid.newpage()
text <- element_text(size = 20, color="black") #change the size of the axes
theme_set(theme_bw()) 

scaleFUN <- function(x) sprintf("%.1f", x)

allcellsizes.sum$size <- factor(allcellsizes.sum$size, levels=c("large", "med", "small"), labels =c ("Large", "Medium", "Small"))

allcellsizes$size <- factor(allcellsizes$size, levels=c("large", "med", "small"), labels =c ("Large", "Medium", "Small"))

allcellsizes.fit.combdata$size <- factor(allcellsizes.fit.combdata$size, levels=c("large", "med", "small"), labels =c ("Large", "Medium", "Small"))


#bw
resize.win(6,6)

ggplot(data=allcellsizes.sum, aes(x=T.factor, y=mean, shape=size)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5, size=1) +
  geom_smooth(data=allcellsizes.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr), color="black", 
              method="lm", stat="identity", alpha=0.2)+ 
  scale_shape_discrete(name="Size") +
  labs(list(x = "Time (min)", y = "Normalized cell count"))+ 
  theme(axis.text=element_text(size=20, colour="black"), axis.title.y=element_text(size=20, vjust=1.5), 
        axis.title.x=element_text(size=20, vjust=-0.5),
        plot.title = element_text(size =24), legend.position="bottom", legend.title=element_blank(),
        strip.text.x = text, strip.text.y = text, legend.text=text, panel.spacing=unit (0.5, "lines"),
        panel.grid.major = element_blank(), panel.spacing.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) 



#for microscale and paper supplementary figure

resize.win(6,7)

ggplot(data=allcellsizes.sum, aes(x=T.factor, y=mean, shape=size, color=size)) + geom_point(size=5)+ 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1) + 
  geom_smooth(data=allcellsizes.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, fill=size), 
              method="lm", stat="identity", alpha=0.1)+ 
  scale_colour_discrete(name="Size") +
  scale_shape_discrete (name="Size") +
  scale_fill_discrete(name="Size") + 
  labs(list(x = "Time (min)", y = "Normalized cell count"))+ 
  theme(axis.text = text,  
        plot.title = element_text(size =24, hjust=0.5), axis.title=element_text (size=20, face="bold"),  legend.position="bottom", 
        strip.text.x = text, strip.text.y = text, legend.title=element_blank(), legend.text=text, panel.spacing=unit (0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))

