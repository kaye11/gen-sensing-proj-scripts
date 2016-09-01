#STEP 1
dir<-"d:/Karen's/PhD/R program/General sensing proj/csv files/dSi old data diff starvation days/"


#get dir/file names
ff<-do.call(rbind, lapply(dir, function(x) {
  ff<-list.files(x, "\\.xls$", include.dirs = FALSE, full.names = TRUE)
  data.frame(dir=basename(x), file=basename(ff), 
             fullpath=ff, stringsAsFactors=F)
}))

#read into data.frame

source("readstack.R") #run the code readstack.R
library(plyr)
merged.data= read.stack(ff$fullpath, extra=list(file=ff$file))

##data table challenge
library (data.table)

##direct input from trackmate 
NT=data.table(merged.data, key="Label")
speed.data=NT[, list(day=file, V=TRACK_MEAN_SPEED, Vlog=log(TRACK_MEAN_SPEED+1)), by=c("Label")]

library(ggplot2)
library(grLabel)
library(gtable)
library(ggthemes)
library(grLabelExtra)
library(mgcv)
library(ggplot2)
library(grid)

source("AED.R")
source("vif.R")
source("summarySE.R")
source("resizewin.R")

resize.win(9,6)

#plotting

qplot(day,V, color = day, data = speed.data,  geom = "boxplot")
qplot(day,Vlog, color = day, data = speed.data,  geom = "boxplot")


speed.sumV <- summarySE(speed.data, measurevar="V", groupvars=c("day"))
speed.sumVlog <- summarySE(speed.data, measurevar="Vlog", groupvars=c("day"))

#make another variable
speed.data$ID=as.factor(paste(speed.data$Label, speed.data$day, sep = "-"))


exp=as.data.frame(data.table(cbind(day=as.numeric(as.factor(speed.data$day)), ID=as.numeric(as.factor(speed.data$ID)))))

cor(exp, method = "spearman")
pairs(exp, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

vif_func(in_frame=exp,thresh=5,trace=T)

#boxplots
op=par(mfrow=c(1,2))
boxplot(Vlog~day, data=speed.data)
boxplot(Vlog~ID, data=speed.data)

library(lawstat)
levene.test(speed.data$Vlog, group=speed.data$day, location="mean") #unequal

#nlme
Form <- formula (Vlog ~ day)
si.gls<- gls(Form, data=speed.data)

#aov
si.aov <- aov(Form, data=speed.data)


#nlme model
si1.lme <- lme (Form, random = ~1|ID, method="REML", data=speed.data)

#si2.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|ID), method="REML", data=speed.data)

si3.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|day), method="REML", data=speed.data)


si4.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|day), correlation=corAR1 (), method="REML", data=speed.data) #best

anova(si.gls, si1.lme, si3.lme, si4.lme) 

#best is si3.lme

library(multcomp)

summary(glht(si3.lme, linfct=mcp(day="Tukey", covariate_average=TRUE)))

speed.data$starveday <- factor(speed.data$day, levels=c("0day.xls", "1day.xls", "3days.xls", "7days.xls"), 
                                  labels=c("Day0", "Day1", "Day3", "Day7"))

speed.sumV <- summarySE(speed.data, measurevar="V", groupvars=c("starveday"))


#let's plot this! #plotted that same colors are not statistically significant to each other

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 

ggplot(data=speed.sumV, aes(x=starveday, y=V, fill=starveday)) + geom_bar(stat="identity")+ 
  geom_errorbar(aes(ymin=V-se, ymax=V+se), width=.2, size=1) +
  scale_fill_manual(values = c(Day0="lightgray", Day1="lightgray", Day3="#999999", Day7 ="black"))+
  labs(list(x = "Length of starvation", y = "Mean cell speed (µm/s)"))+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))

