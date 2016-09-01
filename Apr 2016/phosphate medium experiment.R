#STEP 1
dirdata<-"d:/Karen's/PhD/R program/General sensing proj/csv files/phosphate motility data/"


#get dir/file names
ff<-do.call(rbind, lapply(dirdata, function(x) {
  ff<-list.files(x, "\\.xls$", include.dirs = FALSE, full.names = TRUE)
  data.frame(dir=basename(x), file=basename(ff), 
             fullpath=ff, stringsAsFactors=F)
}))

#read into data.frame

source("readstack.R") #run the code readstack.R
library(plyr)
merged.data= read.stack(ff$fullpath, extra=list(file=ff$file))

library(plyr)

#read into data.frame

pattern <- "(\\w{1})(\\d{1})(_)(\\w+)(.xls)"

well <- gsub(pattern,'\\1\\2',merged.data$file)
treatment <- gsub(pattern,'\\4',merged.data$file)

final.data=cbind(merged.data, well, treatment)
final.data=subset(merged.data2, merged.data2$NUMBER_SPOTS>10, )

##data table challenge
library (data.table)

##direct input from trackmate 
NT=data.table(final.data, key="Label")
speed.data=NT[, list(treatment=treatment, well=well, V=TRACK_MEAN_SPEED, Vlog=log(TRACK_MEAN_SPEED+1)), by=c("Label")]

library(ggplot2)
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

qplot(treatment,V, color = treatment, data = speed.data,  geom = "boxplot")
qplot(treatment,Vlog, color = treatment, data = speed.data,  geom = "boxplot")

speed.sumV <- summarySE(speed.data, measurevar="V", groupvars=c("treatment"))
speed.sumVlog <- summarySE(speed.data, measurevar="Vlog", groupvars=c("treatment"))

#make another variable
speed.data$ID=as.factor(paste(speed.data$well, speed.data$Label, sep = "-"))
speed.data$wellcond = as.factor(paste(speed.data$well, speed.data$treatment, sep = "-"))

#test correlation
exp=as.data.frame(data.table(cbind(treatment=as.numeric(as.factor(speed.data$treatment)), wellcond=as.numeric(as.factor(speed.data$wellcond)),
                                   ID=as.numeric(as.factor(speed.data$ID)))))

cor(exp, method = "spearman")
pairs(exp, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

vif_func(in_frame=exp,thresh=5,trace=T)

#boxplots
op=par(mfrow=c(2,2))
boxplot(Vlog~treatment, data=speed.data)
boxplot(Vlog~ID, data=speed.data)
boxplot (Vlog~wellcond, data=speed.data)

library(lawstat)
levene.test(speed.data$Vlog, group=speed.data$treatment, location="mean") #unequal
levene.test(speed.data$Vlog, group=speed.data$ID, location="mean") #unequal
levene.test(speed.data$Vlog, group=speed.data$wellcond, location="mean") #unequal

#transform data and try
speed.data$Vlog2 = log((2*speed.data$V)+1)
speed.data$Vlog3 = log(speed.data$V+3/8)
speed.data$Vsqrt = sqrt(speed.data$V)
speed.data$Vsqrd= (speed.data$V)^2
speed.data$Vlog10 = log((2*speed.data$V)+1)
#verdict: Vlog which is log+1 transformation is still the best

#nlme
Form <- formula (Vlog~ treatment)
pho.gls<- gls(Form, data=speed.data)

#aov
pho.aov <- aov(Form, data=speed.data)


#nlme model
pho1.lme <- lme (Form, random = ~1|ID, method="REML", data=speed.data)

#variance structures, varIdent
pho2.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|well), method="REML", data=speed.data)

pho3.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|wellcond), method="REML", data=speed.data)

pho4.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), method="REML", data=speed.data) #lowest

#pho5.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|ID), method="REML", data=speed.data) 

anova(pho1.lme, pho2.lme, pho3.lme, pho4.lme)


#try other variance structure (DID NOT WORK! huhubels)
pho4a.lme <- lme (Form, random = ~1|ID,  weights = varExp(form=~fitted(.)), method="REML", data=speed.data) 
pho4b.lme <- lme (Form, random = ~1|ID,  weights = varExp(form= ~ treatment),  method="REML", data=speed.data) 
pho4c.lme <- lme (Form, random = ~1|ID,  weights = varPower(form= ~ wellcond|ID),  method="REML", data=speed.data) 
pho4d.lme <- lme (Form, random = ~1|ID,  weights = varConstPower(form= ~ treatment|well),  method="REML", data=speed.data) 


#combine variance structures

#pho6.lme <- lme (Form, random = ~1|ID,  weights=varComb(varIdent(form=~1|well), varIdent(form=~1|treatment)), method="REML", data=speed.data) 

pho7.lme <- lme (Form, random = ~1|ID,  weights=varComb(varIdent(form=~1|wellcond), varIdent(form=~1|treatment)), method="REML", data=speed.data) 

anova(pho1.lme, pho2.lme, pho3.lme, pho4.lme, pho4a.lme, pho7.lme)


#correlation structures
pho8.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1(), method="REML", data=speed.data) 
pho9.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1 (form=~1|ID/treatment), method="REML", data=speed.data) 
pho10.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1 (form=~1|ID/wellcond), method="REML", data=speed.data) 

anova(pho1.lme, pho2.lme, pho3.lme, pho4.lme, pho4a.lme, pho7.lme, pho8.lme, pho9.lme, pho10.lme)
#all 3 are the same

#other correlation structures

cs1 <- corARMA(c(0.1), p = 1, q = 0)

pho11.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation = cs1,
                 method="REML", data=speed.data) 

anova(pho1.lme, pho2.lme, pho3.lme, pho4.lme, pho4a.lme, pho7.lme, pho8.lme, pho9.lme, pho10.lme, pho11.lme)

#lowest AIC= 4226.896 is pho4.lme 


#try gamm

pho0 <- gamm (Vlog~s(treatment, bs="fs"), method="REML", data=speed.data) 
#smoothing splines for factors can are 3 types cr (cubic regression), cs (shrinkage version of cr) and cc (cyclic cubic)
pho1 <- gamm (Vlog~s(treatment, bs="fs", xt="cr"), method="REML", data=speed.data) #best
pho2 <- gamm (Vlog~s(treatment, bs="fs", xt="cs"), method="REML", data=speed.data)
pho2.1 <- gamm (Vlog~s(treatment, bs="fs", xt="cc"), method="REML", data=speed.data)

anova(pho0$lme, pho1$lme, pho2$lme, pho2.1$lme)
#everything is the same actually, AIC=4257.09 so LME is still better




library(multcomp)

summary(glht(pho4.lme, linfct=mcp(treatment="Tukey", covariate_average=TRUE)))
summary(glht(pho8.lme, linfct=mcp(treatment="Tukey", covariate_average=TRUE)))

speed.data$treat2 <- factor(speed.data$treatment, levels=c("before", "afterASW",  "afternoP"))

speed.sumV <- summarySE(speed.data, measurevar="V", groupvars=c("treat2"))

#let's plot this! #plotted that same colors are not statistically significant to each other

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 

ggplot(data=speed.sumV, aes(x=treat2, y=V, fill=treat2)) + geom_bar(stat="identity")+ 
  geom_errorbar(aes(ymin=V-se, ymax=V+se), width=.2, size=1) +
  scale_fill_manual(values = c(afternoP="lightgray", before="lightgray", afterASW ="black"))+
  labs(list(y = "Mean cell speed (µm/s)"))+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_blank(),
        plot.title = element_text(size =20, face="bold"), axis.text.x=element_blank(),  legend.position="none",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (2, "lines"),
        panel.grid.major = element_blank(),axis.ticks=element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))+
  annotate(geom = "text", x=1, y=-0.15, label="before dPI \naddition", size=6.8, fontface=2)+
  annotate(geom = "text", x=2, y=-0.15, label="+dPI", size=6.8, fontface=2)+
  annotate(geom = "text", x=3, y=-0.15, label="blank addition", size=6.8, fontface=2)


