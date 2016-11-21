#STEP 1
dirdata<-"d:/Karen's/PhD/R program/General sensing proj/csv files/phosphate motility data no dSi/"


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

pattern <- "(\\w{1})(\\d{1})(-)(\\w+)(-)(\\w+)(-)(\\d+)(.xls)"

well <- gsub(pattern,'\\1\\2',merged.data$file)
wellvid <- gsub(pattern,'\\1\\2\\7\\8',merged.data$file)
treatment <- gsub(pattern,'\\4\\5\\6',merged.data$file)

merged.data2=cbind(merged.data, well, wellvid, treatment)

#drop tracks with less than 5 spots
final.data=subset(merged.data2, merged.data2$NUMBER_SPOTS>29, )

##data table challenge
library (data.table)

##direct input from trackmate 
NT=data.table(final.data, key="Label")
speed.data=NT[, list(treatment=treatment, well=well, wellvid=wellvid, V=TRACK_MEAN_SPEED, Vlog=log(TRACK_MEAN_SPEED+1)), by=c("Label")]

library(ggplot2)
library(gtable)
library(ggthemes)
library(mgcv)
library(ggplot2)
library(grid)

source("AED.R")
source("vif.R")
source("summarySE.R")
source("resizewin.R")

resize.win(12,9)

#plotting

qplot(treatment,V, color = treatment, data = speed.data,  geom = "boxplot")
qplot(treatment,Vlog, color = treatment, data = speed.data,  geom = "boxplot")

speed.sumV <- summarySE(speed.data, measurevar="V", groupvars=c("treatment"))
speed.sumVlog <- summarySE(speed.data, measurevar="Vlog", groupvars=c("treatment"))

#make another variable
speed.data$ID=as.factor(paste(speed.data$wellvid, speed.data$Label, sep = "-"))
speed.data$wellcond = as.factor(paste(speed.data$wellvid, speed.data$treatment, sep = "-"))

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


#nlme
Form <- formula (Vlog~ treatment)
pmot.gls<- gls(Form, data=speed.data)

#aov
pmot.aov <- aov(Form, data=speed.data)


#nlme model
pmot1.lme <- lme (Form, random = ~1|ID, method="REML", data=speed.data)

#variance structures, varIdent
pmot2.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|well), method="REML", data=speed.data)

#pmot3.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|wellcond), method="REML", data=speed.data)

pmot4.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), method="REML", data=speed.data) #lowest

#pmot5.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|ID), method="REML", data=speed.data) 

anova(pmot1.lme, pmot2.lme, pmot4.lme)


#try other variance structure (DID NOT WORK! huhubels)
pmot4a.lme <- lme (Form, random = ~1|ID,  weights = varExp(form=~fitted(.)), method="REML", data=speed.data)

#pmot4b.lme <- lme (Form, random = ~1|ID,  weights = varExp(form= ~ treatment),  method="REML", data=speed.data) 
#pmot4c.lme <- lme (Form, random = ~1|ID,  weights = varPower(form= ~ wellcond|ID),  method="REML", data=speed.data) 
#pmot4d.lme <- lme (Form, random = ~1|ID,  weights = varConstPower(form= ~ treatment|well),  method="REML", data=speed.data) 


#combine variance structures

#pmot6.lme <- lme (Form, random = ~1|ID,  weights=varComb(varIdent(form=~1|treatment), varExp(form=~fitted(.))), method="REML", data=speed.data) 

#pmot7.lme <- lme (Form, random = ~1|ID,  weights=varComb(varIdent(form=~1|wellcond), varIdent(form=~1|treatment)), method="REML", data=speed.data) 

anova(pmot1.lme, pmot2.lme, pmot4.lme, pmot4a.lme)


#correlation structures
pmot8.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1(), method="REML", data=speed.data) 
pmot8a.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1(form=~1|ID), method="REML", data=speed.data) 
#pmot8b.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1(form=~1|treatment), method="REML", data=speed.data) 
#pmot8c.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1(form=~1|wellcond), method="REML", data=speed.data) 

pmot8d.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1(-0.2, form=~1|ID), method="REML", data=speed.data) 
pmot8e.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1(0.6, form=~1|ID), method="REML", data=speed.data) 



pmot9.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1 (form=~1|ID/treatment), method="REML", data=speed.data) 
pmot10.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1 (form=~1|ID/wellcond), method="REML", data=speed.data) 

anova(pmot1.lme, pmot2.lme, pmot4.lme, pmot4a.lme, pmot8.lme, pmot8a.lme, pmot8d.lme, pmot8e.lme, pmot9.lme, pmot10.lme)

#all 3 are the same

#other correlation structures

cs1 <- corARMA(c(0.2, 0.2), p = 2, q = 0)

pmot11.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation = cs1,
                 method="REML", data=speed.data) 

anova(pmot1.lme, pmot2.lme, pmot4.lme, pmot4a.lme, pmot8.lme, pmot8a.lme, pmot8d.lme, pmot8e.lme, pmot9.lme, pmot10.lme, pmot11.lme)

summary(pmot4.lme)
anova(pmot4.lme)

#lowest AIC is pmot4.lme=1237.500

library(multcomp)
summary(glht(pmot4.lme, linfct=mcp(treatment="Tukey", covariate_average=TRUE)))

#other comparison methods, same results
#summary(glht(pmot4.lme, linfct=mcp(treatment="Tukey", covariate_average=TRUE)), test = adjusted(type = "bonferroni"))

library(lsmeans)
pmot.lsmeans=lsmeans(pmot4.lme,  ~ treatment)
summary(pmot.lsmeans, adjust="bonf")
contrast(pmot.lsmeans, adjust="bonf")
pairs(pmot.lsmeans)

library(multcompView)
cld(pmot.lsmeans, alpha=0.05)


#graphing
speed.sumV$before= c("dP-rich medium", "dP-rich medium", "low-dP medium", "blank addition", "low-dP medium")
speed.sumV$after= c("addASW", "noadd", "addASW", "blankadd", "noadd")

speed.sumV$beforelabel <- factor(speed.sumV$before, levels=c("dP-rich medium", "low-dP medium", "blank addition"),
                                 labels= c(" non-starved  ", " dP-starved  ", " blank addition")) 

speed.sumV <- data.table (speed.sumV)
speed.sumV [, grouping := ifelse(treatment %in% c("nonstarved-motility", "starved-motility"), "before medium exchange",
                                 ifelse(treatment %in% c("nonstarved-addASW", "starved-addASW", "starved-blankadd"), 
                                        "1h after dP/ blank addition", NA))]


speed.sumV$groupre <- factor(speed.sumV$grouping, levels=c("before medium exchange", "1h after dP/ blank addition"))

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 

resize.win(12,9)

ggplot(speed.sumV, aes(x=treatment, y=V, fill = beforelabel)) + 
  geom_bar(stat="identity", position = "dodge", width=0.8, color="black") + 
  geom_errorbar(aes(ymin=V-se, ymax=V+se), size=1, width=0.2, position=position_dodge(0.9))+
  facet_grid(~groupre, switch = "x", scales = "free_x", space = "free_x") +
  scale_fill_grey(start = 0.3, end = 1) + 
  labs(y = expression("Mean cell speed"~("µm"~s^-1)), title="dP medium exchange experiment") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.title.y=text, 
        axis.title.x=element_blank(),
        plot.title = element_text(size =24), axis.text=text,  legend.position="bottom", legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, legend.direction="horizontal", 
        panel.margin=unit (0, "lines"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(1,1,1,1), "cm"))