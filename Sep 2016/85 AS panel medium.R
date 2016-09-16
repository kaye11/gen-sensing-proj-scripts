#STEP 1
dirdata<-"d:/Karen's/PhD/R program/General sensing proj/csv files/motility 85AS panel/"


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

pattern <- "(\\w{1})(\\d{1})(_)(\\w+)(_)(\\w+)(-)(\\d+)(.xls)"

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
levene.test(speed.data$Vlog, group=speed.data$ID, location="mean") #didnt work
levene.test(speed.data$Vlog, group=speed.data$wellcond, location="mean") #unequal


#nlme
Form <- formula (Vlog~ treatment)
mot85AS.gls <- gls(Form, data=speed.data)

#aov
mot85AS.aov <- aov(Form, data=speed.data)


#nlme model
mot85AS1.lme <- lme (Form, random = ~1|ID, method="REML", data=speed.data)

#variance structures, varIdent
mot85AS2.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|well), method="REML", data=speed.data)

#mot85AS3.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|wellcond), method="REML", data=speed.data)

mot85AS4.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), method="REML", data=speed.data) #lowest

#mot85AS5.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|ID), method="REML", data=speed.data) 

anova(mot85AS1.lme, mot85AS2.lme, mot85AS4.lme)


#try other variance structure (DID NOT WORK! huhubels)
#mot85AS4a.lme <- lme (Form, random = ~1|ID,  weights = varExp(form=~fitted(.)), method="REML", data=speed.data)

#mot85AS4b.lme <- lme (Form, random = ~1|ID,  weights = varExp(form= ~ treatment),  method="REML", data=speed.data) 
#mot85AS4c.lme <- lme (Form, random = ~1|ID,  weights = varPower(form= ~ wellcond|ID),  method="REML", data=speed.data) 
#mot85AS4d.lme <- lme (Form, random = ~1|ID,  weights = varConstPower(form= ~ treatment|well),  method="REML", data=speed.data) 

#correlation structures
mot85AS5.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1(), method="REML", data=speed.data) 
mot85AS5a.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1(form=~1|ID), method="REML", data=speed.data) 
#mot85AS5b.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1(form=~1|treatment), method="REML", data=speed.data) 
#mot85AS5c.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1(form=~1|wellcond), method="REML", data=speed.data) 

mot85AS5d.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1(-0.2, form=~1|ID), method="REML", data=speed.data) 
mot85AS5e.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1(0.6, form=~1|ID), method="REML", data=speed.data) 

mot85AS6.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1 (form=~1|ID/treatment), method="REML", data=speed.data) 
mot85AS7.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation=corAR1 (form=~1|ID/wellcond), method="REML", data=speed.data) 

anova(mot85AS1.lme, mot85AS2.lme, mot85AS4.lme, mot85AS5.lme, mot85AS5a.lme, mot85AS5d.lme, mot85AS5e.lme, mot85AS6.lme, mot85AS7.lme) #

#models5-7 are the same

#other correlation structures

cs1 <- corARMA(c(0.2, 0.2), p = 2, q = 0)

mot85AS11.lme <- lme (Form, random = ~1|ID,  weights=varIdent(form=~1|treatment), correlation = cs1,
                 method="REML", data=speed.data) 

anova(mot85AS1.lme, mot85AS2.lme, mot85AS4.lme, mot85AS5.lme, mot85AS5a.lme, mot85AS5d.lme, mot85AS5e.lme, mot85AS6.lme, mot85AS7.lme, mot85AS11.lme) #

#summaries

summary(mot85AS4.lme)
anova(mot85AS4.lme)


#lowest AIC is mot85AS4.lme=875.1697

library(multcomp)
summary(glht(mot85AS4.lme, linfct=mcp(treatment="Tukey", covariate_average=TRUE)))

#other comparison methods, same results
summary(glht(mot85AS4.lme, linfct=mcp(treatment="Tukey", covariate_average=TRUE)), test = adjusted(type = "bonferroni"))

library(lsmeans)
mot85AS.lsmeans=lsmeans(mot85AS4.lme,  ~ treatment)
summary(mot85AS.lsmeans, adjust="bonf")
contrast(mot85AS.lsmeans, adjust="bonf")
pairs(mot85AS.lsmeans)

library(multcompView)
cld(mot85AS.lsmeans, alpha=0.05)

#alternative to mixed model
kruskal.test(Vlog ~ treatment, data = speed.data) 
dunnTest(Vlog ~ treatment, data = speed.data, method="none")  


#graphing

speed.data$treatment2 <- factor(speed.data$treatment, levels=c("ASW_induced", "ASW_notinduced", "Si_ASW", "Si_induced",  "Si_notinduced", "Si_Si"), 
                                                               labels=c("non-starved,\n induced", "non-starved,\n not induced", 
                                                                        "dSi-starved\n+dSi", "dSi-starved,\n induced",
                                                                        "dSi-starved,\n not induced", "dSi-starved\n+blank addition"))


speed.sumV2 <- summarySE(speed.data, measurevar="V", groupvars=c("treatment2"))

speed.sumV2$treatment2 <- factor(speed.sumV2$treatment2, levels = c("non-starved,\n induced", "non-starved,\n not induced",
                                                                          "dSi-starved,\n induced", "dSi-starved,\n not induced", 
                                                                          "dSi-starved\n+dSi", "dSi-starved\n+blank addition"))

qplot(treatment2, V, color = treatment2, data = speed.sumV2,  geom = "boxplot")

#sigbars are dependent on the multcomp results

speed.sumV2$sigbars <- c ("1", "1", "2", "2", "3", "4")

#bwPalette <- c("#000000", "#666666", "#cccccc", "#ffffff")

resize.win(12,9)

grid.newpage()
text <- element_text(size = 15) #change the size of the axes
theme_set(theme_bw()) 

ggplot(speed.sumV2, aes(treatment2, V, fill = sigbars)) + 
  geom_bar(stat="identity", position = "dodge", width=0.8, color="black") + geom_errorbar(aes(ymin=V-se, ymax=V+se), width=0.2, position=position_dodge(0.8))+
  scale_fill_grey(start = 0.3, end = 1) + 
  labs(y = "Speed (µm/s)") +
  theme(axis.text=element_text(size=15), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_blank(),
        plot.title = element_text(size =15, face="bold"), axis.text=text,  legend.position="none", legend.title = element_blank(),
        strip.text.x = element_text(size=15), strip.text.y = text, legend.title=text, legend.text=element_text(size=12), legend.direction="horizontal", 
        panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "cm")) + 
  guides(fill=guide_legend(keywidth=0.2,keyheight=0.2, default.unit="inch"))