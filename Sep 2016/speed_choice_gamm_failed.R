
#GAMM: induced.sum

qplot(as.factor(timemin), Vlog, color = bin, data = induced.sum,  geom = "boxplot") + facet_grid(bin~., scales="free") +
  stat_smooth (method="loess", formula=y~x, size=1, aes(group=1))

exp.induced=as.data.frame(data.table(cbind(wellvid=as.numeric(as.factor(induced.sum$wellvid)),bin=as.numeric(as.factor(induced.sum$bin)), 
                                           T=as.numeric(induced.sum$timemin), ID=as.numeric(as.factor(induced.sum$ID)))))

cor(exp.induced, method = "spearman")

vif_func(in_frame=exp.induced,thresh=5,trace=T)

pairs(exp.induced, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

#boxplots
op=par(mfrow=c(2,2))
boxplot(Vlog~wellvid, data=induced.sum)
boxplot(Vlog~bin, data=induced.sum)
boxplot(Vlog~ID, data=induced.sum)
boxplot (Vlog~timemin, data=induced.sum)

#levene
library(lawstat)
levene.test(induced.sum$Vlog, group=induced.sum$wellvid, location="mean") #significant
levene.test(induced.sum$Vlog, group=induced.sum$ID, location="mean") # significant
levene.test(induced.sum$Vlog, group=induced.sum$timemin, location="mean")  # significant
levene.test(induced.sum$Vlog, group=induced.sum$bin, location="mean") # significant


#gamm
induced.gam <- gamm (Vlog~s(timemin, by=bin, bs="fs"), method="REML", data = induced.sum)

#smoothing splines for factors can are 3 types cr (cubic regression), cs (shrinkage version of cr) and cc (cyclic cubic)
induced.gam1 <- gamm (Vlog~s(timemin, by=bin, bs="fs", xt="cr"), method="REML", data = induced.sum) 
induced.gam2 <- gamm (Vlog~s(timemin, by=bin, bs="fs", xt="cs"), method="REML", data = induced.sum) #BEST
induced.gam2.1 <- gamm (Vlog~s(timemin, by=bin, bs="fs", xt="cc"), method="REML", data = induced.sum) 
anova(induced.gam$lme, induced.gam1$lme, induced.gam2$lme, induced.gam2.1$lme)

#make random factor and correlations
finduced.sum <- Vlog~s(timemin, by=bin, bs="fs", xt="cr")

induced.gam3 <- gamm (finduced.sum, method="REML",  random=list(ID=~1), data = induced.sum) 
induced.gam4 <- gamm (finduced.sum, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|bin/ID), data = induced.sum) #BEST
induced.gam5 <- gamm (finduced.sum, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|wellvid/ID), data = induced.sum) 
induced.gam6 <- gamm (finduced.sum, method="REML", random=list(ID=~1), correlation= corAR1 (), data = induced.sum) #same with model5

anova(induced.gam$lme, induced.gam1$lme, induced.gam2$lme, induced.gam2.1$lme, induced.gam3$lme, induced.gam4$lme, induced.gam5$lme, induced.gam6$lme)
#model4 is the best but 5 or 6 is the correct formulation

#make variance structures using model4
#induced.gam6 <- gamm (finduced.sum, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|bin/ID), weights = varIdent(form=~1| timemin), data = induced.sum) #no convergence

#induced.gam7 <- gamm (finduced.sum, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|bin/ID), weights = varIdent(form=~1| ID), data = induced.sum) #no convergence

induced.gam8 <- gamm (finduced.sum, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|bin/ID), 
             weights = varIdent(form=~1| bin), data = induced.sum) 

#induced.gam9 <- gamm (finduced.sum, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|bin/ID), weights = varIdent(form=~1| wellvid), data = induced.sum) 

induced.gam10 <- gamm (finduced.sum, method="REML", random=list(ID=~1), correlation= corAR1 (), 
                      weights = varIdent(form=~1| bin), data = induced.sum) 

anova(induced.gam$lme, induced.gam1$lme, induced.gam2$lme, induced.gam2.1$lme, induced.gam3$lme, induced.gam4$lme, induced.gam5$lme, induced.gam8$lme,
      induced.gam10$lme)

#model 10 is the best. lme AIC=  4709.980

#Take induced.gam10
with(induced.sum, tsDiagGamm(induced.gam10, timevar=timemin, observed=Vlog)) #visual checking of residuals

plot(induced.gam10$lme)

summary(induced.gam10$gam)
summary(induced.gam10$lme)
anova(induced.gam10$gam)
anova(induced.gam10$lme)

#checking the model produced by GAMM
resize.win(9,12)
op=par(mfrow=c(3,1), mar=c(4.5,4.5,1.5,1.5))
#plot(induced.gam8$gam, cex.lab=1.1, cex.axis=1.1, xlab ="timemin (s)")
plot(induced.gam10$gam, cex.lab=1.5, cex.axis=1.5, xlab ="timemin (min)")

#try lme

induced.lme1 <- lme(Vlog~as.factor(timemin)*bin, method="REML", random=list(ID=~1), correlation= corAR1 (), 
                    weights = varIdent(form=~1| bin), data = induced.sum)

anova(induced.lme1)
summary(induced.lme1)
library(lsmeans)

pairs(lsmeans(induced.lme1, ~bin|timemin))


##switch to lme and do sequential model building
