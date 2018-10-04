#comparing different starvation days (dSi only)

starvedaybase.si <- subset (starvedaybase, starvedaybase$treatment=="dSi")

expsi=as.data.frame(data.table(cbind(Bin=starvedaybase.si$Bin,T=starvedaybase.si$T, starvedays=starvedaybase.si$starvedays, 
                                     ID=starvedaybase.si$reptreatdays)))
cor(expsi, method = "spearman")


vif_func(in_frame=expsi,thresh=5,trace=T)

pairs(expsi, lower.panel = panel.smooth2,  upper.panel = panel.cor, diag.panel = panel.hist)

levene.test(starvedaybase.si$CellsBase, group=starvedaybase.si$reptreatdays, location="mean") #significant
levene.test(starvedaybase.si$CellsBase, group=starvedaybase.si$T, location="mean") # significant
levene.test(starvedaybase.si$CellsBase, group=starvedaybase.si$Bin, location="mean") #significant
levene.test(starvedaybase.si$CellsBase, group=starvedaybase.si$starveday, location="mean") #not significant

#boxplots
op=par(mfrow=c(2,2))
boxplot(CellsBase~Bin, data=starvedaybase.si)
boxplot(CellsBase~reptreatdays, data=starvedaybase.si)
boxplot (CellsBase~T, data=starvedaybase.si)
boxplot (CellsBase~starvedays, data=starvedaybase.si)


#fit a gls
Form <- formula (CellsBase ~ T*starvedays*Bin)
starvedaybase.si.gls<- gls(Form, data=starvedaybase.si)


#####NLME models######


#random structures #winner: starvedaybase.si1.lme
starvedaybase.si1.lme <- lme (Form, random = ~1|reptreatdays, method="REML", data=starvedaybase.si) 

starvedaybase.si2.lme <- lme (Form, random = ~1|treatment/reptreatdays, method="REML", data=starvedaybase.si)

anova(starvedaybase.si.gls, starvedaybase.si1.lme, starvedaybase.si2.lme)


#random structures+correlation structures #all produce the same AIC values even starvedaybase.si5.lme AIC=271.61

starvedaybase.si3.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (form=~1|reptreatdays/treatment), method="REML", data=starvedaybase.si) 

starvedaybase.si4.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), method="REML", data=starvedaybase.si) 

starvedaybase.si5.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (form=~1|reptreatdays/starvedays), method="REML", data=starvedaybase.si) 

anova(starvedaybase.si.gls, starvedaybase.si1.lme, starvedaybase.si2.lme, starvedaybase.si3.lme, starvedaybase.si4.lme, starvedaybase.si5.lme)


#random structures+correlation structures+weights #starvedaybase.si10.lme wins AIC=-100.29

#starvedaybase.si6.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varIdent(form=~1|reptreatdays), method="REML", data=starvedaybase.si) 

starvedaybase.si7.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varIdent(form=~1|starvedays), method="REML", data=starvedaybase.si) 

starvedaybase.si8.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varIdent(form=~1|Bin), method="REML", data=starvedaybase.si) 

anova(starvedaybase.si.gls, starvedaybase.si1.lme, starvedaybase.si2.lme, starvedaybase.si3.lme, starvedaybase.si4.lme, starvedaybase.si5.lme, 
      starvedaybase.si7.lme, starvedaybase.si8.lme)

starvedaybase.si9.lme <- lme (Form, random = ~1|reptreatdays, correlation=corAR1 (), weights=varComb(varIdent(form=~1|Bin), varIdent(form=~1|starvedays)), 
                   method="REML", data=starvedaybase.si)

anova(starvedaybase.si.gls, starvedaybase.si1.lme, starvedaybase.si2.lme, starvedaybase.si3.lme, starvedaybase.si4.lme, starvedaybase.si5.lme, 
      starvedaybase.si7.lme, starvedaybase.si8.lme, stravedaybase.si9.lme)


#let's plot this!

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 
library (AICcmodavg)



#starvedaybase.si fit
starvedaybase.si.sum <- summarySE(starvedaybase.si, measurevar="CellsBase", groupvars=c("T", "starvedays","Bin"))

starvedaybase.si.fit <- as.data.frame(predictSE.lme(starvedaybase.si9.lme, starvedaybase.si, se.fit = TRUE, level = 0,
                                                    print.matrix = FALSE))

starvedaybase.si.fit$upr <- starvedaybase.si.fit$fit + (1.96 * starvedaybase.si.fit$se)
starvedaybase.si.fit$lwr <- starvedaybase.si.fit$fit - (1.96 * starvedaybase.si.fit$se)

starvedaybase.si.fit.combdata <- cbind(starvedaybase.si, starvedaybase.si.fit)

ggplot(data=starvedaybase.si.sum, aes(x=T, y=CellsBase, shape=starvedays, color=starvedays)) + geom_point(size=5)+
  geom_errorbar(aes(ymin=CellsBase-se, ymax=CellsBase+se), width=30, size=1) +
  geom_smooth(data=starvedaybase.si.fit.combdata, size=1,  aes(y=fit, ymin=lwr, ymax=upr, fill=starvedays), method="lm", stat="identity", alpha=0.1)+ 
  scale_shape_discrete (name="starvedays") +
  scale_fill_discrete(name="starvedays") + facet_wrap (Bin~starvedays)+
  labs(list(x = "Time (s)", y = "Normalized cell count", title="dSi only"))+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_blank(), 
        axis.title.x=element_blank(),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="bottom",
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous (breaks=c(200, 400, 600)) 
