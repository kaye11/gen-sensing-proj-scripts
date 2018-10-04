trackdata <- trackdata [! trackdata$time=="0",  ]
trackdata <- trackdata [! trackdata$ID=="A3-004-0", ]
trackdata <- trackdata [! trackdata$ID=="A3-004-469", ]
trackdata <- trackdata [! trackdata$ID=="A3-004-946", ]


qplot(time2, Vlog, color = treatment, data = trackdata,  geom = "boxplot") + facet_grid(treatment~newbin, scales="free") 

#gamm
BA <- gamm (Vlog~s(time, by=treatment, bs="fs"), method="REML", data = trackdata)

#smoothing splines for factors can are 3 types cr (cubic regression), cs (shrinkage version of cr) and cc (cyclic cubic)
BA1 <- gamm (Vlog~s(time, by=treatment, bs="fs", xt="cr"), method="REML", data = trackdata) 
BA2 <- gamm (Vlog~s(time, by=treatment, bs="fs", xt="cs"), method="REML", data = trackdata) 
BA2.1 <- gamm (Vlog~s(time, by=treatment, bs="fs", xt="cc"), method="REML", data = trackdata) #best but use cr because it yields better results with corAR1

anova(BA$lme, BA1$lme, BA2$lme, BA2.1$lme)

#make random factor and correlations
ftrackdata <- Vlog~s(time, by=treatment, bs="fs", xt="cr")

BA3 <- gamm (ftrackdata, method="REML",  random=list(ID=~1), data = trackdata) 
BA4 <- gamm (ftrackdata, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), data = trackdata) #BEST
BA5 <- gamm (ftrackdata, method="REML", random=list(ID=~1), correlation= corAR1 (), data = trackdata) #same with BA4

anova(BA$lme, BA1$lme, BA2$lme, BA2.1$lme, BA3$lme, BA4$lme, BA5$lme)

#make variance structures
#BA6 <- gamm (ftrackdata, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varIdent(form=~1| time), data = trackdata) #no convergence

#BA7 <- gamm (ftrackdata, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), weights = varIdent(form=~1| ID), data = trackdata) #no convergence

BA8 <- gamm (ftrackdata, method="REML", random=list(ID=~1), correlation= corAR1 (form=~1|treatment/ID), 
             weights = varIdent(form=~1| treatment), data = trackdata)