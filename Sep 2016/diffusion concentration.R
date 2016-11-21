library(gdata)
library(data.table)
library(ggplot2)
library(grid)
library(gtable)

ia= 5.7004E-08
ib= 5.52766E-08
ic= 5.61083E-08
imean= 5.61296E-08


bead2=55*0.0001 #in cm
rad=seq(1, 340, 1)*0.0001
D= 10^-5 #cm^2/s
i=ia*10^-6 #in mol/s
Ci=i/4*pi*bead2*D #in mol/cm3
i2=Ci*4*pi*bead2*D


var2=as.data.frame(rad)

var2$mcm=i/(4*pi*D*var2$rad)
var2$mcmsq=i/(4*pi*D*sqrt(var2$rad))
var2$M=var2$mcm*10^3
var2$Msq=var2$mcmsq*10^3
var2$uM=var2$M*10^6
var2$uMsq=var2$Msq*10^6
var2$rad2=var2$rad/0.0001
var2$ts=(2*var2$rad)^2/D
var2$nMsq=var2$Msq*10^9
var2$log=log10(var2$uMsq)
var2$fMsq=var2$nMsq*1000000


var=subset(var2, var2$rad2>55, )

grid.newpage()
text <- element_text(size = 20, face="bold") #change the size of the axes
theme_set(theme_bw())

ggplot (data=var2, aes(x=rad2, y=uMsq))+geom_line(size=2)+ xlab("Distance from bead (µm)")+
  ylab(expression(paste("µM ", dP / bead, sep="")))+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25), 
        plot.title = element_text(size =25, face="bold"), axis.text=text, 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

write.table (var2, "d:/Karen's/PhD/R program/General sensing proj/Processed_data/diffusion/rep3.csv", 
             sep=";", col.names=T, row.names=F)


