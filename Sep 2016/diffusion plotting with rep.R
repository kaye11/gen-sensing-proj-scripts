library(gdata)
library(data.table)
library(ggplot2)
library(grid)
library(gtable)
library(plyr)

'rep1' <- read.csv("D:/Karen's/PhD/R program/General sensing proj/Processed_data/diffusion/rep1.csv", sep=";")
'rep2' <- read.csv("D:/Karen's/PhD/R program/General sensing proj/Processed_data/diffusion/rep2.csv", sep=";")
'rep3' <- read.csv("D:/Karen's/PhD/R program/General sensing proj/Processed_data/diffusion/rep3.csv", sep=";")

all <- rbind((rep1 [, c (7, 8, 9)]), (rep2 [, c (7, 8, 9)]), 
             (rep3 [, c (7, 8, 9)]))

dev.new(width=6, height=9)
source("resizewin.R")
resize.win(6,9)

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw())

ggplot (data=all, aes(x=rad2, y=uMsq))+geom_line(size=2)+ 
  labs(x="Distance from bead (µm)", 
       y="µM dP/bead")+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25,face="bold"), 
        plot.title = element_text(size =25, face="bold"), axis.text=text, 
        axis.title.y = element_text(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

source("summarySE.R")
difall<- summarySE(all, measurevar="uMsq", groupvars=c("rad2"), na.rm=TRUE)

resize.win (8,8)

ggplot(data=difall, aes(x=rad2, y=uMsq)) + geom_line(size=1)+ 
  geom_ribbon(aes(ymin=uMsq-sd, ymax=uMsq+sd), alpha=0.5, size=3) +  
  labs(x="Distance from bead (µm)", 
       y="µm dP/bead")+
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20, vjust=1.5), 
        axis.title.x=element_text(size=20, vjust=-0.5),
        plot.title = element_text(size =24), axis.text=text,  legend.position="bottom", legend.title=element_blank(),
        legend.key.width=unit(2.5,"cm"),legend.key.height=unit(0.8,"cm"), 
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))  


#for poster

ggplot(data=difall, aes(x=rad2, y=uMsq)) + geom_line(size=1)+ 
  geom_ribbon(aes(ymin=uMsq-sd, ymax=uMsq+sd), alpha=0.5, size=3) +  
  labs(x="Distance from bead (µm)", 
       y="µm dP/bead")+ 
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.position="none", legend.direction="horizontal",
        legend.title=element_blank(),legend.key.width=unit(1.8,"cm"),legend.key.height=unit(0.8,"cm"),  
        strip.text.x = text, strip.text.y = text, legend.title=text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm"))