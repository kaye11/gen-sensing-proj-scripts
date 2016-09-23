induced_DPR <- read.csv("D:/Karen's/PhD/R program/General sensing proj/Processed_data/RMS data_choice/induced_DPR .csv", sep=";")
notinduced_DPR <- read.csv("D:/Karen's/PhD/R program/General sensing proj/Processed_data/RMS data_choice/notinduced_DPR .csv", sep=";")
induced_dSi <- read.csv("D:/Karen's/PhD/R program/General sensing proj/Processed_data/RMS data_choice/induced_dSi .csv", sep=";")
notinduced_dSi <- read.csv("D:/Karen's/PhD/R program/General sensing proj/Processed_data/RMS data_choice/notinduced_dSi .csv", sep=";")

induced_DPR$induction = as.factor("induced")
induced_dSi$induction = as.factor("induced")
notinduced_DPR$induction = as.factor("not induced")
notinduced_dSi$induction = as.factor("not induced")

induced_DPR$bead = as.factor("Diproline bead")
notinduced_DPR$bead = as.factor("Diproline bead")
induced_dSi$bead = as.factor("dSi bead")
notinduced_dSi$bead = as.factor("dSi bead")

library(ggplot2)
library(grid)
library(ggthemes)
library(gridExtra)
library(reshape2)

dev.new(width=6, height=9)
source("resizewin.R")

grid.newpage()
text <- element_text(size = 18, face="bold") #change the size of the axes
theme_set(theme_bw()) 

##RMS plotting overlayed
RMS <- rbind (induced_dSi, notinduced_dSi, induced_DPR, notinduced_DPR)

RMS$timemin <- RMS$time/60


#plot

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 

scaleFUN <- function(x) sprintf("%.1f", x)

resize.win(6,9)

ggplot(RMS, aes(x=timemin, y = MF, linetype=bead)) + geom_line(size=2)+ facet_grid(induction~.) +
  labs(list(x = "Time (s)", y = "RMS (µm)")) + 
  scale_linetype_manual(values=c("dashed", "solid")) +
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20,face="bold", vjust=1.5),
        axis.title.x=element_text(size=20,face="bold", vjust=-0.5),
        plot.title = element_text(size =20, face="bold"), axis.text=text,  legend.title=element_blank(), legend.text=text, 
        legend.key.width=unit(2.5,"cm"),legend.key.height=unit(0.8,"cm"), legend.position="bottom", 
        strip.text = text, legend.title=text, legend.text=text, panel.margin=unit (1.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_x_continuous (breaks=c(0, 2, 4, 6, 8, 10))+
  scale_y_continuous(labels=scaleFUN)


