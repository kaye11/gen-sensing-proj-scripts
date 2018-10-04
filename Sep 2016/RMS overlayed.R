library(readr)

induced_DPR <- read_delim("D:/Karen's/PhD/R program/General sensing proj/Processed_data/RMS data_choice/induced_DPR .csv", 
                          ";", escape_double = FALSE, trim_ws = TRUE)
notinduced_DPR <- read_delim("D:/Karen's/PhD/R program/General sensing proj/Processed_data/RMS data_choice/notinduced_DPR .csv", 
                          ";", escape_double = FALSE, trim_ws = TRUE)

induced_dSi <- read_delim("D:/Karen's/PhD/R program/General sensing proj/Processed_data/RMS data_choice/induced_dSi .csv", 
                          ";", escape_double = FALSE, trim_ws = TRUE)
notinduced_dSi <- read_delim("D:/Karen's/PhD/R program/General sensing proj/Processed_data/RMS data_choice/notinduced_dSi .csv", 
                          ";", escape_double = FALSE, trim_ws = TRUE)

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

#dev.new(width=6, height=9)
source("resizewin.R")

grid.newpage()
text <- element_text(size = 18, face="bold") #change the size of the axes
theme_set(theme_bw()) 

##RMS plotting overlayed
RMS <- rbind (induced_dSi, notinduced_dSi, induced_DPR, notinduced_DPR)

RMS$timemin <- RMS$time/60
RMS$timeminsq <- sqrt(RMS$timemin)
RMS$MSD <- RMS$MF^2

#plot

grid.newpage()
text <- element_text(size = 20) #change the size of the axes
theme_set(theme_bw()) 

scaleFUN <- function(x) sprintf("%.1f", x)

RMS$bead2 <- factor(RMS$bead, levels=c("dSi bead", "Diproline bead"), 
                                    labels =c (" dSi bead  ", " Diproline bead  "))

#resize.win(12,9)

ggplot(RMS, aes(x=time, y = MF, linetype=bead2)) + geom_line(size=2)+ facet_grid(~induction) +
  labs(list(x = "Time (s)", y = "RMS (µm)", title="Differences in the root mean square travelled between beads")) + 
  scale_linetype_manual(values=c("dashed", "solid")) +
  theme(axis.text=element_text(size=20), axis.title.y=element_text(size=20, vjust=1.5), 
        axis.title.x=element_text(size=20, vjust=-0.5),
        plot.title = element_text(size =24), legend.position="bottom", legend.title=element_blank(),
        legend.key.width=unit(2.5,"cm"),legend.key.height=unit(0.8,"cm"), 
        strip.text.x = text, strip.text.y = text, legend.text=text, panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(),panel.margin.y = unit(1, "lines"), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_y_continuous(labels=scaleFUN)


