
speed.sumV <- summarySE(speed.data, measurevar="V", groupvars=c("treatment"))

speed.sumV <- data.table (speed.sumV)
speed.sumV [, grouping := ifelse(treatment %in% c("nonstarved-motility", "starved-motility"), "before medium \nexchange",
                                 ifelse(treatment %in% c("nonstarved-addASW", "starved-addASW", "starved-addminusP", "starved-addminusSi"), 
                                        "1h after medium exchange", NA))]


speed.sumV$grouping2 <- c ("dP replete medium", "dP replete medium", "dP deplete medium", "blank addition", "dP deplete medium+ASW-dSi", 
                           "dP deplete medium")

speed.sumV$groupre <- factor(speed.sumV$grouping, levels=c("before medium \nexchange", "1h after medium exchange"))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

grid.newpage()
text <- element_text(size = 15 #change the size of the axes
theme_set(theme_bw()) 

ggplot(speed.sumV, aes(treatment, V, fill = grouping2)) + 
  geom_bar(stat="identity", position = "dodge", width=0.8) + geom_errorbar(aes(ymin=V-se, ymax=V+se), width=0.2, position=position_dodge(0.8))+
  facet_grid(.~groupre, scale="free_x", space="free") +
  scale_fill_manual(values = cbPalette) + scale_x_discrete(breaks=NULL) +
  labs(y = "Speed (µm/s)") +
  theme(axis.text=element_text(size=15), axis.title.y=element_text(size=20,face="bold", vjust=1.5), 
        axis.title.x=element_blank(),
        plot.title = element_text(size =15, face="bold"), axis.text=text,  legend.position="bottom", legend.title = element_blank(),
        strip.text.x = element_text(size=15), strip.text.y = text, legend.title=text, legend.text=element_text(size=12), legend.direction="horizontal", 
        panel.margin=unit (0.5, "lines"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(1,1,1,1), "cm")) + 
  guides(fill=guide_legend(keywidth=0.2,keyheight=0.2, default.unit="inch"))

