##source this! name the file into s1
library (data.table)
library (ggplot2)
library(grid)

##direct input from trackmate
#tracking is every 2s so T is multiplied by 3 and V is divided by 2
NT=data.table(s1, key="TRACK_ID")
NT=NT[, list(A=TRACK_ID, X=POSITION_X, Y=POSITION_Y, T=FRAME), by=c("TRACK_ID")]
s1a= NT[order(A, T)]
s1=s1a[,list(T=T*2, X=X, Y=Y, V=c(0, (sqrt(diff(X)^2+diff(Y)^2))/2)), by=c("A")]


##parameter computations
beadX<-as.numeric(readline("X position of the Bead?"))
beadY<-as.numeric(readline("Y position of the Bead?"))
beadSize<-as.numeric(readline("radius of the Bead?"))

library(data.table)
library(zoo)

deg<-180/pi

t1=s1

NT<-data.table(t1, key=c("A"))
NT[, V := c(0, V[2:(.N-1)], 0), by = A]

t1=NT[,list(T=T, X=X, Y=Y, V=c(0, sqrt(diff(X)^2+diff(Y)^2)), GD=cumsum(V), ND=sqrt((X-X[1])^2 + (Y-Y[1])^2)), 
      by=c("A")]

t1[,NGDR:=ND/GD]
#t1[,ED:=sqrt((X-X[.N])^2 + (Y-Y[.N])^2), by=A]


t1[,a:=c(NA,(X[2:(.N-1)]),NA),by=A]
t1[,b:=c(NA,(Y[2:(.N-1)]),NA),by=A]
t1[,c:=c(NA,(X[1:(.N-2)]),NA),by=A]
t1[,d:=c(NA,(Y[1:(.N-2)]),NA),by=A]

t1[, scalar:=(a-c)*(a-beadX)+(b-d)*(b-beadY)]
t1[, det:=(a-c)*(a-beadX)-(b-d)*(b-beadY)]
t1[, angle:= atan2(det, scalar)*deg]
t1[, angle2:=c(NA, (na.locf(angle [1:(.N-1)])), NA), by=A]
t1[, angs:=sin(angle2*pi/180)]

t1[,a:=c(NA,(V[2:(.N)])),by=A]
t1[,b:=c(NA,(V[1:(.N-1)])),by=A]

t1[, Accel:=a-b]
t1$a=NULL
t1$b=NULL
t1$c=NULL
t1$d=NULL
t1$scalar=NULL
t1$det=NULL
t1$angle=NULL

#t1[, mag:= sqrt((a-c)^2+(b-d)^2)*sqrt((a-beadX)^2+(b-beadY)^2)]
#t1[, div:= scalar/mag]
#t1[, angleman:= acos(div)*deg]
#t1[, angle:=acos(((a-c)*(a-beadX)+(b-d)*(b-beadY))/(sqrt((a-c)^2+(b-d)^2)*sqrt((a-beadX)^2+(b-beadY)^2)))*deg]


#t1[, angc:=cos(angle2*pi/180)]
#t1[, angatan:=atan2(angs, angc)]

##saving data
VN<- readline("What data did you analyse?")
t1$video <- paste(VN)
Vid<-paste ("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate automatic/raw track data/",VN,".csv")
write.table(t1, Vid, sep=";", col.names=T, row.names=F)


##plotting
#plot is saved but make an evernote copy as well

motplot <- qplot(X, Y, data=t1, color = factor(A), group = factor(A))+ geom_point(size=0.5)+
  geom_path(lineend="round", size=1, arrow=arrow(angle=45, ends="last", type= "closed", length=unit (0.10, "inches")), 
            mapping=aes(group=factor(A))) + scale_y_reverse()+theme_classic() + 
  labs(list(title=paste(VN)))+
  theme(plot.title=element_text(size=15, face="bold"), axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        legend.position="none", plot.title = element_text(size = 15))+
  annotate("path",x = beadX + beadSize*cos(seq(0,2*pi,length.out=100)), y=beadY + beadSize*sin(seq(0,2*pi,length.out=100)), 
           color="black", size=2)

motplotfile <- paste ("d:/Karen's/PhD/R program/General sensing proj/tiff files/bead motility plots/automatic tracking/",VN,".tiff")

ggsave(motplot, filename=motplotfile, dpi=800, width=6, height=6, units="in")

motplot
