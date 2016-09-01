##for checking max. no of tracks per sec

#control
control$A=as.numeric(control$A)
control.count=ddply(control,~T,summarise,DPR=(length(unique((ID)))))

#P bead
Pbead$A=as.numeric(Pbead$A)
Pbead.count=ddply(Pbead,~T,summarise,DPR=(length(unique((ID)))))

#Si bead
Sibead$A=as.numeric(Sibead$A)
Sibead.count=ddply(Sibead,~T,summarise,DPR=(length(unique((ID)))))

#data binding (DPR=A, HLB=HLB)
counttracks=cbind(control.count, Pbead.count, Sibead.count)

ggplot(counttracks, aes(T)) + geom_line(aes(y = DPR, colour = "DPR"),  DPRze=1) + 
  geom_line(aes(y = HLB, colour = "HLBtrol"), DPRze=1)+ 
  labs(list(x = "Time (s)", y = "Number of Tracks")) 


